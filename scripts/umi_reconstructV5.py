#!/usr/bin/env python3
"""
UMI Reads 重构工具 v5.1
========================================
核心逻辑：
  - STARsolo BAM 中 CB/UB 已是校正+collapsing后结果，直接以 (CB, GX, UB) 为唯一键分组
  - read级别过滤：NH==1（唯一比对）& GX!='-'（比对到基因）& proper_pair
  - 每个 (CB, GX, UB) 组的所有reads做CIGAR重构（M/N/D区间合并）
  - 按染色体并行，多线程加速

使用方法：
  python umi_reconstructV5.py \
      -i input.bam \
      -o output_prefix \
      [--threads 10] \
      [--max_umi_reads 500] \
      [--single_end]

输出文件：
  <output_prefix>.tsv             每行一个重构UMI，含完整元信息，tab分隔
  <output_prefix>_summary.txt     统计摘要，可直接 cat 查看

输出字段（tsv列）：
  cb, ub, gx, gn, chrom, strand,
  n_reads,           该UMI组原始read数（PCR扩增深度参考）
  new_cigar_start,   重构后起始位置（1-based）
  new_cigar,         重构后CIGAR（M/N/D）
  n_conflict         M与N区间冲突碱基数（M优先）
"""

import pysam
import argparse
import sys
import re
import os
import time
import pandas as pd
from collections import defaultdict
from typing import List, Tuple, Dict, Optional
from joblib import Parallel, delayed


# ─────────────────────────────────────────────
# CIGAR 解析与区间提取
# ─────────────────────────────────────────────

def parse_cigar(cigar_str: str) -> List[Tuple[int, str]]:
    return [(int(l), op) for l, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar_str)]


def extract_intervals(start_0: int, cigar_ops: List[Tuple[int, str]]):
    """
    展开CIGAR为三类区间（0-based，左闭右开）：
      m_ivs : M/=/X（外显子匹配）
      n_ivs : N（intron跳跃）
      d_ivs : D（deletion）
    """
    m_ivs, n_ivs, d_ivs = [], [], []
    pos = start_0
    for length, op in cigar_ops:
        if op in ('M', '=', 'X'):
            m_ivs.append((pos, pos + length))
            pos += length
        elif op == 'N':
            n_ivs.append((pos, pos + length))
            pos += length
        elif op == 'D':
            d_ivs.append((pos, pos + length))
            pos += length
    return m_ivs, n_ivs, d_ivs


# ─────────────────────────────────────────────
# 区间集合运算
# ─────────────────────────────────────────────

def union_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    if not intervals:
        return []
    ivs = sorted(intervals)
    merged = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [tuple(x) for x in merged]


def subtract_intervals(base: List[Tuple[int, int]],
                       minus: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    result = list(base)
    for ms, me in sorted(minus):
        new_result = []
        for s, e in result:
            if me <= s or ms >= e:
                new_result.append((s, e))
            else:
                if s < ms:
                    new_result.append((s, ms))
                if me < e:
                    new_result.append((me, e))
        result = new_result
    return result


def intersect_size(a: List[Tuple[int, int]],
                   b: List[Tuple[int, int]]) -> int:
    total = 0
    i, j = 0, 0
    a, b = sorted(a), sorted(b)
    while i < len(a) and j < len(b):
        lo = max(a[i][0], b[j][0])
        hi = min(a[i][1], b[j][1])
        if lo < hi:
            total += hi - lo
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return total


def intervals_between(ref_union: List[Tuple[int, int]],
                      skip_union: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """ref和skip都不覆盖的空隙 → D区间"""
    if not ref_union:
        return []
    all_covered = union_intervals(ref_union + skip_union)
    ref_min = ref_union[0][0]
    ref_max = ref_union[-1][1]
    gaps = []
    prev = ref_min
    for s, e in all_covered:
        if s > prev and prev >= ref_min:
            gaps.append((prev, min(s, ref_max)))
        prev = max(prev, e)
    if prev < ref_max:
        gaps.append((prev, ref_max))
    return [(s, e) for s, e in gaps if s < e]


# ─────────────────────────────────────────────
# 新 CIGAR 生成
# ─────────────────────────────────────────────

def build_cigar(ref_union: List[Tuple[int, int]],
                skip_union: List[Tuple[int, int]],
                del_union: List[Tuple[int, int]]) -> Tuple[int, str, int]:
    """
    三类区间 → 新CIGAR字符串 + 起始位点（1-based）
    冲突处理：M优先，从skip中减去与ref重叠部分
    """
    n_conflict = intersect_size(ref_union, skip_union)
    if n_conflict > 0:
        skip_union = subtract_intervals(skip_union, ref_union)

    del_union = intervals_between(ref_union, skip_union)

    labeled = (
        [(s, e, 'M') for s, e in ref_union] +
        [(s, e, 'N') for s, e in skip_union] +
        [(s, e, 'D') for s, e in del_union]
    )
    if not labeled:
        return 0, '', 0

    labeled.sort(key=lambda x: x[0])
    new_start_1 = labeled[0][0] + 1
    cigar = ''.join(f"{e - s}{op}" for s, e, op in labeled if e - s > 0)
    return new_start_1, cigar, n_conflict


# ─────────────────────────────────────────────
# 辅助：读取BAM标签
# ─────────────────────────────────────────────

def get_tag(read: pysam.AlignedSegment, tag: str, default: str = '-') -> str:
    try:
        return str(read.get_tag(tag))
    except KeyError:
        return default


# ─────────────────────────────────────────────
# UMI 组重构
# ─────────────────────────────────────────────

def reconstruct_umi_group(reads_info: List[Dict]) -> Optional[Dict]:
    """
    对一个 (CB, GX, UB) 组的所有reads做CIGAR重构。
    返回重构结果dict，或None（丢弃）。
    """
    if not reads_info:
        return None

    has_r1    = any(r['is_read1'] for r in reads_info)
    is_paired = reads_info[0].get('is_paired', True)
    if is_paired and not has_r1:
        return None

    cb = reads_info[0]['cb']
    ub = reads_info[0]['ub']
    gx = reads_info[0]['gx']
    gn = next((r['gn'] for r in reads_info if r['gn'] != '-'), '-')

    # 染色体多数投票
    chrom_count: Dict[str, int] = {}
    for r in reads_info:
        chrom_count[r['chrom']] = chrom_count.get(r['chrom'], 0) + 1
    chrom = max(chrom_count, key=chrom_count.get)
    reads_info = [r for r in reads_info if r['chrom'] == chrom]

    all_m, all_n, all_d = [], [], []
    for r in reads_info:
        all_m.extend(r['m_ivs'])
        all_n.extend(r['n_ivs'])
        all_d.extend(r['d_ivs'])

    ref_union  = union_intervals(all_m)
    skip_union = union_intervals(all_n)
    del_union  = union_intervals(all_d)

    if not ref_union:
        return None

    new_start_1, new_cigar, n_conflict = build_cigar(ref_union, skip_union, del_union)
    if not new_cigar:
        return None

    r1_strands = [r['strand'] for r in reads_info if r['is_read1']]
    if r1_strands:
        strand = max(set(r1_strands), key=r1_strands.count)
    else:
        all_strands = [r['strand'] for r in reads_info]
        strand = max(set(all_strands), key=all_strands.count)

    return {
        'cb':              cb,
        'ub':              ub,
        'gx':              gx,
        'gn':              gn,
        'chrom':           chrom,
        'strand':          strand,
        'n_reads':         len(reads_info),
        'new_cigar_start': new_start_1,
        'new_cigar':       new_cigar,
        'n_conflict':      n_conflict,
    }


# ─────────────────────────────────────────────
# 单条染色体处理（多线程worker）
# ─────────────────────────────────────────────

OUTPUT_COLS = ['cb', 'ub', 'gx', 'gn', 'chrom', 'strand',
               'n_reads', 'new_cigar_start', 'new_cigar', 'n_conflict']


def process_chrom(bam_path:      str,
                  chrom:         str,
                  max_umi_reads: int,
                  single_end:    bool) -> Tuple[List[Dict], Dict]:
    """处理单条染色体，返回 (results列表, stats字典)"""
    results = []
    stats = {
        'chrom':            chrom,
        'n_reads_total':    0,
        'skip_multimapper': 0,
        'skip_no_gene':     0,
        'skip_no_umi':      0,
        'skip_no_cb':       0,
        'skip_unmapped':    0,
        'skip_improper':    0,
        'skip_secondary':   0,
        'skip_large_umi':   0,
        'skip_no_r1':       0,
        'skip_no_ref':      0,
        'n_umi_out':        0,
    }

    umi_groups: Dict[str, List[Dict]] = defaultdict(list)

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for read in bam.fetch(chrom):
            stats['n_reads_total'] += 1

            if read.is_unmapped or read.cigarstring is None:
                stats['skip_unmapped'] += 1
                continue
            if read.flag & 0x900:
                stats['skip_secondary'] += 1
                continue
            if not single_end and not read.is_proper_pair:
                stats['skip_improper'] += 1
                continue

            try:
                nh = read.get_tag('NH')
            except KeyError:
                nh = 1
            if nh != 1:
                stats['skip_multimapper'] += 1
                continue

            try:
                gx = str(read.get_tag('GX'))
            except KeyError:
                gx = '-'
            if gx == '-' or gx == '':
                stats['skip_no_gene'] += 1
                continue

            try:
                ub = str(read.get_tag('UB'))
            except KeyError:
                ub = '-'
            if ub == '-' or ub == '':
                stats['skip_no_umi'] += 1
                continue

            try:
                cb = str(read.get_tag('CB'))
            except KeyError:
                cb = '-'
            if cb == '-' or cb == '':
                stats['skip_no_cb'] += 1
                continue

            cigar_ops = parse_cigar(read.cigarstring)
            start_0   = read.reference_start
            m_ivs, n_ivs, d_ivs = extract_intervals(start_0, cigar_ops)

            umi_key = f"{cb}__{gx}__{ub}"
            umi_groups[umi_key].append({
                'is_read1':  read.is_read1,
                'is_paired': bool(read.flag & 0x1),
                'chrom':     chrom,
                'strand':    '-' if read.is_reverse else '+',
                'cb':        cb,
                'ub':        ub,
                'gx':        gx,
                'gn':        get_tag(read, 'GN'),
                'm_ivs':     m_ivs,
                'n_ivs':     n_ivs,
                'd_ivs':     d_ivs,
            })

    # 重构所有UMI组
    for umi_key, reads_list in umi_groups.items():
        if len(reads_list) > max_umi_reads:
            stats['skip_large_umi'] += 1
            continue

        result = reconstruct_umi_group(reads_list)
        if result is None:
            has_r1 = any(r['is_read1'] for r in reads_list)
            if not has_r1:
                stats['skip_no_r1'] += 1
            else:
                stats['skip_no_ref'] += 1
            continue

        results.append(result)
        stats['n_umi_out'] += 1

    print(f"  [{chrom}] reads={stats['n_reads_total']:,}  "
          f"UMI输出={stats['n_umi_out']:,}  "
          f"过滤(无基因)={stats['skip_no_gene']:,}  "
          f"过滤(多重比对)={stats['skip_multimapper']:,}",
          flush=True)

    return results, stats


# ─────────────────────────────────────────────
# 主流程
# ─────────────────────────────────────────────

def process_bam(bam_path:      str,
                output_prefix: str,
                threads:       int  = 10,
                max_umi_reads: int  = 500,
                single_end:    bool = False):

    out_tsv     = output_prefix + '.tsv'
    out_summary = output_prefix + '_summary.txt'

    t0 = time.time()

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        chroms = [sq['SN'] for sq in bam.header['SQ']]

    print(f"[INFO] 染色体数 : {len(chroms)}", flush=True)
    print(f"[INFO] 线程数   : {threads}", flush=True)
    print(f"[INFO] 输入BAM  : {bam_path}", flush=True)
    print(f"[INFO] 输出前缀 : {output_prefix}", flush=True)
    print(flush=True)

    # 按染色体并行
    outputs = Parallel(n_jobs=1, verbose=0, backend='threading')(
        delayed(process_chrom)(bam_path, chrom, max_umi_reads, single_end)
        for chrom in chroms
    )

    # 汇总
    all_results = []
    total_stats: Dict = defaultdict(int)
    chrom_umi_counts = {}

    for results, stats in outputs:
        all_results.extend(results)
        for k, v in stats.items():
            if k != 'chrom':
                total_stats[k] += v
        chrom_umi_counts[stats['chrom']] = stats['n_umi_out']

    # 写 tsv
    if all_results:
        df = pd.DataFrame(all_results, columns=OUTPUT_COLS)
        df.to_csv(out_tsv, sep='\t', index=False)
        print(f"\n[INFO] 已保存: {out_tsv}  ({len(df):,} 行)", flush=True)
    else:
        print("\n[WARNING] 没有任何重构结果，请检查过滤参数或BAM标签", flush=True)

    elapsed = time.time() - t0

    # 写 summary.txt
    with open(out_summary, 'w') as f:
        def w(line=''):
            print(line, file=f)
            print(line, flush=True)

        w('=' * 58)
        w('  UMI 重构结果统计  (umi_reconstructV5.1)')
        w('=' * 58)
        w(f'  输入BAM      : {bam_path}')
        w(f'  输出TSV      : {out_tsv}')
        w(f'  线程数       : {threads}')
        w(f'  运行时间     : {elapsed:.1f} s  ({elapsed/60:.1f} min)')
        w()
        w('── Read 级别过滤 ──────────────────────────────────')
        w(f'  输入reads总数          : {total_stats["n_reads_total"]:>12,}')
        w(f'  过滤-未比对            : {total_stats["skip_unmapped"]:>12,}')
        w(f'  过滤-次要/补充比对     : {total_stats["skip_secondary"]:>12,}')
        w(f'  过滤-非proper pair     : {total_stats["skip_improper"]:>12,}')
        w(f'  过滤-多重比对(NH>1)    : {total_stats["skip_multimapper"]:>12,}')
        w(f'  过滤-未比对到基因(GX-) : {total_stats["skip_no_gene"]:>12,}')
        w(f'  过滤-无有效UMI         : {total_stats["skip_no_umi"]:>12,}')
        w(f'  过滤-无有效CB          : {total_stats["skip_no_cb"]:>12,}')
        w()
        w('── UMI 组级别过滤 ─────────────────────────────────')
        w(f'  丢弃-超大组(>{max_umi_reads}reads)   : {total_stats["skip_large_umi"]:>12,}')
        w(f'  丢弃-无Read1的组       : {total_stats["skip_no_r1"]:>12,}')
        w(f'  丢弃-无M区间的组       : {total_stats["skip_no_ref"]:>12,}')
        w()
        w('── 输出 ───────────────────────────────────────────')
        w(f'  输出UMI总数            : {total_stats["n_umi_out"]:>12,}')
        w()
        w('── 各染色体 UMI 数 ────────────────────────────────')
        for chrom in chroms:
            cnt = chrom_umi_counts.get(chrom, 0)
            if cnt > 0:
                w(f'  {chrom:<22}: {cnt:>10,}')
        w('=' * 58)

    print(f"\n[INFO] 已保存统计: {out_summary}", flush=True)


# ─────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='UMI reads重构工具 v5.1 — STARsolo BAM，按染色体多线程并行，输出TSV')
    parser.add_argument('-i', '--input',     required=True,
                        help='输入BAM文件（坐标排序，需有CB/UB/GX/NH标签）')
    parser.add_argument('-o', '--output',    required=True,
                        help='输出文件前缀（不含扩展名）。'
                             '生成 <prefix>.tsv 和 <prefix>_summary.txt')
    parser.add_argument('--threads',         type=int, default=10,
                        help='并行线程数，按染色体并行（默认10）')
    parser.add_argument('--max_umi_reads',   type=int, default=500,
                        help='同一UMI组最多允许reads数，超出丢弃（默认500）')
    parser.add_argument('--single_end',      action='store_true',
                        help='单端模式：关闭proper pair过滤（默认双端）')
    args = parser.parse_args()

    process_bam(
        bam_path=args.input,
        output_prefix=args.output,
        threads=args.threads,
        max_umi_reads=args.max_umi_reads,
        single_end=args.single_end,
    )


if __name__ == '__main__':
    main()
