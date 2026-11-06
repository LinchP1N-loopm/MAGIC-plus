"""
collect_editing_sites.py  —  适配 STARsolo / Patho-DBiT (CB-tag 版本)

与原始脚本的核心差异:
  原始: spot = read_name.rsplit('|')[-1].replace('_', 'x')
        依赖 read name 末尾携带 row_col 坐标

  本版: mpileup 输出的 QNAME 格式为普通 read name
        barcode 在 BAM 的 CB tag 里，但 mpileup --output-QNAME 只输出 read name
        
  解决方案:
        在运行 mpileup 之前，先用 samtools 将每条 read 的 CB tag 写入 read name
        (见配套 shell 脚本 RNA_editing_pipeline_adapted.sh)
        改造后 read name 末尾追加 |CB序列，例如:
          AE01-...:R021C006:0049:8865|GTACAGCGTCTACTGC
        本脚本再从 read name 取出 CB，查 cb_to_pixel 字典得到像素坐标
"""

import os
import sys
import re
import argparse
from collections import defaultdict as dd
import cgranges as cr

min_total_read_count = 10
min_edited_read_count = 0


# ── 新增：读取 CB → pixel 映射 ────────────────────────────────────────
def get_cb_to_pixel(cb_to_pixel_fn):
    """
    读取 cb_to_pixel.tsv
    格式: AACGCCTCAATCTCGA\t139x61
    返回: dict {CB_seq -> pixel_str}
    """
    cb_to_pixel = {}
    with open(cb_to_pixel_fn) as f:
        for line in f:
            ele = line.strip().split('\t')
            if len(ele) < 2:
                continue
            cb_to_pixel[ele[0]] = ele[1]
    sys.stderr.write(f'[INFO] Loaded {len(cb_to_pixel)} CB→pixel mappings\n')
    return cb_to_pixel


def get_spot_to_clu(spot_to_clu_fn, clu_to_name_tsv):
    spot_to_clu = dd(lambda: 'NotClustered')
    clu_to_name = dict()
    with open(spot_to_clu_fn) as spot_fp, open(clu_to_name_tsv) as clu_fp:
        for line in clu_fp:
            ele = line.strip().rsplit()
            if len(ele) < 2:
                continue
            clu, name = ele[0], ele[1]
            clu_to_name[clu] = name
        for line in spot_fp:
            ele = line.strip().rsplit()
            if len(ele) < 2:
                continue
            spot, clu = ele[0], ele[1]
            if clu not in clu_to_name:
                name = 'NoClusterName'
            else:
                name = clu_to_name[clu]
            spot_to_clu[spot] = name
    return spot_to_clu


def get_gene_cr(anno_gtf):
    gene_cr = cr.cgranges()
    gene_idx = dict()
    n_gene = 0
    with open(anno_gtf) as fp:
        for line in fp:
            if line.startswith('#'):
                continue
            ele = line.strip().rsplit('\t')
            if ele[2] != 'gene':
                continue
            chrom, start, end = ele[0], int(ele[3]), int(ele[4])
            gene_id, gene_name = '', ''
            if 'gene_id' in ele[8]:
                gene_id = ele[8][ele[8].index('gene_id') + 9:].split('\"')[0]
            if 'gene_name' in ele[8]:
                gene_name = ele[8][ele[8].index('gene_name') + 11:].split('\"')[0]
            if gene_name != '':
                gene_id = gene_name
            if gene_id != '':
                gene_cr.add(chrom, start-1, end, n_gene)
                gene_idx[n_gene] = gene_id
                n_gene += 1
    gene_cr.index()
    sys.stderr.write(f'[INFO] Loaded {n_gene} genes from GTF\n')
    return gene_cr, gene_idx


def collect_editing_count(editing_count_fn, min_read_count, min_edit_count,
                          site_wise_spot_to_editing, gene_wise_spot_to_editing,
                          site_to_gene,
                          gene_cr, gene_idx,
                          cb_to_pixel):          # ← 新增参数
    n_line = 0
    n_skip_depth = 0
    n_skip_no_cb = 0
    n_site = 0

    with open(editing_count_fn) as f:
        for line in f:
            n_line += 1
            ele = line.strip().rsplit('\t')
            if len(ele) < 7:
                continue
            chrom, pos, ref, total_count1, bases, base_quality, read_names_raw = ele
            if total_count1 == '0':
                continue

            ref, bases = ref.upper(), bases.upper()
            if ref == 'A':
                target = 'G'
            elif ref == 'T':
                target = 'C'
            else:
                continue  # 不是 A-to-I 或 C-to-U 位点，跳过

            # read_names 字段: 逗号分隔
            # 改造后格式: original_readname|CB序列
            read_names = read_names_raw.rsplit(',')

            if len(bases) != len(read_names):
                sys.stderr.write(f'[WARN] bases/read_names length mismatch at {chrom}:{pos}\n')
                continue

            ref_count = bases.count('.')
            alt_count = bases.count(target)
            total_count = ref_count + alt_count
            if total_count < min_read_count or alt_count < min_edit_count:
                n_skip_depth += 1
                continue

            genes = [gene_idx[r[2]] for r in gene_cr.overlap(chrom, int(pos)-1, int(pos))]
            site = chrom + ':' + pos
            n_site += 1

            for base, read_name in zip(bases, read_names):
                # ── 核心改动：从 read name 末尾取 CB ──────────────────
                # 改造后 read name 格式: original_name|GTACAGCGTCTACTGC
                cb = read_name.rsplit('|')[-1]
                if cb not in cb_to_pixel:
                    n_skip_no_cb += 1
                    continue
                spot = cb_to_pixel[cb]   # e.g. "139x61"
                # ────────────────────────────────────────────────────────

                site_wise_spot_to_editing[spot][site]['total_count'] += 1
                for gene in genes:
                    gene_wise_spot_to_editing[spot][gene]['total_count'] += 1
                if base == target:
                    site_wise_spot_to_editing[spot][site]['ed_count'] += 1
                    for gene in genes:
                        gene_wise_spot_to_editing[spot][gene]['ed_count'] += 1
                        site_to_gene[site] = gene

    sys.stderr.write(f'[INFO] Parsed {n_line} mpileup lines\n')
    sys.stderr.write(f'[INFO] Sites passing depth filter: {n_site}\n')
    sys.stderr.write(f'[INFO] Reads skipped (depth filter): {n_skip_depth}\n')
    sys.stderr.write(f'[INFO] Reads skipped (CB not in whitelist): {n_skip_no_cb}\n')
    return site_wise_spot_to_editing, gene_wise_spot_to_editing


def write_spot_stats(out_dir, site_wise_spot_to_editing, gene_wise_spot_to_editing, spot_to_clu):
    ed_site_count_fn = os.path.join(out_dir, 'spot_site_editing_count.txt')
    ed_site_ratio_fn = os.path.join(out_dir, 'spot_site_editing_ratio.txt')
    ed_gene_count_fn = os.path.join(out_dir, 'spot_gene_editing_count.txt')
    ed_gene_ratio_fn = os.path.join(out_dir, 'spot_gene_editing_ratio.txt')
    stat_fn = os.path.join(out_dir, 'spot_editing_stat.txt')
    all_sites = set()
    all_genes = set()
    spot_to_stat = dd(lambda: dd(lambda: 0.0))

    for spot in site_wise_spot_to_editing:
        all_sites.update(site_wise_spot_to_editing[spot].keys())
        all_genes.update(gene_wise_spot_to_editing[spot].keys())

    all_sites = sorted(all_sites)
    all_genes = sorted(all_genes)

    with open(ed_site_count_fn, 'w') as ed_site_cnt_fp, \
         open(ed_site_ratio_fn, 'w') as ed_site_ratio_fp, \
         open(ed_gene_count_fn, 'w') as ed_gene_cnt_fp, \
         open(ed_gene_ratio_fn, 'w') as ed_gene_ratio_fp:

        ed_site_cnt_fp.write('spot\t{}\n'.format('\t'.join(all_sites)))
        ed_site_ratio_fp.write('spot\t{}\n'.format('\t'.join(all_sites)))
        ed_gene_cnt_fp.write('spot\t{}\n'.format('\t'.join(all_genes)))
        ed_gene_ratio_fp.write('spot\t{}\n'.format('\t'.join(all_genes)))

        for spot in sorted(site_wise_spot_to_editing.keys()):
            ed_site_cnt_fp.write('{}'.format(spot))
            ed_site_ratio_fp.write('{}'.format(spot))
            for site in all_sites:
                edit_count = site_wise_spot_to_editing[spot][site]['ed_count']
                total_count = site_wise_spot_to_editing[spot][site]['total_count']
                edit_ratio = 0 if total_count == 0 else edit_count / total_count
                ed_site_cnt_fp.write('\t{}'.format(int(edit_count)))
                ed_site_ratio_fp.write('\t{:.4f}'.format(edit_ratio))
                spot_to_stat[spot]['n_edit_site'] += 1 if edit_count > 0 else 0
                spot_to_stat[spot]['n_edit_read'] += edit_count
                spot_to_stat[spot]['n_total_read'] += total_count
            ed_site_cnt_fp.write('\n')
            ed_site_ratio_fp.write('\n')

            ed_gene_cnt_fp.write('{}'.format(spot))
            ed_gene_ratio_fp.write('{}'.format(spot))
            for gene in all_genes:
                edit_count = gene_wise_spot_to_editing[spot][gene]['ed_count']
                total_count = gene_wise_spot_to_editing[spot][gene]['total_count']
                edit_ratio = 0 if total_count == 0 else edit_count / total_count
                ed_gene_cnt_fp.write('\t{}'.format(int(edit_count)))
                ed_gene_ratio_fp.write('\t{:.4f}'.format(edit_ratio))
                spot_to_stat[spot]['n_edit_gene'] += 1 if edit_count > 0 else 0
            ed_gene_cnt_fp.write('\n')
            ed_gene_ratio_fp.write('\n')

    with open(stat_fn, 'w') as fp:
        fp.write('spot\tn_edit_sites\tn_edit_genes\tn_edit_reads\tedit_ratio\n')
        for spot in sorted(spot_to_stat.keys()):
            n_edit_read = spot_to_stat[spot]['n_edit_read']
            n_total_read = spot_to_stat[spot]['n_total_read']
            edit_ratio = 0 if n_total_read == 0 else n_edit_read / n_total_read
            fp.write('{}\t{}\t{}\t{}\t{:.4f}\n'.format(
                spot,
                int(spot_to_stat[spot]['n_edit_site']),
                int(spot_to_stat[spot]['n_edit_gene']),
                int(n_edit_read),
                edit_ratio))


def write_clu_stats(out_dir, site_wise_spot_to_editing, gene_wise_spot_to_editing, spot_to_clu):
    ed_site_count_fn = os.path.join(out_dir, 'clu_site_editing_count.txt')
    ed_site_ratio_fn = os.path.join(out_dir, 'clu_site_editing_ratio.txt')
    ed_gene_count_fn = os.path.join(out_dir, 'clu_gene_editing_count.txt')
    ed_gene_ratio_fn = os.path.join(out_dir, 'clu_gene_editing_ratio.txt')
    stat_fn = os.path.join(out_dir, 'clu_editing_stat.txt')
    all_sites = set()
    all_genes = set()
    all_clus = set(spot_to_clu.values())
    clu_to_stat = dd(lambda: dd(lambda: 0.0))

    for spot in site_wise_spot_to_editing:
        all_sites.update(site_wise_spot_to_editing[spot].keys())
        all_genes.update(gene_wise_spot_to_editing[spot].keys())

    all_sites = sorted(all_sites)
    all_genes = sorted(all_genes)

    site_wise_clu = dd(lambda: dd(lambda: dd(lambda: 0.0)))
    gene_wise_clu = dd(lambda: dd(lambda: dd(lambda: 0.0)))

    for spot in site_wise_spot_to_editing:
        clu = spot_to_clu[spot]
        for site in site_wise_spot_to_editing[spot]:
            site_wise_clu[clu][site]['ed_count']    += site_wise_spot_to_editing[spot][site]['ed_count']
            site_wise_clu[clu][site]['total_count'] += site_wise_spot_to_editing[spot][site]['total_count']
        for gene in gene_wise_spot_to_editing[spot]:
            gene_wise_clu[clu][gene]['ed_count']    += gene_wise_spot_to_editing[spot][gene]['ed_count']
            gene_wise_clu[clu][gene]['total_count'] += gene_wise_spot_to_editing[spot][gene]['total_count']

    with open(ed_site_count_fn, 'w') as site_cnt_fp, \
         open(ed_site_ratio_fn, 'w') as site_ratio_fp, \
         open(ed_gene_count_fn, 'w') as gene_cnt_fp, \
         open(ed_gene_ratio_fn, 'w') as gene_ratio_fp:

        site_cnt_fp.write('clu\t{}\n'.format('\t'.join(all_sites)))
        site_ratio_fp.write('clu\t{}\n'.format('\t'.join(all_sites)))
        gene_cnt_fp.write('clu\t{}\n'.format('\t'.join(all_genes)))
        gene_ratio_fp.write('clu\t{}\n'.format('\t'.join(all_genes)))

        for clu in sorted(all_clus):
            site_cnt_fp.write('{}'.format(clu))
            site_ratio_fp.write('{}'.format(clu))
            for site in all_sites:
                ed = site_wise_clu[clu][site]['ed_count']
                tot = site_wise_clu[clu][site]['total_count']
                ratio = 0 if tot == 0 else ed / tot
                site_cnt_fp.write('\t{}'.format(int(ed)))
                site_ratio_fp.write('\t{:.4f}'.format(ratio))
                clu_to_stat[clu]['n_edit_site']  += 1 if ed > 0 else 0
                clu_to_stat[clu]['n_edit_read']  += ed
                clu_to_stat[clu]['n_total_read'] += tot
            site_cnt_fp.write('\n')
            site_ratio_fp.write('\n')

            gene_cnt_fp.write('{}'.format(clu))
            gene_ratio_fp.write('{}'.format(clu))
            for gene in all_genes:
                ed = gene_wise_clu[clu][gene]['ed_count']
                tot = gene_wise_clu[clu][gene]['total_count']
                ratio = 0 if tot == 0 else ed / tot
                gene_cnt_fp.write('\t{}'.format(int(ed)))
                gene_ratio_fp.write('\t{:.4f}'.format(ratio))
                clu_to_stat[clu]['n_edit_gene'] += 1 if ed > 0 else 0
            gene_cnt_fp.write('\n')
            gene_ratio_fp.write('\n')

    with open(stat_fn, 'w') as fp:
        fp.write('clu\tn_edit_sites\tn_edit_genes\tn_edit_reads\tedit_ratio\n')
        for clu in sorted(clu_to_stat.keys()):
            n_ed  = clu_to_stat[clu]['n_edit_read']
            n_tot = clu_to_stat[clu]['n_total_read']
            ratio = 0 if n_tot == 0 else n_ed / n_tot
            fp.write('{}\t{}\t{}\t{}\t{:.4f}\n'.format(
                clu,
                int(clu_to_stat[clu]['n_edit_site']),
                int(clu_to_stat[clu]['n_edit_gene']),
                int(n_ed), ratio))


def write_site_stats(out_dir, site_to_gene, site_wise_spot_to_editing):
    site_stat_fn = os.path.join(out_dir, 'site_editing_count_ratio.txt')
    site_to_count = dd(lambda: dd(lambda: 0))
    for spot in site_wise_spot_to_editing:
        for site in site_wise_spot_to_editing[spot]:
            site_to_count[site]['total_count'] += site_wise_spot_to_editing[spot][site]['total_count']
            site_to_count[site]['ed_count']    += site_wise_spot_to_editing[spot][site]['ed_count']
    with open(site_stat_fn, 'w') as out_fp:
        out_fp.write('site\ttotal_count\ted_count\tedit_ratio\n')
        for site in sorted(site_to_count.keys()):
            ed  = site_to_count[site]['ed_count']
            tot = site_to_count[site]['total_count']
            ratio = 0 if tot == 0 else ed / tot
            gene = site_to_gene.get(site, 'NA')
            out_fp.write('{}:{}\t{}\t{}\t{:.4f}\n'.format(
                gene, site, int(tot), int(ed), ratio))


def parser_argv():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('editing_site',        metavar='edit_tsv',           help='mpileup 输出 tsv')
    parser.add_argument('cb_to_pixel_tsv',     metavar='cb_to_pixel.tsv',    help='CB序列→像素坐标映射')
    parser.add_argument('pixel_to_clu_ID_tsv', metavar='pixel_to_clu_ID.tsv',help='像素坐标→cluster ID')
    parser.add_argument('clu_ID_to_region_tsv',metavar='clu_to_region.tsv',  help='cluster ID→脑区名')
    parser.add_argument('anno_gtf',            metavar='anno.gtf',           help='GTF 注释文件')
    parser.add_argument('output_dir',          metavar='out_dir',            help='输出目录')
    opt = parser.add_argument_group('optional arguments')
    opt.add_argument('-m', '--min-read-count', type=int, default=min_total_read_count,
                     help='最小总 read 深度')
    opt.add_argument('-c', '--min-edit-count', type=int, default=min_edited_read_count,
                     help='最小编辑 read 数')
    return parser.parse_args()


if __name__ == '__main__':
    args = parser_argv()
    edit_site_tsv      = args.editing_site
    cb_to_pixel_tsv    = args.cb_to_pixel_tsv
    spot_to_clu_tsv    = args.pixel_to_clu_ID_tsv
    clu_to_name_tsv    = args.clu_ID_to_region_tsv
    anno_gtf           = args.anno_gtf
    out_dir            = args.output_dir
    min_read_count     = args.min_read_count
    min_edit_count     = args.min_edit_count

    os.makedirs(out_dir, exist_ok=True)

    cb_to_pixel  = get_cb_to_pixel(cb_to_pixel_tsv)
    spot_to_clu  = get_spot_to_clu(spot_to_clu_tsv, clu_to_name_tsv)
    gene_cr, gene_idx = get_gene_cr(anno_gtf)

    site_wise_spot_to_editing = dd(lambda: dd(lambda: dd(lambda: 0.0)))
    gene_wise_spot_to_editing = dd(lambda: dd(lambda: dd(lambda: 0.0)))
    site_to_gene = dict()

    collect_editing_count(
        edit_site_tsv, min_read_count, min_edit_count,
        site_wise_spot_to_editing, gene_wise_spot_to_editing,
        site_to_gene, gene_cr, gene_idx,
        cb_to_pixel)

    write_spot_stats(out_dir, site_wise_spot_to_editing, gene_wise_spot_to_editing, spot_to_clu)
    write_clu_stats( out_dir, site_wise_spot_to_editing, gene_wise_spot_to_editing, spot_to_clu)
    write_site_stats(out_dir, site_to_gene, site_wise_spot_to_editing)
    sys.stderr.write('[INFO] Done.\n')
