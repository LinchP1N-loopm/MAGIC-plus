#!/usr/bin/env python3
# Trackplot 绘制 Regions pseudobulk BAM coverage + junction 弧线。
# 用法：python 02_plot_sashimi_junctions.py Apoe --region CA1 --samples G21-4,G21-1,G21-2
# 特殊模式：python 02_plot_sashimi_junctions.py Apoe --tracks G21-4:CA1,G21-1:CA3


import argparse
import csv
import os
import re
import subprocess
import sys
from functools import lru_cache
PROJECT_DIR = '/histor/zhao/projects/lupeng/zhaolab/ZJJ/Fig4_sashimi'
GTF = '/histor/zhao/pangkun/database/ref_mm10_M23_release98/mm10.gtf'
COORDS_DIR = f'{PROJECT_DIR}/gene_coords'
GTF_DIR = f'{PROJECT_DIR}/tracks_gtf'
BAM_BASE = '/histor/zhao/pangkun/Zhaolab/Lupeng/Zhujunjie/20260705_Fig3/rmats_works/pseudobulk_bams'
MERGED_BAM_DIR = f'{PROJECT_DIR}/bigwig'
DEFAULT_OUTDIR = f'{PROJECT_DIR}/sashimi_junctions'
SAMPLE_INFO = {
    'G21-4': {'group': 'Young', 'color': "#7faba2"},
    'G21-7': {'group': 'Young', 'color': '#7faba2'},
    'G21-1': {'group': 'Aged', 'color': '#ac5136'},
    'G21-8': {'group': 'Aged', 'color': '#ac5136'},
    'G21-2': {'group': 'AD', 'color': '#B39DDB'},
    'G21-5': {'group': 'AD', 'color': '#B39DDB'},
}
DEFAULT_SAMPLE_ORDER = ['G21-4', 'G21-7', 'G21-1', 'G21-8', 'G21-2', 'G21-5']
VALID_REGIONS = ['CA1', 'CA3', 'CA1CA3', 'CTX_1', 'CTX_2', 'CTX_3', 'DG', 'FT', 'HIP', 'HY', 'STR', 'TH', 'VS']
WIDTH = 8
HEIGHT = 1.5
DPI = 300
FONT_SIZE = 8
FONT = 'Arial'
FMT = 'pdf'
INTRON_SCALE = 0.5
EXON_SCALE = 1.0
ANNOTATION_SCALE = 0.35
DISTANCE_RATIO = 0.02
REF_COLOR = 'black'
N_Y_TICKS = 3
FILL_STEP = 'post'
JUNCTION_THRESHOLD = 10

def bam_path(sample, region):
    base = MERGED_BAM_DIR if region == 'CA1CA3' else BAM_BASE
    return os.path.join(base, sample, f'{sample}_{region}_pseudobulk.bam')

def parse_args():
    p = argparse.ArgumentParser(description='使用 Regions pseudobulk BAM 绘制 Trackplot junction 图')
    p.add_argument('genes', nargs='*', help='基因名，可输入多个，例如 Apoe Gfap')
    p.add_argument('--gene_info', default=None, help='基因列表 TSV，需含 Gene_Name、gene_name 或 gene 列')
    p.add_argument('--coords', default=None, help='自定义坐标 TSV，需含 gene、chr、plot_start、plot_end')
    p.add_argument('--region', default=None, choices=VALID_REGIONS, help='统一脑区，例如 CA1、CA3、DG、CTX_1')
    p.add_argument('--samples', default=None, help='选择样本，逗号分隔；输入顺序即绘图顺序。例如 G21-4,G21-1,G21-2。若不指定，则使用默认6样本顺序。')
    p.add_argument('--tracks', default=None, help='高级模式：直接指定 sample:region，逗号分隔。例如 G21-4:CA1,G21-1:CA3,G21-2:DG')
    p.add_argument('--outdir', default=DEFAULT_OUTDIR, help=f'输出目录，默认 {DEFAULT_OUTDIR}')
    p.add_argument('--fmt', default=FMT, choices=['pdf', 'svg', 'png'], help=f'输出格式，默认 {FMT}')
    p.add_argument('--plot-region', default=None, help='直接指定区域，例如 chr2:157559078-157563522:+')
    p.add_argument('--plot-start', type=int, default=None, help='覆盖默认坐标起点')
    p.add_argument('--plot-end', type=int, default=None, help='覆盖默认坐标终点')
    p.add_argument('--flank', type=int, default=0, help='最终坐标两侧扩展 bp，默认 0')
    p.add_argument('--width', type=int, default=WIDTH, help=f'整张图宽度；trackplot要求整数，默认 {WIDTH}')
    p.add_argument('--height', type=float, default=HEIGHT, help=f'每个 subplot 高度，默认 {HEIGHT}')
    p.add_argument('--dpi', type=int, default=DPI, help=f'输出分辨率，默认 {DPI}')
    p.add_argument('--font-size', type=int, default=FONT_SIZE, help=f'字体大小，默认 {FONT_SIZE}')
    p.add_argument('--font', default=FONT, help=f'字体，默认 {FONT}')
    p.add_argument('--intron-scale', type=float, default=INTRON_SCALE, help='内含子横向缩放；默认0.5。长内含子可尝试 0.1、0.05、0.02')
    p.add_argument('--exon-scale', type=float, default=EXON_SCALE, help='外显子横向缩放，默认1；可用1.5或2')
    p.add_argument('--annotation-scale', type=float, default=ANNOTATION_SCALE, help='底部 transcript annotation 高度比例')
    p.add_argument('--distance-ratio', type=float, default=DISTANCE_RATIO, help='transcript 标签与结构线之间距离')
    p.add_argument('--ref-color', default=REF_COLOR, help='annotation结构颜色')
    p.add_argument('--transcripts-to-show', default=None, help='只显示指定转录本，多个名称用逗号分隔')
    p.add_argument('--junction-threshold', type=int, default=JUNCTION_THRESHOLD, help='junction 最少支持 reads 数，默认 10；设为 0 显示全部')
    p.add_argument('--library-type', default='frf', help='BAM 建库链方向；默认按 Trackplot 文档示例 frf，请按实际文库确认')
    p.add_argument('--same-y', action='store_true', help='所有 coverage 共用相同 y 轴范围')
    p.add_argument('--n-y-ticks', type=int, default=N_Y_TICKS, help=f'y轴刻度数量，默认 {N_Y_TICKS}')
    p.add_argument('--fill-step', default=FILL_STEP, choices=['pre', 'post', 'mid'], help=f'coverage填充阶梯方式，默认 {FILL_STEP}')
    p.add_argument('--title', default=None, help='自定义标题；默认使用基因名')
    p.add_argument('--no-title', action='store_true', help='不显示标题')
    p.add_argument('--no-region-text', action='store_true', help='不显示底部染色体坐标文字；仅在当前trackplot支持时生效')
    p.add_argument('--process', type=int, default=1, help='trackplot进程数，默认1')
    return p.parse_args()

# 选择样本及脑区
def select_tracks(args):
    if args.tracks and (args.region or args.samples):
        raise ValueError('--tracks 不能与 --region/--samples 同时使用')
    tracks = []
    if args.tracks:
        for item in args.tracks.split(','):
            item = item.strip()
            if not item:
                continue
            if ':' not in item:
                raise ValueError(f'track格式错误: {item}；应为 sample:region')
            sample, region = [x.strip() for x in item.split(':', 1)]
            if sample not in SAMPLE_INFO:
                raise ValueError(f'未知样本: {sample}')
            if region not in VALID_REGIONS:
                raise ValueError(f'未知脑区: {region}')
            tracks.append({'sample': sample, 'region': region, 'label': f'{sample} {region}', 'color': SAMPLE_INFO[sample]['color']})
        if not tracks:
            raise ValueError('--tracks 没有读取到有效 track')
        return tracks
    if args.region is None:
        raise ValueError('请指定 --region，例如 --region CA1；或者使用 --tracks sample:region')
    if args.samples:
        samples = [x.strip() for x in args.samples.split(',') if x.strip()]
    else:
        samples = DEFAULT_SAMPLE_ORDER.copy()
    if not samples:
        raise ValueError('没有读取到任何样本')
    for sample in samples:
        if sample not in SAMPLE_INFO:
            raise ValueError(f'未知样本: {sample}')
        tracks.append({'sample': sample, 'region': args.region, 'label': sample, 'color': SAMPLE_INFO[sample]['color']})
    return tracks

# 读取基因列表
def load_gene_info(fp):
    out = []
    with open(fp) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            gene = (row.get('Gene_Name') or row.get('gene_name') or row.get('gene') or '').strip()
            if gene:
                out.append(gene)
    return out

# 读取基因坐标
def load_coords_file(fp):
    out = []
    with open(fp) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            item = {'name': row['gene'], 'chr': row['chr'], 'plot_start': int(row['plot_start']), 'plot_end': int(row['plot_end'])}
            strand = (row.get('strand') or '').strip()
            if strand in {'+', '-'}:
                item['strand'] = strand
            out.append(item)
    return out
REGION_RE = re.compile('^(?P<chrom>[^:]+):(?P<start>\\d+)-(?P<end>\\d+)(?::(?P<strand>[+-]))?$')

def parse_plot_region(value):
    match = REGION_RE.match(value.strip())
    if not match:
        raise ValueError('plot-region格式错误，应为 chr:start-end 或 chr:start-end:+/-')
    start = int(match.group('start'))
    end = int(match.group('end'))
    if start >= end:
        raise ValueError('plot-region中 start 必须小于 end')
    return {'chr': match.group('chrom'), 'plot_start': start, 'plot_end': end, 'strand': match.group('strand')}

# 确定绘图范围
def make_region(gene, args):
    chrom = gene['chr']
    start = gene['plot_start']
    end = gene['plot_end']
    strand = gene.get('strand')
    if args.plot_region:
        custom = parse_plot_region(args.plot_region)
        chrom = custom['chr']
        start = custom['plot_start']
        end = custom['plot_end']
        strand = custom.get('strand')
    if args.plot_start is not None:
        start = args.plot_start
    if args.plot_end is not None:
        end = args.plot_end
    start = max(1, start - args.flank)
    end = end + args.flank
    if start >= end:
        raise ValueError(f'最终绘图坐标不合法：{chrom}:{start}-{end}')
    region = f'{chrom}:{start}-{end}'
    if strand in {'+', '-'}:
        region += f':{strand}'
    return region

# 基因注释
def gene_gtf(gene_name):
    """
    优先读取：

      PROJECT_DIR/tracks_gtf/<gene>.gtf

    不存在时，从 mm10.gtf 自动提取，
    然后 bgzip + tabix。
    """
    os.makedirs(GTF_DIR, exist_ok=True)
    raw = f'{GTF_DIR}/{gene_name}.gtf'
    sorted_gtf = f'{GTF_DIR}/{gene_name}.sorted.gtf'
    gz = f'{GTF_DIR}/{gene_name}.gtf.gz'
    tbi = f'{gz}.tbi'
    if not os.path.exists(raw) or os.path.getsize(raw) == 0:
        print(f'    extracting GTF for {gene_name} ...')
        cmd = f'''grep -w 'gene_name "{gene_name}"' "{GTF}" | grep -v 'transcript_type "retained_intron"' > "{raw}"'''
        subprocess.run(cmd, shell=True, check=True)
        if not os.path.exists(raw) or os.path.getsize(raw) == 0:
            raise RuntimeError(f'未从参考GTF中找到基因：{gene_name}')
    need_repack = not os.path.exists(gz) or not os.path.exists(tbi) or os.path.getmtime(raw) > os.path.getmtime(gz)
    if need_repack:
        print(f'    bgzip + tabix for {gene_name} ...')
        subprocess.run(f'(grep "^#" "{raw}" || true; grep -v "^#" "{raw}" | sort -k1,1 -k4,4n) > "{sorted_gtf}"', shell=True, check=True)
        subprocess.run(f'bgzip -f -c "{sorted_gtf}" > "{gz}"', shell=True, check=True)
        subprocess.run(f'tabix -f -p gff "{gz}"', shell=True, check=True)
    else:
        print(f'    using cached {gene_name}.gtf.gz')
    return gz

# 生成 BAM 配置
def write_density_tsv(tsv_path, tracks, library_type):
    """
    Trackplot BAM density TSV:

      路径    bam    标签    颜色    建库类型
    """
    with open(tsv_path, 'w') as fh:
        for track in tracks:
            path = bam_path(track['sample'], track['region'])
            fh.write(f"{path}\tbam\t{track['label']}\t{track['color']}\t{library_type}\n")

@lru_cache(maxsize=1)
def trackplot_help():
    result = subprocess.run(['trackplot', '--help'], capture_output=True, text=True)
    return (result.stdout or '') + '\n' + (result.stderr or '')

def trackplot_supports(option):
    return option in trackplot_help()

def add_supported_flag(cmd, option, value=None):
    if not trackplot_supports(option):
        print(f'    [WARN] 当前 trackplot 不支持 {option}，已自动跳过')
        return
    cmd.append(option)
    if value is not None:
        cmd.append(str(value))

def track_tag(tracks):
    regions = {x['region'] for x in tracks}
    if len(regions) == 1:
        region = tracks[0]['region']
        samples = '_'.join((x['sample'] for x in tracks))
        return f'{region}_{samples}'
    return '_'.join((f"{x['sample']}-{x['region']}" for x in tracks))

# 调用 Trackplot 绘图
def plot_gene(gene, outdir, fmt, args, tracks):
    name = gene['name']
    try:
        region = make_region(gene, args)
    except ValueError as exc:
        print(f'  [ERROR] {name}: {exc}')
        return False
    os.makedirs(outdir, exist_ok=True)
    tag = track_tag(tracks)
    density_tsv = os.path.join(outdir, f'{name}_{tag}_density.tsv')
    write_density_tsv(density_tsv, tracks, args.library_type)
    region_tag = re.sub('[^A-Za-z0-9_.-]+', '-', region)
    out_file = os.path.join(outdir, f'{name}_{tag}_junction_intron{args.intron_scale}_{region_tag}.{fmt}')
    try:
        gtf_file = gene_gtf(name)
    except (subprocess.CalledProcessError, RuntimeError) as exc:
        print(f'  [ERROR] GTF prep failed for {name}: {exc}')
        return False
    cmd = ['trackplot', '-e', region, '-r', gtf_file, '--density', density_tsv, '-o', out_file,
           '--show-junction-num', '-t', str(args.junction_threshold), '--dpi', str(args.dpi), '--width', str(int(args.width)),
           '--height', str(args.height), '--font-size', str(args.font_size), '--font', args.font,
           '--intron-scale', str(args.intron_scale), '--exon-scale', str(args.exon_scale),
           '--annotation-scale', str(args.annotation_scale), '--distance-ratio', str(args.distance_ratio),
           '--ref-color', args.ref_color, '--fill-step', args.fill_step, '-p', str(args.process)]
    if args.same_y:
        cmd.append('--same-y')
    add_supported_flag(cmd, '--n-y-ticks', args.n_y_ticks)
    if args.transcripts_to_show:
        add_supported_flag(cmd, '--transcripts-to-show', args.transcripts_to_show)
    if args.no_region_text:
        add_supported_flag(cmd, '--no-region-text')
    if not args.no_title:
        title = args.title if args.title is not None else name
        cmd.extend(['--title', title])
    print(f'\n  plotting {name} ({region})')
    print(f'    tracks          : {len(tracks)}')
    for x in tracks:
        print(f"      {x['sample']} | {x['region']}")
    print(f'    width × height  : {args.width} × {args.height}')
    print(f'    intron-scale    : {args.intron_scale}')
    print(f'    exon-scale      : {args.exon_scale}')
    print(f'    annotation      : {args.annotation_scale}')
    print(f'    same-y          : {args.same_y}')
    print(f'    density TSV     : {density_tsv}')
    print(f'    output          : {out_file}')
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f'  [ERROR] {name}')
        print(result.stderr[-3000:])
        print('\n  [CMD]')
        print('  ' + ' '.join(cmd))
        return False
    print(f'  ✓ {out_file}')
    return True

# 选择基因
def load_genes(args):
    if args.coords:
        return load_coords_file(args.coords)
    if args.gene_info:
        gene_names = load_gene_info(args.gene_info)
    elif args.genes:
        gene_names = args.genes
    else:
        if not os.path.isdir(COORDS_DIR):
            return []
        gene_names = sorted((filename[:-4] for filename in os.listdir(COORDS_DIR) if filename.endswith('.tsv')))
    genes = []
    for gene_name in gene_names:
        coord_file = os.path.join(COORDS_DIR, f'{gene_name}.tsv')
        if not os.path.exists(coord_file):
            print(f'  [SKIP] {gene_name}: coords文件不存在：{coord_file}')
            continue
        loaded = load_coords_file(coord_file)
        if loaded:
            genes.append(loaded[0])
    if args.gene_info:
        order = load_gene_info(args.gene_info)
        rank = {gene: i for i, gene in enumerate(order)}
        wanted = set(order)
        genes = sorted([x for x in genes if x['name'] in wanted], key=lambda x: rank.get(x['name'], 999999))
    elif args.genes:
        rank = {gene: i for i, gene in enumerate(args.genes)}
        wanted = set(args.genes)
        genes = sorted([x for x in genes if x['name'] in wanted], key=lambda x: rank.get(x['name'], 999999))
    return genes

# 主程序
def main():
    args = parse_args()
    try:
        tracks = select_tracks(args)
    except ValueError as exc:
        print(f'[ERROR] {exc}')
        sys.exit(1)
    print('Tracks to plot:')
    for i, x in enumerate(tracks, 1):
        print(f"  {i:>2}. {x['sample']} | {x['region']} | {SAMPLE_INFO[x['sample']]['group']}")
    print('\nChecking BAM files ...')
    missing = []
    for x in tracks:
        path = bam_path(x['sample'], x['region'])
        exists = os.path.exists(path)
        print(f"  {('✓' if exists else '✗')} {x['sample']} | {x['region']}  ({path})")
        if not exists:
            missing.append(path)
    if missing:
        print('\n[ERROR] Missing BAM，请先生成以下文件：')
        for path in missing:
            print(f'  {path}')
        sys.exit(1)
    genes = load_genes(args)
    if not genes:
        print('[ERROR] 没有找到可绘制的基因')
        print(f'请先运行 00_get_coords_regions.sh，坐标目录为：{COORDS_DIR}')
        sys.exit(1)
    if args.plot_region and len(genes) > 1:
        print('[ERROR] --plot-region 只适合一次绘制一个基因；多个基因请分别运行')
        sys.exit(1)
    print(f"\nGenes to plot ({len(genes)}): {[x['name'] for x in genes]}")
    success = []
    failed = []
    for gene in genes:
        ok = plot_gene(gene=gene, outdir=args.outdir, fmt=args.fmt, args=args, tracks=tracks)
        if ok:
            success.append(gene['name'])
        else:
            failed.append(gene['name'])
    print(f'\ndone — success: {len(success)}  failed: {len(failed)}')
    if failed:
        print(f'failed: {failed}')
        sys.exit(1)
if __name__ == '__main__':
    main()
