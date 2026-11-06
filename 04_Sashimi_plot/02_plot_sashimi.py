#!/usr/bin/env python3
"""
02_plot_sashimi_v8_1.py

在 v7 基础上增加：
1. 内含子压缩参数：--intron-scale / --exon-scale / --fixed-intron-bp
2. 论文友好输出：默认 pdf，保留矢量；可选 svg/png
3. 注释轨道优化：--annotation-scale / --distance-ratio / --ref-color
4. 版面优化：--width / --height / --font-size / --same-y / --hide-y-label
5. 转录本筛选：--transcripts-to-show / --choose-primary / --remove-empty
6. 默认绘制 junction 弧线及支持数，形成真正的 Sashimi plot
7. 可直接覆盖绘图区间：--plot-region / --plot-start / --plot-end / --flank
8. junction 布局增强：--junctions-on-top / --included-junctions / --sites

Usage:
    # 推荐：普通长基因，内含子按比例压缩
    python 02_plot_sashimi_v8.py Myl6 -r Cortex_2,VS,DG --fmt pdf --intron-scale 0.02 --exon-scale 1

    # 如果 trackplot >= 0.5.3，可把每个 intron 固定显示成 250 bp
    python 02_plot_sashimi_v8.py Myl6 -r Cortex_2,VS,DG --fixed-intron-bp 250 --fmt pdf

    # 只显示指定 transcript
    python 02_plot_sashimi_v8.py Myl6 -r Cortex_2,VS,DG \
        --transcripts-to-show Myl6-206,Myl6-201 --intron-scale 0.02 --fmt pdf

    # V8 默认显示 junction 弧线和数字；可直接指定绘图范围
    python 02_plot_sashimi_v8.py Myl6 -r TH,DG,VS \
        --plot-region chr10:128490499-128494000:- \
        --junction-threshold 5 --fmt pdf
"""

import os
import sys
import csv
import subprocess
import argparse

# ── paths ──────────────────────────────────────────────────────────────────────
PROJECT_DIR = "/path/to/project/20260427_AS_analysis/Fig2_trackplot"
GTF         = "/path/to/project/database/ref_mm10_M23_release98/mm10.gtf"
COORDS_DIR  = f"{PROJECT_DIR}/gene_coords"
GTF_DIR     = f"{PROJECT_DIR}/tracks_gtf"

BAM_DIR    = ("/path/to/project/20260427_AS_analysis"
              "/Fig2_rmats/ASclu_results/pseudobulk_bams/G21-4-0504")
BAM_PREFIX = "G21-4-0504_"
BAM_SUFFIX = "_pseudobulk.bam"

# ── 短名 → 标准名（BAM 路径用标准名，palette 用标准名）─────────────────────────
ALIAS = {
    "EPI":  "Epithalamus",
    "FT":   "Fiber_tracts",
    "HY":   "Hypothalamus",
    "STR":  "Striatum",
    "TH":   "Thalamus",

    "CA1":      "CA1",
    "CA3":      "CA3",
    "Cortex_1": "Cortex_1",
    "Cortex_2": "Cortex_2",
    "Cortex_3": "Cortex_3",
    "DG":       "DG",
    "VS":       "VS",
}

BAM_NAME = {
    "CA1":          "CA1",
    "CA3":          "CA3",
    "Cortex_1":     "Cortex_1",
    "Cortex_2":     "Cortex_2",
    "Cortex_3":     "Cortex_3",
    "DG":           "DG",
    "Epithalamus":  "Epithalamus",
    "Fiber_tracts": "Fiber_tracts",
    "Hypothalamus": "Hypothalamus",
    "Striatum":     "Striatum",
    "Thalamus":     "Thalamus",
    "VS":           "VS",
}

CLUSTER_PALETTE = {
    "CA1":          "#E64B35",
    "CA3":          "#4DBBD5",
    "Cortex_1":     "#00A087",
    "Cortex_2":     "#3C5488",
    "Cortex_3":     "#F39B7F",
    "DG":           "#8491B4",
    "Epithalamus":  "#91D1C2",
    "Fiber_tracts": "#725AC1",
    "Hypothalamus": "#7E6148",
    "Striatum":     "#B09C85",
    "Thalamus":     "#CCA26E",
    "VS":           "#666666",
    "Meninges":     "#E0E0E0",
    "unknown":      "#BDBDBD",
}

_DEFAULT_REGIONS = ["Cortex_2", "VS", "DG"]

# ── figure defaults：论文友好配置 ──────────────────────────────────────────────
WIDTH              = 12
HEIGHT             = 1.8
DPI                = 300
FONTSIZE           = 8
FMT                = "pdf"

NORMALIZE_FORMAT   = "rpkm"

# 真正的 Sashimi plot：默认显示 junction 弧线和支持 reads 数
# 5 是较平衡的起始阈值；低表达事件可改为 1~3，高覆盖数据可改为 10
JUNCTION_THRESHOLD = 5
SHOW_JUNCTION_NUM  = True
JUNCTIONS_ON_TOP   = True
SHOW_MEAN_JUNCTION_NUM = False
FILL_STEP          = "post"
N_Y_TICKS          = 3
LIBRARY_TYPE       = "fru"   # fru=unstranded；frf/frs=链特异文库

# 长基因推荐先用比例压缩，兼容所有版本 trackplot
INTRON_SCALE       = "0.5"
EXON_SCALE         = "1"

# 增大 annotation track，有助于 transcript 箭头/结构在论文 PDF 里更清楚
ANNOTATION_SCALE   = "0.35"
DISTANCE_RATIO     = "0.02"
REF_COLOR          = "black"


# ── args ───────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(
        description="trackplot sashimi v8：默认显示 junction 弧线和支持数，可覆盖绘图范围",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    p.add_argument("genes", nargs="*", help="基因名，可多个")
    p.add_argument("--gene_info", default=None, help="基因列表 txt（含 Gene_Name/gene_name 列）")
    p.add_argument("--outdir", default=f"{PROJECT_DIR}/sashimi_output_v8", help="输出目录")
    p.add_argument("--coords", default=None, help="指定 gene_coords tsv 路径")
    p.add_argument(
        "--plot-region", default="",
        help=("直接覆盖绘图范围，格式 chr:start-end 或 chr:start-end:+/-。"
              "例如 chr10:128490499-128494000:-；仅建议单基因运行时使用。"),
    )
    p.add_argument("--plot-start", default=None, type=int,
                   help="覆盖 coords 中的绘图起点；必须与 --plot-end 同时使用。")
    p.add_argument("--plot-end", default=None, type=int,
                   help="覆盖 coords 中的绘图终点；必须与 --plot-start 同时使用。")
    p.add_argument("--flank", default=0, type=int,
                   help="在 coords 原范围两侧各扩展 N bp，例如 --flank 500。")

    p.add_argument(
        "-r", "--regions",
        default=None,
        help=(
            "要绘制的脑区，逗号分隔，支持短名。\n"
            "例: -r Cortex_1,CA1,EPI\n"
            "可用: CA1 CA3 Cortex_1 Cortex_2 Cortex_3 DG VS\n"
            "      EPI(Epithalamus) FT(Fiber_tracts) HY(Hypothalamus)\n"
            "      STR(Striatum) TH(Thalamus)\n"
            "不传则使用默认脑区: " + ",".join(_DEFAULT_REGIONS)
        ),
    )

    # output / layout
    p.add_argument("--fmt", default=FMT, choices=["pdf", "svg", "png"], help=f"输出格式，默认 {FMT}")
    p.add_argument("--width", default=WIDTH, type=int, help=f"输出宽度（trackplot 仅接受整数），默认 {WIDTH}")
    p.add_argument("--height", default=HEIGHT, type=float, help=f"每个 track 高度，默认 {HEIGHT}")
    p.add_argument("--dpi", default=DPI, type=int, help=f"分辨率，默认 {DPI}")
    p.add_argument("--font-size", default=FONTSIZE, type=int, help=f"字体大小，默认 {FONTSIZE}")
    p.add_argument("--font", default="Arial", help="字体，默认 Arial")
    p.add_argument("--hide-y-label", action="store_true", help="隐藏 y 轴标签")
    p.add_argument("--same-y", action="store_true", help="不同 density track 使用同一 y 轴范围")
    p.add_argument("--raster", action="store_true", help="对部分复杂图层 rasterize；PDF/SVG 文件会更小，但不再完全矢量")
    p.add_argument("--n-y-ticks", default=N_Y_TICKS, type=int,
                   help=f"每个 density 轨道的 y 轴刻度数量，默认 {N_Y_TICKS}")
    p.add_argument("--no-title", action="store_true", help="不显示主标题（新版 trackplot 支持）")
    p.add_argument("--no-region-text", action="store_true",
                   help="不显示最上方 event 坐标文字（新版 trackplot 支持）")

    # reference annotation
    p.add_argument("--intron-scale", default=INTRON_SCALE,
                   help=("内含子缩放比例。长基因推荐 0.01~0.05；"
                         "例如 0.02 表示 intron 显示长度压缩到 2%%。"))
    p.add_argument("--exon-scale", default=EXON_SCALE,
                   help="外显子缩放比例，默认 1；如果外显子太细可试 2~5。")
    p.add_argument("--fixed-intron-bp", default=None, type=int,
                   help=("如果 trackplot >=0.5.3，可把每个 intron 固定显示为 N bp；"
                         "例如 --fixed-intron-bp 250。设置后会覆盖 --intron-scale。"))
    p.add_argument("--annotation-scale", default=ANNOTATION_SCALE,
                   help=f"annotation 轨道占比，默认 {ANNOTATION_SCALE}；箭头/转录本不清楚可调到 0.4~0.6")
    p.add_argument("--distance-ratio", default=DISTANCE_RATIO,
                   help=f"transcript label 与线之间距离，默认 {DISTANCE_RATIO}")
    p.add_argument("--ref-color", default=REF_COLOR,
                   help=f"exon/annotation 颜色，默认 {REF_COLOR}")
    p.add_argument("--remove-empty", action="store_true",
                   help="移除该区域内没有内容或无关的 transcript，减少空白")
    p.add_argument("--choose-primary", action="store_true",
                   help="只画 primary transcript")
    p.add_argument("--transcripts-to-show", default="",
                   help="只显示指定 transcript name/id，逗号分隔，如 Myl6-206,Myl6-201")
    p.add_argument("--show-id", action="store_true", help="显示 gene/transcript id")
    p.add_argument("--show-exon-id", action="store_true", help="显示 exon id")
    p.add_argument("--no-gene", action="store_true", help="不在 transcript id 旁边显示 gene id/name")

    # density / junction
    p.add_argument("--normalize-format", default=NORMALIZE_FORMAT,
                   choices=["count", "cpm", "rpkm"], help=f"coverage 归一化，默认 {NORMALIZE_FORMAT}")
    p.add_argument("--junction-threshold", default=JUNCTION_THRESHOLD, type=int,
                   help=f"junction 最低支持数，默认 {JUNCTION_THRESHOLD}；低表达事件可试 1~3")
    p.set_defaults(show_junction_num=SHOW_JUNCTION_NUM)
    p.add_argument("--show-junction-num", dest="show_junction_num", action="store_true",
                   help="显示 junction 支持 reads 数（V8 默认开启）")
    p.add_argument("--hide-junction-num", dest="show_junction_num", action="store_false",
                   help="只显示 junction 弧线，不显示数字")
    p.add_argument("--show-mean-junction-num", action="store_true",
                   default=SHOW_MEAN_JUNCTION_NUM,
                   help="聚合多个 BAM 到同一 label 时显示平均 junction count；新版 trackplot 支持")
    p.set_defaults(junctions_on_top=JUNCTIONS_ON_TOP)
    p.add_argument("--junctions-on-top", dest="junctions_on_top", action="store_true",
                   help="把 junction 弧线放到 coverage 上方（V8 默认开启；新版 trackplot 支持）")
    p.add_argument("--junctions-inline", dest="junctions_on_top", action="store_false",
                   help="不强制把 junction 弧线置顶，使用 trackplot 默认布局")
    p.add_argument("--included-junctions", default="",
                   help="指定需要保留的 junction，逗号分隔，如 chr1:1-100,chr1:200-300")
    p.add_argument("--fill-step", default=FILL_STEP, choices=["pre", "post", "mid"],
                   help=f"coverage 填充步进方式，默认 {FILL_STEP}")
    p.add_argument("--sites", default="",
                   help="额外绘制垂直指示线，逗号分隔坐标，如 128491723,128492717")
    p.add_argument("--density-by-strand", action="store_true",
                   help="按链方向画 density；启用时 --junctions-on-top 不生效")
    p.add_argument("--library", default=LIBRARY_TYPE, choices=["fru", "frf", "frs"],
                   help=f"BAM 文库方向：fru=unstranded，frf/frs=链特异；默认 {LIBRARY_TYPE}")
    p.add_argument("--log", default=None, choices=["0", "2", "10", "zscore"],
                   help="y 轴变换：0/2/10/zscore")

    # event strand
    p.add_argument("--strand", default="auto", choices=["auto", "+", "-", "none"],
                   help=("event 区域链方向。auto 会读取 coords 里的 strand 列；"
                         "coords 没有 strand 列则不加链方向。也可强制 + 或 -。"))

    # performance
    p.add_argument("-p", "--process", default=1, type=int, help="trackplot 使用 CPU 数，默认 1")

    return p.parse_args()


def resolve_regions(regions_str):
    raw = [r.strip() for r in regions_str.split(",") if r.strip()]
    resolved, unknown = [], []
    for r in raw:
        if r in ALIAS:
            resolved.append(ALIAS[r])
        else:
            unknown.append(r)
    if unknown:
        print(f"[ERROR] 未知脑区: {unknown}")
        print(f"  可用: {sorted(ALIAS.keys())}")
        sys.exit(1)
    return resolved


# ── helpers ────────────────────────────────────────────────────────────────────
def load_gene_info(fp):
    out = []
    with open(fp) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            g = (row.get("Gene_Name") or row.get("gene_name") or "").strip()
            if g:
                out.append(g)
    return out


def load_coords_file(fp):
    out = []
    with open(fp) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            item = {
                "name": row["gene"],
                "chr": row["chr"],
                "plot_start": int(row["plot_start"]),
                "plot_end": int(row["plot_end"]),
            }
            if "strand" in row and row["strand"] in {"+", "-"}:
                item["strand"] = row["strand"]
            out.append(item)
    return out


def gene_gtf(gene_name):
    """
    返回 bgzip+tabix 的 .gtf.gz 路径。

    优先级：
      1. 优先读 <gene>.gtf，方便手动删减转录本
      2. 若 .gtf 比 .gtf.gz 新，自动重新压缩
      3. 若 .gtf 不存在，从全局 GTF 提取 gene_name
    """
    os.makedirs(GTF_DIR, exist_ok=True)
    raw        = f"{GTF_DIR}/{gene_name}.gtf"
    sorted_gtf = f"{GTF_DIR}/{gene_name}.sorted.gtf"
    gz         = f"{GTF_DIR}/{gene_name}.gtf.gz"
    tbi        = f"{gz}.tbi"

    if not os.path.exists(raw) or os.path.getsize(raw) == 0:
        print(f"    extracting GTF for {gene_name} ...")
        cmd_extract = (
            f'grep -w \'gene_name "{gene_name}"\' "{GTF}" '
            f'| grep -v \'transcript_type "retained_intron"\' > "{raw}"'
        )
        subprocess.run(cmd_extract, shell=True, check=True)

    need_repack = (
        not os.path.exists(gz)
        or not os.path.exists(tbi)
        or os.path.getmtime(raw) > os.path.getmtime(gz)
    )

    if need_repack:
        print(f"    bgzip + tabix for {gene_name} ...")
        cmd_sort = (
            f'(grep "^#" "{raw}" || true; '
            f'grep -v "^#" "{raw}" | sort -k1,1 -k4,4n) > "{sorted_gtf}"'
        )
        subprocess.run(cmd_sort, shell=True, check=True)
        subprocess.run(f'bgzip -f "{sorted_gtf}" -o "{gz}"', shell=True, check=True)
        subprocess.run(f'tabix -p gff "{gz}"', shell=True, check=True)
    else:
        print(f"    using cached {gene_name}.gtf.gz")

    return gz


def bam_path(std_name):
    return os.path.join(BAM_DIR, f"{BAM_PREFIX}{BAM_NAME[std_name]}{BAM_SUFFIX}")


def write_density_tsv(tsv_path, samples, library="fru"):
    with open(tsv_path, "w") as f:
        for s in samples:
            f.write(f"{bam_path(s['name'])}\tbam\t{s['name']}\t{s['color']}\t{library}\n")



def parse_plot_region(region_text):
    """解析 chr:start-end[:strand]。"""
    import re

    pattern = re.compile(r"^([^:]+):(\d+)-(\d+)(?::([+-]))?$")
    match = pattern.match(region_text.strip())
    if not match:
        raise ValueError(
            "--plot-region 格式错误，应为 chr:start-end 或 chr:start-end:+/-，"
            "例如 chr10:128490499-128494000:-"
        )

    chrom, start, end, strand = match.groups()
    start, end = int(start), int(end)
    if start >= end:
        raise ValueError(f"绘图范围必须满足 start < end，当前为 {start}-{end}")

    result = {"chr": chrom, "plot_start": start, "plot_end": end}
    if strand in {"+", "-"}:
        result["strand"] = strand
    return result


def apply_range_overrides(genes, args):
    """应用 --plot-region、--plot-start/--plot-end 或 --flank。"""
    if args.plot_region and (args.plot_start is not None or args.plot_end is not None):
        raise ValueError("--plot-region 不能与 --plot-start/--plot-end 同时使用")

    if (args.plot_start is None) != (args.plot_end is None):
        raise ValueError("--plot-start 和 --plot-end 必须同时提供")

    explicit_range = bool(args.plot_region) or args.plot_start is not None
    if explicit_range and len(genes) != 1:
        raise ValueError("自定义绝对绘图范围时请一次只输入一个基因")

    updated = [dict(g) for g in genes]

    if args.plot_region:
        parsed = parse_plot_region(args.plot_region)
        updated[0].update(parsed)
    elif args.plot_start is not None:
        if args.plot_start >= args.plot_end:
            raise ValueError("--plot-start 必须小于 --plot-end")
        updated[0]["plot_start"] = args.plot_start
        updated[0]["plot_end"] = args.plot_end

    if args.flank < 0:
        raise ValueError("--flank 不能为负数")
    if args.flank:
        for gene in updated:
            gene["plot_start"] = max(0, gene["plot_start"] - args.flank)
            gene["plot_end"] += args.flank

    return updated


_TRACKPLOT_HELP = None


def trackplot_supports(flag):
    """检查当前环境的 trackplot --help 是否支持某个新参数。"""
    global _TRACKPLOT_HELP
    if _TRACKPLOT_HELP is None:
        try:
            result = subprocess.run(
                ["trackplot", "--help"], capture_output=True, text=True, check=False
            )
            _TRACKPLOT_HELP = (result.stdout or "") + "\n" + (result.stderr or "")
        except OSError:
            _TRACKPLOT_HELP = ""
    return flag in _TRACKPLOT_HELP


def add_supported_flag(cmd, condition, flag):
    if not condition:
        return
    if trackplot_supports(flag):
        cmd.append(flag)
    else:
        print(f"    [WARN] 当前 trackplot 不支持 {flag}，已自动跳过")


def add_supported_value(cmd, value, flag):
    if value is None or str(value) == "":
        return
    if trackplot_supports(flag):
        cmd.extend([flag, str(value)])
    else:
        print(f"    [WARN] 当前 trackplot 不支持 {flag}，已自动跳过")


def format_region(gene, strand_mode="auto"):
    region = f"{gene['chr']}:{gene['plot_start']}-{gene['plot_end']}"

    if strand_mode in {"+", "-"}:
        return f"{region}:{strand_mode}"

    if strand_mode == "auto" and gene.get("strand") in {"+", "-"}:
        return f"{region}:{gene['strand']}"

    return region


def add_flag(cmd, condition, flag):
    if condition:
        cmd.append(flag)


def add_value(cmd, value, flag):
    if value is not None and str(value) != "":
        cmd.extend([flag, str(value)])


# ── plot ───────────────────────────────────────────────────────────────────────
def plot_gene(gene, outdir, samples, regions_tag, args):
    name   = gene["name"]
    region = format_region(gene, args.strand)

    os.makedirs(outdir, exist_ok=True)
    density_tsv = os.path.join(outdir, f"{name}_density.tsv")
    write_density_tsv(density_tsv, samples, args.library)

    intron_tag = f"fixed{args.fixed_intron_bp}bp" if args.fixed_intron_bp else f"intron{args.intron_scale}"
    range_tag = f"_{gene['chr']}-{gene['plot_start']}-{gene['plot_end']}"
    tx_tag = ""
    if args.transcripts_to_show:
        tx_tag = "_tx-" + args.transcripts_to_show.replace(",", "-").replace(" ", "")

    out_file = os.path.join(
        outdir,
        f"{name}_{args.normalize_format}_{regions_tag}_{intron_tag}{range_tag}{tx_tag}.{args.fmt}"
    )

    try:
        gtf_file = gene_gtf(name)
    except subprocess.CalledProcessError as e:
        print(f"  [ERROR] GTF prep failed for {name}: {e}")
        return False

    intron_scale = str(args.fixed_intron_bp) if args.fixed_intron_bp else str(args.intron_scale)

    cmd = [
        "trackplot",
        "-e", region,
        "-r", gtf_file,
        "--density", density_tsv,
        "--normalize-format", args.normalize_format,
        "-t", str(args.junction_threshold),
        "--fill-step", args.fill_step,
        "-o", out_file,
        "--dpi", str(args.dpi),
        "--width", str(int(args.width)),
        "--height", str(args.height),
        "--font-size", str(args.font_size),
        "--font", args.font,
        "--title", f"{name}  [{args.normalize_format.upper()} | {regions_tag}]",
        "--intron-scale", intron_scale,
        "--exon-scale", str(args.exon_scale),
        "--annotation-scale", str(args.annotation_scale),
        "--distance-ratio", str(args.distance_ratio),
        "--ref-color", args.ref_color,
        "-p", str(args.process),
    ]

    add_flag(cmd, args.show_junction_num, "--show-junction-num")
    add_supported_flag(cmd, args.show_mean_junction_num, "--show-mean-junction-num")
    add_supported_flag(
        cmd,
        args.junctions_on_top and not args.density_by_strand,
        "--junctions-on-top",
    )
    add_flag(cmd, args.same_y, "--same-y")
    add_flag(cmd, args.hide_y_label, "--hide-y-label")
    add_flag(cmd, args.raster, "--raster")
    add_flag(cmd, args.remove_empty, "--remove-empty")
    add_flag(cmd, args.choose_primary, "--choose-primary")
    add_flag(cmd, args.show_id, "--show-id")
    add_flag(cmd, args.show_exon_id, "--show-exon-id")
    add_flag(cmd, args.no_gene, "--no-gene")
    add_flag(cmd, args.density_by_strand, "--density-by-strand")

    add_value(cmd, args.transcripts_to_show, "--transcripts-to-show")
    add_value(cmd, args.included_junctions, "--included-junctions")
    add_value(cmd, args.sites, "--sites")
    add_value(cmd, args.log, "--log")
    add_supported_value(cmd, args.n_y_ticks, "--n-y-ticks")
    add_supported_flag(cmd, args.no_title, "--no-title")
    add_supported_flag(cmd, args.no_region_text, "--no-region-text")

    print(f"  plotting {name}  ({region})")
    print(f"    regions       : {regions_tag}")
    print(f"    intron-scale  : {intron_scale}")
    print(f"    exon-scale    : {args.exon_scale}")
    print(f"    annotation    : {args.annotation_scale}")
    print(f"    junction min  : {args.junction_threshold}")
    print(f"    junction nums : {args.show_junction_num}")
    print(f"    on top        : {args.junctions_on_top and not args.density_by_strand}")
    print(f"    output        : {out_file}")

    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  [ERROR] {name}\n{r.stderr[-2000:]}")
        print("  [CMD]")
        print("  " + " ".join(cmd))
        return False

    print(f"  ✓ {out_file}")
    return True


# ── main ───────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()

    if args.regions:
        active_std = resolve_regions(args.regions)
    else:
        active_std = [ALIAS.get(r, r) for r in _DEFAULT_REGIONS]

    samples = [{"name": n, "color": CLUSTER_PALETTE[n]} for n in active_std]
    regions_tag = args.regions.replace(",", "-") if args.regions else "-".join(_DEFAULT_REGIONS)

    print(f"Regions: {active_std}\n")
    print("Checking BAM files ...")
    missing = []
    for s in samples:
        bp = bam_path(s["name"])
        bai1 = f"{bp}.bai"
        bai2 = bp.replace(".bam", ".bai")
        exists = os.path.exists(bp)
        indexed = os.path.exists(bai1) or os.path.exists(bai2)

        print(f"  {'✓' if exists else '✗'} {s['name']}  ({bp})")
        if exists and not indexed:
            print(f"    [WARN] BAM index not found: {bai1} or {bai2}")
            print(f"           若 trackplot 报错，请先运行: samtools index {bp}")

        if not exists:
            missing.append(bp)

    if missing:
        print("\n[ERROR] Missing BAM files:")
        for m in missing:
            print(f"  {m}")
        sys.exit(1)

    if args.coords:
        genes = load_coords_file(args.coords)
    else:
        gene_names = (
            load_gene_info(args.gene_info) if args.gene_info
            else args.genes if args.genes
            else sorted(f[:-4] for f in os.listdir(COORDS_DIR) if f.endswith(".tsv"))
        )

        genes = []
        for gname in gene_names:
            per = os.path.join(COORDS_DIR, f"{gname}.tsv")
            if not os.path.exists(per):
                print(f"  [SKIP] {gname}: coords not found")
                continue
            loaded = load_coords_file(per)
            if loaded:
                genes.append(loaded[0])

    if args.gene_info:
        order = load_gene_info(args.gene_info)
        wanted = set(order)
        idx = {n: i for i, n in enumerate(order)}
        genes = sorted(
            [g for g in genes if g["name"] in wanted],
            key=lambda g: idx.get(g["name"], 9999)
        )
    elif args.genes and not args.coords:
        wanted = set(args.genes)
        genes = [g for g in genes if g["name"] in wanted]
        if not genes:
            print("[ERROR] genes not found in coords")
            sys.exit(1)

    try:
        genes = apply_range_overrides(genes, args)
    except ValueError as e:
        print(f"[ERROR] {e}")
        sys.exit(1)

    print(f"\ngenes to plot ({len(genes)}): {[g['name'] for g in genes]}\n")

    ok, fail = [], []
    for g in genes:
        if plot_gene(g, args.outdir, samples, regions_tag, args):
            ok.append(g["name"])
        else:
            fail.append(g["name"])

    print(f"\ndone — success: {len(ok)}  failed: {len(fail)}")
    if fail:
        print(f"failed: {fail}")


if __name__ == "__main__":
    main()
