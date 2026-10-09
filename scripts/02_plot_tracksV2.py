#!/usr/bin/env python3
"""
02_plot_tracksV2.py
Usage:
    python 02_plot_tracksV2.py Penk Gfap Ttr     # specific genes
    python 02_plot_tracksV2.py                    # all genes in gene_coords.tsv
"""

import os, sys, csv, subprocess

# ── paths ──────────────────────────────────────────────────────────────────────
GTF        = "/histor/zhao/pangkun/database/ref_mm10_M23_release98/mm10.gtf"
COORDS_TSV = "./gene_coords.tsv"
BIGWIG_DIR = "./bigwig_data"
OUTDIR     = "./tracks_output"
INI_DIR    = "./tracks_ini"
GTF_DIR    = "./tracks_gtf"

for d in (OUTDIR, INI_DIR, GTF_DIR):
    os.makedirs(d, exist_ok=True)

# ── samples ────────────────────────────────────────────────────────────────────
SAMPLES = [
    {"name": "MAGIC_plus", "color": "#C62828"},
    {"name": "MAGIC-seq",  "color": "#2979FF"},
]

# ── figure settings ────────────────────────────────────────────────────────────
WIDTH     = 24    # cm
DPI       = 300
H_BW      = 2.5   # cm, bigwig track height
H_GENE    = 0.7   # cm, gene structure height
FONTSIZE  = 7     # pt, applied everywhere except sample title
# Note: pyGenomeTracks uses matplotlib fonts.
# Arial requires Arial to be installed on the system.
# If unavailable, matplotlib will fall back to DejaVu Sans silently.
FONT      = "Arial"


# ── extract per-gene GTF (grep once, cached) ──────────────────────────────────
def gene_gtf(gene_name):
    out = f"{GTF_DIR}/{gene_name}.gtf"
    if not os.path.exists(out) or os.path.getsize(out) == 0:
        subprocess.run(
            f'grep -w \'gene_name "{gene_name}"\' "{GTF}" > "{out}"',
            shell=True, check=True
        )
    return out


# ── write ini and run pyGenomeTracks ──────────────────────────────────────────
def plot_gene(gene):
    name   = gene["name"]
    chrom  = gene["chr"]
    start  = gene["plot_start"]
    end    = gene["plot_end"]
    region = f"{chrom}:{start}-{end}"

    ini = f"{INI_DIR}/{name}.ini"
    with open(ini, "w") as f:

        # x-axis top — fontsize controls the bp coordinate label size
        f.write( "[x-axis]\n")
        f.write( "where = top\n")
        f.write(f"fontsize = {FONTSIZE}\n\n")

        for s in SAMPLES:
            f.write(f"[{s['name']}]\n")
            f.write(f"file      = {BIGWIG_DIR}/{s['name']}.bw\n")
            f.write( "file_type = bigwig\n")
            f.write(f"title     = {s['name']}\n")
            f.write(f"height    = {H_BW}\n")
            f.write(f"color     = {s['color']}\n")
            f.write( "min_value = 0\n")
            f.write( "type      = fill\n")
            f.write( "alpha     = 0.88\n")
            f.write( "nans_to_zeros = true\n")
            # f.write( "transform = log1p\n")
            # f.write( "y_axis_values = original\n")
            f.write(f"fontsize  = {FONTSIZE}\n")
            f.write( "\n[spacer]\nheight = 0.2\n\n")

        f.write( "[genes]\n")
        f.write(f"file      = {gene_gtf(name)}\n")
        f.write( "file_type = gtf\n")
        f.write(f"title     = {name}\n")
        f.write(f"height    = {H_GENE}\n")
        f.write( "style     = UCSC\n")
        f.write( "display   = collapsed\n")
        f.write( "merge_transcripts = true\n")
        f.write( "prefered_name     = gene_name\n")
        f.write( "color             = black\n")
        f.write( "color_arrow       = black\n")
        f.write( "color_backbone    = black\n")
        f.write( "labels            = false\n")
        f.write(f"fontsize          = {FONTSIZE}\n")
        f.write( "arrow_interval    = 10\n")
        f.write( "arrowhead_fraction = 0.01\n")
        f.write( "arrowhead_included = true\n")
        f.write( "\n[x-axis]\n")
        f.write( "where = bottom\n")
        f.write(f"fontsize = {FONTSIZE}\n")

    out_pdf = f"{OUTDIR}/{name}_tracks.pdf"
    cmd = [
        "pyGenomeTracks",
        "--tracks",   ini,
        "--region",   region,
        "--width",    str(WIDTH),
        "--dpi",      str(DPI),
        "--fontSize", str(FONTSIZE),  # global fallback font size
        "-o",         out_pdf,
    ]
    print(f"  plotting {name}  ({region})")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  [ERROR]\n{r.stderr[-1000:]}")
        return False
    print(f"  ✓ {out_pdf}")
    return True


# ── main ──────────────────────────────────────────────────────────────────────
def main():
    # Set Arial globally via matplotlib before pyGenomeTracks runs.
    # This works because pyGenomeTracks uses matplotlib internally.
    try:
        import matplotlib
        matplotlib.rcParams["font.family"] = "sans-serif"
        matplotlib.rcParams["font.sans-serif"] = [FONT, "DejaVu Sans"]
        print(f"font set to: {FONT} (fallback: DejaVu Sans)")
    except Exception as e:
        print(f"[WARNING] could not set font: {e}")

    genes = []
    with open(COORDS_TSV) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            genes.append({
                "name":       row["gene"],
                "chr":        row["chr"],
                "plot_start": int(row["plot_start"]),
                "plot_end":   int(row["plot_end"]),
            })

    if sys.argv[1:]:
        wanted = {a.upper() for a in sys.argv[1:]}
        genes  = [g for g in genes if g["name"].upper() in wanted]
        if not genes:
            print("[ERROR] none of the requested genes found in gene_coords.tsv")
            sys.exit(1)

    print(f"genes to plot: {[g['name'] for g in genes]}\n")

    ok, fail = [], []
    for g in genes:
        (ok if plot_gene(g) else fail).append(g["name"])

    print(f"\ndone — success: {len(ok)}  failed: {len(fail)}")
    if fail:
        print(f"failed: {fail}")


if __name__ == "__main__":
    main()
