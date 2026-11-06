#!/bin/bash
set -euo pipefail

# Parameters
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"
COORD_SCRIPT="${SCRIPT_DIR}/00_get_coords_regions.sh"
PLOT_SCRIPT="${SCRIPT_DIR}/02_plot_sashimi_junctions.py"
OUT_DIR="${PROJECT_DIR}/sashimi_junctions/AS_gene_region"
INTRON_SCALE=0.05
JUNCTION_THRESHOLD=10

AD_EVENTS=("Tmem214|TH" "Ctss|CA3" "Eno2|TH" "Bin1|VS" "Mbp|CA1" "Ncam1|FT" "Mag|FT" "Cd74|CTX_3")
AGED_EVENTS=("Atp5j2|CA3" "Ncam1|FT" "Bin1|VS" "Clcn4|TH" "Mag|FT" "Mbp|CA1" "C4b|CTX_3")

# 生成坐标
bash "${COORD_SCRIPT}" Tmem214 Ctss Eno2 Bin1 Mbp Ncam1 Mag Cd74 Atp5j2 Clcn4 C4b

# 绘图
run_events() {
    local comparison="$1" samples="$2"
    shift 2
    local outdir="${OUT_DIR}/${comparison}"
    mkdir -p "${outdir}"

    for event in "$@"; do
        IFS='|' read -r gene region <<< "${event}"
        python "${PLOT_SCRIPT}" "${gene}" --region "${region}" --samples "${samples}" \
            --intron-scale "${INTRON_SCALE}" --junction-threshold "${JUNCTION_THRESHOLD}" \
            --outdir "${outdir}"
    done
}

run_events "AD_vs_Young" "G21-4,G21-2" "${AD_EVENTS[@]}"
run_events "Aged_vs_Young" "G21-4,G21-1" "${AGED_EVENTS[@]}"

echo "Done: ${OUT_DIR}"
