#!/bin/bash
set -euo pipefail

# Parameters
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"
PLOT_SCRIPT="${SCRIPT_DIR}/02_plot_sashimi_junctions.py"
COORD_SCRIPT="${SCRIPT_DIR}/00_get_coords_regions.sh"
OUT_DIR="${PROJECT_DIR}/sashimi_junctions/selected_events"
JUNCTION_THRESHOLD=10
INTRON_SCALE=0.05
LIBRARY_TYPE="frf"  # 请按实际建库方向确认

# 检查文件
[[ -f "${PLOT_SCRIPT}" ]] || { echo "[ERROR] Missing: ${PLOT_SCRIPT}"; exit 1; }
[[ -f "${COORD_SCRIPT}" ]] || { echo "[ERROR] Missing: ${COORD_SCRIPT}"; exit 1; }
mkdir -p "${OUT_DIR}"

# 生成基因坐标
bash "${COORD_SCRIPT}" Pabpc4 Vps29 Pkm Bin1 U2af1 Cnbp

# 绘图
python "${PLOT_SCRIPT}" Pabpc4 --region CTX_3 --samples G21-4,G21-1 --plot-region chr4:123294000-123299800:+ --transcripts-to-show Pabpc4-202,Pabpc4-204 --junction-threshold "${JUNCTION_THRESHOLD}" --intron-scale "${INTRON_SCALE}" --library-type "${LIBRARY_TYPE}" --outdir "${OUT_DIR}"

python "${PLOT_SCRIPT}" Vps29 --region CA1 --samples G21-4,G21-1,G21-2 --plot-region chr5:122353300-122360800:+ --transcripts-to-show Vps29-204,Vps29-208 --junction-threshold "${JUNCTION_THRESHOLD}" --intron-scale "${INTRON_SCALE}" --library-type "${LIBRARY_TYPE}" --outdir "${OUT_DIR}"

python "${PLOT_SCRIPT}" Pkm --region FT --samples G21-4,G21-1,G21-2 --plot-region chr9:59669100-59680400:+ --transcripts-to-show Pkm-202,Pkm-201 --junction-threshold "${JUNCTION_THRESHOLD}" --intron-scale "${INTRON_SCALE}" --library-type "${LIBRARY_TYPE}" --outdir "${OUT_DIR}"

python "${PLOT_SCRIPT}" Bin1 --region CTX_2 --samples G21-4,G21-2 --transcripts-to-show Bin1-201,Bin1-205 --junction-threshold "${JUNCTION_THRESHOLD}" --intron-scale "${INTRON_SCALE}" --library-type "${LIBRARY_TYPE}" --outdir "${OUT_DIR}"

python "${PLOT_SCRIPT}" U2af1 --region DG --samples G21-4,G21-2 --transcripts-to-show U2af1-201,U2af1-202 --junction-threshold "${JUNCTION_THRESHOLD}" --intron-scale "${INTRON_SCALE}" --library-type "${LIBRARY_TYPE}" --outdir "${OUT_DIR}"

python "${PLOT_SCRIPT}" Cnbp --region DG --samples G21-4,G21-2 --transcripts-to-show Cnbp-208,Cnbp-203 --junction-threshold "${JUNCTION_THRESHOLD}" --intron-scale "${INTRON_SCALE}" --library-type "${LIBRARY_TYPE}" --outdir "${OUT_DIR}"

echo "Done: ${OUT_DIR}"
