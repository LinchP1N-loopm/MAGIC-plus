#!/bin/bash
# 00_get_coords_regions.sh
# Purpose : Look up gene coordinates from GTF by gene name,
#           then write TSV files used by sashimi/track plotting.
#
# Usage:
#   bash 00_get_coords_regions.sh Apoe Gfap
#   bash 00_get_coords_regions.sh --file gene_info/all_genes.txt
#   bash 00_get_coords_regions.sh --outfile gene_coords/all_coords.tsv Apoe Gfap
set -euo pipefail
# 路径与参数
PROJECT_DIR="/histor/zhao/projects/lupeng/zhaolab/ZJJ/Fig4_sashimi"
GTF="/histor/zhao/pangkun/database/ref_mm10_M23_release98/mm10.gtf"
COORDS_DIR="${PROJECT_DIR}/gene_coords"
OUTFILE=""
PADDING=1000
# 可选：为特定基因手动指定绘图范围
# 格式：CUSTOM["Gene"]="start,end"
declare -A CUSTOM
# CUSTOM["Mbp"]="82532000,82586000"
# 参数解析
INPUT_FILE=""
GENES=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --file)
            INPUT_FILE="$2"
            shift 2
            ;;
        --outfile)
            OUTFILE="$2"
            shift 2
            ;;
        --padding)
            PADDING="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage:"
            echo "  bash $0 Apoe Gfap"
            echo "  bash $0 --file gene_info/all_genes.txt"
            echo "  bash $0 --outfile gene_coords/all_coords.tsv Apoe Gfap"
            echo "  bash $0 --padding 500 Apoe"
            exit 0
            ;;
        *)
            GENES+=("$1")
            shift
            ;;
    esac
done
# 从文件读取基因名
#
# 默认兼容原脚本：
#   多列表格优先取第4列
#   单列文件直接读取第1列
if [[ -n "${INPUT_FILE}" ]]; then
    if [[ ! -f "${INPUT_FILE}" ]]; then
        echo "[ERROR] File not found: ${INPUT_FILE}"
        exit 1
    fi
    echo "Reading gene names from: ${INPUT_FILE}"
    mapfile -t GENES < <(
        awk -F'\t' '
            NR==1 {
                if (NF==1 && $1 !~ /^(gene|gene_name|Gene_Name)$/) print $1
                next
            }
            /^#/ {
                next
            }
            NF>=4 && $4!="" {
                print $4
                next
            }
            NF>=1 && $1!="" {
                print $1
            }
        ' "${INPUT_FILE}" | sort -u
    )
    echo "  → ${#GENES[@]} unique genes loaded"
fi
if [[ ${#GENES[@]} -eq 0 ]]; then
    echo "[ERROR] No genes supplied."
    echo "Usage: bash $0 Apoe Gfap"
    echo "       bash $0 --file gene_info/all_genes.txt"
    exit 1
fi
# Ensembl ID → gene_name
resolve_id() {
    local query="$1"
    if [[ "${query}" =~ ^ENSMUSG ]]; then
        local name
        name=$(grep -m1 "gene_id \"${query}\"" "${GTF}" \
            | grep -oP 'gene_name "\K[^"]+' \
            | head -1 || true)
        echo "${name:-${query}}"
    else
        echo "${query}"
    fi
}
# 输出准备
mkdir -p "${COORDS_DIR}"
if [[ -n "${OUTFILE}" ]]; then
    mkdir -p "$(dirname "${OUTFILE}")"
    echo -e \
        "gene\tchr\tgene_start\tgene_end\tstrand\tplot_start\tplot_end" \
        > "${OUTFILE}"
    echo "Output mode: single file → ${OUTFILE}"
else
    echo "Output mode: per-gene → ${COORDS_DIR}/<gene>.tsv"
fi
echo "GTF     : ${GTF}"
echo "Project : ${PROJECT_DIR}"
echo "Padding : ${PADDING} bp"
echo "Genes   : ${#GENES[@]}"
echo ""
printf "  %-18s %-7s %-13s %-13s %-8s %-13s %-13s %s\n" \
    "gene" \
    "chr" \
    "gene_start" \
    "gene_end" \
    "strand" \
    "plot_start" \
    "plot_end" \
    "note"
echo "  $(printf '%0.s─' {1..90})"
FOUND=0
NOTFOUND=0
# 查询并输出坐标
for raw_gene in "${GENES[@]}"; do
    gene=$(resolve_id "${raw_gene}")
    line=$(grep -w "gene_name \"${gene}\"" "${GTF}" \
        | awk '$3=="gene"' \
        | head -1 || true)
    if [[ -z "${line}" ]]; then
        echo "  [NOT FOUND] ${raw_gene}"
        NOTFOUND=$((NOTFOUND + 1))
        continue
    fi
    chr=$(echo "${line}" | awk '{print $1}')
    gstart=$(echo "${line}" | awk '{print $4}')
    gend=$(echo "${line}" | awk '{print $5}')
    strand=$(echo "${line}" | awk '{print $7}')
    if [[ -n "${CUSTOM[$gene]+_}" ]]; then
        plot_start=$(echo "${CUSTOM[$gene]}" | cut -d, -f1)
        plot_end=$(echo "${CUSTOM[$gene]}" | cut -d, -f2)
        note="[custom region]"
    else
        plot_start=$((gstart - PADDING))
        plot_end=$((gend + PADDING))
        [[ "${plot_start}" -lt 1 ]] && plot_start=1
        note=""
    fi
    printf "  %-18s %-7s %-13s %-13s %-8s %-13s %-13s %s\n" \
        "${gene}" \
        "${chr}" \
        "${gstart}" \
        "${gend}" \
        "${strand}" \
        "${plot_start}" \
        "${plot_end}" \
        "${note}"
    if [[ -n "${OUTFILE}" ]]; then
        echo -e \
            "${gene}\t${chr}\t${gstart}\t${gend}\t${strand}\t${plot_start}\t${plot_end}" \
            >> "${OUTFILE}"
    else
        local_out="${COORDS_DIR}/${gene}.tsv"
        echo -e \
            "gene\tchr\tgene_start\tgene_end\tstrand\tplot_start\tplot_end" \
            > "${local_out}"
        echo -e \
            "${gene}\t${chr}\t${gstart}\t${gend}\t${strand}\t${plot_start}\t${plot_end}" \
            >> "${local_out}"
    fi
    FOUND=$((FOUND + 1))
done
echo ""
echo "Found: ${FOUND}   Not found: ${NOTFOUND}"
if [[ -n "${OUTFILE}" ]]; then
    echo "→ written to: ${OUTFILE}"
else
    echo "→ written to: ${COORDS_DIR}/<gene>.tsv (${FOUND} files)"
fi
echo ""
echo "Next:"
echo "python 02_plot_sashimi_regions.py <gene> --region CA1 --samples G21-4,G21-1,G21-2"
