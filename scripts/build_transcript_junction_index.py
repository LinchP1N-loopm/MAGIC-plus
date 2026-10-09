#!/usr/bin/env python3
import argparse
import os
import pickle
from collections import defaultdict


DEFAULT_GTF = "/histor/zhao/pangkun/database/ref_mm10_M23_release98/mm10.gtf"
DEFAULT_OUTPUT = "/histor/zhao/projects/lupeng/zhaolab/ZJJ/Reconstruct/transcript_junction_index.pkl"
VALID_CHROM = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}


def parse_attr(attr_str, key):
    for attr in attr_str.split(";"):
        attr = attr.strip()
        if attr.startswith(key + " "):
            return attr.split(" ", 1)[1].strip().strip('"')
    return ""


def build_index(gtf_path):
    transcript_exons = {}
    print(f"解析 GTF: {gtf_path}")

    with open(gtf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in ("transcript", "exon") or fields[0] not in VALID_CHROM:
                continue

            feature, chrom, strand, attrs = fields[2], fields[0], fields[6], fields[8]
            transcript_id = parse_attr(attrs, "transcript_id")
            transcript_name = parse_attr(attrs, "transcript_name") or transcript_id
            gene_name = parse_attr(attrs, "gene_name")

            if feature == "transcript":
                transcript_exons[transcript_name] = {
                    "transcript_id": transcript_id,
                    "gene": gene_name,
                    "chrom": chrom,
                    "strand": strand,
                    "exons": [],
                }
            elif transcript_name in transcript_exons:
                start = int(fields[3]) - 1
                end = int(fields[4])
                transcript_exons[transcript_name]["exons"].append((start, end))

    index = {}
    n_single_exon = 0
    for transcript_name, info in transcript_exons.items():
        exons = sorted(info["exons"])
        if len(exons) < 2:
            n_single_exon += 1
            continue
        index[transcript_name] = {
            "transcript_id": info["transcript_id"],
            "gene": info["gene"],
            "chrom": info["chrom"],
            "strand": info["strand"],
            "junctions": {(exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1)},
            "n_exons": len(exons),
        }

    print(f"转录本总数: {len(transcript_exons):,}")
    print(f"单exon转录本（跳过）: {n_single_exon:,}")
    print(f"有junction的转录本: {len(index):,}")
    print(f"含transcript_id: {sum(bool(x['transcript_id']) for x in index.values()):,}")
    return index


def build_lookup(index):
    lookup = defaultdict(list)
    for transcript_name, info in index.items():
        for junction in info["junctions"]:
            lookup[(info["chrom"], junction)].append(transcript_name)
    return dict(lookup)


def main():
    parser = argparse.ArgumentParser(description="Build transcript junction index from GTF")
    parser.add_argument("--gtf", default=DEFAULT_GTF)
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    index = build_index(args.gtf)
    output = {"index": index, "lookup": build_lookup(index)}
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "wb") as handle:
        pickle.dump(output, handle, protocol=pickle.HIGHEST_PROTOCOL)

    print(f"反向索引条目: {len(output['lookup']):,}")
    print(f"保存: {args.output}")


if __name__ == "__main__":
    main()
