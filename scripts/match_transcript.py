#!/usr/bin/env python3
import argparse
import collections
import os
import pickle
import re


_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def get_junctions_from_cigar(start_1based, cigar_str):
    pos = int(start_1based) - 1
    junctions = []
    for length_text, op in _CIGAR_RE.findall(str(cigar_str)):
        length = int(length_text)
        if op in ("M", "=", "X", "D"):
            pos += length
        elif op == "N":
            junctions.append((pos, pos + length))
            pos += length
    return set(junctions)


def match_transcripts(chrom, read_junctions, transcript_index, lookup):
    candidates = set()
    for junction in read_junctions:
        candidates.update(lookup.get((chrom, junction), []))
    return sorted(t for t in candidates if transcript_index[t]["chrom"] == chrom and read_junctions.issubset(transcript_index[t]["junctions"]))


def get_transcript_id(info):
    for key in ("transcript_id", "ensembl_transcript_id", "transcriptId", "id"):
        if info.get(key):
            return str(info[key])
    return ""


def main():
    parser = argparse.ArgumentParser(description="Build a spot-by-transcript count matrix from reconstructed UMI CIGARs")
    parser.add_argument("--tsv", required=True)
    parser.add_argument("--index", required=True)
    parser.add_argument("--adata", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--obs-cols", default="AS_cluster,brain_cluster,leiden")
    parser.add_argument("--min-spots", type=int, default=2)
    parser.add_argument("--chunk-size", type=int, default=500000)
    args = parser.parse_args()

    import anndata as ad
    import pandas as pd
    from scipy.sparse import csr_matrix
    from tqdm import tqdm

    with open(args.index, "rb") as handle:
        data = pickle.load(handle)
    transcript_index = data["index"]
    lookup = data["lookup"]

    reference = ad.read_h5ad(args.adata)
    if "in_tissue" in reference.obs.columns:
        reference = reference[reference.obs["in_tissue"].astype(str) == "1"].copy()

    obs_cols = [x.strip() for x in args.obs_cols.split(",") if x.strip()]
    missing = [x for x in obs_cols if x not in reference.obs.columns]
    if missing:
        raise KeyError(f"Missing obs columns: {missing}")

    obs_names = reference.obs_names.astype(str).tolist()
    valid_barcodes = set(obs_names)
    obs_to_index = {name: i for i, name in enumerate(obs_names)}
    valid_chroms = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}
    counts = collections.Counter()
    stats = collections.Counter()

    reader = pd.read_csv(args.tsv, sep="\t", usecols=["cb", "gn", "chrom", "strand", "new_cigar_start", "new_cigar"], chunksize=args.chunk_size, dtype={"cb": str, "gn": str, "chrom": str, "strand": str, "new_cigar_start": "Int64", "new_cigar": str})
    for chunk in tqdm(reader, desc="Processing TSV", unit="chunk"):
        stats["total"] += len(chunk)
        has_junction = chunk["new_cigar"].str.contains("N", na=False)
        stats["no_junction"] += int((~has_junction).sum())
        chunk = chunk[has_junction & chunk["chrom"].isin(valid_chroms)]
        valid_cb = chunk["cb"].isin(valid_barcodes)
        stats["invalid_barcode"] += int((~valid_cb).sum())
        for row in chunk[valid_cb].itertuples(index=False):
            junctions = get_junctions_from_cigar(row.new_cigar_start, row.new_cigar)
            matched = match_transcripts(row.chrom, junctions, transcript_index, lookup)
            if not matched:
                stats["no_match"] += 1
                continue
            stats["matched"] += 1
            stats["ambiguous"] += len(matched) > 1
            for transcript in matched:
                counts[(row.cb, transcript)] += 1

    transcripts = sorted({transcript for _, transcript in counts})
    transcript_to_index = {name: i for i, name in enumerate(transcripts)}
    rows, cols, values = [], [], []
    for (barcode, transcript), count in counts.items():
        rows.append(obs_to_index[barcode])
        cols.append(transcript_to_index[transcript])
        values.append(count)

    matrix = csr_matrix((values, (rows, cols)), shape=(len(obs_names), len(transcripts)))
    keep = (matrix > 0).sum(axis=0).A1 >= args.min_spots
    matrix = matrix[:, keep]
    transcripts = [name for name, selected in zip(transcripts, keep) if selected]

    var = pd.DataFrame({
        "transcript_name": transcripts,
        "gene_name": [transcript_index.get(t, {}).get("gene", "") for t in transcripts],
        "chrom": [transcript_index.get(t, {}).get("chrom", "") for t in transcripts],
        "strand": [transcript_index.get(t, {}).get("strand", "") for t in transcripts],
        "n_junctions": [len(transcript_index.get(t, {}).get("junctions", set())) for t in transcripts],
        "transcript_id": [get_transcript_id(transcript_index.get(t, {})) for t in transcripts],
    }, index=transcripts)

    result = ad.AnnData(X=matrix, obs=reference.obs[obs_cols].copy(), var=var)
    result.layers["counts"] = result.X.copy()
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    result.write_h5ad(args.output, compression="gzip")

    print(f"Saved: {args.output}")
    print(f"Matrix: {result.n_obs} spots x {result.n_vars} transcripts")
    print(f"Total UMI: {stats['total']:,}; matched: {stats['matched']:,}; ambiguous: {stats['ambiguous']:,}")


if __name__ == "__main__":
    main()
