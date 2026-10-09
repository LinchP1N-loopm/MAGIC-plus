<p align="center">
  <img src="assets/github_repo_icon.png" alt="MAGIC-plus" width="200">
</p>

# MAGIC-plus

Custom scripts and analysis notebooks for the MAGIC-plus study, covering sequencing data processing, transcript reconstruction, alternative splicing, RNA editing, and spatial transcriptomic analysis.

## Repository contents

| Directory | Contents |
| --- | --- |
| `scripts/` | Sequencing processing, UMI-based reconstruction, transcript matching, rMATS analyses, RNA editing, and coverage/sashimi plotting scripts |
| `Figure4/` | Figure 4 notebooks for differential expression, alternative-splicing enrichment and networks, and spatial RNA/PSI integration |
| `analysis/` | Combined RNA and alternative-splicing analysis notebook |
| `data/` | Supplied mouse reference annotation table (`TABLE1_mm10_v3.txt`) |
| `assets/` | Repository image |

See [file_index.tsv](file_index.tsv) for the complete file list. Figure assignments follow the supplied filenames; additional main and supplementary figure folders can be added when their scripts are available.

## Requirements

These scripts were developed for Linux and Jupyter environments. Shell and PBS scripts include study-specific paths, sample names, environment names, and cluster settings. Update these before running. Reuse the appropriate mamba/micromamba environments where available.

Dependencies vary by workflow:

- Sequencing and QC: STAR/STARsolo, fastp, mawk, pigz, samtools, st_pipeline, and RSeQC.
- Alternative splicing and reconstruction: rMATS, pysam, pandas, and joblib.
- RNA editing: GATK, samtools, cgranges, anndata, and pandas.
- Notebook analysis: Jupyter, NumPy, pandas, SciPy, Scanpy, AnnData, scikit-learn, PyTorch, SpatialGlue, Matplotlib, Seaborn, NetworkX, statsmodels, GSEApy, and tqdm.
- Coverage and junction plots: deepTools, pyGenomeTracks, and Trackplot.

This is an initial code archive. A complete, version-pinned environment has not yet been supplied, and the workflows have not been independently rerun as part of repository preparation.

## Usage

1. Download or clone this repository.
2. Prepare the inputs required by the selected workflow: FASTQ/BAM files, genome references and annotations, barcode/coordinate files, processed AnnData objects, or rMATS outputs.
3. Edit the input/output paths, sample identifiers, reference paths, and tool/environment settings in the relevant scripts or notebook cells. For PBS jobs, also update queue, resource, and working-directory directives.
4. Run the chosen notebook in Jupyter, or submit the appropriate job on a PBS cluster. These files represent separate workflows and should not be run together in alphabetical order.

Examples of entry points:

| Workflow | Files |
| --- | --- |
| Sequencing processing | `scripts/St_pipeline.sh`, `scripts/MAGIC-seq1.RD.sh` |
| UMI reconstruction and transcript matching | `scripts/umi_reconstructV5.py`, `scripts/build_transcript_junction_index.py`, `scripts/match_transcript.py` |
| Alternative splicing | `scripts/step3_control_vs_AD_regions.pbs`, `scripts/step3_run_rmats_G21-4_rm_dup_ALL.sh` |
| RNA editing | `scripts/01_export_maps.pbs`, `scripts/03_pipeline_G21-4_0515.pbs` |
| Coverage and sashimi plots | `scripts/00_get_coords_regions.sh`, `scripts/02_plot_sashimi_junctions.py`, `scripts/02_plot_tracksV2.py` |
| Figure 4 analysis | Notebooks in `Figure4/` |

The original analysis code is retained. Notebook outputs and execution counts are cleared to keep the repository compact. Source filenames, including Chinese filenames, are preserved. Moving files into this layout does not automatically update the original absolute paths.

## Data and citation

Raw sequencing data and large processed datasets are not included. Data accession numbers and the final manuscript citation will be added when available. The provenance and version of the supplied annotation table should be documented before the final manuscript release.

Related method repository: [MAGIC-seq](https://github.com/bioinfo-biols/MAGIC-seq).

## Contact

Please use the repository's GitHub Issues for questions about this code.
