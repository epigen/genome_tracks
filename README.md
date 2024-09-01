[![DOI](https://zenodo.org/badge/438573546.svg)](https://zenodo.org/doi/10.5281/zenodo.10849097)

# Genome Browser Track Visualization Workflow 
A [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow for easy visualization of genome browser tracks of aligned/mapped BAM files (e.g., RNA-seq, ATAC-seq, scRNA-seq, ...) powered by the wrapper [gtracks](https://gitlab.com/salk-tm/gtracks) for the package [pyGenomeTracks](https://github.com/deeptools/pyGenomeTracks) and [IGV-reports](https://github.com/igvteam/igv-reports).

This workflow adheres to the module specifications of [MR.PARETO](https://github.com/epigen/mr.pareto), an effort to augment research by modularizing (biomedical) data science. For more details, instructions and modules check out the project's repository. Please consider starring and sharing modules that are interesting or useful to you, this helps me in prioritizing my efforts!

**If you use this workflow in a publication, please don't forget to give credits to the authors by citing this DOI [10.5281/zenodo.10849097](https://zenodo.org/doi/10.5281/zenodo.10849097).**

![Workflow Rulegraph](./workflow/dags/rulegraph.svg)

Table of contents
----------------
  * [Authors](#authors)
  * [Software](#software)
  * [Methods](#methods)
  * [Features](#features)
  * [Usage](#usage)
  * [Configuration](#configuration)
  * [Examples](#examples)
  * [Links](#links)
  * [Genome Browser Tracks](#genome-browser-tracks)
  * [Resources](#resources)
  * [Publications](#publications)

# Authors
- [Stephan Reichl](https://github.com/sreichl)
- [Christoph Bock](https://github.com/chrbock)

# Software
This project wouldn't be possible without the following software and their dependencies:

| Software | Reference (DOI) |
| :---: | :---: |
| deeptools | https://doi.org/10.1093/nar/gkw257 |
| gtracks | https://gitlab.com/salk-tm/gtracks |
| igv-reports | https://github.com/igvteam/igv-reports |
| pygenometracks | https://doi.org/10.1093/bioinformatics/btaa692 |
| samtools | https://doi.org/10.1093/bioinformatics/btp352 |
| sinto | https://github.com/timoast/sinto |

# Methods
This is a template for the Methods section of a scientific publication and is intended to serve as a starting point. Only retain paragraphs relevant to your analysis. References [ref] to the respective publications are curated in the software table above. Versions (ver) have to be read out from the respective conda environment specifications (workflow/envs/\*.yaml file) or post execution in the result directory (/envs/genome_tracks/\*.yaml). Parameters that have to be adapted depending on the data or workflow configurations are denoted in squared brackets e.g., [X].

__(optional) Single-cell preprocessing.__ Each single-cell BAM file was split into [group]-wise BAM files according to it's cell barcode metadata using filterbarcodes from the command line tool sinto (ver) [ref].

__Processing.__ Aligned (filtered, and indexed) BAM files were merged by [group] using samtools (ver) [ref]. Each merged BAM file's coverage was determined for dowmnstream analysis and visualization using bamCoverage from the command line tool deepTools (ver) [ref] and saved in the bigWig format. Finally, we extracted coordinates, extended start and end by [base_buffer] bases, and number of isoforms of all relevant genes/genomic regions [gene_list] from the 12 column BED file genome [genome] annotation [genome_bed].

__Visualization.__ Visualizations for each relevant gene/genomic region and [category] were generated by using the generated bigWig coverage files and vertically stacking genome browser tracks with their annotation at the [x_axis] and each track scaled by [y_max] reads. The plotting was performed using the python wrapper gtracks (ver) [ref] for the package pyGenomeTracks (ver) [ref]. Additionally, an interactive self-contained IGV-report containing all merged samples and gene/genomic regions of interest was generated using igv-reports (ver) [ref]. Finally, a UCSC genome browser track hub was created for online sharing and inspection using [UCSC Genome Browser](https://genome.ucsc.edu/). Both the plotted tracks and the UCSC genome browser tracks were color coded according to [group].

**The processing and visualizations described here were performed using a publicly available Snakemake [ver] (ref) workflow [[10.5281/zenodo.10849097](https://zenodo.org/doi/10.5281/zenodo.10849097)].**

# Features
The workflow performs the following steps to produce the outlined results (`genome_tracks/`).

- Processing
  0. (optional) Single-cell BAM files are split into groups according to their metadata file using [sinto::filterbarcodes](https://timoast.github.io/sinto/basic_usage.html#filter-cell-barcodes-from-bam-file) and saved in folders named by the md5 hash of the input BAM file path (`sc_bams/{md5hash}/{group}.bam`). Downstream they are processed and visualized the same as bulk samples.
  1. BAM files of the same group are merged and indexed using [samtools::merge](https://www.htslib.org/doc/samtools-merge.html). (`merged_bams/{group}.bam`)
  2. A bigWig file per merged BAM file is generated using [deepTools::bamCoverage](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html). (`bigWigs/{group}.bw`)
  3. Annotations per gene are retrieved from the 12-column BED file (not necessary for genomic regions).
      - coordinates from the 12-column BED file are extracted and extended at start/end by the parameter base_buffer.
      - the number of isoforms i.e. number of lines in the BED file is determined (__only for genes, for genomic regions it is hardcoded to 1__) to plot below the tracks.
      - genes that were not found in the provided genome BED file are reported (`genes_not_found.csv`).
  4. A BED file of all genes and genomic regions is generated for the IGV-report (`genes.bed`).
- Visualization
  - One plot including all groups visualizing the coverage (i.e., bigWigs) per gene/genomic region with the determined gene-parameters i.e, coordinates and gene-rows is generated (`tracks/{gene|region}.{svg|pdf|png}`). The track height (`ymax`) and colors (`track_colors`) are configurable within the `gene_list` or config file, respectively. Default color is black.
  - An interactive self-contained IGV-report from merged BAM files with all genes/regions is generated for inspection for sharing (`igv-report.html`).
  - A UCSC genome browser track hub for all bigWigs is set up (`bigWigs/`), color coded according to configuration. See detailed instructions for usage and sharing below.

# Usage
Here are some tips for the usage of this workflow:
- Start with the 1-5 most interesting genes (e.g., marker genes for cell types as quality control or the most differentially expressed between conditions) and few/relevant samples for a test run.
- Set y-max to auto (i.e., 0 in the `gene_list` CSV table) for the first run to get a feeling for the magnitudes in your samples/groups and adapt to the highest peaks afterward to make the tracks comparable within a gene/region of interest.
- Splitting single-cell BAM files, merging BAM files and generating bigWig files take the longest, but are performed only once (multithreaded 4x configured threads). The plot generation for different genes/genomic regions afterward is very fast. Therefore, it is recommended to get the workflow going with a few samples and then increase the number of samples and genes/regions afterward.
- For simple quality control (QC) of all/individual samples just provide a unique sample name in the `group` column of the annnotation file. Then no samples will be merged into groups, but only renamed (useful before sharing).

This workflow is written with Snakemake and its usage is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog?usage=epigen/genome_tracks).

# Configuration
Detailed specifications can be found here [./config/README.md](./config/README.md)

# Examples
--- COMING SOON ---

Runtime examples for different data modalities:
- 78 ATAC-seq samples in 31 groups took 22 minutes with max 4 cores and 4GB memory.
- 64 RNA-seq samples in 31 groups took 16 minutes with max 4 cores and 4GB memory.
- 2 10x genomics 5' scRNA-seq samples/reactions each ~10k cells split in 2 small subset groups took 31 minutes with max 4 cores and 8GB memory.

# Genome Browser Tracks
The `bigWigs` directory contains the read coverage per sample/group in bigWig format (`{group}.bw`) for visual inspection of each sample e.g., during QC or group e.g., comparison of conditions. Below are instructions for two different approaches (online/local).

## UCSC Genome Browser Track Hub (online)
0. Requirement: web server.
1. Copy (or symlink) the `bigWigs` directory to an externally __accessible__ location on your web server (=`web_server_location`).
2. Create a UCSC Genome Browser hyperlink
    - the general formula is: ucsc_url + genome + web_server_location + bigWigs/hub.txt
    - concretely: `http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=` + {genome} + `&hubUrl=` + {web_server_location} + `bigWigs/hub.txt`
    - mm10 example: [http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=mm10&hubUrl=https://medical-epigenomics.org/data/genome_tracks/mm10test/hub/hub.txt](http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=mm10&hubUrl=https://medical-epigenomics.org/data/genome_tracks/mm10test/hub/hub.txt)
    - hg38 example: [http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg38&hubUrl=https://medical-epigenomics.org/data/genome_tracks/hg38test/hub/hub.txt](http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg38&hubUrl=https://medical-epigenomics.org/data/genome_tracks/hg38test/hub/hub.txt)
3. Share the link with the world e.g., collaborators or upon publication of your data.

**A new feature (2024-08-30)** allows users to download all visible data in the current region directly from our tracks display. This facilitates reproducibility when writing reports or publications as data can update and change over time. This feature can be found in the blue bar menu by going to **Downloads > Download Current Track Data**. The resulting pop-up dialogue box (see screenshot below) can configure the exact tracks to download from all visible tracks, as well as the file name and the output format (JSON, csv, tsv).
![UCSC_download](https://github.com/user-attachments/assets/174ad904-52a2-458f-a1e7-387a0ff95007)

## IGV: Integrative Genomics Viewer (local/offline)
0. Requirement: [IGV Desktop application](https://igv.org/doc/desktop/).
1. Open IGV.
2. Select genome.
3. Drag and drop all/selected bigWig files from the `bigWigs` directory directly into the IGV application.

# Links
- [GitHub Repository](https://github.com/epigen/genome_tracks/)
- [GitHub Page](https://epigen.github.io/genome_tracks/)
- [Zenodo Repository](https://zenodo.org/doi/10.5281/zenodo.10849097)
- [Snakemake Workflow Catalog Entry](https://snakemake.github.io/snakemake-workflow-catalog?usage=epigen/genome_tracks)

# Resources
- Recommended compatible [MR.PARETO](https://github.com/epigen/mr.pareto) modules:
  - for upstream processing (before)
    - [ATAC-seq Data Processing & Quantification Pipeline](https://github.com/epigen/atacseq_pipeline) for processing, quantification and annotation of ATAC-seq samples.
- [UCSC Genome Browser annotation track database](https://genome.ucsc.edu/cgi-bin/hgTables)
    - recommended source for the required 12 column BED file annotation of the respective genome.
    - e.g., for mm10 from UCSC as gzip https://genome.ucsc.edu/cgi-bin/hgTables assembly:mm10 -> track:NCBI RefSeq -> table:refFlat; output format: BED

# Publications
The following publications successfully used this module for their analyses.
- ...
