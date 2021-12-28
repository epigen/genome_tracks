# Genome Browser Tracks Pipeline
Workflow for easy genome browser track plotting using the packages [gtracks](https://github.com/anthony-aylward/gtracks) and [pyGenomeTracks](https://github.com/deeptools/pyGenomeTracks).

From aligned BAM files to genome browser tracks.

![rulegraph](https://github.com/sreichl/genome_tracks/blob/main/workflow/rulegraph.svg)

# Input
- project config
  - paths to 
    - sample annotations
    - bed file
    - results folder (merged_bams, bigwigs, tracks)
    - gene_list
  - gtracks parameters
  - cluster parameters
- annotation file
  - BAM file paths
  - group (of samples to be merged)
  - category (groups to be in one plot)
- 12 column BED file annotation 
  - eg for mm10 from UCSC as gzip https://genome.ucsc.edu/cgi-bin/hgTables assembly:mm10 -> track:NCBI RefSeq -> table:refFlat; output format: BED
- BAM files
- list of genes to plot OR genomic regions (here gene-rows=1 is hardcoded) as text file

# Executed Steps
1. merge and index BAM files with samtools by group (samtools)
2. generate bigwigs per merged bam file (bamCoverage)
4. get information per gene from BED file
    - coordinates from BED file & extend left & right by parameter base_buffer
    - number of isoforms ie number of lines in the BED file (only for genes, for regions it is hardcoded to 1)
3. make one plot per category of bigwigs and gene with the before determined gene-parameters (gene-rows, gene-heights, coordinates)

# Parameters (for project config file)
- genome size
- gtracks
    - y max value (number or '' for auto)
    - x axis position
    - width of plot in cm (default: 40)
- base buffer
- cluster parameters
