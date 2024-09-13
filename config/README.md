# Configuration

You need one configuration file and one annotation file to run the complete workflow.  If in doubt read the comments in the config and/or try the default values.

- project configuration (`config/config.yaml`): different for every project/dataset and configures the analyses to be performed.
- sample annotation (sample_annotation): CSV file consisting of two mandatory columns
    - bam: absolute path to mapped/aligned (and filtered) BAM file (*.bam).
    - group: 
      - bulk samples: a string describing a group of samples to be merged (e.g., RNA_untreated), which are visualized together as one genome track. Note: multiple rows can have the same group and the respective BAM files will be [merged using samtools](https://www.htslib.org/doc/samtools-merge.html).
      - single-cell samples: path to a tab-separated metadata table wihtout header (TSV) where the first column are cell barcodes (CB) and the second a group variable (e.g., untreated/treated or cluster_1). Note: each single-cell BAM file is split into multiple BAM files according to the metadata file using [sinto](https://timoast.github.io/sinto/basic_usage.html#filter-cell-barcodes-from-bam-file), and subsequenytly merged by group across all samples using samtools and downstream processed and visualized the same as bulk samples. It is important that the single-cell metadata.tsv is sample specific, because simply by chance different samples could have the same cell barcodes (CB).

Set workflow-specific `resources` or command line arguments (CLI) in the workflow profile `workflow/profiles/default.config.yaml`, which supersedes global Snakemake profiles.