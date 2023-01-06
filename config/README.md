# Configuration

You need one configuration file and one annotation file to run the complete workflow. You can use the provided example as starting point. Always use absolute paths. If in doubt read the comments in the config and/or try the default values.

- project configuration (config/config.yaml): different for every project/dataset and configures the analyses to be performed
- sample annotation (sample_annotation): CSV file consisting of three mandatory columns
    - bam: absolute path to mapped/aligned (and filtered) BAM file (*.bam)
    - group: a string describing a group of samples to be merged (e.g., RNA_untreated), which are visualized together as one genome track. Note: multiple rows can have the same group and the respective BAM files will be merged using samtools.
    - category: a string describing a category of groups (e.g., RNA), which tracks are stacked vertically together in one plot.
