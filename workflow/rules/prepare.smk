# merge .bam files of the same group according to annotation using samtools
rule merge_bams:
    input:
        get_bams,
    output:
        merged_bam = os.path.join(config["result_path"], module_name, 'merged_bams','{group}.bam'),
    resources:
        mem_mb = config.get("mem", "4000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/pygenometracks.yaml",
    log:
        os.path.join("logs","rules","merge_bams_{group}.log"),
    params:
#         bams = lambda w: annot.loc[annot['group']=="{}".format(w.group),'bam'].to_list(),
        # cluster parameters
        partition = config.get("partition"),
    shell:
        """
        samtools merge {output.merged_bam} {input}
        
        samtools index -b {output.merged_bam}
        """
    
# generate a bigWig file per group using bamCoverage
rule make_bigwigs:
    input:
        merged_bam = os.path.join(config["result_path"], module_name, 'merged_bams','{group}.bam'),
    output:
        bigwig = os.path.join(config["result_path"], module_name, 'bigwigs','{group}.bw'),
    resources:
        mem_mb=config.get("mem", "4000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/pygenometracks.yaml",
    log:
        os.path.join("logs","rules","make_bigwigs_{group}.log"),
    params:
        # bamCoverage parameters
        extendReads =  lambda w: "--extendReads 175" if "ATAC" in "{}".format(w.group) else " ",
        genome_size = genome_size,
        # cluster parameters
        partition = config.get("partition"),
    shell:
        """
        bamCoverage --bam {input.merged_bam} \
            -p max --binSize 10  --normalizeUsing RPGC \
            --effectiveGenomeSize {params.genome_size} {params.extendReads} \
            -o "{output.bigwig}" > "{output.bigwig}.log" 2>&1;
        """