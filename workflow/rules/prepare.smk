# merge .bam files of the same group according to annotation using samtools
rule merge_bams:
    input:
        get_bams,
    output:
        merged_bam = os.path.join(result_path, 'merged_bams','{group}.bam'),
    resources:
        mem_mb = config.get("mem", "4000"),
    threads: 4*config.get("threads", 1)
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
        samtools merge -@ {threads} {output.merged_bam} {input}
        
        samtools index -@ {threads} -b {output.merged_bam}
        """
    
# generate a bigWig file per group using bamCoverage
rule coverage:
    input:
        bam = os.path.join(result_path, 'merged_bams','{group}.bam'),
    output:
        bigWig = os.path.join(result_path, 'bigWigs','{group}.bw'),
    resources:
        mem_mb=config.get("mem", "4000"),
    threads: 4*config.get("threads", 1)
    conda:
        "../envs/pygenometracks.yaml",
    log:
        os.path.join("logs","rules","coverage_{group}.log"),
    params:
        # bamCoverage parameters
        bamCoverage_parameters = config["bamCoverage_parameters"],
        # cluster parameters
        partition = config.get("partition"),
    shell:
        """
        bamCoverage --bam {input.bam} \
            {params.bamCoverage_parameters} \
            -o {output.bigWig} > {output.bigWig}.log 2>&1;
        """