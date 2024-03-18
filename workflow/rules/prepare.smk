# create BED file of genes/genomic regions
rule make_bed:
    output:
        bed_file = os.path.join(result_path,'genes.bed'),
    resources:
        mem_mb = config.get("mem", "1000"),
    threads: config.get("threads", 1)
    log:
        os.path.join("logs","rules","make_bed.log"),
    params:
        # cluster parameters
        partition = config.get("partition"),
    run:
        # Convert DataFrame to .BED format (chr, start, end, name)
        bed_df = gene_annot_df.reset_index().rename(columns={'index': 'name'})
        bed_df = bed_df[['chr', 'start', 'end', 'name']]
        bed_df.to_csv(output.bed_file, sep='\t', header=False, index=False)

# in case of single cell BAM files, create split BAM files by metadata using sinto
# https://timoast.github.io/sinto/basic_usage.html#filter-cell-barcodes-from-bam-file
rule split_sc_bam:
    input:
        get_sc_bam,
    output:
#         directory(os.path.join(result_path, 'sc_bams','{sample}')),
        sc_bams = expand(os.path.join(result_path, 'sc_bams','{{sample}}','{group}.bam'), group = sc_groups),
#         sc_bams = expand(os.path.join(result_path, 'sc_bams','{{sample}}','{group}.bam'), group = lambda w:sc_file_group_dict["{}".format(w.sample)]["groups"]),
#         sc_bams = os.path.join(result_path, 'sc_bams','{sample}__{group}.bam'),
    resources:
        mem_mb = config.get("mem", "4000"),
    threads: 4*config.get("threads", 1)
    conda:
        "../envs/sinto.yaml",
    log:
        os.path.join("logs","rules","split_sc_bam_{sample}.log"),
    params:
        # cluster parameters
        partition = config.get("partition"),
    shell:
        """
        sinto filterbarcodes -b {input[0]} -c {input[1]} --outdir $(dirname {output[0]}) -p {threads}
        
        # generate empty bam file if group is not present in sample
        for file in {output}; do
          if [ ! -f "$file" ]; then
            touch "$file"
          fi
        done
        """

# # merge .bam files of the same group according to annotation using samtools
# rule merge_sc_bams:
#     input:
#         aggregate_group_bams,
#     output:
#         merged_bam = os.path.join(result_path, 'merged_bams','{group}.bam'),
#     resources:
#         mem_mb = config.get("mem", "4000"),
#     threads: 4*config.get("threads", 1)
#     conda:
#         "../envs/pygenometracks.yaml",
#     log:
#         os.path.join("logs","rules","merge_sc_bams_{group}.log"),
#     params:
# #         bams = lambda w: annot.loc[annot['group']=="{}".format(w.group),'bam'].to_list(),
#         # cluster parameters
#         partition = config.get("partition"),
#     shell:
#         """
#         samtools merge -@ {threads} {output.merged_bam} {input}
        
#         samtools index -@ {threads} -b {output.merged_bam}
#         """

# merge .bam files of the same group according to annotation using samtools
# https://www.htslib.org/doc/samtools-merge.html
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
# https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html
rule coverage:
    input:
#         get_prepared_bam,
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