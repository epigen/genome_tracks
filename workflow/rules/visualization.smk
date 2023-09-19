
# plot genome tracks of groups in the same category together using pygenometracks wrapper gtracks
rule plot_tracks:
    input:
        get_bigwigs,
    output:
        genome_track = report(os.path.join(config["result_path"], module_name, 'tracks','{category}_{gene}.'+config["file_type"]), 
                              caption="../report/genome_tracks.rst", 
                              category="{}_{}".format(config["project_name"], module_name),
                              subcategory="{category}",
                              labels={
                                  "data": "{gene}",
                                  "type": "genome track",
                                  "misc": "ymax {}".format(ymax),
                              }),
    resources:
        mem_mb=config.get("mem", "4000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/pygenometracks.yaml",
    log:
        os.path.join("logs","rules","plot_tracks_{category}_{gene}.log"),
    params:
        # gtracks parameters
        gene = lambda w: "{}".format(w.gene),
        genome_bed = config['genome_bed'],
        ymax = lambda w: "--max {}".format(ymax) if ymax!='' else " ",
        xaxis = xaxis,
        coordinates = lambda w: "{}:{}-{}".format(gene_annot_df.loc[w.gene,'chr'], gene_annot_df.loc[w.gene,'start'], gene_annot_df.loc[w.gene,'end']),
        # eg chr14:103052047-103053094
        gene_rows = lambda w: "{}".format(gene_annot_df.loc[w.gene,'count'].astype(np.int64)),
        width = "{}".format(width),
        # cluster parameters
        partition = config.get("partition"),
    shell:
        """
        export GTRACKS_GENES_PATH={params.genome_bed}
        
        gtracks {params.coordinates} \
            {input} \
            {output.genome_track} \
            --genes {params.genome_bed} {params.ymax} \
            --gene-rows {params.gene_rows} \
            --genes-height {params.gene_rows} \
            --x-axis {params.xaxis} \
            --width {params.width}
        """