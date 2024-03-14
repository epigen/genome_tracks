
# prepare a UCSC genome browser track hub
rule ucsc_hub:
    input:
        bigWig_files = expand(os.path.join(result_path, 'bigWigs','{group}.bw'), group=sorted(annot['group'].unique())),
    output:
        bigWig_symlinks = expand(os.path.join(result_path, "bigWigs", config["genome"], "{group}.bw"), group=sorted(annot['group'].unique())),
        genomes_file = os.path.join(result_path, "bigWigs", "genomes.txt"),
        hub_file = os.path.join(result_path, "bigWigs", "hub.txt"),
        trackdb_file = os.path.join(result_path, "bigWigs", config["genome"], "trackDb.txt"),
    params:
        # cluster parameters
        partition=config.get("partition"),
    resources:
        mem_mb=config.get("mem", "1000"),
    threads: config.get("threads", 1)
    log:
        "logs/rules/ucsc_hub.log"
    run:
        # create bigWig symlinks
        for i in range(len(input.bigWig_files)):
            os.symlink(os.path.join('../',os.path.basename(input.bigWig_files[i])), output.bigWig_symlinks[i])

        # create genomes.txt
        with open(output.genomes_file, 'w') as gf:
            genomes_text = f'genome {config["genome"]}\ntrackDb {config["genome"]}/trackDb.txt\n'
            gf.write(genomes_text)

        # create hub file
        with open(output.hub_file, 'w') as hf:
            hub_text = [f'hub {config["project_name"]}',
                        f'shortLabel {config["project_name"]}',
                        f'longLabel {config["project_name"]}',
                        'genomesFile genomes.txt',
                        f'email {config["email"]}\n',]
            hf.write('\n'.join(hub_text))

        # create trackdb file
        with open(output.trackdb_file, 'w') as tf:
            colors = ['166,206,227', '31,120,180', '51,160,44', '251,154,153', '227,26,28',
                              '253,191,111', '255,127,0', '202,178,214', '106,61,154', '177,89,40']
            
            track_db = ['track {}'.format(config["project_name"]),
                        'type bigWig', 'compositeTrack on', 'autoScale on', 'maxHeightPixels 32:32:8',
                        'shortLabel {}'.format(config["project_name"][:8]),
                        'longLabel {}'.format(config["project_name"]),
                        'visibility full',
                        '', '']
            for group in sorted(annot['group'].unique()):
                track_color = '255,40,0'
                
#                 if config["annot_columns"][0]!="":
                color_hash = hash(annot.loc[annot['group']==group,'category'][0]) #hash(samples[sample_name][config["annot_columns"][0]])
                track_color = colors[color_hash % len(colors)]
                
                track = ['track {}'.format(group),
                         'shortLabel {}'.format(group),
                         'longLabel {}'.format(group),
                         'bigDataUrl {}.bigWig'.format(group),
                         'parent {} on'.format(config["project_name"]),
                         'type bigWig', 'windowingFunction mean',
                         'color {}'.format(track_color),
                         '', '']
                
                track_db += track

            tf.write('\n'.join(track_db))


# plot genome tracks of groups in the same category together using pyGenomeTracks wrapper gtracks
rule plot_tracks:
    input:
        get_bigWigs,
    output:
        genome_track = report(os.path.join(result_path, 'tracks','{category}_{gene}.'+config["file_type"]), 
                              caption="../report/genome_tracks.rst", 
                              category="{}_{}".format(config["project_name"], module_name),
                              subcategory="{category}",
                              labels={
                                  "data": "{gene}",
                                  "type": "genome track",
                                  "misc": "ymax {}".format(lambda w: gene_annot_df.loc["{}".format(w.gene),"ymax"]),
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
        ymax = lambda w: "--max {}".format(gene_annot_df.loc["{}".format(w.gene), "ymax"]) if gene_annot_df.loc["{}".format(w.gene),"ymax"]!=0 else " ",
        xaxis = config['x_axis'],
        coordinates = lambda w: "{}:{}-{}".format(gene_annot_df.loc[w.gene,'chr'], gene_annot_df.loc[w.gene,'start'], gene_annot_df.loc[w.gene,'end']),
        gene_rows = lambda w: "{}".format(gene_annot_df.loc[w.gene,'count'].astype(np.int64)),
        width = "{}".format(config['width']),
        colors = get_colors,
        # cluster parameters
        partition = config.get("partition"),
    shell:
        """
        export GTRACKS_GENES_PATH={params.genome_bed}
        
        gtracks {params.coordinates} \
            {input} \
            {output.genome_track} \
            --genes {params.genome_bed} \
            {params.ymax} \
            --gene-rows {params.gene_rows} \
            --genes-height {params.gene_rows} \
            --x-axis {params.xaxis} \
            --width {params.width} \
            --color-palette {params.colors}
        """