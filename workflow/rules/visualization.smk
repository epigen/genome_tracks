
# prepare a UCSC genome browser track hub
rule ucsc_hub:
    input:
        bigWig_files = expand(os.path.join(result_path, 'bigWigs','{group}.bw'), group=plot_groups),
    output:
        bigWig_symlinks = expand(os.path.join(result_path, "bigWigs", config["genome"], "{group}.bw"), group=plot_groups),
        genomes_file = os.path.join(result_path, "bigWigs", "genomes.txt"),
        hub_file = os.path.join(result_path, "bigWigs", "hub.txt"),
        trackdb_file = os.path.join(result_path, "bigWigs", config["genome"], "trackDb.txt"),
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
            
            track_db = ['track {}'.format(config["project_name"]),
                        'type bigWig', 'compositeTrack on', 'autoScale on', 'maxHeightPixels 32:32:8',
                        'shortLabel {}'.format(config["project_name"][:8]),
                        'longLabel {}'.format(config["project_name"]),
                        'visibility full',
                        '', '']
            for group in plot_groups:
                
                hex_color = config["track_colors"][group] if group in config["track_colors"] else "#000000"
                track_color = tuple(int(hex_color[i:i+2], 16) for i in (1, 3, 5)) # convert to RGB
                
                track = ['track {}'.format(group),
                         'shortLabel {}'.format(group),
                         'longLabel {}'.format(group),
                         'bigDataUrl {}.bw'.format(group),
                         'parent {} on'.format(config["project_name"]),
                         'type bigWig', 'windowingFunction mean',
                         'color {}'.format(','.join(map(str, track_color))),
                         '', '']
                
                track_db += track

            tf.write('\n'.join(track_db))


# plot genome tracks of groups using pyGenomeTracks wrapper gtracks
rule plot_tracks:
    input:
        bigWigs = expand(os.path.join(result_path, 'bigWigs','{group}.bw'), group=plot_groups),
#         get_bigWigs,
    output:
        genome_track = report(os.path.join(result_path, 'tracks','{gene}.'+config["file_type"]), 
                              caption="../report/genome_tracks.rst", 
                              category="{}_{}".format(config["project_name"], module_name),
                              subcategory="genome tracks",
                              labels={
                                  "data": "{gene}",
                                  "type": "pyGenomeTrack",
                                  "misc": "ymax {}".format(lambda w: gene_annot_df.loc["{}".format(w.gene),"ymax"]),
                              }),
    resources:
        mem_mb=config.get("mem", "4000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/pygenometracks.yaml",
    log:
        os.path.join("logs","rules","plot_tracks_{gene}.log"),
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
    shell:
        """
        export GTRACKS_GENES_PATH={params.genome_bed}
        
        gtracks {params.coordinates} \
            {input.bigWigs} \
            {output.genome_track} \
            --genes {params.genome_bed} \
            {params.ymax} \
            --gene-rows {params.gene_rows} \
            --genes-height {params.gene_rows} \
            --x-axis {params.xaxis} \
            --width {params.width} \
            --color-palette {params.colors}
        """

# create interactive IGV-report
# https://github.com/igvteam/igv-reports
rule igv_report:
    input:
        bed = os.path.join(result_path,'genes.bed'),
        tracks = expand(os.path.join(result_path, 'merged_bams','{group}.bam'), group=plot_groups),
    output:
        igv_report = report(os.path.join(result_path, "igv-report.html"), 
                              caption="../report/igv_report.rst", 
                              category="{}_{}".format(config["project_name"], module_name),
                              subcategory="IGV report",
                              labels={
                                  "data": "IGV report",
                                  "type": "HTML",
                                  "misc": "interactive",
                              }),
    resources:
        mem_mb = config.get("mem", "4000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/igv_reports.yaml",
    log:
        os.path.join("logs","rules","igv_report.log"),
    params:
        # igv-reports parameters
        genome = config["genome"],
    shell:
        """
        create_report {input.bed} \
            --genome {params.genome} \
            --tracks {input.tracks} \
            --output {output.igv_report}
            
            
        # replace 'Variants' with 'Genes and genomic regions of interest' in the HTML
        sed 's/<label for="collapsible" class="lbl-toggle">Variants<\/label>/<label for="collapsible" class="lbl-toggle">Genes and genomic regions<\/label>/g' {output.igv_report} > {output.igv_report}.tmp && mv {output.igv_report}.tmp {output.igv_report}
        """