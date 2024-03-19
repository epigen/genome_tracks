##### utility functions #####

def get_sc_bam(wildcards):
    return [sc_file_group_dict[wildcards.sample]["bam"], sc_file_group_dict[wildcards.sample]["metadata"]]    

def get_bams(wildcards):
    if wildcards.group in sc_groups:
        samples = [key for key, value in sc_file_group_dict.items() if wildcards.group in value['groups']]
        return expand(os.path.join(result_path, 'sc_bams', "{sample}", "{}.bam".format(wildcards.group)), sample=samples)
    else:
        return annot.loc[annot['group']==wildcards.group,'bam'].to_list()

# def get_bigWigs(wildcards):
#     return expand(os.path.join(result_path, 'bigWigs','{group}.bw'),group=sort(plot_groups))

def parse_gene(gene):
    count = 0
    
    with gzip.open(config['genome_bed'], 'rt') as f:
        for line in f:
            parsed_line = line.split()
            if parsed_line[3] == gene:
                count = count+1
                if count==1:
                    chrom, start, end = parsed_line[:3]
                else:
                    tmp_chrom, tmp_start, tmp_end = parsed_line[:3]
                    if int(tmp_start)<int(start):
                        start = tmp_start
                    if int(tmp_end)>int(end):
                        end = tmp_end
        if count==0:
            return -1
    return chrom, int(start)-config['base_buffer'], int(end)+config['base_buffer'], count

def parse_region(region):
    chrom, start, end = region.replace('-', ':').split(':')
    return chrom, int(start), int(end), 1

def get_colors(wildcards, input):
    # extract group information and order from inputs
    groups = [os.path.basename(path).replace('.bw', '') for path in input]

    colors = [config["track_colors"][group] if group in config["track_colors"] else "#000000" for group in groups]
    colors_str = "' '".join(colors)
    colors_str = "'"+colors_str+"'"
    
    return colors_str

