##### utility functions #####

def get_sc_bam(wildcards):
#     sample, group = wildcards.sample_group.split('__')
#     return [sc_file_group_dict[sample]["bam"], sc_file_group_dict[sample]["metadata"]]
    return [sc_file_group_dict[wildcards.sample]["bam"], sc_file_group_dict[wildcards.sample]["metadata"]]

# def aggregate_group_bams(wildcards):
#     checkpoint_output = checkpoints.split_sc_bam.get(**wildcards).output[0]
#     # Assuming the output directory contains the split BAM files named '{group}.bam'
#     bam_files = glob.glob(os.path.join(checkpoint_output, f"{wildcards.group}.bam"))
#     print(bam_files)
#     return bam_files
    

def get_bams(wildcards):
    if wildcards.group in sc_groups:
        samples = [key for key, value in sc_file_group_dict.items() if wildcards.group in value['groups']]
        return expand(os.path.join(result_path, 'sc_bams', "{sample}", "{}.bam".format(wildcards.group)), sample=samples)
#         return [os.path.join(result_path, 'sc_bams', sample, "{}.bam".format(wildcards.group)) for sample in samples]
        
#         # find all samples with wildcards.group in ["groups"] and generate BAM file list sc_bams/{key}__{group}
#         checkpoint_output = checkpoints.split_sc_bam.get(**wildcards).output[0]
#         # Assuming the output directory contains the split BAM files named '{group}.bam'
#         bam_files = glob.glob(os.path.join(checkpoint_output, wildcards.sample, f"{wildcards.group}.bam"))
#         print(bam_files)
#         return bam_files
    else:
        return annot.loc[annot['group']==wildcards.group,'bam'].to_list()

# def get_merged_bam(wildcards):
#     if wildcards.group in sc_groups:
#         return os.path.join(result_path, 'sc_bams','{group}.bam')
#     else:
#         return os.path.join(result_path, 'merged_bams','{group}.bam')

def get_bigWigs(wildcards):
#     if wildcards.category=='ALL':
    return expand(os.path.join(result_path, 'bigWigs','{group}.bw'),group=plot_groups)
#     else:
#         return expand(os.path.join(result_path, 'bigWigs','{group}.bw'),group=sorted(annot.loc[annot['category']==wildcards.category,'group'].unique()))

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
#     if wildcards.category=='ALL':
    # extract group information and order from inputs
    groups = [os.path.basename(path).replace('.bw', '') for path in input]
    # map groups to categories
#     annot_reduced = annot[['group', 'category']].drop_duplicates()
#     categories = annot_reduced[annot_reduced['group'].isin(groups)]['category'].tolist()
    # map categories to colors & make color parameter string
#     colors = [config["track_colors"][group] for group in groups if group in config["track_colors"] else "#000000"]
    colors = [config["track_colors"][group] if group in config["track_colors"] else "#000000" for group in groups]
    colors_str = "' '".join(colors)
    colors_str = "'"+colors_str+"'"
    return colors_str
#     else:
#         return ("'"+config["track_colors"][wildcards.category]+"' ") * len(input)

