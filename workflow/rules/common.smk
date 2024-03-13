##### utility functions #####
def get_bams(wildcards):
    return annot.loc[annot['group']==wildcards.group,'bam'].to_list()

def get_bigwigs(wildcards):
    if wildcards.category=='ALL':
        return expand(os.path.join(config["result_path"], module_name, 'bigwigs','{group}.bw'),group=sorted(annot['group'].unique()))
    else:
        return expand(os.path.join(config["result_path"], module_name, 'bigwigs','{group}.bw'),group=sorted(annot.loc[annot['category']==wildcards.category,'group'].unique()))

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
    return chrom, int(start)-base_buffer, int(end)+base_buffer, count, gene_dict[gene]["ymax"]

def parse_region(region):
    chrom, start, end = region.replace('-', ':').split(':')
    return chrom, int(start), int(end), 1