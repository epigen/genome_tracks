# rule for used conda environments to document the exact versions and builds of the used software        
rule env_export:
    output:
        report(os.path.join(config["result_path"],'envs', module_name,'{env}.yaml'),
                      caption="../report/software.rst", 
                      category="Software", 
                      subcategory="{}_{}".format(config["project_name"], module_name)
                     ),
    conda:
        "../envs/{env}.yaml"
    resources:
        mem_mb = config.get("mem", "4000"),
    threads: config.get("threads", 1)
    log:
        os.path.join("logs","rules","env_{env}.log"),
    params:
        partition = config.get("partition"),
    shell:
        """
        conda env export > {output}
        """
        
# add configuration files to report        
rule config_export:
    output:
        configs = report(os.path.join(config["result_path"],'configs', module_name,'{}_config.yaml'.format(config["project_name"])), 
                         caption="../report/configs.rst", 
                         category="Configuration", 
                         subcategory="{}_{}".format(config["project_name"], module_name)
                        )
    resources:
        mem_mb = 1000, #config.get("mem", "1000"),
    threads: config.get("threads", 1)
    log:
        os.path.join("logs","rules","config_export.log"),
    params:
        partition=config.get("partition"),
    run:
        with open(output["configs"], 'w') as outfile:
            json.dump(config, outfile)
        
# export used annotation file for documentation and reproducibility
rule annot_export:
    input:
        config["sample_annotation"],
    output:
        annot = report(os.path.join(config["result_path"],'configs', module_name,'{}_annot.csv'.format(config["project_name"])), 
                         caption="../report/configs.rst", 
                         category="Configuration", 
                         subcategory="{}_{}".format(config["project_name"], module_name)
                        )
    resources:
        mem_mb = 1000, #config.get("mem_small", "16000"),
    threads: config.get("threads", 1)
    log:
        os.path.join("logs","rules","annot_export.log"),
    params:
        partition=config.get("partition"),
    shell:
        """
        cp {input} {output}
        """

# export used gene list for documentation and reproducibility
rule gene_list_export:
    input:
        config["gene_list"],
    output:
        gene_list = report(os.path.join(config["result_path"],'configs', module_name,'gene_list.txt'), 
                            caption="../report/gene_list.rst", 
                            category="Configuration", 
                            subcategory="{}_{}".format(config["project_name"], module_name)
                           ),
    resources:
        mem_mb = 1000, #config.get("mem_small", "16000"),config.get("mem", "16000"),
    threads: config.get("threads", 1)
    log:
        os.path.join("logs","rules","gene_list_export.log"),
    params:
        partition=config.get("partition"),
    shell:
        """
        cp {input} {output}
        """

        

        