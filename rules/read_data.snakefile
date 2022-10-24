#Read in pooled data
rule read_input:
    input:
        pool_rc = config["pool_rc"],
        pool_info = config["pool_info"],
        snp_list = config["snp_list"]
    output:
        RData = "data/pools.RData"
    params:
        snp_num = config["snp_num"]
    conda:
        "env_configs/R.yaml"
    shell:
        '''
        Rscript scripts/load_pool_data.R --pool_rc {input.pool_rc} --info {input.pool_info} --snps {input.snp_list} --snp_num {params.snp_num} --out {output.RData}
        '''