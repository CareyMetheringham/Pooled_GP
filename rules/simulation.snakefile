#Simulate GP data
rule simulate:
    param:
        n_pop = config["sim_pops"]
        n_ind = config["sim_ind"]
        n_site = config["sim_sites"]
        h2 = config["sim_h2"]
        maf = config["sim_MAF"]
        threshold = config["sim_threshold"]
    output:
        out1 = "data/training_sim.RData"
        out2 = "data/test_sim.RData"
    conda:
        "env_configs/R.yaml"
    shell:
        '''
        Rscript scripts/simulate.R --p {param.n_pop} --i {param.n_ind} \
        --s {param.n_site} --h {param.h2} --m {param.maf} --t {param.threshold} \
        --o1 {output.out1} --o2 {output.out2}
        '''