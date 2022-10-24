#Simulate GP data
rule simulate:
    params:
        n_pop = config["sim_pops"],
        n_ind = config["sim_ind"],
        n_site = config["sim_sites"],
        h2 = config["sim_h2"],
        maf = config["sim_MAF"],
        threshold = config["sim_threshold"]
    output:
        out1 = "data/training_sim.RData",
        out2 = "data/test_sim.RData"
    conda:
        "env_configs/R.yaml"
    shell:
        '''
        Rscript scripts/simulate.R --n_pop {params.n_pop} --n_ind {params.n_ind} --n_site {params.n_site} --h2 {params.h2} --maf {params.maf} --threshold {params.threshold} --out1 {output.out1} --out2 {output.out2}
        '''