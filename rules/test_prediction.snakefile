#Test prediction of breeding values
rule mixed_solve:
    input:
        test_info = config["ind_info"],
        test_gt = config["gt"],
        test_fix = config["fix"],
        effect_sizes = "data/output/effect_sizes.txt" 
    output:
        out = "data/output/effect_sizes.txt"     
    conda:
        "env_configs/R.yaml"
    shell:
        '''
        Rscript scripts/test_prediction.R 
        '''