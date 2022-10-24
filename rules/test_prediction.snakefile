#Test prediction of breeding values
rule get_gebv:
    input:
        test_info = config["ind_info"],
        test_gt = config["gt"],
        effect_sizes = "data/output/effect_sizes.txt" 
    output:
        gebv = "data/output/gebv.txt"     
    conda:
        "env_configs/R.yaml"
    shell:
        '''
        Rscript scripts/test_prediction.R --es {input.effect_sizes} --info {input.test_info} --gt {input.test_gt} --out {output.gebv}
        '''