#Run mixed solve on read counts
rule mixed_solve:
    input:
        data = "data/pools.RData"
    output:
        out = "data/output/effect_sizes.txt"     
    log:
        "logs/mixed_solve.log"
    params:
        reps = config["rrblup_reps"]
    conda:
        "env_configs/R.yaml"
    shell:
        '''
        Rscript scripts/predict.R --data {input.data} --reps {params.reps} --out {output.out}
        '''