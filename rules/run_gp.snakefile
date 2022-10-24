#Read in pooled data
rule mixed_solve:
    input:
        data = 
    output:
        out = "data/output/effect_sizes.txt"
    params:
        dir = "data/reference"
    log:
        "logs/mixed_solve.log"
    conda:
        "env_configs/R.yaml"
    shell:
        '''
        Rscript scripts/predict.R --data {input.data}
        '''
        
#Run mixed solve on read counts
rule mixed_solve:
    input:
        data = 
    output:
        out = "data/output/effect_sizes.txt"
    params:
        dir = "data/reference"
    log:
        "logs/mixed_solve.log"
    conda:
        "env_configs/R.yaml"
    shell:
        '''
        Rscript scripts/predict.R --data {input.data}
        '''