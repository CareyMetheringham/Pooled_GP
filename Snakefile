configfile: "config.yaml"

include: 'rules/simulation.snakefile',
include: 'rules/read_data.snakefile'
include: 'rules/run_gp.snakefile'

#rule sim_object:
#    input:
#        sim_training_obj = "data/training_sim.RData" 

rule run_gp:
    input:
        RData = "data/pools.RData",
        effect_sizes = "data/output/effect_sizes.txt"
