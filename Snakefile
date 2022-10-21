configfile: "config.yaml"

include: 'rules/run_model.snakefile'
include: 'rules/simulation.snakefile'

rule sim_object:
    input:
        sim_training_obj = "data/training_sim.RData"
