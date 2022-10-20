configfile: "config.yaml"

include: 'rules/qc.snakefile'

rule deseq:
    input:
        tx2gene = 'data/reference/tx2gene.csv',
        gene_counts = "data/quantification/gene.counts.table",
        transcript_counts = "data/quantification/transcript.counts.table",
        diff_ex_sex = "data/expression/diff_ex_sex"
         
rule quant:
    input:
        expand("data/quantification/{sample}_{lane}/quant.sf", sample=config["samples"], lane=config["lanes"]),

rule qc:
    input: 
       report = "data/fastqc/report_quality_control.html"
       