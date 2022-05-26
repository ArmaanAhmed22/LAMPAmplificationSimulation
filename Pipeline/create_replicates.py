import os
file = snakemake.wildcards[0]
replicate_number = snakemake.wildcards[1]

path = ''.join(snakemake.input[0].split("/")[:-1])
os.system(f"mkdir -r {path}")

os.system(f"cp {snakemake.input[0]} {snakemake.output[0]}")
os.system(f"cp {snakemake.input[1]} {snakemake.output[1]}")
