rule create_replicates:
    input:
        "Input/{file}/primers.json",
        "Input/{file}/conditions.json"
    output:
        "Input/{file}/replicates/R{number}/primers.json",
        "Input/{file}/replicates/R{number}/conditions.json"
    script:
        "Pipeline/create_replicates.py"

rule simulate_amplification:
    input:
        "Input/{file}/primers.json",
        "Input/{file}/conditions.json",
    output:
        "Output/{file}/out.png",
        "Output/{file}/out.csv"
    script:
        "Pipeline/simulate_amplification.py"

