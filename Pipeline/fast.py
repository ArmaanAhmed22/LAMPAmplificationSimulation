from numpy import number
from reaction_mechanizer.pathway.reaction import ReactionMechanism, SimpleStep
from reaction_mechanizer.drawing.mechanism_reaction_visualizer import ReactionVisualizer
import seaborn as sns
import matplotlib.pyplot as plt
import json
import pandas as pd
STEPS = 5000
"""
FROM: 10.1002/bit.20555


P: Sense primer
S: Single-stranded template sense
P': Anti-sense primer
S': Single-stranded template antisense
SP: Single-stranded template-primer complex
ESP: Single stranded-template-primer-enzyme complex
N: Nucleotide
D_i: incomplete double-stranded DNA with i nucleotides added to the SP complex
D: DNA target
"""

with open(snakemake.input[0], "r") as h1, open(snakemake.input[1], "r") as h2:
    primers = json.load(h1)
    conditions = json.load(h2)


def extend(label:int, mechanism_list:list[SimpleStep], conditions:dict, annotation:str=""):
    S="S"+annotation
    P="P"+annotation
    D="D"+annotation
    if label == 1:
        mechanism_list.append(SimpleStep.str_to_step(f"E{S}{P}+N->E{S}{P}N"))
        mechanism_list[-1].set_rate_constant(kf=10**10)
        mechanism_list.append(SimpleStep.str_to_step(f"E{S}{P}N->E{D}1"))
        mechanism_list[-1].set_rate_constant(kf=conditions["enzyme"]["k_nucleotide_addition"])
    else:
        mechanism_list.append(SimpleStep.str_to_step(f"E{D}{label-1}+N->E{D}{label-1}N"))
        mechanism_list[-1].set_rate_constant(kf=10**10)
        mechanism_list.append(SimpleStep.str_to_step(f"E{D}{label-1}N->E{D}{label}"))
        mechanism_list[-1].set_rate_constant(kf=conditions["enzyme"]["k_nucleotide_addition"])

def extend_quick(number_bases, mechanism_list:list[SimpleStep], conditions:dict, annotation:str=""):
    S="S"+annotation
    P="P"+annotation
    D="D"+annotation
    mechanism_list.append(SimpleStep.str_to_step(f"E{S}{P}+{number_bases}N->E{D}{number_bases}"))
    mechanism_list[-1].set_rate_constant(kf=conditions["enzyme"]["k_nucleotide_addition"]/number_bases)


def mechanism_basis(conditions: dict) -> list[SimpleStep]:
    mechanism_list:list[SimpleStep] = []
    #Creation of primer to target complex
    k_annealing_primer = conditions["primer"]["k_annealing"]
    k_dissociation_primer = conditions["primer"]["k_dissociation"]
    mechanism_list.append(SimpleStep.str_to_step("D+P->SP+S'")) #Sense primer forming dna-primer complex
    mechanism_list[-1].set_rate_constant(kf=k_annealing_primer, kr=k_dissociation_primer)
    mechanism_list.append(SimpleStep.str_to_step("D+P'->S'P'+S")) #Antisense primer forming dna-primer complex
    mechanism_list[-1].set_rate_constant(kf=k_annealing_primer, kr=k_dissociation_primer)
    mechanism_list.append(SimpleStep.str_to_step("S+P->SP"))
    mechanism_list[-1].set_rate_constant(kf=k_annealing_primer, kr=k_dissociation_primer)
    mechanism_list.append(SimpleStep.str_to_step("S'+P'->S'P'"))
    mechanism_list[-1].set_rate_constant(kf=k_annealing_primer, kr=k_dissociation_primer)
    #------------------------------------------------------------
    #Enzyme joining
    k_enzyme_SP_complex_formation = conditions["enzyme"]["k_enzyme_SP_complex_formation"]
    k_enzyme_SP_complex_dissociation = conditions["enzyme"]["k_enzyme_SP_complex_dissociation"]
    mechanism_list.append(SimpleStep.str_to_step("E+SP->ESP")) #Enzyme joins in
    mechanism_list[-1].set_rate_constant(kf=k_enzyme_SP_complex_formation, kr=k_enzyme_SP_complex_dissociation)
    mechanism_list.append(SimpleStep.str_to_step("E+S'P'->ES'P'")) #Enzyme joins in
    mechanism_list[-1].set_rate_constant(kf=k_enzyme_SP_complex_formation, kr=k_enzyme_SP_complex_dissociation)
    #------------------------------------------------------------------------
    
    return mechanism_list

def get_fill_size(cur_primer_set) -> tuple[int, int]:
    need_to_fill_size_forward = cur_primer_set["PB"]["position"] + cur_primer_set["PB"]["length"] - (cur_primer_set["PF"]["position"] + cur_primer_set["PF"]["length"])
    """
                       PRIMER_B
                       ||||||||
    GGGGGGGGGGGGGGGGGGGGGGGGGGG
    |||||||||||||||||||||||||||
    PRIMER_F*******************
    """
    need_to_fill_size_backward = cur_primer_set["PB"]["position"] - cur_primer_set["PF"]["position"]
    
    return (need_to_fill_size_forward, need_to_fill_size_backward)
    

def primer_specific(mechanism_list, cur_primer_set, conditions) -> list[SimpleStep]:
    need_to_fill_size_forward, need_to_fill_size_backward = get_fill_size(cur_primer_set)
    extend_quick(need_to_fill_size_forward, mechanism_list, conditions)
    mechanism_list.append(SimpleStep.str_to_step(f"ED{need_to_fill_size_forward}->E+D"))
    mechanism_list[-1].set_rate_constant(kf=conditions["enzyme"]["k_enzyme_SP_complex_dissociation"])

    extend_quick(need_to_fill_size_forward, mechanism_list, conditions, annotation="'")
    mechanism_list.append(SimpleStep.str_to_step(f"ED'{need_to_fill_size_backward}->E+D"))
    mechanism_list[-1].set_rate_constant(kf=conditions["enzyme"]["k_enzyme_SP_complex_dissociation"])
    #------------------------------------------------------------------------------
    return mechanism_list

def generate_mechanism(cur_primer_set, conditions):
    cur_mechanism_list = mechanism_basis(conditions)
    mechanism_list = primer_specific(cur_mechanism_list, cur_primer_set, conditions)
    
    mechanism = ReactionMechanism(mechanism_list)
    all_species = {species for step in mechanism_list for chosen in (step.reactants, step.products) for species in chosen.keys()}
    #print(all_species)
    vis = ReactionVisualizer(mechanism)
    params = {
        "initial_state": 
            {
                "P": conditions["primer"]["concentration"],
                "P'": conditions["primer"]["concentration"],
                "N": conditions["nucleotide"]["concentration"],
                "E": conditions["enzyme"]["concentration"],
                "D": conditions["DNA"]["concentration"]
            },
        "time_end": 3600,
        "number_steps":STEPS
    }
    params["initial_state"] = {species:params["initial_state"].get(species, 0) for species in all_species}
    #print(params["initial_state"])
    df = vis.progress_reaction(**params)
    df.to_csv("real_out.csv")
    output = df[["Time", "D"]]
    D_concentration = output["D"]
    print(all_species)
    need_to_fill_size_forward, need_to_fill_size_backward = get_fill_size(cur_primer_set)
    D_concentration+=df[f"ED{need_to_fill_size_forward}"]
    D_concentration+=df[f"ED'{need_to_fill_size_backward}"]

    D_concentration += 1/2*df["S'"] + 1/2*df["S"] + 1/2*df["SP"] + 1/2*df["S'P'"] + 1/2*df["ESP"] + 1/2*df["ES'P'"]
    output["D"] = D_concentration
    return output

def main():
    output = pd.DataFrame({"Time":[], "ID":[], "DNA":[]})
    for i, cur_primer_set in enumerate(primers["primer_sets"]):
        cur_output = generate_mechanism(cur_primer_set, conditions)
        cur_df = pd.DataFrame({"Time":[], "ID":[], "DNA":[]})
        cur_df["Time"] = cur_output["Time"]
        cur_df["ID"] = [i]*STEPS
        cur_df["DNA"] = cur_output["D"]
        
        output = pd.concat([output, cur_df], ignore_index=True)
        
    sns.lineplot(x="Time", y="DNA", data=output, hue="ID")
    output.to_csv(snakemake.output[1])
    plt.savefig(snakemake.output[0], dpi=600, bbox_inches="tight")

main()