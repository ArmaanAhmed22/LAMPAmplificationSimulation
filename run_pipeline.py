import argparse
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def parse():
    parser = argparse.ArgumentParser(description="Processor for AmplificationSimulator Tool")
    parser.add_argument("-c", "--cores", metavar="C", type=int, nargs=1, default=[1], help="Number of CPUs to use")
    parser.add_argument("-r", "--replicates", metavar="R", type=int, nargs=1, default=[1], help="Number of replicates to run")
    args = parser.parse_args()
    return args.cores[0], args.replicates[0]

number_cores, number_replicates = parse()
literal_replicate_string = ' '.join([f"Output/Multiplex1/replicates/R"+str(i)+"/out.png" for i in range(1, number_replicates+1)])
os.system(f"pyenv exec snakemake {literal_replicate_string} -c{number_cores} -f")

def stitch(replicates):
    number_steps = None
    master_dfs = []
    id_dict = {}
    for i in range(1, replicates+1):
        cur_df = pd.read_csv(f"Output/Multiplex1/replicates/R{str(i)}/out.csv")
        if number_steps == None:
            number_steps = cur_df.shape[0]
        ids = cur_df["ID"].unique()
        if len(id_dict.items()) == 0:
            for id in ids:
                id_dict[id] = None
        
        for id in ids:
            prev = None
            is_good = True
            for row in cur_df[cur_df["ID"] == id].itertuples():
                if prev == None:
                    prev = row.DNA
                    
                growth = (row.DNA - prev) / prev
                
                if growth < -0.01:
                    is_good = False
                    print(f"{prev=}")
                    print(f"{row.DNA=}")
                    break
                prev = row.DNA
            if is_good:
                id_dict[id] = i-1
        master_dfs.append(cur_df)
    
    for key, val in id_dict.items():
        if val == None:
            raise Exception(f"No replicate was found with a good simulation for ID {key}")
    final_df = pd.DataFrame({"Time": [], "DNA": [], "ID": []})
    for key, val in id_dict.items():
        cur_df = master_dfs[val]
        cur_df = cur_df[cur_df["ID"] == key]
        final_df = pd.concat([final_df, cur_df])
    sns.lineplot(x="Time", y="DNA", hue="ID", data=final_df)
    plt.savefig("Output/Multiplex1/out.png", dpi=600, bbox_inches="tight")
    final_df.to_csv("Output/Multiplex1/out.csv")
    
stitch(number_replicates)
os.system("rm -R Input/Multiplex1/replicates")