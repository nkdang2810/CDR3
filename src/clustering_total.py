
from multiprocessing import Pool
from pathlib import Path
import pandas as pd
from tqdm import tqdm
from clustering_similar import clustering
import subprocess
from itertools import groupby
import re

CLUSTER_COVERAGE = 0.03
SIMILARITY_THRESHOLD = 0.8


def clustering(args):
    """
    Use CDHIT to cluster junction_aa of same v_call and j_call records
    """
    df, temp_path = args

    name = f"{df.head(1).index[0]}_{df.head(1)['v_call'].values[0]}_{df.head(1)['j_call'].values[0]}"
    new_path = temp_path.joinpath(name)
    new_path.mkdir(parents=True, exist_ok=True)

    if len(df) > 1:
        with open(new_path.joinpath(f"{name}.fasta"), "w") as file:
            for idx, row in df.iterrows():
                file.write(
                    f">{idx}@{name}_{row['junction_aa']}\n")
                file.write(f"{row['junction_aa']}\n")

        subprocess.run(["cd-hit", "-g", "1", "-c", "0.8", "-M",
                        "0", "-l", "4", "-d", "0", "-T", "0", "-i", str(new_path.joinpath(f"{name}.fasta")), "-o", str(new_path.joinpath(name))], stdout=subprocess.DEVNULL)

        with open(new_path.joinpath(f"{name}.clstr"), "r") as file:
            groups = [list(group) for key, group in groupby(
                file, lambda line: line.startswith(">Cluster")) if not key]
            ids = [re.findall(r">(.*)\@{1}", "".join(g)) for g in groups]
        return ids
    else:
        return df.index.tolist()


def get_top(data, size, out_path):
    temp = data.groupby(['v_call', 'j_call', 'cluster']).agg(
        {"sample_id": lambda x: x.nunique()})
    temp['percent'] = temp['sample_id']/size
    temp = temp.reset_index(level=2)
    tops = temp.merge(data.groupby(['v_call', 'j_call'])[
        "junction_aa"].first(), how="left", left_index=True, right_index=True)
    tops = tops.sort_values("sample_id", ascending=False)
    tops = tops.reset_index()

    top_clusters = tops[tops['percent'] > 0.05]['cluster'].values
    tops = tops[['v_call', 'j_call', 'junction_aa', 'percent']]
    tops.to_csv(out_path.joinpath("summary.tsv"), sep='\t', index=False)

    # return only the most common clusters
    return top_clusters


cancers = ["LUSC", "LUAD"]

source_path = Path(f"/rsrch4/scratch/mol_cgenesis/nkdang/CDR3/trust4")
results_path = Path(
    f"/rsrch4/scratch/mol_cgenesis/nkdang/CDR3/results/heavy")

results_all_path = results_path.joinpath("all")
results_all_path.mkdir(parents=True, exist_ok=True)

all_data = []
all_size = 0
for cancer in cancers:
    cancer_size = len(list(source_path.joinpath(cancer).glob('**/*_airr.tsv')))
    cancer_data = pd.read_csv(results_path.joinpath(
        cancer).joinpath("clusters.tsv"), sep='\t')

    all_data.append(cancer_data)
    all_size += cancer_size

print(f"Analyzing all cancer types")

all_cancer_data = pd.concat(all_data, ignore_index=True)
all_cancer_data = all_cancer_data.drop("cluster", axis=1)

print(f"Number of patients {all_size}")

all_cancer_data['junction_aa_length'] = all_cancer_data['junction_aa'].apply(
    lambda x: len(x))
all_groups = all_cancer_data.groupby(
    ['v_call', 'j_call', 'junction_aa_length'])

# multi-processing
with Pool(processes=None) as pool:

    # have your pool map the file names to dataframes
    indices_list = list(tqdm(
        pool.imap(clustering, [(i, results_all_path.joinpath("temp")) for _, i in all_groups]), total=len(all_groups)))


# # WRITE CLUSTERS TO FILES
print(f"Asigning clusters ...")
cluster_indices = [0] * len(all_cancer_data)
cluster_count = 1
for ls in tqdm(indices_list, total=len(indices_list)):
    if any(isinstance(x, list) for x in ls):
        for i in ls:
            for j in i:
                cluster_indices[int(j)] = cluster_count
            cluster_count += 1
    else:
        cluster_indices[int(ls[0])] = cluster_count
    cluster_count += 1

all_cancer_data['cluster'] = cluster_indices
all_cancer_data = all_cancer_data.drop("junction_aa_length", axis=1)
all_cancer_data.to_csv(results_all_path.joinpath(
    "clusters.tsv"), sep='\t', index=False)

# CALCULATE THE COVERAGE
temp = all_cancer_data.groupby('cluster').agg(
    {"sample_id": lambda x: x.nunique()})
temp = temp.rename(columns={'sample_id': 'sample_count'})
temp['percent'] = temp['sample_count']/cancer_size
# temp = temp.reset_index(level=2)
tops = temp.merge(all_cancer_data.groupby('cluster')[
    ["v_call", "j_call", "junction_aa"]].first(), how="left", left_index=True, right_index=True)
tops = tops.sort_values("sample_count", ascending=False)
tops = tops.reset_index()

tops = tops[tops['percent'] > CLUSTER_COVERAGE]
top_clusters = tops['cluster'].values

tops = tops[['v_call', 'j_call', 'junction_aa', 'sample_count', 'percent']]
tops.to_csv(results_all_path.joinpath("summary.tsv"), sep='\t', index=False)

print(f"Writing clusters ...")
for c in tqdm(top_clusters, total=len(top_clusters)):
    temp = all_cancer_data[all_cancer_data['cluster'] == c]
    temp = temp.drop("cluster", axis=1)
    temp = temp.reset_index(drop=True)
    temp.to_csv(results_all_path.joinpath(
        f"{temp['v_call'][0]}_{temp['j_call'][0]}_{temp['junction_aa'][0]}.tsv"), sep='\t', index=False)
print(f"Finish.")
