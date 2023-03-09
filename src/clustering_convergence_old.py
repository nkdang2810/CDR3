from pathlib import Path
import logging
from p_tqdm import p_map
import pandas as pd
from collections import Counter
from itertools import combinations
import networkx as nx


class ClusteringConver():
    """
    Clustering based on v_gene, j_gene and the amino acid groups of junction_aa
    """

    def __init__(self, similar_cluster_path, patient_min,
                 blosum, outdir) -> None:

        logging.info("CONVERGENCE JUNCTION-REGION CLUSTERING")
        self.similar_cluster_path = similar_cluster_path
        self.blosum = pd.read_csv(blosum)

        self.patient_min = patient_min

        self.clustering_group_path = outdir.joinpath(
            "clustering_convergence")
        self.clustering_group_path.mkdir(parents=True, exist_ok=True)
        self.cluster_path = self.clustering_group_path.joinpath("clusters")
        self.cluster_path.mkdir(parents=True, exist_ok=True)

    def __cluster_helper__(self, filename: Path):
        df = pd.read_csv(filename, dtype={"sample_id": str,
                                          "v_call": str,
                                          "j_call": str,
                                          "junction_aa": str,
                                          "patient_id": str,
                                          "consensus_count": int,
                                          "frequency": float,
                                          "row_index": str})

        def check_group(s1, s2):
            for i, _ in enumerate(s1):
                if self.blosum.loc[s1[i], s2[i]] <= 0:
                    return False
            return True

        G = nx.Graph()
        G.add_edges_from([(i, j) for i, j in list(combinations(
            df['junction_aa'].unique(), 2)) if check_group(i, j)])

        clusters = []
        while G.nodes():
            cliques = sorted(list(nx.find_cliques(G)),
                             key=len, reverse=True)
            if len(cliques[0]) > 1:
                clusters.append(cliques[0])
            G.remove_nodes_from(cliques[0])

        df_summary = pd.DataFrame()
        for c in clusters:
            df_out = df[df['junction_aa'].isin(c)]
            patient_uniq = df_out['patient_id'].nunique()
            if patient_uniq >= self.patient_min:
                cluster_path = self.cluster_path.joinpath(
                    str(patient_uniq))
                cluster_path.mkdir(parents=True, exist_ok=True)

                name = f"{df_out.head(1)['v_call'].values[0]}_{df_out.head(1)['j_call'].values[0]}_{df_out.head(1)['junction_aa'].values[0]}"
                df_out.to_csv(cluster_path.joinpath(
                    name+".csv"), index=False)

                df_out_sum = df_out.groupby(['v_call', 'j_call']).aggregate(
                    {"consensus_count": 'sum', "frequency": "sum", "patient_id": lambda x: "_".join(set(x))}).reset_index()
                df_out_sum['junction_aa'] = Counter(
                    df_out['junction_aa']).most_common()[0][0]
                df_summary = pd.concat([df_summary, df_out_sum])
        return df_summary

    def cluster(self):
        """
        Assign cluster for each records based on the similarity threshold and amino acid groups
        """
        if not self.clustering_group_path.joinpath("summary.csv").is_file():
            logging.info(
                f"Clustering based on v_call, j_call, similar junction_aa fields and group of amino acids ...")

            # write clusters using multi-processing with tqdm wrapper
            logging.info(f"Assigning clusters ...")
            total_cases = len(
                list(self.similar_cluster_path.glob('**/*.csv')))
            summary_list = p_map(self.__cluster_helper__, self.similar_cluster_path.glob(
                '**/*.csv'), total=total_cases)
            total_clusters = len(list(self.cluster_path.glob('**/*.csv')))
            logging.info(
                f"Total clusters based on v_gene, j_gene and convergence junction_aa: {total_clusters:,d}")

            # generate summary file
            logging.info("Writing summary file ...")
            summary = pd.concat(summary_list, ignore_index=True)
            if not summary.empty:
                summary['patient_count'] = summary['patient_id'].str.count(
                    "_")+1
                summary = summary.sort_values(
                    by=['patient_count'], ascending=False)
                summary = summary.round(5)
                summary = summary[['v_call', 'j_call', 'junction_aa',
                                   'consensus_count', 'frequency', 'patient_count', 'patient_id']]
                summary.to_csv(self.clustering_group_path.joinpath(
                    "summary.csv"), index=False)
            else:
                logging.info("There's no convergence VDJ cluster.")
        logging.info(f"Finish clustering.\n")
        return self.cluster_path
