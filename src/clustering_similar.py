import subprocess
import pandas as pd
import logging
from p_tqdm import p_map
from collections import Counter

import numpy as np
from scipy.spatial.distance import pdist, squareform
from itertools import combinations
import networkx as nx


class ClusteringSimilar():
    """
    Use exhaustive clustering algorithm to cluster junction_aa of same v_call, j_call records and same length of junctions aa
    """

    def __init__(self, data_clean_path, patient_min, similarity_threshold, outdir) -> None:
        logging.info("SIMILAR JUNCTION-REGION CLUSTERING")
        self.data_clean_path = data_clean_path
        self.patient_min = patient_min

        self.clustering_similar_path = outdir.joinpath(
            "clustering_similar")
        self.clustering_similar_path.mkdir(parents=True, exist_ok=True)
        # self.temp_path = self.clustering_similar_path.joinpath("temp")
        # self.temp_path.mkdir(parents=True, exist_ok=True)
        self.cluster_path = self.clustering_similar_path.joinpath("clusters")
        self.cluster_path.mkdir(parents=True, exist_ok=True)

        self.similarity_threshold = similarity_threshold

    def __cluster_helper__(self, df):
        """
        Assign cluster for each records based on the similarity threshold using clique finding algorithm,
        one value can be assigned to one cluster
        """
        def hamstring_similarity(x, y):
            return sum([c1 == c2 for c1, c2 in zip(list(x), list(y))])/len(x)

        G = nx.Graph()
        G.add_edges_from([(i, j) for i, j in list(combinations(df.index, 2)) if (hamstring_similarity(df.loc[i]['junction_aa'], df.loc[j]
                         ['junction_aa']) >= self.similarity_threshold)])
        clusters = []
        while G.nodes():
            cliques = sorted(list(nx.find_cliques(G)),
                             key=len, reverse=True)
            biggest_cluster = cliques[0]
            if len(biggest_cluster) > 1:
                clusters.append(biggest_cluster)
            G.remove_nodes_from(biggest_cluster)

        df_summary = pd.DataFrame()
        for c in clusters:
            df_out = df.loc[c]
            patient_uniq = df_out['patient_id'].nunique()
            if patient_uniq >= self.patient_min:
                cluster_path = self.cluster_path.joinpath(
                    str(patient_uniq))
                cluster_path.mkdir(parents=True, exist_ok=True)

                name = f"{df_out.head(1)['v_call'].values[0]}_{df_out.head(1)['j_call'].values[0]}_{df_out.head(1)['junction_aa'].values[0]}"
                df_out = df_out.drop(
                    ['junction_aa_length'], axis=1)
                df_out.to_csv(cluster_path.joinpath(
                    name+".csv"), index=False)

                df_out_sum = df_out.groupby(['v_call', 'j_call']).aggregate(
                    {"consensus_count": 'mean', "patient_id": lambda x: "_".join(set(x))}).reset_index()
                df_out_sum['junction_aa'] = Counter(
                    df_out['junction_aa']).most_common()[0][0]
                df_summary = pd.concat([df_summary, df_out_sum])

        if not df_summary.empty:
            return df_summary

    def cluster(self):
        """
        Assign cluster for each records based on the similarity threshold
        """
        if not self.clustering_similar_path.joinpath("summary.csv").is_file():
            logging.info(
                f"Clustering based on v_call, j_call and similar junction_aa fields ...")
            df = pd.read_feather(self.data_clean_path)
            df['junction_aa_length'] = df['junction_aa'].apply(
                lambda x: len(x))

            groups = [i for _, i in df.groupby(
                ['v_call', 'j_call', 'junction_aa_length']) if i['patient_id'].nunique() >= self.patient_min]
            logging.info(
                f"Total clusters based on v_gene, j_gene and length of junction_aa: {len(groups):,d}")

            # write clusters using multi-processing with tqdm wrapper
            logging.info(f"Assigning clusters ...")
            summary_list = p_map(self.__cluster_helper__, groups)
            total_clusters = len(list(self.cluster_path.glob('**/*.csv')))
            logging.info(
                f"Total clusters based on v_gene, j_gene and similar junction_aa: {total_clusters:,d}")

            # generate summary file
            logging.info("Writing summary file ...")
            if len(summary_list) != 0:
                summary = pd.concat(summary_list, ignore_index=True)
                summary['patient_count'] = summary['patient_id'].str.count(
                    "_")+1
                summary = summary.sort_values(
                    by=['patient_count'], ascending=False)
                summary = summary.round(5)
                summary = summary[['v_call', 'j_call', 'junction_aa',
                                   'consensus_count', 'patient_count', 'patient_id']]
                summary.to_csv(self.clustering_similar_path.joinpath(
                    "summary.csv"), index=False)
            else:
                logging.info("There's no similar VDJ cluster.")
        logging.info(f"Reports are at {self.clustering_similar_path}")
        logging.info(f"Finish clustering.\n")
