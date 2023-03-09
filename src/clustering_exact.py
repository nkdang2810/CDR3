import pandas as pd
import logging
from p_tqdm import p_map


class ClusteringExact():
    """
    Clustering based on v_gene, j_gene and the exact match of junction_aa
    """

    def __init__(self, data_clean_path, patient_min, outdir) -> None:
        logging.info("EXACT JUNCTION-REGION CLUSTERING")
        self.data_clean_path = data_clean_path
        self.patient_min = patient_min

        self.clustering_exact_path = outdir.joinpath(
            "clustering_exact")
        self.clustering_exact_path.mkdir(parents=True, exist_ok=True)
        self.cluster_path = self.clustering_exact_path.joinpath("clusters")
        self.cluster_path.mkdir(parents=True, exist_ok=True)

    def __cluster_helper__(self, df):
        out = df.groupby(['v_call', 'j_call', 'junction_aa']).aggregate(
            {"consensus_count": "mean", "patient_id": lambda x: "_".join(set(x))}).reset_index()
        out['patient_count'] = out['patient_id'].str.count("_")+1

        out_path = self.cluster_path.joinpath(
            str(out['patient_count'].values[0]))
        out_path.mkdir(parents=True, exist_ok=True)

        df.to_csv(out_path.joinpath(
            f"{out['v_call'].values[0]}_{out['j_call'].values[0]}_{out['junction_aa'].values[0]}.csv"), index=False)
        return out

    def cluster(self):
        if not self.clustering_exact_path.joinpath("summary.csv").is_file() or True:
            logging.info(
                f"Clustering based on v_call, j_call and exact match of junction_aa ...")
            df = pd.read_feather(self.data_clean_path)
            total_patient = df['patient_id'].nunique()
            logging.info(f"Number of patients: {total_patient:,d}")

            logging.info(
                f"Total clusters based on v_gene, j_gene and junction_aa: {(len(df.groupby(['v_call', 'j_call', 'junction_aa']))):,d}")

            # cluster_patientid = df.groupby('cluster_exact')[
            #     "patient_id"].count().reset_index().set_index('cluster_exact')
            # cluster_patientid = cluster_patientid[cluster_patientid['patient_id']
            #                                       >= self.patient_min]
            # df_filtered = df[df['cluster_exact'].isin(cluster_patientid.index)]

            groups = [i for _, i in df.groupby(['v_call', 'j_call', 'junction_aa']) if (
                len(i) > 1) and (i['patient_id'].nunique() >= self.patient_min)]
            logging.info(
                f"Total clusters having minimum of {self.patient_min} members: {len(groups):,d}")

            # write the clusters
            logging.info(f"Assigning clusters ...")
            summary_list = p_map(self.__cluster_helper__, groups)
            summary = pd.concat(summary_list, ignore_index=True)
            if not summary.empty:
                logging.info(f"Writing summary file ...")
                summary = summary.sort_values(
                    by=['patient_count'], ascending=False)
                summary = summary[['v_call', 'j_call', 'junction_aa',
                                   'consensus_count', 'patient_count', 'patient_id']]
                summary = summary.round(5)
                summary.to_csv(self.clustering_exact_path.joinpath(
                    "summary.csv"), index=False)
            else:
                logging.info("There's no exact VDJ cluster.")
        logging.info(f"Reports are at {self.clustering_exact_path}")
        logging.info(f"Finish clustering.\n")
