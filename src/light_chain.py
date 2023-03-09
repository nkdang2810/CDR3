import requests
import pandas as pd
from pathlib import Path
import json
from tqdm import tqdm
import numpy as np
import plotly.express as px
import logging
from sklearn.cluster import DBSCAN
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import random


class LightChainAnalysis():
    """
    Calculate the matrix score of a cluster of heavy chain with all other light chains in one cancer type
    Select the best candidate light chains based on the CUSTOMIZED CORRELATION and FREQUENCY SCORE between light chain and heavy chain
    * CUSTOMIZED CORRELATION: group the nearby score using DBSCAN with min distance of 0.001 and min sample of 1. If there're only 2 data points left, the correlation score is -1
    * FREQUENCY SCORE: square of distance of heavychain frequency and lightchain frequecy
    Filters:
    - The frequency of lightchain must be 
    Parameters:
        - 
    """

    def __init__(self, group_cluster_path, data_heavy_clean_path, data_light_clean_path, cancer) -> None:
        self.group_cluster_path = group_cluster_path
        self.data_heavy_clean_path = data_heavy_clean_path
        self.data_light_clean_path = data_light_clean_path
        self.cancer = cancer

        self.min_distance = 0.3
        self.color_list = ['#1f77b4',  # muted blue
                           '#ff7f0e',  # safety orange
                           '#2ca02c',  # cooked asparagus green
                           '#d62728',  # brick red
                           '#9467bd',  # muted purple
                           '#8c564b',  # chestnut brown
                           '#e377c2',  # raspberry yogurt pink
                           '#7f7f7f',  # middle gray
                           '#bcbd22',  # curry yellow-green
                           '#17becf'   # blue-teal
                           ]

    def __draw_correlation__(self, path, light_chain, heavy_chain):

        ls1 = self.F_L_high_expression.loc[light_chain].reset_index(drop=True)
        ls2 = pd.Series(self.F_Ht)
        cor = self.customized_correlation(ls1, ls2)

        if cor > 0:
            fig = make_subplots(rows=1, cols=2, subplot_titles=(
                f"Original data - Cov: {ls1.corr(ls2)}", f"Merged data - Cov: {cor}"))

            points = np.array([[x, y] for x, y in zip(ls2, ls1)])
            clustering = DBSCAN(eps=0.001, min_samples=1).fit(points)
            labels = clustering.labels_.astype(str)
            colors = random.sample(self.color_list, len(set(labels)))

            for idx, g in enumerate(set(labels)):
                indices = np.where(labels == g)
                x, y = np.mean(points[indices], axis=0)
                fig.add_trace(go.Scatter(x=[x], y=[y], mode="markers", marker_color=f'{colors[idx]}'), row=1, col=2
                              )
                fig.add_trace(
                    go.Scatter(x=ls2.iloc[indices], y=ls1.iloc[indices], mode="markers", marker_color=f'{colors[idx]}'), row=1, col=1
                )

            m = np.max(points)
            fig.update_yaxes(range=[-0.0005, m+0.0005])
            fig.update_xaxes(range=[-0.0005, m+0.0005])
            fig.update_layout(title_text=heavy_chain+"   "+"_".join(light_chain),
                              height=1000, width=2000, showlegend=False)

            viz_path = Path(str(path).replace(
                "heavy", "light")).with_suffix("")
            viz_path.mkdir(parents=True, exist_ok=True)
            fig.write_html(viz_path.joinpath("_".join(light_chain)+".html"))

    def customized_correlation(self, ls1: np.array, ls2: np.array) -> float:
        """
        Input: list1 and list2
        Output: correlation value after grouping nearby points 
        """
        points = np.array([[x, y] for x, y in zip(ls1, ls2)])
        clustering = DBSCAN(eps=0.001, min_samples=1).fit(points)
        labels = clustering.labels_
        ls1_new = []
        ls2_new = []
        for g in set(labels):
            indices = np.where(labels == g)
            group = points[indices]
            x, y = np.mean(group, axis=0)
            ls1_new.append(x)
            ls2_new.append(y)

        # the correlation is meaningless if there're only 2 data points
        if len(ls1_new) < 3:
            return -1
        return pd.Series(ls1_new).corr(pd.Series(ls2_new))

    def execute(self):
        """
        Find corresponding light chains for each heavy chain cluster
        """

        light_chains = pd.read_feather(self.data_light_clean_path)
        heavy_chains = pd.read_feather(self.data_heavy_clean_path)

        patients_light_count = pd.DataFrame(light_chains.value_counts(
            "patient_id"), columns=['total_expression'])
        patients_heavy_count = pd.DataFrame(heavy_chains.value_counts(
            "patient_id"), columns=['total_expression'])

        light_groups = pd.DataFrame(light_chains.value_counts(
            ['v_call', 'j_call', 'junction_aa', 'patient_id']), columns=['light_chain_expression'])
        light_groups = light_groups.reset_index()

        # p_map(self.__cluster_helper__, groups)
        results = pd.DataFrame(
            columns=['score', 'correlation', 'distance', 'light', 'heavy', 'patient_count'])
        for f in tqdm(self.group_cluster_path.glob("**/*.csv"), total=len(list(self.group_cluster_path.glob("**/*.csv")))):
            heavy_cluster = pd.read_csv(f)
            heavy_cluster = heavy_cluster.astype({"patient_id": str})
            heavy_chain = f.stem

            heavy_expression = pd.DataFrame(heavy_cluster.value_counts(
                ['patient_id']).sort_index(), columns=['cluster_expression'])
            heavy_expression.index = heavy_expression.index.get_level_values(0)
            heavy_expression = heavy_expression.merge(
                patients_heavy_count, how='left', left_index=True, right_index=True)
            heavy_expression['percent'] = heavy_expression['cluster_expression'] / \
                heavy_expression['total_expression']
            patients = heavy_expression.index.tolist()

            F_Ht = heavy_expression['percent'].values

            # FREQUENCY LIGHT CHAIN MATRIX
            groups = light_groups[light_groups['patient_id'].isin(
                patients)]
            groups = groups.sort_values(by=["patient_id"])
            groups = groups.merge(patients_light_count, how='left',
                                  left_on='patient_id', right_on='patient_id')
            groups['percent'] = groups['light_chain_expression'] / \
                groups['total_expression']

            F_L = pd.crosstab(index=[groups["v_call"], groups["j_call"], groups['junction_aa']],
                              columns=groups['patient_id'], values=groups['percent'], aggfunc='sum').fillna(0)
            # remove low signal light chain
            F_L_high_expression = F_L[F_L.sum(axis=1) > F_Ht.sum()/5]

            # additional filters? (half number of values must be non-zero number)
            F_L_high_expression = F_L_high_expression[(F_L_high_expression < 0.00001).astype(
                int).sum(axis=1) < len(F_L_high_expression.columns)/2]

            F_L_high_expression["correlation"] = np.array([row.reset_index(drop=True).corr(
                pd.Series(F_Ht)) for _, row in F_L_high_expression.iterrows()])
            F_L_high_expression["fined_correlation"] = np.array(
                [self.customized_correlation(row, F_Ht) for _, row in F_L_high_expression.iterrows()])

            F_L_high_expression_high_correlation = F_L_high_expression[
                F_L_high_expression['fined_correlation'] > 0]
            score = pd.DataFrame(np.power(F_L_high_expression_high_correlation.drop(
                ['correlation', "fined_correlation"], axis=1)-F_Ht, 2).sum(axis=1), columns=['score'])
            cov_match = []
            for _, row in F_L_high_expression_high_correlation.iterrows():
                cov_match.append(row['fined_correlation'])
            score['correlation'] = cov_match
            score['distance'] = np.power(
                np.power(score['score'], 2) + np.power(score['correlation']-1, 2), 1/2)

            # distance must be lesser than {}
            score = score[score['distance'] < self.min_distance]
            score = score.sort_values('distance')
            if not score.empty:
                outdir = Path(str(f).replace("heavy", "light"))
                outdir.parent.mkdir(parents=True, exist_ok=True)
                score.to_csv(outdir)

                # score viz
                fig = px.scatter(score, x="score", y="correlation", width=800, height=800, hover_data=[
                    score.index.values], title=f"Scoring matrix: {self.cancer} - {heavy_chain}")
                fig.update_yaxes(range=[-0.01, 1.01])
                fig.update_traces(
                    hovertemplate='Light chain: %{customdata[0]}<br>Score: %{x}<br>Correlation: %{y}')
                fig.update_yaxes(range=[0.45, 1.05])
                fig.update_xaxes(range=[-0.0001, 0.001])
                fig.write_html(outdir.with_suffix(".html"))

            # correlation viz
            self.F_L_high_expression = F_L_high_expression
            self.F_Ht = F_Ht
            for light, row in score.iterrows():
                self.__draw_correlation__(f, light, heavy_chain)

            # report
            score['light'] = ["_".join(i) for i in score.index]
            score['heavy'] = [heavy_chain] * len(score)
            score['patient_count'] = [f.parent.stem] * len(score)
            results = pd.concat([results, score], axis=0, ignore_index=True)
        logging.info("Writing summary file ...")
        results = results.sort_values("distance")
        results.to_csv(Path(str(f).replace(
            "heavy", "light")).parent.parent.parent.parent.joinpath('summary.csv'), index=False)
        logging.info("Finish.")
