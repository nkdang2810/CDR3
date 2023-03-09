import argparse
import networkx as nx
from pathlib import Path
import pandas as pd
import numpy as np
from itertools import combinations
from tqdm import tqdm
import plotly.graph_objects as go


def main(args):
    blosum = pd.read_csv(
        "/rsrch4/home/mol_cgenesis/nkdang/CDR3/reference/BLOSUM62.csv")
    blosum = blosum.astype(int)

    def check_group(s1, s2):
        for i, _ in enumerate(s1):
            if blosum.loc[s1[i], s2[i]] <= 0:
                return False
        return True

    def hamstring_similarity(x, y):
        return sum([c1 == c2 for c1, c2 in zip(list(x), list(y))])/len(x)

    # from plotly.colors import n_colors
    # cluster_size = sorted([int(i.stem) for i in list(Path("/rsrch4/home/mol_cgenesis/EMC_BIC_rsrch4/nkdang/CDR3/results/SKCM/heavy/clustering_convergence/clusters").glob("*"))])
    # cluster_size_mapper = dict(zip(cluster_size, n_colors('rgb(255, 200, 200)', 'rgb(255, 0, 0)',len(cluster_size) , colortype='rgb')))
    # cluster_size_mapper[0] = 'rgb(0,0,0)'
    # cluster_size_mapper

    in_dir = Path(args.input)
    cancer = in_dir.parent.parent.name
    cluster_type = in_dir.name
    in_dir = in_dir.joinpath("clusters")
    G = nx.Graph()

    cluster_counter = 1

    for file in tqdm(in_dir.glob("**/*.csv"), total=len(list(in_dir.glob("**/*.csv")))):
        temp = pd.read_csv(file)
        # print(temp)
        # if temp.empty:
        #     continue
        # temp = temp.drop(['sample_id', 'row_index'], axis=1)
        temp['patient_count'] = file.parent.name
        temp['vdj_count'] = len(temp)
        temp['name'] = temp['patient_id']+"_"+temp["v_call"] + \
            "_"+temp["j_call"]+"_"+temp["junction_aa"]
        temp['cluster'] = cluster_counter
        G.add_edges_from([(i, j) for i, j in list(combinations(temp['name'], 2))], cluster=str(
            cluster_counter), cluster_size=file.parent.stem)
        nx.set_node_attributes(G, temp.set_index('name').to_dict('index'))
        cluster_counter += 1
    filters = [lambda x: not G.has_edge(x[0][0], x[1][0]),
               lambda x: len(x[0][1]['junction_aa']) == len(
                   x[1][1]['junction_aa']),
               lambda x: x[0][1]['v_call'] == x[1][1]['v_call'],
               lambda x: x[0][1]['j_call'] == x[1][1]['j_call'],
               lambda x: check_group(
                   x[0][1]['junction_aa'], x[1][1]['junction_aa']),
               lambda x: hamstring_similarity(x[0][1]['junction_aa'], x[1][1]['junction_aa']) >= 0.8]
    result = combinations(list(G.nodes(data=True)), 2)
    for f in filters:
        result = filter(f, result)
    result = list(result)
    G.add_edges_from([(i[0][0], i[1][0]) for i in result], cluster='0')

    # nx.set_node_attributes(G, {k:{"pos":v} for k,v in nx.shell_layout(G).items()})
    nx.set_node_attributes(
        G, {k: {"pos": v} for k, v in nx.nx_agraph.graphviz_layout(G, prog="circo").items()})
    # nx.set_node_attributes(G, {k:{"pos":v} for k,v in nx.kamada_kawai_layout(G).items()})
    edge_remove = [(i, j)
                   for i, j, a in G.edges(data=True) if a['cluster'] == '0']
    for e in edge_remove:
        G.remove_edge(e[0], e[1])

    def color_mapper(size):
        if size < 5:
            return 'rgb(0,0,0)'
        elif size < 10:
            return 'rgb(0,0,200)'
        elif size < 20:
            return 'rgb(0,200,0)'
        else:
            return 'rgb(200,0,0)'

    node_trace = []
    for node, att in G.nodes(data=True):
        x_pos, y_pos = att['pos']
        node_text = f"{att['v_call']}<br>{att['j_call']}<br>{att['junction_aa']}<br>Normalized count: {att['normalized_count']} (multiplied by 1e3)<br>Patient: {att['patient_id']}<br>Patient count in cluster: {att['patient_count']}<br>VDJs count: {att['vdj_count']}"
        node_trace.append(dict(type='scatter',
                               x=[x_pos],
                               y=[y_pos],
                               mode='markers',
                               hoverinfo='text',
                               text=node_text,
                               marker=dict(colorscale='Reds',
                                           color=color_mapper(
                                               int(att['patient_count'])),
                                           size=np.maximum(
                                               1, np.log(att['normalized_count']*100)),
                                           line=dict(width=1,
                                                     color='rgb(255,255,255)'))))

    edge_trace = []
    for node1, node2, att in G.edges(data=True):
        x_pos = [G.nodes[node1]['pos'][0], G.nodes[node2]['pos'][0]]
        y_pos = [G.nodes[node1]['pos'][1], G.nodes[node2]['pos'][1]]
        line_dash = "solid" if G.nodes[node1]['junction_aa'] == G.nodes[node2]['junction_aa'] else "dot"
        line_width = 2 if G.nodes[node1]['junction_aa'] == G.nodes[node2]['junction_aa'] else 1
        edge_trace.append(dict(type='scatter',
                               x=x_pos,
                               y=y_pos,
                               mode='lines',
                               line=dict(width=line_width,
                                         color=color_mapper(
                                             int(att['cluster_size'])),
                                         dash=line_dash)))
    layout = go.Layout(height=1600,
                       title=f'<br>{cluster_type} of {cancer}',
                       titlefont_size=16,
                       showlegend=False,
                       hovermode='closest',
                       margin=dict(b=20, l=5, r=5, t=40),
                       xaxis=dict(showgrid=False, zeroline=False,
                                  showticklabels=False),
                       yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))

    fig = go.Figure(data=edge_trace+node_trace, layout=layout)
    viz_path = in_dir.parent.parent.parent.parent.parent.joinpath(
        "viz").joinpath(cluster_type)
    viz_path.mkdir(parents=True, exist_ok=True)
    fig.write_html(viz_path.joinpath(f"{cancer}.html"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='CDR3 analysis')
    parser.add_argument('input',
                        type=str,
                        help='input directory path')
    args = parser.parse_args()
    main(args)
