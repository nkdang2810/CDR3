# CDR3 project
## Overview
Finding the common VDJs between patients from TRUST4 output

python /rsrch4/home/mol_cgenesis/nkdang/CDR3/src/main.py

---
## Input
-i: path to TRUST4 output directory. The directory structure is /{sample_id}/*_airr.tsv

-g: path to metadata file. The metadata file contains the sample_id and patient_id columns (Make sure the sample_id and patient_id columns are same size and one-one relations)

-oL path to output directory

---
## Output
Output directory: /rsrch4/home/mol_cgenesis/EMC_BIC_rsrch4/nkdang/CDR3/results/{cancer}

There're 3 types of heavy-chain clustering: exact, similar and convergence. In each subdirectory, there's a summary file

The visualization directory: /rsrch4/home/mol_cgenesis/EMC_BIC_rsrch4/nkdang/CDR3/viz

---

## Methods
1. Preprocessing
    - Merge all _airr.tsv from the input directory. Remove the subgroup in V gene and J gene columns. Get the patient_id from metadata file. Calculate the `normalized count = consensus_count / total_count * 1000`
    - Keep columns sample_id, row_index, v_call, j_call, junction_aa, normalized_count for next steps
    - Separate `heavy_chains` and `light_chains` based on the V gene name. Light chains may be lamba or kappa.
    - Aggregate rows having same sample_id, v_call, j_call, junction_aa, patient_id. Merged rows have the **concatenated row_index** and **sum of normalized_count**
    - Three output categories: heavy_chains, light_chains and others
    /rsrch4/home/mol_cgenesis/EMC_BIC_rsrch4/nkdang/CDR3/results/{cancer}/heavy
    /rsrch4/home/mol_cgenesis/EMC_BIC_rsrch4/nkdang/CDR3/results/{cancer}/light
    /rsrch4/home/mol_cgenesis/EMC_BIC_rsrch4/nkdang/CDR3/results/{cancer}/others
2. Heavy-chain exact VDJs clustering
    - Cluster VDJs based on the exact matches of `V gene`, `J gene` and `Junction AA`.
    - Generate the summary file
    /rsrch4/home/mol_cgenesis/EMC_BIC_rsrch4/nkdang/CDR3/results/{cancer}/heavy/clustering_exact
3. Heavy-chain similar VDJs clustering (params: similarity_threshold)
    - Group VDJs based on the exact matches of `V gene`, `J gene` and `length of Junction AA`.
    - For each group, generate a graph G(N,E). Each node represent a unique VDJ (unique combination of V gene, J gene and Junction AA). Each edge represent a similarity of 2 Junction AA `above or equal the similarity_threshold`. Apply Bron–Kerbosch algorithm to find all maximum cliques in G. Record the biggest clique and remove all of its nodes from G. Repeat the algorithm until G is empty.
    - Generate the summary file
    /rsrch4/home/mol_cgenesis/EMC_BIC_rsrch4/nkdang/CDR3/results/{cancer}/heavy/clustering_similar
4. Heavy-chain convergence VDJs clustering
    - Similar to similar clustering method except the edge definition. Each edge represent a similarity of 2 Junction AA `above or equal the similarity_threshold` AND `same amino acid group for ALL amino acid positions`. Same amino acid group means the `BLOSUM62` score of these 2 amino acid is positive.
5. Heavy-chain clusters visualization
    - The network visualization are constructed by `networkx` and `plotlt`. Using layout `circo` from `graphviz`.
    - Size of node represents the normalized count
    - Solid line represents exact junction_aa, otherwise dot-line. 
    - Red color means clusters having patient count >=20, blue: from 10 to 20, green: from 5 to 10, black: < 5
    - Same v gene, j gene and junction_aa length VDJs tends to sit close together
    /rsrch4/home/mol_cgenesis/EMC_BIC_rsrch4/nkdang/CDR3/viz/{cancer}
6. Light-chain detection