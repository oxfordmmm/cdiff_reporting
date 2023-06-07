from argparse import ArgumentParser
import os, sys, time
import json
import numpy as np
import pandas as pd
import networkx as nx
import scipy.cluster.hierarchy as shc
from scipy.spatial.distance import squareform
from Bio import Phylo
import matplotlib.pyplot as plt

def read_json(dir: str, ids: list) -> dict:
    data = {}

    sys.stdout.write("Reading in JSON files\n")
    start = time.time()
    
    for filename in os.listdir(dir):
        if filename.endswith("cgmlst.json"):
            filebase = filename.replace("_cgmlst.json", "")
            if not filebase in ids:
                continue

            filepath = os.path.join(dir, filename)
            with open(filepath, "r") as file:
                try:
                    json_data = json.load(file)
                    data[filebase] = json_data
                except json.JSONDecodeError:
                    print(f"Error parsing JSON file: {filename}")
    
    end = time.time()
    sys.stdout.write("Seconds to read in files: {}\n".format(end - start))

    return data

def calculate_distances(cgmlst_profiles: dict) -> list:
    #list of bad genes to exclude
    excludeList = ['CD630_17960', 'CD630_23930', 'CD630_25030', 'CD630_08260', 'CD630_09010', 'CD630_12750', 'CD630_22750', 'CD630_27680', 'CD630_11800', 'CD630_12080', 'CD630_15400', 'CD630_16950', 'CD630_20730', 'CD630_32240', 'CD630_34030']
	
    samples= list(cgmlst_profiles.keys())[:50]
    distances = []

    sys.stdout.write("Comparing profiles\n")
    start = time.time()

    for i in range(0, len(samples)):
        profile1 = cgmlst_profiles[samples[i]]
        alleles1 = [profile1['alleles'][k] for k in profile1['alleles'].keys() if k not in excludeList]
        for j in range(i, len(samples)):
            profile2 = cgmlst_profiles[samples[j]]
            alleles2 = [profile2['alleles'][k] for k in profile2['alleles'].keys() if k not in excludeList]
	    
            compared = len([True for a1,a2 in zip(alleles1,alleles2) if a1 != "" and a2 != ""])
            diff = len([True for a1,a2 in zip(alleles1,alleles2) if a1 != "" and a2 != "" and a1 != a2])
	    
            distances.append([samples[i], samples[j], compared, diff])
	
    end = time.time()
    sys.stdout.write("Seconds to compare all pairs: {}\n".format(end - start))

    return distances

def to_newick_iter(node, newick: list, parentdist: float, leaf_names: list) -> list:
    if node.is_leaf():
        return newick + [f'{leaf_names[node.id]}:{parentdist - node.dist}']

    if len(newick) > 0:
        newick.append(f'):{parentdist - node.dist}')
    else:
        newick.append(');')
    newick = to_newick_iter(node.get_left(), newick, node.dist, leaf_names)
    newick.append(',')
    newick = to_newick_iter(node.get_right(), newick, node.dist, leaf_names)
    newick.append('(')
    return newick

def to_newick(linkage, leaf_names: list):
    T = shc.to_tree(linkage)
    newick_list = to_newick_iter(T, [], T.dist, leaf_names)
    return ''.join(newick_list[::-1])

def create_clusters(samples, distances, cutoff: int):
    G = nx.Graph()
    for s in samples:
        G.add_node(s)
    
    for s, t, compared, diff in distances:
        if s != t and diff <= cutoff and compared > 1000:
            G.add_edge(s, t, weight=diff)
    
    components = nx.connected_components(G)
    distances_df = pd.DataFrame(distances, columns=['sample_1', 'sample_2', 'compared', 'differences'])

    clusters_dict = {}
    clusters = [] # will be list of (cluster, sample). 0 means no cluster
    cluster_i = 1

    # Now need to form 
    newicks = []
    for component in components:
        component = list(component)
        if len(component) <= 1:
            for sample in component:
                clusters.append([0, sample])
            continue

        clusters_dict[cluster_i] = component
        for sample in component:
            clusters.append([cluster_i, sample])
        cluster_i += 1
        comp_df = distances_df.query('sample_1 in @component').query('sample_2 in @component')
        indexes = {sample:index for index,sample in enumerate(component)}


        dist_matrix = np.zeros((len(component), len(component)), dtype=int)
        for _, row in comp_df.iterrows():
            u = row['sample_1']
            v = row['sample_2']
            diff = row['differences']
            dist_matrix[indexes[u], indexes[v]] = diff
            dist_matrix[indexes[v], indexes[u]] = diff

        condensed_dist = squareform(dist_matrix)
        Z = shc.linkage(condensed_dist, 'single', optimal_ordering=True)

        newick = to_newick(Z, component)
        newicks.append(newick)
    
    return newicks, pd.DataFrame(clusters, columns=['cluster_number', 'id']), clusters_dict

if __name__ == '__main__':
    parser = ArgumentParser(description='Produce cgmlst distance matrix and cluster trees.')
    parser.add_argument('-f', '--fasta_list_tsv', required=True,
                        help='tsv with overall consensus fasta list')
    parser.add_argument('-s', '--samples_json_dir', required=True,
                        help='Directory containing cgmlst json files')
    parser.add_argument('-o', '--output_dir', required=True,
                        help='Ouput directory')
    parser.add_argument('-c', '--cutoff', required=False,
                        default=20,
                        help='cgmlst distance to use as cutoff for clusters')
    args = parser.parse_args()
    consensus_fastas =  pd.read_csv(args.fasta_list_tsv, sep='\t', names=['id', 'path'])
    samples_dir = args.samples_json_dir
    output_dir = args.output_dir
    cutoff = int(args.cutoff)

    consensus_ids = list(consensus_fastas['id'])
    print("consensus_ids")
    print(consensus_ids)
    profiles = read_json(samples_dir, consensus_ids)
    distances = calculate_distances(profiles)
    newicks, clusters, clusters_dict = create_clusters(profiles.keys(), distances, cutoff)

    distances_df = pd.DataFrame(distances, columns=['sample_1', 'sample_2', 'compared', 'differences'])
    distances_df.to_csv(f"{output_dir}/distances.csv", index=False)
    clusters.sort_values('cluster_number').to_csv(f"{output_dir}/clusters.tsv",
                                           index=False,
                                           sep='\t')

    for cluster, newick in enumerate(newicks):
        f = open(f"{output_dir}/cluster_{cluster+1}_cgmlst_scaled.newick", 'w')
        f.write(newick)
        f.close()
    
    for cluster, ids in clusters_dict.items():
        consensus_fastas.query('id in @ids').to_csv(f"{output_dir}/cluster_{cluster}_consensus.tsv", sep='\t', header=False, index=False)

