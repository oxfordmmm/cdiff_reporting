#Usage: python3 draw_all_cluster_trees.py path/to/cfml/newick/files/ (in cluster_ml of the output stem from the .ini file)
import os
from Bio import Phylo
import matplotlib.pyplot as plt

# Set the directory path containing the Newick files
dir_path = "/mnt/scratch/soft/runListCompare_rev/out_rlc_seq3_rerun/cluster_ml/"

# Loop through all the files in the directory with the .newick extension (cfml labelled tree)
for file_name in os.listdir(dir_path):
    if file_name.endswith("_iqtree_scaled.tree"):
        # Load the Newick file
        file_path = os.path.join(dir_path, file_name)
        tree = Phylo.read(file_path, "newick")

        # Draw the tree using Phylo module
        tree.rooted = True
        Phylo.draw(tree)
        #Phylo.draw(tree, branch_labels=lambda c: c.branch_length)

        # Save the tree as a PNG file
        file_base = os.path.splitext(file_name)[0]
        png_file = file_base + ".png"
        png_path = os.path.join(dir_path, png_file)
        plt.savefig(png_path)

        # Clear the plot to prepare for the next tree
        plt.clf()
