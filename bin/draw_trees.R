
library(tidyverse)
library(ggtree)
library(treeio)

args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]

tree_files <- list.files(dir, pattern = "\\.newick$", full.names = T)

get_node_depth <- function(tree, n) {
  to_search = c(n)
  depths = c(0)
  max_d = 0
  
  while (length(to_search) > 0) {
    m <- tail(to_search, 1)
    to_search <- head(to_search, -1)
    current_depth <- tail(depths, 1)
    depths <- head(depths, -1)
    
    this_branch_length <- tree %>% 
      filter(node == m) %>% 
      .$branch.length %>% 
      sum(na.rm = T)
    
    current_depth <- current_depth + this_branch_length
    max_d <- max(max_d, current_depth)
    children <- child(tree, m)$node
    to_search <- append(to_search, children)
    
    if (length(children) > 0) {
      for (i in 1:length(children)) {
        depths <- append(depths, current_depth)
      }
    }
  }
  
  return(max_d)
}

get_depth <- function(tree_tb) {
  root <- tree_tb %>% 
    filter(parent == node) %>% 
    .$node
  depth <- get_node_depth(tree_tb, root)
  return(depth)
}

# Scales tree so that total depth is 1 (or target_length)
normalise_tree <- function(tree_tb, target_length = 1) {
  depth <- get_depth(tree_tb)

  if(depth==0) {
    stop("Can't normalise when tree depth is 0")
  }
  
  tree_tb %>% 
    mutate(branch.length = branch.length * target_length / depth) %>% 
    return()
}

draw_tree <- function(newick_file) {
  print(sprintf("drawing %s", newick_file))
  tree <- read.newick(newick_file)
  n_leaves <- Ntip(tree)
  height <- min(3 + n_leaves, 18)
  text_size = if_else(n_leaves > 20, 3, 8)
  tree <- as_tibble(tree)
  
  depth <- get_depth(tree)
  print(sprintf("depth %s", depth))
  # tree <- normalise_tree(tree)

  scale_width <- 1
  scale_label <- 'distance'
  if(str_detect(newick_file, 'iqtree_scaled') | str_detect(newick_file, 'cf_scaled') | str_detect(newick_file, 'phyml_scaled')) {
    scale_width <- 10
    scale_label <- 'snp distance'
  }

  if(str_detect(newick_file, 'cgmlst_scaled')) {
    scale_width <- 1
    scale_label <- 'cgmlst distance'
  }
  
  gg_tr <- ggtree(as.treedata(tree)) +
    geom_tiplab(size=text_size) +
    # theme_tree2() +
    geom_treescale(width = scale_width) +
    xlab(scale_label) +
    theme(axis.title.x = element_text(size=18))

  if(depth < 7) {
    gg_tr <- gg_tr + scale_x_continuous(limits=c(NA, 10))
  }
  else{
    gg_tr <- gg_tr + hexpand(.4, direction = 1)
  }
  
  ggsave(str_replace(newick_file, '.newick', '.png'), plot=gg_tr, width=12, height=height, dpi=300)
  return(gg_tr)
}
draw_tree <- Vectorize(draw_tree)

trees <- draw_tree(tree_files)
