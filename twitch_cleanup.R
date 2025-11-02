# Code for cleaning data from the raw version available at: https://snap.stanford.edu/data/twitch_gamers.html

library(tidyverse)
library(igraph)

twitch_el <- read_csv("large_twitch_edges.csv")
twitch_el <- twitch_el |> 
  rename(node1 = numeric_id_1, node2 = numeric_id_2) |> 
  mutate(node1 = node1 + 1, node2 = node2 + 1) |> 
  filter(node1 < node2)

twitch_features <- read_csv("large_twitch_features.csv")
twitch_features <- twitch_features |> 
  rename(node = numeric_id) |> 
  mutate(node = node + 1)

# clean to 10983 nodes
cutoff_views <- quantile(twitch_features$views, 0.5)
cutoff_life_time <- quantile(twitch_features$life_time, 0.5)

# remove EN speakers and dead accounts
# and nodes with views < cutoff views and life_time < cutoff life_time
twitch_sub <- twitch_features |> 
  filter(language != 'EN') |> 
  filter(
    dead_account == 0 & 
      views >= cutoff_views & 
      life_time >= cutoff_life_time
  )

nrow(twitch_sub)
table(twitch_sub$language)

twitch_nodes <- twitch_sub$node
twitch_el_sub <- twitch_el |> 
  filter(node1 %in% twitch_nodes & node2 %in% twitch_nodes)

# largest connected component

twitch.igraph <- graph_from_data_frame(twitch_el_sub, 
                                       directed = FALSE, vertices = twitch_nodes)
twitch.comp <- components(twitch.igraph)
largest.comp <- which.max(twitch.comp$csize)
node.id <- V(twitch.igraph)[twitch.comp$membership == largest.comp]$name

twitch.igraph.lcc <- induced_subgraph(twitch.igraph, node.id)

twitch.adj.lcc <- as_adjacency_matrix(twitch.igraph.lcc, sparse = TRUE)
twitch.el.lcc <- as_data_frame(twitch.igraph.lcc, what = "edges")

twitch.lang <- twitch_sub |> 
  filter(node %in% as.integer(node.id)) |> 
  pull(language)

twitch.comm <- factor(twitch.lang) |> 
  as.integer()

twitch.lang.family <- case_when(
  twitch.lang %in% c('DA', 'DE', 'FI', 'NL', 'NO', 'SV') ~ 'Germanic',
  twitch.lang %in% c('FR', 'IT') ~ 'French',
  twitch.lang %in% c('ES', 'PT') ~ 'Spanish',
  twitch.lang %in% c('JA', 'KO', 'TH', 'TR', 'ZH') ~ 'Asian',
  .default = 'E.EU'
)

twitch.comm.family <- factor(twitch.lang.family) |> 
  as.integer()

saveRDS(list(adjacency = twitch.adj.lcc, classes = twitch.comm.family), "twitch.RData")
