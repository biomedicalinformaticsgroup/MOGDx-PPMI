library(SNFtool)
library(igraph)
library(data.table)
source('~/Year2/MOGDx/R/preprocess_functions.R')

setwd('~/Year2/MOGDx/')

dataset <- 'PPMI'
project <- 'All'
trait <- c('CONCOHORT_DEFINITION')
TimeStep <- 'V08'
index_col <- 'PATNO'

# The list of modalities
modalities <- c( 'mRNA' , 'miRNA' , 'DNAm' , 'SNP' , 'Clinical' )

# Initialize an empty list to store sublists
mod_list <- list()

for (comb_length in 2:length(modalities)) { 
  
  len_mod_list <- length(mod_list)
  # Get all combinations without repetition
  combinations <- combn(modalities, comb_length)
  
  # Convert the matrix of combinations into a list of lists
  for (i in (len_mod_list+1):(len_mod_list + ncol(combinations))) {
    sublist <- c(combinations[, i-len_mod_list])
    mod_list[[i]] <- sublist
  }
}

for (sub_mod_list in mod_list) {
  colnames <- c('PATNO' ,  'race' , 'GENDER' , 'AGE_AT_VISIT' , trait)
  datMeta <- t(data.frame( row.names = colnames))
  for (mod in sub_mod_list) {
    print(mod)
    datMeta <- rbind(datMeta , read.csv(paste0('./data/',dataset,'/raw/',project,'/',TimeStep,'/output/',TimeStep,'/datMeta_',mod,'.csv') , row.names = 1)[ , colnames])
  }
  datMeta <- datMeta[!(duplicated(datMeta)),]
  rownames(datMeta) <- datMeta[[index_col]]
  print(dim(datMeta))
  
  all_idx <- c()
  g_list <- list()
  for (net in list.files('./Network/')) {
    if ((unlist(strsplit(net , '_'))[2] %in% sub_mod_list)&&(unlist(strsplit(net , '_'))[1] == TimeStep)) {
      print(net)
      net_graph <- read.csv(paste0('./Network/SNF/',net) , row.names = 1)
      patients <- unique(data.frame(id = c(net_graph$from_name , net_graph$to_name) ,
                                    class = c(net_graph$from_class , net_graph$to_class)))
      relation <- data.frame(from = net_graph$from_name , 
                             to = net_graph$to_name )
      
      g_net <- graph_from_data_frame(relation , directed = FALSE , vertices = patients)
      g_net <- simplify(g_net, remove.multiple=TRUE, remove.loops=TRUE)
      
      g_list[[net]] <- g_net
      all_idx <- unique(append(all_idx,V(g_net)$name))
    }
  }
  
  # This for loop extracts the adjacency (similarity/affinity) matrix from each graph.
  adjacency_graphs <- list()
  for (graph_names in names(g_list)) {
    
    missing_idx <- setdiff(all_idx , V(g_list[[graph_names]])$name)
    g_list[[graph_names]] <- add_vertices(g_list[[graph_names]] , length(missing_idx) , name = missing_idx)
    
    graph_adj <- as.matrix(as_adjacency_matrix(g_list[[graph_names]]))[all_idx,all_idx]
    
    adjacency_graphs[[graph_names]] <- graph_adj
    
  }
  
  ## First, set all the parameters:
  K = 15;		# number of neighbors, usually (10~30)
  T = 10; 	# Number of Iterations, usually (10~20)
  
  #change this to similarity matrix
  W = SNF(adjacency_graphs, K , T)
  W <- W - diag(0.5 , dim(W)[1]) 
  
  g <- snf.to.graph(W , datMeta , trait , all_idx , sub_mod_list)
  
  print(length(V(g)))
  write.csv(as_long_data_frame(g) , file = paste0('./data/',project,'/raw/',project,'/',TimeStep,'/output/',TimeStep,'/',TimeStep,'_',paste0(sub_mod_list , collapse = '_'),'_graph.csv'))
}

