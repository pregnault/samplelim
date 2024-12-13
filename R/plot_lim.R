plot.lim <- function(lim, seed=NULL, ...) {
  # Check if 'Flowmatrix' exists in the object 'lim'
  if (!exists("Flowmatrix", where = lim)) {
    stop("'Flowmatrix' does not exist in used 'lim' object.")
  }
  
  # Check if 'Flowmatrix' is a matrix
  if (!is.matrix(lim$Flowmatrix)) {
    stop("'Flowmatrix' from 'lim' object is not a matrix. ")
  }
  
  # Check if the dimensions of the matrix are equal to the expected dimensions
  # As dimensions depends of NComponents and NExternal we first check if they exist and are integers
  if (!exists("NComponents", where = lim)){
    stop("'NComponents' does not exist in used 'lim' object.")
  }
  
  if (!is.integer(lim$NComponents)){
    stop("'NComponents' from 'lim' object must be an integer.")
  }
  
  if (!exists("NExternal", where = lim)){
    stop("'NExternal' does not exist in used 'lim' object.")
  }
  
  if (!is.integer(lim$NExternal)){
    stop("'NExternal' from 'lim' object must be an integer.")
  }
  
  # Calculate the number of nodes
  N_nodes <- lim$NComponents + lim$NExternal
  
  expected_dim <- c(N_nodes, N_nodes)
  if (all(dim(lim$Flowmatrix) != expected_dim)) {
    stop("Flowmatrix from 'lim' object does not match the expected dimension.")
  }

# Set the random seed is specified
  if (!is.null(seed)) {
    set.seed(seed) 
  }

  adjacency_matrix <- as.matrix(lim$Flowmatrix > 0)
  model_graph <- graph_from_adjacency_matrix(adjacency_matrix)
  plot.igraph(model_graph, ... )

}
