#' Plot a metabolic network
#'
#' @param x A list of class `lim` containing a component `Flowmatrix`.
#' @param seed The seed for the PRNG of the R session
#' @param ... Additionnal arguments to be passed to the function \code{plot.igraph()} from the package \code{\\{igraph\\}}. 
#'
#' @return Returns NULL, invisibly.
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph plot.igraph
#' @method plot lim
#' @export
#'
#' @examples
#' DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
#' BOWF <- df2lim(DF)
#' plot(BOWF)
plot.lim <- function(x, seed=NULL, ...) {
  # Check if 'Flowmatrix' exists in the object 'lim'
  if (!exists("Flowmatrix", where = x)) {
    stop("'Flowmatrix' does not exist in used 'lim' object.")
  }
  
  # Check if 'Flowmatrix' is a matrix
  if (!is.matrix(x$Flowmatrix)) {
    stop("'Flowmatrix' from 'lim' object is not a matrix. ")
  }
  
  # Check if the dimensions of the matrix are equal to the expected dimensions
  # As dimensions depends of NComponents and NExternal we first check if they exist and are integers
  if (!exists("NComponents", where = x)){
    stop("'NComponents' does not exist in used 'lim' object.")
  }
  
  if (!is.integer(x$NComponents)){
    stop("'NComponents' from 'lim' object must be an integer.")
  }
  
  if (!exists("NExternal", where = x)){
    stop("'NExternal' does not exist in used 'lim' object.")
  }
  
  if (!is.integer(x$NExternal)){
    stop("'NExternal' from 'lim' object must be an integer.")
  }
  
  # Calculate the number of nodes
  N_nodes <- x$NComponents + x$NExternal
  
  expected_dim <- c(N_nodes, N_nodes)
  if (all(dim(x$Flowmatrix) != expected_dim)) {
    stop("Flowmatrix from 'lim' object does not match the expected dimension.")
  }

# Set the random seed is specified
  if (!is.null(seed)) {
    set.seed(seed) 
  }

  adjacency_matrix <- as.matrix(x$Flowmatrix > 0)
  model_graph <- graph_from_adjacency_matrix(adjacency_matrix)
  plot.igraph(model_graph, ...)

}
