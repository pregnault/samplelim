
#' Imports Declaration File as an lim object
#'
#' The function \code{df2lim()} turns a declaration file with linear inverse model specifications and creates the corresponding lim object.
#'
#' @param filename A character string specifying the path to the declaration file to import.
#'
#' @return An lim object containing, among other components, the matrices A, B, G and H
#' that characterize the equality and inequality constraints defining the polytope 
#' associated to the LIM.
#' @export
#'
#' @examples
#' # Create a lim object from a Description file
#' DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
#' BOWF <- df2lim(DF)
df2lim<-function(filename){
  return(Setup.liminput(Read(filename)))
}