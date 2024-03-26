

#' Determine ranges of a polytope
#'
#' The functions \code{pol.ranges()} and \code{lim.ranges()} compute the theoretical ranges of the polytope along each dimension of a given polytope \eqn{\mathcal{P}= \{ x \in \mathbb{R}^n: Ax = B, Gx \geq H \}}.  
#'
#' @param A A matrix corresponding to \code{A} in the description of the polytope \eqn{\mathcal{P}}. 
#' @param B A numeric vector corresponding to \code{B} in the description of the polytope \eqn{\mathcal{P}}. 
#' @param G A matrix corresponding to \code{G} in the description of the polytope \eqn{\mathcal{P}}. 
#' @param H A numeric vector corresponding to \code{H} in the description of the polytope \eqn{\mathcal{P}}. 
#'
#' @return A \eqn{n \times 3} matrix, with \eqn{n} the dimension of the polytope. The three columns of the matrix correspond respectively
#' to the minimum, maximum and range along each dimension of the polytope.
#' @importFrom Rglpk Rglpk_solve_LP
#' @export
#'
#' @rdname pol.ranges
#' @examples
#' # Create a lim object from a Description file
#' DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
#' BOWF <- df2lim(DF)
#' pol.ranges(A = BOWF$A, B = BOWF$B, G = BOWF$G, H = BOWF$H)
pol.ranges <- function(A=NULL,B=NULL,G,H)   {

  if (is.data.frame(A)) A <- as.matrix(A)
  if (is.data.frame(G)) G <- as.matrix(G)
  if (is.vector(A)) A <- t(A)
  if (is.vector(G)) G <- t(G)
  
  if (is.null(G)){
    stop("G is NULL, the polytope has 0 dimensions.")
  }
  
  if (nrow(G)!=length(H)){
    stop("G and H have incompatible dimensions.")
  }
  
  
  
  Nconstr<-nrow(G)
  if (is.null(A)){
    Neq=0
  }else{
    Neq<-nrow(A)
  }
  
  if(is.null(B) & Neq !=0){
    stop("A and B have incompatible dimensions.")
  }
  
  if (Neq !=length(B)){
    stop("A and B have incompatible dimensions.")
  }
  
  d<-ncol(G)
  ranges<-matrix(nrow=d,ncol=3)
  colnames(ranges)<-c("min","max","range")
  constraints<-rbind(G,A)
  rhs_vec<-c(H,B)
  constraints_direction<-c(rep(">=",Nconstr),rep("==",Neq))
  bounds <- list(lower = list(ind = c(1:d), val = rep(-Inf,d)),
                 upper = list(ind = c(1:d), val = rep(Inf ,d)))
  for (i in 1:d){
    obj<-rep(0,d)
    obj[i]=1
    ranges[i,1]=Rglpk_solve_LP(obj=obj,mat=constraints,dir=constraints_direction,rhs=rhs_vec,bounds=bounds,max=FALSE)$optimum
    ranges[i,2]=Rglpk_solve_LP(obj=obj,mat=constraints,dir=constraints_direction,rhs=rhs_vec,bounds=bounds,max=TRUE)$optimum
    ranges[i,3]=ranges[i,2]-ranges[i,1]
  }
  return(round(ranges,digits=8))
}


#' @param lim A list with four components \code{A}, \code{B}, \code{G} and \code{H} representing
#' the polytope.
#' @export
#' @rdname pol.ranges
lim.ranges<- function(lim){
  ranges<-pol.ranges(A=lim$A,B=lim$B,G=lim$G,H=lim$H)
  if(!is.null(lim$Unknowns)&&length(lim$Unknowns)==nrow(ranges)){
    rownames(ranges)<-lim$Unknowns
  }
  return(ranges)
}

