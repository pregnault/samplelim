
rpol <- function(A=NULL,B=NULL,G=NULL,H=NULL, W=1, iter=3000, type="mirror", jmp=NULL, x0=NULL){
  
  automatedjump <- function(G,H,scale=10)   {
    ranges<-poly_ranges(G=G,H=H)[,3]
    return(ranges/scale)
    
  }
  

  ## conversions vectors to matrices and checks
  if (is.data.frame(A)) A <- as.matrix(A)
  if (is.data.frame(G)) G <- as.matrix(G)
  if (is.vector(A)) A <- t(A)
  if (is.vector(G)) G <- t(G)
  
  if(is.null(G)){
    stop("G is null, the polytope is not defined")
  }
  
  ## full dimensionnal polytope
  
  if(is.null(A)){
    x0<-0
    Z<-diag(ncol(G))
    g=G
    h=as.numeric(H)
  }else{
    full_pol <-full_dim_poly(A,B,G,H)
    x0<-full_pol$x0
    Z<-full_pol$Z
    g=full_pol$G
    h=as.numeric(full_pol$H)
      
    }
  
  P <- Hpolytope(A = -g, b = -h)
  #exploration
  
  if (type=="Biw"){
    if (is.null(jmp)) {
      random_walk <- list("walk"="BiW")
    }else{
      random_walk <- list("walk"="BiW",L=jmp)
    }
  }
  #mirror
  else if (type=="mirror"){
    if (is.null(jmp)) {
      
      random_walk <- list("walk"="mirror","jump"=automatedjump(g,h))
    }else{
      random_walk <- list("walk"="mirror","jump"=jmp)
    }
  }
  else{
    stop("The walk type is not valid")
    
  }
  res_redspace<-as.matrix(sample_points(P,n=W*iter,random_walk = random_walk))
  
  res<-x0+Z%*%res_redspace
  x<-t(res)
  x<-x[c(1:iter)*W,]
  
  
  xnames <- colnames(A)
  if (is.null(xnames)) xnames <- colnames(G)
  colnames (x) <- xnames
  return(x)
  
}


#' Uniformly sampling high-dimensional polytopes
#'
#' The function \code{rlim()} provides implementations of two MCMC methods 
#' for uniformly sampling high-dimensional polytopes: the so-called mirror and 
#' Billard walks (MiW and BiW).
#'
#' @param lim A \code{lim} object representing the polytope to be sampled; see the documentation of \code{\link{lim}} object.
#' @param W An integer giving the thinning parameter for the MCMC algorithm. Default is 1 (no thinning).
#' @param iter An integer value giving the number of sampled points in the polytope. 
#' @param type A character string specifying the MCMC algorithm to be used for sampling, 
# 'whether \code{"MiW"} for the Mirror Walk or \code{"BiW"} for the Billard Walk. See \emph{Details} Section below.
#' @param jmp A numeric value, the jump rate parameter of the MCMC algorithm to be used. If \code{NULL}, a default value is computed from \env{lim}. See \emph{Details} Section below.
#' @param tol A numeric value specifying the tolerance for numeric computations.
#' @param x0 A numeric vector giving the coordinates of a point inside the polytope, used as starting point in the MCMC algorithm.
#'
#' @return A \code{iter}*\code{p} matrix whose rows are the coordinates of the \code{iter} points sampled in the polytope.
#' @export
#'
#' @details
#' ADD.
#'
#' @examples
rlim<- function(lim, W=1L, iter=3000L, type="mirror", jmp=NULL,
                tol=sqrt(.Machine$double.eps), x0=NULL){

return(rpol(A=lim$A,B=lim$B,G=lim$G,H=lim$H, W=W, iter=iter, type=type, jmp=jmp, x0=x0))   }
