#' Uniformly sampling high-dimensional polytopes
#'
#' The function \code{rlim()} provides access to implementations of two MCMC methods 
#' for uniformly sampling high-dimensional polytopes: the so-called mirror and 
#' Billard walks (MiW and BiW).
#'
#' @param lim A list with four components \code{A}, \code{B}, \code{G} and \code{H} representing the full-dimensional polytope to be sampled; see the section \emph{Details} below.
#' @param thin An integer giving the thinning parameter for the MCMC algorithm. Default is 1 (no thinning).
#' @param burn An integer giving the burning period for the MCMC algorithm. Default is 0 (no burn-in).
#' @param nsamp An integer; the length of the output sample of points into the polytope.
#' @param nsamp An integer value giving the number of sampled points in the polytope. 
#' @param type A character string specifying the MCMC algorithm to be used for sampling, 
# 'whether \code{"MiW"} for the Mirror Walk or \code{"BiW"} for the Billard Walk. See \emph{Details} Section below.
#' @param jmp A numeric value, the jump length parameter of the MCMC algorithm to be used. If \code{NULL}, a default value is computed from \env{lim}. See the section \emph{Details} below.
#' @param tol A numeric value specifying the tolerance for numeric computations.
#' @param starting_point A numeric vector giving the coordinates of a point inside the polytope, used as starting point in the MCMC algorithm.
#'
#' @return A \code{nsamp}*\code{p} matrix whose rows are the coordinates of the \code{nsamp} points sampled in the polytope.
#' The dimension \code{p} is:  
#'  * the dimension of the inner space the polytope relies into (i.e., the number of flows), if Hpol is NULL;
#'  * the dimension of Hpol is Hpol is an object of class Hpolytope.
#' @export
#'
#' @details
#' ADD.
#' 
#' @examples
#' # Create a lim object from a Description file
#' DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
#' BOWF <- df2lim(DF)
#' # Manage sampling in this LIM with rlim()
#' rlim(lim = BOWF, nsamp = 20, seed = 123)
#' # Alternatively, you can use a Hpolytope object as input
#' # It allows using volesti input inside rlim
#' A = matrix(c(-1,0,0,-1,1,1), ncol=2, nrow=3, byrow=TRUE)
#' b = c(0,0,1)
#' P = Hpolytope(A = A, b = b)
#' rlim(lim=NULL, Hpol = P, nsamp = 20, seed = 123)
#' # Samples can be thinned and burnt-in 
#' rlim(lim= BOWF, nsamp = 20, burn = 50, thin =10, seed = 123)
rlim<- function(lim, Hpol = NULL,
                type = "MiW", jmp = NULL, scale = 10,
                thin = 1, burn = 0, nsamp = 3000, nsim = NULL,
                starting_point = NULL, 
                tol = sqrt(.Machine$double.eps), 
                seed = NULL){

  automatedjump <- function(G,H,scale)   {
    ranges<-pol.ranges(G=G,H=H)[,3]
    return(ranges/scale)
    
  }
  
  if (!is.null(lim) && !is.null(lim$G) && !is.null(lim$H)){
  # Extracting from lim, the matrices of equalities and inequalities defining the polytope.
  A <- lim$A
  B <- lim$B
  G <- lim$G
  H <- lim$H
  
  ## conversions vectors to matrices and checks
  if (is.data.frame(A)) A <- as.matrix(A)
  if (is.data.frame(G)) G <- as.matrix(G)
  if (is.vector(A)) A <- t(A)
  if (is.vector(G)) G <- t(G)
  
  if(is.null(G)){
    stop("G is null, the polytope is not defined")
  }
  }
  # If Hpol is NULL (default), the reduced polytope is computed, as a Hpolytope object.
  if (is.null(Hpol)){
    if(is.null(A)){
      x0<-0
      Z<-diag(ncol(G))
      g=G
      h=as.numeric(H)
    }else{
      full_pol <-lim.redpol(lim)
      x0<-full_pol$x0
      Z<-full_pol$Z
      g=full_pol$G
      h=as.numeric(full_pol$H)
      
    }
    # P: reduced polytope
    P <- Hpolytope(A = -g, b = -h)
  }
  # If Hpol is a Hpolytope object, then lim is ignored and 
  # sampling is performed directly into Hpol.
  if (!is.null(Hpol) && class(Hpol) == "Hpolytope") {
    warning("Argument lim is ignored as Hpol has been set manually.  
            Sampling is performed directly into Hpol.")
    P <- Hpol
    g <- -Hpol@A
    h <- -Hpol@b
  }
  #exploration
  if (type=="BiW"){
    if (is.null(jmp)) {
      random_walk <- list("walk"="BiW")
    }else{
      if (!is.numeric(jmp) | length(jpm) != 1) {
        stop("jmp has to be a single numeric value for type='BiW'.")
      }
      random_walk <- list("walk"="BiW",L=jmp)
    }
  }
  #mirror
  else if (type=="MiW"){
    if (is.null(jmp)) {
      
      random_walk <- list("walk"="mirror","jump"=automatedjump(g,h,scale=scale))
    }else{
      random_walk <- list("walk"="mirror","jump"=jmp)
    }
  }
  else{
    stop("The walk type is not valid")
    
  }
  if (!is.null(thin)){
    random_walk<-c(random_walk,list("walk_length"=thin))
  }
  if (!is.null(burn)){
    random_walk<-c(random_walk,list("nburns"=burn))
  }
  if (!is.null(starting_point)){
    random_walk<-c(random_walk,list("starting_point"=starting_point))
  }
  if (!is.null(seed)){
    random_walk<-c(random_walk,list("seed"=seed))
  }
  if ((!is.null(nsim)) && is.numeric(nsim) && length(nsim) == 1L) {
    warning("nsamp is ignored as nsim has been set manually.
            nsamp is set to (nsim-burn)/thin.")
    nsamp <- (nsim-burn)/thin
  }
  if (!is.null(nsim) && !is.numeric(nsim)){
    warning("nsim is ignored as it is not numeric.")
  }
  res_redspace<-as.matrix(sample_points(P,n=nsamp,random_walk = random_walk))
  
  if (is.null(Hpol)) {
    res<-x0+Z%*%res_redspace
    x<-t(res)
    
    
    xnames <- colnames(A)
    if (is.null(xnames)) xnames <- colnames(G)
    colnames (x) <- xnames
  }
  if (!is.null(Hpol) && class(Hpol) == "Hpolytope") {
    x <- t(res_redspace)
  }
  return(x)
  
}
