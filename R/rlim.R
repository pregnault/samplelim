#' Uniformly sampling high-dimensional polytopes
#'
#' The function \code{rlim()} provides access to implementations of two MCMC methods 
#' for uniformly sampling high-dimensional polytopes: the so-called mirror and 
#' billard walks (MiW and BiW).
#'
#' @param lim A list with four components \code{A}, \code{B}, \code{G} and \code{H} representing
#' the polytope to be sampled; 
#' see the section \emph{Details} below.
#' @param Hpol A Hpolytope object representing the reduced polytope to be sampled. If \code{NULL}
#' (the default), it is computed from \env{lim}.
#' \env{Hpol} is used only if \env{lim} is set to \code{NULL}.
#' @param type A character string specifying the MCMC algorithm to be used for sampling, 
#' @param jmp A numeric value (if \env{type} is \code{"BiW"}) or a numeric vector (if \env{type} is
#' \code{"MiW"}) giving the jump rate parameter of the MCMC algorithm to be used. 
#' If \code{NULL} (the default), it is automatically computed from \env{lim}; 
#' see the section \emph{Details} below.
#' @param thin An integer giving the thinning parameter for the MCMC algorithm.
#' Default is 1 (no thinning).
#' @param scale A numeric value. See the section \emph{Details} below for details. 
#' @param burn An integer giving the burning period for the MCMC algorithm. 
#' Default is 0 (no burn-in).
#' @param nsamp An integer; the length of the output sample of points into the polytope. 
#' Default is 3000.
#' @param nsim An integer value giving the number of sampled points in the polytope. 
#' If not \code{NULL} (the default value), then \env{nsamp} is set to \code{(nsim-burn)/thin}.
# 'whether \code{"MiW"} for the Mirror Walk or \code{"BiW"} for the Billard Walk. 
#' See \emph{Details} Section below.
#' @param starting_point A numeric vector giving the coordinates of a point inside the polytope, 
#' used as starting point in the MCMC algorithm.
#' @param tol A numeric value specifying the tolerance for numeric computations.
#' @param seed An integer used to set the seed of the PRNG.
#'
#' @return A \code{nsamp}*\code{p} matrix whose rows are the coordinates of the \code{nsamp} points
#' sampled in the polytope.
#' The dimension \code{p} is:  
#' \itemize{
#'   \item the dimension of the inner space the polytope relies into (i.e., the number of flows), 
#'   if Hpol is NULL;
#'   \item the dimension of Hpol is Hpol is an object of class Hpolytope.
#'  }
#' @export
#'
#' @details
#' A polytope \eqn{\mathcal{P}} is a convex subset or \eqn{\mathbb{R}^n}, \eqn{n \in \mathbb{N^*}},
#' defined as the intersection of hyper-planes and half-spaces. 
#' More precisely, \eqn{\mathcal{P} = \{ x \in \mathbb{R}^n: Ax = B, Gx \geq H \}}, 
#' where \eqn{A} is an \eqn{m\times n} matrix, with \eqn{m \leq n}, \eqn{B \in \mathbb{R}^m},
#' \eqn{G} is a \eqn{k \times n} matrix and \eqn{H \in \mathbb{R}^k}.
#' For polytopes involved in Linear Inverse Models (LIM), \eqn{n} is 
#' the number of flows in the metabolic network. 
#' The \eqn{k} inequality constraints insures that \eqn{\mathcal{P}} is bounded. 
#' The polytope to be sampled is given to \code{rlim()} whether by:
#' \itemize{
#'   \item creating a list with four elements A, B, G, H, giving the set of equality
#'   and inequality constraints for the polytope, and using it for \code{lim} argument. 
#'   This list can be created manually or by importing a LIM declaration file thanks to 
#'   the function \code{\link{df2lim}()}. Or by
#'   \item creating an object of class Hypolytope (originating from package \pkg{volesti}) 
#'   thanks to the function \code{Hpolytope()}; see \link{Hpolytope-class}. 
#' This Hypolytope object is passed to argument \code{Hpol}, and requires argument \code{lim} 
#' to be set to \code{NULL} for being processed. 
#' }
#' 
#' The main practical difference between lim and Hpolytope objects is that Hypolytope allows 
#' for inequality constraints only, hence the polytope is of full dimension \eqn{n}. 
#'  
#' 
#' 
#' Two MCMC methods for sampling a polytope are implemented. 
#' The MiW has been introduced by \cite{Van Oevelen et al. (2010)}, 
#' while the BiW has been introduced by \cite{Polyak and Gryazina (2014)}.
#' Both are reflective Hamiltonian algorithms: each new point of the sample 
#' is generated from the previous one by choosing randomly a direction and moving forward 
#' this direction for a random length. 
#' If the end point of this segment belongs to the polytope, then it is added to the sample; 
#' if not, then, the trajectory is reflected on borders of polytope until 
#' the end point belongs to the polytope.
#' BiW and MiW differ in the probability distributions used for generating 
#' the length of the trajectory; see \cite{Girardin at al. (2023)} for details.
#' Still, these distributions depend on a parameter called (in both cases) 
#' the jump rate of the trajectory.
#' 
#' The jump rate has a strong influence on the performance of the algorithm: 
#' too small, a very large number of points \env{nsim} would be needed to 
#' efficiently sample the polytope. 
#' Too large, the trajectories reflect an important number of times before resulting 
#' in a point inside the polytope, making the process very slow.
#' If \env{jmp} is \code{NULL} (default), it is computed automatically, 
#' depending the \env{type} (MiW or BiW). Precisely,
#' \itemize{
#'  \item if \env{type} is \code{"MiW"}, then \env{jmp} is a random vector of length 
#'  equal to the dimension of the reduced polytope. 
#'  It is equal to the ranges of the reduced polytope divided by \env{scale} 
#'  (default for \env{scale} is 10);
#'  \item if \env{type} is \code{"BiW"}, the \env{jmp} is a single value 
#'  equal to the radius of the largest ball included in the reduced polytope.
#' }
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
#' # Samples can be thinned and/or burnt-in 
#' sample1 <- rlim(lim= BOWF, nsamp = 20, burn = 50, thin =10, seed = 123)
#' # Default behaviour is to set the size nsamp of the returned sample, 
#' # taking into account burn-in and thinning.
#' # Alternatively, the pilot sample size (i.e., the number of points actually
#' # sampled before burn-in and thinning) nsim can be set
#' sample2 <- rlim(lim = BOWF, nsim = 250, burn = 50, thin = 10, seed = 123)
#' all(sample1 == sample2)
#' @seealso \code{\link{df2lim}()} for creating an lim object from a declaration file, 
#' \link{Hpolytope-class} for class Hpolytope.
#' @references {
#' D. Van Oevelen, K. Van den Meersche, F. J. R. Meysman, K. Soetaert, 
#' J. J. Middelburg and A. F. Vézina, 
#' \emph{Quantifying Food Web Flows Using Linear Inverse Models},
#'  Ecosystems \strong{13}, 32-45 (2010).
#' 
#' B.T. Polyak and E.N. Gryazina, 
#' \emph{Billiard walk - a new sampling algorithm for control and optimization}, 
#' IFAC Proceedings Volumes, \strong{47(3)}, 6123-6128 (2014).
#' 
#' V. Girardin, T. Grente, N. Niquil and P. Regnault,
#' \emph{Comparing and updating R packages of MCMC Algorithms for 
#' Linear Inverse Modeling of Metabolic Networks},
#' hal:    (2023)
#' }
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
  if (!is.null(Hpol) && inherits(Hpol, "Hpolytope")) {
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
      if (!is.numeric(jmp) | length(jmp) != 1) {
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
  if (!is.null(Hpol) && inherits(Hpol, "Hpolytope")) {
    x <- t(res_redspace)
  }
  return(x)
  
}
