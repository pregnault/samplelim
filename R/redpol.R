#' Projection of full polytope into the reduced polytope
#'
#' The function \code{lim.redpol()} takes as input a polytope  and returns its projection into the non-empty reduced polytope.
#' Precisely, taking the polytope \eqn{\mathcal{P}= \{ x \in \mathbb{R}^n: Ax = B, Gx \geq H \}} as input, the function returns 
#' \itemize{
#' \item the matrix \code{Z}, basis of the right null space of A and \eqn{x_0} a particular solution of \eqn{P}, used for the reduction;
#' \item the matrix \eqn{G'=GZ} and the vector \eqn{H'=H-Gx_0} describing the reduced polytope \eqn{\mathcal{P'}= \{ x \in \mathbb{R}^{n-k}: G'x \geq H'\}} with \eqn{k=\mathtt{rank}(A)}.
#'  }
#'  
#'  The function \code{full2red()} (resp. \code{red2full()}) turns a sample of points inside the full polytope \eqn{\mathcal{P}} (resp. inside the reduced polytope \eqn{\mathcal{P'}})
#'  into the sample of corresponding points inside the reduced polytope \eqn{\mathcal{P'}} (resp. the full polytope \eqn{\mathcal{P}}).
#'
#'
#' @param lim A list with four components \code{A}, \code{B}, \code{G} and \code{H} representing
#' the polytope to be reduced.
#' @param test A boolean if equal \code{TRUE} checks for additional equalities hidden in inequalities.
#'
#' @return A list with four components; namely:
#' \itemize{
#'   \item \code{G}
#'   \item \code{H}
#'   \item \code{x0}
#'   \item \code{Z}
#' }
#' 
#' @rdname lim.redpol
#' @importFrom lsei lsei
#' @importFrom MASS Null
#' @export
#'
#' @examples
#' DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
#' BOWF <- df2lim(DF)
#' BOWFred <- lim.redpol(BOWF)
#' str(BOWFred, max.length = 1)
lim.redpol <-function(lim,test=TRUE){
  A=lim$A
  B=lim$B
  G=lim$G
  H=lim$H
  tol=sqrt(.Machine$double.eps) #the smallest positive floating-point number x such that 1 + x != 1 on the current machine
  ## 0. Setup problem
  
  if (is.null(A)){
    stop("no equalities found")
  }
  
  NUnknowns<-ncol(A)
  
  ## additional checks for equalities, hidden in inequalities... (Karline S.)
  if (test && !is.null(G))   {
    xr <- pol.ranges(A,B,G,H)
    ii <- which (xr[,3]==0)
    if (length(ii)>0) { # if they exist: add regular equalities !
      dia <- diag(nrow=nrow(xr))
      A  <- rbind(A,dia[ii,])
      B  <- c(B,xr[ii,1])
  }}
  
  
  ## find a particular solution x0
  D<-matrix(data=0, nrow=1,ncol=NUnknowns)
  d<-0
  dvec  <- crossprod(D, d)
  Dmat<-crossprod(D,D)
  diag(Dmat)<-diag(Dmat)+1e-11
  sol <- lsei(a=Dmat,b=dvec,c=A, d=B,e=G,f=H)
  x0<-round(sol,digits=8)
  
  Z <- Null(t(A)); Z[abs(Z)<tol] <- 0  #x=x0+Zq ; AZ=0
  
  k <- ncol(Z)
  
  #Projection de G et H sur l'espace réduit
  g <- G%*%Z
  h <- H-G%*%x0                                            
  g[abs(g)<tol] <- 0
  h[abs(h)<tol] <- 0
  h<-as.numeric(h)
  
  return(list("G"=g,"H"=h,"x0"=x0,"Z"=Z))
  
}


#' @rdname lim.redpol
#' @param sample  A matrix where each row corresponds to a point inside either the full or the reduced polytope.
#' @param x0 A numeric vector of size \eqn{n} corresponding to the particular solution used during the reduction.
#' @param Z The matrix used during the reduction of the polytope (returned by \code{lim.redpol()}).
#' @export

red2full<- function(sample,x0,Z){
  res<-x0+Z%*%t(sample)
  x<-t(res)
  
  return(x)
}

#' @rdname lim.redpol
#' @export

full2red<- function(sample,x0,Z){
  res<-solve(t(Z)%*%Z)%*%t(Z)%*%(t(sample)-x0)
  return(t(res))

}
