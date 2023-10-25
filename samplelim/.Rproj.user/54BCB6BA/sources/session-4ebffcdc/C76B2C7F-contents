
#' @importFrom lsei lsei
#' @importFrom MASS Null


full_dim_poly <-function(A,B,G,H,test=TRUE){
  tol=sqrt(.Machine$double.eps) #the smallest positive floating-point number x such that 1 + x != 1 on the current machine
  ## 0. Setup problem
  
  if (is.null(A)){
    stop("no equalities found")
  }
  
  NUnknowns<-ncol(A)
  
  ## additional checks for equalities, hidden in inequalities... (Karline S.)
  if (test && !is.null(G))   {
    xr <- poly_ranges(A,B,G,H)
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

lim_full_dim_poly<-function(lim,test=TRUE){
  return(full_dim_poly(A=lim$A,B=lim$B,G=lim$G,H=lim$H,test=test))
}

transform_sample_from_redspace <- function(sample,x0,Z){
  xnames=colnames(sample)
  res<-x0+Z%*%t(as.matrix(sample))
  x<-t(res)
  colnames (x) <- xnames
  
  return(x)
}