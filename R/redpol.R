
#' @importFrom lsei lsei
#' @importFrom MASS Null


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



red2full<- function(sample,x0,Z){
  res<-x0+Z%*%t(sample)
  x<-t(res)
  
  return(x)
}

full2red<- function(sample,x0,Z){
  res<-solve(t(Z)%*%Z)%*%t(Z)%*%(t(sample)-x0)
  return(t(res))

}
