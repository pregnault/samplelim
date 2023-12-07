
#' @importFrom Rglpk Rglpk_solve_LP

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


lim.ranges<- function(lim){
  return(pol.ranges(A=lim$A,B=lim$B,G=lim$G,H=lim$H))
}

