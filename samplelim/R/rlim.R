
rpol <- function(A=NULL,B=NULL,G=NULL,H=NULL, walk_length=NULL, nburns=NULL, iter=3000, type="mirror", jmp=NULL, starting_point=NULL,seed=NULL){
  
  automatedjump <- function(G,H,scale=10)   {
    ranges<-pol.ranges(G=G,H=H)[,3]
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
    full_pol <-redpol(A,B,G,H)
    x0<-full_pol$x0
    Z<-full_pol$Z
    g=full_pol$G
    h=as.numeric(full_pol$H)
      
    }
  
  P <- Hpolytope(A = -g, b = -h)
  #exploration
  
  if (type=="BiW"){
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
  if (!is.null(walk_length)){
    random_walk<-c(random_walk,list("walk_length"=walk_length))
  }
  if (!is.null(nburns)){
    random_walk<-c(random_walk,list("nburns"=nburns))
  }
  if (!is.null(starting_point)){
    random_walk<-c(random_walk,list("starting_point"=starting_point))
  }
  if (!is.null(seed)){
    random_walk<-c(random_walk,list("seed"=seed))
  }
  
  res_redspace<-as.matrix(sample_points(P,n=iter,random_walk = random_walk))
  
  res<-x0+Z%*%res_redspace
  x<-t(res)
  
  
  xnames <- colnames(A)
  if (is.null(xnames)) xnames <- colnames(G)
  colnames (x) <- xnames
  return(x)
  
}


rlim<- function(lim, walk_length=NULL, nburns=NULL, iter=3000, type="mirror", jmp=NULL,
                tol=sqrt(.Machine$double.eps), starting_point=NULL,seed=NULL){
return(rpol(A=lim$A,B=lim$B,G=lim$G,H=lim$H, walk_length=walk_length, nburns=nburns, iter=iter, type=type, jmp=jmp, starting_point=starting_point,seed=seed))   }
