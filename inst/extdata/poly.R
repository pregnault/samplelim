library("LIM")
library("SampleLIM")
library("GGally")
library(MASS)
library("tidyverse")


get_Poly_from_LIM <-function(lim){
  tol=sqrt(.Machine$double.eps) #the smallest positive floating-point number x such that 1 + x != 1 on the current machine
  ## 0. Setup problem
  
  Ncomp  <- lim$NComponents
  Nx     <- lim$NUnknowns
  
  #Equalities A*X=B
  A    <- lim$A[]
  B    <- lim$B[]
  
  # Inequalities G*X>=H
  G <- lim$G
  H <- lim$H
  
  
  ## find a particular solution x0
  
  l <- lsei(A=NULL,B=NULL,E=A,F=B,G=G,H=H)
  if (l$residualNorm>1e-6)
    stop("no particular solution found;incompatible constraints")
  else
    x0 <- l$X
  lx <- length(x0)
  
  
  Z <- Null(t(A)); Z[abs(Z)<tol] <- 0  #x=x0+Zq ; AZ=0
  
  k <- ncol(Z)
  
  #Projection de G et H sur l'espace réduit
  g <- G%*%Z
  h <- H-G%*%x0                                            
  g[abs(g)<tol] <- 0
  h[abs(h)<tol] <- 0
  h<-as.numeric(h)
  
  return(list("A"=-g,"b"=-h,"x0"=x0,"Z"=Z))
  
}

transform_sample_from_redspace <- function(sample,x0,Z){
  
  return(x0+Z%*%as.matrix(sample))
}

automatedjump <- function(g,h,g.scale=5){
  k=ncol(g)
  if (is.null(g)) s <- rep(NA,k)
  else  {
    q_ranges <- xranges(E=NULL,F=NULL,-g,-h)
    s <- abs(q_ranges[,1]-q_ranges[,2])/g.scale
    s[which(s>1e20)]<-1e-2
  }
  if (any (is.na(s)))  {
    if (all(is.na(s)))  {
      warning(" problem is unbounded - all jump lengths are set to 1")
      s[] <- 1
    } else {
      warning(" problem is unbounded - some jump lengths are set arbitrarily")
      s[is.na(s)] <- mean(s,na.rm=T)*100
    }
  }
  return(s)
}


# Construction de l'objet LIM à partir du declaration file

filename="DeclarationFileBOWF.txt"
web.lim<- Setup(filename)


# Utilisation de limsolve pour générer des échantillons /!\ Prend du temps sur des gros problèmes
res_xsample<-Xsample(web.lim,iter=10,jmp=0.05)

# Construction du polytope de solutions sur l'espace réduit, renvoie un objet contenant la matrice A et le vecteur b
# définissant le polytope Ax<b. Ainsi que la matrice de passage Z et une solution particulière x0. 
# Pour tout point q de l'espace réduit, le point x correspondant dans l'espace de départ est x=x0+Zq

pol<-get_Poly_from_LIM(web.lim)

#Construit un vecteur de variances à partir des ranges du polytope /!\ Pas fiable sur les gros problèmes, à investiguer

jump<-automatedjump(g=pol$A,h=pol$b,g.scale = 5)

#Construit le polytope avec pour équation Ax<b

P <- Hpolytope(A = pol$A, b = pol$b)


#Sampling des points dans le polytope réduit avec volesti mirroir

res_redspace<-sample_points(P,n=100,random_walk = list("walk"="mirror","jump"=0.05))

#Transformation des solutions dans l'espace réduit en solutions de l'espace de départ 

res<-transform_sample_from_redspace(res_redspace,pol$x0,pol$Z)

# Construction du dataset correspondant aux solutions pour les afficher

rownames(res) = web.lim$Unknowns
res<-as_tibble(t(res))

# Affichage des solutions projetées sur 2 fluxs avec leur covariance

ggpairs(res,columns=5:7, lower=list(continuous=wrap("points",alpha=0.3,size=0.1),combo=wrap("dot",alpha=0.4,size=0.2)),title =paste(filename,"Mirroir C"))
