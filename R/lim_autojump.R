
lim.autojump <- function(lim,scale=10)   {
  full_pol<-redpol(lim$A,lim$B,lim$G,lim$H)
  ranges<-pol.ranges(G=full_pol$G,H=full_pol$H)[,3]
  return(ranges/scale)
}