
lim.autojump <- function(lim,scale=10)   {
  full_pol<-lim.redpol(lim)
  ranges<-pol.ranges(G=full_pol$G,H=full_pol$H)[,3]
  return(ranges/scale)
}