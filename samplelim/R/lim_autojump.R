
lim_autojump <- function(lim,scale=10)   {
  full_pol<-full_dim_poly(lim$A,lim$B,lim$G,lim$H)
  ranges<-poly_ranges(G=full_pol$G,H=full_pol$H)[,3]
  return(ranges/scale)
}