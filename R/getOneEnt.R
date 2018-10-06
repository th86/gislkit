getOneEnt <- function(x, bin=6, so=3){
  n <- length(x)
  if(so >= bin){stop("spline order must be less than bin")}

     out <- .C("entropy1R", x = as.double(x), n=as.integer(n), bin=as.integer(bin), so = as.integer(so), eOut = 0)

  e <- out$eOut
  return (e)

}