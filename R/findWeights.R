getWeights <- function(x, bin=6, so=3){
  n <- length(x)
  if(so >= bin){stop("spline order must be less than bin")}

    out <- .C("findWeightsR", x = as.double(x), n=as.integer(n), bin=as.integer(bin), so = as.integer(so), eOut = matrix(0, n, bin))

  e <- out$eOut
  return (e)

}