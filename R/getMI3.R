getMI3 <- function(x, y,z , bin=6, so=3, rankBased=FALSE, normalize=TRUE, negateMI = FALSE){
  n <- length(x)
  if(length(y) != n){stop("legnth of two vectors are different!")}
  if(so >= bin){stop("spline order must be less than bin")}
 
  out <- .C("mi3", x = as.double(x), y = as.double(y), z = as.double(z), n = as.integer(n), bin = as.integer(bin), so = as.integer(so), miOut = as.double(0) )
  mi <- out$miOut
  return (mi)
  
}
