rangestat <-
function(n,N){
  
  rangescdf <- matrix(0,1,N)
  rangespdf <- matrix(0,1,N)
  rangescdf[,1] <- rangespdf[,1] <- (1/N)^(n-1)
  for(i in 1:(N-1)){
    rangescdf[,(i+1)] <- rangescdf[,i] + (N-i)*(((i+1)/N)^n-2*(i/N)^n+((i-1)/N)^n)
    rangespdf[,(i+1)] <- (N-i)*(((i+1)/N)^n-2*(i/N)^n+((i-1)/N)^n)
  }
  colnames(rangescdf) <- 0:(N-1)
  colnames(rangespdf) <- 0:(N-1)
  ranges <- list(rangescdf,rangespdf)
  names(ranges) <- c("cdf","pdf")
  ranges
}
