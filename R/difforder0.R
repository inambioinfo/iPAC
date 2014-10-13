difforder0 <-
function(N,n,i,j){
  diffs <- matrix(0,1,N)
  diffs[,1] <- 1-pbinom(j-1,n,prob=1/N)
  diffs[,N] <- pbinom(i-1,n,prob=(N-1)/N)
  for(x in 2:(N-1)){
    Mp <- c((x-1)/N,1/N,1-x/N)
    for(w in 0:(i-1)){
      for(u in 0:(n-j)){
        mr <- c(w,n-w-u,u)
        diffs[,x] <- diffs[,x] + dmultinom(mr,prob=Mp)
      }
    }
  }
  diffs
}
