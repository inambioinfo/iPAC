difforder1 <-
function(N,n,i,j,r=1){
  diffs <- matrix(0,1,N-r)
  Mp <- c(1/N,1/N,1-2/N)
  for(u in 0:(n-j)){
    for(w in 0:(j-i-1)){
      mr <- c(i+w,j-i-w+u,n-j-u)
      diffs[,1] <- dmultinom(mr,prob=Mp)+ diffs[,1]
    }
  }
  if(N-r > 1){
    for(x in 2:(N-r)){
      Mp <- c((x-1)/N,1/N,1/N,1-(x+1)/N)
      for(u in 0:(i-1)){
        for(w in 0:(j-i-1)){
          for(qs in 0:(n-j)){
            mr <- c(i-1-u,u+1+w,j-i-w+qs,n-j-qs)
            diffs[,x] <- diffs[,x] + dmultinom(mr,prob=Mp)
          }
        }
      }
    }
  }
  diffs
}
