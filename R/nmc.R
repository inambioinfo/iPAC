nmc <-
function(x,alpha=0.05,multtest=c("Bonferroni","BH","None")){
  
  N <- dim(x)[2] #Columns of the matrix
  n <- length(which(x==1)) #gets the total mutation
  newx <- apply(x,2,sum) # sums down over columns
  mp <- which(newx != 0) # gets the columns which are not equal to 0 - eg those that have mutations
  mp2 <- cumsum(newx[mp]) #sums across the columns
  pvalues2 <- matrix(0,length(mp),2)      #same position
  temp <- difforder0(N,sum(newx),1,mp2[1])
  pvalues2[1,1] <- mp2[1]
  pvalues2[1,2] <- sum(temp)
  if(length(mp) > 1){
    for(i in 2:length(mp)){
      if((mp2[i-1]+1) != mp2[i]){
        temp <- difforder0(N,sum(newx),(mp2[i-1]+1),mp2[i])         }
      else {
        temp <- 1        #if one position has one number
      }
      pvalues2[i,1] <- newx[mp[i]]
      pvalues2[i,2] <- sum(temp)
    }
  }
  pvalues2 <- cbind(mp,mp,pvalues2)
  
  multtest <- match.arg(multtest)
  
  #Bonferroni correction
  if(multtest == "Bonferroni"){
    criterion <- alpha/(length(mp) + (length(mp)*(length(mp)-1)/2))
  }
  else{
    criterion <- alpha
  }
  
  if(n > 1){
    if(length(mp) > 1){
      pvalues <- matrix(0,(length(mp)*(length(mp)-1)/2),2)
      tmp <- rep(mp[1],length(mp)-1)
      tmp2 <- mp[-1]
      tmp3 <- rep(1,length(mp)-1)
      tmp4 <- mp2[-1]
      for(i in 2:length(mp)){
        tmp <- c(tmp,rep(mp[i],length(mp)-i))
        tmp2 <- c(tmp2,mp[-(1:i)])
        tmp3 <- c(tmp3,rep(mp2[i-1]+1,length(mp)-i))
        tmp4 <- c(tmp4,mp2[-(1:i)])
      }
      for(i in 1:dim(pvalues)[1]){
        if(((tmp2[i]-tmp[i]) == (max(tmp2)-min(tmp)))){
          temp <- rangestat(sum(newx),N)
          temp <- temp$cdf[(max(tmp2)-min(tmp))+1]
          pvalues[i,1] <- sum(newx[tmp[i]:tmp2[i]])
          pvalues[i,2] <- sum(temp)
        }
        else if((tmp2[i]-tmp[i]) == 1){
          temp <- difforder1(N,sum(newx),tmp3[i],tmp4[i])
          if((sum(temp) != 0) & (sum(temp) <= criterion)){
            temp2 <- difforder0(N,sum(newx),tmp3[i],tmp4[i])
            temp <- c(temp,temp2)
          }
          pvalues[i,1] <- sum(newx[tmp[i]:tmp2[i]])
          pvalues[i,2] <- sum(temp)
        }
        
        else{
          pvalues[i,1] <- sum(newx[tmp[i]:tmp2[i]])
          pvalues[i,2] <- pbeta(q=abs(tmp[i]-tmp2[i])/N, shape1=tmp4[i]-tmp3[i], shape2=n-tmp4[i]+tmp3[i]+1)
        }
      }
      pvalues <- cbind(tmp,tmp2,pvalues)
      pvalues <- rbind(pvalues,pvalues2)
      pvalues <- pvalues[order(pvalues[,4]),]
      pvalues <- cbind((pvalues[,2]-pvalues[,1]+1),pvalues)
      
      colnames(pvalues) <- c("cluster_size","start","end","number","p_value")
      
      if(multtest == "BH"){
        ap <- mt.rawp2adjp(pvalues[,5],"BH")
        pvalues[,5] <- ap$adjp[order(ap$index),2]
        
        if(length(which(as.numeric(pvalues[,5]) <= criterion)) >= 1){
          pvalues <- pvalues[which(as.numeric(pvalues[,5]) <= criterion),]       
        }
        else {
          pvalues <- NULL
        }
        pvalues
      }
      else{
        if(length(which(as.numeric(pvalues[,5]) <= criterion)) >= 1){
          pvalues <- pvalues[which(as.numeric(pvalues[,5]) <= criterion),]       
        }
        else {
          pvalues <- NULL
        }
        pvalues
      }
    }
    else{
      pvalues <- cbind((pvalues2[,2]-pvalues2[,1]+1),pvalues2)
      colnames(pvalues) <- c("cluster_size","start","end","number","p_value")
      
      if(multtest == "BH"){
        ap <- mt.rawp2adjp(pvalues[,5],"BH")
        pvalues[,5] <- ap$adjp[order(ap$index),2]
        
        if(length(which(as.numeric(pvalues[,5]) <= criterion)) >= 1){
          pvalues <- pvalues[which(as.numeric(pvalues[,5]) <= criterion),]       
        }
        else {
          pvalues <- NULL
        }
        pvalues
      }
      else{
        if(length(which(as.numeric(pvalues[,5]) <= criterion)) >= 1){
          pvalues <- pvalues[which(as.numeric(pvalues[,5]) <= criterion),]      
        }
        else {
          pvalues <- NULL
        }
        pvalues
      }
    }
  }
  else{
    pvalues <- cbind((pvalues2[,2]-pvalues2[,1]+1),pvalues2)
    colnames(pvalues) <- c("cluster_size","start","end","number","p_value")
    
    if(multtest == "BH"){
      ap <- mt.rawp2adjp(pvalues[,5],"BH")
      pvalues[,5] <- ap$adjp[order(ap$index),2]
      
      if(length(which(as.numeric(pvalues[,5]) <= criterion)) >= 1){
        pvalues <- pvalues[which(as.numeric(pvalues[,5]) <= criterion),]       
      }
      else {
        pvalues <- NULL
      }
      pvalues
    }
    else{
      if(length(which(as.numeric(pvalues[,5]) <= criterion)) >= 1){
        pvalues <- pvalues[which(as.numeric(pvalues[,5]) <= criterion),]      
      }
      else {
        pvalues <- NULL
      }
      pvalues
    }
  }
  
}
