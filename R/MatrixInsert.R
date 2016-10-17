MatrixInsert <-
function(missing.values.matrix, results.matrix){
    missing.values.matrix <- cbind(missing.values.matrix, missing.values.matrix[,2]-missing.values.matrix[,1]+1)
    new.results <- rbind(results.matrix)
    
    if(dim(missing.values.matrix)[1]!=0){
      for(i in 1 : dim(new.results)[1]){
        for(j in 1: dim(missing.values.matrix)[1]){
          if(new.results[i,1]>= missing.values.matrix[j,1]){
            new.results[i,1]<- new.results[i,1]+missing.values.matrix[j,3]
          }
          #For the right hand side there is a second condition to check if you end up in a missing section as well. 
          if(new.results[i,2]>=missing.values.matrix[j,2] || ( (new.results[i,2]>=missing.values.matrix[j,1]) && (new.results[i,2] <= missing.values.matrix[j,2]) )){
            new.results[i,2] <- new.results[i,2]+missing.values.matrix[j,3]
          }
        }
      }
    }
    return(new.results)
  }
