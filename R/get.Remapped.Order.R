get.Remapped.Order <- function(mutation.data, position.data,method = "MDS", 
                                OriginX = min(position.data[,4]), OriginY=min(position.data[,5]), OriginZ=min(position.data[,6])){
  
  #mutation.data.culled removes the columns from the mutation matrix that do not correspond to positional entries
  mutation.data.culled<-mutation.data[,position.data[,2]]
  
  if(method == "Linear"){  
    #Gets the Remapped Protein and Plots if required
    remapped<-LinearRemapperMain(position.data,OriginX, OriginY, OriginZ, mutation.data.culled)
    aa.pos <- as.numeric(substring(colnames(remapped[[1]]), 2))
    remapped.culled.order <-as.numeric(aa.pos[remapped[[2]]])
  }else if(method == "MDS"){
    remapped<-MDSRemapperMain(position.data,mutation.data.culled)
    aa.pos <- as.numeric(substring(colnames(remapped[[1]]), 2))
    remapped.culled.order <-as.numeric(aa.pos[remapped[[2]]])
  }
  return(remapped.culled.order)
}