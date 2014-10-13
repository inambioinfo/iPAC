MDSRemapperMain <-
function(positions, mutation.data){
  ordering <- calcorder(positions)   
  remapping<- sortpositions(ordering,mutation.data)
  return(remapping)   
}
