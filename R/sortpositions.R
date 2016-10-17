sortpositions <-
function(ordering,mutation.data){
  protein.ordering<- order(ordering) # gets the new protein ordering
  remapped.mutation.data <- mutation.data #instantialize sorted protein data
  
  #sorts the protein to the new allocation
  for(i in 1:dim(mutation.data)[2]){ 
    remapped.mutation.data[,i]<- mutation.data[,protein.ordering[i]]
  }
  
  #this returns a list of 3 elements. The first element is the reordered matrix of mutation data. The second is the numerical ordering of the proteins. The third
  #is the original ordering passed in.
  return(list(remapped.mutation.data, protein.ordering, ordering))
}
