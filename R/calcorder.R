calcorder <-
function(positions){
  
  #creates the dist object and then runs the MDS scale command
  distance.matrix<- dist(positions[,4:6], method ="euclidean", p  = 2)
  distance.vector <- cmdscale(distance.matrix, k = 1, eig = FALSE, add = FALSE, x.ret= FALSE)
  return(distance.vector)
}
