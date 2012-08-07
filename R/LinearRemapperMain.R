LinearRemapperMain <-
function(positions, OriginX, OriginY, OriginZ, mutation.data){
  distance.from.origin <- calcdistance(positions, OriginX, OriginY, OriginZ)   
  remapping<- sortpositions(distance.from.origin, mutation.data)
  return(remapping)   
}
