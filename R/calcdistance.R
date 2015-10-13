calcdistance <-
function(positions,OriginX, OriginY, OriginZ){
  distance.vector<- sqrt((positions[,4] - OriginX)^2 + (positions[,5]-OriginY)^2 + (positions[,6]-OriginZ)^2)
  return(distance.vector)
}
