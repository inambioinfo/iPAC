plotMDS <-
function(positions,Graph.Output.Path=NULL, Graph.File.Name=NULL, Graph.Title, ordering, Show.Graph){
  
  if(!is.null(Graph.Output.Path) || !is.null(Graph.File.Name)){
    filename <- paste(Graph.Output.Path,Graph.File.Name, sep ="")
  }
  
  if(Show.Graph == "N"){
    pdf(filename)
  }
  
  
  #Create the main plot with all the points
  answer<-scatterplot3d(x = positions[,4], y = positions[,5], z = positions[,6],type = 'o', color = 'blue',xlab = "x-axis", ylab = "y-axis", zlab = "z-axis",main = Graph.Title) 
  answer$points3d(positions[,4], positions[,5], positions[,6],col='red')

  x.axis.length <- max(positions[,4])-min(positions[,4]) #gets the x-axis length from smallest to largest
  x.axis.min <- min(positions[,4]) #gets the starting point of the x-axis
  mds.ordering.length<-max(ordering) - min(ordering) #gets the ordering length
  
  
  #The line calculation is as follows: p1<- (distance in ordering - min(ordering))/ordering.length
  #p1 is a percentage of how far along the ordering length (in the MDS scale) the amino acid is. Thus if p1 = 40% that means that the AA is 40%
  #of the way to the final ordering position.
  #we then calculate p1* length(x-axis) + min(x.axis). Specifically, we calculate how far along the AA is in physical space and offset it by the starting position.
  
  #this plots the line to each amino acid.
  for(i in 1:dim(positions)[1]){
    answer$points3d(c( ((ordering[i,] - min(ordering))/mds.ordering.length) * x.axis.length + x.axis.min, positions[i,4]), c(min(positions[,5]),positions[i,5]), c(min(positions[,6]),positions[i,6]), type = 'l', col = i, pch = 8, lty = 2)
  }
  
  if(Show.Graph == "N"){
    dev.off()
  }
  else if(Show.Graph == "Y" && (!is.null(Graph.Output.Path) && !is.null(Graph.File.Name))){
    dev.copy2pdf(file = filename)
  }
  


  
}
