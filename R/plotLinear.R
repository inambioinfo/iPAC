plotLinear <-
function(OriginX, OriginY, OriginZ, positions, Graph.File.Name, Graph.Title, Graph.Output.Path, Show.Graph){
  if(!is.null(Graph.Output.Path) && !is.null(Graph.File.Name)){
    filename <- paste(Graph.Output.Path,Graph.File.Name, sep ="")
  }
  
  if(Show.Graph == "N"){
    pdf(filename)
  }
  
  #Create the main plot with all the points
  answer<-scatterplot3d(x = positions[,4], y = positions[,5], z = positions[,6],type = 'o', color = 'blue',xlab = "x-axis", ylab = "y-axis", zlab = "z-axis",main = Graph.Title) 
  answer$points3d(positions[,4], positions[,5], positions[,6],col='red')
  
  #this plots the line to each amino acid.
  for(i in 1:dim(positions)[1]){
    answer$points3d(c(OriginX,positions[i,4]), c(OriginY,positions[i,5]), c(OriginZ,positions[i,6]), type = 'l', col = "green", pch = 8, lty = 2)
  }

  if(Show.Graph == "N"){
    dev.off()
  }
  else if(Show.Graph == "Y" && !is.null(Graph.Output.Path)&& !is.null(Graph.File.Name)){
    dev.copy2pdf(file = filename)
  }
 

  
}
