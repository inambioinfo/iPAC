Plot.Protein.Linear <- function(path, colCount, cex=0.5, height = 1,width = 1, color.palette = "heat"){
  total<-length(path)
  if(color.palette == "heat"){
    pal <- rev(heat.colors(total))
  }else if(color.palette == "gray"){
    pal <- rev(gray(seq(1/2, 1, length=total)))
  }else if(color.palette == "topo"){
    pal <- rev(topo.colors(total))
  }else if(color.palette == "cm"){
    pal <- rev(cm.colors(total))
  }
  
  positions.accounted <- sort(path)  
 
  colCount <- colCount # number per row
  rowCount <- ceiling(total/colCount)
  
  par(xpd = TRUE)
  plot( c(1,colCount), c(0,rowCount), type="n", ylab="", xlab="",
        axes=FALSE, ylim=c(rowCount,0))
  title("Protein Reordering")
  
  for (j in 0:(rowCount-1))
  {
    base <- j*colCount
    remaining <- total - base
    RowSize <- ifelse(remaining < colCount, remaining, colCount)
    
    for(i in 1:RowSize){
      if(base + i <= total){
        current.pos <- positions.accounted[base+i]
        current.dist <- which(path == (current.pos))
        rect(i-width/2,j-height/2,i+width/2, j+height/2, border = "black", col = pal[current.dist])
        text(x = i, y =j, labels = paste(current.pos), cex =cex, col = SetTextContrastColor(pal[current.dist]))
      }
      
    }
  }
  par(xpd=FALSE)
}