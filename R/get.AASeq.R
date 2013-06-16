get.AASeq <-
function(filelocation){
  filedata <- readLines(filelocation, n=-1)
  AACount = 1
  AASeq <- vector()
  for(i in 1:length(filedata)){
    if(substring(filedata[i],1,1)!=">" & filedata[i]!=""){
      AASeq<-c(AASeq,substring(filedata[i],seq(from=1, to =nchar(filedata[i]),by=1),seq(from=1, to =nchar(filedata[i]),by=1)))
    }
  }
   if(AASeq[length(AASeq)]=="*"){
    return.value <- AASeq[1:length(AASeq)-1]
  }else{
    return.value <- AASeq[1:length(AASeq)]
  }
  return (return.value)
}
