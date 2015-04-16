ClusterFind <-
function(mutation.data, position.data,method = "MDS", alpha = 0.05, MultComp = "Bonferroni",Include.Culled = "Y", Include.Full = "Y",
                        create.map ="Y", Show.Graph = "Y", Graph.Output.Path=NULL,Graph.File.Name = "Map.pdf", Graph.Title = "Mapping",  
                        OriginX = min(position.data[,4]), OriginY=min(position.data[,5]), OriginZ=min(position.data[,6])){
  
  #mutation.data.culled removes the columns from the mutation matrix that do not correspond to positional entries
  mutation.data.culled<-mutation.data[,position.data[,2]]
  
  #If positions (x_i, x_i+n_i) are missing, skipped.position.left.vector will keep track of the x_i, and skipped.position.right.vector will keep track of the x_i+n_i
  LHS <- numeric()
  RHS <- numeric()
  
  missing.positions<-setdiff(c(1:dim(mutation.data)[2]), position.data[,2])
  
  LHS.found <- FALSE
  for(i in 1:dim(mutation.data)[2]){
    if(is.element(i, missing.positions) && LHS.found != TRUE){ #This block finds the left hand side
      LHS <- append(LHS, i)
      LHS.found <- TRUE
    }else if(!(is.element(i, missing.positions)) && LHS.found == TRUE){ #This block finds the right hand side
      RHS<- append(RHS, i-1)
      LHS.found  <- FALSE
    }else if(i == dim(mutation.data)[2] && is.element(i, missing.positions)){ #This block handles the case if the final aa is not included
      RHS <- append(RHS, i)
      LHS.found <- FALSE
    }
      
  }

  missing.values.matrix <- cbind(LHS, RHS)

  
  ################################################################################
  ###Begin Protein Positional Remapper
  ################################################################################
  
  if(method == "Linear"){  
    #Gets the Remapped Protein and Plots if required
    remapped<-LinearRemapperMain(position.data,OriginX, OriginY, OriginZ, mutation.data.culled)
    if(create.map == "Y"){
      plotLinear(OriginX=OriginX, OriginY=OriginY, OriginZ=OriginZ,positions=position.data, Graph.File.Name = Graph.File.Name, Graph.Title = Graph.Title,
                 Graph.Output.Path=Graph.Output.Path, Show.Graph = Show.Graph)
      }
  }else if(method == "MDS"){
    remapped<-MDSRemapperMain(position.data,mutation.data.culled) 
    if(create.map == "Y"){
      plotMDS(positions = position.data,Graph.Output.Path=Graph.Output.Path, Graph.File.Name = Graph.File.Name, Graph.Title = Graph.Title, ordering = remapped[[3]],Show.Graph = Show.Graph) 
    }
  }
  

  ################################################################################
  ###Begin The NMC Algorithm Here
  ################################################################################
  
  
  print("Running Remapped")
  if(inherits(try(ClusteringResults <- rbind(nmc(remapped[[1]],alpha=alpha,multtest=MultComp))),"try-error")==TRUE){
    ClusteringResults <- "ERROR! Error While running clustering on the remapped protein! Verify that your mutation data matrix is not all 0's in the remaining positions!"
  }
  
  #The original NMC algorithm returns a vector if there is only one element in the cluster. This causes the 
  #call to the dimensions of the matrix in 'MatrixInsert' to fail. As such, we call an rbind to FORCE it to 
  #be returned a matrix
  if(Include.Full == "Y"){
    print("Running Full")   
    if(inherits(try(OriginalPaperResult <- rbind(nmc(mutation.data, alpha=alpha,multtest = MultComp))),"try-error")==TRUE){
      OriginalPaperResult <- "ERROR! Error While running clustering on the culled linear protein! Verify that your mutation data matrix is not all 0's in the remaining positions!"
    }
  }else{
    OriginalPaperResult<- NULL
  }
  
  if(Include.Culled =="Y"){
    print("Running Culled")
    if(inherits(try(OriginalPaperResultCulled<- rbind(nmc(mutation.data.culled, alpha=alpha,multtest = MultComp))),"try-error")==TRUE){
      OriginalPaperResultCulled <- "ERROR! Error While running clustering on the full linear protein! Verify that your mutation data matrix is not all 0's!"
    }
    
  }else{
    OriginalPaperResultCulled<- NULL
  }
  
  
  #The code below adjusts for the missing AA if you remove some columns
  if(!is.null(OriginalPaperResultCulled) && substring(OriginalPaperResultCulled,1,6)!="ERROR!"){
    OriginalPaperAdjustedClusters <- MatrixInsert(missing.values.matrix, OriginalPaperResultCulled[,2:3])
    OriginalPaperResultCulled[,2:3]<- OriginalPaperAdjustedClusters
  }
  
  ################################################################################
  ###Map Clusters to Original Protein Sequence Then Calls MatrixInsert To Calculate The Correct Size
  ################################################################################
  Unmapped.Clustering.Results <- ClusteringResults
  if(!is.null(ClusteringResults) && substring(ClusteringResults,1,6)!="ERROR!"){
    
    #This part maps it to the original but does not account for "jumps". Thus if the position.data starts at 12, and the first clustering result is at 12-13, the output
    #would be 1 -2. Adjustment is needed to fill in the gaps
    for(i in 1: dim(ClusteringResults)[1]){
      Unmapped.Clustering.Results[i,2] <- remapped[[2]][ClusteringResults[i,2]]
      Unmapped.Clustering.Results[i,3] <- remapped[[2]][ClusteringResults[i,3]]
        if(Unmapped.Clustering.Results[i,2] > Unmapped.Clustering.Results[i,3] ){
          temp <- Unmapped.Clustering.Results[i,2]
          Unmapped.Clustering.Results[i,2]<- Unmapped.Clustering.Results[i,3]
          Unmapped.Clustering.Results[i,3]<- temp
        }
      Unmapped.Clustering.Results[i,4]<-sum(mutation.data.culled[,Unmapped.Clustering.Results[i,2]:Unmapped.Clustering.Results[i,3]])
    }
    
  
    
  AdjustedClusters <- MatrixInsert(missing.values.matrix, Unmapped.Clustering.Results[,2:3])
  Unmapped.Clustering.Results[,1]<- AdjustedClusters[,2]-AdjustedClusters[,1]+1
  Unmapped.Clustering.Results[,2:3] <- AdjustedClusters
  }
  
  ############################################################################
  ###Output Results
  ################################################################################
  final.result<-list(Unmapped.Clustering.Results,OriginalPaperResultCulled,OriginalPaperResult,missing.values.matrix)
  names(final.result)<- c("Remapped","OriginalCulled","Original","MissingPositions")
  return(final.result)
}
