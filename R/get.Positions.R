get.Positions <-
function(CIF.File.Location, Fasta.File.Location, chain.required = "A",RequiredModelNum =NULL){
  PDBXPoly<- get.CIFSequence(CIF.File.Location,chain.required,RequiredModelNum =NULL)
  PDBXPoly$Sequence<- PDBXPoly$Sequence[which(PDBXPoly$Sequence$AuthSeqNum!="?"),]
  PDBXPoly$Diffs.Matrix<- PDBXPoly$Diffs.Matrix[which(PDBXPoly$Diffs.Matrix$Remark!= "EXPRESSION TAG"),] #Removes Expression Tags since no longer needed
  StartAlign <-PDBXPoly$Align.Canonical[1]
  mapping <- match(PDBXPoly$Sequence$Seq_id,PDBXPoly$AtomMatrix$SeqId)
  
  PDBXPoly$AtomMatrix<-PDBXPoly$AtomMatrix[mapping,]
  
  aligned.positions<-PDBXPoly$Sequence$Count +StartAlign-1 #-1 to account for the nomenclature that data frames start at index 1
  
  if(length(aligned.positions)==dim(PDBXPoly$AtomMatrix)[1]){
    PDBXPoly$AtomMatrix<-cbind(PDBXPoly$AtomMatrix, aligned.positions)
    if(length(PDBXPoly$Diffs.Matrix$Seq_num)!=0){
      Diffs.Mapping<-match(PDBXPoly$Diffs.Matrix$Seq_num, PDBXPoly$AtomMatrix$SeqId)
      
      #THe line below removes all the difference which map to NA, eg they are outside the sequence that we are pulling
      PDBXPoly$Diffs.Matrix$Seq_num[which(Diffs.Mapping!= "NA")] <- Diffs.Mapping[which(Diffs.Mapping!="NA")]
      PDBXPoly$Diffs.Matrix<- PDBXPoly$Diffs.Matrix[which(Diffs.Mapping!="NA",arr.ind=TRUE),]
      PDBXPoly$Diffs.Matrix<-cbind(PDBXPoly$Diffs.Matrix, PDBXPoly$AtomMatrix[PDBXPoly$Diffs.Matrix$Seq_num,]$aligned.positions)
      PDBXPoly$Diffs.Matrix<-PDBXPoly$Diffs.Matrix[,c(1,2,5,4)]
      names(PDBXPoly$Diffs.Matrix)<-c("PDB.Residue","Canonical.Residue","Canonical.Num","Remark")
    }
  }else{
    print("Misalignment Occurred!")
    print("Manual Alignment Required!")
    PDBXPoly$AtomMAtrix <-NULL
  }
  
  
  
  #VerifyStep
  AAseq.ref<- get.AASeq(Fasta.File.Location)[PDBXPoly$AtomMatrix$aligned.positions]
  AAseq.PDB<- sapply(PDBXPoly$AtomMatrix$Residue,get.SingleLetterCode,USE.NAMES=FALSE)
  try(diffs <- sum(!AAseq.ref == AAseq.PDB))
  mismatched <- which(AAseq.ref!=AAseq.PDB)
  
  if(is.na(diffs)){
    Result = "Failure"
    Mismatch <- NULL
  }
  else if(diffs !=0){ #Diffs between Canonical and PDB Exists. Checks if explained with PDB file.
    Mismatch<-data.frame(AAseq.PDB[which(AAseq.ref!=AAseq.PDB)],
                         AAseq.ref[which(AAseq.ref!=AAseq.PDB)],
                         PDBXPoly$AtomMatrix[which(AAseq.ref!=AAseq.PDB),7])
    colnames(Mismatch)<-c("PDB.Residue","Canonical.Residue","Canonical.Num")
    
    if(diffs!= length(PDBXPoly$Diffs.Matrix$Canonical.Num)){#Number of inconsistencies does not match up with external count. Automatic Failure
      Result = "Failure"
      
    }
    else if(Mismatch$Canonical.Num==PDBXPoly$Diffs.Matrix$Canonical.Num && Mismatch$PDB.Residue ==PDBXPoly$Diffs.Matrix$PDB.Residue &&
      Mismatch$Canonical.Residue ==PDBXPoly$Diffs.Matrix$Canonical.Residue){ #Diffs Explained within PDB file. OK to use.
      Result = "OK"
    }else{
      Result = "Failure" #Diffs Not Explained.
    }
    
  }
  else{ #No Diffs between Canonical and PDB
    Mismatch <- NULL
    PDBXPoly$Diffs.Matrix<-NULL
    Result = "OK"
  }
  
  Final.Position.Matrix<- data.frame(PDBXPoly$AtomMatrix$Residue, as.numeric(PDBXPoly$AtomMatrix$aligned.positions), 
                                     PDBXPoly$AtomMatrix$SideChain, as.numeric(PDBXPoly$AtomMatrix$XCoord), 
                                     as.numeric(PDBXPoly$AtomMatrix$YCoord), as.numeric(PDBXPoly$AtomMatrix$ZCoord))
  colnames(Final.Position.Matrix)<-c("Residue","Can.Count","SideChain","XCoord","YCoord","ZCoord")
  
  return.value <- list(Positions = Final.Position.Matrix, External.Mismatch=Mismatch, PDB.Mismatch=PDBXPoly$Diffs.Matrix,Result=Result)
  return(return.value)
}
