get.AlignedPositions <-
function(CIF.File.Location, Fasta.File.Location, chain.required = "A", RequiredModelNum =NULL,
                               patternQuality = PhredQuality(22L), subjectQuality = PhredQuality(22L),
                               type = "global-local", substitutionMatrix =NULL, fuzzyMatrix = NULL,
                               gapOpening = -10, gapExtension = -4, scoreOnly = FALSE){
  
  #The two lines below prevent a warning of occurring in the build. They have no effect on the code.
  subject <-NULL
  pattern <-NULL
	
  AtomSequence<- get.CIFSequence.Aligned(CIF.File.Location,chain.required, RequiredModelNum)
  AtomSequence<-cbind(AtomSequence, seq(from =1, to = dim(AtomSequence)[1],by=1))
  colnames(AtomSequence)=c(colnames(AtomSequence)[1:6],"Ali.Count")
 
  #gets the canonical and extracted protein and performs the alignment
  canonical.protein<-paste(get.AASeq(Fasta.File.Location),collapse="")
  protein.extracted<- paste(sapply(AtomSequence$Residue,get.SingleLetterCode,USE.NAMES=FALSE), collapse ="")
  alignment <-pairwiseAlignment(pattern = canonical.protein, subject=protein.extracted, patternQuality=patternQuality,
                               subjectQuality=subjectQuality,type = type, substitutionMatrix= substitutionMatrix,
                                fuzzyMatrix=fuzzyMatrix,gapOpening=gapOpening,gapExtension=gapExtension,
                                scoreOnly=scoreOnly)
								 
  #Gets the beginning position of alignment in the Atom Sequence (eg, cuts off the top of AtomSequence if not in alignment)
  StartSubjectPosition <- start(subject(alignment))
  AtomSequence<-AtomSequence[StartSubjectPosition:dim(AtomSequence)[1],]
  
  #gets the full string list for the pattern and subject  
  PatternString<-strsplit(toString(pattern(alignment)),split="")
  AlignedString<-strsplit(toString(subject(alignment)),split="")
  
 
    
  #makes a matrix to put them side by side and gets the start position
  AlignedMatrix<-data.frame(PatternString[[1]],AlignedString[[1]])
  colnames(AlignedMatrix)<-c("Can.Res","Ali.Res")
  StartPosition <- start(pattern(alignment))
  
  #attaches a columnn with counts
  AlignedMatrix<-cbind(AlignedMatrix, rep(0,dim(AlignedMatrix)[1]), rep(0,dim(AlignedMatrix)[1]))
  colnames(AlignedMatrix)<-c(colnames(AlignedMatrix)[1:2],"Can.Count","Ali.Count")
  Can.Count <-StartPosition-1
  Ali.Count <-StartSubjectPosition-1
  AlignedMatrix$Can.Count<-as.numeric(AlignedMatrix$Can.Count)
  AlignedMatrix$Ali.Count<-as.numeric(AlignedMatrix$Ali.Count)
  
  
  for(i in 1:dim(AlignedMatrix)[1]){
    if(as.character(AlignedMatrix$Can.Res[i]) != "-"){
      Can.Count <- Can.Count +1 
      AlignedMatrix$Can.Count[i]<-Can.Count
      AtomSequenceAligned <- NULL
    }
    else{
      AlignedMatrix$Can.Count[i]<- -1
    }
  
    if(as.character(AlignedMatrix$Ali.Res[i]) != "-"){
      Ali.Count <- Ali.Count + 1
      AlignedMatrix$Ali.Count[i]<-Ali.Count
    }
    else{
      AlignedMatrix$Ali.Count[i] <- -1
    }
  }
  
  
  #gets all the positions where sequences are aligned but with a substitution. 
  #gets the positions where the canonical protein is matched to gaps in the subject.
  #gets the positions where the subject protein is matched to gaps in the pattern
  MisMatches <- intersect(intersect(which(AlignedMatrix$Ali.Count!= -1),which(AlignedMatrix$Can.Count != -1)),which(as.character(AlignedMatrix$Can.Res)!=as.character(AlignedMatrix$Ali.Res)))
  Missing.Positions.in.Subject <-  which(AlignedMatrix$Ali.Count == -1)
  Missing.Positions.in.Pattern <- which(AlignedMatrix$Can.Count == -1)

    
  #gets the rows where the subject pattern is matched to (the canonical protein or a gap in the canonical protein)
  Subject.Positions.Matched.to.Pattern <- which(AlignedMatrix$Ali.Count != -1)
  Aligned.Trimmed.Matrix <- AlignedMatrix[Subject.Positions.Matched.to.Pattern,]


  #shortens the atom sequence if we didn't align all of them and some are excluded.
  if(dim(Aligned.Trimmed.Matrix)[1]<=dim(AtomSequence)[1]){ 
    AtomSequence <- AtomSequence[1: dim(Aligned.Trimmed.Matrix)[1],]
  }
  
  #At this point we should have Our Trimmed Matrix to be as long as the our AtomSequence.
  if(dim(Aligned.Trimmed.Matrix)[1] == dim(AtomSequence)[1]){
 
    AtomSequenceAligned <- AtomSequence
  
    #We save the full sequence into AtomSequence$Ali.Count to make it the Can.Count
    #At this point AtomSequence$Can.Count can have -1 in it, in the case the subject has elements that the pattern does not
    AtomSequenceAligned$Ali.Count<- Aligned.Trimmed.Matrix$Can.Count
    colnames(AtomSequenceAligned) <- c(colnames(AtomSequence)[1:6],"Can.Count")
   
  
    #Removes from the sequence positions where the canonical protein is not referenced
    AtomSequenceAligned <- AtomSequenceAligned[which(AtomSequenceAligned$Can.Count !=-1),]
    Aligned.Trimmed.Matrix<-Aligned.Trimmed.Matrix[which(Aligned.Trimmed.Matrix$Can.Count!=-1),]

    #At this point, we should again have the same length
    if(dim(Aligned.Trimmed.Matrix)[1] == dim(AtomSequenceAligned)[1]){
    
      #Gets which Canonical Counts in the Aligned Trimmed Matrix match up (thus we exclude those that don't match)
      sequence.matches <- Aligned.Trimmed.Matrix$Can.Count[which(as.character(Aligned.Trimmed.Matrix$Can.Res)==as.character(Aligned.Trimmed.Matrix$Ali.Res))]
      AtomSequenceAligned <- AtomSequenceAligned[match(sequence.matches, AtomSequenceAligned$Can.Count),]
    
    
      #some housekeeping  -- rearange the matrix into the same order as the main program expects
      AtomSequenceAligned<- AtomSequenceAligned[,c(2,7,3,4,5,6)]
  
      #Check Alignment
      extracted.letters <- sapply(AtomSequenceAligned$Residue,get.SingleLetterCode,USE.NAMES=FALSE)
      canonical.letters <- unlist(strsplit(canonical.protein, split=""))[AtomSequenceAligned$Can.Count] 
         
      diff.count = sum( extracted.letters != canonical.letters [AtomSequence$AtomCount])
      diff.positions = which(extracted.letters!=canonical.letters[AtomSequence$AtomCount])
 
      if(diff.count == 0){
       Result = "OK"
       diff.positions <- NULL
      }
      else if(diff.count != 0){
        Result = "Extraction Failure"
      }
    }
    else{
      print("Sequences of different length after GAP removal in Canonical Sequence")
      AlignedMatrix = NULL
      diff.count = NULL
      diff.positions = NULL
      Alignment.Result = alignment
      Result = "Error! Sequence Lengths Not Equal"  
    }
   }else{
    print("Sequences of different length after GAP removal in Subject Sequence")
    AlignedMatrix = NULL
    diff.count = NULL
    diff.positions = NULL
    Alignment.Result = alignment
    Result = "Error! Sequence Lengths Not Equal"
  
  }

  return.result <- list(Positions = AtomSequenceAligned, Diff.Count = diff.count, Diff.Positions = diff.positions, Alignment.Result = alignment, 
                        Result = Result)
  return(return.result)
}
