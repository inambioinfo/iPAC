get.CIFSequence <-
function(sourcefile, chain_required="A",RequiredModelNum = NULL){
  filedata<- readLines(sourcefile, n=-1)
  AACount<- 0
  AtomCount <- 0
  DiffsCount <-0
  N <- 10^5
  PositionMatrix <- data.frame(Count = rep(NA, N), Seq_id = rep(NA,N), Residue = rep(NA,N), PDBNum=rep(NA,N), AuthSeqNum = rep(NA,N))
  AtomMatrix <- data.frame(SeqId = rep(NA,N), Residue = rep(NA,N), SideChain = rep(NA, N), XCoord = rep(NA,N), YCoord = rep(NA,N), ZCoord = rep(NA,N))
  pdbx_poly_found = FALSE
  struct_ref_found = FALSE
  atom_site_found = FALSE
  diffs_found = FALSE
  Diffs_Matrix <- data.frame(PDB.Residue = rep(NA,N),Canonical.Residue = rep(NA,N),Seq_num = rep(NA,N), Remark = rep(NA,N))
  Diffs_Entries <-vector()
  PDBX_Entries <- vector()
  Struct_Entries <- vector()
  Atom_Site_Entries <-vector()
  ModelNumFound = FALSE
  
  ExpressionTagPositions <- vector()
  
  Align.Canonical.Begin =0
  Align.Canonical.End = 0
  Alignment.Canonical = 0
  Align.PDB.Begin = 0
  Align.PDB.End = 0
  Alignment.PDB =0

  for(i in 1:length(filedata)){
    if(substring(filedata[i],1,16) =="_struct_ref_seq."){
      struct_ref_found= TRUE
      items<- scan(text=filedata[i], what="character",quiet = TRUE)
      if(length(items)==1){
        Struct_Entries<- c(Struct_Entries, gdata::trim(filedata[i]))
      }else if(length(items)==2){
        if(items[1] == "_struct_ref_seq.db_align_beg"){
          Align.Canonical.Begin = as.numeric(items[2])
          Struct_Entries<-c(Struct_Entries, items[1])
        }
        else if(items[1]=="_struct_ref_seq.db_align_end"){
          Align.Canonical.End = as.numeric(items[2])
          Struct_Entries<-c(Struct_Entries, items[1])
        }
        else if(items[1]=="_struct_ref_seq.pdbx_auth_seq_align_beg"){
          Align.PDB.Begin = as.numeric(items[2])
          Struct_Entries<-c(Struct_Entries, items[1])
        }
        else if(items[1]=="_struct_ref_seq.pdbx_auth_seq_align_end"){
          Align.PDB.End = as.numeric(items[2])
          Struct_Entries<-c(Struct_Entries, items[1])
        }
        Alignment.Canonical = c(as.numeric(Align.Canonical.Begin),as.numeric(Align.Canonical.End))
        Alignment.PDB =c(as.numeric(Align.PDB.Begin), as.numeric(Align.PDB.End))
      }
      
    }
    else if(substring(filedata[i],1,22)=="_pdbx_poly_seq_scheme."){
      pdbx_poly_found = TRUE
      PDBX_Entries<- c(PDBX_Entries, gdata::trim(filedata[i]))
    }
    else if(substring(filedata[i],1,11)=="_atom_site."){
      atom_site_found = TRUE
      Atom_Site_Entries <- c(Atom_Site_Entries,gdata::trim(filedata[i]))
    }
    else if(substring(filedata[i],1,20)=="_struct_ref_seq_dif."){
      diffs_found = TRUE
      Diffs_Entries<- c(Diffs_Entries,gdata::trim(filedata[i]))
      
    }
    else if(substring(filedata[i],1,1)=="#"){ #We mark both as false to avoid an extra if condition
      pdbx_poly_found = FALSE
      struct_ref_found = FALSE
      atom_site_found = FALSE
      diffs_found = FALSE
    }
    
    else if(diffs_found == TRUE && substring(filedata[i],1,20)!="_struct_ref_seq_dif."){
      items<- scan(text=filedata[i], quote="'", what="character",quiet = TRUE)
      if(items[which(Diffs_Entries == "_struct_ref_seq_dif.pdbx_pdb_strand_id")] == chain_required){
        DiffsCount <- DiffsCount +1
        Diffs_Matrix[DiffsCount,]<-c(get.SingleLetterCode(items[which(Diffs_Entries == "_struct_ref_seq_dif.mon_id")]),
                                     get.SingleLetterCode(items[which(Diffs_Entries == "_struct_ref_seq_dif.db_mon_id")]),
                                     items[which(Diffs_Entries == "_struct_ref_seq_dif.seq_num")],
                                     items[which(Diffs_Entries =="_struct_ref_seq_dif.details")])
          if(items[which(Diffs_Entries =="_struct_ref_seq_dif.details")] == "EXPRESSION TAG"){
            ExpressionTagPositions <- c(ExpressionTagPositions,as.numeric(items[which(Diffs_Entries == "_struct_ref_seq_dif.seq_num")]) )
          }
      }
      
    }
    
    else if(struct_ref_found == TRUE && substring(filedata[i],1,16) !="_struct_ref_seq."){
      items<- unlist(strsplit(filedata[i]," +"))
      if(items[which(Struct_Entries=="_struct_ref_seq.pdbx_strand_id")] == chain_required){
        Align.Canonical.Begin = items[which(Struct_Entries=="_struct_ref_seq.db_align_beg")]
        Align.Canonical.End = items[which(Struct_Entries=="_struct_ref_seq.db_align_end")]
        Alignment.Canonical = c(as.numeric(Align.Canonical.Begin),as.numeric(Align.Canonical.End))
        Align.PDB.Begin = items[which(Struct_Entries=="_struct_ref_seq.pdbx_auth_seq_align_beg")]
        Align.PDB.End = items[which(Struct_Entries=="_struct_ref_seq.pdbx_auth_seq_align_end")]
        Alignment.PDB =c(as.numeric(Align.PDB.Begin), as.numeric(Align.PDB.End))
      }
    }
    else if(atom_site_found && substring(filedata[i],1,11)!="_atom_site."){
      items<-unlist(strsplit(filedata[i]," +"))
      
      if(ModelNumFound == FALSE){
        first_model_found = as.numeric(items[which(Atom_Site_Entries == "_atom_site.pdbx_PDB_model_num")])
        if(is.null(RequiredModelNum)){
          RequiredModelNum = first_model_found
        }
        ModelNumFound = TRUE
      }
      
      if(items[which(Atom_Site_Entries == "_atom_site.auth_asym_id")] == chain_required &&
         items[which(Atom_Site_Entries == "_atom_site.auth_atom_id")] == "CA" &&
         as.numeric(items[which(Atom_Site_Entries == "_atom_site.pdbx_PDB_model_num")]) == RequiredModelNum){
          altID = items[which(Atom_Site_Entries == "_atom_site.label_alt_id")]
          if(altID == "." || altID =="A"){
            AtomCount <- AtomCount+1
            AtomMatrix[AtomCount,] <- c(items[which(Atom_Site_Entries == "_atom_site.label_seq_id")],
                                        items[which(Atom_Site_Entries == "_atom_site.label_comp_id")],
                                        items[which(Atom_Site_Entries == "_atom_site.auth_asym_id")],
                                        items[which(Atom_Site_Entries == "_atom_site.Cartn_x")],
                                        items[which(Atom_Site_Entries == "_atom_site.Cartn_y")],
                                        items[which(Atom_Site_Entries == "_atom_site.Cartn_z")])
          }
        
      }
      
    }
    else if(pdbx_poly_found == TRUE && substring(filedata[i],1,22)!="_pdbx_poly_seq_scheme."){
      items<- unlist(strsplit(filedata[i]," +"))
      strandIDCol = which(PDBX_Entries=="_pdbx_poly_seq_scheme.pdb_strand_id")
      PDBSequenceCol = which(PDBX_Entries=="_pdbx_poly_seq_scheme.pdb_seq_num")
      ResidueNameCol = which(PDBX_Entries == "_pdbx_poly_seq_scheme.mon_id")
      SeqIdCol = which(PDBX_Entries=="_pdbx_poly_seq_scheme.seq_id") 
      MissingCol = which(PDBX_Entries =="_pdbx_poly_seq_scheme.auth_seq_num")
      
      if(items[strandIDCol]==chain_required){
        
        #Second Condition Cuts out expression tags since they should not influence AACount
        if(as.numeric(items[PDBSequenceCol]) >= Align.PDB.Begin && !(as.numeric(items[SeqIdCol]) %in% ExpressionTagPositions)){
            if(as.numeric(items[PDBSequenceCol])>0){
              AACount = AACount+1
              #cat(sprintf("Row %s: Data: %s: %s, %s,\n", AACount,filedata[i],items[ResidueNameCol],items[PDBSequenceCol] ))
              PositionMatrix[AACount,]<- c(AACount, Seq_id =items[SeqIdCol],items[ResidueNameCol],items[PDBSequenceCol],items[MissingCol])
            }
        }
      }
      
    }

  }
  
  if(DiffsCount>0){
    Diffs_Matrix<-Diffs_Matrix[1:DiffsCount,]
    Diffs_Matrix[,3]<-as.numeric(Diffs_Matrix[,3])
  }else{
    Diffs_Matrix<-NULL
  }
  
  AtomMatrix<-AtomMatrix[1:AtomCount,]
  AtomMatrix[,1]<-as.numeric(AtomMatrix[,1])
  AtomMatrix[,4]<-as.numeric(AtomMatrix[,4])
  AtomMatrix[,5]<-as.numeric(AtomMatrix[,5])
  AtomMatrix[,6]<-as.numeric(AtomMatrix[,6])

  PositionMatrix<-PositionMatrix[1:AACount,]
  PositionMatrix[,1]<- as.numeric(PositionMatrix[,1])
  PositionMatrix[,2]<- as.numeric(PositionMatrix[,2])
  PositionMatrix[,4]<- as.numeric(PositionMatrix[,4])
  result<- list(PositionMatrix,AtomMatrix,Alignment.Canonical, Alignment.PDB, Diffs_Matrix)
  names(result)<-c("Sequence","AtomMatrix","Align.Canonical", "Align.PDB", "Diffs.Matrix")
  return(result)
}
