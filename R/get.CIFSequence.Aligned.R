get.CIFSequence.Aligned <-
function(sourcefile, chain_required="A", RequiredModelNum = NULL){
  filedata<- readLines(sourcefile, n=-1)
  AtomCount <- 0
  N <- 10^5
  AtomMatrix <- data.frame(SeqId = rep(NA,N), Residue = rep(NA,N), SideChain = rep(NA, N), XCoord = rep(NA,N), YCoord = rep(NA,N), ZCoord = rep(NA,N))
  atom_site_found = FALSE
  Atom_Site_Entries <-vector()
  ModelNumFound = FALSE

  for(i in 1:length(filedata)){
    if(substring(filedata[i],1,11)=="_atom_site."){
      atom_site_found = TRUE
      Atom_Site_Entries <- c(Atom_Site_Entries,scan(text=filedata[i], what="character",quiet = TRUE))
    }
    else if(substring(filedata[i],1,1)=="#"){ #We mark both as false to avoid an extra if condition
      atom_site_found = FALSE
    }
    
    else if(atom_site_found && substring(filedata[i],1,11)!="_atom_site." && substring(filedata[i],1,4)=="ATOM"){
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
  }
  

  AtomMatrix<-AtomMatrix[1:AtomCount,]
  AtomMatrix[,1]<-as.numeric(AtomMatrix[,1])
  AtomMatrix[,4]<-as.numeric(AtomMatrix[,4])
  AtomMatrix[,5]<-as.numeric(AtomMatrix[,5])
  AtomMatrix[,6]<-as.numeric(AtomMatrix[,6])

 
  return(AtomMatrix)
}
