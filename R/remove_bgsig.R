# Wrap for subtracting the background from the target signal to estimate the KO-associated signature
Wrap_KOSig <- function(MutCatalogue,bg_column,ko_column,sampling_number, start_num,boundary,outputname){
  
  KOSig <- RemoveBackground_vector_single(MutCatalogue[,bg_column], MutCatalogue[,ko_column],sampling_number, start_num,boundary)
  KOSig$MutationType <- MutCatalogue[,"MutationType"]
  write.table(KOSig,paste0(outputname,".txt"),sep = "\t",col.names = T, row.names = F, quote = F)
  plotCountbasis(KOSig,1,6,9,paste0(outputname,".pdf"))
  plotPercentagebasis(KOSig,1,6,9,paste0(outputname,"_percentage.pdf"))
  
}

