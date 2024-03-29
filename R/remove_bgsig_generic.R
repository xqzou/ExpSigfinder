
#' Wrap for subtracting the background to estimate the experimentally generated signature
#' @param MutCatalogue 96-channel mutational catalog
#' @param bg_column the column name of the background catalog
#' @param ko_column the column names of the experiment catalog
#' @param sampling_number number of bootstrapping samples
#' @param start_num mutation burden to start
#' @param boundary Range of signature around centroid (default = 2)
#' @param channel_number number of channels in user-defined catalogue
#' @param outputname output file name
#' @export
Wrap_KOSig_channel <- function(MutCatalogue,bg_column,ko_column,sampling_number, start_num,boundary,channel_number,outputname){
  
  KOSig <- RemoveBackground_vector_single_channel(MutCatalogue[,bg_column], MutCatalogue[,ko_column],sampling_number, start_num,boundary,channel_number)
  KOSig$MutationType <- MutCatalogue[,"MutationType"]
  utils::write.table(KOSig,paste0(outputname,".txt"),sep = "\t",col.names = T, row.names = F, quote = F)
  
}


#' Return single vector
#' @param background_profile background mutational catalog
#' @param sig_profile experimentally generated mutational catalog
#' @param sampling_number number of bootstrapping samples
#' @param start_num mutation burden to start
#' @param boundary Range of signature around centroid (default = 2)
#' @param channel_number number of channels in user-defined catalogue
#' @return data.frame including background signature and experiment signature
#' @export 
RemoveBackground_vector_single_channel <- function(background_profile, sig_profile,sampling_number, start_num,boundary,channel_number){
  
  # Remove weak mutation types in sig_profile
  removeWeakMutationTypes <- 0.01
  genomesOriginal <- as.data.frame(sig_profile)
  Totalmutations <- sum(sig_profile)
  removeMutations_max <- removeWeakMutationTypes * Totalmutations
  removeIdx <- which(cumsum(rowSums(genomesOriginal[,1:dim(genomesOriginal)[2],drop = FALSE])[order(rowSums(genomesOriginal[,1:dim(genomesOriginal)[2],drop = FALSE]))])<= min(removeMutations_max,4))
  mutationTypesToRemoveSet <- order(rowSums(genomesOriginal[,1:dim(genomesOriginal)[2],drop = FALSE]))[removeIdx]
  genomesReducted <- as.data.frame(genomesOriginal[-mutationTypesToRemoveSet,])
  reducedMutationtypes <- dim(genomesReducted)[1]
  
  # bootstrap sig profile
  centroid_sig <- rowMeans(genomesReducted[,1:dim(genomesReducted)[2],drop = FALSE])
  RepSig <- matrix(rep(centroid_sig,sampling_number),ncol = sampling_number)
  Sig_bootstraps <- bootstrapGenomesfun2(RepSig,sum(centroid_sig))
  
  i=start_num
  reachLimit <- FALSE
  diff_all_save <- NULL
  while(i<sum(centroid_sig) & !reachLimit){
    # Remove the same weak mutation types in background profile
    backgroundReducted <- background_profile[-mutationTypesToRemoveSet]
    RepControl <- matrix(rep(backgroundReducted,sampling_number),ncol = sampling_number)
    bg_bootstraps <- bootstrapGenomesfun2(RepControl,i)
    
    # Range of background
    centroid_background <- rowMeans(bg_bootstraps)
    sd_background <- apply(bg_bootstraps,1,sd)
    boundary_background <- centroid_background+boundary*sd_background
    
    # Range of sig
    centroid_sig <- rowMeans(Sig_bootstraps)
    sd_sig <- apply(Sig_bootstraps,1,sd)
    boundary_sig <- centroid_sig+boundary*sd_sig
    
    
    diff_all_boundary <- boundary_sig-centroid_background
    diff_all <- centroid_sig-centroid_background
    
    if(length(which(diff_all_boundary<0))>0){
      reachLimit <- TRUE
    }
    
    if(length(which(diff_all_boundary<0))==0){
      diff_all_boundary_save <- diff_all_boundary
      diff_all_save <- diff_all
      diff_all_save[which(diff_all_save<0)] <- 0
    }
    
    
    
    i = i+1
    
  }
  
  if(length(diff_all_save)==0) {
    stop("You need to reduce the start_number! ", 
         "Exiting...",call.=FALSE)
  }
  # Add Weak mutations for experiment expoure
  exposure <- rep(0,channel_number)
  origArrayIndex <- 1
  for(i in 1:channel_number){
    if(! i %in% mutationTypesToRemoveSet){
      exposure[i] <- diff_all_save[origArrayIndex]
      origArrayIndex=origArrayIndex+1
    }
  }
  
  # Add Weak mutations for background
  background_exposure<- rep(0,channel_number)
  origArrayIndex <- 1
  for(i in 1:channel_number){
    if(! i %in% mutationTypesToRemoveSet){
      background_exposure[i] <- centroid_background[origArrayIndex]
      origArrayIndex=origArrayIndex+1
    }
  }
  
  
  return(data.frame("KO_exposure"=exposure,"background_exposure"=background_exposure))
  
  
}
