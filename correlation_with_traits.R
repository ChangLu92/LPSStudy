library(WGCNA)

correlation_with_traits <- function(MEs, traits,  threthold_corr, threthold_pv ) {
  
  
  ind <- sapply(traits, is.factor)
  traits[ind] <- lapply(traits[ind], function(x) {a <- as.numeric(x); return(a)})
  
  nSamples = apply(traits,2,function(x) sum(!is.na(x)))
  names(nSamples) <- colnames(traits)

  moduleTraitCor = cor(MEs, traits, use = "pairwise.complete.obs")
  x=1:ncol(traits)
  moduleTraitPvalue = as.matrix(sapply(x,function(i){
    corPvalueStudent(moduleTraitCor[,i], nSamples[i])
  } 
  ))
  colnames(moduleTraitPvalue) = colnames(moduleTraitCor) 
  
  cor_logic <- apply(moduleTraitCor,2,function(x) any(abs(x)>threthold_corr))
  pv_logic <- apply(moduleTraitPvalue,2,function(x) any(x<threthold_pv))

  if(ncol(traits)!= 1){
    traits <- traits[,which(cor_logic&pv_logic==TRUE)]
    moduleTraitCor <- moduleTraitCor[,which(cor_logic&pv_logic)]
    moduleTraitPvalue <- moduleTraitPvalue[,which(cor_logic&pv_logic)]
  }
  my_list <- list("moduleTraitCor" = moduleTraitCor, "moduleTraitPvalue" = moduleTraitPvalue)
  return(my_list)
  
}
