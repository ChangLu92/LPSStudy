library(WGCNA)

correlation_with_traits <- function(MEs, traits, cex.lab,cex.text,wi,he,
                                    threthold_corr, threthold_pv, filename ) {
  
  # nGenes = ncol(datExpr);
  
  ind <- sapply(traits, is.factor)
  traits[ind] <- lapply(traits[ind], function(x) {a <- as.numeric(x); return(a)})
  
  
  nSamples = apply(traits,2,function(x) sum(!is.na(x)))
  names(nSamples) <- colnames(traits)

  moduleTraitCor = cor(MEs, traits, use = "pairwise.complete.obs")
  x=1:ncol(traits)
  moduleTraitPvalue = as.matrix(sapply(x,function(i){
    corPvalueStudent(moduleTraitCor[,i], nSamples[i])
    
  } ))
  
  
  cor_logic <- apply(moduleTraitCor,2,function(x) any(abs(x)>threthold_corr))
  pv_logic <- apply(moduleTraitPvalue,2,function(x) any(x<threthold_pv))
  

  if(ncol(traits)!= 1){
    traits <- traits[,which(cor_logic&pv_logic==TRUE)]
    moduleTraitCor <- moduleTraitCor[,which(cor_logic&pv_logic)]
    moduleTraitPvalue <- moduleTraitPvalue[,which(cor_logic&pv_logic)]
  }
  
  
  
  
  sizeGrWindow(10,6)
  pdf(file = paste0("../results/",filename, ".pdf"), wi = wi, he = he)
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(8, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(traits),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = cex.text,
                 cex.lab =cex.lab,
                 zlim = c(-1,1),
                 main = paste(filename))
  
  dev.off()
  
  my_list <- list("moduleTraitCor" = moduleTraitCor, "moduleTraitPvalue" = moduleTraitPvalue)
  return(my_list)
  
}
