options(stringsAsFactors = F)

# library(PharmacoGx)
library(CoreGx)
library(piano)
library(dplyr)
library(parallel) 
library(fgsea)
library(tibble)

# gse70138_info = list()
# gse70138_info$inst_info <- read.delim(  "../data/referencefile/GSE70138_Broad_LINCS_inst_info.txt", sep="\t", stringsAsFactors=F)
# gse70138_info$pert_info <- read.delim("../data/referencefile/GSE70138_Broad_LINCS_pert_info.txt", sep="\t", stringsAsFactors=F)
# gse70138_info$sig_info <- read.delim("../data/referencefile/GSE70138_Broad_LINCS_sig_info.txt", sep="\t", stringsAsFactors=F)
# gse70138_info$sig_metrics <- read.delim("../data/referencefile/GSE70138_Broad_LINCS_sig_metrics.txt", sep="\t", stringsAsFactors=F)
# 
# 
# 
# 
# 
# gse70138_info$cell_info <- read.delim("../data/referencefile/GSE92742_Broad_LINCS_cell_info.txt", sep="\t", stringsAsFactors=F)
# gse70138_info$gene_info <- read.delim("../data/referencefile/GSE92742_Broad_LINCS_gene_info.txt", sep="\t", stringsAsFactors=F)
# 

# gse92742_info = list()
# gse92742_info$inst_info <- read.delim("../data/referencefile/GSE92742_Broad_LINCS_inst_info.txt", sep="\t", stringsAsFactors=F)
# gse92742_info$pert_info <- read.delim("../data/referencefile/GSE92742_Broad_LINCS_pert_info.txt", sep="\t", stringsAsFactors=F)
# gse92742_info$pert_info_metrics <- read.delim("./referencefile/GSE92742_Broad_LINCS_pert_metrics.txt", sep="\t", stringsAsFactors=F)
# gse92742_info$sig_info <- read.delim("./referencefile/GSE92742_Broad_LINCS_sig_info.txt", sep="\t", stringsAsFactors=F)
# gse92742_info$sig_metrics <- read.delim("./referencefile/GSE92742_Broad_LINCS_sig_metrics.txt", sep="\t", stringsAsFactors=F)


# gse92742_info$cell_info <- read.delim("./referencefile/GSE92742_Broad_LINCS_cell_info.txt", sep="\t", stringsAsFactors=F)
# gse92742_info$gene_info <- read.delim("./referencefile/GSE92742_Broad_LINCS_gene_info.txt", sep="\t", stringsAsFactors=F)
# gse92742_info$gene_info_landmark <- read.delim("./referencefile/GSE92742_Broad_LINCS_gene_info_delta_landmark.txt", sep="\t", stringsAsFactors=F)


# geneID_drugs = gse92742_info$gene_info


loadGSE <- function(GSEname, filepath){
  gse_info = list()
  gse_info$inst_info <- read.delim(paste0(filepath,GSEname, "_Broad_LINCS_inst_info.txt") , sep="\t", stringsAsFactors=F)
  gse_info$pert_info <- read.delim(paste0(filepath,GSEname, "_Broad_LINCS_pert_info.txt"), sep="\t", stringsAsFactors=F)
  # gse_info$pert_info_metrics <- read.delim(paste0(filepath,GSEname, "_Broad_LINCS_pert_metrics.txt"), sep="\t", stringsAsFactors=F)
  gse_info$sig_info <- read.delim(paste0(filepath,GSEname, "_Broad_LINCS_sig_info.txt"), sep="\t", stringsAsFactors=F)
  gse_info$sig_metrics <- read.delim(paste0(filepath,GSEname, "_Broad_LINCS_sig_metrics.txt"), sep="\t", stringsAsFactors=F)
  
  
  gse_info$cell_info <- read.delim(paste0(filepath,"GSE92742_Broad_LINCS_cell_info.txt"), sep="\t", stringsAsFactors=F)
  gse_info$gene_info <- read.delim(paste0(filepath,"GSE92742_Broad_LINCS_gene_info.txt"), sep="\t", stringsAsFactors=F)
  # gse_info$gene_info_landmark <- read.delim(paste0(filepath,"GSE92742_Broad_LINCS_gene_info_delta_landmark.txt"), sep="\t", stringsAsFactors=F)
   return(gse_info)
  }


loadGCTXData = function(GSE, cell_line, filepath)
{
  library(cmapR)
  if(GSE == "GSE70138")
  {
    GSE_info = gse70138_info
    ds_path <- paste0(filepath,"GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx")
  }
  if(GSE == "GSE92742")
  {
    GSE_info = gse92742_info
    ds_path <- paste0(filepath,"GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx")
  }
  # gene_list = gene_list[gene_list$SYMBOL %in% GSE_info$gene_info$pr_gene_symbol,]
  # gene_info = GSE_info$gene_info[match(gene_list$SYMBOL, GSE_info$gene_info$pr_gene_symbol),]
  gene_info = GSE_info$gene_info
  
  idx <- which(GSE_info$sig_info$cell_id==cell_line)
  sig_ids <- GSE_info$sig_info$sig_id[idx]
  drug_ds <- parse_gctx(ds_path, cid=sig_ids, rid = as.character(GSE_info$gene_info$pr_gene_id))
  drug_matrix = list()
  # drug_matrix$target = gene_list
  drug_matrix$gene_info_drug = gene_info
  drug_matrix$drug_matrix = drug_ds
  return(drug_matrix)
}


# my_ds_rank_by_column <- rank_gct(gse92742_drug_data$drug_matrix, dim="col")
# m <- mat(gse92742_drug_data$drug_matrix)
# # plot z-score vs rank for the first 25 genes (rows)
# ranked_m <- mat(my_ds_rank_by_column)
# plot(ranked_m[1:25, ],
#      m[1:25, ],
#      xlab="rank",
#      ylab="differential expression score",
#      main="score vs. rank")
# rm(m,ranked_m,my_ds_rank_by_column)


getDrugMatrix = function(cell_line, filepath)
{
  # results$target = gene_list[gene_list$SYMBOL %in% gse70138_info$gene_info$pr_gene_symbol,]
  results = list()
  if(cell_line %in% unique(gse92742_info$sig_info$cell_id))
  {
    gse92742_drug_data = loadGCTXData("GSE92742", cell_line, filepath)
  }
  if(cell_line %in% unique(gse70138_info$sig_info$cell_id))
  {
    gse70138_drug_data = loadGCTXData("GSE70138", cell_line, filepath)
  }
  
  if(exists("gse92742_drug_data") & exists("gse70138_drug_data"))
  {
    results$drug_FC_matrix = cbind(gse92742_drug_data$drug_matrix@mat, gse70138_drug_data$drug_matrix@mat)
    results$drugs_info = rbind(gse92742_info$sig_info[,c(1:5,8,11,12)], gse70138_info$sig_info)
    results$gene_info = gse92742_drug_data$gene_info_drug
  }
  if(!exists("gse92742_drug_data") & exists("gse70138_drug_data"))
  {
    print("only gse70138")
    results$drug_FC_matrix = gse70138_drug_data$drug_matrix@mat
    results$drugs_info = gse70138_info$sig_info
    results$gene_info = gse70138_drug_data$gene_info_drug
  }
  if(exists("gse92742_drug_data") & !exists("gse70138_drug_data"))
  {
    print("only gse92742")
    results$drug_FC_matrix = gse92742_drug_data$drug_matrix@mat
    results$drugs_info = gse92742_info$sig_info[,c(1:5,8,11,12)]
    results$gene_info = gse92742_drug_data$gene_info_drug
  }
  results$drugs_info = results$drugs_info[results$drugs_info$cell_id == cell_line,]
  results$drugs_info$pert_idose = gsub(results$drugs_info$pert_idose, pattern = "ÂµM", replacement = "µM")
  # results$target_list_up = target_list_up
  # results$target_list_down = target_list_down
  return(results)
}



getConnectivityScore = function(drugs, target_list, method = 'fgsea')
{
  druginfo = drugs$drugs_info[-c(1)]
  rownames(druginfo) = drugs$drugs_info$sig_id
  if(method == 'fgsea'){
    res <- apply(drugs$drug_FC_matrix, 2,  function(x, y){ 
      return(connectivityScore(x=x, 
                               y=y, 
                               method= method, nperm=1000))},y=target_list)}
  else if(method == 'gwc'){
    res <- apply(drugs$drug_FC_matrix, 2,  function(x, y){ 
      return(connectivityScore(x=cbind(x,rep(1, times= length(x))), 
                               y=y, 
                               method= method, nperm=1000))},y=target_list)}
  
  else{ stop('PLEASE INPUT ONE OF THE METHODS: fgsea or gwc !')}
  
  
  rownames(res) <- c("Connectivity", "PValue")
  res <- t(res)
  res <- res[order(res[,1], decreasing=TRUE),]
  results <- merge(res, druginfo, by="row.names", all=TRUE) 
  results <- results[order(results$Connectivity), ]
  return(results)
}


          
getEnrichedScores = function(drugs, targets, ncore = 5)
{ 
  druginfo = drugs$drugs_info[-c(1)]
  rownames(druginfo) = drugs$drugs_info$sig_id
  
  gset <- cbind(gene = rownames(targets), set = ifelse(as.numeric(targets[, 1]) >= 0, "UP", "DOWN"))
  gset <- piano::loadGSC(gset)
  n=ncol(drugs$drug_FC_matrix)
  
  # res.up = res.down = data.frame(matrix(nrow = n,ncol=8)) 
  # for (i in 1:n) {
  #   fgseaMultilevelRes <- fgseaMultilevel(gset$gsc, stats=drugs$drug_FC_matrix[,i], maxSize=500)
  #   res.up[i,] = fgseaMultilevelRes %>% dplyr::filter(pathway == 'UP')
  #   res.down[i,] = fgseaMultilevelRes %>% dplyr::filter(pathway == 'DOWN')
  # }
  
  res <- lapply(1:n,  function(x) {
    fgseaMultilevelRes <- fgsea::fgseaMultilevel(gset$gsc, stats=drugs$drug_FC_matrix[,x], maxSize=500)
    res.up = fgseaMultilevelRes %>% dplyr::filter(pathway == 'UP')
    res.down = fgseaMultilevelRes %>% dplyr::filter(pathway == 'DOWN')
    return(list(res.up=res.up,res.down=res.down))})
  
  # cl <- makeCluster(ncore)
  # clusterEvalQ(cl, c(library(dplyr), library(fgsea)))
  # clusterExport(cl, c('gset','drugs','n'))
  # res <- parLapply(cl, 1:n,  function(x) {
  #   fgseaMultilevelRes <- fgsea::fgseaMultilevel(gset$gsc, stats=drugs$drug_FC_matrix[,x], maxSize=500)
  #   res.up = fgseaMultilevelRes %>% dplyr::filter(pathway == 'UP')
  #   res.down = fgseaMultilevelRes %>% dplyr::filter(pathway == 'DOWN')
  #   return(list(res.up=res.up,res.down=res.down))})
  # stopCluster(cl)

  res.down = do.call(rbind, lapply(res, function(x) x$res.down ) )
  res.up = do.call(rbind, lapply(res, function(x) x$res.up ) )
  
  res.down$row.names = res.up$row.names = colnames(drugs$drug_FC_matrix)
  druginfo = druginfo %>% tibble::rownames_to_column(var = 'row.names' )
  
  res.down = res.down %>% left_join(druginfo, by = 'row.names')
  res.up = res.up %>% left_join(druginfo, by = 'row.names')
  
  res = rbind(res.down,res.up)
  # colnames(res.down) =colnames(res.up) = colnames(fgseaMultilevelRes)
 
  # rownames(res.down) <- colnames(drugs$drug_FC_matrix)
  # res.down <- merge(res.down, druginfo, by="row.names", all=TRUE) 
  # res.up <- merge(res.up, druginfo, by="row.names", all=TRUE) 
  
  return(res)
}


  




