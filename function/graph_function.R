
library(tibble)
library(WGCNA)
library(dplyr)
library(igraph)
library(RobustRankAggreg)

#pearson correlation and adjusted pvalue
corr.mat <- function(X,Y,adjmethod = 'fdr'){
  pv = matrix(nrow = ncol(X), ncol = ncol(Y))
  rho = matrix(nrow = ncol(X), ncol = ncol(Y))
  rownames(pv) = rownames(rho) = colnames(X)
  colnames(pv) = colnames(rho) = colnames(Y)
  
  for (i in 1:ncol(X)) {
    for(j in 1:ncol(Y)){
      corre = stats::cor.test(x = X[,i] ,y = Y[,j], na.action =  na.omit)
      pv[i,j] = corre$p.value
      rho[i,j] = corre$estimate
    }
  }
  adjpv = apply(pv,2 , function(x) p.adjust(x, method = adjmethod))
  return(list('rho' = rho,'pv' = pv,'adjpv'=adjpv))
}


generateSubnetwork <- function (datExpr, 
                   module, 
                   color, 
                   omitColors = "grey", 
                   power = 3, 
                   type = "unsigned",
                   t = 0.1){
  
  modules = names(table(module))
  if (!is.na(omitColors)[1]){modules = modules[!is.element(modules, omitColors)]} 
  if (is.null(colnames(datExpr))) {
    colnames(datExpr) = 1:dim(datExpr)[2]
    isIndex = TRUE
  }
  adj = WGCNA::adjacency(datExpr[, names(module)[module %in% color] ], power = power, type = type)
  
  # a <- sort(rowSums(adj),decreasing = TRUE)
  # b = as.data.frame(a) %>% rownames_to_column(var ='genes') %>% dplyr::rename(sumweight = 'a')

  adjm= adj
  adjm[adj<t] <- 0
  g <- graph_from_adjacency_matrix(adjm, weighted=TRUE, mode = "undirected")
  g <- igraph::simplify(g,remove.multiple = TRUE)
  return(g)
  
  
}

# calcalculate Hubgenes using aggregateRanks (RRA)
calcHubgenes <- function(vertexdf, genename, network){
  vertexdf$eigenvector <- evcent(network)$vector
  vertexdf$degree <- igraph::degree(network) 
  vertexdf$betweenness <- betweenness(network, directed = F) 
  vertexdf$closeness <- closeness(network)
  vertexdf$page_rank <- page_rank(network, directed = TRUE)$vector
  glist <-  lapply(vertexdf %>% select(eigenvector:page_rank),  function(x){
    names(x) <- vertexdf[,genename]
    return(names(sort(x, decreasing = TRUE)))
  })
  r = rankMatrix(glist)
  ranking_pv <- aggregateRanks(glist, exact = TRUE)
  return(ranking_pv)
}



# chooseTopHubgenes <- function (datExpr, module, color, vertexinfo, omitColors = "grey", power = 2, type = "unsigned",top = 0.1, t = 0.5) 
# {
#   modules = names(table(module))
#   if (!is.na(omitColors)[1]){modules = modules[!is.element(modules, omitColors)]} 
#   if (is.null(colnames(datExpr))) {
#     colnames(datExpr) = 1:dim(datExpr)[2]
#     isIndex = TRUE
#   }
#   adj = adjacency(datExpr[, names(module)[module %in% color] ], power = power, type = type)
#   a <- sort(rowSums(adj),decreasing = TRUE)
#   
#   hub = a[1:round(top*length(a))]
#   b = as.data.frame(a) %>% rownames_to_column(var ='genes') %>% dplyr::rename(sumweight = 'a')
#   
#   b$hubgenes[b$genes %in% names(hub)] = 'yes'
#   b$hubgenes[is.na(b$hubgenes)] = 'no'
#   
#   adjm= adj
#   adjm[adj<t] <- 0
#   g <- graph_from_adjacency_matrix(adjm, weighted=TRUE, mode = "undirected")
#   c = vertexinfo  %>% left_join(b,by='genes') 
#   c$direction[c$logFC > 0] = 1
#   c$direction[c$logFC < 0] = -1
#   vertex_attr(g) <- as.list(c)
#   g <- igraph::simplify(g,remove.multiple = TRUE)
#   # hub = names(a)[1:round(top*length(a))]
#   return(list(hub,g))
# }


DorotheatfGraph <- function(modulegenes, modulevertexdf, targetonly = TRUE){
  
  if(targetonly == TRUE){
    SubAracnenetworkR <- NULL
  }else{
    SubAracnenetworkR <- regulonsALL[regulonsALL$tf %in% modulegenes,]
  }
  
  SubAracnenetworkT <- regulonsALL[regulonsALL$target %in% modulegenes,]
  SubAracnenetworkX <- dplyr::bind_rows(SubAracnenetworkR,SubAracnenetworkT) %>% dplyr::distinct()
  SubAracnenetworkVector <- c(SubAracnenetworkX$tf,SubAracnenetworkX$target)
  SubAracnenetworkVectoruni <- unique(SubAracnenetworkVector)
  
  univetricesubDoro<- corrall2[corrall2$hgnc_symbol %in% SubAracnenetworkVectoruni,]
  univetricesubDoro  <-  univetricesubDoro %>% relocate(hgnc_symbol, .before = genes) 
  univetricesubDoro$TF[univetricesubDoro$hgnc_symbol %in% SubAracnenetworkX$tf] <- 'Yes'
  univetricesubDoro$TF[!(univetricesubDoro$hgnc_symbol %in% SubAracnenetworkX$tf)] <- 'No'
  
  univetricesubDoro = univetricesubDoro %>% left_join(modulevertexdf, by = 'genes')
  
  # univetricesubDoro$WGCNAweight = as.numeric(univetricesubDoro$WGCNAweight)
  # univetricesubDoro$WGCNAweight = (univetricesubDoro$WGCNAweight-min(univetricesubDoro$WGCNAweight))/(max(univetricesubDoro$WGCNAweight)-min(univetricesubDoro$WGCNAweight))
  
  SubAracnenetworkX = SubAracnenetworkX  %>% relocate(target, .before = confidence)
  SubAracnenetworkGraph <- graph_from_data_frame(SubAracnenetworkX,directed=TRUE,vertices = univetricesubDoro)
  E(SubAracnenetworkGraph)$weight = E(SubAracnenetworkGraph)$mor
  SubAracnenetworkGraphS <- igraph::simplify(SubAracnenetworkGraph,remove.multiple = TRUE,remove.loops = TRUE)
  return(list(SubAracnenetworkGraphS,univetricesubDoro))
  
}


moduletfGraph <- function(targetgenedf,
                          # modulevertexdf,
                          regulons,
                          vertattrdf
                          # tw = 10,
                          # td = 3
                          ){
  
  #SubAracnenetworkT <- ARACNEnet[ARACNEnet$Target %in% modulegenes,]
  
  a = sapply(regulons, function(x) any(names(x$tfmode) %in% targetgenedf$genes ))
  regulons_df_ara <- regulons[a]
  SubAracnenetworkT =  purrr::map_df(regulons_df_ara, .f = function(i) {
    tf_target = i$tfmode %>%
      enframe(name = "target", value="mor") %>%
      mutate(likelihood = i$likelihood)
  },
  .id = "tf")
  
  SubAracnenetworkX = SubAracnenetworkT[SubAracnenetworkT$target %in% targetgenedf$genes,]
  SubAracnenetworkVector <- c(SubAracnenetworkX$tf,SubAracnenetworkX$target)
  SubAracnenetworkVectoruni <- unique(SubAracnenetworkVector)
  univetriceDFSubAracnenetwork<- vertattrdf[vertattrdf$genes %in% SubAracnenetworkVectoruni, ]
  univetriceDFSubAracnenetwork$TF[univetriceDFSubAracnenetwork$genes %in% SubAracnenetworkX$tf] <- 'Yes'
  univetriceDFSubAracnenetwork$TF[!(univetriceDFSubAracnenetwork$genes %in% SubAracnenetworkX$tf)] <- 'No'
  univetriceDFSubAracnenetwork = univetriceDFSubAracnenetwork %>% left_join(targetgenedf, by = 'genes')
  
  
  
  # de = table(SubAracnenetworkT$tf)
  # wei = sapply(unique(SubAracnenetworkT$tf), function(x) sum(SubAracnenetworkT$MI[SubAracnenetworkT$Regulator == x]) ) 
  # names(wei) = unique(SubAracnenetworkT$Regulator)
  
  # de_wei = cbind(de,wei)
  
  # SubAracnenetworkX = SubAracnenetworkT[SubAracnenetworkT$Regulator %in% intersect(names(wei)[wei>tw], names(de)[de>td]),]
  # SubAracnenetworkVector <- c(SubAracnenetworkX$Regulator,SubAracnenetworkX$Target)
  # SubAracnenetworkVectoruni <- unique(SubAracnenetworkVector)
  
  # univetriceDFSubAracnenetwork<- vertattrdf[vertattrdf$genes %in% SubAracnenetworkVectoruni, ]
  # univetriceDFSubAracnenetwork$TF[univetriceDFSubAracnenetwork$genes %in% SubAracnenetworkX$Regulator] <- 'Yes'
  # univetriceDFSubAracnenetwork$TF[!(univetriceDFSubAracnenetwork$genes %in% SubAracnenetworkX$Regulator)] <- 'No'
  # univetriceDFSubAracnenetwork = univetriceDFSubAracnenetwork %>% left_join(modulevertexdf, by = 'genes')
  
  #univetriceDFSubAracnenetwork$comnodes[univetriceDFSubAracnenetwork$TF== 'Yes']  <- 'TF'
  #univetriceDFSubAracnenetwork$comnodes[univetriceDFSubAracnenetwork$TF== 'No' & univetriceDFSubAracnenetwork$hubgenes == 'yes']  <- 'Hub'
 # univetriceDFSubAracnenetwork$comnodes[univetriceDFSubAracnenetwork$TF== 'No' & univetriceDFSubAracnenetwork$hubgenes == 'no']  <- 'Others'
  
  SubAracnenetworkGraph <- graph_from_data_frame(SubAracnenetworkX,directed=TRUE,vertices = univetriceDFSubAracnenetwork)
  
  # E(SubAracnenetworkGraph)$weight = E(SubAracnenetworkGraph)$MI
  # E(SubAracnenetworkGraph)$mor = SubAracnenetworkGraph$mor
  SubAracnenetworkGraphS <- igraph::simplify(SubAracnenetworkGraph,remove.multiple = TRUE,remove.loops = TRUE,edge.attr.comb = "median")
  return(SubAracnenetworkGraphS)
}



# PLOT tf activity per module
moduletfactvityplot <- function(name, TFinfo, t = 5){
  
  # I=TFinfo$degree[rank(-TFinfo$degree) == 5]
  I = t
  p1 = ggplot(TFinfo, aes(x = TF_LogFC.NES, y =TF_BPD.cor)) + 
    geom_point(aes(size = degree)) + 
    theme_bw() + 
    stat_poly_eq(formula = y ~ x,
                 label.x = "right",
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~")),
                 parse = TRUE) +geom_smooth(method=lm)+
    geom_label_repel(data = TFinfo %>% filter(degree>=I ), aes(label = hgnc_symbol) ) #+ggtitle( paste(name, collapse="&") )
  # p2 = ggplot(TFinfo, aes(x = TF_LogFC.NES, y =TF_BPS.cor)) + 
  #   geom_point(aes(size = degree )) + 
  #   theme_bw() + 
  #   stat_poly_eq(formula = y ~ x,
  #                label.x = "right",
  #                eq.with.lhs = "italic(hat(y))~`=`~",
  #                aes(label = paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~")),
  #                parse = TRUE) +geom_smooth(method=lm)+
  #   geom_label_repel(data = TFinfo %>%  filter(degree>=I ), aes(label = hgnc_symbol) )
  
  jpeg(paste0('../results/',paste(name, collapse=" & "),'_tfDOROTHEA.jpg'), width = 800, height = 400,res = 100)
  # print(ggarrange(p1, p2, ncol = 2, nrow = 1, 
  #                 common.legend = TRUE, 
  #                 labels= paste(name, collapse=" & "),
  #                 legend="right",
  #                 vjust = 1,hjust = -0.2) )
  print(p)
  dev.off()
}