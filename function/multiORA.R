#ORA
library(piano)
library(ggnewscale)
library(stringr)

GOFILE = 'G:/CTMM/CTMM by Chang/codes/AnalysisAfterClusteringProject/referencefile/c5.go.v7.4.symbols.gmt'
pathways = "G:/CTMM/CTMM by Chang/codes/AnalysisAfterClusteringProject/referencefile/c2.cp.v7.2.symbols.gmt"


multiORA = function(geneset,
                  NESname = NULL,
                  PVname = NULL,
                  genename = NULL,
                  cases, 
                  bg,
                  pvcutoff = 0.05){
  
  if(is.null(PVname)){
    selectedgeneset <- switch(  
      cases[1],  
      up= geneset[geneset[,NESname]>0, genename],
      dn= geneset[geneset[,NESname]<0, genename],
      all = geneset[, genename])
  }else{
    selectedgeneset <- switch(  
      cases[1],  
      up= geneset[(geneset[,NESname]>0 & geneset[,PVname]< pvcutoff), genename],
      dn= geneset[(geneset[,NESname]<0 & geneset[,PVname]< pvcutoff), genename],
      all = geneset[, genename])
  }
  
  pathwayfilelist <- switch(  
    cases[2],  
    pathway= pathways,
    GO= GOFILE
  ) 
  
  candidate_genes <- unique(selectedgeneset)
  sig_pathways <- runGSAhyper(genes = candidate_genes, universe = unique(bg) , gsc = loadGSC(pathwayfilelist),pcutoff = 1)
  sig_pathways_df <- as.data.frame(sig_pathways$resTab)  %>% tibble::rownames_to_column(var = "pathway")
  sig_pathways_df <- data.frame(t(apply(sig_pathways_df, 1, function(r){
    aux = unlist(strsplit( sub("_",";", r["pathway"]), ";" ))
    r["pathway"] = gsub("_", " ", aux[2])
    return(c(r, "source" = aux[1]))})))
  colnames(sig_pathways_df)[length(sig_pathways_df)] = 'source'
  
  #data for plotting
  PathwaysSelect <- sig_pathways_df %>%
    dplyr::rename(pvalue = `p.value`, AdjPvalu = `Adjusted.p.value`) %>% 
    dplyr::mutate(pathway = as.factor(pathway))
  # return(list(PathwaysSelect, sig_pathways$gsc))
  return(PathwaysSelect)
}


multiGSEA = function(geneset,
                    NESname = NULL,
                    genename,
                    cases){
  
  pathwayfilelist <- switch(  
    cases,  
    pathway= pathways,
    GO= GOFILE
  ) 
  selectedgeneset <- geneset[,NESname]
  names(selectedgeneset) <- geneset[,genename]
  
  sig_pathways <- runGSA(geneLevelStats = selectedgeneset,  
                         gsc = loadGSC(pathwayfilelist),
                         adjMethod = 'fdr',
                         geneSetStat="wilcoxon")

  sig_pathways_df <- as.data.frame(sig_pathways$resTab)  %>% tibble::rownames_to_column(var = "pathway")
  sig_pathways_df <- data.frame(t(apply(sig_pathways_df, 1, function(r){
    aux = unlist(strsplit( sub("_",";", r["pathway"]), ";" ))
    r["pathway"] = gsub("_", " ", aux[2])
    return(c(r, "source" = aux[1]))})))
  colnames(sig_pathways_df)[length(sig_pathways_df)] = 'source'
  
  #data for plotting
  PathwaysSelect <- sig_pathways_df %>%
    dplyr::rename(pvalue = `p.value`, AdjPvalu = `Adjusted.p.value`) %>% 
    dplyr::mutate(pathway = as.factor(pathway))
  return(PathwaysSelect)
  
}




reorganisedf <- function(PathwayLists){
  
  for (i in 1:length(PathwayLists)) {
    pl = PathwayLists[[i]]
    modulename = names(PathwayLists)[[i]]
    
    df <- do.call('rbind',lapply(1:length(pl),function(x) pl[[x]]))
    modulenames <- rep(modulename, times = nrow(df))
    class <- do.call('c', lapply(1:length(pl),
                                   function(x){
                                     rep(paste0(modulename,'_',names(pl)[x]), 
                                         times = nrow(pl[[x]]))}))
    if(i == 1){
      dfall = df
      classall = class
      modules = modulenames
    }else{
      dfall = rbind(dfall,df) 
      classall = c(classall,class)
      modules = c(modules,modulenames)
    }
    
  }
  
  sign = vector(mode="character",length = length(classall))
  sign[str_detect(classall,'up' )] = 'activated'
  sign[str_detect(classall,'dn' )] = 'suppressed'
  sign[str_detect(classall,'all' )] = 'all'
  # unique(sign)
  
  PathwayAll_df<- dfall %>%
    tibble::add_column(class = factor(classall), modules = modules ,.sign = sign)%>%  
    dplyr::rename(Count = `Significant..in.gene.set.`, OtherintheGeneset = `Non.significant..in.gene.set.`) 
  PathwayAll_df[,2:7] <- sapply(dfall[,2:7], as.numeric)
  
  
  
  # fedf <- do.call('rbind',lapply(1:length(femalePathwayList),function(x) femalePathwayList[[x]]))
  # genderfe <- rep(modulename1, times = nrow(fedf))
  # classfe <- do.call('c', lapply(1:length(femalePathwayList),
  #                                function(x){
  #                                  rep(paste0('female_',names(femalePathwayList)[x]), 
  #                                      times = nrow(femalePathwayList[[x]]))}))
  # 
  # 
  # madf <- do.call('rbind',lapply(1:length(malePathwayList),function(x) malePathwayList[[x]]))
  # genderma <- rep(modulename2, times = nrow(madf))
  # classma <- do.call('c', lapply(1:length(malePathwayList),
  #                                function(x){
  #                                  rep(paste0('male_',names(femalePathwayList)[x]), 
  #                                             times = nrow(femalePathwayList[[x]]))})) 
  # 
  # df = rbind(fedf,madf) 
  # class = c(classfe,classma)
  
  # sign = vector(mode="character",length = length(class))
  # sign[str_detect(class,'up' )] = 'activated'
  # sign[str_detect(class,'dn' )] = 'suppressed'
  # # unique(sign)
  # 
  # PathwayAll_df<- df %>%
  #   add_column(class = factor(class), gender = factor(c(genderfe,genderma)),.sign = sign)%>%  
  #   dplyr::rename(Count = `Significant..in.gene.set.`, OtherintheGeneset = `Non.significant..in.gene.set.`) 
  # PathwayAll_df[,2:7] <- sapply(df[,2:7], as.numeric)
  # # dplyr::mutate(GeneRatio = Count / (Count+OtherintheGeneset)) 
  
  return(PathwayAll_df)
}


showCompEnrichDotplot_UD <- function(PathwayAll_df,
                                  pathwayname = 'GO',
                                  split = FALSE,
                                  top = 10,
                                  pcutoff = 0.05,
                                  filename = NULL,
                                  wid = 1200,
                                  break_by = 2,
                                  height = 1000){
  
  if(pathwayname == 'GO'){
    PathwayAll_df_filter = PathwayAll_df %>% 
      dplyr::arrange(AdjPvalu) %>%
      dplyr::filter(AdjPvalu < pcutoff & source %in% c("GOBP" , "GOMF","GOCC") )
  }else if(pathwayname == 'pathway'){
    PathwayAll_df_filter = PathwayAll_df %>% 
      dplyr::arrange(AdjPvalu) %>%
      dplyr::filter(AdjPvalu < pcutoff & !(source %in% c("GOBP" , "GOMF","GOCC")))}
  
  
  clProf.reshape.df <- do.call("rbind", lapply(unique(as.character(PathwayAll_df_filter$class)), function(x, top){
    a <- PathwayAll_df_filter %>%
      dplyr::filter(class %in% x) %>% 
      dplyr::arrange(AdjPvalu) %>%
      mutate(logpv = -log10(AdjPvalu), class = sub('_([^_]*)$', '', class) )
      
    # a$pathway <- factor(a$pathway, levels =a$pathway[order(a$class)])
    a <- a[1:min(c(top,nrow(a))),]
  },top = top))
  
  
  clProf.reshape.df$pathway = sapply( strwrap(clProf.reshape.df$pathway, 100, simplify=FALSE), paste, collapse="\n" )
  
  # if(floor(min(clProf.reshape.df$logadjpv))>0){
  #   brs = seq(0 , floor(max(clProf.reshape.df$logadjpv)) ,by=break_by)
  #   labels = seq(0 , floor(max(clProf.reshape.df$logadjpv)) ,by=break_by)
  # }else if(floor(max(clProf.reshape.df$logadjpv))<0){
  #   brs = seq(floor(min(clProf.reshape.df$logadjpv)) ,0, by=break_by)
  #   labels = seq(floor(min(clProf.reshape.df$logadjpv)) ,0, by=break_by)
  # }else{
  # brs = seq( floor(min(clProf.reshape.df$logadjpv))  , floor(max(clProf.reshape.df$logadjpv)) ,by=break_by)
  # labels = abs(seq( floor(min(clProf.reshape.df$logadjpv))  , floor(max(clProf.reshape.df$logadjpv)) ,by=break_by))
  # }
  brs = seq( 0 , floor(max(abs(clProf.reshape.df$logpv)) ),by=break_by)
  labels = abs(brs)
  
  
  p<- ggplot()+
    geom_point(data = clProf.reshape.df , aes(x = class, y = pathway, size = Count,color =logpv)) +
    geom_point(data = clProf.reshape.df , aes(x = class, y = pathway , size = Count), shape = 1,colour = "black")+
    scale_color_gradient("-log(Adj.pv)",
                          # colours = c("blue","White","red"),
                          # values = c(0, -min(clProf.reshape.df$logadjpv)/(max(clProf.reshape.df$logadjpv)-min(clProf.reshape.df$logadjpv))   ,1),
                          low = "White",high = "red", space = "Lab" ,
                          breaks= brs ,
                          labels= labels,
                          # labels = abs(seq( floor(min(clProf.reshape.df$logadjpv))  , floor(max(clProf.reshape.df$logadjpv)) ,by=2)),
                          limits=c(0,max(abs(clProf.reshape.df$logpv))),
                          guide=guide_colorbar(nbin = 500, barheight = 9))+
    theme(legend.text=element_text(size=rel(0.9)), axis.text.y = element_text(size=rel(1)), axis.text.x = element_text(size=rel(1.2)) )+
    ylab('')+xlab('')
  
  
  if(split){p = p+ facet_wrap(~source, ncol = length(unique(clProf.reshape.df$source)))+ theme(strip.text = element_text(size = 14))}
  if(is.null(filename)){
    p
  }else{
    jpeg(filename, width = wid, height = height, res = 150)
    print(p)
    dev.off()
  }
  
}



showCompEnrichDotplot <- function(PathwayAll_df,
                                  pathwayname = 'GO',
                                  split = FALSE,
                                  top = 10,
                                  pcutoff = 0.05,
                                  filename = NULL,
                                  wid = 1200,
                                  tophighlight = 3,
                                  break_by = 2,
                                  height = 1000){
  
  if(pathwayname == 'GO'){
    PathwayAll_df_filter = PathwayAll_df %>% 
      dplyr::arrange(AdjPvalu) %>%
      dplyr::filter(AdjPvalu < pcutoff & source %in% c("GOBP" , "GOMF","GOCC") )
  }else if(pathwayname == 'pathway'){
    PathwayAll_df_filter = PathwayAll_df %>% 
      dplyr::arrange(AdjPvalu) %>%
      dplyr::filter(AdjPvalu < pcutoff & !(source %in% c("GOBP" , "GOMF","GOCC")))}
  
  # if(ifClusterProfile){
  #   label = c('female(suppressed)'='female\ndn',
  #                         'male(suppressed)'='male\ndn',
  #                         "female(activated)"= "female\nup",
  #                         "male(activated)" = "male\nup")
  # }
  
  # }else if(pathwayname == 'ClusterProfile'){
  #   PathwayAll_df_filter = PathwayAll_df %>% 
  #     dplyr::arrange(AdjPvalu) %>%
  #     dplyr::filter(AdjPvalu < pcutoff)
  #    label = c('female(suppressed)'='female\ndn',
  #             'male(suppressed)'='male\ndn',   
  #             "female(activated)"= "female\nup", 
  #             "male(activated)" = "male\nup")
    
  # }
  
  
   clProf.reshape.df <- do.call("rbind", lapply(unique(as.character(PathwayAll_df_filter$class)), function(x, top){
     a <- PathwayAll_df_filter %>%
       dplyr::filter(class %in% x) %>% 
       dplyr::arrange(AdjPvalu) %>%
       mutate(logpv = -log10(AdjPvalu),)
     # a$pathway <- factor(a$pathway, levels =a$pathway[order(a$class)])
     
     
     a <- a[1:min(c(top,nrow(a))),]
   },top = top))
   
   
   
   
   # updata <- clProf.reshape.df %>% dplyr::filter(.sign == 'activated')
   # updata$pathway <- factor(updata$pathway, levels = updata$pathway[order(updata$class)])
   
   # downdata <- clProf.reshape.df %>% dplyr::filter(.sign == 'suppressed') 
   # downdata$pathway <- factor(downdata$pathway, levels = downdata$pathway[order(downdata$AdjPvalu)])
   
   
   clProf.reshape.df$logadjpv[clProf.reshape.df$.sign=='suppressed'] = log10(clProf.reshape.df$AdjPvalu[clProf.reshape.df$.sign=='suppressed'])
   clProf.reshape.df$logadjpv[clProf.reshape.df$.sign=='activated'] = -log10(clProf.reshape.df$AdjPvalu[clProf.reshape.df$.sign=='activated'])
   
   clProf.reshape.df = clProf.reshape.df[clProf.reshape.df$.sign!='all',]
   
   clProf.reshape.df$pathway = sapply( strwrap(clProf.reshape.df$pathway, 50, simplify=FALSE), paste, collapse="\n" )
   
   # if(floor(min(clProf.reshape.df$logadjpv))>0){
   #   brs = seq(0 , floor(max(clProf.reshape.df$logadjpv)) ,by=break_by)
   #   labels = seq(0 , floor(max(clProf.reshape.df$logadjpv)) ,by=break_by)
   # }else if(floor(max(clProf.reshape.df$logadjpv))<0){
   #   brs = seq(floor(min(clProf.reshape.df$logadjpv)) ,0, by=break_by)
   #   labels = seq(floor(min(clProf.reshape.df$logadjpv)) ,0, by=break_by)
   # }else{
   # brs = seq( floor(min(clProf.reshape.df$logadjpv))  , floor(max(clProf.reshape.df$logadjpv)) ,by=break_by)
   # labels = abs(seq( floor(min(clProf.reshape.df$logadjpv))  , floor(max(clProf.reshape.df$logadjpv)) ,by=break_by))
   # }
   
    brs1 = seq(0, floor(max(abs(clProf.reshape.df$logadjpv)) ),by=break_by)
    brs = sort(unique(c(-brs1,brs1))) 
    labels = abs(brs)
     
     clProf.reshape.df$pathway <- factor(clProf.reshape.df$pathway, levels= unique(clProf.reshape.df$pathway))

     facelist = rep('plain', times = nrow(clProf.reshape.df))
     tophl = clProf.reshape.df %>% dplyr::top_n(tophighlight, wt= abs(logadjpv)) %>% dplyr::select(pathway)
     facelist[clProf.reshape.df$pathway %in% tophl$pathway] = 'bold' 
   
   p<- ggplot()+
     geom_point(data = clProf.reshape.df , aes(x = modules, y = pathway, size = Count,color =logadjpv)) +
     geom_point(data = clProf.reshape.df , aes(x = modules, y = pathway , size = Count), shape = 1,colour = "black")+
     scale_color_gradient2("-log(Adj.pv)",
                           # colours = c("blue","White","red"),
                           # values = c(0, -min(clProf.reshape.df$logadjpv)/(max(clProf.reshape.df$logadjpv)-min(clProf.reshape.df$logadjpv))   ,1),
                           midpoint = 0,
                           low = "blue", mid = "White",high = "red", space = "Lab" ,
                           
                           breaks= brs ,
                           labels= labels,
                           
                           # labels = abs(seq( floor(min(clProf.reshape.df$logadjpv))  , floor(max(clProf.reshape.df$logadjpv)) ,by=2)),
                           limits=c(-max(abs(clProf.reshape.df$logadjpv)),max(abs(clProf.reshape.df$logadjpv))),
                           guide=guide_colorbar(nbin = 500, barheight = 8))+
     theme(legend.text=element_text(size=rel(0.8)), 
           axis.text.y = element_text(size=rel(1.1), face = facelist), 
           axis.text.x = element_text(size=rel(1.2)) )+
     ylab('')+xlab('')
   
   

   
   
# if(nrow(downdata)!=0){
#   p <- ggplot()+
#     geom_point(data = downdata , aes(x = gender, y = as.character(AdjPvalu), size = Count,color = -log10(AdjPvalu))) + 
#     geom_point(data = downdata , aes(x = gender, y = as.character(AdjPvalu) , size = Count), shape = 1,colour = "black")+
#     # scale_x_discrete(labels = pathway)+
#     scale_color_gradient2("DOWN\n-log(Adj.pv)", 
#                           midpoint = -log10(downdata$AdjPvalu[1])/2, 
#                           low = "white", mid = "lightblue", high = "blue", space = "Lab" ,
#                           breaks=seq(1,-log10(downdata$AdjPvalu[1]),by=2),
#                           limits=c(1, -log10(downdata$AdjPvalu[1])),
#                           guide=guide_colorbar(nbin = 500, barheight = 7))+ new_scale_color()
# }else{p <- ggplot()}
#  
#    
#    # ggplot()+
#    #   +  geom_point(data = clProf.reshape.df , aes(x = gender, y = pathway, size = Count,color = -log10(AdjPvalu))) + 
#    #   + geom_point(data = clProf.reshape.df , aes(x = gender, y = pathway , size = Count), shape = 1,colour = "black")  
#    
#    
#      
# if(nrow(updata)!=0){
#   p <- p +
#     geom_point(data = updata , aes(x = gender, y = as.character(AdjPvalu), size = Count, color = -log10(AdjPvalu)))+ 
#     geom_point(data = updata , aes(x = gender, y = as.character(AdjPvalu), size = Count), shape = 1,colour = "black")+
#     scale_color_gradient2("UP\n-log(Adj.pv)", 
#                           midpoint = -log10(updata$AdjPvalu[1])/2, 
#                           low = "White", mid = "red",high = "darkred", space = "Lab" ,
#                           breaks=seq(1,-log10(updata$AdjPvalu[1]),by=1),
#                           limits=c(1, -log10(updata$AdjPvalu[1])),
#                           guide=guide_colorbar(nbin = 500, barheight = 7))+
#     theme(legend.text=element_text(size=rel(0.7)))+
#     ylab('')+xlab('')
# }

   
   # p <- ggplot()+
   #   geom_point(data = updata, aes(x = gender, y = pathway, size = Count,color = -log10(AdjPvalu))) +
   #   geom_point(data = updata, aes(x = gender, y = pathway, size = Count), shape = 1,colour = "black")+
   #   scale_color_gradient2("UP\n-log(Adj.pv)",
   #                         midpoint = -log10(updata$AdjPvalu[1])/2,
   #                         low = "White", mid = "red",high = "darkred", space = "Lab" ,
   #                         breaks=seq(1,-log10(updata$AdjPvalu[1]),by=1),
   #                         limits=c(1, -log10(updata$AdjPvalu[1])),
   #                         guide=guide_colorbar(nbin = 500, barheight = 7))+
   #   new_scale_color() +
   #   geom_point(data = downdata , aes(x = gender, y = pathway, size = Count, color = -log10(AdjPvalu)))+
   #   geom_point(data = downdata , aes(x = gender, y = pathway, size = Count), shape = 1,colour = "black")+
   #   scale_color_gradient2("DOWN\n-log(Adj.pv)",
   #                         midpoint = -log10(downdata$AdjPvalu[1])/2,
   #                         low = "white", mid = "lightblue", high = "blue", space = "Lab" ,
   #                         breaks=seq(1,-log10(downdata$AdjPvalu[1]),by=1),
   #                         limits=c(1, -log10(downdata$AdjPvalu[1])),
   #                         guide=guide_colorbar(nbin = 500, barheight = 7, reverse = TRUE))+
   #   # guides(color = guide_colourbar(reverse = TRUE))+
   #   theme(legend.text=element_text(size=rel(0.7)), axis.text.y = element_text(size=rel(1.2)) )+
   #   ylab('')+xlab('')
   
   
   
   
# p <- ggplot()+
#     geom_point(data = updata, aes(x = gender, y = pathway, size = Count,color = -log10(AdjPvalu))) +
#     geom_point(data = updata, aes(x = gender, y = pathway, size = Count), shape = 1,colour = "black")+
#     scale_color_gradient2("UP\n-log(Adj.pv)",
#                         midpoint = -log10(updata$AdjPvalu[1])/2,
#                         low = "White", mid = "red",high = "darkred", space = "Lab" ,
#                         breaks=seq(1,-log10(updata$AdjPvalu[1]),by=1),
#                         limits=c(1, -log10(updata$AdjPvalu[1])),
#                         guide=guide_colorbar(nbin = 500, barheight = 7))+
#     new_scale_color() +
#     geom_point(data = downdata , aes(x = gender, y = pathway, size = Count, color = -log10(AdjPvalu)))+
#     geom_point(data = downdata , aes(x = gender, y = pathway, size = Count), shape = 1,colour = "black")+
#     scale_color_gradient2("DOWN\n-log(Adj.pv)",
#                         midpoint = -log10(downdata$AdjPvalu[1])/2,
#                         low = "white", mid = "lightblue", high = "blue", space = "Lab" ,
#                         breaks=seq(1,-log10(downdata$AdjPvalu[1]),by=1),
#                         limits=c(1, -log10(downdata$AdjPvalu[1])),
#                         guide=guide_colorbar(nbin = 500, barheight = 7, reverse = TRUE))+
#     # guides(color = guide_colourbar(reverse = TRUE))+
#     theme(legend.text=element_text(size=rel(0.7)), axis.text.y = element_text(size=rel(1.2)) )+
#     ylab('')+xlab('')


  if(split){p = p+ facet_wrap(~source, ncol = length(unique(clProf.reshape.df$source)))+ theme(strip.text = element_text(size = 14))}
  if(is.null(filename)){
    p
  }else{
    jpeg(filename, width = wid, height = height, res = 150)
    print(p)
    dev.off()
  }
  
}



showALLCompEnrichDotplot <- function(PathwayAll_df,
                                  pathwayname = 'GO',
                                  split = FALSE,
                                  top = 10,
                                  pcutoff = 0.05,
                                  filename = NULL,
                                  wid = 1200,
                                  break_by = 2,
                                  height = 1000){
  
  if(pathwayname == 'GO'){
    PathwayAll_df_filter = PathwayAll_df %>% 
      dplyr::arrange(AdjPvalu) %>%
      dplyr::filter(AdjPvalu < pcutoff & source %in% c("GOBP" , "GOMF","GOCC") )
  }else if(pathwayname == 'pathway'){
    PathwayAll_df_filter = PathwayAll_df %>% 
      dplyr::arrange(AdjPvalu) %>%
      dplyr::filter(AdjPvalu < pcutoff & !(source %in% c("GOBP" , "GOMF","GOCC")))}

  
  
  clProf.reshape.df <- do.call("rbind", lapply(unique(as.character(PathwayAll_df_filter$class)), function(x, top){
    a <- PathwayAll_df_filter %>%
      dplyr::filter(class %in% x) %>% 
      dplyr::arrange(AdjPvalu) %>%
      mutate(logpv = -log10(AdjPvalu))
    # a$pathway <- factor(a$pathway, levels =a$pathway[order(a$class)])
    a <- a[1:min(c(top,nrow(a))),]
  },top = top))
  
  
  # updata <- clProf.reshape.df %>% dplyr::filter(.sign == 'activated')
  # downdata <- clProf.reshape.df %>% dplyr::filter(.sign == 'suppressed') 
  alldata <- clProf.reshape.df %>% dplyr::filter(.sign == 'all') 
  
  alldata$pathway = sapply( strwrap(alldata$pathway, 100, simplify=FALSE), paste, collapse="\n" )
  
  p<- ggplot()+
    geom_point(data = alldata , aes(x = modules, y = pathway, size = Count,color =logpv)) +
    geom_point(data = alldata , aes(x = modules, y = pathway , size = Count), shape = 1,colour = "black")+
    scale_color_gradient2("-log(Adj.pv)",
                          low= "White",high = "red", space = "Lab" ,
                          breaks= seq( 0, floor(max(alldata$logpv)) ,by=break_by),
                          labels= abs(seq( 0, floor(max(alldata$logpv)) ,by=break_by)),
                          guide=guide_colorbar(nbin = 500, barheight = 15))+
    theme(legend.text=element_text(size=rel(0.9)), axis.text.y = element_text(size=rel(1)), axis.text.x = element_text(size=rel(1.2)) )+
    ylab('')+xlab('')
  
  
  
  if(split){p = p+ facet_wrap(~source, ncol = length(unique(clProf.reshape.df$source)))+ theme(strip.text = element_text(size = 14))}
  if(is.null(filename)){
    p
  }else{
    jpeg(filename, width = wid, height = height, res = 150)
    print(p)
    dev.off()
  }
  
}


  


  
  



# 
# #------------------------
# TFPathwayAll_df_filter = TFPathwayAll_df %>% 
#   dplyr::arrange(AdjPvalu) %>%
#   dplyr::filter(AdjPvalu < 0.05 & !(source %in% c("GOBP" , "GOMF","GOCC")))
# clProf.reshape.df <- do.call("rbind", lapply(unique(as.character(TFPathwayAll_df_filter$class)), function(x, top){
#   a <- TFPathwayAll_df_filter %>%
#     dplyr::filter(class %in% x) %>% 
#     dplyr::arrange(AdjPvalu)
#   a <- a[1:min(c(top,nrow(a))),]
# },top = 10))
# 
# jpeg(paste0("./results/limma+batch+age/CVD_vs_Control/Pathway TF0.1.jpg"), width = 1500, height = 800, res = 150)
# ggplot(clProf.reshape.df, aes(x = class, y = pathway, size = Count, color=-log10(AdjPvalu)))+
#   geom_point() +
#   scale_color_gradient2(midpoint = 8, low = "blue", mid = "white",
#                         high = "red", space = "Lab" ,breaks=seq(0,max(-log10(clProf.reshape.df$AdjPvalu)),by=1),guide=guide_colorbar(nbin = 500, barheight = 15))+
#   facet_wrap(~source, ncol = length(unique(clProf.reshape.df$source)))+
#   theme(legend.text=element_text(size=rel(0.7)))+
#   scale_x_discrete(labels=c('female_dn_pathway'='female\ndn',
#                             'male_dn_pathway'='male\ndn',   
#                             "female_up_pathway"= "female\nup", 
#                             "male_up_pathway" = "male\nup"))+
#   ylab('')+xlab('')
# dev.off()




