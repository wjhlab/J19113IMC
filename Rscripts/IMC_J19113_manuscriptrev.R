.libPaths(paste0(work,"/Libraries"))


rm(list = ls())

####READ and CLUSTER FUNCTIONS####
returnfcs <- function(FDR_cutoff=.05,
                      metaDataFile=NULL,
                      panelDataFile=NULL,
                      dataDirectory=NULL,
                      shape_conditions=NULL,
                      color_conditions=NULL){
  ## This function generates an fcs file, subtype_markers, colors and shapes for clustering 
  require(scales);require(readxl);require(dplyr);require(flowCore)
  ## Directory and metadatafile checking
  if(!dir.exists(dataDirectory)) {stop('ERR: cannot find data directory')}
  if(!file.exists(metaDataFile)) {stop('ERR: cannot find metadata.xlsx or .csv file')}
  ## Read-in metadata and clean
  ifelse(grepl(metaDataFile,pattern='.xls'),md <- read_excel(metaDataFile),md <- read.csv(metaDataFile,header = TRUE))#must be in xl format or csv
  md$file_name <- factor(md$file_name)
  md$File_order <- factor(md$File_order)
  md$Case <- factor(md$Case)
  md$Tumor <- factor(md$Tumor)
  md$Timepoint <- factor(md$Timepoint)
  md$TimeGroup <- factor(md$TimeGroup)
  
  rownames(md) = md$sample_id
  md$sample_id <- md$sample_id
  
  ## Make sure all files in metadata present in datadirectory
  if(!all(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])){
    print(paste('ERR: not all filenames in metadata present in data folder - missing',md$file_name[!which(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])],'Subsetting...'))
    md <- md[-c(!which(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])),]
  }
  
  ## Read fcs into fcs_raw
  fcs_raw <- read.flowSet(paste0(dataDirectory,"/",md$file_name), transformation = FALSE, truncate_max_range = FALSE)
  panel <- read_xlsx(panelDataFile)
  head(data.frame(panel))
  panel$Parameter <- gsub('-', '_', panel$Parameter)
  
  
  ## Export out the parameter/panel data from the flowFrame to edit
  ## use panel$Antigen to fix description in panel_fcs
  ## use metal+isotope as mapping between panel from xlsx and panel from the fcs files
  
  panel_fcs <- pData(parameters(fcs_raw[[1]]))
  panel_fcs$desc <- gsub('-', '_', panel_fcs$desc)
  panel_fcs$desc[is.na(panel_fcs$desc)] <- paste0('NA_',which(is.na(panel_fcs$desc)))  
  
  rownames(panel_fcs) = panel_fcs$name
  
  ## Replace desc with revised Name
  panel_fcs[panel$Parameter,]$desc<-panel$Name
  
  ## Replace parameter data in flowSet with edits
  pData(parameters(fcs_raw[[1]])) <- panel_fcs
  
  ## Assign objects to marker lists
  subtype_markers <- panel$Name[panel$Subtype == 1]
  functional_markers <- panel$Name[panel$Functional == 1]
  otherparameters <- panel$Name[panel$Other ==1]
  cluster_by <- panel$Name[panel$Cluster == 1]
  
  ## Check marker lists
  if(!all(subtype_markers %in% panel_fcs$desc)){stop('ERR: Not all subtype_markers in panel_fcs$desc')}
  if(!all(functional_markers %in% panel_fcs$desc)){stop('ERR: Not all functional_markers in panel_fcs$desc')}
  if(!all(otherparameters %in% panel_fcs$desc)){stop('ERR: Not all otherparameters in panel_fcs$desc')}
  if(!all(cluster_by %in% panel_fcs$desc)){stop('ERR: Not all cluster markers in panel_fcs$desc')}
  
  fcs <- fsApply(fcs_raw, function(x, cofactor = 0.8){
    colnames(x) <- panel_fcs$desc
    expr <- exprs(x)
    expr <- asinh(expr[, union(subtype_markers,functional_markers)] / cofactor)
    exprs(x) <- expr
    x
  })
  
  ## Save out the original rownames from the parameter list from the fcs flowFrame
  panellist <- rownames(pData(parameters(fcs[[1]])))
  
  ## Create dummy list to save all expression data
  exprTr_list<-c()
  
  ## Save arc sin transformed expression + spatial data for each flowFrame
  for(i in 1:length(md$file_name)){
    
    ## Expression data is combined with spatial parameter data
    
    exprRaw<-exprs(fcs_raw[[i]])
    
    colnames(exprRaw)<-panel_fcs$desc
    
    expr<-cbind(exprs(fcs[[i]])[, union(subtype_markers,functional_markers)],exprRaw[,otherparameters])
    
    ## Combine other (spatial) data with the protein data
    colnames(expr)<-c(colnames(exprs(fcs[[i]])),otherparameters)
    
    ## Filter out any event that is 95th percentile for BOTH CD29 and CD45 (antibody aggregates)
    ##expr<-expr[expr[,"CD29"] < quantile(expr[,"CD29"], probs=0.95) & expr[,"CD45"] < quantile(expr[,"CD45"], probs=0.95),]
    
    ## Add to list
    exprTr_list[[i]]<-expr
    
  }
  
  ## Create a flowSet based on the list of expression data
  fcs1<-flowSet(sapply(exprTr_list,flowFrame))
  
  ## Change parameter rownames
  panel_fcs1 <- pData(parameters(fcs1[[1]]))
  rownames(pData(parameters(fcs1[[1]]))) <- rownames(panel_fcs[panel_fcs$desc %in% pData(parameters(fcs1[[1]]))$name,])
  
  
  
  ###to scale every flowframe
  ## Save out the original rownames from the parameter list from the fcs flowFrame
  panellist <- rownames(pData(parameters(fcs[[1]])))
  
  ## Create dummy list to save all expression data
  exprTr_list<-c()
  
  ## Save arc sin transformed expression + spatial data for each flowFrame
  for(i in 1:length(md$file_name)){
    
    ## Expression data is combined with spatial parameter data
    
    expr<-exprs(fcs[[i]])
    
    expr<-t(scale(t(expr)))
    
    ## Add to list
    exprTr_list[[i]]<-expr
    
  }
  
  ## Create a flowSet based on the list of expression data
  fcs2<-flowSet(sapply(exprTr_list,flowFrame))
  
  
  
  
  ## Get sample ids
  sample_ids <- rep(md$sample_id, fsApply(fcs1, nrow))
  
  ## Return: 
  ## fcs (only marker expressions arcsin transformed), 
  ## fcs1 (arcsin transformed + spatial parameters), 
  ## fcs2 (scaled arcsin expr per flowframe)
  ## fcsraw (all raw data), and all marker/parameter lists
  return(list('fcs'=fcs,'fcs1'=fcs1,'fcs2'=fcs2,'fcsraw'=fcs_raw,'subtype_markers'=subtype_markers,'functional_markers'=functional_markers,'otherparameters'=otherparameters,'cluster_by'=cluster_by,'sample_ids'=sample_ids,'meta_data'=md))
}



clusterfcs <- function(fcs=output$fcs,
                       cluster_by = output$cluster_by,
                       seed=1234,plottitle='consensus_plots',
                       scaleoption, scaled.center, scaled.scale,
                       numclusters=40){
  ## Cell population identification with FlowSOM and ConsensusClusterPlus
  require(dplyr);require(FlowSOM);require(ConsensusClusterPlus)
  set.seed(seed)
  som <- ReadInput(fcs, 
                   transform = FALSE, 
                   scale = scaleoption,
                   scaled.center = scaled.center,
                   scaled.scale = scaled.scale) %>% BuildSOM(colsToUse = cluster_by)
  
  ## Get the cell clustering into 100 SOM codes
  cell_clustering_som <- som$map$mapping[,1]
  
  ## Metaclustering into numclusters with ConsensusClusterPlus
  codes <- som$map$codes
  mc <- ConsensusClusterPlus(t(codes), maxK = numclusters, reps = 100,
                             pItem = 0.9, pFeature = 1, title = plottitle, 
                             plot = "png", clusterAlg = "hc", 
                             innerLinkage = "average", finalLinkage = "average",
                             distance = "euclidean", seed = 1234)
  
  ## Get cluster ids for each cell
  code_clustering <- mc[[numclusters]]$consensusClass#metaclusters consensus
  cell_clustering <- code_clustering[cell_clustering_som]#cell clustering from som
  return(list('code_clustering'=code_clustering,'cell_clustering'=cell_clustering,'metaclusters'=mc))
}


####DIAGNOSTIC HEATMAP FUNCTIONS ####
plot_clustering_heatmap_wrapper <- function(fcs, cell_clustering, nclusters,
                                            number_clusters, cluster_merging = NULL, 
                                            cluster_by=output$cluster_by,
                                            clusterMergeFile=NULL,
                                            fileName = 'clusteringheatmap.pdf'){
  require(matrixStats);require(dplyr);require(RColorBrewer);require(pheatmap);require(readxl);require(flowCore);require(scales);require(pals);require(ComplexHeatmap)
  ## Will output the heatmap object and print it 
  color_clusters = kovesi.rainbow_bgyrm_35_85_c69(number_clusters)
  
  
  #get expression
  expr <- fsApply(fcs, exprs);expr <-expr[,cluster_by]
  ## Scale expression of all markers to values between 0 and 1
  rng <- colQuantiles(expr, probs = c(0.01, 0.99))
  expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0; expr01[expr01 > 1] <- 1;expr01 <-expr01[,cluster_by]
  
  
  ## Calculate the mean expression##################################################
  
  expr_mean <- data.frame(expr, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  expr01_mean <- data.frame(expr01, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  
  ## Calculate cluster frequencies
  
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_mean[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
  rownames(expr_heat) <- expr01_mean$cell_clustering
  
  ## Annotation for the original clusters
  
  annotation_row <- data.frame(Cluster = factor(expr01_mean$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  
  ## Annotation for the merged clusters
  
  if(!is.null(clusterMergeFile)){
    ifelse(grepl(clusterMergeFile,pattern='.xls'),cluster_merging <- read_excel(clusterMergeFile),cluster_merging <- read.csv(clusterMergeFile,header = TRUE))
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$Merged <- cluster_merging$new_cluster
    color_clusters2 <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters2) <- levels(cluster_merging$new_cluster)
    annotation_colors$Merged <- color_clusters2
  }
  
  ## Colors for the heatmap
  
  legend_breaks = seq(from = 0, to = 1, by = 0.2)
  
  clusternames<-sort(unique(cluster_merging$new_cluster))
  
  colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_merging$new_cluster)))
  
  names(colorassigned)<-clusternames
  
  color_list = list(clusters=colorassigned)
  
  color_list_byoriginal = colorassigned[match(cluster_merging$new_cluster,names(colorassigned))]
  
  cp<-rowAnnotation(clusters=cluster_merging$new_cluster,
                    col=color_list,
                    gp = gpar(col = "white", lwd = .5),
                    prop=anno_barplot(
                      clustering_prop, 
                      gp = gpar(fill=color_list_byoriginal, col=F),
                      border = F,
                      bar_width = 0.75, 
                      width = unit(2,"cm")))
  
  q <- Heatmap(expr_heat, name="scaled",
               col=rev(brewer.rdbu(100)),
               row_order = cluster_merging[order(cluster_merging$new_cluster),]$original_cluster,
               cluster_columns = T,
               cluster_rows = T,
               border = NA,
               rect_gp = gpar(col = "white", lwd = .5),
               right_annotation = cp,
               show_row_names = T,
               row_names_gp = gpar(fontsize=7),
               column_names_gp = gpar(fontsize=10),
               heatmap_legend_param = list(at=seq(from = 0, to = 1, by = 0.2)),
               width = unit(10, "cm"))
  
  print('Colors:')
  print(color_clusters)
  
  pdf(fileName, width=8, height=6) 
  
  return(q)
  
  dev.off() 
  
}


plot_clustering_heatmap_wrapper2 <- function(fcs, cell_clustering, nclusters=40,
                                             color_clusters=colorassigned,
                                             colorbar=rev(brewer.rdbu(100)),
                                             subtype_markers,
                                             fileName = 'clusteringheatmap.pdf'){
  require(matrixStats);require(dplyr);require(RColorBrewer);require(pheatmap);require(readxl);require(flowCore);require(scales)
  ## Will output the heatmap object and print it 
  #get expression
  expr <- fsApply(fcs, exprs);expr <-expr[,subtype_markers]
  ## Scale expression of all markers to values between 0 and 1
  rng <- colQuantiles(expr, probs = c(0.01, 0.99))
  expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0; expr01[expr01 > 1] <- 1;expr01 <-expr01[,subtype_markers]
  ## Calculate the mean expression##################################################
  pdf(fileName, width=8, height=11) 
  expr_mean <- data.frame(expr, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  expr01_mean <- data.frame(expr01, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  
  ## Calculate cluster frequencies
  
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_mean[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
  rownames(expr_heat) <- expr01_mean$cell_clustering
  
  ## Colors for the heatmap
  
  #legend_breaks = seq(from = 0, to = 1, by = 0.2)
  #labels_row <- expr01_mean$cell_clustering
  
  labels_row <- paste0(expr01_mean$cell_clustering, " ")
  
  ## Annotation for the original clusters
  
  annotation_row <- data.frame(Cluster = factor(expr01_mean$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  
  p <- pheatmap(expr_heat, 
                color = colorbar, 
                cluster_cols = F,
                cluster_rows = F, 
                labels_row = labels_row,
                #scale="row",
                display_numbers = F, 
                number_color = "black",
                fontsize = 9, fontsize_number = 6,  
                #legend_breaks = legend_breaks,
                annotation_row = annotation_row, 
                annotation_colors = annotation_colors,
                cellwidth = 8,
                cellheight = 8,
                border_color = "white",
                annotation_legend = F
  )
  dev.off() 
  print('Colors:')
  print(color_clusters)
  print(p);return(p)
}


plot_clustering_heatmap_wrapper3 <- function(fcs, cell_clustering, nclusters,
                                             cluster_by=output$cluster_by,
                                             colorassigned,
                                             sampleno,
                                             clusterMergeFile=NULL,
                                             fileName = 'clusteringheatmap.pdf'){
  require(matrixStats);require(dplyr);require(RColorBrewer);require(pheatmap);require(readxl);require(flowCore);require(scales);require(pals);require(ComplexHeatmap)
  ## Will output the heatmap object and print it 
  
  number_clusters <- length(unique(cell_clustering[which(output$sample_ids==output$meta_data$sample_id[sampleno])]))
  
  color_clusters = kovesi.rainbow_bgyrm_35_85_c69(number_clusters)
  
  #get expression
  expr <- exprs(fcs);expr <-expr[,cluster_by]
  ## Scale expression of all markers to values between 0 and 1
  rng <- colQuantiles(expr, probs = c(0.01, 0.99))
  expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0; expr01[expr01 > 1] <- 1;expr01 <-expr01[,cluster_by]
  
  ## Calculate the mean expression##################################################
  
  cell_clustering <- cell_clustering[which(output$sample_ids==output$meta_data$sample_id[sampleno])]
  
  expr_mean <- data.frame(expr, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  expr01_mean <- data.frame(expr01, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  
  ## Calculate cluster frequencies
  
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_mean[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
  rownames(expr_heat) <- expr01_mean$cell_clustering
  
  ## Annotation for the original clusters
  
  annotation_row <- data.frame(Cluster = factor(expr01_mean$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  
  ## Annotation for the merged clusters
  
  if(!is.null(clusterMergeFile)){
    ifelse(grepl(clusterMergeFile,pattern='.xls'),cluster_merging1 <- read_excel(clusterMergeFile),cluster_merging1 <- read.csv(clusterMergeFile,header = TRUE))
    cluster_merging1 <- cluster_merging1[cluster_merging1$original_cluster %in% unique(cell_clustering),]
    cluster_merging1$new_cluster <- factor(cluster_merging1$new_cluster)
    annotation_row$Merged <- cluster_merging1$new_cluster
    color_clusters2 <- colorassigned[c(annotation_row$Merged)]
    names(color_clusters2) <- levels(cluster_merging1$new_cluster)
    annotation_colors$Merged <- color_clusters2
  }
  
  ## Colors for the heatmap
  
  legend_breaks = seq(from = 0, to = 1, by = 0.2)
  
  clusternames<-sort(unique(cluster_merging1$new_cluster))
  
  colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_merging1$new_cluster)))
  
  names(colorassigned)<-clusternames
  
  color_list = list(clusters=colorassigned)
  
  color_list_byoriginal = colorassigned[match(cluster_merging1$new_cluster,names(colorassigned))]
  
  cp<-rowAnnotation(clusters=cluster_merging1$new_cluster,
                    col=color_list,
                    gp = gpar(col = "white", lwd = .5),
                    prop=anno_barplot(
                      clustering_prop, 
                      gp = gpar(fill=color_list_byoriginal, col=F),
                      border = F,
                      bar_width = 0.75, 
                      width = unit(2,"cm")))
  
  q <- Heatmap(expr_heat,
               right_annotation = cp)
  
  r <- Heatmap(expr_heat, name="scaled",
               col=rev(brewer.rdbu(100)),
               #row_order = cluster_merging[order(cluster_merging$new_cluster),]$original_cluster,
               cluster_columns = T,
               cluster_rows = T,
               border = NA,
               rect_gp = gpar(col = "white", lwd = .5),
               right_annotation = cp,
               show_row_names = T,
               row_names_gp = gpar(fontsize=7),
               column_names_gp = gpar(fontsize=10),
               heatmap_legend_param = list(at=seq(from = 0, to = 1, by = 0.2)),
               width = unit(10, "cm"))
  
  print('Colors:')
  print(color_clusters)
  
  pdf(fileName, width=8, height=6) 
  
  return(q)
  
  dev.off() 
  
}



plot_clustering_heatmap_wrapper4 <- function(fcs, cell_clustering, nclusters=40,
                                             color_clusters=colorassigned,
                                             colorbar=rev(brewer.rdbu(100)),
                                             rowcut=rowcut,
                                             colcut=colcut,
                                             subtype_markers,
                                             fileName = 'clusteringheatmap.pdf'){
  require(matrixStats);require(dplyr);require(RColorBrewer);require(pheatmap);require(readxl);require(flowCore);require(scales)
  ## Will output the heatmap object and print it 
  #get expression
  expr <- fsApply(fcs, exprs);expr <-expr[,subtype_markers]
  ## Scale expression of all markers to values between 0 and 1
  rng <- colQuantiles(expr, probs = c(0.1, 0.9))
  expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0; expr01[expr01 > 1] <- 1;expr01 <-expr01[,subtype_markers]
  ## Calculate the mean expression##################################################
  pdf(fileName, width=8, height=11) 
  expr_mean <- data.frame(expr, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  expr01_mean <- data.frame(expr01, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  
  ## Calculate cluster frequencies
  
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_mean[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
  rownames(expr_heat) <- expr01_mean$cell_clustering
  
  ## Colors for the heatmap
  
  #legend_breaks = seq(from = 0, to = 1, by = 0.2)
  #labels_row <- expr01_mean$cell_clustering
  
  labels_row <- paste0(expr01_mean$cell_clustering, " ")
  
  ## Annotation for the original clusters
  
  annotation_row <- data.frame(Cluster = factor(expr01_mean$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  
  p <- pheatmap(expr_heat, 
                color = colorbar, 
                cluster_cols = F,
                cluster_rows = F, 
                labels_row = labels_row,
                #scale="row",
                display_numbers = F, 
                number_color = "black",
                fontsize = 9, fontsize_number = 6,  
                #legend_breaks = legend_breaks,
                annotation_row = annotation_row, 
                annotation_colors = annotation_colors,
                cellwidth = 10,
                cellheight = 10,
                border_color = NA,
                annotation_legend = F,
                gaps_col = colcut,
                gaps_row = rowcut
  )
  dev.off() 
  print('Colors:')
  print(color_clusters)
  print(p);return(p)
}



####REQUIRED LIBRARIES####
library(reshape2)
library(pals)
library(ggplot2)
library(Hmisc)
library(ComplexHeatmap)
library(ggiraphExtra)
library(ggvoronoi)
library(ggpubr)
library(multcomp)
library(sf)
library(clusterSim)
library(circlize)
library(RColorBrewer)
library(stringr)
library(igraph)
library(readxl)
library(dplyr)
library(packcircles)
library(gridExtra)
library(limma)
library(qgraph)
library(flowCore)
library(basetheme)
library(pheatmap)
library(premessa)
library(spatstat)


#======================
####RUNNING DATA####
#======================


####DATA LOADING####

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
work<-getwd()



## If loading from prior run
output<-readRDS('backup_output.rds')


## Read (skip if previously ran)

output <- returnfcs(metaDataFile = paste0(work,"/Config/J19113_metadata.xlsx"),
                    panelDataFile = paste0(work,"/Config/J19113_panel.xlsx"),
                    dataDirectory = paste0(work,"/Data"))


## Set up levels

samplevels=c("P01_t1_1","P01_t1_2","P01_t3_1","P01_t3_2","P01_t3_3",
             "P04_t3_1","P04_t3_3",
             "P05_t1_1","P05_t3_1","P05_t3_2","P05_t3_3",
             "P06_t1_1","P06_t1_2","P06_t3_4","P06_t3_5","P06_t3_6",
             "P07_t1_1","P07_t1_2","P07_t1_3","P07_t3_1","P07_t3_2","P07_t3_3",
             "P08_t1_1","P08_t1_2","P08_t1_3","P08_t1_4","P08_t3_1","P08_t3_2",
             "P09_t1_1","P09_t1_2","P09_t1_3","P09_t1_4","P09_t2_1","P09_t2_2","P09_t2_3",
             "P010_t1_1","P010_t2_1","P010_t2_2","P010_t2_3",
             "P011_t1_3","P011_t2_1","P011_t2_3",
             "P012_t1_1","P012_t1_2","P012_t1_3",
             "P013_t1_1","P013_t1_2","P013_t1_3","P013_t2_1","P013_t2_2",
             "P015_t3_4","P015_t3_4_1","P015_t3_4_2","P015_t3_5","P015_t3_6",
             "P018_t1_1","P018_t1_2","P018_t1_3","P018_t3_1",
             "P019_t1_1","P019_t1_2","P019_t3_1","P019_t3_2")

timegrouplevels=c("PreTx","OnTx")

timepointlevels=c("Baseline","Cycle1","Cycle2")

tumorlevels=c("LUNG","LIVER", "OTHER","LIVER_EXC")

caselevels=c("P01", "P04", "P05", "P06", "P07", "P08", "P09", "P010", "P011", "P012", "P013", "P015", "P018", "P019")

####DIAGNOSTICS####
## Spot check - number of cells per sample
cell_table <- table(output$sample_ids)
ggdf <- data.frame(sample_id = names(cell_table), 
                   cell_counts = as.numeric(cell_table))
ggdf$case <- factor(output$meta_data$Case[match(ggdf$sample_id,output$meta_data$sample_id)], levels = caselevels)
ggdf$tumor <- factor(output$meta_data$Tumor[match(ggdf$sample_id,output$meta_data$sample_id)], levels=tumorlevels)
ggdf$timepoint <- factor(output$meta_data$Timepoint[match(ggdf$sample_id,output$meta_data$sample_id)], levels=timepointlevels)
ggdf$timegroup <- factor(output$meta_data$TimeGroup[match(ggdf$sample_id,output$meta_data$sample_id)], levels=timegrouplevels)

ggp<-ggplot(ggdf, aes(x = sample_id, y = cell_counts, fill = sample_id)) + 
  geom_bar(stat = 'identity') + 
  #geom_text(aes(label = cell_counts), angle = 45, hjust = 0.5, vjust = -0.5, size = 2) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=5)) +  
  scale_x_discrete(drop = FALSE)
pdf('Diagnostics_cellcounts.pdf',width=6, height=4);ggp; dev.off()

## Multi-dimensional scaling plot to show similarities between samples
## Get the mean marker expression per sample
expr_mean_sample_tbl <- data.frame(sample_id = output$sample_ids, fsApply(output$fcs,exprs)) %>%
  group_by(sample_id) %>%  summarize_all(funs(mean))
expr_mean_sample <- t(expr_mean_sample_tbl[, -1])
colnames(expr_mean_sample) <- expr_mean_sample_tbl$sample_id
mds <- plotMDS(expr_mean_sample, plot = FALSE)
ggdf <- data.frame(MDS1 = mds$x, MDS2 = mds$y,
                   sample_id = colnames(expr_mean_sample))
ggdf$case <- factor(output$meta_data$Case[match(ggdf$sample_id,output$meta_data$sample_id)], levels = caselevels)
ggdf$tumor <- factor(output$meta_data$Tumor[match(ggdf$sample_id,output$meta_data$sample_id)], levels=tumorlevels)
ggdf$timepoint <- factor(output$meta_data$Timepoint[match(ggdf$sample_id,output$meta_data$sample_id)], levels=timepointlevels)
ggdf$timegroup <- factor(output$meta_data$TimeGroup[match(ggdf$sample_id,output$meta_data$sample_id)], levels=timegrouplevels)
ggp<-ggplot(ggdf, aes(x = MDS1, y = MDS2, color = tumor, shape = timegroup)) +
  geom_point(size = 2.5) +
  #geom_text(aes(label = patient_id)) +
  theme_bw()+
  theme(plot.background = element_rect(fill="black"),
        panel.background = element_rect(fill="black"),
        panel.grid = element_blank(),
        axis.line = element_line(color="white"),
        axis.text = element_text(color="white"),
        axis.title = element_text(color="white"),
        legend.background = element_rect(fill="black"),
        legend.key = element_rect(fill="white"),
        legend.text = element_text(color="white"))
pdf('Diagnostics_MDS_sample.pdf',width=6, height=6);ggp; dev.off()

## Colors for the heatmap
color <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)
pdf('Diagnostics_Heatmap.pdf',width=8, height=8)
pheatmap(expr_mean_sample_tbl[,c(2:36,39)], color = color, display_numbers = TRUE,
         number_color = "black", fontsize_number = 3, 
         clustering_method = "average")
dev.off()


####CLUSTERING####

##Revised loading depending on the diagnostics if needed

output <- returnfcs(metaDataFile = paste0(work,"/Config/J19113_metadata.xlsx"),
                    panelDataFile = paste0(work,"/Config/J19113_panel.xlsx"),
                    dataDirectory = paste0(work,"/Data"))


##Clustering

output[(length(output)+1):(length(output)+3)] <- clusterfcs(fcs=output$fcs, numclusters=50, scaleoption = T, scaled.center = T, scaled.scale = T) 
#output$fcs uses just arcsin transformed data
#scaleoption scales the dataset across the channels so that the channels with the highest intensities will not dominate the clustering

names(output)[(length(output)-2):(length(output))] <- c('code_clustering','cell_clustering','metaclusters')


####ANNOTATIONS AND VISUALIZATION OF CLUSTERS####

## Load merge file
## Assign colors
clusterMergeFile = paste0(work,"/Config/J19113_merge.xlsx") #create dummy merger numbers prior to annotation
cluster_merging <- read_excel(clusterMergeFile)

clusterlevels=c(1:50)

colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_merging$new_cluster)))

clusternames<-clusterlevels
names(colorassigned)<-clusternames
mm1 <- match(output$cell_clustering, cluster_merging$original_cluster)
cell_clustering1m <- cluster_merging$new_cluster[mm1]
output$cell_clustering1m <- cell_clustering1m

## Metacluster heatmaps
plot_clustering_heatmap_wrapper(fcs=output$fcs,
                                number_clusters = length(unique(output$cell_clustering1m)),
                                cell_clustering = output$cell_clustering, 
                                cluster_by = output$cluster_by,
                                clusterMergeFile = clusterMergeFile,
                                fileName = 'Clusteringheatmap_all_fcs.pdf'); dev.off()

clusterMergeFile = paste0(work,"/Config/J19113_merge.xlsx") #create dummy merger numbers prior to annotation
cluster_merging <- read_excel(clusterMergeFile)

clusterlevels=c("B",
                "DC",
                "Gran",
                "Mac_I",
                "Mac_II",
                "Mac_III",
                "Mac_IV",
                "Myeloid",
                "NK",
                "Stroma_I",
                "Stroma_II",
                "Tc",
                "Th",
                "Treg",
                "Tumor",
                "UNS")

colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_merging$new_cluster)))
clusternames<-clusterlevels
names(colorassigned)<-clusternames
mm1 <- match(output$cell_clustering, cluster_merging$original_cluster)
cell_clustering1m <- cluster_merging$new_cluster[mm1]
output$cell_clustering1m <- cell_clustering1m

plot_clustering_heatmap_wrapper2(fcs=output$fcs,
                                 colorbar = kovesi.diverging_bwr_40_95_c42(100),
                                 subtype_markers = c("CD45", "CD45RA", "CD45RO","DCLAMP", "FOXP3", "CD3", "CD4", "CD8", "GZMB", "CD20", "CD57", "CD15", "CD14", "CD16", "CD68", "CD163", "CD206", "SMA", "COL", "PDPN", "PDGFRA", "PanCK", "ECAD"),
                                 color_clusters = colorassigned,
                                 cell_clustering = factor(output$cell_clustering1m, levels=clusterlevels), 
                                 fileName = 'Clusteringheatmap_merged_fcs.pdf');dev.off()


####for SU2C presentation####
#to highlight CXCL12
clusterlevels2=c("DC","Gran", #2
                 "Mac_I","Mac_II","Mac_III","Mac_IV","Myeloid", #7
                 "B","Tc","Th","Treg", #11
                 "Stroma_I","Stroma_II", #13
                 "Tumor", #14
                 "UNS", "NK")

wheretorowcut=c(2,7,11,13,14) #for cluster gaps

wheretocolcut=c(6,13,17,19) #for marker gaps

#heatmap scaled 10 to 90
plot_clustering_heatmap_wrapper4(fcs=output$fcs2,
                                 colorbar = rev(brewer.rdylbu(100)),
                                 subtype_markers = c("FOXP3", "CD3", "CD4", "CD8", "GZMB", "CD20",#6
                                                     "DCLAMP","CD15", "CD14", "CD16", "CD68", "CD163", "CD206", #13
                                                     "SMA", "COL", "PDPN", "PDGFRA", #17
                                                     "PanCK", "ECAD", #19
                                                     "CXCL12"),
                                 color_clusters = colorassigned,
                                 rowcut = wheretorowcut,
                                 colcut = wheretocolcut,
                                 cell_clustering = factor(output$cell_clustering1m, levels=clusterlevels2), 
                                 fileName = 'Clusteringheatmap_merged_forfigure1.pdf');dev.off()


## Remove any samples after clustering
samplestoexclude=c("P015_t3_4_1","P015_t3_4_2") #duplicate with "P015_t3_4"; 
include_inds <- which(output$sample_ids %nin% samplestoexclude)

## Save output list
saveRDS(output, file="backup_output.rds")


####ABUNDANCE PLOTS####
## Exclude non-liver samples
nonliversamples=c("P01_t1_1",
                  "P01_t1_2",
                  "P01_t3_1",
                  "P01_t3_2",
                  "P01_t3_3",
                  "P06_t1_1",
                  "P06_t1_2",
                  "P06_t3_4",
                  "P06_t3_5",
                  "P06_t3_6",
                  "P08_t1_1",
                  "P08_t1_2",
                  "P08_t1_3",
                  "P08_t1_4",
                  "P08_t3_1",
                  "P08_t3_2",
                  "P010_t1_1",
                  "P010_t2_1",
                  "P010_t2_2",
                  "P010_t2_3") #P01 is liver but was excluded due to poor staining

## Proportion calculations
counts_table <- table(output$cell_clustering1m[include_inds], output$sample_ids[include_inds])
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
counts <- as.data.frame.matrix(counts_table)
props <- as.data.frame.matrix(props_table)

#Density calculations
areas <- tibble(sample_id=colnames(counts))
raw_areas <- read_xlsx(paste0(work,'/Config/J19113_area.xlsx'))
areas$TotalArea <- raw_areas$TotalArea[match(areas$sample_id, raw_areas$sample_id)]
densities <- t(t(counts)/areas$TotalArea)

write.csv(counts,'Results_counts.csv')
write.csv(props,'Results_props.csv')
write.csv(densities, 'Results_densities.csv')

## Set up the data frame for proportional plotting
ggdf <- melt(data.frame(cluster = rownames(props),props, check.names = FALSE),
             id.vars = "cluster", value.name = "proportion", 
             variable.name = "sample_id")
ggdf$sample_id <- factor(ggdf$sample_id, levels=samplevels)
ggdf$case <- factor(output$meta_data$Case[match(ggdf$sample_id,output$meta_data$sample_id)], levels = caselevels)
ggdf$tumor <- factor(output$meta_data$Tumor[match(ggdf$sample_id,output$meta_data$sample_id)], levels=tumorlevels)
ggdf$timepoint <- factor(output$meta_data$Timepoint[match(ggdf$sample_id,output$meta_data$sample_id)], levels=timepointlevels)
ggdf$timegroup <- factor(output$meta_data$TimeGroup[match(ggdf$sample_id,output$meta_data$sample_id)], levels=timegrouplevels)


## Set up dataframe for densities
ggdfd <- melt(data.frame(cluster = rownames(densities),densities, check.names = FALSE),
              id.vars = "cluster", value.name = "densities", 
              variable.name = "sample_id")
ggdfd$sample_id <- factor(ggdfd$sample_id, levels=samplevels)
ggdfd$case <- factor(output$meta_data$Case[match(ggdfd$sample_id,output$meta_data$sample_id)], levels = caselevels)
ggdfd$tumor <- factor(output$meta_data$Tumor[match(ggdfd$sample_id,output$meta_data$sample_id)], levels=tumorlevels)
ggdfd$timepoint <- factor(output$meta_data$Timepoint[match(ggdfd$sample_id,output$meta_data$sample_id)], levels=timepointlevels)
ggdfd$timegroup <- factor(output$meta_data$TimeGroup[match(ggdfd$sample_id,output$meta_data$sample_id)], levels=timegrouplevels)


## Set up ggdf's to include liver samples ONLY - proportion dataframe and density dataframe respectively
ggdf_liver<-ggdf[ggdf$sample_id %nin% nonliversamples,]
ggdfd_liver<-ggdfd[ggdfd$sample_id %nin% nonliversamples,]


#plot box plots

#% CELLS
ggp2<-ggplot(ggdf_liver,aes(x=case,y=proportion,fill=case,shape=timepoint))+
  geom_boxplot(outlier.shape=NA, lwd=0.25, color="black")+
  geom_jitter(width=0.2, size=2)+
  facet_wrap(~cluster,ncol=7,scales="free")+
  ylab("% of Cells")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=12, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=7, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white")  
  )

pdf('Abundance_box2_merged.pdf',width=32,height=16)
ggp2

dev.off()

#for densities

ggp2<-ggplot(ggdfd_liver,aes(x=case,y=densities,fill=case))+
  geom_boxplot(outlier.shape=NA, lwd=0.25, color="black")+
  geom_jitter(width=0.2, size=1)+
  facet_wrap(~cluster,ncol=4,scales="free")+
  ylab("# of Cells / mm^2")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=12, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=12, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=11),
        legend.key = element_rect(fill="white")  
  )

pdf('Abundance_box2_densities_merged.pdf',width=12,height=9)
ggp2

dev.off()


## Stacked bars
bp <- ggplot(ggdf, aes(x = sample_id, y = proportion, fill=cluster, order=cluster)) +
  geom_bar(stat = "identity", position="fill", width=0.8) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, color="black", size=6),
        axis.text.y = element_text(color="black"),
        axis.ticks = element_line(size=0.25),
        strip.text=element_text(size=8),
        strip.background = element_rect(fill=NA, color=NA),
        legend.title = element_blank(),
        legend.key.size = unit(.75,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white"),
  ) +
  ylab("% of Total Cells")+
  scale_fill_manual(values = colorassigned,
                    breaks = clusterlevels,
                    labels = clusterlevels)+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(ncol=2))
pdf('Abundance_stackedbar.pdf', width=6, height=6); bp; dev.off()


####LINE PLOTS##############
## Line plots

## Line Plot - Mean Case Density
linegraph_data <- ggdfd_liver %>% group_by(cluster, case, tumor, timepoint, timegroup) %>% summarize_at(vars(densities), funs(mean))
## Line Plot - subset for cases with matched timepoints
linegraph_data_pairedonly <-linegraph_data[linegraph_data$case %nin% c("P04","P012","P015"),]
## Line Plot - include a merged "Mac" cluster w/ all Mac subsets
linegraph_data_pairedonly_mac <- linegraph_data_pairedonly[linegraph_data_pairedonly$cluster %in% c("Mac_I","Mac_II","Mac_III","Mac_IV"),]
linegraph_data_cleaned <- linegraph_data_pairedonly[,c(1,2,5,6)]
colorassigned2<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(ggdfd_liver$case)))

#merged mac clusters line plot
merged_mac_lp <- aggregate(linegraph_data_pairedonly_mac$densities, by=list(case=linegraph_data_pairedonly_mac$case, timegroup=linegraph_data_pairedonly_mac$timegroup), FUN=sum)
colnames(merged_mac_lp)[3] <- "densities"
merged_mac_lp <-cbind(cluster = "Mac", merged_mac_lp)

linegraph_data_cleaned <- rbind(linegraph_data_cleaned, merged_mac_lp)


lp <- ggplot (data=linegraph_data_cleaned, aes(x=timegroup, y = densities, group = case, colour = case)) +
  geom_line(size = 0.75)+
  geom_point(size = 2) +
  facet_wrap(~cluster, scales = "free")+
  labs(
    x = "Treatment Time Group",
    y = "# over Tissue Area (mm²)",
    title = "Cluster Line Plots",
    group = "Case"
  )

## Line Plot Aesthetics

lp_styled <- lp + theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 18, colour = "black"),
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    panel.grid.major = element_line(colour = "grey70", size = 0.2),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(), #element_rect(fill = "black", color = "black", size = 1),
    strip.text = element_text(face = "bold", size = 12,color = "black"),
    legend.background = element_rect(fill = "white", colour = "white"),
  )
pdf('Line_Plot_merged.pdf',width=9.5,height=9);lp_styled;
dev.off()


lp <- ggplot (data=linegraph_data_cleaned, aes(x=timegroup, y = densities, group = case, colour = case)) +
  geom_line(size = 0.75)+
  geom_point(size = 2) +
  facet_wrap(~cluster, scales = "free")+
  labs(
    x = "Treatment Time Group",
    y = "# over Tissue Area (mm²)",
    title = "Cluster Line Plots",
    group = "Case"
  ) +
  scale_color_manual(name = "Case", values = colorassigned2)+
  scale_y_sqrt(expand = expansion(mult = c(0.05, 0.15)))+
  stat_compare_means(comparisons = list(c("PreTx", "OnTx")),
                     paired=T,
                     method = "t.test",
                     label = "p.format",
                     label.x = 1.4,
                     label.y.npc="top",
                     hide.ns = T,
                     size=5)

## Line Plot Aesthetics

lp_styled <- lp + theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 18, colour = "black"),
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    panel.grid.major = element_line(colour = "grey70", size = 0.2),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(), #element_rect(fill = "black", color = "black", size = 1),
    strip.text = element_text(face = "bold", size = 12,color = "black"),
    legend.background = element_rect(fill = "white", colour = "white"),
  )
pdf('Line_Plot_merged_statspvalue.pdf',width=9.5,height=9);lp_styled;
dev.off()


##################################

#cluster heatmap for macrophage markers and clusters only
markerlist_mac = c("CD163","CD206","HLADR","CD86")
exprtbl <- data.frame(fsApply(output$fcs,exprs)[, markerlist_mac], sample_id = output$sample_ids, cluster = output$cell_clustering1m)
exprtbl <- exprtbl[exprtbl$sample_id %nin% nonliversamples,]
exprtbl <- exprtbl %>% group_by(cluster) %>% summarise_at(markerlist_mac, mean, na.rm=TRUE)
exprmtx <- as.matrix(exprtbl[,2:(length(markerlist_mac)+1)])
rownames(exprmtx) <- unlist(exprtbl[,1])

pdf("Clusterheatmap_focused_Mac.pdf", width=5, height=5)
pheatmap(exprmtx[clusterlevels[str_detect(clusterlevels, "Mac")],], 
         scale="column",
         cluster_rows = F,
         cellwidth=10,
         cellheight=10)
dev.off()

####DISTANCE RELATIONSHIP####


##colors
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(output$cell_clustering1m)))
hex <- hue_pal()(9)
colorassignedbroad <- c(rep(hex[1],1), #B - lym
                        rep(hex[2],2), #DC, Gran - myl
                        rep(hex[2],4), #Macs - myl
                        rep(hex[2],1), #Myeloid - myl
                        rep(hex[1],1), #NK - lym
                        rep(hex[3],2), #Str
                        rep(hex[1],3), #T cells - lym
                        rep(hex[4],1), #tumor
                        rep(hex[5],1)) #other
clusternames<-clusterlevels
names(colorassigned)<-clusternames

#cell type legend 
allcelltypes<-clusterlevels
legendctype<-as.data.frame(cbind(paste0("ctype",1:length(allcelltypes)),allcelltypes))
legendctype$maintype<-1
legendctype$maintype[str_detect(legendctype$allcelltypes,"T")]<-"Lym"
legendctype$maintype[str_detect(legendctype$allcelltypes,"B")]<-"Lym"
legendctype$maintype[str_detect(legendctype$allcelltypes,"NK")]<-"Lym"
legendctype$maintype[str_detect(legendctype$allcelltypes,"DC")]<-"Myl"
legendctype$maintype[str_detect(legendctype$allcelltypes,"M")]<-"Myl"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Gran")]<-"Myl"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Str")]<-"Stroma"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Tumor")]<-"Tumor"
legendctype$maintype[str_detect(legendctype$allcelltypes,"UNS")]<-"Other"

####DISTANCE RELATIONSHIPS####

counts_table <- table(output$cell_clustering1m, output$sample_ids)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
counts <- as.data.frame.matrix(counts_table)
props <- as.data.frame.matrix(props_table)

ggdft <- melt(data.frame(cluster = rownames(counts), counts, check.names = FALSE),
              id.vars = "cluster", value.name = "counts", 
              variable.name = "sample_id")
ggdft$sample_id <- factor(ggdft$sample_id, levels=samplevels)
ggdft$cluster <- factor(ggdft$cluster, levels=clusterlevels)
ggdft$case <- factor(output$meta_data$Case[match(ggdft$sample_id,output$meta_data$sample_id)], levels=caselevels)
ggdft$tumor <- factor(output$meta_data$Tumor[match(ggdft$sample_id,output$meta_data$sample_id)], levels=tumorlevels)
ggdft$timepoint <- factor(output$meta_data$Timepoint[match(ggdft$sample_id,output$meta_data$sample_id)], levels=timepointlevels)
ggdft$timegroup <- factor(output$meta_data$TimeGroup[match(ggdft$sample_id,output$meta_data$sample_id)], levels=timegrouplevels)

ggdft_liver <- ggdft[ggdft$tumor=="LIVER",]

totalcounts<- ggdft_liver %>% group_by(cluster, case, tumor, timegroup) %>% summarize_at(vars(counts),funs(sum))

write.csv(totalcounts,"Totalcounts.csv")

totalcounts_timegrouptype <- ggdft_liver %>% group_by(cluster, timegroup) %>% summarize_at(vars(counts),funs(sum))

write.csv(totalcounts_timegrouptype, "Totalcounts_timegrouptype.csv")

#percentage of each respective total

totalcounts_pretx <- totalcounts_timegrouptype[totalcounts_timegrouptype$timegroup=="PreTx",]
totalcounts_ontx <- totalcounts_timegrouptype[totalcounts_timegrouptype$timegroup=="OnTx",]

totalpretx<-sum(totalcounts_pretx$counts)
totalontx<-sum(totalcounts_ontx$counts)

pct_pretx <- totalcounts_pretx$counts/totalpretx*100
names(pct_pretx)<- totalcounts_pretx$cluster
pct_ontx <- totalcounts_ontx$counts/totalontx*100
names(pct_ontx)<- totalcounts_ontx$cluster


ggdf_pctpretx <- melt(pct_pretx);ggdf_pctpretx$timegroup<-"PreTx"
ggdf_pctontx <- melt(pct_ontx);ggdf_pctontx$timegroup<-"OnTx"


ggdf_pct<-rbind(ggdf_pctpretx,ggdf_pctontx)
ggdf_pct$cluster<-c(rownames(ggdf_pctpretx),rownames(ggdf_pctontx))
rownames(ggdf_pct)<-1:nrow(ggdf_pct)
ggdf_pct$timegroup<-factor(ggdf_pct$timegroup, levels=c("PreTx","OnTx"))

bp <- ggplot(ggdf_pct, aes(x = timegroup, y = value, fill=cluster, order=cluster)) +
  geom_bar(stat = "identity", position="fill", width=0.8) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, color="black", size=6),
        axis.text.y = element_text(color="black"),
        axis.ticks = element_line(size=0.25),
        strip.text=element_text(size=8),
        strip.background = element_rect(fill=NA, color=NA),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(.75,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white"),
  ) +
  ylab("% of Total Cells")+
  scale_fill_manual(values = colorassigned,
                    breaks = clusterlevels,
                    labels = clusterlevels)+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(ncol=2))
pdf('Abundance_stackedbar_timegroup_liveronly.pdf', width=3, height=4); bp; dev.off()


#identify which cell types in either of the data subsets are very rare
#also include UNS cluster = ctype16
exclude_pretx<-which(pct_pretx<5/totalpretx*100)
exclude_pretx<-legendctype$V1[match(names(exclude_pretx), legendctype$allcelltypes)]
exclude_pretx<-union(exclude_pretx, "ctype16")
exclude_ontx<-which(pct_ontx<65/totalontx*100)
exclude_ontx<-legendctype$V1[match(names(exclude_ontx), legendctype$allcelltypes)]
exclude_ontx<-union(exclude_ontx, "ctype16")



#get out expression levels, X, and Y coords

expr <- fsApply(output$fcs1, exprs) #create expression matrix

expr0<-data.frame(expr[,c(union(output$subtype_markers,output$functional_markers),"Object Id","XMin","YMin")],
                  cluster=output$cell_clustering1m,
                  sample_id=output$sample_ids)

expr0$case <- factor(output$meta_data$Case[match(expr0$sample_id,output$meta_data$sample_id)], levels=caselevels)
expr0$tumor <- factor(output$meta_data$Tumor[match(expr0$sample_id,output$meta_data$sample_id)], levels=tumorlevels)
expr0$timepoint <- factor(output$meta_data$Timepoint[match(expr0$sample_id,output$meta_data$sample_id)], levels=timepointlevels)
expr0$timegroup <- factor(output$meta_data$TimeGroup[match(expr0$sample_id,output$meta_data$sample_id)], levels=timegrouplevels)


expr0<-expr0[expr0$cluster!="UNS",]
expr0_liver<- expr0[expr0$tumor=="LIVER",]
expr0_pre <-expr0_liver[expr0_liver$timegroup=="PreTx",]
expr0_on <-expr0_liver[expr0_liver$timegroup=="OnTx",]



###CREATE DISTANCE MATRICES FOR PROGRESSION CRITERIA IN LUM or TNC ONLY

PreTxid<-unique(expr0_pre$sample_id)
OnTxid<-unique(expr0_on$sample_id)

##Liver Pre-Tx ONLY====
expr_pre<-c()

for(k in 1:length(PreTxid)){
  expr_k<-expr0_pre[expr0_pre$sample_id==PreTxid[k],] 
  
  #create placer cols for shortest distance to each cell type
  dummy<-matrix(nrow=nrow(expr_k),ncol=length(allcelltypes)) 
  colnames(dummy)<-legendctype$allcelltypes
  dummy<-as.data.frame(dummy)
  expr_k<-data.frame(expr_k,dummy)
  
  X <- ppp(x=expr_k$XMin,y=expr_k$YMin,marks = as.factor(expr_k$cluster), window = owin(xrange = range(expr_k$XMin), yrange = range(expr_k$YMin)),checkdup = F)
  distByCelltype <- nndist(X, by=marks(X))
  
  for (i in 1:(ncol(distByCelltype))){
    d <- which(colnames(expr_k) == (colnames(distByCelltype)[i]))
    expr_k[d] <- distByCelltype[,i]
  }
  expr_pre <- rbind(expr_pre, expr_k)
}
colnames(expr_pre)<-c(colnames(expr0_pre),legendctype$V1)
expr_pre_m<-as.matrix(expr_pre[,colnames(expr_pre)[str_detect(colnames(expr_pre),"ctype")]])
expr_pre$ctype_no<-paste0("ctype",match(expr_pre$cluster,legendctype$allcelltypes))
rownames(expr_pre_m)<-expr_pre$ctype_no
expr_pre_m[is.infinite(expr_pre_m)]<-NA

expr_pre_combined <- cbind(expr0_pre,expr_pre_m)

expr_pre<-expr_pre_m
expr_pre[expr_pre == 0] <- NA

saveRDS(expr_pre,'backup_dist_expr_pre.rds')

#load previously saved distance matrix

expr_pre<- readRDS('backup_dist_expr_pre.rds')

#create expression + distance data frames

expr_pre_combined <- cbind(expr0_pre,expr_pre)

#improving robustness of distance relationships in the dataset by removing cell types without relationships or very rare cell types that would be overrepresented (sampling bias)

#remove cell type columns where there are no cells

expr_prerm0<-expr_pre[,colnames(expr_pre) %nin% colnames(expr_pre)[colSums(expr_pre, na.rm=T)==0]]

#remove rows/columns where there are very rare Cell types

expr_prerm<-expr_prerm0[rownames(expr_prerm0) %nin% exclude_pretx, colnames(expr_prerm0) %nin% exclude_pretx]


##Liver On-Tx ONLY====


expr_on<-c()

for(k in 1:length(OnTxid)){
  expr_k<-expr0_on[expr0_on$sample_id==OnTxid[k],] 
  
  #create placer cols for shortest distance to each cell type
  dummy<-matrix(nrow=nrow(expr_k),ncol=length(allcelltypes)) 
  colnames(dummy)<-legendctype$allcelltypes
  dummy<-as.data.frame(dummy)
  expr_k<-data.frame(expr_k,dummy)
  
  X <- ppp(x=expr_k$XMin,y=expr_k$YMin,marks = as.factor(expr_k$cluster), window = owin(xrange = range(expr_k$XMin), yrange = range(expr_k$YMin)),checkdup = F)
  distByCelltype <- nndist(X, by=marks(X))
  
  for (i in 1:(ncol(distByCelltype))){
    d <- which(colnames(expr_k) == (colnames(distByCelltype)[i]))
    expr_k[d] <- distByCelltype[,i]
  }
  expr_on <- rbind(expr_on, expr_k)
}
colnames(expr_on)<-c(colnames(expr0_on),legendctype$V1)
expr_on_m<-as.matrix(expr_on[,colnames(expr_on)[str_detect(colnames(expr_on),"ctype")]])
expr_on$ctype_no<-paste0("ctype",match(expr_on$cluster,legendctype$allcelltypes))
rownames(expr_on_m)<-expr_on$ctype_no
expr_on_m[is.infinite(expr_on_m)]<-NA

expr_on_combined <- cbind(expr0_on,expr_on_m)

expr_on<-expr_on_m
expr_on[expr_on == 0] <- NA

saveRDS(expr_on,'backup_dist_expr_on.rds')

#load previously saved distance matrix

expr_on<- readRDS('backup_dist_expr_on.rds')

#create expression + distance data frames
expr_on_combined <- cbind(expr0_on,expr_on)

#improving robustness of distance relationships in the dataset by removing cell types without relationships or very rare cell types that would be overrepresented (sampling bias)

#remove cell type columns where there are no cells

expr_onrm0<-expr_on[,colnames(expr_on) %nin% colnames(expr_on)[colSums(expr_on, na.rm=T)==0]]

#remove rows/columns where there are very rare Cell types

expr_onrm<-expr_onrm0[rownames(expr_onrm0) %nin% exclude_ontx, colnames(expr_onrm0) %nin% exclude_ontx]

##################################



####NETWORK VISUALIZATION####
###by a distance object (of matrix) of shortest distance between cell types

#for liver_pretx
mat_liver_pre=aggregate(x=expr_prerm, by=list(rownames(expr_prerm)), FUN=mean, na.rm=T)
groupnames<-mat_liver_pre$Group.1
mat_liver_pre<-as.matrix(mat_liver_pre[,2:ncol((mat_liver_pre))])
rownames(mat_liver_pre)<-str_sub(groupnames,6,) # add the rownames
colnames(mat_liver_pre)<-str_sub(colnames(mat_liver_pre),6,) #simplify colnames
mat_liver_preex<-mat_liver_pre[rownames(mat_liver_pre)[c(rownames(mat_liver_pre) %nin% "16")], colnames(mat_liver_pre)[c(colnames(mat_liver_pre) %nin% "16")]] #excluding Other subtypes
dist_liver_pre<-mat_liver_preex

#for liver_ontx
mat_liver_on=aggregate(x=expr_onrm, by=list(rownames(expr_onrm)), FUN=mean, na.rm=T)
groupnames<-mat_liver_on$Group.1
mat_liver_on<-as.matrix(mat_liver_on[,2:ncol((mat_liver_on))])
rownames(mat_liver_on)<-str_sub(groupnames,6,) # add the rownames
colnames(mat_liver_on)<-str_sub(colnames(mat_liver_on),6,) #simplify colnames
mat_liver_onex<-mat_liver_on[rownames(mat_liver_on)[c(rownames(mat_liver_on) %nin% "16")], colnames(mat_liver_on)[c(colnames(mat_liver_on) %nin% "16")]] #excluding Other subtypes
dist_liver_on<-mat_liver_onex

#color of legends revised
legendctypeex<-legendctype[legendctype$maintype!="Other",] #excluding Other subtypes
cols <- cbbPalette
colorlist <- cols[as.numeric(as.factor(legendctypeex$maintype))]


#generate visualizations of distance relationships using network qgraph

dev.off()




####Network Viz######
#for Liver in Pre-Tx
pdf("Distance_Relationships_Liver_PreTx.pdf",width=8,height=6)

xx<-dist_liver_pre
yy<-xx[as.character(sort(as.numeric(rownames(xx)))),#make sure all rows and cols are ordered
       as.character(sort(as.numeric(colnames(xx))))]

colorlist_liver_pre<-colorlist[as.numeric(colnames(mat_liver_pre))] #set color
names(colorlist_liver_pre)<-as.character(colnames(mat_liver_pre)) #set color reference names

g<-qgraph(1/yy, DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_pretx[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/yy,
          title="Pre-Tx Liver 1/yy (from each cluster node)",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_liver_pre,1/dist_liver_on, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=T, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlist_liver_pre[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=levels(as.factor(legendctype$maintype))[levels(as.factor(legendctype$maintype)) %nin% "Other"], 
       col = cols , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

g<-qgraph(1/t(yy), DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_pretx[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/t(yy), 
          title="Pre-Tx Liver 1/t(yy) (to each cluster node)",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_liver_pre,1/dist_liver_on, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=T, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlist_liver_pre[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=levels(as.factor(legendctype$maintype))[levels(as.factor(legendctype$maintype)) %nin% "Other"], 
       col = cols , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)


dev.off()


#for Liver in On-Tx
pdf("Distance_Relationships_Liver_OnTx.pdf",width=8,height=6)

xx<-dist_liver_on
yy<-xx[as.character(sort(as.numeric(rownames(xx)))),#make sure all rows and cols are ordered
       as.character(sort(as.numeric(colnames(xx))))]

colorlist_liver_on<-colorlist[as.numeric(colnames(mat_liver_on))] #set color
names(colorlist_liver_on)<-as.character(colnames(mat_liver_on)) #set color reference names

g<-qgraph(1/yy, DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_ontx[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/yy,
          title="On-Tx Liver 1/yy (from each cluster node)",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_liver_pre,1/dist_liver_on, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=T, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlist_liver_on[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=levels(as.factor(legendctype$maintype))[levels(as.factor(legendctype$maintype)) %nin% "Other"], 
       col = cols , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

g<-qgraph(1/t(yy), DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_ontx[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/t(yy), 
          title="On-Tx Liver 1/t(yy) (to each cluster node)",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_liver_pre,1/dist_liver_on, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=T, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlist_liver_on[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=levels(as.factor(legendctype$maintype))[levels(as.factor(legendctype$maintype)) %nin% "Other"], 
       col = cols , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)


dev.off()


###DIFF ANALYSIS OF DISTANCES: from Tc (ctype12)//SS ==============

dist_ttest<-function(group1=expr_pre,
                     group2=expr_on,
                     ctype12="ctype12",
                     ctype1="ctype1"){
  c1toc2_g1 = as.vector(group1[rownames(group1)==ctype12,ctype1])
  c1toc2_g2 = as.vector(group2[rownames(group2)==ctype12,ctype1])
  res<-t.test(x=c1toc2_g1,
                   y=c1toc2_g2,
                   paired=F)
  medians<-c(median(c1toc2_g1, na.rm=T),median(c1toc2_g2, na.rm=T))
  names(medians)<-c("PreTx_median","OnTx_median")
  estimate<-res$estimate
  names(estimate)<-c("PreTx_mean","OnTx_mean")
  diff_med<-medians["PreTx_median"]-medians["OnTx_median"]
  diff_mean<-estimate["PreTx_mean"]-estimate["OnTx_mean"]
  combine<-cbind(t(estimate), t(medians), diff_med=diff_med, diff_mean=diff_mean, pvalue=res$p.value)
  rownames(combine)<-paste(ctype12,ctype1,sep="_")
  print(combine)
}


res0<-data.frame(PreTx_mean=NA, OnTx_mean=NA, PreTx_median=NA, OnTx_median=NA, diff_med=NA, diff_mean=NA, pvalue=NA)
for(i in 1:length(legendctype$V1)){
  res_save<-dist_ttest(group1=expr_prerm0, group2=expr_onrm0, ctype12="ctype12", ctype1=as.character(legendctype$V1[i])) 
  res0<-rbind(res0,res_save)
} #omitted ctype16 (UA) but the loop goes from 1:16 - "ctype16" step will error out as subscript out of bounds
res_ctype1<-res0[2:nrow(res0),]  

res_ctype1$padj<-p.adjust(res_ctype1$pvalue, method = "BH")

write.csv(res_ctype1,"Results_distances_Tc.csv") 
write.csv(legendctype,"Legendsforctype.csv")

###VIOLIN PLOTS OF DISTANCES
#From Tc (ctype12)
df_Pre<-melt(expr_pre)
df_Pre$timegroup <- "PreTx"
df_On<-melt(expr_on)
df_On$timegroup <- "OnTx"
df_combined<-rbind(df_Pre,df_On)
colnames(df_combined)<-c("index","target","distance", "timegroup")
df_combined<-df_combined[!is.na(df_combined$distance),]

df_combined_c1ind<-df_combined[df_combined$index=="ctype12",]


pdf('Distance_fromTc_violinplots.pdf',width=3,height=3)
df_combined_c1ind$timegroup <- as.character(df_combined_c1ind$timegroup)
df_combined_c1ind$timegroup <- factor(df_combined_c1ind$timegroup, levels=c("PreTx", "OnTx"))
#Treg cells
ggplot(df_combined_c1ind[df_combined_c1ind$target=="ctype14",], aes(x=timegroup, y=distance, fill=timegroup))+
  ylim(0,1150)+ggtitle('Treg')+
  geom_violin()+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(color="black"))

#Th cells
ggplot(df_combined_c1ind[df_combined_c1ind$target=="ctype13",], aes(x=timegroup, y=distance, fill=timegroup))+
  ylim(0,950)+ggtitle('Th')+
  geom_violin()+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(color="black"))

#DC
ggplot(df_combined_c1ind[df_combined_c1ind$target=="ctype2",], aes(x=timegroup, y=distance, fill=timegroup))+
  ylim(0,850)+ggtitle('DC')+
  geom_violin()+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(color="black"))

#Gran
ggplot(df_combined_c1ind[df_combined_c1ind$target=="ctype3",], aes(x=timegroup, y=distance, fill=timegroup))+
  ylim(0,550)+ggtitle('Gran')+
  geom_violin()+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(color="black"))

#Tumor
ggplot(df_combined_c1ind[df_combined_c1ind$target=="ctype15",], aes(x=timegroup, y=distance, fill=timegroup))+
  ylim(0,550)+ggtitle('Tumor')+
  geom_violin()+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(color="black"))

#Mac_Merged
ggplot(df_combined_c1ind[df_combined_c1ind$target %in% c("ctype4", "ctype5", "ctype6", "ctype7"),], aes(x=timegroup, y=distance, fill=timegroup))+
  ylim(0,1200)+ggtitle('Mac_Merged')+
  geom_violin()+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(color="black"))

#Mac_I
ggplot(df_combined_c1ind[df_combined_c1ind$target=="ctype4",], aes(x=timegroup, y=distance, fill=timegroup))+
  ylim(0,900)+ggtitle('Mac_I')+
  geom_violin()+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(color="black"))

#Mac_II
ggplot(df_combined_c1ind[df_combined_c1ind$target=="ctype5",], aes(x=timegroup, y=distance, fill=timegroup))+
  ylim(0,500)+ggtitle('Mac_II')+
  geom_violin()+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(color="black"))

#Mac_III
ggplot(df_combined_c1ind[df_combined_c1ind$target=="ctype6",], aes(x=timegroup, y=distance, fill=timegroup))+
  ylim(0,900)+ggtitle('Mac_III')+
  geom_violin()+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(color="black"))

#Mac_IV
ggplot(df_combined_c1ind[df_combined_c1ind$target=="ctype7",], aes(x=timegroup, y=distance, fill=timegroup))+
  ylim(0,1200)+ggtitle('Mac_IV')+
  geom_violin()+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(color="black"))
dev.off()


#From Th
dist_ttest<-function(group1=expr_pre,
                     group2=expr_on,
                     ctype13="ctype13",
                     ctype1="ctype1"){
  c1toc2_g1 = as.vector(group1[rownames(group1)==ctype13,ctype1])
  c1toc2_g2 = as.vector(group2[rownames(group2)==ctype13,ctype1])
  res<-t.test(x=c1toc2_g1,
              y=c1toc2_g2,
              paired=F)
  medians<-c(median(c1toc2_g1, na.rm=T),median(c1toc2_g2, na.rm=T))
  names(medians)<-c("PreTx_median","OnTx_median")
  estimate<-res$estimate
  names(estimate)<-c("PreTx_mean","OnTx_mean")
  diff_med<-medians["PreTx_median"]-medians["OnTx_median"]
  diff_mean<-estimate["PreTx_mean"]-estimate["OnTx_mean"]
  combine<-cbind(t(estimate), t(medians), diff_med=diff_med, diff_mean=diff_mean, pvalue=res$p.value)
  rownames(combine)<-paste(ctype13,ctype1,sep="_")
  print(combine)
}


res0<-data.frame(PreTx_mean=NA, OnTx_mean=NA, PreTx_median=NA, OnTx_median=NA, diff_med=NA, diff_mean=NA, pvalue=NA)
for(i in 1:length(legendctype$V1)){
  res_save<-dist_ttest(group1=expr_prerm0, group2=expr_onrm0, ctype13="ctype13", ctype1=as.character(legendctype$V1[i])) 
  res0<-rbind(res0,res_save)
} #omitted ctype16 (UA) but the loop goes from 1:16 - "ctype16" step will error out as subscript out of bounds
res_ctype1<-res0[2:nrow(res0),]  

res_ctype1$padj<-p.adjust(res_ctype1$pvalue, method = "BH")

write.csv(res_ctype1,"Results_distances_Th.csv") 
write.csv(legendctype,"Legendsforctype.csv")

#From Th
df_Pre<-melt(expr_pre)
df_Pre$timegroup <- "PreTx"
df_On<-melt(expr_on)
df_On$timegroup <- "OnTx"
df_combined<-rbind(df_Pre,df_On)
colnames(df_combined)<-c("index","target","distance", "timegroup")
df_combined<-df_combined[!is.na(df_combined$distance),]

df_combined_c1ind<-df_combined[df_combined$index=="ctype13",]


pdf('Distance_fromTh_violinplots.pdf',width=3,height=3)
df_combined_c1ind$timegroup <- as.character(df_combined_c1ind$timegroup)
df_combined_c1ind$timegroup <- factor(df_combined_c1ind$timegroup, levels=c("PreTx", "OnTx"))
#Treg cells
ggplot(df_combined_c1ind[df_combined_c1ind$target=="ctype14",], aes(x=timegroup, y=distance, fill=timegroup))+
  ylim(0,850)+ggtitle('Treg')+
  geom_violin()+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(color="black"))

#Tc cells
ggplot(df_combined_c1ind[df_combined_c1ind$target=="ctype12",], aes(x=timegroup, y=distance, fill=timegroup))+
  ylim(0,350)+ggtitle('Tc')+
  geom_violin()+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(color="black"))

#DC
ggplot(df_combined_c1ind[df_combined_c1ind$target=="ctype2",], aes(x=timegroup, y=distance, fill=timegroup))+
  ylim(0,750)+ggtitle('DC')+
  geom_violin()+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(color="black"))

#Gran
ggplot(df_combined_c1ind[df_combined_c1ind$target=="ctype3",], aes(x=timegroup, y=distance, fill=timegroup))+
  ylim(0,350)+ggtitle('Gran')+
  geom_violin()+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(color="black"))

#Tumor
ggplot(df_combined_c1ind[df_combined_c1ind$target=="ctype15",], aes(x=timegroup, y=distance, fill=timegroup))+
  ylim(0,200)+ggtitle('Tumor')+
  geom_violin()+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(color="black"))

#Mac_Merged
ggplot(df_combined_c1ind[df_combined_c1ind$target %in% c("ctype4", "ctype5", "ctype6", "ctype7"),], aes(x=timegroup, y=distance, fill=timegroup))+
  ylim(0,1200)+ggtitle('Mac_Merged')+
  geom_violin()+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(color="black"))

#Mac_I
ggplot(df_combined_c1ind[df_combined_c1ind$target=="ctype4",], aes(x=timegroup, y=distance, fill=timegroup))+
  ylim(0,750)+ggtitle('Mac_I')+
  geom_violin()+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(color="black"))

#Mac_II
ggplot(df_combined_c1ind[df_combined_c1ind$target=="ctype5",], aes(x=timegroup, y=distance, fill=timegroup))+
  ylim(0,450)+ggtitle('Mac_II')+
  geom_violin()+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(color="black"))

#Mac_III
ggplot(df_combined_c1ind[df_combined_c1ind$target=="ctype6",], aes(x=timegroup, y=distance, fill=timegroup))+
  ylim(0,700)+ggtitle('Mac_III')+
  geom_violin()+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(color="black"))

#Mac_IV
ggplot(df_combined_c1ind[df_combined_c1ind$target=="ctype7",], aes(x=timegroup, y=distance, fill=timegroup))+
  ylim(0,1200)+ggtitle('Mac_IV')+
  geom_violin()+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(color="black"))

dev.off()

