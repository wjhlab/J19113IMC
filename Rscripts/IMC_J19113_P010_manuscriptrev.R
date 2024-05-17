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
  
  panel_fcs <- pData(flowCore::parameters(fcs_raw[[1]]))
  panel_fcs$desc <- gsub('-', '_', panel_fcs$desc)
  panel_fcs$desc[is.na(panel_fcs$desc)] <- paste0('NA_',which(is.na(panel_fcs$desc)))  
  
  rownames(panel_fcs) = panel_fcs$name
  
  ## Replace desc with revised Name
  panel_fcs[panel$Parameter,]$desc<-panel$Name
  
  ## Replace parameter data in flowSet with edits
  pData(flowCore::parameters(fcs_raw[[1]])) <- panel_fcs
  
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
  panellist <- rownames(pData(flowCore::parameters(fcs[[1]])))
  
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
  panel_fcs1 <- pData(flowCore::parameters(fcs1[[1]]))
  rownames(pData(flowCore::parameters(fcs1[[1]]))) <- rownames(panel_fcs[panel_fcs$desc %in% pData(flowCore::parameters(fcs1[[1]]))$name,])
  
  
  
  ###to scale every flowframe
  ## Save out the original rownames from the parameter list from the fcs flowFrame
  panellist <- rownames(pData(flowCore::parameters(fcs[[1]])))
  
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
output<-readRDS('backup_output_P010.rds')


## Read (skip if previously ran)

output <- returnfcs(metaDataFile = paste0(work,"/Config/J19113_P010_metadata.xlsx"),
                    panelDataFile = paste0(work,"/Config/J19113_P010_panel.xlsx"),
                    dataDirectory = paste0(work,"/Data"))


## Set up levels

samplevels=c("P010_t1_1","P010_t1_ablation2_1","P010_t1_ablation2_2",
             "P010_t2_1","P010_t2_2","P010_t2_3")

timegrouplevels=c("PreTx","OnTx")

timepointlevels=c("Baseline","Cycle1")

tumorlevels=c("ADRENAL")

caselevels=("P010")

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

output <- returnfcs(metaDataFile = paste0(work,"/Config/J19113_P010_metadata.xlsx"),
                    panelDataFile = paste0(work,"/Config/J19113_P010_panel.xlsx"),
                    dataDirectory = paste0(work,"/Data"))


##Clustering

output[(length(output)+1):(length(output)+3)] <- clusterfcs(fcs=output$fcs, numclusters=40, scaleoption = T, scaled.center = T, scaled.scale = T) 
#output$fcs uses just arcsin transformed data
#scaleoption scales the dataset across the channels so that the channels with the highest intensities will not dominate the clustering

names(output)[(length(output)-2):(length(output))] <- c('code_clustering','cell_clustering','metaclusters')


####ANNOTATIONS AND VISUALIZATION OF CLUSTERS####

## Load merge file
## Assign colors
clusterMergeFile = paste0(work,"/Config/J19113_P010_merge.xlsx") #create dummy merger numbers prior to annotation
cluster_merging <- read_excel(clusterMergeFile)

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
                                fileName = 'Clusteringheatmap_all_fcs_P010.pdf'); dev.off()

clusterMergeFile = paste0(work,"/Config/J19113_P010_merge.xlsx") #create dummy merger numbers prior to annotation
cluster_merging <- read_excel(clusterMergeFile)

clusterlevels=c("DC",
                "Gran",
                "Mac",
                "NK",
                "Stroma",
                "Tc_I",
                "Tc_II",
                "Th",
                "Tumor",
                "Tumor_Prolif",
                "UA")

colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_merging$new_cluster)))
clusternames<-clusterlevels
names(colorassigned)<-clusternames
mm1 <- match(output$cell_clustering, cluster_merging$original_cluster)
cell_clustering1m <- cluster_merging$new_cluster[mm1]
output$cell_clustering1m <- cell_clustering1m

plot_clustering_heatmap_wrapper2(fcs=output$fcs,
                                 colorbar = kovesi.diverging_bwr_40_95_c42(100),
                                 subtype_markers = c("CD45", "CD45RA", "CD45RO","DCLAMP", "FOXP3", "CD3", "CD4", "CD8", "GZMB", "CD20", "CD57", "CD15", "CD14", "CD16", "CD68", "CD163", "CD206", "SMA", "COL", "PDPN", "PDGFRA", "PanCK", "ECAD", "KI67"),
                                 color_clusters = colorassigned,
                                 cell_clustering = factor(output$cell_clustering1m, levels=clusterlevels), 
                                 fileName = 'Clusteringheatmap_merged_fcs_P010.pdf');dev.off()



## Save output list
saveRDS(output, file="backup_output_P010.rds")

#clusterheatmap for figure
clusterlevels2=c("DC","Gran", #2
                 "Mac", #3
                 "NK", "Th", "Tc_I", "Tc_II", #7
                 "Stroma", #8
                 "Tumor","Tumor_Prolif", #10
                 "UA")

wheretorowcut=c(2,3,7,8,10) #for cluster gaps

wheretocolcut=c(5,12,16,19) #for marker gaps

#heatmap scaled 10 to 90
plot_clustering_heatmap_wrapper4(fcs=output$fcs2,
                                 colorbar = rev(brewer.rdylbu(100)),
                                 subtype_markers = c("CD3", "CD4", "CD8", "GZMB", "CD57",#5
                                                     "DCLAMP","CD15", "CD14", "CD16", "CD68", "CD163", "CD206", #12
                                                     "SMA", "COL", "PDPN", "PDGFRA", #16
                                                     "PanCK", "ECAD", "KI67", #19
                                                     "CXCL12"),
                                 color_clusters = colorassigned,
                                 rowcut = wheretorowcut,
                                 colcut = wheretocolcut,
                                 cell_clustering = factor(output$cell_clustering1m, levels=clusterlevels2), 
                                 fileName = 'Clusteringheatmap_merged_forfigure.pdf');dev.off()


####ABUNDANCE PLOTS####

## Proportion calculations
counts_table <- table(output$cell_clustering1m, output$sample_ids)
counts <- as.data.frame.matrix(counts_table)


#Density calculations
areas <- tibble(sample_id=colnames(counts))
raw_areas <- read_xlsx(paste0(work,'/Config/J19113_P010_area.xlsx'))
areas$TotalArea <- raw_areas$TotalArea[match(areas$sample_id, raw_areas$sample_id)]
densities <- t(t(counts)/areas$TotalArea)

write.csv(counts,'Results_counts.csv')
write.csv(props,'Results_props.csv')
write.csv(densities, 'Results_densities.csv')

#Using counts, sum each cell type for PreTx and sum each cell type for OnTx
counts$PreTx <- rowSums(counts[,1:3])
counts$OnTx <- rowSums(counts[,4:6])
tg_counts <- counts[,7:8]

#calculate the sum of all PreTx area and all OnTx area in J19113_P010_timegrouparea.xlsx based on J19113_P010_area.xlsx
tg_areas <- read_xlsx(paste0(work,'/Config/J19113_P010_timegrouparea.xlsx'))
tg_densities <- t(t(tg_counts)/tg_areas$Area)

## Set up dataframe for timegroup densities
ggdfd <- melt(data.frame(cluster = rownames(tg_densities),tg_densities, check.names = FALSE),
              id.vars = "cluster", value.name = "tg_densities", 
              variable.name = "timegroup")

ggdfdclean <- ggdfd[ggdfd$cluster!="UNS",]


#plot box plots

#for densities
ggp2<-ggplot(ggdfd,aes(x=timegroup,y=tg_densities,fill=timegroup))+
  geom_bar(stat='identity', lwd=0.25, color="black")+
  facet_wrap(~cluster,ncol=4,scales="free")+
  ylab("# of Cells / Total Timegroup Area (mm^2)")+
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

pdf('Densities_box_TG_P010.pdf',width=5,height=5)
ggp2

dev.off()


