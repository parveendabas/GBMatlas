rm(list=ls()) # clear workspace

library(Seurat)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(grid)
library(gtable)
library(data.table)
library(cowplot)
library(plyr)
library(tidyr)
library(pheatmap)
library(tidyverse)
library(dplyr)
library(gridBase)
library(ggrepel)
library(ggplotify)
library(viridis)
library(SummarizedExperiment)

source("/Users/kumarpa/Desktop/Work/Jax/Scripts/Scripts_Kyuson/Functions_GBMatlas.R")

GroupName <- "Patient"
ClusterGroup <- "Cluster_08_09_10_11"
RDSname <- paste0("Seurat_Formatted_Normalized_USE_GBM_",ClusterGroup)
pkWD <- "/Users/kumarpa/Desktop/Work/Jax/Kyuson/Manuscript_Draft"
RDSdir <- paste0("/Users/kumarpa/Desktop/Work/Jax/Kyuson/Bill_Pipeline_Data/Data_Input/",ClusterGroup)
ColToUse <- "cluster.midres"
Suffix <- "MyeloidCells"
## Change Figure2 as well
NameInPdf.main <- paste0(ColToUse,"_Based_",Suffix)
IdentToSubsetColName="cluster.midres"
IdentToSubset="NA"
#downsampleHeatmap=25000


FDR=0.1
FoldChangeCutoff=1.4

  setwd(pkWD)
  plotWD <- paste(getwd(),paste0("For_Github_Main_Figures"),sep="/"); print(plotWD)
  dir.create(file.path(getwd(),paste0("For_Github_Main_Figures")), showWarnings = FALSE)

  setwd(plotWD)
  plotWD1 <- paste(getwd(),paste0("Figure2_",Suffix),sep="/"); print(plotWD1)
  dir.create(file.path(getwd(),paste0("Figure2_",Suffix)), showWarnings = FALSE)
  
  ClusOrder.main <- c("Cluster01", "Cluster02", "Cluster03", "Cluster04", "Cluster05", "Cluster06", "Cluster07", "Cluster08", "Cluster09", "Cluster10", "Cluster11", "Cluster12", "Cluster13", "Cluster14"); ClusOrder.main
  cbPalette.Cluster <- c("darkblue", "#0073ff" ,"#009FFF" ,"#00EEFF", "#FF0073", "#FF00D3" ,"#C900FF" ,"#FF8C00" ,"#FF5900", "#FF2500" ,"darkred" ,"#5FFF00", "#AFFF00" ,"#FFFF00")
  MClusOrder <- c("MC1", "MC2", "MC3", "MC4", "MC5", "MC6", "MC7"); MClusOrder
  MClusOrderPalette <- c("darkred", "#008000", "blue", "#00ff00", "#ff8c00", "red",  "#ffff00")
  CTOrder <- c("s-macr1", "a-microglia", "APC", "r-microglia", "s-macr2", "dividing-mac", "M1-mac"); CTOrder
  cbPalette.CT <- c("darkred", "#008000", "blue", "#00ff00", "#ff8c00", "red",  "#ffff00")
  GroupOrder.temp = c("CNSTM-068", "CNSTM-070", "CNSTM-081", "CNSTM-096")
  GroupPalette.temp <- c("#0039d1","#c20f00","purple","#26cc00")
  FragOrder = c("CNSTM-068-A", "CNSTM-068-B", "CNSTM-068-C", "CNSTM-068-D", "CNSTM-070-A", "CNSTM-070-C", "CNSTM-070-D", "CNSTM-070-F", "CNSTM-081-A", "CNSTM-081-B", "CNSTM-081-C", "CNSTM-081-D", "CNSTM-096-1", "CNSTM-096-2", "CNSTM-096-4", "CNSTM-096-5")
  FragPalette <- c("#0000FF", "#0045FF", "#008AFF" ,"#00CFFF", "#FF9A00", "#FF6E00", "#FF1400" ,"darkred" ,"#FF0052", "#FF00A6", "#FF00F9" ,"#B000FF","darkgreen","#30FF00","#B9FF00", "#FFFF00")
  PredictionOrder <- c("diploid", "aneuploid")
  PredictionPalette <- c("#ff7f50", "#009FFF")
  
  
  setwd(RDSdir)
  FullInfo <- read.table(file = paste0("Full_Information_",ClusterGroup,"_Cells.txt"), header = T, sep = "\t"); head(FullInfo); dim(FullInfo)
  
  setwd(RDSdir)
  SCdata.main <- readRDS(paste0(RDSname,".rds"))
  SCdata.main
  head(SCdata.main@meta.data)
  dim(GetAssayData(object = SCdata.main, slot = "scale.data")) 
  colnames(SCdata.main@meta.data)[colnames(SCdata.main@meta.data) == 'cluster.mid-res'] <- 'cluster.midres'
  SCdata.main$Custom.Cluster <- gsub(" ", "",SCdata.main$Custom.Cluster)
  table(SCdata.main@meta.data[,ColToUse])
  table(as.character(as.numeric(SCdata.main@meta.data[,ColToUse])+1))
  SCdata.main@meta.data[,ColToUse] <- as.character(as.numeric(SCdata.main@meta.data[,ColToUse]) + 1)
  SCdata.main@meta.data[,ColToUse] <- paste0("MC",SCdata.main@meta.data[,ColToUse])
  table(SCdata.main@meta.data[,ColToUse])
  SCdata.main@meta.data$Patient <- FullInfo$Patient[match(rownames(SCdata.main@meta.data), rownames(FullInfo))]
  SCdata.main@meta.data$Fragment <- FullInfo$Fragment[match(rownames(SCdata.main@meta.data), rownames(FullInfo))]
  SCdata.main@meta.data$Custom.Cluster <- FullInfo$Custom.Cluster[match(rownames(SCdata.main@meta.data), rownames(FullInfo))]
  SCdata.main@meta.data$prediction <- FullInfo$prediction[match(rownames(SCdata.main@meta.data), rownames(FullInfo))]
  SCdata.main$prediction <- gsub("tumor", "aneuploid",SCdata.main$prediction)
  SCdata.main$prediction <- factor(SCdata.main$prediction, levels = PredictionOrder)
  Idents(SCdata.main) <- ColToUse
  DefaultAssay(SCdata.main) <- "RNA"
  
  
  
  
  ########################################################################################################################
  ################################*********************************#######################################################
  RUNFigure2="YES"
  if(RUNFigure2=="YES"){
  ### Figure2 Clsuter and heatmaps  #"#787878", "#9A6324", "#ee1289", 
  ToUsePallete <- MClusOrderPalette
  ToUseOrder <- MClusOrder
  ToUseCol <- ColToUse
  
  MarkerDir <- paste0("/Users/kumarpa/Desktop/Work/Jax/Kyuson/Bill_Pipeline_Data/DGEs/DGEs_",ToUseCol,"_Based_",Suffix)
  setwd(MarkerDir)
  MarkerGenes <- read.table(file = paste0("DEGs_Heatmap_",ToUseCol,"_Based_",Suffix,".txt"), header = T, sep = "\t"); head(MarkerGenes); dim(MarkerGenes)
  
  
  ######## **************************************************  Panel a: HEATMAP  ****************************************************
  ### Main Figures
  RUNFigure2a="YES"
  if(RUNFigure2a=="YES"){
  markers <- MarkerGenes
  markers <- markers[markers$p_val_adj < FDR,]; dim(markers)
  markers <- markers[!rownames(markers) %like% "^RP[SL]",]; dim(markers)
  markers <- markers[!rownames(markers) %like% "^RP11",]; dim(markers)
  markers <- markers[!rownames(markers) %like% "^MT",]; dim(markers)
  markers <- markers[!rownames(markers) %like% "^RNA",]; dim(markers)
  topnumber=30
  top <- markers %>% group_by(cluster) %>% top_n(n = topnumber, wt = avg_logFC); dim(top)
  topFDR <- markers %>% dplyr::arrange(p_val_adj, desc(avg_logFC)) %>% dplyr::group_by(cluster) %>% dplyr::slice(1:topnumber); dim(topFDR)
  setwd(MarkerDir)
  temp <- as.data.frame(topFDR)
  temp.unique <- temp[!duplicated(temp$gene),]
  write.table(temp, file = paste0("DEGs_Heatmap_",ToUseCol,"_Based_",Suffix,"_Top",topnumber,"_Genes.txt"),quote=F,sep="\t")
  write.table(temp.unique, file = paste0("DEGs_Heatmap_",ToUseCol,"_Based_",Suffix,"_Top",topnumber,"_Genes_UNIQUE.txt"),quote=F,sep="\t")
  SCdata.main$Temp <- "Temp"
  Idents(SCdata.main) <- "Temp"
  #SCdata.temp.Heatmap <- subset(SCdata.main, downsample=downsampleHeatmap)
  SCdata.temp.Heatmap <- SCdata.main
  Idents(SCdata.main) <- ToUseCol
  Idents(SCdata.main) <- factor(Idents(SCdata.main), levels= ToUseOrder)
  SCdata.temp.Heatmap <- ScaleData(object = SCdata.temp.Heatmap, verbose = FALSE, features = markers$gene, scale.max = 2)
  dtype="scale.data"
  nor.exp <- GetAssayData(object = SCdata.temp.Heatmap, slot = dtype); print(dim(nor.exp))
  UseGenes <- intersect(markers$gene, rownames(nor.exp)); length(UseGenes)
  nor.exp <- nor.exp[UseGenes,,drop=FALSE]; dim(nor.exp)
  table(SCdata.temp.Heatmap@meta.data[,ToUseCol], SCdata.temp.Heatmap@meta.data[,GroupName])
  #meta.data.plot <- SCdata.temp.Heatmap@meta.data[,c(GroupName, "Fragment", ColToUse)]
  meta.data.plot <- SCdata.temp.Heatmap@meta.data[,c(GroupName, ColToUse)]
  meta.data.plot[,ColToUse] <- factor(meta.data.plot[,ColToUse], levels = ToUseOrder)
  meta.data.plot[,GroupName] <- factor(meta.data.plot[,GroupName], levels = GroupOrder.temp)
  meta.data.plot <- meta.data.plot[order(factor(meta.data.plot[,ColToUse], levels = ToUseOrder), factor(meta.data.plot[,GroupName], levels = GroupOrder.temp)),]; head(meta.data.plot)
  meta.data.plot[,ColToUse] <- as.character(meta.data.plot[,ColToUse])
  nor.exp <- nor.exp[,rownames(meta.data.plot), drop=FALSE]
  colnames(meta.data.plot) <- c("Patient", "Cluster")
  head(meta.data.plot)
  print(dim(nor.exp))
  
  Cluster.FULL =  MClusOrderPalette[1:length(MClusOrder)]; names(Cluster.FULL) <- MClusOrder; Cluster.FULL
  Fragment.FULL = FragPalette[1:length(FragOrder)]; names(Fragment.FULL) <- FragOrder; Fragment.FULL
  prediction.FULL = PredictionPalette[1:length(PredictionOrder)]; names(prediction.FULL) <- PredictionOrder; prediction.FULL
  Patient.FULL = GroupPalette.temp[1:length(GroupOrder.temp)]; names(Patient.FULL) <- GroupOrder.temp; Patient.FULL
  
  ann_colors = list(
    Cluster = Cluster.FULL[names(Cluster.FULL) %in% as.character(unique(meta.data.plot$Cluster))],
    Fragment = Fragment.FULL[names(Fragment.FULL) %in% as.character(unique(meta.data.plot$Fragment))],
    prediction = prediction.FULL[names(prediction.FULL) %in% as.character(unique(meta.data.plot$prediction))],
    Patient = Patient.FULL[names(Patient.FULL) %in% as.character(unique(meta.data.plot$Patient))]
  )
  
  ### Lighter Shade of blue near Zero
  colors <- c(seq(-2,2,by=0.01))
  my_palette <- c(colorRampPalette(colors = c("darkblue", "#a7c5f2", "#e6f0f5", "gray97", "darksalmon", "orangered3", "darkred")) (n = length(colors)))
  
  Temp <- pheatmap(nor.exp[as.character(unique(topFDR$gene)),],annotation_col=meta.data.plot,show_colnames=F,show_rownames=T,  color=my_palette, breaks=colors, 
                   annotation_colors = ann_colors, cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 10); dev.off()
  #as.grob(pheatmap(nor.exp[as.character(top$gene),],annotation_col=meta.data.plot,show_colnames=F,show_rownames=T,  color=my_palette, breaks=colors, 
  #                       annotation_colors = ann_colors, cluster_rows = FALSE, cluster_cols = FALSE, main = paste0("Heatmap All markers"), fontsize = 13)); dev.off()
  
  sigGenes = c("FN1", "S100A4", "IBSP", "CCL3", "CCL4", "CD83", "IL1B", "B2M", "CD74", "HLA-DRA", "CD14", "APOE", "TMEM119", "BHLHE41", "SORL1", "CX3CR1", "MIF", "MKI67", "IER3", "NFKBIZ")
  ##heatmap <- Temp$gtable
  ##new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  ##new.label$label
  
  p4 <- Gene.Labels.pheatmap(Temp, kept.labels = sigGenes, repel.degree = 0)
  
  
  
  setwd(plotWD1)
  pdf(paste0("Figure2_",ToUseCol,"_",Suffix,"_Panel_a_HEATMAP_UPLOAD.pdf"), height = 5, width = 8)
  print(plot_grid(p4))
  dev.off()
  }
  
  
  ######## **************************************************  Panel b: UMAP  ****************************************************
  ### Main Figures
  RUNFigure2b="YES"
  if(RUNFigure2b=="YES"){
    Idents(SCdata.main) <- ToUseCol
    Idents(SCdata.main) <- factor(Idents(SCdata.main), levels= ToUseOrder)
    r2.1 <- DimPlot(SCdata.main, reduction = "umap", cols = ToUsePallete, label = F, label.size = 1.5, pt.size = 0.01) + 
      theme(axis.title.x=element_blank(), axis.title.y = element_blank(), legend.key.size = unit(0.4,"line"), legend.text=element_text(size=13)) +
      guides(color = guide_legend(override.aes = list(size = 1.2)))
    
    setwd(plotWD1)
    pdf(paste0("Figure2_",ToUseCol,"_UMAP_",Suffix,"_Panel_b_Cluster_UPLOAD.pdf"), height = 3, width = 5)
    print(plot_grid(r2.1))
    dev.off()
  }
  
  
  
  
  }
  