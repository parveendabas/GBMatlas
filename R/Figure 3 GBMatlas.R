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

GroupName <- "Patient"
ClusterGroup <- "Cluster_07"
RDSname <- paste0("Seurat_Formatted_Normalized_USE_GBM_",ClusterGroup)
pkWD <- "/Users/kumarpa/Desktop/Work/Jax/Kyuson/Manuscript_Draft"
RDSdir <- paste0("/Users/kumarpa/Desktop/Work/Jax/Kyuson/Bill_Pipeline_Data/Data_Input/",ClusterGroup)
ColToUse <- "cluster.R2"
Suffix <- "TCells"
## Change Figure4 as well
NameInPdf.main <- paste0(ColToUse,"_Based_",Suffix)
IdentToSubsetColName="cluster.R2"
IdentToSubset="NA"
##downsampleHeatmap=5000

FDR=0.1
FoldChangeCutoff=1.4

source("/Users/kumarpa/Desktop/Work/Jax/Scripts/Scripts_Kyuson/Functions_Kyuson_Project_GitHub.R")

setwd(pkWD)
plotWD <- paste(getwd(),paste0("For_Github_Main_Figures"),sep="/"); print(plotWD)
dir.create(file.path(getwd(),paste0("For_Github_Main_Figures")), showWarnings = FALSE)


  setwd(plotWD)
  plotWD1 <- paste(getwd(),paste0("Figure3_",Suffix),sep="/"); print(plotWD1)
  dir.create(file.path(getwd(),paste0("Figure3_",Suffix)), showWarnings = FALSE)
  
  #Custom.Cluster = c(`Cluster01` = "#0000ee", `Cluster02` = "#27408B", `Cluster03` = "#56B4E9",  `Cluster14` = "#00ffff", `Cluster04` = "#00ff00", `Cluster05` = "#DEB887", `Cluster06` = "#008000", `Cluster07` = "#bf3eff", `Cluster08` = "#ff0000", `Cluster09` = "#8b0000", `Cluster10` = "#CD5C5C", `Cluster11` = "#ff80bf", `Cluster12` = "#F0E442", `Cluster13` = "#ffff00"),
  ClusOrder.main <- c("Cluster01", "Cluster02", "Cluster03", "Cluster04", "Cluster05", "Cluster06", "Cluster07", "Cluster08", "Cluster09", "Cluster10", "Cluster11", "Cluster12", "Cluster13", "Cluster14"); ClusOrder.main
  cbPalette.Cluster <- c("darkblue", "#0073ff" ,"#009FFF" ,"#00EEFF", "#FF0073", "#FF00D3" ,"#C900FF" ,"#FF8C00" ,"#FF5900", "#FF2500" ,"darkred" ,"#5FFF00", "#AFFF00" ,"#FFFF00")
  TClusOrder <- c("TC1", "TC2", "TC3", "TC4", "TC5", "TC6", "TC7", "TC8", "TC9"); TClusOrder
  TClusOrderPalette <- c("#008000", "darkred", "#ffff00", "#00ff00", "#C900FF" ,"#00EEFF", "#0073ff", "#FF2500", "#FF00D3")
  CTOrder <- c("CD4 T", "NK1", "Dividing T", "Tregs", "Resting T", "CD8 Trm", "CD8 T", "NK2", "CD8 Tmem"); CTOrder
  cbPalette.CT <- c("#008000", "darkred", "#ffff00", "#00ff00", "#C900FF" ,"#00EEFF", "#0073ff", "#FF2500", "#FF00D3")
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
  colnames(SCdata.main@meta.data)[colnames(SCdata.main@meta.data) == 'cluster.R2'] <- 'cluster.R2'
  SCdata.main$Custom.Cluster <- gsub(" ", "",SCdata.main$Custom.Cluster)
  table(SCdata.main@meta.data[,ColToUse])
  table(as.character(as.numeric(SCdata.main@meta.data[,ColToUse])))
  SCdata.main@meta.data[,ColToUse] <- as.character(as.numeric(SCdata.main@meta.data[,ColToUse]))
  SCdata.main@meta.data[,ColToUse] <- paste0("TC",SCdata.main@meta.data[,ColToUse])
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
  RUNFigure3="YES"
  if(RUNFigure3=="YES"){
  ### Figure3 Clsuter and heatmaps  #"#787878", "#9A6324", "#ee1289", 
  ToUsePallete <- TClusOrderPalette
  ToUseOrder <- TClusOrder
  ToUseCol <- ColToUse
  
  MarkerDir <- paste0("/Users/kumarpa/Desktop/Work/Jax/Kyuson/Bill_Pipeline_Data/DGEs/DGEs_",ToUseCol,"_Based_",Suffix)
  setwd(MarkerDir)
  MarkerGenes <- read.table(file = paste0("DEGs_Heatmap_",ToUseCol,"_Based_",Suffix,".txt"), header = T, sep = "\t"); head(MarkerGenes); dim(MarkerGenes)
  
  
  ######## **************************************************  Panel a: HEATMAP  ****************************************************
  ### Main Figures
  RUNFigure3a="YES"
  if(RUNFigure3a=="YES"){
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
    write.table(temp, file = paste0("DEGs_Heatmap_",ToUseCol,"_Based_",Suffix,"_Top",topnumber,"_Genes.txt"),quote=F,sep="\t")
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
    check <- intersect(GroupOrder.temp, unique(SCdata.temp.Heatmap@meta.data[,GroupName])); check
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
    
    Cluster.FULL =  TClusOrderPalette[1:length(TClusOrder)]; names(Cluster.FULL) <- TClusOrder; Cluster.FULL
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
    
    sigGenes = c("IL7R", "MKI67", "FOXP3", "IL2RA", "CTLA4", "TIGIT", "S100A4", "GZMK", "HLA-DRB1", "CD2", "IFI6", "GNLY", "NKG7", "CD3E", "CD8A", "HOPX")
    ##heatmap <- Temp$gtable
    ##new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
    ##new.label$label
    
    p4 <- Gene.Labels.pheatmap(Temp, kept.labels = sigGenes, repel.degree = 0)
    
    setwd(plotWD1)
    pdf(paste0("Figure3_",ToUseCol,"_",Suffix,"_Panel_a_HEATMAP_UPLOAD.pdf"), height = 4, width = 7.5)
    print(plot_grid(p4))
    dev.off()
  }
  
  
  
  ######## **************************************************  Panel b: UMAP  ****************************************************
  ### Main Figures
  RUNFigure3b="YES"
  if(RUNFigure3b=="YES"){
  Idents(SCdata.main) <- ToUseCol
  Idents(SCdata.main) <- factor(Idents(SCdata.main), levels= ToUseOrder)
  
  
  r2.1 <- DimPlot(SCdata.main, reduction = "umap", cols = ToUsePallete, label = F, label.size = 3, pt.size = 1.8) + 
    theme(axis.title.x=element_blank(), axis.title.y = element_blank(), legend.key.size = unit(0.45,"line"), legend.text=element_text(size=13)) +
    guides(color = guide_legend(override.aes = list(size = 1.2)))
  
  setwd(plotWD1)
  pdf(paste0("Figure3_",ToUseCol,"_",Suffix,"_Panel_b_Cluster_UMAP_UPLOAD.pdf"), height = 3, width = 5.5)
  print(plot_grid(r2.1))
  dev.off()
  }
  
  
  
  
  
  
  RUNFigure3c="YES"
  if(RUNFigure3c=="YES"){
    setwd(MarkerDir)
    Idents(SCdata.main) <- ToUseCol
    Idents(SCdata.main) <- factor(Idents(SCdata.main), levels= ToUseOrder)
    
    TregGenes <- c("FOXP3", "IL2RA", "TIGIT", "CTLA4", "S100A4")
    PlotGenes <- data.frame(Genes=TregGenes, CT=rep("CT",5))
    GeneCol="Genes"; GroupCol="CT"
    row.object <- Voilin_Plot_OneGenePerLine(SCdata.main, PlotGenes, GeneCol, GroupCol, cbPalette.CT, 1) 
    #p1 <- plot_grid(plotlist=row.object, nrow = length(row.object), labels=unique(PlotGenes[,GroupCol]), label_colour="darkred", hjust = -0.15)
    setwd(plotWD1)
    pdf(paste0("Figure3_",ToUseCol,"_",Suffix,"_Panel_c_Violin_UPLOAD.pdf"), height = 7, width = 5)
    print(plot_grid(plotlist=row.object, nrow = length(row.object)))
    dev.off()
    
  }
  
  
  ######## **************************************************  Panel c: Violin  ****************************************************
  ### Main Figures
  RUNFigure3d="YES"
  if(RUNFigure3d=="YES"){
    setwd(plotWD1)
    Idents(SCdata.main) <- ToUseCol
    Idents(SCdata.main) <- factor(Idents(SCdata.main), levels= ToUseOrder)
    pdf(paste0("Figure3_",ToUseCol,"_",Suffix,"_Panel_d_Violin_Tcell_UPLOAD.pdf"), height = 3, width = 4.3)
    print(VlnPlot(SCdata.main, features = c("S100A4"), cols = ToUsePallete, pt.size = 0.0)  + 
            theme(axis.title.x=element_blank(), axis.title.y = element_blank(), axis.text=element_text(size=14),
                  legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = 16), legend.position = 'none') +
            stat_summary(fun=median, geom="point", size=1, color="black"))
    dev.off()
  }
  
  
  
  
  }
  ################################*********************************#######################################################
  ########################################################################################################################
  
  