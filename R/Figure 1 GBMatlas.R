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
ClusterGroup <- "PatientALL"
RDSname <- paste0("GBM_",ClusterGroup)
pkWD <- "/Users/kumarpa/Desktop/GBMatlas_Code"
RDSdir <- paste0("/Users/kumarpa/Desktop/GBMatlas_Code/data/")
ColToUse <- "Custom.Cluster"
Suffix <- "PatientALL"
NameInPdf.main <- paste0(ColToUse,"_Based_",Suffix)
IdentToSubsetColName="Custom.Cluster"
downsampleHeatmap=25000
FDR=0.1

source("/Users/kumarpa/Desktop/GBMatlas_Code/Functions_GBMatlas.R")

setwd(pkWD)
plotWD <- paste(getwd(),paste0("GBMatlas_Main_Figures"),sep="/"); print(plotWD)
dir.create(file.path(getwd(),paste0("GBMatlas_Main_Figures")), showWarnings = FALSE)

setwd(plotWD)
plotWD1 <- paste(getwd(),paste0("Figure1_",Suffix),sep="/"); print(plotWD1)
dir.create(file.path(getwd(),paste0("Figure1_",Suffix)), showWarnings = FALSE)

## Assign Color Palettes
ClusOrder.main <- c("Cluster01", "Cluster02", "Cluster03", "Cluster04", "Cluster05", "Cluster06", "Cluster07", "Cluster08", "Cluster09", "Cluster10", "Cluster11", "Cluster12", "Cluster13", "Cluster14"); ClusOrder.main
ClusOrder <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"); ClusOrder
ClusterPalette <- c("darkblue", "#0073ff" ,"#009FFF" ,"#00EEFF", "#FF0073", "#FF00D3" ,"#C900FF" ,"#FF8C00" ,"#FF5900", "#FF2500" ,"darkred" ,"#5FFF00", "#AFFF00" ,"#FFFF00")
CTOrder <- c("Glioma", "Pericyte", "Tcells", "Myeloid", "Oligo", "Mix"); CTOrder
CTPalette <- c("darkblue", "#FF00D3", "#C900FF", "darkred", "#5FFF00", "grey")
GroupOrder = c("CNSTM-068", "CNSTM-070", "CNSTM-081", "CNSTM-096")
GroupPalette <- c("#0039d1","#c20f00","purple","#26cc00")
FragOrder = c("CNSTM-068-A", "CNSTM-068-B", "CNSTM-068-C", "CNSTM-068-D", "CNSTM-070-A", "CNSTM-070-C", "CNSTM-070-D", "CNSTM-070-F", "CNSTM-081-A", "CNSTM-081-B", "CNSTM-081-C", "CNSTM-081-D", "CNSTM-096-1", "CNSTM-096-2", "CNSTM-096-4", "CNSTM-096-5")
FragPalette <- c("#0000FF", "#0045FF", "#008AFF" ,"#00CFFF", "#FF9A00", "#FF6E00", "#FF1400" ,"darkred" ,"#FF0052", "#FF00A6", "#FF00F9" ,"#B000FF","darkgreen","#30FF00","#B9FF00", "#FFFF00")
PredictionOrder <- c("diploid", "aneuploid")
PredictionPalette <- c("#ff7f50", "#009FFF")

## Load Data
setwd(RDSdir)
SCdata.main <- readRDS(paste0(RDSname,".rds"))
#saveRDS(SCdata.main, paste0(RDSname,".rds"))
DefaultAssay(SCdata.main) <- "RNA"


########################################################################################################################
################################*********************************#######################################################
RUNFigure1="YES"
if(RUNFigure1=="YES"){
  ### Figure1 Clsuter and heatmaps  #"#787878", "#9A6324", "#ee1289", 
  ToUsePallete <- ClusterPalette
  ToUseOrder <- ClusOrder.main
  ToUseCol <- ColToUse
  
  setwd(RDSdir)
  MarkerGenes <- read.table(file = paste0("DEGs_",Suffix,".txt"), header = T, sep = "\t"); head(MarkerGenes); dim(MarkerGenes)
  
  ######## **************************************************  Panel b: Heatmap  ****************************************************
  ### Main Figures
  # Panel b: Heatmap
  RUNFigure1b="YES"
  if(RUNFigure1b=="YES"){
  markers <- MarkerGenes[MarkerGenes$p_val_adj < FDR,]; dim(markers)
  
  topnumber=30
  topFDR <- markers %>% dplyr::arrange(p_val_adj, desc(avg_logFC)) %>% dplyr::group_by(cluster) %>% dplyr::slice(1:topnumber); dim(topFDR)
  
  SCdata.main$Temp <- "Temp"
  Idents(SCdata.main) <- "Temp"
  SCdata.temp.Heatmap <- subset(SCdata.main, downsample=downsampleHeatmap)
  Idents(SCdata.main) <- ToUseCol
  Idents(SCdata.main) <- factor(Idents(SCdata.main), levels= ToUseOrder)
  SCdata.temp.Heatmap <- ScaleData(object = SCdata.temp.Heatmap, verbose = FALSE, features = markers$gene, scale.max = 2)
  nor.exp <- GetAssayData(object = SCdata.temp.Heatmap, slot = "scale.data"); print(dim(nor.exp))
  UseGenes <- intersect(markers$gene, rownames(nor.exp)); length(UseGenes)
  nor.exp <- nor.exp[UseGenes,,drop=FALSE]; dim(nor.exp)
  meta.data.plot <- SCdata.temp.Heatmap@meta.data[,c(GroupName, ColToUse)]
  meta.data.plot[,ColToUse] <- factor(meta.data.plot[,ColToUse], levels = ToUseOrder)
  meta.data.plot[,GroupName] <- factor(meta.data.plot[,GroupName], levels = GroupOrder)
  meta.data.plot <- meta.data.plot[order(factor(meta.data.plot[,ColToUse], levels = ToUseOrder), factor(meta.data.plot[,GroupName], levels = GroupOrder)),]; head(meta.data.plot)
  meta.data.plot[,ColToUse] <- as.character(meta.data.plot[,ColToUse])
  nor.exp <- nor.exp[,rownames(meta.data.plot), drop=FALSE]
  colnames(meta.data.plot) <- c("Patient", "Cluster")
  
  Custom.Cluster.FULL =  ClusterPalette[1:length(ClusOrder.main)]; names(Custom.Cluster.FULL) <- ClusOrder.main; Custom.Cluster.FULL
  Fragment.FULL = FragPalette[1:length(FragPalette)]; names(Fragment.FULL) <- FragPalette; Fragment.FULL
  prediction.FULL = PredictionPalette[1:length(PredictionOrder)]; names(prediction.FULL) <- PredictionOrder; prediction.FULL
  Patient.FULL = GroupPalette[1:length(GroupOrder)]; names(Patient.FULL) <- GroupOrder; Patient.FULL
  
  ann_colors = list(
    Cluster = Custom.Cluster.FULL[names(Custom.Cluster.FULL) %in% as.character(unique(meta.data.plot$Cluster))],
    Fragment = Fragment.FULL[names(Fragment.FULL) %in% as.character(unique(meta.data.plot$Fragment))],
    prediction = prediction.FULL[names(prediction.FULL) %in% as.character(unique(meta.data.plot$prediction))],
    Patient = Patient.FULL[names(Patient.FULL) %in% as.character(unique(meta.data.plot$Patient))]
  )
  
  
  ### Lighter Shade of blue near Zero
  colors <- c(seq(-2,2,by=0.01))
  my_palette <- c(colorRampPalette(colors = c("darkblue", "#a7c5f2", "#e6f0f5", "gray97", "darksalmon", "orangered3", "darkred")) (n = length(colors)))
  
  Temp <- pheatmap(nor.exp[as.character(unique(topFDR$gene)),],annotation_col=meta.data.plot,show_colnames=F,show_rownames=T,  color=my_palette, breaks=colors, 
                   annotation_colors = ann_colors, cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 10, fontsize_row=9); dev.off()
  
  sigGenes = c("PDGFRA", "OLIG1", "GFAP", "SOX2", "ACTA2", "CD3G", "S100A9", "S100A4", "ITGAM", "CD68", "APOE", "HLA-DRA", "P2RY12", "CCL3", "MBP")
  
  p4 <- Gene.Labels.pheatmap(Temp, kept.labels = sigGenes, repel.degree = 0)
  
  setwd(plotWD1)
  pdf(paste0("Figure1","_",Suffix,"_Panel_b.pdf"), height = 8, width = 6)
  print(plot_grid(p4))
  dev.off()
  rm(SCdata.temp.Heatmap)
  }
  
  
  ######## **************************************************  Panel c: UMAP  ****************************************************
  ### Main Figures
  
  RUNFigure1c="YES"
  if(RUNFigure1c=="YES"){
  Idents(SCdata.main) <- ToUseCol
  Idents(SCdata.main) <- factor(Idents(SCdata.main), levels= ToUseOrder)
  
  setwd(plotWD1)
  pdf(paste0("Figure1_",Suffix,"_Panel_c.1.pdf"), height = 2, width = 3)
  print(DimPlot(SCdata.main, reduction = "umap", cols = ToUsePallete, label = F, label.size = 1.5, pt.size = 0.01) + 
    theme(axis.title.x=element_blank(), axis.title.y = element_blank(), legend.key.size = unit(0.4,"line"), legend.text=element_text(size=9)) +
    guides(color = guide_legend(override.aes = list(size = 1.2))))
  dev.off()
  
  CellToPlot <- rownames(SCdata.main@meta.data[SCdata.main@meta.data$prediction %in% c("aneuploid", "diploid"),])
  setwd(plotWD1)
  pdf(paste0("Figure1_",Suffix,"_Panel_c.2.pdf"), height = 2, width = 3)
  print(DimPlot(SCdata.main, reduction = "umap", cells = CellToPlot, cols = PredictionPalette, label = F, label.size = 4, pt.size = 0.01, group.by = "prediction") + 
    theme(axis.title.x=element_blank(), axis.title.y = element_blank(), legend.key.size = unit(0.5,"line"), legend.text=element_text(size=11)) +
    guides(color = guide_legend(override.aes = list(size = 2))))
  dev.off()
  
  
  Idents(SCdata.main) <- "CT"
  Idents(SCdata.main) <- factor(Idents(SCdata.main), levels= CTOrder)
  setwd(plotWD1)
  pdf(paste0("Figure1_",Suffix,"_Panel_c.3.pdf"), height = 2, width = 3)
  print(DimPlot(SCdata.main, reduction = "umap", cols = CTPalette[CTOrder %in% unique((SCdata.main$CT))], label = F, label.size = 6, pt.size = 0.01) + 
    theme(axis.title.x=element_blank(), axis.title.y = element_blank(), legend.key.size = unit(0.5,"line"), legend.text=element_text(size=11)) +
    guides(color = guide_legend(override.aes = list(size = 2))))
  dev.off()
  
  
  ### Avoid visual biases from default ordering by shuffling the points
  SE <- as.SingleCellExperiment(SCdata.main)
  red.dim <- SingleCellExperiment::reducedDim(SE, "UMAP");
  plot.data <- data.frame(X=red.dim[, 1], Y=red.dim[, 2], row.names=colnames(SE));
  plot.data$ColorBy <- colData(SE)[, "Patient"];
  plot.data[["ColorBy"]] <- factor(plot.data[["ColorBy"]]);
  set.seed(60025);
  plot.data <- plot.data[sample(nrow(plot.data)),,drop=FALSE];
  setwd(plotWD1)
  pdf(paste0("Figure1_",Suffix,"_Panel_c.4.pdf"), height = 2, width = 3)
  print(ggplot() + geom_point(aes(x=X, y=Y, color=ColorBy), alpha=1, plot.data, size=0.1) + labs(x=NULL, y=NULL, color="Patient") + scale_color_manual(values=GroupPalette, na.value='grey50', drop=FALSE) +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(), axis.title.y = element_blank(), text = element_text(size=14), legend.key.size = unit(0.7,"line"), legend.text=element_text(size=8), legend.title = element_text(size=10)) +
    guides(color = guide_legend(override.aes = list(size = 1.2))))
  dev.off()
  rm(SE); rm(red.dim); rm(plot.data)
  }
  
  
  
  
  ######## **************************************************  Panel d: DotPlot  ****************************************************
  ### Main Figures
  
  RUNFigure1d="YES"
  if(RUNFigure1d=="YES"){
    Idents(SCdata.main) <- ToUseCol
    Idents(SCdata.main) <- factor(Idents(SCdata.main), levels= ToUseOrder)
    PlotGenes <- c("OLIG2", "PDGFRA", "SOX2", "GFAP", "PTPRZ1", "NES", "VEGFA", "VIM", "TGFBI", "PDGFRB", "PTPRC", "CD8A", "ITGAM", "CD14", "MSR1", "CD68", "MBP")
    setwd(plotWD1)
    pdf(file = paste0("Figure1_",Suffix,"_Panel_d.pdf"), height = 5, width = 9)
    print(DotPlot(SCdata.main, features = PlotGenes, cols= c("gray80", "red"))  + 
            theme(axis.title.x=element_blank(), axis.title.y = element_blank(), axis.text=element_text(size=13), legend.text=element_text(size=13), legend.title=element_text(size=15),
                  legend.key.size = unit(0.4, "cm")) + RotatedAxis() + scale_colour_viridis_c(option = "plasma"))
    dev.off()
    
  }
  
  
  ######## **************************************************  Panel e: Barplot  ****************************************************
  ### Main Figures
  
  RUNFigure1e="YES"
  if(RUNFigure1e=="YES"){
    column.col=c("#0039d1","#c20f00","purple","#26cc00")
    names(column.col)=c("CNSTM-068", "CNSTM-070", "CNSTM-081", "CNSTM-096")
    column.col2= c("darkblue", "#0073ff" ,"#009FFF" ,"#00EEFF", "#FF0073", "#FF00D3" ,"#C900FF" ,"#FF8C00" ,"#FF5900", "#FF2500" ,"darkred" ,"#5FFF00", "#AFFF00" ,"#73C2B5")
    names(column.col2)=c("Cluster01", "Cluster02", "Cluster03", "Cluster04", "Cluster05", "Cluster06", "Cluster07", "Cluster08", "Cluster09", "Cluster10", "Cluster11", "Cluster12", "Cluster13", "Cluster14")
    
    setwd(plotWD1)
    pdf(file = paste0("Figure1_",Suffix,"_Panel_e.pdf"), height = 5, width = 6)
    print(ggplot(SCdata.main@meta.data, aes(Patient, fill=Custom.Cluster))+geom_bar(stat="count",position="fill",colour = "black")+  scale_fill_manual(values = column.col2)+
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
                  axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=10),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=10),
                  axis.text=element_text(size=6.5, face="bold"), axis.title=element_text(size=12,face="bold"),
                  plot.margin = unit(c(0.1, 0.5,0.2, 0.5),"cm"))+labs(y ="Composition (percentage of cells)", x= NULL)+scale_y_continuous(expand = c(0,0)))
    print(ggplot(SCdata.main@meta.data, aes(Custom.Cluster, fill=Custom.Cluster))+geom_bar(stat="count",colour = "black")+  scale_fill_manual(values = column.col2)+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
            axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=10),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=10),
            axis.text=element_text(size=6.5, face="bold"), axis.title=element_text(size=12,face="bold"),
            plot.margin = unit(c(0.1, 0.5,0.2, 0.5),"cm"))+labs(y ="number of cells", x= NULL)+ scale_y_continuous(expand=c(0,0),trans ="sqrt",limits=c(0,20000),breaks=c(1,250,1000,2250,4000,6250,9000,12250,16000)))
    dev.off()
    
    
  }
  
  
}
