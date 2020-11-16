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
RDSname <- paste0("Seurat_Formatted_Normalized_USE_GBM_",ClusterGroup)
pkWD <- "/Users/kumarpa/Desktop/Work/Jax/Kyuson/Manuscript_Draft"
RDSdir <- paste0("/Users/kumarpa/Desktop/Work/Jax/Kyuson/Bill_Pipeline_Data/Data_Input/",ClusterGroup)
ColToUse <- "Custom.Cluster"
Suffix <- "PatientALL"
## Change Figure1 as well
NameInPdf.main <- paste0(ColToUse,"_Based_",Suffix)
IdentToSubsetColName="Custom.Cluster"
IdentToSubset="NA"
downsampleHeatmap=25000

FDR=0.1
FoldChangeCutoff=1.4

source("/Users/kumarpa/Desktop/Work/Jax/Scripts/Scripts_Kyuson/Functions_Kyuson_Project_GitHub.R")

setwd(pkWD)
plotWD <- paste(getwd(),paste0("For_Github_Main_Figures"),sep="/"); print(plotWD)
dir.create(file.path(getwd(),paste0("For_Github_Main_Figures")), showWarnings = FALSE)

setwd(plotWD)
plotWD1 <- paste(getwd(),paste0("Figure1_",Suffix),sep="/"); print(plotWD1)
dir.create(file.path(getwd(),paste0("Figure1_",Suffix)), showWarnings = FALSE)

#Custom.Cluster = c(`Cluster01` = "#0000ee", `Cluster02` = "#27408B", `Cluster03` = "#56B4E9",  `Cluster14` = "#00ffff", `Cluster04` = "#00ff00", `Cluster05` = "#DEB887", `Cluster06` = "#008000", `Cluster07` = "#bf3eff", `Cluster08` = "#ff0000", `Cluster09` = "#8b0000", `Cluster10` = "#CD5C5C", `Cluster11` = "#ff80bf", `Cluster12` = "#F0E442", `Cluster13` = "#ffff00"),
ClusOrder.main <- c("Cluster01", "Cluster02", "Cluster03", "Cluster04", "Cluster05", "Cluster06", "Cluster07", "Cluster08", "Cluster09", "Cluster10", "Cluster11", "Cluster12", "Cluster13", "Cluster14"); ClusOrder.main
ClusOrder <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"); ClusOrder
cbPalette.Cluster <- c("darkblue", "#0073ff" ,"#009FFF" ,"#00EEFF", "#FF0073", "#FF00D3" ,"#C900FF" ,"#FF8C00" ,"#FF5900", "#FF2500" ,"darkred" ,"#5FFF00", "#AFFF00" ,"#FFFF00")
#, "#CAFF70", "#b5d4ff", "#000000", "#787878", "#999999", "#D3D3D3"
# Has 20 colors cbPalette.Cluster <- c("#0000ee", "#27408B", "#56B4E9", "#ff0000", "#CD5C5C", "#DEB887", "#008000", "#E69F00", "#00ffff", "#F0E442", "#bf3eff", "#00ff00", "#ff80bf", "#ffff00", "#CAFF70", "#b5d4ff", "#000000", "#787878", "#999999", "#D3D3D3")
#cbPalette.CT <- c("#0000ee", "#ff0000", "#008000", "#E69F00", "#00ffff", "#F0E442", "#bf3eff", "#00ff00", "#ff80bf", "#00EE00", "#FF0000", "#787878", "#999999", "#D3D3D3", "#e066ff", "#b452cd")
#cbPalette.CT <- c("#00008b", "#87ceff", "#009acd", "#8470ff", "#cd5b45", "#ff7f50", "#ff0000", "#8b0000", "#cdad00", "#FFD700", "#eee685", "#FFFF00", "#008b00", "#00ff00", "#6b8e23", "#20B2AA")
CTOrder <- c("Glioma", "Pericyte", "Tcells", "Myeloid", "Oligo", "Mix"); CTOrder
cbPalette.CT <- c("darkblue", "#FF00D3", "#C900FF", "darkred", "#5FFF00", "grey")
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
SCdata.main$Custom.Cluster <- gsub(" ", "",SCdata.main$Custom.Cluster)
SCdata.main$prediction <- gsub("tumor", "aneuploid",SCdata.main$prediction)
SCdata.main$prediction <- gsub("normal", "diploid",SCdata.main$prediction)
SCdata.main$prediction <- factor(SCdata.main$prediction, levels = PredictionOrder)
Idents(SCdata.main) <- ColToUse
DefaultAssay(SCdata.main) <- "RNA"


########################################################################################################################
################################*********************************#######################################################
RUNFigure1="YES"
if(RUNFigure1=="YES"){
  ### Figure1 Clsuter and heatmaps  #"#787878", "#9A6324", "#ee1289", 
  ToUsePallete <- cbPalette.Cluster
  ToUseOrder <- ClusOrder.main
  ToUseCol <- ColToUse
  
  MarkerDir <- paste0("/Users/kumarpa/Desktop/Work/Jax/Kyuson/Bill_Pipeline_Data/DGEs/DGEs_",ToUseCol,"_Based_",Suffix)
  setwd(MarkerDir)
  MarkerGenes <- read.table(file = paste0("DEGs_Heatmap_",ToUseCol,"_Based_",Suffix,".txt"), header = T, sep = "\t"); head(MarkerGenes); dim(MarkerGenes)
  
  
  
  ######## **************************************************  Panel b: Heatmap  ****************************************************
  ### Main Figures
  # Panel b: Heatmap
  RUNFigure1b="YES"
  if(RUNFigure1b=="YES"){
  markers <- MarkerGenes
  #markers <- Remove_Genes_Rp_mt_Rna_Human(markers)
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
  SCdata.temp.Heatmap <- subset(SCdata.main, downsample=downsampleHeatmap)
  #SCdata.temp.Heatmap <- subset(SCdata.main, downsample=10)
  Idents(SCdata.main) <- ToUseCol
  Idents(SCdata.main) <- factor(Idents(SCdata.main), levels= ToUseOrder)
  SCdata.temp.Heatmap <- ScaleData(object = SCdata.temp.Heatmap, verbose = FALSE, features = markers$gene, scale.max = 2)
  dtype="scale.data"
  nor.exp <- GetAssayData(object = SCdata.temp.Heatmap, slot = dtype); print(dim(nor.exp))
  UseGenes <- intersect(markers$gene, rownames(nor.exp)); length(UseGenes)
  nor.exp <- nor.exp[UseGenes,,drop=FALSE]; dim(nor.exp)
  table(SCdata.temp.Heatmap@meta.data[,ToUseCol], SCdata.temp.Heatmap@meta.data[,GroupName])
  meta.data.plot <- SCdata.temp.Heatmap@meta.data[,c(GroupName, ColToUse)]
  meta.data.plot[,ColToUse] <- factor(meta.data.plot[,ColToUse], levels = ToUseOrder)
  meta.data.plot[,GroupName] <- factor(meta.data.plot[,GroupName], levels = GroupOrder.temp)
  meta.data.plot <- meta.data.plot[order(factor(meta.data.plot[,ColToUse], levels = ToUseOrder), factor(meta.data.plot[,GroupName], levels = GroupOrder.temp)),]; head(meta.data.plot)
  meta.data.plot[,ColToUse] <- as.character(meta.data.plot[,ColToUse])
  nor.exp <- nor.exp[,rownames(meta.data.plot), drop=FALSE]
  colnames(meta.data.plot) <- c("Patient", "Cluster")
  print(dim(nor.exp))
  
  Custom.Cluster.FULL =  cbPalette.Cluster[1:length(ClusOrder.main)]; names(Custom.Cluster.FULL) <- ClusOrder.main; Custom.Cluster.FULL
  Fragment.FULL = FragPalette[1:length(FragPalette)]; names(Fragment.FULL) <- FragPalette; Fragment.FULL
  prediction.FULL = PredictionPalette[1:length(PredictionOrder)]; names(prediction.FULL) <- PredictionOrder; prediction.FULL
  Patient.FULL = GroupPalette.temp[1:length(GroupOrder.temp)]; names(Patient.FULL) <- GroupOrder.temp; Patient.FULL
  
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
  #as.grob(pheatmap(nor.exp[as.character(top$gene),],annotation_col=meta.data.plot,show_colnames=F,show_rownames=T,  color=my_palette, breaks=colors, 
  #                       annotation_colors = ann_colors, cluster_rows = FALSE, cluster_cols = FALSE, main = paste0("Heatmap All markers"), fontsize = 13)); dev.off()
  
  sigGenes = c("PDGFRA", "OLIG1", "GFAP", "SOX2", "ACTA2", "CD3G", "S100A9", "S100A4", "ITGAM", "CD68", "APOE", "HLA-DRA", "P2RY12", "CCL3", "MBP")
  ##heatmap <- Temp$gtable
  ##new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  ##new.label$label
  
  p4 <- Gene.Labels.pheatmap(Temp, kept.labels = sigGenes, repel.degree = 0)
  
  setwd(plotWD1)
  pdf(paste0("Figure1_",ToUseCol,"_",Suffix,"_Panel_b_HEATMAP_UPLOAD.pdf"), height = 8, width = 6)
  print(plot_grid(p4))
  dev.off()
  
  }
  
  
  ######## **************************************************  Panel c: UMAP  ****************************************************
  ### Main Figures
  
  RUNFigure1c="YES"
  if(RUNFigure1c=="YES"){
  Idents(SCdata.main) <- ToUseCol
  Idents(SCdata.main) <- factor(Idents(SCdata.main), levels= ToUseOrder)
  r2.1 <- DimPlot(SCdata.main, reduction = "umap", cols = ToUsePallete, label = F, label.size = 1.5, pt.size = 0.01) + 
    theme(axis.title.x=element_blank(), axis.title.y = element_blank(), legend.key.size = unit(0.4,"line"), legend.text=element_text(size=9)) +
    guides(color = guide_legend(override.aes = list(size = 1.2)))
  
  setwd(plotWD1)
  pdf(paste0("Figure1_",ToUseCol,"_UMAP_",Suffix,"_Panel_c.1_Cluster_UPLOAD.pdf"), height = 2, width = 3)
  print(plot_grid(r2.1))
  dev.off()
  
  #r2.2 <- DimPlot(SCdata.main, reduction = "umap", cols = GroupPalette.temp, label = F, label.size = 4, pt.size = 0.01, group.by = "Patient") + 
  #  theme(axis.title.x=element_blank(), axis.title.y = element_blank(), legend.key.size = unit(0.5,"line"), legend.text=element_text(size=9)) +
  #  guides(color = guide_legend(override.aes = list(size = 2)))
  
  CellToPlot <- rownames(SCdata.main@meta.data[SCdata.main@meta.data$prediction %in% c("aneuploid", "diploid"),])
  r2.3 <- DimPlot(SCdata.main, reduction = "umap", cells = CellToPlot, cols = PredictionPalette, label = F, label.size = 4, pt.size = 0.01, group.by = "prediction") + 
    theme(axis.title.x=element_blank(), axis.title.y = element_blank(), legend.key.size = unit(0.5,"line"), legend.text=element_text(size=11)) +
    guides(color = guide_legend(override.aes = list(size = 2)))
  
  setwd(plotWD1)
  pdf(paste0("Figure1_",ToUseCol,"_UMAP_",Suffix,"_Panel_c.2_Status_UPLOAD.pdf"), height = 2, width = 3)
  print(plot_grid(r2.3))
  dev.off()
  
  
  SCdata.main@meta.data$CT <- "FILL"
  table(SCdata.main@meta.data$CT, SCdata.main@meta.data$Custom.Cluster)
  SCdata.main@meta.data[SCdata.main@meta.data$Custom.Cluster %in% c("Cluster01", "Cluster02", "Cluster03", "Cluster04"),"CT"] <- "Glioma"
  SCdata.main@meta.data[SCdata.main@meta.data$Custom.Cluster %in% "Cluster05","CT"] <- "Mix"
  SCdata.main@meta.data[SCdata.main@meta.data$Custom.Cluster %in% "Cluster06","CT"] <- "Pericyte"
  SCdata.main@meta.data[SCdata.main@meta.data$Custom.Cluster %in% "Cluster07","CT"] <- "Tcells"
  SCdata.main@meta.data[SCdata.main@meta.data$Custom.Cluster %in% c("Cluster08", "Cluster09", "Cluster10", "Cluster11"),"CT"] <- "Myeloid"
  SCdata.main@meta.data[SCdata.main@meta.data$Custom.Cluster %in% c("Cluster12", "Cluster13"),"CT"] <- "Oligo"
  SCdata.main@meta.data[SCdata.main@meta.data$Custom.Cluster %in% c("Cluster14"),"CT"] <- "Glioma"
  table(SCdata.main@meta.data$CT, SCdata.main@meta.data$Custom.Cluster)
  Idents(SCdata.main) <- "CT"
  Idents(SCdata.main) <- factor(Idents(SCdata.main), levels= CTOrder)
  r2.4 <- DimPlot(SCdata.main, reduction = "umap", cols = cbPalette.CT[CTOrder %in% unique((SCdata.main$CT))], label = F, label.size = 6, pt.size = 0.01) + 
    theme(axis.title.x=element_blank(), axis.title.y = element_blank(), legend.key.size = unit(0.5,"line"), legend.text=element_text(size=11)) +
    guides(color = guide_legend(override.aes = list(size = 2)))
  
  setwd(plotWD1)
  pdf(paste0("Figure1_",ToUseCol,"_UMAP_",Suffix,"_Panel_c.3_CT_UPLOAD.pdf"), height = 2, width = 3)
  print(plot_grid(r2.4))
  dev.off()
  
  
  ##Code for new dimplot--> from iSEE package
  SE <- as.SingleCellExperiment(SCdata.main)
  red.dim <- SingleCellExperiment::reducedDim(SE, "UMAP");
  plot.data <- data.frame(X=red.dim[, 1], Y=red.dim[, 2], row.names=colnames(SE));
  
  plot.data$ColorBy <- colData(SE)[, "Patient"];
  plot.data[["ColorBy"]] <- factor(plot.data[["ColorBy"]]);
  # Avoid visual biases from default ordering by shuffling the points
  set.seed(60025);
  plot.data <- plot.data[sample(nrow(plot.data)),,drop=FALSE];
  r2.2 <- ggplot() + geom_point(aes(x=X, y=Y, color=ColorBy), alpha=1, plot.data, size=0.1) + labs(x=NULL, y=NULL, color="Patient") + scale_color_manual(values=GroupPalette.temp, na.value='grey50', drop=FALSE) +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(), axis.title.y = element_blank(), text = element_text(size=14), legend.key.size = unit(0.7,"line"), legend.text=element_text(size=8), legend.title = element_text(size=10)) +
    guides(color = guide_legend(override.aes = list(size = 1.2)))
  
  setwd(plotWD1)
  pdf(paste0("Figure1_",ToUseCol,"_UMAP_",Suffix,"_Panel_c.4_Patient_UPLOAD.pdf"), height = 2, width = 3)
  print(plot_grid(r2.2))
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
    pdf(file = paste0("Figure1_",ToUseCol,"_",Suffix,"_Panel_d_DOTPLOT_UPLOAD.pdf"), height = 5, width = 9)
    print(DotPlot(SCdata.main, features = PlotGenes, cols= c("gray80", "red"))  + 
            theme(axis.title.x=element_blank(), axis.title.y = element_blank(), axis.text=element_text(size=13), legend.text=element_text(size=13), legend.title=element_text(size=15),
                  legend.key.size = unit(0.4, "cm")) + RotatedAxis() + scale_colour_viridis_c(option = "plasma"))
    dev.off()
    
  }
  
  
}
