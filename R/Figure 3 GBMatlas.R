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
RDSname <- paste0("GBM_",ClusterGroup)
pkWD <- "/Users/kumarpa/Desktop/GBMatlas_Code"
RDSdir <- paste0("/Users/kumarpa/Desktop/GBMatlas_Code/data/")
ColToUse <- "cluster.R2"
Suffix <- "TCells"
## Change Figure4 as well
NameInPdf.main <- paste0(ColToUse,"_Based_",Suffix)
IdentToSubsetColName="cluster.R2"
FDR=0.1

source("/Users/kumarpa/Desktop/GBMatlas_Code/Functions_GBMatlas.R")

setwd(pkWD)
plotWD <- paste(getwd(),paste0("GBMatlas_Main_Figures"),sep="/"); print(plotWD)
dir.create(file.path(getwd(),paste0("GBMatlas_Main_Figures")), showWarnings = FALSE)


  setwd(plotWD)
  plotWD1 <- paste(getwd(),paste0("Figure3_",Suffix),sep="/"); print(plotWD1)
  dir.create(file.path(getwd(),paste0("Figure3_",Suffix)), showWarnings = FALSE)
  
  ClusOrder.main <- c("Cluster01", "Cluster02", "Cluster03", "Cluster04", "Cluster05", "Cluster06", "Cluster07", "Cluster08", "Cluster09", "Cluster10", "Cluster11", "Cluster12", "Cluster13", "Cluster14"); ClusOrder.main
  ClusterPalette <- c("darkblue", "#0073ff" ,"#009FFF" ,"#00EEFF", "#FF0073", "#FF00D3" ,"#C900FF" ,"#FF8C00" ,"#FF5900", "#FF2500" ,"darkred" ,"#5FFF00", "#AFFF00" ,"#FFFF00")
  TClusOrder <- c("TC1", "TC2", "TC3", "TC4", "TC5", "TC6", "TC7", "TC8", "TC9"); TClusOrder
  TClusOrderPalette <- c("#008000", "darkred", "#ffff00", "#00ff00", "#C900FF" ,"#00EEFF", "#0073ff", "#FF2500", "#FF00D3")
  CTOrder <- c("CD4 T", "NK1", "Dividing T", "Tregs", "Resting T", "CD8 Trm", "CD8 T", "NK2", "CD8 Tmem"); CTOrder
  cbPalette.CT <- c("#008000", "darkred", "#ffff00", "#00ff00", "#C900FF" ,"#00EEFF", "#0073ff", "#FF2500", "#FF00D3")
  GroupOrder = c("CNSTM-068", "CNSTM-070", "CNSTM-081", "CNSTM-096")
  GroupPalette <- c("#0039d1","#c20f00","purple","#26cc00")
  FragOrder = c("CNSTM-068-A", "CNSTM-068-B", "CNSTM-068-C", "CNSTM-068-D", "CNSTM-070-A", "CNSTM-070-C", "CNSTM-070-D", "CNSTM-070-F", "CNSTM-081-A", "CNSTM-081-B", "CNSTM-081-C", "CNSTM-081-D", "CNSTM-096-1", "CNSTM-096-2", "CNSTM-096-4", "CNSTM-096-5")
  FragPalette <- c("#0000FF", "#0045FF", "#008AFF" ,"#00CFFF", "#FF9A00", "#FF6E00", "#FF1400" ,"darkred" ,"#FF0052", "#FF00A6", "#FF00F9" ,"#B000FF","darkgreen","#30FF00","#B9FF00", "#FFFF00")
  PredictionOrder <- c("diploid", "aneuploid")
  PredictionPalette <- c("#ff7f50", "#009FFF")
  
  
  setwd(RDSdir)
  SCdata.main <- readRDS(paste0(RDSname,".rds"))
  #saveRDS(SCdata.main, paste0(RDSname,".rds"))
  DefaultAssay(SCdata.main) <- "RNA"
  
  
  
  
  ########################################################################################################################
  ################################*********************************#######################################################
  RUNFigure3="YES"
  if(RUNFigure3=="YES"){
  ### Figure3 Clsuter and heatmaps  #"#787878", "#9A6324", "#ee1289", 
  ToUsePallete <- TClusOrderPalette
  ToUseOrder <- TClusOrder
  ToUseCol <- ColToUse
  
  setwd(RDSdir)
  MarkerGenes <- read.table(file = paste0("DEGs_",Suffix,".txt"), header = T, sep = "\t"); head(MarkerGenes); dim(MarkerGenes)
  
  
  ######## **************************************************  Panel a: HEATMAP  ****************************************************
  ### Main Figures
  RUNFigure3a="YES"
  if(RUNFigure3a=="YES"){
    markers <- MarkerGenes[MarkerGenes$p_val_adj < FDR,]; dim(markers)
    
    library(ComplexHeatmap)
    
    HumanTcells.sce=as.SingleCellExperiment(SCdata.main)
    plot.data<-as.data.frame(assay(HumanTcells.sce, "logcounts"))
    plot.data<-plot.data[markers$gene,]
    plot.data <- plot.data - rowMeans(plot.data)
    plot.data=na.omit(plot.data)
    column_annot <-SCdata.main@meta.data[,c("cluster.R2","Patient"),drop=F]
    colnames(column_annot) <- c("Cluster","Patient")
    column_annot=with(column_annot, column_annot[order(Cluster), , drop=F])
    plot.data<-plot.data[,row.names(column_annot)]
    column.col=GroupPalette
    names(column.col)=c("CNSTM-068", "CNSTM-070", "CNSTM-081", "CNSTM-096")
    column.col2= TClusOrderPalette
    names(column.col2)=c("TC1", "TC2", "TC3", "TC4", "TC5", "TC6", "TC7", "TC8","TC9")
    
    column.colors=list()
    column.colors[["Patient"]]<-column.col
    column.colors[["Cluster"]]<-column.col2
    Patient=as.matrix(column_annot[,c("Patient"),drop=F])
    Cluster=as.matrix(column_annot[,c("Cluster"),drop=F])
    colanno <- columnAnnotation(df=column_annot,
                                show_annotation_name =T,show_legend = F,col=column.colors)
    genes=c('IL7R','GNLY','HOPX','NKG7','MKI67','FOXP3',"IL2RA","CTLA4","TIGIT","S100A4","GZMK","CD8A","HLA-DRB1","CD2","IFI6")
    
    rows=rowAnnotation(sel = anno_mark(at = match(genes,row.names(plot.data)), labels = genes,labels_gp =gpar(col = "black", fontsize = 8)))
    
    col<- circlize::colorRamp2(breaks = c(-2, 0, 2), colors = c("#007dd1", "white", "#ab3000"))
    #features=top30.markers.all.HumanTcells$gene
    HM=Heatmap(name="logcounts",as.matrix(plot.data),cluster_rows =F,cluster_columns = F,top_annotation = colanno, right_annotation = rows,row_names_gp = gpar(fontsize=5),
               col = col,show_column_names= F,show_row_names=F,border = F,show_heatmap_legend = T,use_raster = T)        
    lgd=Legend(title = "logcounts", at=  c(-2,0, 2), col=col)
    lgd1=Legend(labels = levels(as.factor(column_annot$Cluster)),title="Cluster",legend_gp = gpar(fill=column.col2,fontsize=20))
    lgd2=Legend(labels = levels(as.factor(Patient)),title="Patient",legend_gp = gpar(fill=column.col,fontsize=20))
    
    setwd(plotWD1)
    pdf(paste0("Figure3_",Suffix,"_Panel_a.pdf"),width=7,height=4)
    draw(HM,heatmap_legend_list = list(lgd1, lgd2), heatmap_legend_side = "right")
    dev.off()
    
    
  }
  
  
  
  ######## **************************************************  Panel b: UMAP  ****************************************************
  ### Main Figures
  RUNFigure3b="YES"
  if(RUNFigure3b=="YES"){
  Idents(SCdata.main) <- ToUseCol
  Idents(SCdata.main) <- factor(Idents(SCdata.main), levels= ToUseOrder)
  
  setwd(plotWD1)
  pdf(paste0("Figure3_",Suffix,"_Panel_b.pdf"), height = 3, width = 5.5)
  
  print(DimPlot(SCdata.main, reduction = "umap", cols = ToUsePallete, label = F, label.size = 3, pt.size = 1.8) + 
    theme(axis.title.x=element_blank(), axis.title.y = element_blank(), legend.key.size = unit(0.45,"line"), legend.text=element_text(size=13)) +
    guides(color = guide_legend(override.aes = list(size = 1.2))))
  dev.off()
  }
  
  
  
  
  
  
  RUNFigure3c="YES"
  if(RUNFigure3c=="YES"){
    Idents(SCdata.main) <- ToUseCol
    Idents(SCdata.main) <- factor(Idents(SCdata.main), levels= ToUseOrder)
    
    TregGenes <- c("FOXP3", "IL2RA", "TIGIT", "CTLA4", "S100A4")
    PlotGenes <- data.frame(Genes=TregGenes, CT=rep("CT",5))
    GeneCol="Genes"; GroupCol="CT"
    row.object <- Voilin_Plot_OneGenePerLine(SCdata.main, PlotGenes, GeneCol, GroupCol, cbPalette.CT, 1) 
    #p1 <- plot_grid(plotlist=row.object, nrow = length(row.object), labels=unique(PlotGenes[,GroupCol]), label_colour="darkred", hjust = -0.15)
    setwd(plotWD1)
    pdf(paste0("Figure3_",Suffix,"_Panel_c.pdf"), height = 7, width = 5)
    print(plot_grid(plotlist=row.object, nrow = length(row.object)))
    dev.off()
    
  }
  
  
  RUNFigure3e="YES"
  if(RUNFigure3e=="YES"){
    
    
    library(survival)
    library(survminer)
    #to make figure 3e
    #load data
    setwd(RDSdir)
    CGGA.EXP <- read.table("CGGA.mRNAseq_325.RSEM-genes.20200506.txt",as.is=T,header = T,row.names=1)
    AnnotationCGGA <- read.csv("2020-08-20_CGGA_pheno.csv",as.is=T,header = T,row.names=1)
    CGGA.GBM=filter(AnnotationCGGA,Histology=="GBM")
    #extract S100A4 expression and merge with annotation
    data <- CGGA.EXP["S100A4",]
    t.data=t(data)
    CGGA.dat=merge(t.data,AnnotationCGGA,by.x=0, by.y=0)
    #get KM plots
    CGGA.dat$OS<-as.numeric(CGGA.dat$survival)
    CGGA.dat$status<-as.numeric(CGGA.dat$status)
    surv_object <- Surv(time = CGGA.dat$OS, event = CGGA.dat$status)
    #hist(CGGA.dat$S100A4)
    median(CGGA.dat$S100A4)
    CGGA.dat <-CGGA.dat %>% mutate(S100A4.Levels = ifelse(S100A4 >=14.91, "High", "Low"))
    CGGA.dat$S100A4.Levels <-factor(CGGA.dat$S100A4.Levels)
    fit1 <- survfit(surv_object ~ S100A4.Levels, data = CGGA.dat)
    summary(fit1)
    P1=ggsurvplot(fit1, data = CGGA.dat, pval = TRUE,pval.coord = c(50, 1),pval.size=6,legend="top",palette = c("Red","Blue"),risk.table = F,font.main = c(12, "bold"),legend.title = "Expression",
                  title="S100A4-AllGliomas-CGGA", legend.labs = c("High", "Low"),font.x = c(14, "bold"),font.y = c(14, "bold"),font.tickslab = c(12, "plain"))
    #get GBM patients survival only
    Counts=as.data.frame(CGGA.EXP)
    data <- Counts["S100A4",]
    t.data=t(data)
    CGGA.dat=merge(t.data,CGGA.GBM,by.x=0, by.y=0)
    CGGA.dat$OS<-as.numeric(CGGA.dat$survival)
    CGGA.dat$status<-as.numeric(CGGA.dat$status)
    surv_object <- Surv(time = CGGA.dat$OS, event = CGGA.dat$status)
    #hist(CGGA.dat$S100A4)
    median(CGGA.dat$S100A4)
    CGGA.dat <-CGGA.dat %>% mutate(S100A4.Levels = ifelse(S100A4 >=38.39, "High", "Low"))
    CGGA.dat$S100A4.Levels <-factor(CGGA.dat$S100A4.Levels)
    fit1 <- survfit(surv_object ~ S100A4.Levels, data = CGGA.dat)
    summary(fit1)
    P2=ggsurvplot(fit1, data = CGGA.dat, pval = TRUE,pval.coord = c(50, 1),pval.size=6,legend="top",palette = c("Red","Blue"),risk.table = F,font.main = c(12, "bold"),legend.title = "Expression",
                  title="S100A4-GBM only-CGGA", legend.labs = c("High", "Low"),font.x = c(14, "bold"),font.y = c(14, "bold"),font.tickslab = c(12, "plain"))
    P=list(P1,P2)
    
    setwd(plotWD1)
    pdf(file =paste0("Figure3_",Suffix,"_Panel_e.pdf"), height = 3.5, width =14,onefile = T)
    arrange_ggsurvplots(P,print = TRUE, ncol = 4,nrow = 1,surv.plot.height = NULL,risk.table.height = NULL,ncensor.plot.height = NULL)
    dev.off()
    
    
    
  }
  
  
  
  
  }
  ################################*********************************#######################################################
  ########################################################################################################################
  
  