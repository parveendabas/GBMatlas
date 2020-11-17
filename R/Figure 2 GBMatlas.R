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
library(stringr)
library(GOplot)
library(ggpubr)
library(msigdbr)


source("/Users/kumarpa/Desktop/GBMatlas_Code/Functions_GBMatlas.R")

GroupName <- "Patient"
ClusterGroup <- "Cluster_08_09_10_11"
RDSname <- paste0("GBM_",ClusterGroup)
pkWD <- "/Users/kumarpa/Desktop/GBMatlas_Code"
RDSdir <- paste0("/Users/kumarpa/Desktop/GBMatlas_Code/data/")
ColToUse <- "cluster.midres"
Suffix <- "MyeloidCells"
## Change Figure2 as well
NameInPdf.main <- paste0(ColToUse,"_Based_",Suffix)
IdentToSubsetColName="cluster.midres"
IdentToSubset="NA"
FDR=0.1

source("/Users/kumarpa/Desktop/GBMatlas_Code/Functions_GBMatlas.R")

  setwd(pkWD)
  plotWD <- paste(getwd(),paste0("GBMatlas_Main_Figures"),sep="/"); print(plotWD)
  dir.create(file.path(getwd(),paste0("GBMatlas_Main_Figures")), showWarnings = FALSE)

  setwd(plotWD)
  plotWD1 <- paste(getwd(),paste0("Figure2_",Suffix),sep="/"); print(plotWD1)
  dir.create(file.path(getwd(),paste0("Figure2_",Suffix)), showWarnings = FALSE)
  
  ClusOrder.main <- c("Cluster01", "Cluster02", "Cluster03", "Cluster04", "Cluster05", "Cluster06", "Cluster07", "Cluster08", "Cluster09", "Cluster10", "Cluster11", "Cluster12", "Cluster13", "Cluster14"); ClusOrder.main
  ClusterPalette <- c("darkblue", "#0073ff" ,"#009FFF" ,"#00EEFF", "#FF0073", "#FF00D3" ,"#C900FF" ,"#FF8C00" ,"#FF5900", "#FF2500" ,"darkred" ,"#5FFF00", "#AFFF00" ,"#FFFF00")
  MClusOrder <- c("MC1", "MC2", "MC3", "MC4", "MC5", "MC6", "MC7"); MClusOrder
  MClusOrderPalette <- c("darkred", "#008000", "blue", "#00ff00", "#ff8c00", "red",  "#ffff00")
  CTOrder <- c("s-macr1", "a-microglia", "APC", "r-microglia", "s-macr2", "dividing-mac", "M1-mac"); CTOrder
  CTPalette <- c("darkred", "#008000", "blue", "#00ff00", "#ff8c00", "red",  "#ffff00")
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
  RUNFigure2="YES"
  if(RUNFigure2=="YES"){
  
  ToUsePallete <- MClusOrderPalette
  ToUseOrder <- MClusOrder
  ToUseCol <- ColToUse
  
  setwd(RDSdir)
  MarkerGenes <- read.table(file = paste0("DEGs_",Suffix,".txt"), header = T, sep = "\t"); head(MarkerGenes); dim(MarkerGenes)
  
  
  ######## **************************************************  Panel a: HEATMAP  ****************************************************
  ### Main Figures
  RUNFigure2a="YES"
  if(RUNFigure2a=="YES"){
  markers <- MarkerGenes[MarkerGenes$p_val_adj < FDR,]; dim(markers)
  
  topnumber=30
  topFDR <- markers %>% dplyr::arrange(p_val_adj, desc(avg_logFC)) %>% dplyr::group_by(cluster) %>% dplyr::slice(1:topnumber); dim(topFDR)
  
  SCdata.main$Temp <- "Temp"
  Idents(SCdata.main) <- "Temp"
  SCdata.temp.Heatmap <- SCdata.main
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
  
  
  Cluster.FULL =  MClusOrderPalette[1:length(MClusOrder)]; names(Cluster.FULL) <- MClusOrder; Cluster.FULL
  Fragment.FULL = FragPalette[1:length(FragOrder)]; names(Fragment.FULL) <- FragOrder; Fragment.FULL
  prediction.FULL = PredictionPalette[1:length(PredictionOrder)]; names(prediction.FULL) <- PredictionOrder; prediction.FULL
  Patient.FULL = GroupPalette[1:length(GroupOrder)]; names(Patient.FULL) <- GroupOrder; Patient.FULL
  
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
  
  sigGenes = c("FN1", "S100A4", "IBSP", "CCL3", "CCL4", "CD83", "IL1B", "B2M", "CD74", "HLA-DRA", "CD14", "APOE", "TMEM119", "BHLHE41", "SORL1", "CX3CR1", "MIF", "MKI67", "IER3", "NFKBIZ")
  
  p4 <- Gene.Labels.pheatmap(Temp, kept.labels = sigGenes, repel.degree = 0)
  
  setwd(plotWD1)
  pdf(paste0("Figure2_",Suffix,"_Panel_a.pdf"), height = 5, width = 8)
  print(plot_grid(p4))
  dev.off()
  rm(SCdata.temp.Heatmap)
  }
  
  
  ######## **************************************************  Panel b: UMAP  ****************************************************
  ### Main Figures
  RUNFigure2b="YES"
  if(RUNFigure2b=="YES"){
    Idents(SCdata.main) <- ToUseCol
    Idents(SCdata.main) <- factor(Idents(SCdata.main), levels= ToUseOrder)
    setwd(plotWD1)
    pdf(paste0("Figure2_",Suffix,"_Panel_b.pdf"), height = 3, width = 5)
    print(DimPlot(SCdata.main, reduction = "umap", cols = ToUsePallete, label = F, label.size = 1.5, pt.size = 0.01) + 
      theme(axis.title.x=element_blank(), axis.title.y = element_blank(), legend.key.size = unit(0.4,"line"), legend.text=element_text(size=13)) +
      guides(color = guide_legend(override.aes = list(size = 1.2))))
    dev.off()
  }
  
  
  ######## **************************************************  Panel 3d: UMAP  ****************************************************
  ### Main Figures
  
  RUNFigure3d="YES"
  if(RUNFigure3d=="YES"){
    setwd(plotWD1)
    pdf(paste0("Figure3_",Suffix,"_Panel_d.pdf"), height = 3, width = 4)
    Idents(SCdata.main) <- ToUseCol
    Idents(SCdata.main) <- factor(Idents(SCdata.main), levels= ToUseOrder)
    print(VlnPlot(SCdata.main, features = c("S100A4"), cols = ToUsePallete, pt.size = 0.0)  + 
            theme(axis.title.x=element_blank(), axis.title.y = element_blank(), axis.text=element_text(size=13), legend.text=element_text(size=14), legend.title=element_text(size=14),
                  legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = 16), legend.position = 'none') + 
            stat_summary(fun=median, geom="point", size=1, color="black"))
    dev.off()
  }
  
  
  
  ######## **************************************************  Panel c: circular dendrograms  ****************************************************
  ### Main Figures
  RUNFigure2c="YES"
  if(RUNFigure2c=="YES"){
    
    Idents(SCdata.main) <- ToUseCol
    Idents(SCdata.main) <- factor(Idents(SCdata.main), levels= ToUseOrder)
    #Myeloid.Markers= FindAllMarkers(SCdata.main,test.use = "MAST",only.pos = F)
    setwd(RDSdir)
    Myeloid.Markers <- read.csv(file = paste0("MyeloidmarkersallMAST.csv"), header = T, row.names = 1); head(Myeloid.Markers)
    names(Myeloid.Markers)[1] <- "P.Value"
    names(Myeloid.Markers)[2] <- "logFC"
    names(Myeloid.Markers)[5] <- "adj.P.Val"
    names(Myeloid.Markers)[7] <- "ID"
    head(Myeloid.Markers)
    
    ##switch to directory containing gprofiler files... this file is generated using the g:profiler website (https://biit.cs.ut.ee/gprofiler/gost) by uploading the DE genes for each cluster
    setwd(RDSdir)
    filess <- list.files(pattern="gprofiler.csv"); filess
    
    for (file in filess){
      X=read.csv(paste(file))
      NAME= str_replace(paste(file),"gprofiler.csv","")
      assign(NAME, X)
    }
    
    Final=list()
    for (i in 1:7){
      #i=1
      DEgenes=Myeloid.Markers%>%dplyr::filter(cluster==i-1)
      MC=paste0("MC",i)
      circ<-circle_dat(eval(parse(text=MC)),DEgenes)
      process = c("canonical glycolysis","cellular response to hypoxia","Antigen processing-Cross presentation","cellular response to interferon-gamma","cellular response to tumor necrosis factor","NF-kappa B signaling pathway","oxidative phosphorylation","M Phase")
      genes=as.data.frame(DEgenes[,c("ID","logFC")])
      chord=chord_dat(circ,process = process,genes = genes)
      Chord1=chord[,colSums(chord)>0]
      term.col<-c("#73A27E","#73C2B5","#008F91","#194D94","#7D0202","#CC2929","#F47C54","#F9A953")
      names(term.col)=process
      P1=GOCluster(chord,process = process,lfc.col =c("darkred","white","darkblue"),term.col = term.col,lfc.min = -1,lfc.max = 1)+guides(size=guide_legend("GO Terms",title.position ="top",order = 20,ncol=2,byrow=F,override.aes=list(shape=22,fill=term.col,size = 4,linetype = 0)))+
        theme(legend.position='right',legend.background = element_rect(fill='transparent'),legend.spacing.x = unit(0.1, 'cm'),legend.spacing.y = unit(0.1, 'cm'),legend.title = element_text(size=10,face="bold"),
              legend.box='horizontal',legend.box.just = "center",legend.key.height = unit(0.3,'cm'),legend.key.width =unit(0.5,'cm'),legend.direction='horizontal',legend.text = element_text(size=9,face="bold"),plot.margin = unit(c(0,0,0,0),'cm'))+ggtitle(paste0("MC",i))+
        theme(plot.title = element_text(hjust = 0.5,vjust = -5,face = "bold",size = 10),panel.spacing = unit(c(0,0,0,0),'cm'),plot.margin = unit(c(0,-1,-0.5,-1),'cm'),axis.text = element_blank(),axis.ticks.length=unit(0,'cm'))+labs(x=NULL, y=NULL)
      Final[[i]]=ggplotGrob(P1)
      NAME=paste0("Graphs.MC",i)
      assign(NAME, P1)
    }
    g_legend<-function(a.gplot){
      tmp <- ggplot_gtable(ggplot_build(a.gplot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend <- tmp$grobs[[leg]]
      return(legend)}
    
    mylegend<-g_legend(Graphs.MC1)
    lay <- rbind(c(1,2,3,4),
                 c(5,6,7,NA),
                 c(8,8,8,8))
    #png("Figure2c.png",width=6.5,height=5,units = "in",res = 1200,type = "quartz",pointsize = 10)
    setwd(plotWD1)
    pdf(paste0("Figure2_",Suffix,"_Panel_c.pdf"), height = 6, width = 8)
    print(grid.arrange(Graphs.MC1+ theme(legend.position="none"),Graphs.MC2+ theme(legend.position="none"),
                 Graphs.MC3+ theme(legend.position="none"),Graphs.MC4+ theme(legend.position="none"),
                 Graphs.MC5+ theme(legend.position="none"),Graphs.MC6+ theme(legend.position="none"),
                 Graphs.MC7+ theme(legend.position="none"),mylegend, ncol = 4, nrow = 4,layout_matrix=lay))
    dev.off()
    
    
  }
  
  
  
  ######## **************************************************  Panel c: Butterfly Plot  ****************************************************
  ### Main Figures
  RUNFigure2d="YES"
  if(RUNFigure2d=="YES"){
    
    library(scrabble)
    source("/Users/kumarpa/Desktop/GBMatlas_Code/to make suva's butterfly.R")
    
    Idents(SCdata.main) <- ToUseCol
    Idents(SCdata.main) <- factor(Idents(SCdata.main), levels= ToUseOrder)
    Counts <- as.data.frame(GetAssayData(object = SCdata.main, assay = "RNA", slot = "data"))
    #get count matrix
    CellInfoHM=SCdata.main@meta.data
    #get Hallmark gene sets
    m_df_H<- msigdbr(species = "Homo sapiens", category = "H")
    Filtered.genesets=filter(m_df_H,gs_name=="HALLMARK_OXIDATIVE_PHOSPHORYLATION"|gs_name=="HALLMARK_TNFA_SIGNALING_VIA_NFKB"|gs_name=="HALLMARK_INTERFERON_GAMMA_RESPONSE"|gs_name=="HALLMARK_HYPOXIA")
    geneset.all <-   Filtered.genesets[Filtered.genesets$gene_symbol %in% rownames(Counts),]
    geneset <- geneset.all[!duplicated(geneset.all$gene_symbol),]
    genesetlist=geneset.all%>% split(x = .$gene_symbol, f = .$gs_name)
    #get score
    sc=score(mat=Counts,
             groups=genesetlist,
             binmat = NULL,
             bins = NULL,
             controls = NULL,
             bin.control = F,
             center = T,
             nbin = 30,
             n = 100,
             replace = T)
    #get coordinates
    h=hierarchy (sc, quadrants = c("HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_HYPOXIA"), log.scale =T)
    #change names and add colors in metadata
    CellInfoHM$Cluster <- CellInfoHM$cluster.midres
    table(CellInfoHM$Cluster)
    
    Cluster <- CellInfoHM$Cluster
    names(Cluster) <- colnames(x = SCdata.main)
    SCdata.main <- AddMetaData(object = SCdata.main,metadata = Cluster,col.name = 'Cluster')
    
    CellInfoHM$Clustercolor <- "FILL"
    CellInfoHM$Clustercolor[CellInfoHM$`cluster.midres`== "MC1"] <- "darkred"
    CellInfoHM$Clustercolor[CellInfoHM$`cluster.midres`== "MC2"]  <- "darkgreen"
    CellInfoHM$Clustercolor[CellInfoHM$`cluster.midres`== "MC3"] <- "blue"
    CellInfoHM$Clustercolor[CellInfoHM$`cluster.midres`== "MC4"]  <- "green"
    CellInfoHM$Clustercolor[CellInfoHM$`cluster.midres`== "MC5"] <- "orange"
    CellInfoHM$Clustercolor[CellInfoHM$`cluster.midres`== "MC6"] <- "red"
    CellInfoHM$Clustercolor[CellInfoHM$`cluster.midres`== "MC7"] <- "yellow"
    table(CellInfoHM$Clustercolor)
    
    Clustercolor <- CellInfoHM$Clustercolor
    names(Clustercolor) <- colnames(x = SCdata.main)
    SCdata.main <- AddMetaData(object = SCdata.main,metadata = Clustercolor,col.name = 'Clustercolor')
    #figure 2d
    #Make the big plot--> all clusters
    groups=SCdata.main@meta.data[,c("Clustercolor","Cluster")]
    matrix=merge(h,groups,by.x=0,by.y=0,all.x = T)
    matrix$Clustercolor[is.na(matrix$Clustercolor)] <- "gray"
    matrix$Cluster[is.na(matrix$Cluster)] <- "Other"
    row.names(matrix)=matrix$Row.names
    x=matrix$Clustercolor
    y=matrix$Cluster
    col=x[!duplicated(x)]
    names(col)=y[!duplicated(y)]
    matrix=matrix[,-1]
    title="All Myeloid Clusters"
    p0 <- ggplot(matrix, aes(x = X,
                             y =Y,color=factor(Cluster)))+geom_point()+geom_point(data = subset(matrix, Cluster !="Other"))+ 
      scale_color_manual(values=col,aesthetics = "colour", breaks = waiver()) +labs(x=NULL,y=NULL)+theme(legend.position = "none")+
      scale_x_continuous(expand = c(0, 0), limits = c(-1,1)) + scale_y_continuous(expand = c(0, 0), limits = c(-0.8,0.8))+
      theme(panel.background = element_rect(fill = "white",colour = "white"),axis.ticks.x=element_blank(),axis.text.x=element_blank(),
            axis.ticks.y=element_blank(), axis.text.y=element_blank())+
      ggtitle(title)+theme(plot.title =element_text(size=25,face="bold") )+ 
      theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+
      geom_hline(yintercept=0, color = "black", size=0.5)+
      geom_vline(xintercept=0, color = "black", size=0.5)+
      annotate("rect", xmin = -1, xmax = -0.3, ymin = 0.6, ymax = 0.8, fill= "black")  + 
      annotate("text",x = -0.65, y = 0.7,label = "IFNG-res",color="white",fontface="bold",size=9)+ 
      annotate("rect", xmin = 0.3, xmax = 1, ymin = 0.6, ymax = 0.8, fill= "black")  +
      annotate("text",x = 0.65, y = 0.7,label = "Hypoxia",color="white",fontface="bold",size=9)+
      annotate("rect", xmin = -1, xmax = -0.3, ymin = -0.8, ymax = -0.6, fill= "black")  + 
      annotate("text",x = -0.65, y = -0.7,label = "Ox-Phos",color="white",fontface="bold",size=9)+ 
      annotate("rect", xmin = 0.3, xmax =1, ymin = -0.8, ymax = -0.6, fill= "black")  +
      annotate("text",x = 0.65, y = -0.7,label = "TNFA-sig",color="white",fontface="bold",size=9)  
    P0=ggplotGrob(p0)
    ##make the small plots--> one per cluster
    Final <- list()
    Final[[1]]=ggplotGrob(p0)
    for (i in 1:7) {
      ClusterMD=SCdata.main@meta.data[SCdata.main@meta.data$Cluster==paste0("MC",i),]
      groups=ClusterMD[,c("Clustercolor","Cluster")]
      title=paste0("MC",i)
      matrix=merge(h,groups,by.x=0,by.y=0,all.x = T)
      matrix$Clustercolor[is.na(matrix$Clustercolor)] <- "gray"
      matrix$Cluster[is.na(matrix$Cluster)] <- "Other"
      row.names(matrix)=matrix$Row.names
      x=matrix$Clustercolor
      y=matrix$Cluster
      col=x[!duplicated(x)]
      names(col)=y[!duplicated(y)]
      matrix=matrix[,-1]
      P=ggplot(matrix, aes(x = X,
                           y =Y,color=factor(Cluster)))+geom_point()+geom_point(data = subset(matrix, Cluster !="Other"))+ 
        scale_color_manual(values=col,aesthetics = "colour", breaks = waiver()) +labs(x=NULL,y=NULL)+theme(legend.position = "none")+
        scale_x_continuous(expand = c(0, 0), limits = c(-1,1)) + scale_y_continuous(expand = c(0, 0), limits = c(-0.8,0.8))+
        theme(panel.background = element_rect(fill = "white",colour = "white"),axis.ticks.x=element_blank(),axis.text.x=element_blank(),
              axis.ticks.y=element_blank(), axis.text.y=element_blank())+
        ggtitle(title)+theme(plot.title =element_text(size=25,face="bold") )+ 
        theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+
        geom_hline(yintercept=0, color = "black", size=0.5)+
        geom_vline(xintercept=0, color = "black", size=0.5)+
        annotate("rect", xmin = -1, xmax = -0.3, ymin = 0.6, ymax = 0.8, fill= "black")  + 
        annotate("text",x = -0.65, y = 0.7,label = "IFNG-res",color="white",fontface="bold",size=6)+ 
        annotate("rect", xmin = 0.3, xmax = 1, ymin = 0.6, ymax = 0.8, fill= "black")  +
        annotate("text",x = 0.65, y = 0.7,label = "Hypoxia",color="white",fontface="bold",size=6)+
        annotate("rect", xmin = -1, xmax = -0.3, ymin = -0.8, ymax = -0.6, fill= "black")  + 
        annotate("text",x = -0.65, y = -0.7,label = "Ox-Phos",color="white",fontface="bold",size=6)+ 
        annotate("rect", xmin = 0.3, xmax =1, ymin = -0.8, ymax = -0.6, fill= "black")  +
        annotate("text",x = 0.65, y = -0.7,label = "TNFA-sig",color="white",fontface="bold",size=6)  
      Final[[i+1]] = ggplotGrob(P)
    }
    Final[[1]]=ggplotGrob(p0)
    
    setwd(plotWD1)
    pdf(paste0("Figure2_",Suffix,"_Panel_d.pdf"), height = 6, width = 20)
    grid.arrange(grobs=Final, widths = c(1.5,1,1,1,1),layout_matrix = rbind(c(1, 2, 3,4,5),c(1,6,7,8,NA))) 
    dev.off()
    
    
    
  }
  
  
  
  RUNFigure2EF="YES"
  if(RUNFigure2EF=="YES"){
    
    library(Seurat)
    library(scrabble)
    ##to get figure 2e
    #Load Dataset
    setwd(RDSdir)
    CGGA.EXP <- read.table("CGGA.mRNAseq_325.RSEM-genes.20200506.txt",as.is=T,header = T,row.names=1)
    AnnotationCGGA <- read.csv("2020-08-20_CGGA_pheno.csv",as.is=T,header = T,row.names=1)
    CGGA.GBM=filter(AnnotationCGGA,Histology=="GBM")
    #got markers
    Myeloid.Markers= read.csv("Myeloidmarkers.csv")
    
    #remove genes from signature gene lists that do not exist in expression matrix ##otherwise the code won't work
    Counts=as.data.frame(CGGA.EXP)
    geneset.all <-   Myeloid.Markers[ Myeloid.Markers$gene %in% rownames(Counts),]
    geneset <- geneset.all[!duplicated(geneset.all$gene),]
    genesetlist=geneset.all%>% split(x =.$gene, f = .$cluster,drop = F)
    Scounts=Counts[geneset$gene,]
    #get scores 
    sc=score(mat=Scounts,
             groups=genesetlist,
             binmat = NULL,
             bins = NULL,
             controls = NULL,
             bin.control = F,
             center = T,
             nbin = 30,
             n = 100,
             replace = T)
    #write.csv(sc,"CGGAscores.csv")
    #scores will slightly chage everytime you run them because of random sampling of control genes
    sc=read.csv("CGGAscores.csv",row.names = 1)
    ##Survival Analysis
    
    library(survival)
    library(survminer)
    CGGA.sc=merge(sc,AnnotationCGGA,by.x=0, by.y=0)
    CGGA.sc$OS<-as.numeric(CGGA.sc$survival)
    CGGA.sc$status<-as.numeric(CGGA.sc$status)
    surv_object <- Surv(time = CGGA.sc$OS, event = CGGA.sc$status)
    #for loop to assign "positive" or "negative" to cell scores then draw survival plots
    data.subset=colnames(sc)
    data.subset=sort(data.subset, decreasing = FALSE)
    CGGA.sc2 =CGGA.sc
    Final <- list()
    for( i in 1: length(data.subset))
    {
      YY= data.subset[i]
      OO=paste(YY,"Levels",sep=".")
      CGGA.sc2 <- CGGA.sc2%>% mutate(OO = ifelse(CGGA.sc2[YY]>=0, "Positive", "Negative"))
      CGGA.sc2$OO <-factor(CGGA.sc2$OO)
      fit1 <- survfit(surv_object ~OO, data = CGGA.sc2)
      XX <- ggsurvplot(fit1, data = CGGA.sc2, pval = TRUE,pval.coord = c(50, 1),pval.size=6,legend="top",palette = c("Blue","Red"),risk.table = F,font.main = c(12, "bold"),legend.title = "status",
                       title=YY, legend.labs = c("Not-Enriched", "Enriched"),font.x = c(14, "bold"),font.y = c(14, "bold"),font.tickslab = c(12, "plain"),xlab="Time (Months)")
      Final[[i]] = XX
    }
    
    setwd(plotWD1)
    #pdf(file ="Figure2E.pdf", height = 3.5, width =14,onefile = T)
    pdf(paste0("Figure2_",Suffix,"_Panel_e.pdf"), height = 3.5, width = 14)
    arrange_ggsurvplots(Final,print = TRUE, ncol = 4,nrow = 1,surv.plot.height = NULL,risk.table.height = NULL,ncensor.plot.height = NULL)
    dev.off()
    
    ### to get figure 2f- I subset the matrix and only get the scores for GBM patients
    x=intersect(row.names(CGGA.GBM),colnames(Counts))
    Scounts=Scounts[,x]
    #get scores ###don't forget to  load the functions from ##scrabble package into your environment
    sc=score(mat=Scounts,
             groups=genesetlist,
             binmat = NULL,
             bins = NULL,
             controls = NULL,
             bin.control = F,
             center = T,
             nbin = 30,
             n = 100,
             replace = T)
    #write.csv(sc,"GBMonlyCGGAscores.csv")
    #Since everytime you run the scoring the results are slightly different I provided the scoring file that reproduces the exact figures in the paper
    setwd(RDSdir)
    sc=read.csv("GBMonlyCGGAscores.csv",row.names = 1)
    ##Survival Analysis
    CGGA.sc=merge(sc,AnnotationCGGA,by.x=0, by.y=0)
    CGGA.sc$OS<-as.numeric(CGGA.sc$survival)
    CGGA.sc$status<-as.numeric(CGGA.sc$status)
    surv_object <- Surv(time = CGGA.sc$OS, event = CGGA.sc$status)
    #for loop to assign "positive" or "negative" to cell scores then draw survival plots
    data.subset=colnames(sc)
    data.subset=sort(data.subset, decreasing = FALSE)
    CGGA.sc2 =CGGA.sc
    Final <- list()
    for( i in 1: length(data.subset))
    {
      YY= data.subset[i]
      OO=paste(YY,"Levels",sep=".")
      CGGA.sc2 <- CGGA.sc2%>% mutate(OO = ifelse(CGGA.sc2[YY]>=0, "Positive", "Negative"))
      CGGA.sc2$OO <-factor(CGGA.sc2$OO)
      fit1 <- survfit(surv_object ~OO, data = CGGA.sc2)
      XX <- ggsurvplot(fit1, data = CGGA.sc2, pval = TRUE,pval.coord = c(50, 1),pval.size=6,legend="top",palette = c("Blue","Red"),risk.table = F,font.main = c(12, "bold"),legend.title = "status",
                       title=YY, legend.labs = c("Not-Enriched", "Enriched"),font.x = c(14, "bold"),font.y = c(14, "bold"),font.tickslab = c(12, "plain"),xlab="Time (Months)")
      Final[[i]] = XX
    }
    
    setwd(plotWD1)
    pdf(file =paste0("Figure2_",Suffix,"_Panel_f.pdf"), height = 3.5, width =14,onefile = T)
    arrange_ggsurvplots(Final,print = TRUE, ncol = 4,nrow = 1,surv.plot.height = NULL,risk.table.height = NULL,ncensor.plot.height = NULL)
    dev.off()
    
  }
  
  
  
  
  RUNFigure2G="YES"
  if(RUNFigure2G=="YES"){
  
    library(easypackages)
    
    MyPackages<-c("dplyr","pheatmap","ggplot2","ggpubr","gridExtra","viridis",
                  "grid","lattice","gtools","Biobase","edgeR","limma",
                  "RColorBrewer","biomaRt","GSEABase","genefilter","gridtext",
                  "GSVA","msigdbr","gage","fgsea","Seurat","SummarizedExperiment","ComplexHeatmap")
    libraries(MyPackages)
    
    #load data
    setwd(RDSdir)
    myeloidneftel.subset=readRDS(paste0("nefteletalGBM-allclusters_subset.rds"))
    #saveRDS(myeloidneftel.subset, paste0("nefteletalGBM-allclusters_subset.rds"))
    
    myeloid.metadata=myeloidneftel.subset@meta.data
    
    #heatmap from 2g
    mneftel.sce=as.SingleCellExperiment(myeloidneftel.subset)
    
    setwd(RDSdir)
    MCsignaturegenes=read.csv("Myeloidmarkers.csv")
    top30.Myeloid.Markers=MCsignaturegenes %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
    genes<-top30.Myeloid.Markers[!duplicated(top30.Myeloid.Markers$gene),]
    genes$cnum=as.numeric(gsub("Cluster", "", genes$cluster))
    genes=genes[order(genes$cnum),]
    genes$Cluster=gsub("Cluster", "MC", genes$cluster)
    plot.data<-as.data.frame(assay(mneftel.sce, "logcounts"))
    plot.data<-plot.data[genes$gene,]
    plot.data <- plot.data - rowMeans(plot.data)
    plot.data=na.omit(plot.data)
    col<- circlize::colorRamp2(breaks = c(-4, 0, 4), colors = c("#007dd1", "white", "#ab3000"))
    genes=as.data.frame(genes)
    row.names(genes)=genes$gene
    row_annot=genes[rownames(plot.data),"Cluster",drop=F]
    row_colors=list()
    row_cols <-c("darkred","darkgreen","blue","green","orange","red","yellow")
    names(row_cols)=c("MC1","MC2","MC3","MC4","MC5","MC6","MC7")
    row_colors[["MC"]] <-row_cols
    MC=as.matrix(row_annot[,c("Cluster")])
    rowanno <- rowAnnotation(MC=as.matrix(row_annot[,c("Cluster")]),show_annotation_name = FALSE,col=row_colors,show_legend = F)
    column_annot <-myeloid.metadata[,c("Patient"),drop=F]
    column_annot=with(column_annot, column_annot[order(Cluster), , drop=F])
    plot.data<-plot.data[,row.names(column_annot)]
    column.col=brewer.pal(9,"Set3")
    names(column.col)=c("X102" ,"X105A", "X114", "X115" ,"X118", "X124" ,"X125","X126" ,"X143")
    
    column.colors=list()
    column.colors[["Patient"]]<-column.col
    Patient=as.matrix(column_annot[,c("Patient"),drop=F])
    colanno <- columnAnnotation(df=column_annot,
                                show_annotation_name =F,show_legend = F,col=column.colors)
    
    HM=Heatmap(name="logcounts (centered)",as.matrix(plot.data),cluster_rows = F,cluster_columns = T,left_annotation = rowanno,top_annotation = colanno,
               col = col,show_column_names= FALSE,show_row_names=F,row_title = "Our myeloid cluster signature genes",row_split = MC,
               column_title = gt_render("Neftel *et al.* GBM dataset"),column_km = 7,
               column_dend_reorder = T,border = "gray",show_heatmap_legend = F,use_raster = T)        
    lgd=Legend(title = "logcounts", at=  c(-4,-2,0, 2, 4), col=col)
    lgd1=Legend(labels = levels(as.factor(Patient)),title="Patient",legend_gp = gpar(fill=column.col))
    lgd2=Legend(labels = levels(as.factor(row_annot$Cluster)),title="MC",legend_gp = gpar(fill=c("darkred","darkgreen","blue","green","orange","red","yellow")))
    
    setwd(plotWD1)
    pdf(file =paste0("Figure2_",Suffix,"_Panel_g.pdf"), height = 12, width =8)
    draw(HM,heatmap_legend_list = list(lgd,lgd1, lgd2), heatmap_legend_side = "right",legend_labels_gp =gpar(col = "black", font = 3,face="bold"))
    dev.off()
    
    
  
  }
  
  
  
  }
  
  
  
  