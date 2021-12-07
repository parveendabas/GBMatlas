library(easypackages)
MyPackages<-c("dplyr","ggplot2","ggpubr","gridExtra","viridis","egg",
              "grid","lattice","gtools","Biobase","RColorBrewer","tibble",
              "Seurat","cowplot","patchwork","stringr","ComplexHeatmap","SingleCellExperiment", "ggmin")
libraries(MyPackages)

#Set working directory
Directory="~/Box/Yun lab projects/scRNAseq Project/Human Glioma/ALLSAMPLESNovember2021/"

#dir.create(Directory)

setwd(Directory)

#Change variables

ObjName= "HumanGlioma-AllSamples-"
Subset="Glioma"
resolution= 0.2
dimz=1:30
abb="GC"

RobjDirectory="~/Box/Yun lab projects/scRNAseq Project/Human Glioma/ALLSAMPLESJuly2021/R_Objects/"
#paste0(Directory,"R_Objects/")

OutputDirectory=paste0(Directory,"Output/",Subset,"/")
#dir.create(OutputDirectory)
featuremapDirectory=paste0(OutputDirectory,"Featuremaps/")
#dir.create(featuremapDirectory)

#Load File
SeuratObj=readRDS(paste0(RobjDirectory,"GliomaClusters-11-3-21 patients renamed.rds"))

#Extract Meta Data
CellInfo <- SeuratObj@meta.data

#Assign colors
SampleColors=c(`LGG-03` = "#003F5C", `LGG-04` = "#008796", `ndGBM-01` = "#2F4B7C", 
               `ndGBM-02` = "#227ED4", `ndGBM-03` = "#665191", `ndGBM-04` = "#9E71C7", 
               `ndGBM-05` = "#A05195", `ndGBM-06` = "#E082C1", `ndGBM-07` = "#D45087", 
               `ndGBM-08` = "#F57689", `ndGBM-09` = "#BD3E3E", `ndGBM-10` = "#F95D6A", 
               `ndGBM-11` = "#E07B18", `rGBM-01` = "#FF7C43", `rGBM-02` = "#FFA600", 
               `rGBM-03` = "#F5D256", `rGBM-04` = "#42A665", `rGBM-05` = "#8DB032"
)

ClusterColors=c(GC01 = "#4BA75E", GC02 = "#003F5C", GC03 = "#235A82", GC04 = "#877EAE", 
                GC05 = "#4367B2", GC06 = "#8DB032", GC07 = "#007D8E", GC08 = "#266AB3", 
                GC09 = "#7B5DA5")

FragmentColors=c(`LGG-03` = "#003F5C", `LGG-04-1` = "#005B72", `LGG-04-2` = "#007789", 
                 `LGG-04-3` = "#087B91", `ndGBM-01-A` = "#1B6486", `ndGBM-01-C` = "#2D4C7C", 
                 `ndGBM-01-D` = "#2A5D9C", `ndGBM-01-F` = "#2572BF", `ndGBM-02-1` = "#2D76C9", 
                 `ndGBM-02-2` = "#4764AE", `ndGBM-02-4` = "#625394", `ndGBM-02-5` = "#795CA3", 
                 `ndGBM-03-1` = "#8F68B9", `ndGBM-03-2` = "#9E6CC0", `ndGBM-03-3` = "#9F5FAC", 
                 `ndGBM-04` = "#9F5398", `ndGBM-05` = "#B460A3", `ndGBM-06` = "#CE74B4", 
                 `ndGBM-07` = "#DE7CBA", `ndGBM-08` = "#D968A3", `ndGBM-09` = "#D5548C", 
                 `ndGBM-10` = "#DD5B87", `ndGBM-11-A` = "#EB6A88", `ndGBM-11-B` = "#EF7082", 
                 `ndGBM-11-C` = "#D95A64", `ndGBM-11-D` = "#C34446", `rGBM-01-A` = "#CD464A", 
                 `rGBM-01-B` = "#E5525B", `rGBM-01-C` = "#F75F64", `rGBM-01-D` = "#ED6A43", 
                 `rGBM-02-2` = "#E37623", `rGBM-02-3` = "#E77B23", `rGBM-02-4` = "#F47B34", 
                 `rGBM-02-5` = "#FF7D3F", `rGBM-03-1` = "#FF8E25", `rGBM-03-2` = "#FF9F0A", 
                 `rGBM-03-3` = "#FCB013", `rGBM-04-1` = "#F8C136", `rGBM-04-2` = "#F0D056", 
                 `rGBM-04-3` = "#AABF5C", `rGBM-04-4` = "#63AE62", `rGBM-05-1` = "#51A85A", 
                 `rGBM-05-2` = "#6FAC46", `rGBM-05-3` = "#8DB032")
TypeColors=c(Astrocytoma = "#003F5C", GBM = "#795192", Oligodendroglioma = "#EC5873", 
                   `Recurrent GBM` = "#FFA600")
GradeColors=c(II = "#FFA600", IV = "#003F5C")


#Make figure 2a
Idents(SeuratObj)=CellInfo$Cluster
markers = FindAllMarkers(SeuratObj,logfc.threshold = 0.5,test.use = "wilcox",only.pos = T )
write.csv(markers,paste0(OutputDirectory,ObjName,Subset," cluster markers ","res",resolution,".csv"))
#markers=read.csv("~/Box/Yun lab projects/scRNAseq Project/Human Glioma/ALLSAMPLESJuly2021/Output/Glioma/HumanGlioma-AllSamples-Glioma cluster markers res0.2.csv",row.names = 1)
markers=markers[order(markers$cluster),]
top20.markers <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top20.markers$cluster=as.character(top20.markers$cluster)
top20.markers=top20.markers[order(top20.markers$cluster),]
set.seed(12)
Idents(SeuratObj)=CellInfo$Sample
subset=subset(SeuratObj,downsample=1000)
DimPlot(subset,split.by = "Sample",group.by = "Cluster")
SeuratObj.sce=as.SingleCellExperiment(subset)
plot.data<-as.data.frame(assay(SeuratObj.sce, "logcounts"))
plot.data<-plot.data[top20.markers$gene,]
plot.data <- plot.data - rowMeans(plot.data)
plot.data=na.omit(plot.data)
CellInfoS=subset@meta.data
column_annot <-CellInfoS[,c("Cluster","Patient","Type"),drop=F]
column_annot$Patient = as.factor(as.character(column_annot$Patient))
column_annot=with(column_annot, column_annot[order(Patient), , drop=F])
column_annot=with(column_annot, column_annot[order(Cluster), , drop=F])
plot.data<-plot.data[,row.names(column_annot)]

column.colors=list()
column.colors[["Patient"]]<-SampleColors
column.colors[["Cluster"]]<-ClusterColors
column.colors[["Type"]]<-TypeColors

Patient=as.matrix(column_annot[,c("Patient"),drop=F])
Cluster=as.matrix(column_annot[,c("Cluster"),drop=F])
Type=as.matrix(column_annot[,c("Type"),drop=F])

colanno <- columnAnnotation(df=column_annot,
                            show_annotation_name =T,show_legend = F,col=column.colors)
top2 <- markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
genes= top2$gene
genes=c(genes,"APOD", "PDGFRA", "SPP1", "CHI3L1","CD74", "VEGFA","GFAP", "ATP5E", "FAM64A", "NMU","OLIG2","OLIG1","GFAP", "NES", "PROM1", "VIM", "ERBB2", "TGFBI", "SOX2", "ERBB3", "EGFR", "SOX9", "S100A4", "PTPRZ"," MBP", "B2M", "HLA-DRA", "CCL3")
genes= unique(genes)
rows=rowAnnotation(sel = anno_mark(at = match(genes,row.names(plot.data)), labels = genes,labels_gp =gpar(col = "black", fontsize = 9)))

col<- circlize::colorRamp2(breaks = c(-2, 0, 2), colors = c("#007dd1", "white", "#ab3000"))
HM=Heatmap(name="logcounts",as.matrix(plot.data),cluster_rows =F,cluster_columns = F,top_annotation = colanno, right_annotation = rows,row_names_gp = gpar(fontsize=5),
           col = col,show_column_names= F,show_row_names=F,border = F,show_heatmap_legend = F,use_raster = T)        
lgd=Legend(title = "logcounts", at=  c(-2,0, 2),col_fun = col)
lgd1=Legend(labels = levels(as.factor(column_annot$Cluster)),title="Cluster",legend_gp = gpar(fill=ClusterColors,fontsize=5))
lgd2=Legend(labels = levels(as.factor(column_annot$Patient)),title="Patient",legend_gp = gpar(fill=SampleColors,fontsize=5))
lgd4=Legend(labels = levels(as.factor(column_annot$Type)),title="Type",legend_gp = gpar(fill=TypeColors,fontsize=5))

pdf(paste0(OutputDirectory,"heatmap",ObjName,Subset,"res",resolution,"top20 genes per cluster Fig2a.pdf"),width=6,height=5)
draw(HM,heatmap_legend_list = list( lgd,lgd1,lgd4, lgd2), heatmap_legend_side = "right")
dev.off()


# Fig2b
U1=DimPlot(SeuratObj, reduction = "umap",group.by = "Cluster" ,cols =  ClusterColors,label=T)+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="Clusters")+
  guides(colour = guide_legend(override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=15),ncol = 4))

U2=DimPlot(SeuratObj, reduction = "umap",group.by = "Patient",cols = SampleColors,shuffle = T )+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="Patients")+
  guides(colour = guide_legend(override.aes = list(size=4),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=10),ncol = 3))

U3=DimPlot(SeuratObj, reduction = "umap",group.by = "Type" ,cols =  TypeColors,label=F)+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="Type")+
  guides(colour = guide_legend(override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=10,face = "bold"),ncol = 2))

U4=DimPlot(SeuratObj, reduction = "umap",group.by = "Grade" ,cols =  GradeColors,label=T)+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="Grade")+
  guides(colour = guide_legend(override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=10,face = "bold"),ncol = 3))

plots2=ggarrange(U1,U2,U3,U4,ncol = 2)

pdf(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"all UMAPs Fig2b.pdf"),width = 8,height = 9,family = "ArialMT")
print(plots2)
dev.off()

#Fig2c
library(msigdbr)
library(scrabble)

CellInfo=SeuratObj@meta.data
Idents(SeuratObj)=CellInfo$Cluster
subset=subset(SeuratObj, downsample=5000)

Counts <- as.matrix(GetAssayData(object = subset, assay = "RNA", slot = "counts"))

#Read the Genelist file #from Neftel et al
Subtype.Markers=read.csv("~/Box/Yun lab manuscripts/GBM Single Cells/Nature Communications revision/RCode/SuvaGBMsubtype genelist.csv")

#get scores
subtypelist=Subtype.Markers%>% split(x = .$human_gene_symbol, f = .$gs_name)
sc=score(Counts,
         groups=subtypelist,
         binmat = NULL,
         bins = NULL,
         controls = NULL,
         bin.control = F,
         center = T,
         nbin = 30,
         n = 100,
         replace = T)
#write.csv(sc,paste0(OutputDirectory,"HumanGliomacellstates20000cells.csv"))
sc=read.csv(paste0(OutputDirectory,"HumanGliomacellstates20000cells.csv"),row.names = 1)
#get coordinates
h=hierarchy (sc, quadrants = c("AC","OPC","Mesenchymal","Proneural"), log.scale = T)
# make plots

##per cluster
CellInfo=SeuratObj@meta.data
for( i in 1: length(names(ClusterColors))){
  name=names(ClusterColors)[i]
  color=ClusterColors[i]
  CellInfo$Clustercolor [CellInfo$Cluster== name] <- color
}
Clustercolor <- CellInfo$Clustercolor
names(Clustercolor ) <- row.names(CellInfo)
SeuratObj <- AddMetaData(object = SeuratObj,metadata = Clustercolor ,col.name = 'Clustercolor')

# make plots
xmax=max(h$X)+(max(h$X)*0.2)
xmin=min(h$X)+(min(h$X)*0.2)
ymax=max(h$Y)+(max(h$Y)*0.6)
ymin=min(h$Y)+(min(h$Y)*0.6)

#Make the big plot--> all Clusters
groups=SeuratObj@meta.data[,c("Clustercolor","Cluster")]
matrix=merge(h,groups,by.x=0,by.y=0,all.x = T)
matrix$Clustercolor[is.na(matrix$Clustercolor)] <- "gray"
matrix$Cluster[is.na(matrix$Cluster)] <- "Other"
row.names(matrix)=matrix$Row.names
x=matrix$Clustercolor
y=matrix$Cluster
col=x[!duplicated(x)]
names(col)=y[!duplicated(y)]
matrix=matrix[,-1]
title="All Glioma Clusters"
p0 <- ggplot(matrix, aes(x = X,
                         y =Y,color=factor(Cluster)))+geom_point()+geom_point(data = subset(matrix, Cluster !="Other"))+ 
  scale_color_manual(values=col,aesthetics = "colour", breaks = waiver()) +labs(x=NULL,y=NULL,colour="Cluster")+theme(legend.position = "none")+ 
  theme(legend.text=element_text(size=15),legend.title =element_text(size=15,face="bold") )+ 
  theme(axis.title.x = element_text(hjust = 0, vjust=-2, colour="black",size=10,face="bold"))+ 
  theme(axis.title.y = element_text(hjust = 0, vjust=4, colour="black",size=10,face="bold"))+  
  scale_x_continuous(expand = c(0, 0), limits = c(xmin,xmax)) + scale_y_continuous(expand = c(0, 0), limits = c(ymin,ymax))+
  theme(panel.background = element_rect(fill = "white",colour = "white"),axis.ticks.x=element_blank(),axis.text.x=element_blank(),
        axis.ticks.y=element_blank(), axis.text.y=element_blank())+
  ggtitle(title)+theme(plot.title =element_text(size=25,face="bold") )+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+
  geom_hline(yintercept=0, color = "black", size=0.5)+
  geom_vline(xintercept=0, color = "black", size=0.5)+
  annotate("rect", xmin = xmin, xmax = xmin+5, ymin = ymax-2, ymax = ymax, fill= "black")  + 
  annotate("text",x = xmin+2.5, y = ymax-1,label = "Mes-Like",color="white",fontface="bold",size=7.5)+ 
  annotate("rect", xmin = xmax-5, xmax = xmax, ymin = ymax-2, ymax = ymax, fill= "black")  +
  annotate("text",x = xmax-2.5, y = ymax-1,label = "NPC-Like",color="white",fontface="bold",size=7.5)+
  annotate("rect", xmin = xmin, xmax = xmin+5, ymin = ymin+2, ymax = ymin, fill= "black")  + 
  annotate("text",x = xmin+2.5, y = ymin+1,label = "AC-Like",color="white",fontface="bold",size=7.5)+ 
  annotate("rect", xmin = xmax-5, xmax =xmax, ymin = ymin+2, ymax = ymin, fill= "black")  +
  annotate("text",x = xmax-2.5, y = ymin+1,label = "OPC-Like",color="white",fontface="bold",size=7.5)  
P0=ggplotGrob(p0)
##make the small plots--> one per Cluster
Final <- list()
Final[[1]]=ggplotGrob(p0)
for (i in 1:length(names(ClusterColors))) {
  ClusterMD=SeuratObj@meta.data[SeuratObj@meta.data$Cluster==paste0(names(ClusterColors)[i]),]
  groups=ClusterMD[,c("Clustercolor","Cluster")]
  title=paste0(names(ClusterColors)[i])
  matrix=merge(h,groups,by.x=0,by.y=0,all.x = T)
  matrix$Clustercolor[is.na(matrix$Clustercolor)] <- "gray"
  matrix$Cluster=as.character(matrix$Cluster)
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
    scale_x_continuous(expand = c(0, 0), limits = c(xmin,xmax)) + scale_y_continuous(expand = c(0, 0), limits = c(ymin,ymax))+
    theme(panel.background = element_rect(fill = "white",colour = "white"),axis.ticks.x=element_blank(),axis.text.x=element_blank(),
          axis.ticks.y=element_blank(), axis.text.y=element_blank())+
    ggtitle(title)+theme(plot.title =element_text(size=25,face="bold") )+ 
    theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+
    geom_hline(yintercept=0, color = "black", size=0.5)+
    geom_vline(xintercept=0, color = "black", size=0.5)+
    annotate("rect", xmin = xmin, xmax = xmin+5, ymin = ymax-2, ymax = ymax, fill= "black")  + 
    annotate("text",x = xmin+2.5, y = ymax-1,label = "Mes-Like",color="white",fontface="bold",size=3)+ 
    annotate("rect", xmin = xmax-5, xmax = xmax, ymin = ymax-2, ymax = ymax, fill= "black")  +
    annotate("text",x = xmax-2.5, y = ymax-1,label = "NPC-Like",color="white",fontface="bold",size=3)+
    annotate("rect", xmin = xmin, xmax = xmin+5, ymin = ymin+2, ymax = ymin, fill= "black")  + 
    annotate("text",x = xmin+2.5, y = ymin+1,label = "AC-Like",color="white",fontface="bold",size=3)+ 
    annotate("rect", xmin = xmax-5, xmax =xmax, ymin = ymin+2, ymax = ymin, fill= "black")  +
    annotate("text",x = xmax-2.5, y = ymin+1,label = "OPC-Like",color="white",fontface="bold",size=3)  
  Final[[i+1]] = ggplotGrob(P)
}
Final[[1]]=ggplotGrob(p0)
numofplots= length(Final)


pdf(file =paste0(OutputDirectory,ObjName,Subset," subtype butterfly by Cluster8-6.pdf"), height = 8, width =12,onefile = T)
grid.arrange(grobs=Final, widths = c(2,1,1,1),layout_matrix = rbind(c(1, 2, 3,4),c(1,5,6,7),c(1,8,9,10))) 
dev.off()



#Figure 2d

library(presto)
library(msigdbr)
library(dplyr)
library(fgsea)
library(tibble)
library(pheatmap)
  
m_df_H<- msigdbr(species = "Homo sapiens", category = cat) 
Annot.pathway2=as.data.frame(levels(as.factor(m_df_H$gs_name)))
names(Annot.pathway2)="pathway"
fgsea_sets<- m_df_H %>% split(x = .$gene_symbol, f = .$gs_name)
Clusters=levels(as.factor(Annot$Cluster))
#run test
Clusters.genes <- wilcoxauc(SeuratObj , 'Cluster')
dplyr::count(Clusters.genes, group)
objname="Glioma"
name="Hallmark"

  for (i in 1: length(Clusters)){
    X=  Clusters[i]
    print(X)
    Genes<-Clusters.genes %>%
      dplyr::filter(group == Clusters[[i]]) %>%
      arrange(desc(auc)) %>% 
      dplyr::select(feature, auc)
    ranks<- deframe(Genes)
    head(ranks)
    
    fgseaRes<- fgsea(fgsea_sets, stats = ranks,minSize=10)
    fgseaResTidy <- fgseaRes %>%
      as_tibble() %>%
      arrange(desc(NES))
    fgseaResTidy %>% 
      dplyr::select(-leadingEdge, -ES, -log2err) %>% 
      arrange(padj) %>% 
      head()
    fgseaResTidy$Enrichment = ifelse(fgseaResTidy$NES > 0, "Up-regulated", "Down-regulated")
    Filtidy<-fgseaResTidy %>% filter(padj < 0.05) 
    filtRes = rbind(Filtidy %>% filter(Enrichment == "Up-regulated") %>% head(n= 10),
                    Filtidy %>% filter(Enrichment == "Down-regulated") %>% head(n= 10))
    Y<-fgseaResTidy
    hmc=as.data.frame(Y)
    hmc=apply(Y,2,as.character) 
    #write.csv(hmc,paste(name, X,".csv"))
    names(Y)[names(Y)=="NES"]=paste(X)
    Annot.pathway<-Y[,c("pathway",paste(X))]
    Annot.pathway2<-merge(Annot.pathway, Annot.pathway2, by.x="pathway", by.y="pathway")
    
    
    
  }
  ##make a heatmap
  rownames(Annot.pathway2)=Annot.pathway2$pathway
  Annot.pathway2=Annot.pathway2[,-1]
  positions=Clusters
  
  Annot.pathway2=Annot.pathway2[,positions]
  Annot.pathway2[is.na(Annot.pathway2)]=0
  pathnames=c(row.names(Annot.pathway2 %>% top_n(n = 2, wt = eval(parse(text=names(Annot.pathway2)[1])))),row.names(Annot.pathway2 %>% top_n(n = -2, wt = eval(parse(text=names(Annot.pathway2)[1])))))
  for (i in 2: length(names(Annot.pathway2))){
    pathnames1=c(row.names(Annot.pathway2 %>% top_n(n = 2, wt = eval(parse(text=names(Annot.pathway2)[i])))),row.names(Annot.pathway2 %>% top_n(n = -2, wt = eval(parse(text=names(Annot.pathway2)[i])))))
    pathnames=c(pathnames,pathnames1)
  }
  pathnames =pathnames[!duplicated(pathnames)]
  Annot.pathway4=Annot.pathway2[pathnames,]
  Annot.pathway4[is.na(Annot.pathway4)]=0
  redblackblu=CustomPalette(low = "blue", high = "red", mid = "white", k = 1000)
  pheatmap(Annot.pathway4,scale = "none",cluster_rows = T,cluster_cols = F,border_color = "black",treeheight_row = 0,legend_breaks = c(-5,-2,0,2,5),
           fontsize =8, cellwidth = 15,cellheight =10,color= plasma(10),angle_col = 90,
           main = paste("                                          ",name,objname," (NES)"),filename = paste("GSEA",name," (NES) top and bottom 2each-Heatmap Fig2c.pdf"))
  
  
  #Fig 2e
  ## Figure 2e
  
  m_df_H<- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  Filtered.genesets=filter(m_df_H,gs_name=="HALLMARK_HYPOXIA"|gs_name=="HALLMARK_INTERFERON_GAMMA_RESPONSE"|gs_name=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"|gs_name=="HALLMARK_MYC_TARGETS_V1")
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
  h=hierarchy (sc, quadrants = c("HALLMARK_HYPOXIA","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_MYC_TARGETS_V1"), log.scale =T)
  
  #add cluster colors to metadata
  
  CellInfo=SeuratObj@meta.data
  for( i in 1: length(names(ClusterColors))){
    name=names(ClusterColors)[i]
    color=ClusterColors[i]
    CellInfo$Clustercolor [CellInfo$Cluster== name] <- color
  }
  Clustercolor <- CellInfo$Clustercolor
  names(Clustercolor ) <- row.names(CellInfo)
  SeuratObj <- AddMetaData(object = SeuratObj,metadata = Clustercolor ,col.name = 'Clustercolor')
  
  #Make the big plot--> all clusters
  
  
  xmax=max(h$X)+(max(h$X)*0.2)
  xmin=min(h$X)+(min(h$X)*0.2)
  ymax=max(h$Y)+(max(h$Y)*0.6)
  ymin=min(h$Y)+(min(h$Y)*0.6)
  
  #Make the big plot--> all Clusters
  groups=SeuratObj@meta.data[,c("Clustercolor","Cluster")]
  matrix=merge(h,groups,by.x=0,by.y=0,all.x = T)
  matrix$Clustercolor[is.na(matrix$Clustercolor)] <- "gray"
  matrix$Cluster[is.na(matrix$Cluster)] <- "Other"
  row.names(matrix)=matrix$Row.names
  x=matrix$Clustercolor
  y=matrix$Cluster
  col=x[!duplicated(x)]
  names(col)=y[!duplicated(y)]
  matrix=matrix[,-1]
  title="All Glioma Clusters"
  p0 <- ggplot(matrix, aes(x = X,
                           y =Y,color=factor(Cluster)))+geom_point()+geom_point(data = subset(matrix, Cluster !="Other"))+ 
    scale_color_manual(values=col,aesthetics = "colour", breaks = waiver()) +labs(x=NULL,y=NULL,colour="Cluster")+theme(legend.position = "none")+ 
    theme(legend.text=element_text(size=15),legend.title =element_text(size=15,face="bold") )+ 
    theme(axis.title.x = element_text(hjust = 0, vjust=-2, colour="black",size=10,face="bold"))+ 
    theme(axis.title.y = element_text(hjust = 0, vjust=4, colour="black",size=10,face="bold"))+  
    scale_x_continuous(expand = c(0, 0), limits = c(xmin,xmax)) + scale_y_continuous(expand = c(0, 0), limits = c(ymin,ymax))+
    theme(panel.background = element_rect(fill = "white",colour = "white"),axis.ticks.x=element_blank(),axis.text.x=element_blank(),
          axis.ticks.y=element_blank(), axis.text.y=element_blank())+
    ggtitle(title)+theme(plot.title =element_text(size=25,face="bold") )+ 
    theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+
    geom_hline(yintercept=0, color = "black", size=0.5)+
    geom_vline(xintercept=0, color = "black", size=0.5)+
    annotate("rect", xmin = xmin, xmax = xmin+4, ymin = ymax-1, ymax = ymax, fill= "black")  + 
    annotate("text",x = xmin+2, y = ymax-0.5,label = "EMT",color="white",fontface="bold",size=7.5)+ 
    annotate("rect", xmin = xmax-4, xmax = xmax, ymin = ymax-1, ymax = ymax, fill= "black")  +
    annotate("text",x = xmax-2, y = ymax-0.5,label = "MYC-Tar",color="white",fontface="bold",size=7.5)+
    annotate("rect", xmin = xmin, xmax = xmin+4, ymin = ymin+1, ymax = ymin, fill= "black")  + 
    annotate("text",x = xmin+2, y = ymin+0.5,label = "Hypoxia",color="white",fontface="bold",size=7.5)+ 
    annotate("rect", xmin = xmax-4, xmax =xmax, ymin = ymin+1, ymax = ymin, fill= "black")  +
    annotate("text",x = xmax-2, y = ymin+0.5,label = "IFNg-res",color="white",fontface="bold",size=7.5)  
  P0=ggplotGrob(p0)
  ##make the small plots--> one per Cluster
  Final <- list()
  Final[[1]]=ggplotGrob(p0)
  for (i in 1:length(names(ClusterColors))) {
    ClusterMD=SeuratObj@meta.data[SeuratObj@meta.data$Cluster==paste0(names(ClusterColors)[i]),]
    groups=ClusterMD[,c("Clustercolor","Cluster")]
    title=paste0(names(ClusterColors)[i])
    matrix=merge(h,groups,by.x=0,by.y=0,all.x = T)
    matrix$Clustercolor[is.na(matrix$Clustercolor)] <- "gray"
    matrix$Cluster=as.character(matrix$Cluster)
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
      scale_x_continuous(expand = c(0, 0), limits = c(xmin,xmax)) + scale_y_continuous(expand = c(0, 0), limits = c(ymin,ymax))+
      theme(panel.background = element_rect(fill = "white",colour = "white"),axis.ticks.x=element_blank(),axis.text.x=element_blank(),
            axis.ticks.y=element_blank(), axis.text.y=element_blank())+
      ggtitle(title)+theme(plot.title =element_text(size=25,face="bold") )+ 
      theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+
      geom_hline(yintercept=0, color = "black", size=0.5)+
      geom_vline(xintercept=0, color = "black", size=0.5)+
      
      annotate("rect", xmin = xmin, xmax = xmin+4, ymin = ymax-1, ymax = ymax, fill= "black")  + 
      annotate("text",x = xmin+2, y = ymax-0.5,label = "EMT",color="white",fontface="bold",size=3)+ 
      annotate("rect", xmin = xmax-4, xmax = xmax, ymin = ymax-1, ymax = ymax, fill= "black")  +
      annotate("text",x = xmax-2, y = ymax-0.5,label = "MYC-Tar",color="white",fontface="bold",size=3)+
      annotate("rect", xmin = xmin, xmax = xmin+4, ymin = ymin+1, ymax = ymin, fill= "black")  + 
      annotate("text",x = xmin+2, y = ymin+0.5,label = "Hypoxia",color="white",fontface="bold",size=3)+ 
      annotate("rect", xmin = xmax-4, xmax =xmax, ymin = ymin+1, ymax = ymin, fill= "black")  +
      annotate("text",x = xmax-2, y = ymin+0.5,label = "IFNg-res",color="white",fontface="bold",size=3)  
    Final[[i+1]] = ggplotGrob(P)
  }
  Final[[1]]=ggplotGrob(p0)
  length(Final)
  pdf(file =paste0(OutputDirectory,ObjName,Subset," pathway butterfly by Cluster8-6.pdf"), height = 8, width =12,onefile = T)
  grid.arrange(grobs=Final, widths = c(2,1,1,1),layout_matrix = rbind(c(1, 2, 3,4),c(1,5,6,7),c(1,8,9,10))) 
  dev.off()
  ##Fig2f
  Idents(SeuratObj)=CellInfo$Cluster
  subset=subset(SeuratObj, downsample=5000)
  
  Counts <- as.matrix(GetAssayData(object = subset, assay = "RNA", slot = "counts"))
  
  #get pathway scores 
  Filtered.genesets=filter(m_df_H,gs_name=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"|gs_name=="HALLMARK_HYPOXIA"|gs_name=="HALLMARK_INTERFERON_GAMMA_RESPONSE"|gs_name=="HALLMARK_MYC_TARGETS_V1"|gs_name=="HALLMARK_OXIDATIVE_PHOSPHORYLATION"|gs_name=="HALLMARK_TNFA_SIGNALING_VIA_NFKB")
  geneset.all <-   Filtered.genesets[Filtered.genesets$gene_symbol %in% rownames(Counts),]
  geneset <- geneset.all[!duplicated(geneset.all$gene_symbol),]
  genesetlist=geneset.all%>% split(x = .$gene_symbol, f = .$gs_name)
  #get score
  PathwayScore=score(mat=Counts,
           groups=genesetlist,
           binmat = NULL,
           bins = NULL,
           controls = NULL,
           bin.control = F,
           center = T,
           nbin = 30,
           n = 100,
           replace = T)
  head(PathwayScore)
  colnames(PathwayScore)=c("EMT","Hypoxia","IFNg-res","MYC-tar","Ox-Phos","TNFa-sig")
  
  #get subtype cell scores
  Subtype.Markers=read.csv("~/Box/Yun lab manuscripts/GBM Single Cells/Nature Communications revision/RCode/SuvaGBMsubtype genelist.csv")
  subtypelist=Subtype.Markers%>% split(x = .$human_gene_symbol, f = .$gs_name)
  SubtypeScore=score(Counts,
           groups=subtypelist,
           binmat = NULL,
           bins = NULL,
           controls = NULL,
           bin.control = F,
           center = T,
           nbin = 30,
           n = 100,
           replace = T)
  head(SubtypeScore)
  colnames(SubtypeScore)=c("AC","Mes","OPC","NPC")
  
  mergedsc=cbind(PathwayScore,SubtypeScore)
  corrsc=cor(mergedsc)
  corrsctrunk=corrsc[1:6,7:10]
  p.mat <- corrplot::cor.mtest(mergedsc)
  p.mattrunk=p.mat$p
  p.mattrunk=p.mattrunk[1:6,7:10]
  
  pdf(file =paste0(OutputDirectory,ObjName,Subset," correlogram fig2e.pdf"))
  corrplot::corrplot(corrsctrunk,type = "full",method = "pie",p.mat = p.mattrunk,insig = "label_sig",sig.level = 0.05,
                     col = colorRampPalette(c("green","black","red"))(5),tl.col = "black",	cl.ratio= 0.5,pch.col	="black")
  dev.off()
  