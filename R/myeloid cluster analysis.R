library(easypackages)
MyPackages<-c("dplyr","ggplot2","ggpubr","gridExtra","viridis","egg",
              "grid","lattice","gtools","Biobase","RColorBrewer","tibble",
              "Seurat","cowplot","patchwork","stringr","ComplexHeatmap","SingleCellExperiment", "ggmin")
libraries(MyPackages)

#Set working directory
Directory="~/Box/Yun lab projects/scRNAseq Project/Human Glioma/ALLSAMPLESJuly2021//"
RobjDirectory=paste0(Directory,"R_Objects/")
#dir.create(RobjDirectory)

setwd(Directory)

#Load File
SeuratObj=readRDS("~/Box/Kyuson_GBM/2021_GBM_All_Three_Batches_Analysis_USE_THIS/Denovo_DataSets/Myeloid/GBM_ALL_Batches_Correction_subset_Myeloid_minGenes_500_PCA_30_Res_0.1_theta_2_2.rds")
DimPlot(SeuratObj)

ObjName= "HumanGlioma-AllSamples-"
Subset="Myeloid-no Tcells or Astrocytes"
resolution= 0.3
dimz=1:30
abb="MC"

OutputDirectory=paste0(Directory,"Output/",Subset,"/")
dir.create(OutputDirectory)
featuremapDirectory=paste0(OutputDirectory,"Featuremaps/")
dir.create(featuremapDirectory)

CellInfo=SeuratObj@meta.data
CellInfo$Fragment=as.character(CellInfo$Fragment)
levels(as.factor(CellInfo$Fragment))
Patients=levels(as.factor(SeuratObj@meta.data$Sample))
Patientnum=readr::parse_number(Patients)
Patientnum=str_remove(Patientnum,"-")
Patientnum=as.numeric(Patientnum)
for(j in 1:length(Patients)){
  if(Patientnum[j]<10){
    #CellInfo$Sample[CellInfo$Sample==Patients[j]]=str_replace(CellInfo$Sample[CellInfo$Sample==Patients[j]],"-","-0")
    CellInfo$Fragment[CellInfo$Sample==Patients[j]]=str_replace(CellInfo$Fragment[CellInfo$Sample==Patients[j]],"MDAG-","MDAG-0")
  }else{
    CellInfo$Fragment[CellInfo$Sample==Patients[j]]=str_replace(CellInfo$Fragment[CellInfo$Sample==Patients[j]],"MDAG-","MDAG-")

  }
}
CellInfo$Patient=CellInfo$Sample

#Rename Clusters
CellInfo$Cluster=str_replace_all(CellInfo$Cluster,"MC","MC0")
SeuratObj@meta.data=CellInfo


# Add Sample color info
# Sample color
SampleColors= Nour_pal("all")(length(levels(as.factor(SeuratObj$Patient))))
names(SampleColors)=levels(as.factor(SeuratObj$Patient))

#Fragment color
CellInfo=SeuratObj@meta.data
CellInfo$Fragment[which(str_detect(row.names(CellInfo), "57"))] <- "CNSTM-394-3"
SeuratObj@meta.data=CellInfo
FragmentColors= Nour_pal("all",reverse = F)(length(levels(as.factor(SeuratObj@meta.data$Fragment))))
names(FragmentColors)=levels(as.factor(SeuratObj@meta.data$Fragment))

#
CellInfo=SeuratObj@meta.data
Assign= list()
Assign[["GBM"]]=c("CNSTM-070","CNSTM-096","MDAG-04",   "MDAG-06",   "MDAG-07",   "MDAG-09",   "MDAG-10", "CNSTM-394","MDAG-12")
Assign[["Astrocytoma"]]=c("CNSTM-081","MDAG-01","MDAG-11")
Assign[["Recurrent GBM"]]=c("CNSTM-068","CNSTM-375", "CNSTM-379", "CNSTM-390",  "CNSTM-397")
Assign[["Oligodendroglioma"]]=c("MDAG-03")

CellInfo$Type=NA
for(i in 1:length(Assign)){
  CellInfo$Type[CellInfo$ID %in% Assign[[i]]]=names(Assign[i])
}

Assign= list()
Assign[["IV"]]=c("CNSTM-070","CNSTM-096","MDAG-04",   "MDAG-06",   "MDAG-07",   "MDAG-09",   "MDAG-10", "MDAG-12","CNSTM-068","CNSTM-375", "CNSTM-379", "CNSTM-390", "CNSTM-394", "CNSTM-397")
Assign[["II"]]=c("CNSTM-081","MDAG-11","MDAG-03")
Assign[["III"]]=c("MDAG-01")

CellInfo$Grade=NA
for(i in 1:length(Assign)){
  CellInfo$Grade[CellInfo$ID %in% Assign[[i]]]=names(Assign[i])
}
SeuratObj@meta.data=CellInfo

CellInfo=SeuratObj@meta.data
CellInfo$ID=CellInfo$Patient


CellInfo$Patient[which(str_detect(CellInfo$ID, "CNSTM-070"))] <- "ndGBM-01"
CellInfo$Patient[which(str_detect(CellInfo$ID, "CNSTM-096"))] <- "ndGBM-02"
CellInfo$Patient[which(str_detect(CellInfo$ID, "CNSTM-394"))] <- "ndGBM-03"
CellInfo$Patient[which(str_detect(CellInfo$ID, "MDAG-04"))] <- "ndGBM-04"
CellInfo$Patient[which(str_detect(CellInfo$ID, "MDAG-06"))] <- "ndGBM-05"
CellInfo$Patient[which(str_detect(CellInfo$ID, "MDAG-07"))] <- "ndGBM-06"
CellInfo$Patient[which(str_detect(CellInfo$ID, "MDAG-09"))] <- "ndGBM-07"
CellInfo$Patient[which(str_detect(CellInfo$ID, "MDAG-10"))] <- "ndGBM-08"
CellInfo$Patient[which(str_detect(CellInfo$ID, "MDAG-12"))] <- "ndGBM-09"
CellInfo$Patient[which(str_detect(CellInfo$ID, "CNSTM-068"))] <- "rGBM-01"
CellInfo$Patient[which(str_detect(CellInfo$ID, "CNSTM-375"))] <- "rGBM-02"
CellInfo$Patient[which(str_detect(CellInfo$ID, "CNSTM-379"))] <- "rGBM-03"
CellInfo$Patient[which(str_detect(CellInfo$ID, "CNSTM-390"))] <- "rGBM-04"
CellInfo$Patient[which(str_detect(CellInfo$ID, "CNSTM-397"))] <- "rGBM-05"
CellInfo$Patient[which(str_detect(CellInfo$ID, "CNSTM-081"))] <- "LGG-01"
CellInfo$Patient[which(str_detect(CellInfo$ID, "MDAG-01"))] <- "LGG-02"
CellInfo$Patient[which(str_detect(CellInfo$ID, "MDAG-11"))] <- "LGG-03"
CellInfo$Patient[which(str_detect(CellInfo$ID, "MDAG-03"))] <- "LGG-04"

Patient  <- CellInfo$Patient
names(Patient ) <- colnames(x = SeuratObj)
SeuratObj <- AddMetaData(object = SeuratObj,metadata =Patient,col.name = 'Patient')
#---
levels(as.factor(CellInfo$Frag))
id= levels(as.factor(CellInfo$ID))
CellInfo$Frag=NA
for (j in 1:length(id)){
    CellInfo$Frag[CellInfo$ID==id[j]]=str_replace(CellInfo$Fragment[CellInfo$ID==id[j]],id[j],CellInfo$Patient[CellInfo$ID==id[j]])
  }
CellInfo$IDfrag=CellInfo$Fragment
CellInfo$Fragment=CellInfo$Frag
SeuratObj@meta.data=CellInfo

#Type color
TypeColors= Nour_pal("main",reverse = F)(length(levels(as.factor(SeuratObj@meta.data$Type))))
names(TypeColors)=levels(as.factor(SeuratObj@meta.data$Type))

#Grade color
GradeColors= Nour_pal("main",reverse = T)(length(levels(as.factor(SeuratObj@meta.data$Grade))))
names(GradeColors)=levels(as.factor(SeuratObj@meta.data$Grade))
#Run PCA
SeuratObj <- RunPCA(SeuratObj, features = VariableFeatures(object = SeuratObj))
print(SeuratObj[["pca"]], dims = 1:5, nfeatures = 5)
VS=VizDimLoadings(SeuratObj, dims = 1:2, reduction = "pca")
VSD=DimPlot(SeuratObj, reduction = "pca",group.by = "Sample",cols = SampleColors)+theme_min()
pdf(paste0(OutputDirectory,ObjName,Subset," PCA features and plot.pdf"),width = 12,height = 4.5,family = "ArialMT")
grid.arrange(VS[[1]]+theme_min(),VS[[2]]+theme_min(),VSD,nrow=1,widths=c(0.6,0.6,1))
dev.off()
#SeuratObj <- JackStraw(SeuratObj, num.replicate = 100)
#SeuratObj <- ScoreJackStraw(SeuratObj, dims = 1:20)
#JackStrawPlot(SeuratObj, dims = 1:15)
ElbowPlot(SeuratObj,ndims = 25)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

SeuratObj <- CellCycleScoring(SeuratObj, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
DimPlot(SeuratObj,group.by = "Phase")
SeuratObj <- ScaleData(SeuratObj, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(SeuratObj))

#initial clustering
library(harmony)
SeuratObj <- RunHarmony(SeuratObj, c("Fragment","sex","Type"))
SeuratObj <- FindNeighbors(SeuratObj,reduction = "harmony", dims = 1:15)

SeuratObj <- RunUMAP(SeuratObj, reduction = "harmony",dims = 1:15 )

#SeuratObj <- FindNeighbors(SeuratObj, dims = dimz)
resolution=0.25
SeuratObj <- FindClusters(SeuratObj, resolution = 0.16)
DimPlot(SeuratObj,split.by = "Patient",ncol = 6,label = T)

#head(Idents(SeuratObj), 5)
#SeuratObj <- RunUMAP(SeuratObj, dims = dimz)

DimPlot(SeuratObj)

#supplementary fig 5a
CellInfo=SeuratObj@meta.data
#CellInfo$oldclusters=CellInfo$Cluster
#Rename Clusters
clusters=levels(as.factor(SeuratObj@meta.data$seurat_clusters))
for(j in 1:length(clusters)){
  if(j<10){
    CellInfo$Cluster[CellInfo$seurat_clusters==j-1]=paste0(abb,"0",j)
  }else{
    CellInfo$Cluster[CellInfo$seurat_clusters==j-1]=paste0(abb,j)
  }
}
SeuratObj@meta.data=CellInfo
# Cluster color
set.seed(001) # just to make it reproducible

ClusterColors= sample(Nour_pal("all",reverse = T)(length(levels(as.factor(SeuratObj@meta.data$Cluster)))))

names(ClusterColors)=levels(as.factor(SeuratObj@meta.data$Cluster))

DimPlot(SeuratObj, reduction = "umap",group.by = "Fragment",cols =FragmentColors,split.by = "Patient",ncol = 3, pt.size = 0.2)+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 16, y.title = 16,main = 16)+ theme(legend.position = "right",text =element_text(size=15) )+
  guides(colour = guide_legend(title="Clusters", override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=15)))
DimPlot(SeuratObj, reduction = "umap",group.by = "Patient",cols =SampleColors,split.by = "ID",ncol = 3, pt.size = 0.2,label = T)+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 16, y.title = 16,main = 16)+ theme(legend.position = "right",text =element_text(size=15) )+
  guides(colour = guide_legend(title="Clusters", override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=15)))

#CellInfo= CellInfo %>% filter(Cluster!="MC07")
#SeuratObj= subset(SeuratObj, cells=row.names(CellInfo))

# Get number of cells per cluster and per Sample
write.csv(as.matrix(table(SeuratObj@meta.data$Cluster,SeuratObj@meta.data$Patient)),
          file = paste0(OutputDirectory,ObjName,Subset," number of cells per cluster and sample.csv"))
write.csv(as.matrix(table(SeuratObj@meta.data$Cluster,SeuratObj@meta.data$Type)),
          file = paste0(OutputDirectory,ObjName,Subset," number of cells per cluster and Type.csv"))
write.csv(as.matrix(table(SeuratObj@meta.data$Cluster,SeuratObj@meta.data$sex)),
          file = paste0(OutputDirectory,ObjName,Subset," number of cells per cluster and Sex.csv"))
write.csv(as.matrix(table(SeuratObj@meta.data$Cluster,SeuratObj@meta.data$Fragment)),
          file = paste0(OutputDirectory,ObjName,Subset," number of cells per cluster and Fragment.csv"))
write.csv(as.matrix(table(SeuratObj@meta.data$Cluster,SeuratObj@meta.data$Grade)),
          file = paste0(OutputDirectory,ObjName,Subset," number of cells per cluster and Grade.csv"))
#make heatmap
resolution=0.25
Idents(SeuratObj)=CellInfo$Assignment2
markers = FindAllMarkers(SeuratObj,logfc.threshold = 0.25,test.use = "wilcox",only.pos = T )
write.csv(markers,paste0(OutputDirectory,ObjName,Subset," Assignment2 markers ","res",resolution,".csv"))
#markers = read.table(paste0(OutputDirectory,ObjName,Subset," Assignment2 markers ","res",resolution,".csv"))
markers=markers[order(markers$cluster),]
top20.markers <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top20.markers$cluster=as.character(top20.markers$cluster)
top20.markers=top20.markers[order(top20.markers$cluster),]
set.seed(12)
Idents(SeuratObj)=CellInfo$Sample
subset=subset(SeuratObj, downsample=1000)
DimPlot(subset,split.by = "Sample",group.by = "Cluster")
SeuratObj.sce=as.SingleCellExperiment(subset)
plot.data<-as.data.frame(assay(SeuratObj.sce, "logcounts"))
plot.data<-plot.data[top20.markers$gene,]
plot.data <- plot.data - rowMeans(plot.data)
plot.data=na.omit(plot.data)
CellInfoS=subset@meta.data


U1=DimPlot(SeuratObj, reduction = "umap",group.by = "Cluster" ,cols =  column.col2,label=T)+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="Clusters")+
  guides(colour = guide_legend(override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=15),ncol = 5))

pdf(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"Clusters UMAP.pdf"),width = 5.5,height =6,family = "ArialMT")
U1
dev.off()
Genes=c('IL1B','IL6','IL12A','IL23A','CXCL9'
        ,'NOS2','CD80','CD86','TNF','MRC1','IL1R1',
        'CCL17','TGFB1','IGF1','FN1','CCL1','IL10','TNFSF14',
        'MERTK','CXCL13','CD163','VEGFA',"FN1", "S100A4", "IBSP",
        "S100A8","CCL3", "CCL4", "CD83", "IL1B", "IL1A","HLA-DRA",
        "CD74", "CD14", "APOE","B2M", "SPRY1", "BHLHE41", "SORL1",
        "IFNGR1", "MALAT1", "CX3CR1", "BNIP3", "MIF","MKI67", "IER3",
        "NFKBIZ","TMEM119","P2RY12","CX3CR1","SLC2A5","CST3","P2RY13","CCL2","EGR2","CCL4","CD83","EGR3")

pdf(paste0(OutputDirectory,"dotplot ",ObjName,Subset, "by Cluster.pdf"),width=10,height=3.5)
DotPlot(SeuratObj,group.by = "Cluster",dot.scale = 5 ,features = unique(Genes) ,scale = T)+scale_colour_viridis_c(option = "plasma")+
  ggmin::theme_min()+ RotatedAxis()+theme(axis.title = element_blank())
DotPlot(SeuratObj,idents = c("MC01","MC02","MC06","MC07"), features = c("TMEM119","P2RY12","CX3CR1","SLC2A5","CST3","P2RY13","CCL2","EGR2","CCL4","CD83","EGR3"),
        group.by = "Cluster",scale=T)+RotatedAxis()+scale_colour_viridis_c(option = "plasma")+theme(text = element_text(size = 11),axis.title = element_blank())
dev.off()

DotPlot(SeuratObj,group.by = "Cluster" ,features = c('IL1B','IL6','IL12A','IL23A','CXCL9',"NFKBIZ"
                                                        ,'NOS2','CD80','CD86','TNF','MRC1','IL1R1',
                                                        'CCL17','TGFB1','IGF1','FN1','CCL1','IL10','TNFSF14',
                                                        'MERTK','CXCL13','CD163','VEGFA'),
        scale = T)+RotatedAxis()+RotatedAxis()+scale_colour_viridis_c(option = "plasma")
DotPlot(SeuratObj,group.by = "Cluster" ,features = c('ITGAM','ITGAX','CD14','CD33','CD68',"MSR1"
                                                     ,'HLA-DRB1','IFNGR1','CD1C','S100A4','MIF','LYZ',
                                                     'CX3CR1','TMEM119','P2RY12','BHLHE41','SPRY1','NFKBIZ','CCL4',
                                                     'CCL3','ITGA4'),
        scale = T)+RotatedAxis()+RotatedAxis()+scale_colour_viridis_c(option = "plasma")



dev.off()
CellInfo=SeuratObj@meta.data

Assign= list()
Assign[["MDSC/mac"]]=c("MC04")
Assign[["i-microglia"]]=c("MC01")
Assign[["a-microglia"]]=c("MC07")
Assign[["h-microglia"]]=c("MC02")
Assign[["AP-mac"]]=c("MC03")
Assign[["AP-microglia"]]=c("MC06")
Assign[["s-mac"]]=c("MC05")
Assign[["DCs"]]=c("MC08")
Assign[["Proliferating"]]=c("MC09")


CellInfo$Assignment=NA
for(i in 1:length(Assign)){
  CellInfo$Assignment[CellInfo$Cluster %in% Assign[[i]]]=names(Assign[i])
}
SeuratObj@meta.data=CellInfo

Assign= list()
Assign[["MDSCs"]]=c("MDCS?")
Assign[["Microglia"]]=c("a-microglia","h-microglia","mic/mac APC","Proliferating")
Assign[["Macs"]]=c("s-mac1","s-mac2")
Assign[["DCs"]]=c("DCs")
Assign[["Proliferating"]]=c("Proliferating")


CellInfo$Assignment2=NA
for(i in 1:length(Assign)){
  CellInfo$Assignment2[CellInfo$Assignment %in% Assign[[i]]]=names(Assign[i])
}
SeuratObj@meta.data=CellInfo

CellInfo=SeuratObj@meta.data

Assign= list()
Assign[["Microglia"]]=c("MC01","MC02","MC07")
Assign[["Macs"]]=c("MC04","MC05")
Assign[["APCs"]]=c("MC03","MC06","MC08")
Assign[["Proliferating"]]=c("MC09")

CellInfo$Assignment3=NA
for(i in 1:length(Assign)){
  CellInfo$Assignment3[CellInfo$Cluster %in% Assign[[i]]]=names(Assign[i])
}
SeuratObj@meta.data=CellInfo
Assignment3Colors= Nour_pal("all")(length(levels(as.factor(SeuratObj@meta.data$Assignment3))))
names(Assignment3Colors)= levels(as.factor(SeuratObj@meta.data$Assignment3))

AssignmentColors= Nour_pal("all")(length(levels(as.factor(SeuratObj@meta.data$Assignment))))
names(AssignmentColors)= levels(as.factor(SeuratObj@meta.data$Assignment))
U3=DimPlot(SeuratObj, reduction = "umap",group.by = "Assignment" ,cols =  AssignmentColors,label=T)+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="Assignment")+
  guides(colour = guide_legend(override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=10,face = "bold"),ncol = 4))

pdf(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"Assignment UMAP.pdf"),width = 5.5,height =6,family = "ArialMT")
U3
dev.off()
Assignment2Colors= Nour_pal("all")(length(levels(as.factor(SeuratObj@meta.data$Assignment2))))
names(AssignmentColors)= levels(as.factor(SeuratObj@meta.data$Assignment))
U6=DimPlot(SeuratObj, reduction = "umap",group.by = "Assignment2" ,cols =  Assignment2Colors,label=T)+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="Assignment")+
  guides(colour = guide_legend(override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=10,face = "bold"),ncol = 4))

pdf(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"Assignment2 UMAP.pdf"),width = 5.5,height =6,family = "ArialMT")
U6
dev.off()
DimPlot(SeuratObj, reduction = "umap",group.by = "Assignment",split.by = "Cluster" ,cols =  AssignmentColors,label=T)+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="Assignment")+
  guides(colour = guide_legend(override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=10,face = "bold"),ncol = 4))
DimPlot(SeuratObj, reduction = "umap",group.by = "Assignment",split.by = "Patient" ,cols =  AssignmentColors,label=F,ncol = 6)+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="Assignment")+
  guides(colour = guide_legend(override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=10,face = "bold"),ncol = 4))



library(msigdbr)
Idents(SeuratObj)=CellInfo$Cluster
sampled.cells <- sample(x = SeuratObj@active.ident, size = 20000, replace = F)
subset=subset(SeuratObj,cells=names(sampled.cells))

Counts <- as.matrix(GetAssayData(object = subset, assay = "RNA", slot = "counts"))



m_df_H<- msigdbr(species = "Homo sapiens", category = "H")
Filtered.genesets=filter(m_df_H,gs_name=="HALLMARK_MYC_TARGETS_V1"|gs_name=="HALLMARK_TNFA_SIGNALING_VIA_NFKB"|gs_name=="HALLMARK_ALLOGRAFT_REJECTION"|gs_name=="HALLMARK_HYPOXIA")
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
h=hierarchy (sc, quadrants = c("HALLMARK_MYC_TARGETS_V1","HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_ALLOGRAFT_REJECTION","HALLMARK_HYPOXIA"), log.scale =T)
group.cols= ClusterColors
CellInfo=SeuratObj@meta.data
for( i in 1: length(names(column.col2))){
  name=names(column.col2)[i]
  color=column.col2[i]
  CellInfo$Clustercolor [CellInfo$Cluster== name] <- color
}
Clustercolor <- CellInfo$Clustercolor
names(Clustercolor ) <- row.names(CellInfo)
SeuratObj <- AddMetaData(object = SeuratObj,metadata = Clustercolor ,col.name = 'Clustercolor')
#figure 2d
#Make the big plot--> all clusters
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
title="All Myeloid Clusters"

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
title="All Myeloid Clusters"
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
  annotate("text",x = xmin+2, y = ymax-0.5,label = "Allo-rej",color="white",fontface="bold",size=8)+
  annotate("rect", xmin = xmax-4, xmax = xmax, ymin = ymax-1, ymax = ymax, fill= "black")  +
  annotate("text",x = xmax-2, y = ymax-0.5,label = "Hypoxia",color="white",fontface="bold",size=8)+
  annotate("rect", xmin = xmin, xmax = xmin+4, ymin = ymin+1, ymax = ymin, fill= "black")  +
  annotate("text",x = xmin+2, y = ymin+0.5,label = "MYC-tar",color="white",fontface="bold",size=8)+
  annotate("rect", xmin = xmax-4, xmax =xmax, ymin = ymin+1, ymax = ymin, fill= "black")  +
  annotate("text",x = xmax-2, y = ymin+0.5,label = "TNFA-sig",color="white",fontface="bold",size=8)
P0=ggplotGrob(p0)
##make the small plots--> one per Cluster
Final <- list()
Final[[1]]=ggplotGrob(p0)
for (i in 1:length(names(column.col2))) {
  ClusterMD=SeuratObj@meta.data[SeuratObj@meta.data$Cluster==paste0(names(column.col2)[i]),]
  groups=ClusterMD[,c("Clustercolor","Cluster")]
  title=paste0(names(column.col2)[i])
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
    annotate("text",x = xmin+2, y = ymax-0.5,label = "Allo-rej",color="white",fontface="bold",size=4)+
    annotate("rect", xmin = xmax-4, xmax = xmax, ymin = ymax-1, ymax = ymax, fill= "black")  +
    annotate("text",x = xmax-2, y = ymax-0.5,label = "Hypoxia",color="white",fontface="bold",size=4)+
    annotate("rect", xmin = xmin, xmax = xmin+4, ymin = ymin+1, ymax = ymin, fill= "black")  +
    annotate("text",x = xmin+2, y = ymin+0.5,label = "MYC-tar",color="white",fontface="bold",size=4)+
    annotate("rect", xmin = xmax-4, xmax =xmax, ymin = ymin+1, ymax = ymin, fill= "black")  +
    annotate("text",x = xmax-2, y = ymin+0.5,label = "TNFA-sig",color="white",fontface="bold",size=4)
  Final[[i+1]] = ggplotGrob(P)
}
Final[[1]]=ggplotGrob(p0)
numofplots= length(Final)
pdf(file =paste0(OutputDirectory,ObjName,Subset," pathway butterfly by Cluster myc trgets instead of oxphos and allograft rej insted of ifng.pdf"), height = 6, width =19.5,onefile = T)
grid.arrange(grobs=Final, widths = c(2,1,1,1,1,1),layout_matrix = rbind(c(1, 2, 3,4,5,6),c(1,7,8,9,10,NA)))
dev.off()
##
m_df_H<- msigdbr(species = "Homo sapiens", category = "H")
Filtered.genesets=filter(m_df_H,gs_name=="HALLMARK_OXIDATIVE_PHOSPHORYLATION"|gs_name=="HALLMARK_TNFA_SIGNALING_VIA_NFKB"|gs_name=="HALLMARK_INTERFERON_GAMMA_RESPONSE"|gs_name=="HALLMARK_HYPOXIA")
geneset.all <-   Filtered.genesets[Filtered.genesets$gene_symbol %in% rownames(Counts),]
genesetlist=geneset.all%>% split(x = .$gene_symbol, f = .$gs_name)
i=45
for(i in 1:length(genesetlist)){
  genes=genesetlist[[i]]
  plot.data<-as.data.frame(assay(SeuratObj.sce, "logcounts"))
  plot.data<-plot.data[genes,]
  plot.data<-filter(plot.data,rowMeans(plot.data)>0.2)
  plot.data <- plot.data - rowMeans(plot.data)
  plot.data=na.omit(plot.data)
  plot.data<-plot.data[order(rowMeans(plot.data),decreasing=T),]

  column_annot <-CellInfoS[,c("Cluster","Patient","Type"),drop=F]
  column_annot$Patient = as.factor(as.character(column_annot$Patient))
  column_annot=with(column_annot, column_annot[order(Patient), , drop=F])
  column_annot=with(column_annot, column_annot[order(Cluster), , drop=F])
  plot.data<-plot.data[,row.names(column_annot)]
  column.col= SampleColors

  column.col2= ClusterColors
  column.col3= FragmentColors
  column.col4= TypeColors


  column.colors=list()
  column.colors[["Patient"]]<-column.col
  #column.colors[["Fragment"]]<-column.col3
  column.colors[["Cluster"]]<-column.col2
  column.colors[["Type"]]<-column.col4

  Patient=as.matrix(column_annot[,c("Patient"),drop=F])
  Cluster=as.matrix(column_annot[,c("Cluster"),drop=F])
  #Fragment=as.matrix(column_annot[,c("Fragment"),drop=F])
  Type=as.matrix(column_annot[,c("Type"),drop=F])

  col<- circlize::colorRamp2(breaks = c(-2, 0, 2), colors = c("#007dd1", "white", "#ab3000"))
  XX=names(genesetlist[i])
  YY=str_remove(XX, "HALLMARK_")
  HM=Heatmap(name=XX,as.matrix(t(plot.data)),cluster_rows =F,cluster_columns = T,column_title =YY,column_title_gp = gpar(fontsize=18,fontface="bold"),height = 5,width = 7,row_split = Cluster,
             col = col,show_column_names=T,show_row_names=F,show_column_dend = F,row_names_side = "left",column_names_gp =gpar(fontsize=4.5) ,border = "black",show_heatmap_legend = F,use_raster = T)
  assign(paste0("HM",i),HM)
}
colanno <- rowAnnotation(df=column_annot,
                         show_annotation_name =F,show_legend = F,col=column.colors)
lgd=Legend(title = "logcounts", at=  c(-2,0, 2),col_fun = col)
lgd1=Legend(labels = levels(as.factor(column_annot$Cluster)),title="Cluster",legend_gp = gpar(fill=column.col2,fontsize=5))
lgd2=Legend(labels = levels(as.factor(column_annot$Patient)),title="Patient",legend_gp = gpar(fill=column.col,fontsize=5))
#lgd3=Legend(labels = levels(as.factor(column_annot$Fragment)),title="Fragment",legend_gp = gpar(fill=column.col3,fontsize=5))
lgd4=Legend(labels = levels(as.factor(column_annot$Type)),title="Type",legend_gp = gpar(fill=column.col4,fontsize=5))
HM=HM1+HM2+HM3+HM4+colanno
pdf(paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"large Pathway combined heatmap.pdf")),width=25,height=10)
draw(HM,heatmap_legend_list = list(lgd,lgd4, lgd2,lgd1), heatmap_legend_side = "right",legend_labels_gp =gpar(col = "black", fontsize = 20,fontface="bold"))
dev.off()

##

macmarkclassic=read.csv("~/Desktop/paper final figures/Supp Fig 3/macmarkclassic.csv")
geneset.all <-   macmarkclassic[ macmarkclassic$Gene %in% rownames(Counts),]
geneset <- geneset.all[!duplicated(geneset.all$gene),]
genesetlist=geneset.all%>% split(x =.$Gene, f = .$MacType,drop = F)
Scounts=Counts[geneset$Gene,]
#get scores
sc=score(mat=Counts,
         groups=genesetlist,
         binmat = NULL,
         bins = NULL,
         controls = NULL,
         bin.control = F,
         center = F,
         nbin = 30,
         n = 100,
         replace = T)
SC.dat=merge(sc,CellInfo,by.x=0,by.y=0)
#write.csv(SC.dat,"scscoreforclassic mac markers.csv")
colors =ClusterColors
Mactype=levels(as.factor(macmarkclassic$MacType))
Final=list()
for(i in 1:length(Mactype)){
  YY=  Mactype[i]
  XX=ggplot(SC.dat ,aes(Cluster,y=eval(parse(text=YY))))+ geom_boxplot(colour = "black",width=0.5,lwd=0.3)+
    geom_boxplot(colour = "black",fill=colors,width=0.5,lwd=0.3)+
    labs(y ="Metamodule Score", x= NULL,title = YY)+
    theme(plot.title = element_text(face = "bold",color = "black",hjust = 0.5,size=20),
          axis.text.y = element_text(color="black",size=10),
          axis.text.x = element_text(color="black",angle = 25,hjust = 1),
          axis.text=element_text(size=15, face="bold"),
          axis.title=element_text(size=15,face="bold"), plot.margin = unit(c(0.1, 0.5,0.2, 0.5),"cm"))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    stat_summary(fun=median, colour="black", geom="text", show.legend = FALSE,vjust=-0.7, aes( label=round(..y.., digits=2)))
  Final[[i]]=ggplotGrob(XX)
}
pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"classic mac marker scores.pdf")), height = 10, width =20,onefile = T)
grid.arrange(grobs=Final, widths = c(1,1,1),layout_matrix = rbind(c(1, 2, 3),c(4,5,NA)))
dev.off()
###
##to get figure 2e
#Load Dataset
CGGA.EXP <- read.table("~/Downloads/CGGA.mRNAseq_325.RSEM-genes.20200506.txt",as.is=T,header = T,row.names=1)
AnnotationCGGA <- read.csv("~/Downloads/2020-08-20_CGGA_pheno.csv",as.is=T,header = T,row.names=1)
CGGA.GBM=filter(AnnotationCGGA,Histology=="GBM")
#got markers
AllMyeloid.Markers= read.csv("~/Box/Yun lab projects/scRNAseq Project/Human Glioma/ALLSAMPLESJuly2021/Output/Myeloid-no Tcells or Astrocytes/HumanGlioma-AllSamples-Myeloid-no Tcells or Astrocytes cluster markers res0.25.csv")
Myeloid.Markers=AllMyeloid.Markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(Myeloid.Markers,"~/Box/Yun lab manuscripts/GBM Single Cells/Nature Communications revision/RCode/files needed for figure 2e and f/Myeloidmarkers.csv")
#Myeloid.Markers=top20.markers

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
#write.csv(sc,paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"cluster CGGAscores.csv")))
#scores will slightly chage everytime you run them because of random sampling of control genes
#sc=read.csv(paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"assignment CGGAscores.csv")),row.names = 1)
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
  Expression.Level=paste(YY,"Levels",sep=".")
  CGGA.sc2 <- CGGA.sc2%>% mutate(Expression.Level = ifelse(CGGA.sc2[YY]>=0, "Positive", "Negative"))
  CGGA.sc2$Expression.Level <-factor(CGGA.sc2$Expression.Level)
  fit1 <- survfit(surv_object ~Expression.Level, data = CGGA.sc2)
  XX <- ggsurvplot(fit1, data = CGGA.sc2, pval = TRUE,pval.coord = c(50, 1),pval.size=6,legend="top",palette = c("Blue","Red"),risk.table = F,font.main = c(12, "bold"),legend.title = "status",
                   title=YY, legend.labs = c("Not-Enriched", "Enriched"),font.x = c(14, "bold"),font.y = c(14, "bold"),font.tickslab = c(12, "plain"),xlab="Time (Months)")
  Final[[i]] = XX
}
pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"new cluster survival plots CGGA allFigure2E.pdf")), height = 3.5, width =16)
arrange_ggsurvplots(Final,print = TRUE, ncol = 5,nrow = 1,surv.plot.height = NULL,risk.table.height = NULL,ncensor.plot.height = NULL)
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
#write.csv(sc,paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"cluster GBMonlyCGGAscores.csv")))
#Since everytime you run the scoring the results are slightly different I provided the scoring file that reproduces the exact figures in the paper
#sc=read.csv(paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"assignment GBMonlyCGGAscores.csv")),row.names = 1)
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
  Expression.Level=paste(YY,"Levels",sep=".")
  CGGA.sc2 <- CGGA.sc2%>% mutate(Expression.Level = ifelse(CGGA.sc2[YY]>=0, "Positive", "Negative"))
  CGGA.sc2$Expression.Level <-factor(CGGA.sc2$Expression.Level)
  fit1 <- survfit(surv_object ~Expression.Level, data = CGGA.sc2)
  XX <- ggsurvplot(fit1, data = CGGA.sc2, pval = TRUE,pval.coord = c(50, 1),pval.size=6,legend="top",palette = c("Blue","Red"),risk.table = F,font.main = c(12, "bold"),legend.title = "status",
                   title=YY, legend.labs = c("Not-Enriched", "Enriched"),font.x = c(14, "bold"),font.y = c(14, "bold"),font.tickslab = c(12, "plain"),xlab="Time (Months)")
  Final[[i]] = XX
}
pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"new cluster GBM only survivalFigure2F.pdf")), height = 3.5, width =16)
arrange_ggsurvplots(Final,print = TRUE, ncol = 5,nrow = 1,surv.plot.height = NULL,risk.table.height = NULL,ncensor.plot.height = NULL)
dev.off()

#
### to make supplementary figure 4b
#Load Data
TCGA.EXP <- read.table("~/Downloads/Human__TCGA_GBM__UNC__RNAseq__GA_RNA__01_28_2016__BI__Gene__Firehose_RSEM_log2.cct",as.is=T,header = T,row.names=1)
Annotation <- read.table("~/Downloads//Human__TCGA_GBM__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi",as.is=T,header = T,row.names=1)
Annotation =t(Annotation)

Counts=as.data.frame(TCGA.EXP)
geneset.all <-   Myeloid.Markers.surv[ Myeloid.Markers.surv$gene %in% rownames(Counts),]
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
#write.csv(sc,paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"cluster TCGAscores.csv")))

#scores will slightly chage everytime you run them because of random sampling of control genes
scx=sc
#sc=read.csv(paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"cluster TCGAscores.csv")),row.names = 1)
##Survival Analysis
library(survival)
library(survminer)
TCGA.sc=merge(sc,Annotation,by.x=0, by.y=0)
TCGA.sc$OS<-as.numeric(TCGA.sc$overall_survival)
TCGA.sc$status<-as.numeric(TCGA.sc$status)
surv_object <- Surv(time = TCGA.sc$OS, event = TCGA.sc$status)
#for loop to assign "positive" or "negative" to cell scores then draw survival plots
data.subset=colnames(sc)
data.subset=sort(data.subset, decreasing = FALSE)
TCGA.sc2 =TCGA.sc
Final <- list()
for( i in 1: length(data.subset))
{
  YY= data.subset[i]
  Expression.Level=paste(YY,"Levels",sep=".")
  TCGA.sc2 <- TCGA.sc2%>% mutate(Expression.Level = ifelse(TCGA.sc2[YY]>=0, "Positive", "Negative"))
  TCGA.sc2$Expression.Level <-factor(TCGA.sc2$Expression.Level)
  fit1 <- survfit(surv_object ~Expression.Level, data = TCGA.sc2)

  XX <- ggsurvplot(fit1, data = TCGA.sc2, pval = TRUE,pval.coord = c(1000, 1),pval.size=6,legend="top",palette = c("Blue","Red"),risk.table = F,font.main = c(12, "bold"),legend.title = "status",
                   title=YY, legend.labs = c("Not-Enriched", "Enriched"),font.x = c(14, "bold"),font.y = c(14, "bold"),font.tickslab = c(12, "plain"),xlab="Time (Days)")
  Final[[i]] = XX
}
pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"new cluster survival plots TCGA supp4E.pdf")), height = 3.5, width =16)
arrange_ggsurvplots(Final,print = TRUE, ncol = 5,nrow = 1,surv.plot.height = NULL,risk.table.height = NULL,ncensor.plot.height = NULL)
dev.off()

####
TCGA.EXP <- read.table("~/Downloads/Human__TCGA_GBM__UNC__RNAseq__HiSeq_RNA__01_28_2016__BI__Gene__Firehose_RSEM_log2.cct",as.is=T,header = T,row.names=1)
Annotation <- read.table("~/Downloads/2020-08-19_TCGA_GBM_pheno.txt",as.is=T,header = T,row.names=1)

Counts=as.data.frame(TCGA.EXP)
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
#write.csv(sc,paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"cluster TCGA2scores.csv")))


#sc=read.csv("~/Desktop/paper final figures/supp fig 4/TCGA2cellstates.csv",row.names = 1)
##Survival Analysis
library(survival)
library(survminer)
TCGA.sc=merge(sc,Annotation,by.x=0, by.y=0)
TCGA.sc$OS<-as.numeric(TCGA.sc$survival)
TCGA.sc$status<-as.numeric(TCGA.sc$status)
surv_object <- Surv(time = TCGA.sc$OS, event = TCGA.sc$status)
#for loop to assign "positive" or "negative" to cell scores then draw survival plots
data.subset=colnames(sc)
data.subset=sort(data.subset, decreasing = FALSE)
TCGA.sc2 =TCGA.sc
Final <- list()
for( i in 1: length(data.subset))
{
  YY= data.subset[i]
  Expression.Level=paste(YY,"Levels",sep=".")
  TCGA.sc2 <- TCGA.sc2%>% mutate(Expression.Level = ifelse(TCGA.sc2[YY]>=0, "Positive", "Negative"))
  TCGA.sc2$Expression.Level <-factor(TCGA.sc2$Expression.Level)
  fit1 <- survfit(surv_object ~Expression.Level, data = TCGA.sc2)
  XX <- ggsurvplot(fit1, data = TCGA.sc2, pval = TRUE,pval.coord = c(30, 1),pval.size=6,legend="top",palette = c("Blue","Red"),risk.table = F,font.main = c(12, "bold"),legend.title = "status",
                   title=YY, legend.labs = c("Not-Enriched", "Enriched"),font.x = c(14, "bold"),font.y = c(14, "bold"),font.tickslab = c(12, "plain"),xlab="Time (Months)")
  Final[[i]] = XX
}
pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"new cluster survival plots TCGA#2 supp4E.pdf")), height = 3.5, width =16)
arrange_ggsurvplots(Final,print = TRUE, ncol = 5,nrow = 1,surv.plot.height = NULL,risk.table.height = NULL,ncensor.plot.height = NULL)
dev.off()
###
CGGA.EXP <- read.table("~/Downloads/CGGA.mRNAseq_693.RSEM-genes.20200506.txt",as.is=T,header = T,row.names=1)
#AnnotationCGGA <- read.csv("~/Downloads/2020-08-20_CGGA_pheno.csv",as.is=T,header = T,row.names=1)
CGGA.GBM=filter(AnnotationCGGA,Histology=="GBM")

#remove genes from signature gene lists that do not exist in expression matrix ##otherwise the code won't work
Counts=as.data.frame(CGGA.EXP)
geneset.all <-   Myeloid.Markers.surv[ Myeloid.Markers.surv$gene %in% rownames(Counts),]
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
write.csv(sc,"CGGA#2scores.csv")#scores will slightly chage everytime you run them because of random sampling of control genes
sc=read.csv("CGGA#2scores.csv",row.names = 1)
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
rnames=c("Expression.LevelPositive", "GenderFemale ", "GenderMale",
        "GenderMale ", "HistologyAnaplastic Oligoastrocytoma", "HistologyAnaplastic Oligodendrolgioma",
        "HistologyAstrocytoma", "HistologyGBM", "HistologyOligoastrocytoma",
        "HistologyOligodendroglioma", "RecurrenceRecurrent", "MGMT_statusun-methylated",
        "IDH_mutation_statusWildtype")
Summtable=as.data.frame(row.names =rnames,x=rnames  )

for( i in 1: length(data.subset))
{
  YY= data.subset[i]
  Expression.Level=paste(YY,"Levels",sep=".")
  CGGA.sc2 <- CGGA.sc2%>% mutate(Expression.Level = ifelse(CGGA.sc2[YY]>=0, "Positive", "Negative"))
  CGGA.sc2$Expression.Level <-factor(CGGA.sc2$Expression.Level)
  fit1 <- survfit(surv_object ~Expression.Level, data = CGGA.sc2)
  #next three lines sdded by Joshy
  fit1.mv <- coxph(surv_object  ~Expression.Level + Gender + Histology + Recurrence  + MGMT_status +IDH_mutation_status , data = CGGA.sc2)
  mv.res <- summary(fit1.mv)$coefficients
  mv.fname <- paste(OutputDirectory,ObjName,Subset,"res",resolution,data.subset[i],"_CGGA#2-multivariat_results.csv",sep="")
  write.csv(mv.res,mv.fname)
  mv.res=as.data.frame(mv.res)
  Summtable[i]=mv.res$`Pr(>|z|)`
  colnames(Summtable)[i]=YY

  XX <- ggsurvplot(fit1, data = CGGA.sc2, pval = TRUE,pval.coord = c(50, 1),pval.size=6,legend="top",palette = c("Blue","Red"),risk.table = F,font.main = c(12, "bold"),legend.title = "status",
                   title=YY, legend.labs = c("Not-Enriched", "Enriched"),font.x = c(14, "bold"),font.y = c(14, "bold"),font.tickslab = c(12, "plain"),xlab="Time (Months)")
  Final[[i]] = XX
}
write.csv(Summtable,paste(OutputDirectory,ObjName,Subset,"res",resolution,"CGGA2 summary_multivariat_results.csv",sep=""))

pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"new cluster survival plots GCGA#2 ALL GLIOMASsupp4E.pdf")), height = 3.5, width =16)
arrange_ggsurvplots(Final,print = TRUE, ncol = 5,nrow = 1,surv.plot.height = NULL,risk.table.height = NULL,ncensor.plot.height = NULL)
dev.off()
##Figure4c
#Myeloid.Markers= FindAllMarkers(SeuratObj,only.pos = F)
Clusters=levels(as.factor(CellInfo$Cluster))
for (i in 1: length(Clusters)){
  X=  Clusters[i]
  Genes<-SeuratObjClusters.genes %>%
    dplyr::filter(group == X) %>%
    dplyr::filter(padj <= 0.05) %>%
    dplyr::filter(logFC <= -0.1|logFC >= 0.1) %>%
    arrange(desc(auc)) %>%
    dplyr::select(feature)
  Y=gost(query = Genes$feature,organism = "hsapiens",ordered_query = T,multi_query = F,evcodes =T,
         significant = T,correction_method = "fdr" )
  Z=Y[["result"]]
  names(Z)[10] <- "Category"
  names(Z)[3] <- "adj_pval"
  names(Z)[16] <- "Genes"
  names(Z)[9] <- "ID"
  names(Z)[11] <- "Term"

  assign(X, Z)
}
library("gprofiler2")

library(stringr)
library(dplyr)
library(GOplot)
Myeloid.Markers=SeuratObjClusters.genes
names(Myeloid.Markers)[1] <- "ID"
names(Myeloid.Markers)[8] <- "adj.P.Val"
names(Myeloid.Markers)[2] <- "cluster"

##switch to directory containing gprofiler files... this file is generated using the g:profiler website (https://biit.cs.ut.ee/gprofiler/gost) by uploading the DE genes for each cluster
#filess <- list.files(path = "~/Downloads/",pattern="-gProfiler.csv")

for (file in filess){
  X=read.csv(paste0("~/Downloads/",file))
  colnames(X)=c("Category",	"Term",	"ID",	"adj_pval",	"negative_log10_of_adjusted_p_value",	"term_size",	"query_size",
                "intersection_size",	"effective_domain_size",	"Genes")
  NAME= str_replace(paste(file),"-gProfiler.csv","")
  assign(NAME, X)
}

Final=list()
for (i in 1:9){
  MC=paste0("MC0",i)
  print(MC)
  DEgenes=Myeloid.Markers%>%dplyr::filter(cluster==MC)
  df=eval(parse(text=MC))
  names(df)[11] <- "Term"
  names(df)[9] <- "ID"
  circ<-circle_dat(df,DEgenes)
  process = c("canonical glycolysis","cellular response to hypoxia","Antigen processing-Cross presentation","cellular response to interferon-gamma","cellular response to tumor necrosis factor","NF-kappa B signaling pathway","oxidative phosphorylation","M Phase")
  #processx=process[process%in%circ$term]
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
  NAME=paste0("Graphs.MC0",i)
  assign(NAME, P1)
  print(paste(MC,"finished"))
}
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(Graphs.MC01)
lay <- rbind(c(1,2,3),
             c(4,5,6),
             c(7,8,9),
             c(10,10,10))
mylegend$vp$angle=0
png(paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"Figure4c.png")),width=7,height=6.5,units = "in",res = 1200,type = "quartz",pointsize = 10)
grid.arrange(Graphs.MC01+ theme(legend.position="none"),Graphs.MC02+ theme(legend.position="none"),
             Graphs.MC03+ theme(legend.position="none"),Graphs.MC04+ theme(legend.position="none"),
             Graphs.MC05+ theme(legend.position="none"),Graphs.MC06+ theme(legend.position="none"),
             Graphs.MC07+ theme(legend.position="none"),Graphs.MC08+ theme(legend.position="none"),
             Graphs.MC09+ theme(legend.position="none"),mylegend, ncol = 3, nrow = 4,layout_matrix=lay)
dev.off()
Idents(SeuratObj)=SeuratObj$Cluster
for (i in 1:length(c(Genes))){
  Y=Genes[i]
  X=FeaturePlot(SeuratObj,features=Y,order = T,label = T,label.size = 2,cols = Nour_cols(c("darkpurple","lightorange")))+theme_min()
  pdf(paste0(featuremapDirectory,Y," featuremap ",ObjName," ",Subset,".pdf"),width = 4.5,height = 4)
  print(X)
  dev.off()

}
##
pat2=c("LGG-01", "LGG-04", "ndGBM-01", "ndGBM-02", "ndGBM-03", "rGBM-01", "rGBM-02", "rGBM-03", "rGBM-04", "rGBM-05")

for(i in 1:length(pat2)){
  x=ggplot(SeuratObj@meta.data%>%filter(Patient==pat2[i]), aes(Patient, fill=Assignment))+geom_bar(stat="count",position="fill",colour = "white",width =0.7)+
    scale_fill_manual(values = AssignmentColors)+coord_polar("y",start=0) + facet_grid(~ Fragment) + #labs(title=pat2[i])+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_blank(),axis.ticks = element_blank(),
          panel.background = element_blank(), axis.title = element_blank(),strip.background=element_rect(fill="white"),
          strip.text=element_text(face = "bold",color = "black",size=8),plot.margin = unit(c(0,-0.2,-0.2,0),'cm'))
  Final[[i]]=x
}
mylegend<-get_legend(Final[[2]]+theme_min()+ NoAxes()+guides(colour = guide_colourbar(title=Y,title.position = "top",title.theme = element_text(size=10))))

pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"piecharts fragment and assignment.pdf")), height = 6, width =10,onefile = T)
ggpubr::ggarrange(Final[[1]]+NoLegend(),Final[[6]]+NoLegend(),Final[[2]]+NoLegend(),Final[[7]]+NoLegend(),Final[[3]]+NoLegend(),
                  Final[[8]]+NoLegend(),Final[[4]]+NoLegend(),Final[[9]]+NoLegend(),Final[[5]]+NoLegend(),Final[[10]],ncol = 2,nrow = 5,
                  legend.grob = mylegend,legend = "right",widths = c(1,1))
dev.off()
SeuratObj@meta.data$Type2=factor(SeuratObj@meta.data$Type, levels=c("Astrocytoma", "Oligodendroglioma", "GBM", "Recurrent GBM"))
pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"piecharts fragment and assignment.pdf")), height = 6, width =10,onefile = T)

x1=ggplot(SeuratObj@meta.data, aes("", fill=Assignment))+geom_bar(stat="count",position="fill",colour = "white",width =0.7)+
  scale_fill_manual(values = AssignmentColors)+coord_polar("y",start=0) + facet_grid(~ Type2) + #labs(title=pat2[i])+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_blank(),axis.ticks = element_blank(),
        panel.background = element_blank(), axis.title = element_blank(),strip.background=element_rect(fill="white"),
        strip.text=element_text(face = "bold",color = "black",size=8),plot.margin = unit(c(0,-0.2,-0.2,0),'cm'))
x2=ggplot(SeuratObj@meta.data, aes("", fill=Assignment))+geom_bar(stat="count",position="fill",colour = "white",width =0.7)+
  scale_fill_manual(values = AssignmentColors)+coord_polar("y",start=0) + facet_grid(~ sex) + #labs(title=pat2[i])+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_blank(),axis.ticks = element_blank(),
        panel.background = element_blank(), axis.title = element_blank(),strip.background=element_rect(fill="white"),
        strip.text=element_text(face = "bold",color = "black",size=8),plot.margin = unit(c(0,-0.2,-0.2,0),'cm'))
x3=ggplot(SeuratObj@meta.data, aes("", fill=Assignment))+geom_bar(stat="count",position="fill",colour = "white",width =0.7)+
  scale_fill_manual(values = AssignmentColors)+coord_polar("y",start=0) + facet_grid(~ Grade) + #labs(title=pat2[i])+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_blank(),axis.ticks = element_blank(),
        panel.background = element_blank(), axis.title = element_blank(),strip.background=element_rect(fill="white"),
        strip.text=element_text(face = "bold",color = "black",size=8),plot.margin = unit(c(0,-0.2,-0.2,0),'cm'))

mylegend<-get_legend(Final[[2]]+theme_min()+ NoAxes()+guides(colour = guide_colourbar(title=Y,title.position = "top",title.theme = element_text(size=10))))

pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"piecharts type,sexandgrade and assignment.pdf")), height = 4, width =6,onefile = T)
ggpubr::ggarrange(x1+NoLegend(),x2+NoLegend(),x3+NoLegend(),ncol =1,nrow = 3,
                  legend.grob = mylegend,legend = "right")
dev.off()
pat2=c("LGG-01", "LGG-04", "ndGBM-01", "ndGBM-02", "ndGBM-03", "rGBM-01", "rGBM-02", "rGBM-03", "rGBM-04", "rGBM-05")

for(i in 1:length(pat2)){
  x=ggplot(SeuratObj@meta.data%>%filter(Patient==pat2[i]), aes(Patient, fill=Assignment3))+geom_bar(stat="count",position="fill",colour = "white",width =0.7)+
    scale_fill_manual(values = Assignment3Colors)+coord_polar("y",start=0) + facet_grid(~ Fragment) + #labs(title=pat2[i])+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_blank(),axis.ticks = element_blank(),
          panel.background = element_blank(), axis.title = element_blank(),strip.background=element_rect(fill="white"),
          strip.text=element_text(face = "bold",color = "black",size=8),plot.margin = unit(c(0,-0.2,-0.2,0),'cm'))
  Final[[i]]=x
}
mylegend<-get_legend(Final[[2]]+theme_min()+ NoAxes()+guides(colour = guide_colourbar(title=Y,title.position = "top",title.theme = element_text(size=10))))

pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"piecharts fragment and Assignment3.pdf")), height = 6, width =10,onefile = T)
ggpubr::ggarrange(Final[[1]]+NoLegend(),Final[[6]]+NoLegend(),Final[[2]]+NoLegend(),Final[[7]]+NoLegend(),Final[[3]]+NoLegend(),
                  Final[[8]]+NoLegend(),Final[[4]]+NoLegend(),Final[[9]]+NoLegend(),Final[[5]]+NoLegend(),Final[[10]],ncol = 2,nrow = 5,
                  legend.grob = mylegend,legend = "right",widths = c(1,1))
dev.off()

pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"piecharts patient and assignment.pdf")), height = 3, width =6,onefile = T)

ggplot(SeuratObj@meta.data, aes("", fill=Assignment))+geom_bar(stat="count",position="fill",colour = "white",width =0.7)+
  scale_fill_manual(values = AssignmentColors)+coord_polar("y",start=0) + facet_wrap(~ Patient,ncol = 6) + #labs(title=pat2[i])+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_blank(),axis.ticks = element_blank(),
        panel.background = element_blank(), axis.title = element_blank(),strip.background=element_rect(fill="white"),
        strip.text=element_text(face = "bold",color = "black",size=8),plot.margin = unit(c(0,-0.2,-0.2,0),'cm'))
dev.off()

pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"piecharts patient and Cluster.pdf")), height = 3, width =6,onefile = T)

ggplot(SeuratObj@meta.data, aes("", fill=Cluster))+geom_bar(stat="count",position="fill",colour = "white",width =0.7)+
  scale_fill_manual(values = ClusterColors)+coord_polar("y",start=0) + facet_wrap(~ Patient,ncol = 6) + #labs(title=pat2[i])+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_blank(),axis.ticks = element_blank(),
        panel.background = element_blank(), axis.title = element_blank(),strip.background=element_rect(fill="white"),
        strip.text=element_text(face = "bold",color = "black",size=8),plot.margin = unit(c(0,-0.2,-0.2,0),'cm'))
dev.off()

pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"piecharts fragment and Assignment3.pdf")), height = 6, width =10,onefile = T)

x1=ggplot(SeuratObj@meta.data, aes("", fill=Assignment3))+geom_bar(stat="count",position="fill",colour = "white",width =0.7)+
  scale_fill_manual(values = Assignment3Colors)+coord_polar("y",start=0) + facet_grid(~ Type2) + #labs(title=pat2[i])+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_blank(),axis.ticks = element_blank(),
        panel.background = element_blank(), axis.title = element_blank(),strip.background=element_rect(fill="white"),
        strip.text=element_text(face = "bold",color = "black",size=8),plot.margin = unit(c(0,-0.2,-0.2,0),'cm'))
x2=ggplot(SeuratObj@meta.data, aes("", fill=Assignment3))+geom_bar(stat="count",position="fill",colour = "white",width =0.7)+
  scale_fill_manual(values = Assignment3Colors)+coord_polar("y",start=0) + facet_grid(~ sex) + #labs(title=pat2[i])+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_blank(),axis.ticks = element_blank(),
        panel.background = element_blank(), axis.title = element_blank(),strip.background=element_rect(fill="white"),
        strip.text=element_text(face = "bold",color = "black",size=8),plot.margin = unit(c(0,-0.2,-0.2,0),'cm'))
x3=ggplot(SeuratObj@meta.data, aes("", fill=Assignment3))+geom_bar(stat="count",position="fill",colour = "white",width =0.7)+
  scale_fill_manual(values = Assignment3Colors)+coord_polar("y",start=0) + facet_grid(~ Grade) + #labs(title=pat2[i])+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_blank(),axis.ticks = element_blank(),
        panel.background = element_blank(), axis.title = element_blank(),strip.background=element_rect(fill="white"),
        strip.text=element_text(face = "bold",color = "black",size=8),plot.margin = unit(c(0,-0.2,-0.2,0),'cm'))

mylegend<-get_legend(Final[[2]]+theme_min()+ NoAxes()+guides(colour = guide_colourbar(title=Y,title.position = "top",title.theme = element_text(size=10))))

pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"piecharts type,sexandgrade and Assignment3.pdf")), height = 4, width =6,onefile = T)
ggpubr::ggarrange(x1+NoLegend(),x2+NoLegend(),x3+NoLegend(),ncol =1,nrow = 3,
                  legend.grob = mylegend,legend = "right")
dev.off()
saveRDS(SeuratObj,file = paste0(RobjDirectory,"MyeloidClusters-7-20-21-patients renamed.rds"))
SeuratObj=readRDS(paste0(RobjDirectory,"MyeloidClusters-7-20-21-patients renamed.rds"))
VlnPlot(SeuratObj,features = "S100A4",cols=ClusterColors,pt.size = 0)+theme_min2()+NoLegend()+theme(title = element_blank())
