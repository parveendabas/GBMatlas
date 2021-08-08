library(easypackages)
MyPackages<-c("dplyr","ggplot2","ggpubr","gridExtra","viridis","egg",
              "grid","lattice","gtools","Biobase","RColorBrewer","tibble",
              "Seurat","cowplot","patchwork","stringr","ComplexHeatmap","SingleCellExperiment", "ggmin")
libraries(MyPackages)

#Set working directory
#Directory="~/Box/Yun lab projects/scRNAseq Project/Human Glioma/ALLSAMPLESJuly2021//"
Directory="/projects/ebf-lab/PK/Kyuson/GBM/Analysis_GBM_Batch1_2_3_ALL/FINAL_Robjects_from_Nour"
RobjDirectory=paste0(Directory,"R_Objects/")
#dir.create(RobjDirectory)

setwd(Directory)

#Load File
SeuratObj=readRDS("~/Box/Kyuson_GBM/2021_GBM_All_Three_Batches_Analysis_USE_THIS/Denovo_DataSets/Glioma/GBM_ALL_Batches_Correction_subset_Glioma_minGenes_500_PCA_30_Res_0.09_theta_0.5_0.5.rds")
DimPlot(SeuratObj)

ObjName= "HumanGlioma-AllSamples-"
Subset="Glioma"
resolution= 0.2
dimz=1:30
abb="GC"

OutputDirectory=paste0(Directory,"Output/",Subset,"/")
#dir.create(OutputDirectory)
featuremapDirectory=paste0(OutputDirectory,"Featuremaps/")
#dir.create(featuremapDirectory)

CellInfo=SeuratObj@meta.data
Patients=levels(as.factor(SeuratObj@meta.data$Sample))
Patientnum=readr::parse_number(Patients)
Patientnum=str_remove(Patientnum,"-")
Patientnum=as.numeric(Patientnum)
CellInfo$Patient=CellInfo$Sample
for(j in 1:length(Patients)){
  if(Patientnum[j]<10){
    CellInfo$Patient[CellInfo$Sample==Patients[j]]=str_replace(CellInfo$Sample[CellInfo$Sample==Patients[j]],"-","-0")
    CellInfo$Fragment[CellInfo$Sample==Patients[j]]=str_replace(CellInfo$Fragment[CellInfo$Sample==Patients[j]],"MDAG-","MDAG-0")
  }else{
    CellInfo$Fragment[CellInfo$Sample==Patients[j]]=str_replace(CellInfo$Fragment[CellInfo$Sample==Patients[j]],"MDAG-","MDAG-")

  }
}
levels(as.factor(CellInfo$Patient))
levels(as.factor(CellInfo$Fragment))
Assign= list()
Assign[["GBM"]]=c("CNSTM-070","CNSTM-096","MDAG-04",   "MDAG-06",   "MDAG-07",   "MDAG-09",   "MDAG-10", "MDAG-12")
Assign[["Astrocytoma"]]=c("CNSTM-081","MDAG-01","MDAG-11")
Assign[["Recurrent GBM"]]=c("CNSTM-068","CNSTM-375", "CNSTM-379", "CNSTM-390", "CNSTM-394", "CNSTM-397")
Assign[["Oligodendroglioma"]]=c("MDAG-03")

CellInfo$Type=NA
for(i in 1:length(Assign)){
  CellInfo$Type[CellInfo$Patient %in% Assign[[i]]]=names(Assign[i])
}
levels(as.factor(CellInfo$Type))
Assign= list()
Assign[["IV"]]=c("CNSTM-070","CNSTM-096","MDAG-04",   "MDAG-06",   "MDAG-07",   "MDAG-09",   "MDAG-10", "MDAG-12","CNSTM-068","CNSTM-375", "CNSTM-379", "CNSTM-390", "CNSTM-394", "CNSTM-397")
Assign[["II"]]=c("CNSTM-081","MDAG-11","MDAG-03")
Assign[["III"]]=c("MDAG-01")

CellInfo$Grade=NA
for(i in 1:length(Assign)){
  CellInfo$Grade[CellInfo$Patient %in% Assign[[i]]]=names(Assign[i])
}
levels(as.factor(CellInfo$Grade))
SeuratObj@meta.data=CellInfo
CellInfo=SeuratObj@meta.data
Assign= list()
Assign[["Hiseq"]]=c("CNSTM-070","CNSTM-096","CNSTM-081","CNSTM-068")
Assign[["Novoseq"]]=c("CNSTM-375", "CNSTM-379", "CNSTM-390", "CNSTM-394", "CNSTM-397")
Assign[["Novoseq2"]]=c("MDAG-01","MDAG-03","MDAG-04",   "MDAG-06",   "MDAG-07",   "MDAG-09",   "MDAG-10", "MDAG-12","MDAG-11")

CellInfo$Platform=NA
for(i in 1:length(Assign)){
  CellInfo$Platform[CellInfo$Patient %in% Assign[[i]]]=names(Assign[i])
}

Assign= list()
Assign[["Batch1"]]=c("CNSTM-068")
Assign[["Batch1"]]=c("CNSTM-070")
Assign[["Batch1"]]=c("CNSTM-081")
Assign[["Batch1"]]=c("CNSTM-096")
Assign[["Batch5"]]=c("CNSTM-375", "CNSTM-379", "CNSTM-390", "CNSTM-394", "CNSTM-397")
Assign[["Batch6"]]=c("MDAG-01","MDAG-03","MDAG-04",   "MDAG-06",   "MDAG-07",   "MDAG-09",   "MDAG-10", "MDAG-12","MDAG-11")

CellInfo$Batch=NA
for(i in 1:length(Assign)){
  CellInfo$Batch[CellInfo$Patient %in% Assign[[i]]]=names(Assign[i])
}

SeuratObj@meta.data=CellInfo
levels(as.factor(CellInfo$Platform))
DimPlot(SeuratObj,group.by = "Platform")
#Type color
TypeColors= Nour_pal("main",reverse = F)(length(levels(as.factor(SeuratObj@meta.data$Type))))
names(TypeColors)=levels(as.factor(SeuratObj@meta.data$Type))

#Grade color
GradeColors= Nour_pal("main",reverse = T)(length(levels(as.factor(SeuratObj@meta.data$Grade))))
names(GradeColors)=levels(as.factor(SeuratObj@meta.data$Grade))##
##
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

SeuratObj <- CellCycleScoring(SeuratObj, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
DimPlot(SeuratObj,group.by = "Phase")
SeuratObj <- ScaleData(SeuratObj, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(SeuratObj))

#initial clustering
library(harmony)
SeuratObj <- RunHarmony(SeuratObj, c("Patient","sex","Platform"))
ElbowPlot(SeuratObj,ndims = 25,reduction = "harmony")

SeuratObj <- FindNeighbors(SeuratObj,reduction = "harmony", dims = 1:19)

SeuratObj <- RunUMAP(SeuratObj, reduction = "harmony",dims = 1:19 )
DimPlot(SeuratObj)
#SeuratObj <- FindNeighbors(SeuratObj, dims = dimz)
resolution=0.2
SeuratObj <- FindClusters(SeuratObj, resolution = 0.2)
DimPlot(SeuratObj,label = T)
DimPlot(SeuratObj,split.by = "Grade",ncol = 6,label = T)

#head(Idents(SeuratObj), 5)
#SeuratObj <- RunUMAP(SeuratObj, dims = dimz)

DimPlot(SeuratObj)
abb="GC"
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
set.seed(12)
ClusterColors= sample(Nour_pal("cool",reverse = T)(length(levels(as.factor(SeuratObj@meta.data$Cluster)))))
names(ClusterColors)=levels(as.factor(SeuratObj@meta.data$Cluster))

# Get number of cells per cluster and per Sample
write.csv(as.matrix(table(SeuratObj@meta.data$Cluster,SeuratObj@meta.data$Sample)),
          file = paste0(OutputDirectory,ObjName,Subset," number of cells per cluster and sample.csv"))

pat2=c("LGG-01", "LGG-04", "ndGBM-01", "ndGBM-02", "ndGBM-03", "rGBM-01", "rGBM-02", "rGBM-03", "rGBM-04", "rGBM-05")

for(i in 1:length(pat2)){
  x=ggplot(SeuratObj@meta.data%>%filter(Patient==pat2[i]), aes(Patient, fill=Cluster))+geom_bar(stat="count",position="fill",colour = "white",width =0.7)+
    scale_fill_manual(values = ClusterColors)+coord_polar("y",start=0) + facet_grid(~ Fragment) + #labs(title=pat2[i])+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_blank(),axis.ticks = element_blank(),
          panel.background = element_blank(), axis.title = element_blank(),strip.background=element_rect(fill="white"),
          strip.text=element_text(face = "bold",color = "black",size=8),plot.margin = unit(c(0,-0.2,-0.2,0),'cm'))
  Final[[i]]=x
}
mylegend<-get_legend(Final[[2]]+theme_min()+ NoAxes()+guides(colour = guide_colourbar(title=Y,title.position = "top",title.theme = element_text(size=10))))

pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"piecharts fragment and Cluster.pdf")), height = 6, width =10,onefile = T)
ggpubr::ggarrange(Final[[1]]+NoLegend(),Final[[6]]+NoLegend(),Final[[2]]+NoLegend(),Final[[7]]+NoLegend(),Final[[3]]+NoLegend(),
                  Final[[8]]+NoLegend(),Final[[4]]+NoLegend(),Final[[9]]+NoLegend(),Final[[5]]+NoLegend(),Final[[10]],ncol = 2,nrow = 5,
                  legend.grob = mylegend,legend = "right",widths = c(1,1))
dev.off()

SeuratObj@meta.data$Type2=factor(SeuratObj@meta.data$Type, levels=c("Astrocytoma", "Oligodendroglioma", "GBM", "Recurrent GBM"))
pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"piecharts fragment and Cluster.pdf")), height = 6, width =10,onefile = T)

x1=ggplot(SeuratObj@meta.data, aes("", fill=Cluster))+geom_bar(stat="count",position="fill",colour = "white",width =0.7)+
  scale_fill_manual(values = ClusterColors)+coord_polar("y",start=0) + facet_grid(~ Type2) + #labs(title=pat2[i])+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_blank(),axis.ticks = element_blank(),
        panel.background = element_blank(), axis.title = element_blank(),strip.background=element_rect(fill="white"),
        strip.text=element_text(face = "bold",color = "black",size=8),plot.margin = unit(c(0,-0.2,-0.2,0),'cm'))
x2=ggplot(SeuratObj@meta.data, aes("", fill=Cluster))+geom_bar(stat="count",position="fill",colour = "white",width =0.7)+
  scale_fill_manual(values = ClusterColors)+coord_polar("y",start=0) + facet_grid(~ sex) + #labs(title=pat2[i])+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_blank(),axis.ticks = element_blank(),
        panel.background = element_blank(), axis.title = element_blank(),strip.background=element_rect(fill="white"),
        strip.text=element_text(face = "bold",color = "black",size=8),plot.margin = unit(c(0,-0.2,-0.2,0),'cm'))
x3=ggplot(SeuratObj@meta.data, aes("", fill=Cluster))+geom_bar(stat="count",position="fill",colour = "white",width =0.7)+
  scale_fill_manual(values = ClusterColors)+coord_polar("y",start=0) + facet_grid(~ Grade) + #labs(title=pat2[i])+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_blank(),axis.ticks = element_blank(),
        panel.background = element_blank(), axis.title = element_blank(),strip.background=element_rect(fill="white"),
        strip.text=element_text(face = "bold",color = "black",size=8),plot.margin = unit(c(0,-0.2,-0.2,0),'cm'))

mylegend<-get_legend(Final[[2]]+theme_min()+ NoAxes()+guides(colour = guide_colourbar(title=Y,title.position = "top",title.theme = element_text(size=10))))

pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"piecharts type,sexandgrade and Cluster.pdf")), height = 4, width =6,onefile = T)
ggpubr::ggarrange(x1+NoLegend(),x2+NoLegend(),x3+NoLegend(),ncol =1,nrow = 3,
                  legend.grob = mylegend,legend = "right")
dev.off()







CellInfo$Patient=CellInfo$Sample
SeuratObj@meta.data=CellInfo
SeuratObj <- NormalizeData(SeuratObj, normalization.method = "LogNormalize", scale.factor = 10000)
SeuratObj <- FindVariableFeatures(SeuratObj, selection.method = "vst", nfeatures = 2000)
SeuratObj <- ScaleData(SeuratObj, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(SeuratObj))

#Run PCA
SeuratObj <- RunPCA(SeuratObj, features = VariableFeatures(object = SeuratObj))


#initial clustering
library(harmony)
SeuratObj <- RunHarmony(SeuratObj, c("Patient","sex"),theta=c(2,0.5))
SeuratObj <- FindNeighbors(SeuratObj,reduction = "harmony", dims = 1:15)

SeuratObj <- RunUMAP(SeuratObj,reduction = "harmony", dims = 1:15 )
DimPlot(SeuratObj)

#SeuratObj <- FindNeighbors(SeuratObj, dims = dimz)

SeuratObj <- FindClusters(SeuratObj, resolution = 0.3)
DimPlot(SeuratObj,split.by = "Patient",ncol = 6,label = T)

#head(Idents(SeuratObj), 5)
#SeuratObj <- RunUMAP(SeuratObj, dims = dimz)

DimPlot(SeuratObj)

#supplementary fig 5a
CellInfo=SeuratObj@meta.data

#Rename Clusters
clusters=levels(as.factor(SeuratObj@meta.data$RNA_snn_res.0.3))
for(j in 1:length(clusters)){
  if(j<10){
    CellInfo$Cluster[CellInfo$RNA_snn_res.0.3==j-1]=paste0(abb,"0",j)
  }else{
    CellInfo$Cluster[CellInfo$RNA_snn_res.0.3==j-1]=paste0(abb,j)
  }
}
SeuratObj@meta.data=CellInfo

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



#Rename Clusters
#CellInfo$Cluster=str_replace_all(CellInfo$Cluster,"C","GC0")
#SeuratObj@meta.data=CellInfo


# Add Sample color info
# Sample color
SampleColors= Nour_pal("all")(length(levels(as.factor(SeuratObj$Patient))))
names(SampleColors)=levels(as.factor(SeuratObj$Patient))

# Cluster color
set.seed(12)
ClusterColors= sample(Nour_pal("cool",reverse = T)(length(levels(as.factor(SeuratObj@meta.data$Cluster)))))
names(ClusterColors)=levels(as.factor(SeuratObj@meta.data$Cluster))

#Fragment color
CellInfo=SeuratObj@meta.data
CellInfo$Fragment[which(str_detect(row.names(CellInfo), "57"))] <- "CNSTM-394-3"
SeuratObj@meta.data=CellInfo

FragmentColors= Nour_pal("all",reverse = F)(length(levels(as.factor(SeuratObj@meta.data$Fragment))))
names(FragmentColors)=levels(as.factor(SeuratObj@meta.data$Fragment))
#
CellInfo=SeuratObj@meta.data
Assign= list()
Assign[["GBM"]]=c("CNSTM-070","CNSTM-096","MDAG-04",   "MDAG-06",   "MDAG-07",   "MDAG-09",   "MDAG-10", "MDAG-12")
Assign[["Astrocytoma"]]=c("CNSTM-081","MDAG-01","MDAG-11")
Assign[["Recurrent GBM"]]=c("CNSTM-068","CNSTM-375", "CNSTM-379", "CNSTM-390", "CNSTM-394", "CNSTM-397")
Assign[["Oligodendroglioma"]]=c("MDAG-03")

CellInfo$Type=NA
for(i in 1:length(Assign)){
  CellInfo$Type[CellInfo$Patient %in% Assign[[i]]]=names(Assign[i])
}

Assign= list()
Assign[["IV"]]=c("CNSTM-070","CNSTM-096","MDAG-04",   "MDAG-06",   "MDAG-07",   "MDAG-09",   "MDAG-10", "MDAG-12","CNSTM-068","CNSTM-375", "CNSTM-379", "CNSTM-390", "CNSTM-394", "CNSTM-397")
Assign[["II"]]=c("CNSTM-081","MDAG-11","MDAG-03")
Assign[["III"]]=c("MDAG-01")

CellInfo$Grade=NA
for(i in 1:length(Assign)){
  CellInfo$Grade[CellInfo$Patient %in% Assign[[i]]]=names(Assign[i])
}
SeuratObj@meta.data=CellInfo


#Type color
TypeColors= Nour_pal("main",reverse = F)(length(levels(as.factor(SeuratObj@meta.data$Type))))
names(TypeColors)=levels(as.factor(SeuratObj@meta.data$Type))

#Grade color
GradeColors= Nour_pal("main",reverse = T)(length(levels(as.factor(SeuratObj@meta.data$Grade))))
names(GradeColors)=levels(as.factor(SeuratObj@meta.data$Grade))


#make heatmap

Idents(SeuratObj)=CellInfo$Cluster
markers = FindAllMarkers(SeuratObj,logfc.threshold = 0.5,test.use = "wilcox",only.pos = T )
write.csv(markers,paste0(OutputDirectory,ObjName,Subset," cluster markers ","res",resolution,".csv"))
#markers = read.csv(paste0(OutputDirectory,ObjName,Subset," cluster markers ","res",resolution,".csv"))
markers=markers[order(markers$cluster),]
top20.markers <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top20.markers$cluster=as.character(top20.markers$cluster)
top20.markers=top20.markers[order(top20.markers$cluster),]
set.seed(12)
Idents(SeuratObj)=CellInfo$Sample
sampled.cells <- sample(x = SeuratObj@active.ident, size = 5000, replace = F)
subset=subset(SeuratObj,cells=names(sampled.cells))
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
column.col= SampleColors

column.col2= ClusterColors
#column.col3= FragmentColors
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

colanno <- columnAnnotation(df=column_annot,
                            show_annotation_name =T,show_legend = F,col=column.colors)
top2 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
genes= top2$gene
genes=c("APOD", "PDGFRA", "SPP1", "CHI3L1","CD74", "VEGFA","GFAP", "ATP5E", "FAM64A", "NMU","OLIG2","OLIG1","GFAP", "NES", "PROM1", "VIM", "ERBB2", "TGFBI", "SOX2", "ERBB3", "EGFR", "SOX9", "S100A4", "PTPRZ"," MBP", "B2M", "HLA-DRA", "CCL3")

rows=rowAnnotation(sel = anno_mark(at = match(genes,row.names(plot.data)), labels = genes,labels_gp =gpar(col = "black", fontsize = 9)))

col<- circlize::colorRamp2(breaks = c(-2, 0, 2), colors = c("#007dd1", "white", "#ab3000"))
HM=Heatmap(name="logcounts",as.matrix(plot.data),cluster_rows =F,cluster_columns = F,top_annotation = colanno, right_annotation = rows,row_names_gp = gpar(fontsize=5),
           col = col,show_column_names= F,show_row_names=F,border = F,show_heatmap_legend = F,use_raster = T)
lgd=Legend(title = "logcounts", at=  c(-2,0, 2),col_fun = col)
lgd1=Legend(labels = levels(as.factor(column_annot$Cluster)),title="Cluster",legend_gp = gpar(fill=column.col2,fontsize=5))
lgd2=Legend(labels = levels(as.factor(column_annot$Patient)),title="Patient",legend_gp = gpar(fill=column.col,fontsize=5))
#lgd3=Legend(labels = levels(as.factor(column_annot$Fragment)),title="Fragment",legend_gp = gpar(fill=column.col3,fontsize=5))
lgd4=Legend(labels = levels(as.factor(column_annot$Type)),title="Type",legend_gp = gpar(fill=column.col4,fontsize=5))

pdf(paste0(OutputDirectory,"heatmap",ObjName,Subset,"res",resolution,"top20 genes per cluster nofrag-namedit.pdf"),width=6,height=5)
draw(HM,heatmap_legend_list = list( lgd,lgd1,lgd4, lgd2), heatmap_legend_side = "right")
dev.off()

##Umaps fig3b
U1=DimPlot(SeuratObj, reduction = "umap",group.by = "Cluster" ,cols =  column.col2,label=T)+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="Clusters")+
  guides(colour = guide_legend(override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=15),ncol = 4))

pdf(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"Clusters UMAP.pdf"),width = 5.5,height =6,family = "ArialMT")
U1
dev.off()


U2=DimPlot(SeuratObj, reduction = "umap",group.by = "Patient",cols = SampleColors,shuffle = T )+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="Patients")+
  guides(colour = guide_legend(override.aes = list(size=4),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=10),ncol = 3))

pdf(paste0(OutputDirectory,ObjName,Subset,"Samples UMAP.pdf"),width = 5.5,height =6,family = "ArialMT")
U2
dev.off()

U3=DimPlot(SeuratObj, reduction = "umap",group.by = "Type" ,cols =  TypeColors,label=F)+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="Type")+
  guides(colour = guide_legend(override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=10,face = "bold"),ncol = 2))

pdf(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"Type UMAP.pdf"),width = 5.5,height =6,family = "ArialMT")
U3
dev.off()
U4=DimPlot(SeuratObj, reduction = "umap",group.by = "Grade" ,cols =  GradeColors,label=T)+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="Grade")+
  guides(colour = guide_legend(override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=10,face = "bold"),ncol = 3))

pdf(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"Grade UMAP.pdf"),width = 5.5,height =6,family = "ArialMT")
U4
dev.off()

plots2=ggarrange(U1,U2,U3,U4,ncol = 2)

pdf(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"all UMAPs.pdf"),width = 8,height = 9,family = "ArialMT")
print(plots2)
dev.off()




pdf(paste0(OutputDirectory,ObjName,Subset,"Cluster UMAP Iteration by sample.pdf"),width = 12,height =8,family = "ArialMT")
DimPlot(SeuratObj, reduction = "umap",group.by = "Cluster",cols =column.col2,split.by = "Sample",ncol = 3, pt.size = 0.2)+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 16, y.title = 16,main = 16)+ theme(legend.position = "right",text =element_text(size=15) )+
  guides(colour = guide_legend(title="Clusters", override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=15)))
dev.off()

#figure3c
MyPackages<-c("ggplot2","ggpubr","gridExtra","grid","Seurat","scrabble")
libraries(MyPackages)
Idents(SeuratObj)=CellInfo$Sample
sampled.cells <- sample(x = SeuratObj@active.ident, size = 20000, replace = F)
subset=subset(SeuratObj,cells=names(sampled.cells))

Counts <- as.matrix(GetAssayData(object = subset, assay = "RNA", slot = "counts"))



group.cols= ClusterColors

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
##Perpatient
CellInfo=SeuratObj@meta.data
for( i in 1: length(names(column.col))){
  name=names(column.col)[i]
  color=column.col[i]
  CellInfo$Patientcolor [CellInfo$Patient== name] <- color
}
Patientcolor <- CellInfo$Patientcolor
names(Patientcolor ) <- row.names(CellInfo)
SeuratObj <- AddMetaData(object = SeuratObj,metadata = Patientcolor ,col.name = 'Patientcolor')

# make plots
xmax=max(h$X)+(max(h$X)*0.2)
xmin=min(h$X)+(min(h$X)*0.2)
ymax=max(h$Y)+(max(h$Y)*0.6)
ymin=min(h$Y)+(min(h$Y)*0.6)

#Make the big plot--> all Patients
groups=SeuratObj@meta.data[,c("Patientcolor","Patient")]
matrix=merge(h,groups,by.x=0,by.y=0,all.x = T)
matrix$Patientcolor[is.na(matrix$Patientcolor)] <- "gray"
matrix$Patient[is.na(matrix$Patient)] <- "Other"
row.names(matrix)=matrix$Row.names
x=matrix$Patientcolor
y=matrix$Patient
col=x[!duplicated(x)]
names(col)=y[!duplicated(y)]
matrix=matrix[,-1]
title="All Patients"
p0 <- ggplot(matrix, aes(x = X,
                         y =Y,color=factor(Patient)))+geom_point()+geom_point(data = subset(matrix, Patient !="Other"))+
  scale_color_manual(values=col,aesthetics = "colour", breaks = waiver()) +labs(x=NULL,y=NULL,colour="Patient")+theme(legend.position = "none")+
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
  annotate("text",x = xmin+2.5, y = ymax-1,label = "Mes-Like",color="white",fontface="bold",size=9)+
  annotate("rect", xmin = xmax-5, xmax = xmax, ymin = ymax-2, ymax = ymax, fill= "black")  +
  annotate("text",x = xmax-2.5, y = ymax-1,label = "NPC-Like",color="white",fontface="bold",size=9)+
  annotate("rect", xmin = xmin, xmax = xmin+5, ymin = ymin+2, ymax = ymin, fill= "black")  +
  annotate("text",x = xmin+2.5, y = ymin+1,label = "AC-Like",color="white",fontface="bold",size=9)+
  annotate("rect", xmin = xmax-5, xmax =xmax, ymin = ymin+2, ymax = ymin, fill= "black")  +
  annotate("text",x = xmax-2.5, y = ymin+1,label = "OPC-Like",color="white",fontface="bold",size=9)
P0=ggplotGrob(p0)
##make the small plots--> one per Patient
Final <- list()
Final[[1]]=ggplotGrob(p0)
for (i in 1:length(names(column.col))) {
  PatientMD=SeuratObj@meta.data[SeuratObj@meta.data$Patient==paste0(names(column.col)[i]),]
  groups=PatientMD[,c("Patientcolor","Patient")]
  title=paste0(names(column.col)[i])
  matrix=merge(h,groups,by.x=0,by.y=0,all.x = T)
  matrix$Patientcolor[is.na(matrix$Patientcolor)] <- "gray"
  matrix$Patient=as.character(matrix$Patient)
  matrix$Patient[is.na(matrix$Patient)] <- "Other"
  row.names(matrix)=matrix$Row.names
  x=matrix$Patientcolor
  y=matrix$Patient
  col=x[!duplicated(x)]
  names(col)=y[!duplicated(y)]
  matrix=matrix[,-1]
  P=ggplot(matrix, aes(x = X,
                       y =Y,color=factor(Patient)))+geom_point()+geom_point(data = subset(matrix, Patient !="Other"))+
    scale_color_manual(values=col,aesthetics = "colour", breaks = waiver()) +labs(x=NULL,y=NULL)+theme(legend.position = "none")+
    scale_x_continuous(expand = c(0, 0), limits = c(xmin,xmax)) + scale_y_continuous(expand = c(0, 0), limits = c(ymin,ymax))+
    theme(panel.background = element_rect(fill = "white",colour = "white"),axis.ticks.x=element_blank(),axis.text.x=element_blank(),
          axis.ticks.y=element_blank(), axis.text.y=element_blank())+
    ggtitle(title)+theme(plot.title =element_text(size=25,face="bold") )+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+
    geom_hline(yintercept=0, color = "black", size=0.5)+
    geom_vline(xintercept=0, color = "black", size=0.5)+
    annotate("rect", xmin = xmin, xmax = xmin+2, ymin = ymax-1, ymax = ymax, fill= "black")  +
    annotate("rect", xmin = xmin, xmax = xmin+5, ymin = ymax-2, ymax = ymax, fill= "black")  +
    annotate("text",x = xmin+2.5, y = ymax-1,label = "Mes-Like",color="white",fontface="bold",size=4)+
    annotate("rect", xmin = xmax-5, xmax = xmax, ymin = ymax-2, ymax = ymax, fill= "black")  +
    annotate("text",x = xmax-2.5, y = ymax-1,label = "NPC-Like",color="white",fontface="bold",size=4)+
    annotate("rect", xmin = xmin, xmax = xmin+5, ymin = ymin+2, ymax = ymin, fill= "black")  +
    annotate("text",x = xmin+2.5, y = ymin+1,label = "AC-Like",color="white",fontface="bold",size=4)+
    annotate("rect", xmin = xmax-5, xmax =xmax, ymin = ymin+2, ymax = ymin, fill= "black")  +
    annotate("text",x = xmax-2.5, y = ymin+1,label = "OPC-Like",color="white",fontface="bold",size=4)
  Final[[i+1]] = ggplotGrob(P)
}
Final[[1]]=ggplotGrob(p0)
numofplots= length(Final)
pdf(file =paste0(OutputDirectory,ObjName,Subset," subtype butterfly by Patient ordered by type2.pdf"), height = 10, width =20,onefile = T)
grid.arrange(grobs=Final, widths = c(3,1,1,1,1,1),layout_matrix = rbind(c(1, 2, 3,4,NA,5),c(1,6,7,8,9,10),
                                                                        c(1,11,12,13,14,NA),
                                                                        c(1,15,16,17,18,19)))
dev.off()
##per cluster
CellInfo=SeuratObj@meta.data
for( i in 1: length(names(column.col2))){
  name=names(column.col2)[i]
  color=column.col2[i]
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
  ggtitle(title)+theme(plot.title =element_text(size=22,face="bold") )+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+
  geom_hline(yintercept=0, color = "black", size=0.5)+
  geom_vline(xintercept=0, color = "black", size=0.5)+
  annotate("rect", xmin = xmin, xmax = xmin+5, ymin = ymax-2, ymax = ymax, fill= "black")  +
  annotate("text",x = xmin+2.5, y = ymax-1,label = "Mes-Like",color="white",fontface="bold",size=7)+
  annotate("rect", xmin = xmax-5, xmax = xmax, ymin = ymax-2, ymax = ymax, fill= "black")  +
  annotate("text",x = xmax-2.5, y = ymax-1,label = "NPC-Like",color="white",fontface="bold",size=7)+
  annotate("rect", xmin = xmin, xmax = xmin+5, ymin = ymin+2, ymax = ymin, fill= "black")  +
  annotate("text",x = xmin+2.5, y = ymin+1,label = "AC-Like",color="white",fontface="bold",size=7)+
  annotate("rect", xmin = xmax-5, xmax =xmax, ymin = ymin+2, ymax = ymin, fill= "black")  +
  annotate("text",x = xmax-2.5, y = ymin+1,label = "OPC-Like",color="white",fontface="bold",size=7)
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
    ggtitle(title)+theme(plot.title =element_text(size=22,face="bold") )+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+
    geom_hline(yintercept=0, color = "black", size=0.5)+
    geom_vline(xintercept=0, color = "black", size=0.5)+
    annotate("rect", xmin = xmin, xmax = xmin+2, ymin = ymax-1, ymax = ymax, fill= "black")  +
    annotate("rect", xmin = xmin, xmax = xmin+5, ymin = ymax-2, ymax = ymax, fill= "black")  +
    annotate("text",x = xmin+2.5, y = ymax-1,label = "Mes-Like",color="white",fontface="bold",size=4)+
    annotate("rect", xmin = xmax-5, xmax = xmax, ymin = ymax-2, ymax = ymax, fill= "black")  +
    annotate("text",x = xmax-2.5, y = ymax-1,label = "NPC-Like",color="white",fontface="bold",size=4)+
    annotate("rect", xmin = xmin, xmax = xmin+5, ymin = ymin+2, ymax = ymin, fill= "black")  +
    annotate("text",x = xmin+2.5, y = ymin+1,label = "AC-Like",color="white",fontface="bold",size=4)+
    annotate("rect", xmin = xmax-5, xmax =xmax, ymin = ymin+2, ymax = ymin, fill= "black")  +
    annotate("text",x = xmax-2.5, y = ymin+1,label = "OPC-Like",color="white",fontface="bold",size=4)
  Final[[i+1]] = ggplotGrob(P)
}
Final[[1]]=ggplotGrob(p0)
numofplots= length(Final)
pdf(file =paste0(OutputDirectory,ObjName,Subset," subtype butterfly by Cluster realign for figs.pdf"), height = 6, width =18,onefile = T)
grid.arrange(grobs=Final, widths = c(2,1,1,1,1,1),layout_matrix = rbind(c(1, 2, 3,4,5,6),c(1,7,8,9,10,NA)))
dev.off()


##perfragment

CellInfo=SeuratObj@meta.data
for( i in 1: length(names(column.col3))){
  name=names(column.col3)[i]
  color=column.col3[i]
  CellInfo$Fragmentcolor [CellInfo$Fragment== name] <- color
}
Fragmentcolor <- CellInfo$Fragmentcolor
names(Fragmentcolor ) <- row.names(CellInfo)
SeuratObj <- AddMetaData(object = SeuratObj,metadata = Fragmentcolor ,col.name = 'Fragmentcolor')

# make plots
xmax=max(h$X)+(max(h$X)*0.2)
xmin=min(h$X)+(min(h$X)*0.2)
ymax=max(h$Y)+(max(h$Y)*0.6)
ymin=min(h$Y)+(min(h$Y)*0.6)

#Make the big plot--> all Fragments
groups=SeuratObj@meta.data[,c("Fragmentcolor","Fragment")]
matrix=merge(h,groups,by.x=0,by.y=0,all.x = T)
matrix$Fragmentcolor[is.na(matrix$Fragmentcolor)] <- "gray"
matrix$Fragment[is.na(matrix$Fragment)] <- "Other"
row.names(matrix)=matrix$Row.names
x=matrix$Fragmentcolor
y=matrix$Fragment
col=x[!duplicated(x)]
names(col)=y[!duplicated(y)]
matrix=matrix[,-1]
title="All Glioma Fragments"
p0 <- ggplot(matrix, aes(x = X,
                         y =Y,color=factor(Fragment)))+geom_point()+geom_point(data = subset(matrix, Fragment !="Other"))+
  scale_color_manual(values=col,aesthetics = "colour", breaks = waiver()) +labs(x=NULL,y=NULL,colour="Fragment")+theme(legend.position = "none")+
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
  annotate("text",x = xmin+2.5, y = ymax-1,label = "Mes-Like",color="white",fontface="bold",size=9)+
  annotate("rect", xmin = xmax-5, xmax = xmax, ymin = ymax-2, ymax = ymax, fill= "black")  +
  annotate("text",x = xmax-2.5, y = ymax-1,label = "NPC-Like",color="white",fontface="bold",size=9)+
  annotate("rect", xmin = xmin, xmax = xmin+5, ymin = ymin+2, ymax = ymin, fill= "black")  +
  annotate("text",x = xmin+2.5, y = ymin+1,label = "AC-Like",color="white",fontface="bold",size=9)+
  annotate("rect", xmin = xmax-5, xmax =xmax, ymin = ymin+2, ymax = ymin, fill= "black")  +
  annotate("text",x = xmax-2.5, y = ymin+1,label = "OPC-Like",color="white",fontface="bold",size=9)
P0=ggplotGrob(p0)
##make the small plots--> one per Fragment
Final <- list()
Final[[1]]=ggplotGrob(p0)
for (i in 1:length(names(column.col3))) {
  FragmentMD=SeuratObj@meta.data[SeuratObj@meta.data$Fragment==paste0(names(column.col3)[i]),]
  groups=FragmentMD[,c("Fragmentcolor","Fragment")]
  title=paste0(names(column.col3)[i])
  matrix=merge(h,groups,by.x=0,by.y=0,all.x = T)
  matrix$Fragmentcolor[is.na(matrix$Fragmentcolor)] <- "gray"
  matrix$Fragment=as.character(matrix$Fragment)
  matrix$Fragment[is.na(matrix$Fragment)] <- "Other"
  row.names(matrix)=matrix$Row.names
  x=matrix$Fragmentcolor
  y=matrix$Fragment
  col=x[!duplicated(x)]
  names(col)=y[!duplicated(y)]
  matrix=matrix[,-1]
  P=ggplot(matrix, aes(x = X,
                       y =Y,color=factor(Fragment)))+geom_point()+geom_point(data = subset(matrix, Fragment !="Other"))+
    scale_color_manual(values=col,aesthetics = "colour", breaks = waiver()) +labs(x=NULL,y=NULL)+theme(legend.position = "none")+
    scale_x_continuous(expand = c(0, 0), limits = c(xmin,xmax)) + scale_y_continuous(expand = c(0, 0), limits = c(ymin,ymax))+
    theme(panel.background = element_rect(fill = "white",colour = "white"),axis.ticks.x=element_blank(),axis.text.x=element_blank(),
          axis.ticks.y=element_blank(), axis.text.y=element_blank())+
    ggtitle(title)+theme(plot.title =element_text(size=25,face="bold") )+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+
    geom_hline(yintercept=0, color = "black", size=0.5)+
    geom_vline(xintercept=0, color = "black", size=0.5)+
    annotate("rect", xmin = xmin, xmax = xmin+2, ymin = ymax-1, ymax = ymax, fill= "black")  +
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
pdf(file =paste0(OutputDirectory,ObjName,Subset," subtype butterfly by Fragment.pdf"), height = 15, width =35,onefile = T)
grid.arrange(grobs=Final, widths = c(4,1,1,1,1,1,1,1,1,1),layout_matrix = rbind(c(1, 2, 3,4,5,6,7,8,9,10),c(1, 11, 12,13,14,15,16,17,18,19),
                                                                                c(1, 20, 21,22,23,24,25,26,27,NA),c(1, 28,29, 30,31,32,33,34,35,36),
                                                                                c(1, 37,38, 39,40,41,42,43,44,45)))
dev.off()
##per type



CellInfo=SeuratObj@meta.data
for( i in 1: length(names(column.col4))){
  name=names(column.col4)[i]
  color=column.col4[i]
  CellInfo$Typecolor [CellInfo$Type== name] <- color
}
Typecolor <- CellInfo$Typecolor
names(Typecolor ) <- row.names(CellInfo)
SeuratObj <- AddMetaData(object = SeuratObj,metadata = Typecolor ,col.name = 'Typecolor')

# make plots
xmax=max(h$X)+(max(h$X)*0.2)
xmin=min(h$X)+(min(h$X)*0.2)
ymax=max(h$Y)+(max(h$Y)*0.6)
ymin=min(h$Y)+(min(h$Y)*0.6)

#Make the big plot--> all Types
groups=SeuratObj@meta.data[,c("Typecolor","Type")]
matrix=merge(h,groups,by.x=0,by.y=0,all.x = T)
matrix$Typecolor[is.na(matrix$Typecolor)] <- "gray"
matrix$Type[is.na(matrix$Type)] <- "Other"
row.names(matrix)=matrix$Row.names
x=matrix$Typecolor
y=matrix$Type
col=x[!duplicated(x)]
names(col)=y[!duplicated(y)]
matrix=matrix[,-1]
title="All Glioma Types"
p0 <- ggplot(matrix, aes(x = X,
                         y =Y,color=factor(Type)))+geom_point()+geom_point(data = subset(matrix, Type !="Other"))+
  scale_color_manual(values=col,aesthetics = "colour", breaks = waiver()) +labs(x=NULL,y=NULL,colour="Type")+theme(legend.position = "none")+
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
  annotate("text",x = xmin+2.5, y = ymax-1,label = "Mes-Like",color="white",fontface="bold",size=6)+
  annotate("rect", xmin = xmax-5, xmax = xmax, ymin = ymax-2, ymax = ymax, fill= "black")  +
  annotate("text",x = xmax-2.5, y = ymax-1,label = "NPC-Like",color="white",fontface="bold",size=6)+
  annotate("rect", xmin = xmin, xmax = xmin+5, ymin = ymin+2, ymax = ymin, fill= "black")  +
  annotate("text",x = xmin+2.5, y = ymin+1,label = "AC-Like",color="white",fontface="bold",size=6)+
  annotate("rect", xmin = xmax-5, xmax =xmax, ymin = ymin+2, ymax = ymin, fill= "black")  +
  annotate("text",x = xmax-2.5, y = ymin+1,label = "OPC-Like",color="white",fontface="bold",size=6)
P0=ggplotGrob(p0)
##make the small plots--> one per Type
Final <- list()
Final[[1]]=ggplotGrob(p0)
for (i in 1:length(names(column.col4))) {
  TypeMD=SeuratObj@meta.data[SeuratObj@meta.data$Type==paste0(names(column.col4)[i]),]
  groups=TypeMD[,c("Typecolor","Type")]
  title=paste0(names(column.col4)[i])
  matrix=merge(h,groups,by.x=0,by.y=0,all.x = T)
  matrix$Typecolor[is.na(matrix$Typecolor)] <- "gray"
  matrix$Type=as.character(matrix$Type)
  matrix$Type[is.na(matrix$Type)] <- "Other"
  row.names(matrix)=matrix$Row.names
  x=matrix$Typecolor
  y=matrix$Type
  col=x[!duplicated(x)]
  names(col)=y[!duplicated(y)]
  matrix=matrix[,-1]
  P=ggplot(matrix, aes(x = X,
                       y =Y,color=factor(Type)))+geom_point()+geom_point(data = subset(matrix, Type !="Other"))+
    scale_color_manual(values=col,aesthetics = "colour", breaks = waiver()) +labs(x=NULL,y=NULL)+theme(legend.position = "none")+
    scale_x_continuous(expand = c(0, 0), limits = c(xmin,xmax)) + scale_y_continuous(expand = c(0, 0), limits = c(ymin,ymax))+
    theme(panel.background = element_rect(fill = "white",colour = "white"),axis.ticks.x=element_blank(),axis.text.x=element_blank(),
          axis.ticks.y=element_blank(), axis.text.y=element_blank())+
    ggtitle(title)+theme(plot.title =element_text(size=25,face="bold") )+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+
    geom_hline(yintercept=0, color = "black", size=0.5)+
    geom_vline(xintercept=0, color = "black", size=0.5)+
    annotate("rect", xmin = xmin, xmax = xmin+2, ymin = ymax-1, ymax = ymax, fill= "black")  +
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
pdf(file =paste0(OutputDirectory,ObjName,Subset," subtype butterfly by Type.pdf"), height = 6, width =15,onefile = T)
grid.arrange(grobs=Final, widths = c(2,1,1),layout_matrix = rbind(c(1, 2, 3),c(1, 4, 5)))
dev.off()
###
library(msigdbr)
CellInfo=SeuratObj@meta.data
Idents(SeuratObj)=CellInfo$Cluster
#sampled.cells <- sample(x = SeuratObj@active.ident, size = 20000, replace = F)
#subset=subset(SeuratObj,cells=names(sampled.cells))
subset=subset(SeuratObj, downsample=5000)

Counts <- as.matrix(GetAssayData(object = subset, assay = "RNA", slot = "counts"))



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
title="All Glioma Clusters"

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
numofplots= length(Final)
pdf(file =paste0(OutputDirectory,ObjName,Subset," pathway butterfly by Cluster, IFNg, kras, myc targ and EMT.pdf"), height = 9, width =12,onefile = T)
grid.arrange(grobs=Final, widths = c(2,1,1,1),layout_matrix = rbind(c(1, 2, 3,4),c(1,5,6,7),c(1,8,9,10)))
dev.off()
pdf(file =paste0(OutputDirectory,ObjName,Subset," pathway butterfly by Cluster realign for figs.pdf"), height = 6, width =18,onefile = T)
grid.arrange(grobs=Final, widths = c(2,1,1,1,1,1),layout_matrix = rbind(c(1, 2, 3,4,5,6),c(1,7,8,9,10,NA)))
dev.off()


scpath=sc
#scsuva=sc
mergedsc=cbind(scpath,scsuva)
corrsc=cor(mergedsc)
corrsctrunk=corrsc[1:6,7:10]
row.names(corrsctrunk)=c("EMT","Hypoxia","IFNg-res","MYC-tar","Ox-Phos","TNFa-sig")
pheatmap(corrsctrunk)
col<- circlize::colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))

ComplexHeatmap::Heatmap(name="meta-module corr",as.matrix(corrsctrunk),cluster_rows =T,cluster_columns = T,
                        row_names_gp = gpar(fontsize=10),column_names_gp = gpar(fontsize=10),
    show_column_names= T,show_row_names=T,border = F,show_heatmap_legend = T,use_raster = F,col = col)
p.mat <- corrplot::cor.mtest(mergedsc)
p.mattrunk=p.mat$p
p.mattrunk=p.mattrunk[1:6,7:10]
corrsctrunk

corrplot::corrplot(corrsctrunk,type = "full",method = "pie",p.mat = p.mattrunk,insig = "label_sig",sig.level = 0.05,
                   col = colorRampPalette(c("green","black","red"))(5),tl.col = "black",	cl.ratio= 0.5,pch.col	="black"

)

FragmentColors=c(`LGG-01-A` = "#003F5C", `LGG-01-B` = "#005B72", `LGG-01-C` = "#007789", `LGG-01-D` = "#087B91", `LGG-02` = "#1B6486", `LGG-03` = "#2D4C7C", `LGG-04-1` = "#2A5D9C", `LGG-04-2` = "#2572BF", `LGG-04-3` = "#2D76C9", `ndGBM-01-A` = "#4764AE", `ndGBM-01-C` = "#625394", `ndGBM-01-D` = "#795CA3", `ndGBM-01-F` = "#8F68B9", `ndGBM-02-1` = "#9E6CC0", `ndGBM-02-2` = "#9F5FAC", `ndGBM-02-4` = "#9F5398", `ndGBM-02-5` = "#B460A3", `ndGBM-03-1` = "#CE74B4", `ndGBM-03-2` = "#DE7CBA", `ndGBM-03-3` = "#D968A3", `ndGBM-04` = "#D5548C", `ndGBM-05` = "#DD5B87", `ndGBM-06` = "#EB6A88", `ndGBM-07` = "#EF7082", `ndGBM-08` = "#D95A64", `ndGBM-09` = "#C34446", `rGBM-01-A` = "#CD464A", `rGBM-01-B` = "#E5525B", `rGBM-01-C` = "#F75F64", `rGBM-01-D` = "#ED6A43", `rGBM-02-2` = "#E37623", `rGBM-02-3` = "#E77B23", `rGBM-02-4` = "#F47B34", `rGBM-02-5` = "#FF7D3F", `rGBM-03-1` = "#FF8E25", `rGBM-03-2` = "#FF9F0A", `rGBM-03-3` = "#FCB013", `rGBM-04-1` = "#F8C136", `rGBM-04-2` = "#F0D056", `rGBM-04-3` = "#AABF5C", `rGBM-04-4` = "#63AE62", `rGBM-05-1` = "#51A85A", `rGBM-05-2` = "#6FAC46", `rGBM-05-3` = "#8DB032")

CellInfo=SeuratObj@meta.data
for( i in 1: length(names(FragmentColors))){
  name=names(FragmentColors)[i]
  color=FragmentColors[i]
  CellInfo$Fragmentcolor [CellInfo$Fragment== name] <- color
}
Fragmentcolor <- CellInfo$Fragmentcolor
names(Fragmentcolor ) <- row.names(CellInfo)
SeuratObj <- AddMetaData(object = SeuratObj,metadata = Fragmentcolor ,col.name = 'Fragmentcolor')

# make plots
xmax=max(h$X)+(max(h$X)*0.2)
xmin=min(h$X)+(min(h$X)*0.2)
ymax=max(h$Y)+(max(h$Y)*0.6)
ymin=min(h$Y)+(min(h$Y)*0.6)

#Make the big plot--> all Fragments
groups=SeuratObj@meta.data[,c("Fragmentcolor","Fragment")]
matrix=merge(h,groups,by.x=0,by.y=0,all.x = T)
matrix$Fragmentcolor[is.na(matrix$Fragmentcolor)] <- "gray"
matrix$Fragment[is.na(matrix$Fragment)] <- "Other"
row.names(matrix)=matrix$Row.names
x=matrix$Fragmentcolor
y=matrix$Fragment
col=x[!duplicated(x)]
names(col)=y[!duplicated(y)]
matrix=matrix[,-1]
title="All Glioma Fragments"
p0 <- ggplot(matrix, aes(x = X,
                         y =Y,color=factor(Fragment)))+geom_point()+geom_point(data = subset(matrix, Fragment !="Other"))+
  scale_color_manual(values=col,aesthetics = "colour", breaks = waiver()) +labs(x=NULL,y=NULL,colour="Fragment")+theme(legend.position = "none")+
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
  annotate("text",x = xmin+2, y = ymax-0.5,label = "EMT",color="white",fontface="bold",size=8)+
  annotate("rect", xmin = xmax-4, xmax = xmax, ymin = ymax-1, ymax = ymax, fill= "black")  +
  annotate("text",x = xmax-2, y = ymax-0.5,label = "MYC-Tar",color="white",fontface="bold",size=8)+
  annotate("rect", xmin = xmin, xmax = xmin+4, ymin = ymin+1, ymax = ymin, fill= "black")  +
  annotate("text",x = xmin+2, y = ymin+0.5,label = "Hypoxia",color="white",fontface="bold",size=8)+
  annotate("rect", xmin = xmax-4, xmax =xmax, ymin = ymin+1, ymax = ymin, fill= "black")  +
  annotate("text",x = xmax-2, y = ymin+0.5,label = "IFNg-res",color="white",fontface="bold",size=8)
P0=ggplotGrob(p0)
##make the small plots--> one per Fragment
Final <- list()
Final[[1]]=ggplotGrob(p0)
for (i in 1:length(names(FragmentColors))) {
  FragmentMD=SeuratObj@meta.data[SeuratObj@meta.data$Fragment==paste0(names(FragmentColors)[i]),]
  groups=FragmentMD[,c("Fragmentcolor","Fragment")]
  title=paste0(names(FragmentColors)[i])
  matrix=merge(h,groups,by.x=0,by.y=0,all.x = T)
  matrix$Fragmentcolor[is.na(matrix$Fragmentcolor)] <- "gray"
  matrix$Fragment=as.character(matrix$Fragment)
  matrix$Fragment[is.na(matrix$Fragment)] <- "Other"
  row.names(matrix)=matrix$Row.names
  x=matrix$Fragmentcolor
  y=matrix$Fragment
  col=x[!duplicated(x)]
  names(col)=y[!duplicated(y)]
  matrix=matrix[,-1]
  P=ggplot(matrix, aes(x = X,
                       y =Y,color=factor(Fragment)))+geom_point()+geom_point(data = subset(matrix, Fragment !="Other"))+
    scale_color_manual(values=col,aesthetics = "colour", breaks = waiver()) +labs(x=NULL,y=NULL)+theme(legend.position = "none")+
    scale_x_continuous(expand = c(0, 0), limits = c(xmin,xmax)) + scale_y_continuous(expand = c(0, 0), limits = c(ymin,ymax))+
    theme(panel.background = element_rect(fill = "white",colour = "white"),axis.ticks.x=element_blank(),axis.text.x=element_blank(),
          axis.ticks.y=element_blank(), axis.text.y=element_blank())+
    ggtitle(title)+theme(plot.title =element_text(size=25,face="bold") )+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+
    geom_hline(yintercept=0, color = "black", size=0.5)+
    geom_vline(xintercept=0, color = "black", size=0.5)+
    annotate("rect", xmin = xmin, xmax = xmin+4, ymin = ymax-1, ymax = ymax, fill= "black")  +
    annotate("text",x = xmin+2, y = ymax-0.5,label = "EMT",color="white",fontface="bold",size=4)+
    annotate("rect", xmin = xmax-4, xmax = xmax, ymin = ymax-1, ymax = ymax, fill= "black")  +
    annotate("text",x = xmax-2, y = ymax-0.5,label = "MYC-Tar",color="white",fontface="bold",size=4)+
    annotate("rect", xmin = xmin, xmax = xmin+4, ymin = ymin+1, ymax = ymin, fill= "black")  +
    annotate("text",x = xmin+2, y = ymin+0.5,label = "Hypoxia",color="white",fontface="bold",size=4)+
    annotate("rect", xmin = xmax-4, xmax =xmax, ymin = ymin+1, ymax = ymin, fill= "black")  +
    annotate("text",x = xmax-2, y = ymin+0.5,label = "IFNg-res",color="white",fontface="bold",size=4)
  Final[[i+1]] = ggplotGrob(P)
}
Final[[1]]=ggplotGrob(p0)
numofplots= length(Final)
pdf(file =paste0(OutputDirectory,ObjName,Subset," pathway butterfly by Fragment.pdf"), height = 15, width =35,onefile = T)
grid.arrange(grobs=Final, widths = c(4,1,1,1,1,1,1,1,1,1),layout_matrix = rbind(c(1, 2, 3,4,5,6,7,8,9,10),c(1, 11, 12,13,14,15,16,17,18,19),
                                                                                c(1, 20, 21,22,23,24,25,26,27,NA),c(1, 28,29, 30,31,32,33,34,35,36),
                                                                                c(1, 37,38, 39,40,41,42,43,44,45)))
dev.off()
for( i in 1: length(names(column.col))){
  name=names(column.col)[i]
  color=column.col[i]
  CellInfo$Patientcolor [CellInfo$Patient== name] <- color
}
Patientcolor <- CellInfo$Patientcolor
names(Patientcolor ) <- row.names(CellInfo)
SeuratObj <- AddMetaData(object = SeuratObj,metadata = Patientcolor ,col.name = 'Patientcolor')

# make plots
xmax=max(h$X)+(max(h$X)*0.2)
xmin=min(h$X)+(min(h$X)*0.2)
ymax=max(h$Y)+(max(h$Y)*0.6)
ymin=min(h$Y)+(min(h$Y)*0.6)

#Make the big plot--> all Patients
groups=SeuratObj@meta.data[,c("Patientcolor","Patient")]
matrix=merge(h,groups,by.x=0,by.y=0,all.x = T)
matrix$Patientcolor[is.na(matrix$Patientcolor)] <- "gray"
matrix$Patient[is.na(matrix$Patient)] <- "Other"
row.names(matrix)=matrix$Row.names
x=matrix$Patientcolor
y=matrix$Patient
col=x[!duplicated(x)]
names(col)=y[!duplicated(y)]
matrix=matrix[,-1]
title="All Patients"
p0 <- ggplot(matrix, aes(x = X,
                         y =Y,color=factor(Patient)))+geom_point()+geom_point(data = subset(matrix, Patient !="Other"))+
  scale_color_manual(values=col,aesthetics = "colour", breaks = waiver()) +labs(x=NULL,y=NULL,colour="Patient")+theme(legend.position = "none")+
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
  annotate("text",x = xmin+2, y = ymax-0.5,label = "EMT",color="white",fontface="bold",size=9)+
  annotate("rect", xmin = xmax-4, xmax = xmax, ymin = ymax-1, ymax = ymax, fill= "black")  +
  annotate("text",x = xmax-2, y = ymax-0.5,label = "MYC-Tar",color="white",fontface="bold",size=9)+
  annotate("rect", xmin = xmin, xmax = xmin+4, ymin = ymin+1, ymax = ymin, fill= "black")  +
  annotate("text",x = xmin+2, y = ymin+0.5,label = "Hypoxia",color="white",fontface="bold",size=9)+
  annotate("rect", xmin = xmax-4, xmax =xmax, ymin = ymin+1, ymax = ymin, fill= "black")  +
  annotate("text",x = xmax-2, y = ymin+0.5,label = "IFNg-res",color="white",fontface="bold",size=9)
P0=ggplotGrob(p0)
##make the small plots--> one per Patient
Final <- list()
Final[[1]]=ggplotGrob(p0)
for (i in 1:length(names(column.col))) {
  PatientMD=SeuratObj@meta.data[SeuratObj@meta.data$Patient==paste0(names(column.col)[i]),]
  groups=PatientMD[,c("Patientcolor","Patient")]
  title=paste0(names(column.col)[i])
  matrix=merge(h,groups,by.x=0,by.y=0,all.x = T)
  matrix$Patientcolor[is.na(matrix$Patientcolor)] <- "gray"
  matrix$Patient=as.character(matrix$Patient)
  matrix$Patient[is.na(matrix$Patient)] <- "Other"
  row.names(matrix)=matrix$Row.names
  x=matrix$Patientcolor
  y=matrix$Patient
  col=x[!duplicated(x)]
  names(col)=y[!duplicated(y)]
  matrix=matrix[,-1]
  P=ggplot(matrix, aes(x = X,
                       y =Y,color=factor(Patient)))+geom_point()+geom_point(data = subset(matrix, Patient !="Other"))+
    scale_color_manual(values=col,aesthetics = "colour", breaks = waiver()) +labs(x=NULL,y=NULL)+theme(legend.position = "none")+
    scale_x_continuous(expand = c(0, 0), limits = c(xmin,xmax)) + scale_y_continuous(expand = c(0, 0), limits = c(ymin,ymax))+
    theme(panel.background = element_rect(fill = "white",colour = "white"),axis.ticks.x=element_blank(),axis.text.x=element_blank(),
          axis.ticks.y=element_blank(), axis.text.y=element_blank())+
    ggtitle(title)+theme(plot.title =element_text(size=25,face="bold") )+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+
    geom_hline(yintercept=0, color = "black", size=0.5)+
    geom_vline(xintercept=0, color = "black", size=0.5)+
    annotate("rect", xmin = xmin, xmax = xmin+4, ymin = ymax-1, ymax = ymax, fill= "black")  +
    annotate("text",x = xmin+2, y = ymax-0.5,label = "EMT",color="white",fontface="bold",size=4)+
    annotate("rect", xmin = xmax-4, xmax = xmax, ymin = ymax-1, ymax = ymax, fill= "black")  +
    annotate("text",x = xmax-2, y = ymax-0.5,label = "MYC-Tar",color="white",fontface="bold",size=4)+
    annotate("rect", xmin = xmin, xmax = xmin+4, ymin = ymin+1, ymax = ymin, fill= "black")  +
    annotate("text",x = xmin+2, y = ymin+0.5,label = "Hypoxia",color="white",fontface="bold",size=4)+
    annotate("rect", xmin = xmax-4, xmax =xmax, ymin = ymin+1, ymax = ymin, fill= "black")  +
    annotate("text",x = xmax-2, y = ymin+0.5,label = "IFNg-res",color="white",fontface="bold",size=4)
  Final[[i+1]] = ggplotGrob(P)
}
Final[[1]]=ggplotGrob(p0)
numofplots= length(Final)
pdf(file =paste0(OutputDirectory,ObjName,Subset," pathway butterfly by Patient ordered by type2.pdf"), height = 10, width =20,onefile = T)
grid.arrange(grobs=Final, widths = c(3,1,1,1,1,1),layout_matrix = rbind(c(1, 2, 3,4,NA,5),c(1,6,7,8,9,10),
                                                                        c(1,11,12,13,14,NA),
                                                                        c(1,15,16,17,18,19)))
dev.off()


###
saveRDS(SeuratObj,file = paste0(RobjDirectory,"GliomaClusters-7-21-21 patients renamed.rds"))
#SeuratObj=readRDS(paste0(RobjDirectory,"GliomaClusters-7-21-21 patients renamed.rds"))
