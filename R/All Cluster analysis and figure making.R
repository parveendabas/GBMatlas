library(easypackages)
MyPackages<-c("dplyr","ggplot2","ggpubr","gridExtra","viridis","egg",
              "grid","lattice","gtools","Biobase","RColorBrewer",
              "Seurat","cowplot","patchwork","stringr","ComplexHeatmap","SingleCellExperiment", "ggmin")
libraries(MyPackages)


#Set working directory
Directory="~/Box/Yun lab projects/scRNAseq Project/Human Glioma/ALLSAMPLESJuly2021/"
#dir.create(Directory)

setwd(Directory)
RobjDirectory=paste0(Directory,"R_Objects/")
#dir.create(RobjDirectory)
OutDirectory=paste0(Directory,"Output/")
#dir.create(OutDirectory)
OutputDirectory=paste0(OutDirectory,"AllClusters/")
#dir.create(OutputDirectory)
featuremapDirectory=paste0(OutputDirectory,"Featuremaps/")
#dir.create(featuremapDirectory)

#Load File
SeuratObj <- readRDS("~/Box/Kyuson_GBM/2021_GBM_All_Three_Batches_Analysis_USE_THIS/GBM_ALL_Batches_Correction_minGenes_500_PCA_50_Res_0.1_Theta_0.5_0.5_CellTypes_Added_C13_Removed.rds")


#Change variables
ObjName= "HumanGlioma-AllSamples-"
Subset="All_Clusters"
resolution= 0.1



#Extract Meta Data
CellInfo <- SeuratObj@meta.data

CellInfo$Patient=CellInfo$Sample
# Add Sample color info
SampleColors= Nour_pal("all")(length(levels(as.factor(SeuratObj$Patient))))
names(SampleColors)=levels(as.factor(SeuratObj$Patient))
#Fragment color
CellInfo=SeuratObj@meta.data
CellInfo$Fragment[which(str_detect(row.names(CellInfo), "57"))] <- "CNSTM-394-3"
SeuratObj@meta.data=CellInfo
FragmentColors= Nour_pal("all",reverse = F)(length(levels(as.factor(SeuratObj@meta.data$Fragment))))
names(FragmentColors)=levels(as.factor(SeuratObj@meta.data$Fragment))


#TO Make Barplot with statistics
## First find the optimum Y axis limist
Max=max(as.matrix(table(SeuratObj@meta.data$Sample)))
###function to round it up
roundUpNice <- function(x, nice=c(1,1.5,2,2.5,3,4,5,6,7,8,9,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}
Height=roundUpNice(Max)


##Number of cells Barplot
P1=CellInfo %>%   ggplot(aes(x=Sample, fill=Sample)) +scale_fill_manual(values = SampleColors)+
  geom_bar(color="black",width = 0.7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.text = element_text(color="black",size=13,face = "bold"),legend.title = element_text(color="black",size=13,face = "bold"),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15, face="bold"), axis.title=element_text(size=15,face="bold"),
        plot.margin = unit(c(0.5, 0.5,0.2, 0.5),"cm"))+labs(y ="Composition(Number of cells)", x= NULL)+ scale_y_continuous(expand = c(0,0),limits = c(0,Height),breaks = seq(0, Height, by = 2000))
pdf(paste0(OutputDirectory,ObjName,Subset," Barplot number of cells per sample.pdf"),width = 8,height = 5.5,family = "ArialMT")
P1
dev.off()

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

Idents(SeuratObj)=CellInfo$Cluster
#markers = FindAllMarkers(SeuratObj,logfc.threshold = 0.5,test.use = "wilcox",only.pos = T )
#write.csv(markers,paste0(OutputDirectory,ObjName,Subset," cluster markers ","res",resolution,".csv"))
markers = read.table("~/Downloads/DEGs_Heatmap_seurat_clusters_Based_ALLcells.txt")
markers=markers[order(markers$cluster),]
top20.markers <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top20.markers=top20.markers[order(top20.markers$cluster),]
set.seed(12)
Idents(SeuratObj)=CellInfo$Sample
#sampled.cells <- sample(x = SeuratObj@active.ident, size = 5000, replace = F)
subset=subset(SeuratObj, downsample=1000)
DimPlot(subset,split.by = "Sample",group.by = "Cluster")
SeuratObj.sce=as.SingleCellExperiment(subset)
plot.data<-as.data.frame(assay(SeuratObj.sce, "logcounts"))
plot.data<-plot.data[top20.markers$gene,]
plot.data <- plot.data - rowMeans(plot.data)
plot.data=na.omit(plot.data)
CellInfoS=subset@meta.data

column_annot <-CellInfoS[,c("Cluster","Patient"),drop=F]
column_annot$Patient = as.factor(as.character(column_annot$Patient))
column_annot=with(column_annot, column_annot[order(Patient), , drop=F])
column_annot=with(column_annot, column_annot[order(Cluster), , drop=F])
plot.data<-plot.data[,row.names(column_annot)]
column.col= SampleColors

column.col2= Nour_pal("all",reverse = F)(length(levels(as.factor(SeuratObj@meta.data$Cluster))))
names(column.col2)=levels(as.factor(SeuratObj@meta.data$Cluster))
column.col3= Nour_pal("all",reverse = F)(length(levels(as.factor(SeuratObj@meta.data$Fragment))))
names(column.col3)=levels(as.factor(SeuratObj@meta.data$Fragment))


column.colors=list()
column.colors[["Patient"]]<-column.col
#column.colors[["Fragment"]]<-column.col3
column.colors[["Cluster"]]<-column.col2
Patient=as.matrix(column_annot[,c("Patient"),drop=F])
Cluster=as.matrix(column_annot[,c("Cluster"),drop=F])
#Fragment=as.matrix(column_annot[,c("Fragment"),drop=F])

colanno <- columnAnnotation(df=column_annot,
                            show_annotation_name =T,show_legend = F,col=column.colors)
top2 <- markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
genes= top2$gene
rows=rowAnnotation(sel = anno_mark(at = match(genes,row.names(plot.data)), labels = genes,labels_gp =gpar(col = "black", fontsize = 9)))

col<- circlize::colorRamp2(breaks = c(-2, 0, 2), colors = c("#007dd1", "white", "#ab3000"))
HM=Heatmap(name="logcounts",as.matrix(plot.data),cluster_rows =F,cluster_columns = F,top_annotation = colanno, right_annotation = rows,row_names_gp = gpar(fontsize=5),
           col = col,show_column_names= F,show_row_names=F,border = F,show_heatmap_legend = F,use_raster = T)
lgd=Legend(title = "logcounts", at=  c(-2,0, 2),col_fun = col)
lgd1=Legend(labels = levels(as.factor(column_annot$Cluster)),title="Cluster",legend_gp = gpar(fill=column.col2,fontsize=5))
lgd2=Legend(labels = levels(as.factor(column_annot$Patient)),title="Patient",legend_gp = gpar(fill=column.col,fontsize=5))
#lgd3=Legend(labels = levels(as.factor(column_annot$Fragment)),title="Fragment",legend_gp = gpar(fill=column.col3,fontsize=5))

pdf(paste0(OutputDirectory,"heatmap",ObjName,Subset,"res",resolution,"top20 genes per cluster nofragx.pdf"),width=7.5,height=9)
draw(HM,heatmap_legend_list = list( lgd,lgd1, lgd2), heatmap_legend_side = "right")
dev.off()

###
MaxS=max(as.matrix(table(SeuratObj@meta.data$Sample)))
###function to round it up
HeightS=roundUpNice(MaxS)
##Number of cells Barplot
P2=CellInfo %>%   ggplot(aes(x=Sample, fill=Cluster)) +scale_fill_manual(values = column.col2)+
  geom_bar(color="black",width = 0.7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.text = element_text(color="black",size=13,face = "bold"),legend.title = element_text(color="black",size=13,face = "bold"),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15, face="bold"), axis.title=element_text(size=15,face="bold"),
        plot.margin = unit(c(0.5, 0.5,0.2, 0.5),"cm"))+labs(y ="Composition(Number of cells)", x= NULL)+ scale_y_continuous(expand = c(0,0),limits = c(0,HeightS),breaks = seq(0, HeightS, by = 1000))
pdf(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"number of cells per sample and clusters barplot.pdf"),width = 6,height = 5.5,family = "ArialMT")
P2
dev.off()
##
P3=CellInfo %>%   ggplot(aes(x=Cluster, fill=Sample)) +scale_fill_manual(values = SampleColors)+
  geom_bar(color="black",position = "fill",width = 0.7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.text = element_text(color="black",size=13,face = "bold"),legend.title = element_text(color="black",size=13,face = "bold"),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15, face="bold"), axis.title=element_text(size=15,face="bold"),
        plot.margin = unit(c(0.2, 0.5,0.2, 0.5),"cm"))+labs(y ="Composition(percentage of cells)", x= NULL)+ scale_y_continuous(expand = c(0,0),limits = c(0,1),breaks = c(0,0.20,0.40,0.60,0.80,1))
pdf(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"percent of cells per sample and clusters barplot.pdf"),width = 6.5,height = 5.5,family = "ArialMT")
P3
dev.off()
##
MaxC=max(as.matrix(table(SeuratObj@meta.data$Cluster)))
HeightC=roundUpNice(MaxC)

P4=CellInfo %>%   ggplot(aes(x=Cluster, fill=Sample)) +scale_fill_manual(values = SampleColors)+
  geom_bar(color="black",width = 0.7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.text = element_text(color="black",size=13,face = "bold"),legend.title = element_text(color="black",size=13,face = "bold"),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15, face="bold"), axis.title=element_text(size=15,face="bold"),
        plot.margin = unit(c(0.5, 0.5,0.2, 0.5),"cm"))+labs(y ="Composition(Number of cells)", x= NULL)+ scale_y_continuous(expand = c(0,0),limits = c(0,HeightC),breaks = seq(0, HeightC, by = 1000))
pdf(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"number of cells per Cluster and Sample barplot.pdf"),width = 6.5,height = 5.5,family = "ArialMT")
P4
dev.off()

#Umaps
U1=DimPlot(SeuratObj, reduction = "umap",group.by = "Cluster" ,cols =  column.col2,label=T)+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="Clusters")+
  guides(colour = guide_legend(override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=15),ncol = 5))

pdf(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"Clusters UMAP.pdf"),width = 5.5,height =6,family = "ArialMT")
U1
dev.off()


U2=DimPlot(SeuratObj, reduction = "umap",group.by = "Patient",cols = SampleColors )+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="Patients")+
  guides(colour = guide_legend(override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=15),ncol = 4))

pdf(paste0(OutputDirectory,ObjName,Subset,"Samples UMAP.pdf"),width = 5,height =6)
U2
dev.off()

pdf(paste0(OutputDirectory,ObjName,Subset,"Cluster UMAP Iteration by sample.pdf"),width = 20,height =10,family = "ArialMT")
DimPlot(SeuratObj, reduction = "umap",group.by = "Cluster",cols =column.col2,split.by = "Patient",ncol = 6, pt.size = 0.8)+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 16, y.title = 16,main = 16)+ theme(legend.position = "right",text =element_text(size=15) )+
  guides(colour = guide_legend(title="Clusters", override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=15)))
dev.off()

## to do dotplot
Genes= c('EGFR','OLIG2','PDGFRA','GFAP','S100B','VIM','MBP','CSPG4',
         'ACTA2','P2RY12','TMEM119','CX3CR1','PTPRC','ITGAM','CD14',
         'CD68','MSR1','HLA-DRA','ITGAX','CD1C','CD3E','CD4','CD8A')

pdf(paste0(OutputDirectory,"dotplot ",ObjName,Subset, "by Cluster fig2#2.pdf"),width=6,height=3.2)
DotPlot(SeuratObj,group.by = "Cluster",dot.scale = 5 ,features = c("PDGFRA","SOX2","OLIG1","NES","GFAP","S100B","EGFR","MBP","PDGFRB","ACTA2",
                                                                   "PROM1","PECAM1","TEK",
                                                                   "P2RY12","PTPRC","ITGAM","ITGAX","S100A8","CD68","MRC1","NKG7","CD3E",
                                                                   "CD8A","CD4","CD79A","MKI67"),cols = c("blue","red"),scale = T)+scale_colour_viridis_c(option = "plasma")+
  ggmin::theme_min()+ RotatedAxis()+theme(axis.title = element_blank())
dev.off()
pdf(paste0(OutputDirectory,"dotplot ",ObjName,Subset, "by Cluster.pdf"),width=7,height=3.5)
DotPlot(SeuratObj,group.by = "Cluster",dot.scale = 5 ,features = Genes ,scale = T)+scale_colour_viridis_c(option = "plasma")+
  ggmin::theme_min()+ RotatedAxis()+theme(axis.title = element_blank())
dev.off()

top2=top2[order(as.character(top2$cluster)),]
genes=top2$gene
pdf(paste0(OutputDirectory,"dotplot ",ObjName,Subset, "top 2 genes by Cluster.pdf"),width=10,height=3.5)
DotPlot(SeuratObj,group.by = "Cluster",dot.scale = 5 ,features = unique(genes) ,scale = T)+scale_colour_viridis_c(option = "plasma")+
  ggmin::theme_min()+ RotatedAxis()+theme(axis.title = element_blank())
dev.off()

## assign clusters

Assign= list()
Assign[["Glioma"]]=c("C2","C6","C9")
Assign[["TCells"]]=c("C3")
Assign[["Myeloid"]]=c("C1","C4","C7")
Assign[["Oligo"]]=c("C5")
Assign[["Endo"]]=c("C10")
Assign[["Pericytes"]]=c("C8")
Assign[["BCells"]]=c("C11")
Assign[["Other"]]=c("C12")
CellInfo$Assignment=NA
for(i in 1:length(Assign)){
  CellInfo$Assignment[CellInfo$Cluster %in% Assign[[i]]]=names(Assign[i])
}
SeuratObj@meta.data=CellInfo
CellInfo=SeuratObj@meta.data

AssignmentColors= Nour_pal("all",reverse = T)(length(levels(as.factor(SeuratObj@meta.data$Assignment))))
#names(AssignmentColors)= levels(as.factor(SeuratObj@meta.data$Assignment))
names(AssignmentColors)=c("TCells", "BCells","Other" ,  "Myeloid",  "Endo", "Oligo",
  "Pericytes","Glioma")
MaxA=max(as.matrix(table(SeuratObj@meta.data$Assignment)))
HeightA=roundUpNice(MaxA)
P5=CellInfo %>%   ggplot(aes(x=Sample, fill=Assignment)) +scale_fill_manual(values = AssignmentColors)+
  geom_bar(color="black",width = 0.7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.text = element_text(color="black",size=13,face = "bold"),legend.title = element_text(color="black",size=13,face = "bold"),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15, face="bold"), axis.title=element_text(size=15,face="bold"),
        plot.margin = unit(c(0.5, 0.5,0.2, 0.5),"cm"))+labs(y ="Composition(Number of cells)", x= NULL)+ scale_y_continuous(expand = c(0,0),limits = c(0,HeightS),breaks = seq(0, HeightS, by = 1000))
pdf(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"number of cells per sample and Assignment barplot.pdf"),width = 9,height = 5.5,family = "ArialMT")
P5
dev.off()


P6=CellInfo %>%   ggplot(aes(x=Assignment, fill=Sample)) +scale_fill_manual(values = SampleColors)+
  geom_bar(color="black",position = "fill",width = 0.7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.text = element_text(color="black",size=13,face = "bold"),legend.title = element_text(color="black",size=13,face = "bold"),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15, face="bold"), axis.title=element_text(size=15,face="bold"),
        plot.margin = unit(c(0.2, 0.5,0.2, 0.5),"cm"))+labs(y ="Composition(percentage of cells)", x= NULL)+ scale_y_continuous(expand = c(0,0),limits = c(0,1),breaks = c(0,0.20,0.40,0.60,0.80,1))
pdf(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"percent of cells per sample and Assignment barplot.pdf"),width = 6,height = 5.5,family = "ArialMT")
P6
dev.off()

P7=CellInfo %>%   ggplot(aes(x=Assignment, fill=Sample)) +scale_fill_manual(values = SampleColors)+
  geom_bar(color="black",width = 0.7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.text = element_text(color="black",size=13,face = "bold"),legend.title = element_text(color="black",size=13,face = "bold"),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15, face="bold"), axis.title=element_text(size=15,face="bold"),
        plot.margin = unit(c(0.5, 0.5,0.2, 0.5),"cm"))+labs(y ="Composition(Number of cells)", x= NULL)+ scale_y_continuous(expand = c(0,0),limits = c(0,HeightA),breaks = seq(0, HeightA, by = 1000))
pdf(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"percentage of cells per major population and Sample barplot.pdf"),width = 7,height = 5.5,family = "ArialMT")
P7
dev.off()
plots=ggarrange(P1,P2,P3,P4,P5,P6,P7,ncol = 2)

pdf(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"all barplots.pdf"),width = 18,height = 18,family = "ArialMT")
print(plots)
dev.off()


#assignment UMAP
U3=DimPlot(SeuratObj, reduction = "umap",group.by = "Assignment" ,cols =  AssignmentColors,label=T)+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="Assignment")+
  guides(colour = guide_legend(override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=15),ncol = 4))

pdf(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"Assignment UMAP.pdf"),width = 5,height =6)
U3
dev.off()
PredictionPalette <- c("#ff7f50", "#009FFF")
names(PredictionPalette)<- c("diploid", "aneuploid")

U4=DimPlot(SeuratObj, reduction = "umap",group.by = "CopyKatPrediction" ,cols =  PredictionPalette,label=T,
           cells = row.names(CellInfo2%>% filter(CopyKatPrediction
                                                 %in% na.omit(SeuratObj$CopyKatPrediction))))+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="CopyKat Prediction")+
  guides(colour = guide_legend(override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=15),ncol = 4))

pdf(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"copyKat UMAP new.pdf"),width = 5,height =6)
U4
dev.off()

plots2=ggarrange(U1,U2,U3,U4,ncol = 2)

pdf(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"all UMAPs2 updated copykat.pdf"),width = 12,height = 10)
print(plots2)
dev.off()
CellInfo$IDfrag
pdf(paste0(OutputDirectory,ObjName,Subset,"Assignment UMAP Iteration by sample.pdf"),width = 8,height =8,family = "ArialMT")
DimPlot(SeuratObj, reduction = "umap",group.by = "Assignment",cols =AssignmentColors,split.by = "Sample",ncol = 2, pt.size = 0.2)+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 16, y.title = 16,main = 16)+ theme(legend.position = "right",text =element_text(size=15) )+
  guides(colour = guide_legend(title="Assignment", override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=15)))
dev.off()
theme_min2 <- function(base_size = 11, base_family = "") {

  theme_light(base_size = 11, base_family = "") +
    theme(     plot.title = element_text(size = rel(0.9),hjust = 0.5,vjust = 0.5),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               panel.border = element_rect(fill = NA, colour = "grey90", size = 1),
               strip.background = element_rect(fill = NA, colour = NA),
               strip.text.x = element_text(colour = "black", size = rel(1.2)),
               strip.text.y = element_text(colour = "black", size = rel(1.2)),
               title = element_text(size = rel(0.9),hjust = 0.5,vjust = 0.5),
               axis.text = element_text(colour = "black", size = rel(0.8)),
               axis.title = element_blank(),
               legend.title = element_text(colour = "black", size = rel(0.9),hjust = 0.5),
               legend.key.size = unit(0.9, "lines"),
               legend.text = element_text(size = rel(0.7), colour = "black"),
               legend.key = element_rect(colour = NA, fill = NA),
               legend.background = element_rect(colour = NA, fill = NA),
               plot.margin = unit(c(0.1,0,0,-0.2), "lines")
    )
}

#
for (i in 1:length(Genes)){
  Y=Genes[i]
  X=FeaturePlot(SeuratObj,features=Y,split.by = "IDfrag",order = T,label = T,label.size = 2,cols = Nour_cols(c("darkpurple","lightorange")),pt.size =0.1)
  mylegend<-get_legend(X[[2]]+theme_min()+ NoAxes()+guides(colour = guide_colourbar(title=Y,title.position = "top",title.theme = element_text(size=10))))

  Z=ggpubr::ggarrange(ggparagraph(text=" ",  size = 0),X[[1]]+theme_min2()+NoLegend() + NoAxes(),
                      X[[2]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[3]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[4]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[5]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[6]]+theme_min2()+NoLegend()+NoAxes() ,ggparagraph(text=" ",  size = 0),
                      X[[7]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[8]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[9]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[10]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[11]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[12]]+theme_min2()+NoLegend()+NoAxes() ,ggparagraph(text=" ",  size = 0),
                      X[[13]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[14]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[15]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[16]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[17]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[18]]+theme_min2()+NoLegend()+NoAxes() ,ggparagraph(text=" ",  size = 0),
                      X[[19]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[20]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[21]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[22]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[23]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[24]]+theme_min2()+NoLegend()+NoAxes() ,ggparagraph(text=" ",  size = 0),
                      X[[25]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[26]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[27]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[28]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[29]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[30]]+theme_min2()+NoLegend()+NoAxes() ,ggparagraph(text=" ",  size = 0),
                      X[[31]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[32]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[33]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[34]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[35]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[36]]+theme_min2()+NoLegend()+NoAxes() ,ggparagraph(text=" ",  size = 0),
                      X[[37]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[38]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[39]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[40]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[41]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[42]]+theme_min2()+NoLegend()+NoAxes() ,ggparagraph(text=" ",  size = 0),
                      X[[43]]+theme_min2()+NoLegend()+NoAxes() ,
                      X[[44]]+theme_min2()+NoLegend()+NoAxes() ,
                      ncol=7,nrow=8,widths = c(0.03,1,1,1,1,1,1),legend.grob = mylegend,legend = "right")
  pdf(paste0(featuremapDirectory,Y," featuremap ",ObjName," ",Subset,".pdf"),width = 8,height = 2)
  print(Z)
  dev.off()

}

#SeuratObj <- RunPCA(SeuratObj, features = c(s.genes, g2m.genes))
DimPlot(SeuratObj,reduction = "umap",group.by = "Phase")
color=Nour_pal("all")(length(levels(as.factor(SeuratObj@meta.data$Phase))))
names(color)=levels(as.factor(SeuratObj@meta.data$Phase))
pdf(paste0(OutputDirectory,ObjName,Subset,"CellCycle UMAP by Iteration Sample.pdf"),width = 8,height =8,family = "ArialMT")
DimPlot(SeuratObj, reduction = "umap",label = T,repel = T,label.size = 2,cols =color,split.by = "Sample",ncol = 2, pt.size = 0.2)+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+labs(title = "Seurat Cell Cycle scoring")+
  FontSize(x.title = 16, y.title = 16,main = 16)+ theme(legend.position = "right",text =element_text(size=15) )+
  guides(colour = guide_legend(title="Cell Cycle", override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=15)))
dev.off()

U4=DimPlot(SeuratObj, reduction = "umap",group.by = "Phase" ,cols = color,label=T,label.color = c("yellow","black","black"))+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="Cell Cycle")+
  guides(colour = guide_legend(override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=15),ncol = 4))

pdf(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"Cell Cycle UMAP.pdf"),width = 5.5,height =6,family = "ArialMT")
U4
dev.off()

plots2=ggarrange(U1,U2,U4,U3,ncol = 4)

pdf(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"all UMAPs2.pdf"),width = 18,height = 5.0,family = "ArialMT")
print(plots2)
dev.off()


for (i in 1:length(c(Genes))){
  Y=Genes[i]
  X=FeaturePlot(SeuratObj,features=Y,split.by = "Sample",order = T,label = T,label.size = 2,cols = Nour_cols(c("darkpurple","lightorange")),pt.size =0.1)
  mylegend<-get_legend(X[[2]]+theme_min()+ NoAxes()+guides(colour = guide_colourbar(title=Y,title.position = "top",title.theme = element_text(size=10))))

  Z=ggpubr::ggarrange(ggparagraph(text=" ",  size = 0),X[[1]]+theme_min2()+NoLegend() + NoAxes(),
                      X[[2]]+theme_min2()+NoLegend()+NoAxes() ,
                      ncol=3,nrow=1,widths = c(0.03,1,1),legend.grob = mylegend,legend = "right")
  pdf(paste0(featuremapDirectory,Y," featuremap ",ObjName," ",Subset,".pdf"),width = 4.5,height = 2)
  print(Z)
  dev.off()

}
Idents(SeuratObj)=SeuratObj$Cluster
for (i in 1:length(c(Genes))){
  Y=Genes[i]
  X=FeaturePlot(SeuratObj,features=Y,order = T,label = T,label.size = 2,cols = Nour_cols(c("darkpurple","lightorange")))+theme_min()
  pdf(paste0(featuremapDirectory,Y," featuremap ",ObjName," ",Subset,".pdf"),width = 4.5,height = 2)
  print(X)
  dev.off()

}


Y="CD68"

###Automated Assignment



table(CellInfo$Assignment,CellInfo$Cluster)
table(CellInfo2$AutoAssignment,CellInfo2$Cluster)


saveRDS(SeuratObj,file = paste0(RobjDirectory,"Allhuman-7-21-21 patients renamed.rds"))
#SeuratObj=readRDS(paste0(RobjDirectory,"Allhuman-7-21-21 patients renamed.rds"))
Idents(SeuratObj)=SeuratObj$Assignment

FeaturePlot(SeuratObj,features="PDCD1", label=T,order = T,cols = Nour_cols(c("darkpurple","lightorange")),pt.size =0.1)+theme_min()
VlnPlot(SeuratObj,features="PDCD1",group.by = "Patient")
tcs= subset(SeuratObj,cells=WhichCells(SeuratObj,idents = "TCells"))
DimPlot(tcs)
VlnPlot(tcs,features="PDCD1",group.by = "Patient",cols =SampleColors ,log = T)+theme_min()+RotatedAxis()
