library(easypackages)
MyPackages<-c("dplyr","ggplot2","ggpubr","gridExtra","viridis","egg",
              "grid","lattice","gtools","Biobase","RColorBrewer",
              "Seurat","cowplot","patchwork","stringr","ComplexHeatmap","SingleCellExperiment", "ggmin")
libraries(MyPackages)

#Set working directory
Directory="~/Box/Yun lab projects/scRNAseq Project/Human Glioma/ALLSAMPLESJuly2021/"
#dir.create(Directory)

setwd(Directory)


#Change variables
ObjName= "HumanGlioma-AllSamples-"
Subset="All_Clusters"
resolution= 0.1

RobjDirectory=paste0(Directory,"R_Objects/")
#dir.create(RobjDirectory)
OutDirectory=paste0(Directory,"Output/")
#dir.create(OutDirectory)
OutputDirectory=paste0(OutDirectory,"AllClusters/")
#dir.create(OutputDirectory)
featuremapDirectory=paste0(OutputDirectory,"Featuremaps/")
#dir.create(featuremapDirectory)

#Load File
SeuratObj=readRDS(paste0(RobjDirectory,"Allhuman-11-3-21 patients renamed.rds"))

#Extract Meta Data
CellInfo <- SeuratObj@meta.data

#Assign colors
SampleColors=c(`LGG-01` = "#003F5C", `LGG-03` = "#008796", `LGG-04` = "#2F4B7C", 
               `ndGBM-01` = "#227ED4", `ndGBM-02` = "#665191", `ndGBM-03` = "#9E71C7", 
               `ndGBM-04` = "#A05195", `ndGBM-05` = "#E082C1", `ndGBM-06` = "#D45087", 
               `ndGBM-07` = "#F57689", `ndGBM-08` = "#BD3E3E", `ndGBM-09` = "#F95D6A", 
               `ndGBM-10` = "#E07B18", `rGBM-01` = "#FF7C43", `rGBM-02` = "#FFA600", 
               `rGBM-03` = "#F5D256", `rGBM-04` = "#42A665", `rGBM-05` = "#8DB032"
)


ClusterColors=c(C1 = "#003F5C", C2 = "#196687", C3 = "#2879CD", C4 = "#8965B3", 
                C5 = "#AB599D", C6 = "#D75D96", C7 = "#E56674", C8 = "#EE5761", 
                C9 = "#EB7B27", C10 = "#FFA206", C11 = "#A3BD5C", C12 = "#8DB032"
)

FragmentColors=c(`LGG-01-A` = "#003F5C", `LGG-01-B` = "#005B72", `LGG-01-C` = "#007789", 
                 `LGG-01-D` = "#087B91", `LGG-03` = "#1B6486", `LGG-04-1` = "#2D4C7C", 
                 `LGG-04-2` = "#2A5D9C", `LGG-04-3` = "#2572BF", `ndGBM-01-A` = "#2D76C9", 
                 `ndGBM-01-C` = "#4764AE", `ndGBM-01-D` = "#625394", `ndGBM-01-F` = "#795CA3", 
                 `ndGBM-02-1` = "#8F68B9", `ndGBM-02-2` = "#9E6CC0", `ndGBM-02-4` = "#9F5FAC", 
                 `ndGBM-02-5` = "#9F5398", `ndGBM-03-1` = "#B460A3", `ndGBM-03-2` = "#CE74B4", 
                 `ndGBM-03-3` = "#DE7CBA", `ndGBM-04` = "#D968A3", `ndGBM-05` = "#D5548C", 
                 `ndGBM-06` = "#DD5B87", `ndGBM-07` = "#EB6A88", `ndGBM-08` = "#EF7082", 
                 `ndGBM-09` = "#D95A64", `ndGBM-10` = "#C34446", `rGBM-01-A` = "#CD464A", 
                 `rGBM-01-B` = "#E5525B", `rGBM-01-C` = "#F75F64", `rGBM-01-D` = "#ED6A43", 
                 `rGBM-02-2` = "#E37623", `rGBM-02-3` = "#E77B23", `rGBM-02-4` = "#F47B34", 
                 `rGBM-02-5` = "#FF7D3F", `rGBM-03-1` = "#FF8E25", `rGBM-03-2` = "#FF9F0A", 
                 `rGBM-03-3` = "#FCB013", `rGBM-04-1` = "#F8C136", `rGBM-04-2` = "#F0D056", 
                 `rGBM-04-3` = "#AABF5C", `rGBM-04-4` = "#63AE62", `rGBM-05-1` = "#51A85A", 
                 `rGBM-05-2` = "#6FAC46", `rGBM-05-3` = "#8DB032")

AssignmentColors=c(TCells = "#8DB032", BCells = "#F9BF31", Other = "#E47B1E", 
                   Myeloid = "#CD4D53", Endo = "#DC73B0", Oligo = "#966CBF", Pericytes = "#2960A1", 
                   Glioma = "#003F5C")
PredictionColors=c(diploid = "#ff7f50", aneuploid = "#009FFF")

###Figure1b

colortxt=c("white", "white","black","black","black","black","black", "black","black","black","black","black")

U1=DimPlot(SeuratObj, reduction = "umap",group.by = "Cluster" ,cols =  ClusterColors,label=T,label.color = colortxt)+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="Clusters")+
  guides(colour = guide_legend(override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=15),ncol = 5))

U2=DimPlot(SeuratObj, reduction = "umap",group.by = "Patient",cols = SampleColors )+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="Patients")+
  guides(colour = guide_legend(override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=15),ncol = 4))


colortxt2=c("white", "black","black","black", "black","black","black","black")

U3=DimPlot(SeuratObj, reduction = "umap",group.by = "Assignment" ,cols =  AssignmentColors,label=T,label.color = colortxt2)+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="Assignment")+
  guides(colour = guide_legend(override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=15),ncol = 4))

U4=DimPlot(SeuratObj, reduction = "umap",group.by = "CopyKatPrediction" ,cols =  PredictionColors,label=T,
           cells = row.names(CellInfo%>% filter(CopyKatPrediction 
                                                 %in% na.omit(SeuratObj$CopyKatPrediction)))      )+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="CopyKat Prediction")+
  guides(colour = guide_legend(override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=15),ncol = 4))

plots2=ggarrange(U1,U2,U3,U4,ncol = 2)

pdf(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"Figure 1b.pdf"),width = 12,height = 10)
print(plots2)
dev.off()
#Supp Table 2
write.csv(as.matrix(table(SeuratObj@meta.data$Cluster,SeuratObj@meta.data$Patient)),
          file = paste0(OutputDirectory,ObjName,Subset," Supp Table 2 number of cells per cluster and sample.csv"))

#Figure 1c
Idents(SeuratObj)=CellInfo$Cluster
markers = FindAllMarkers(SeuratObj,logfc.threshold = 0.5,test.use = "wilcox",only.pos = T )
write.csv(markers,paste0(OutputDirectory,ObjName,Subset," Supp Table 3 cluster markers ","res",resolution,".csv"))
#markers=read.csv(paste0(OutputDirectory,"HumanGlioma-AllSamples-All_Clusters cluster markers res0.1.csv"),row.names = 1)
#markers$cnum=stringr::str_remove_all(markers$cluster,"C")
#markers$cnum=as.numeric(markers$cnum)
#markers=markers[order(markers$cnum),]
top20.markers <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
#top20.markers=top20.markers[order(top20.markers$cluster),]
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

column.colors=list()
column.colors[["Patient"]]<-SampleColors
column.colors[["Cluster"]]<-ClusterColors
Patient=as.matrix(column_annot[,c("Patient"),drop=F])
Cluster=as.matrix(column_annot[,c("Cluster"),drop=F])

colanno <- columnAnnotation(df=column_annot,
                            show_annotation_name =T,show_legend = F,col=column.colors)
top2 <- markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
genes= top2$gene
rows=rowAnnotation(sel = anno_mark(at = match(genes,row.names(plot.data)), labels = genes,labels_gp =gpar(col = "black", fontsize = 9)))

col<- circlize::colorRamp2(breaks = c(-2, 0, 2), colors = c("#007dd1", "white", "#ab3000"))
HM=Heatmap(name="logcounts",as.matrix(plot.data),cluster_rows =F,cluster_columns = F,top_annotation = colanno, right_annotation = rows,row_names_gp = gpar(fontsize=5),
           col = col,show_column_names= F,show_row_names=F,border = F,show_heatmap_legend = F,use_raster = T)        
lgd=Legend(title = "logcounts", at=  c(-2,0, 2),col_fun = col)
lgd1=Legend(labels = levels(as.factor(column_annot$Cluster)),title="Cluster",legend_gp = gpar(fill=ClusterColors,fontsize=5))
lgd2=Legend(labels = levels(as.factor(column_annot$Patient)),title="Patient",legend_gp = gpar(fill=SampleColors,fontsize=5))
#lgd3=Legend(labels = levels(as.factor(column_annot$Fragment)),title="Fragment",legend_gp = gpar(fill=column.col3,fontsize=5))

pdf(paste0(OutputDirectory,"heatmap",ObjName,Subset,"res",resolution,"top20 genes per cluster Fig1c.pdf"),width=7.5,height=9)
draw(HM,heatmap_legend_list = list( lgd,lgd1, lgd2), heatmap_legend_side = "right")
dev.off()

#Figure 1d
pdf(paste0(OutputDirectory,"dotplot ",ObjName,Subset, "by Cluster fig 1d.pdf"),width=6,height=3.2)
DotPlot(SeuratObj,group.by = "Cluster",dot.scale = 5 ,features = c("PDGFRA","SOX2","OLIG1","NES","GFAP","S100B","EGFR","MBP","PDGFRB","ACTA2",
                                                                   "PROM1","PECAM1","TEK",
                                                                   "P2RY12","PTPRC","ITGAM","ITGAX","S100A8","CD68","MRC1","NKG7","CD3E",
                                                                   "CD8A","CD4","CD79A","MKI67"),cols = c("blue","red"),scale = T)+scale_colour_viridis_c(option = "plasma")+
  ggmin::theme_min()+ RotatedAxis()+theme(axis.title = element_blank())
dev.off()


#Figure 1e
library(forcats)
P1=CellInfo %>%   ggplot(aes(x=Patient, fill=fct_rev(Cluster))) +scale_fill_manual(values = ClusterColors)+
  geom_bar(color="black",position = "fill",width = 0.7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.text = element_text(color="black",size=13,face = "bold"),legend.title = element_text(color="black",size=13,face = "bold"),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15, face="bold"), axis.title=element_text(size=15,face="bold"), 
        plot.margin = unit(c(0.2, 0.5,0.2, 0.5),"cm"))+labs(y ="Composition(percentage of cells)", x= NULL)+ 
  scale_y_continuous(expand = c(0,0),limits = c(0,1),breaks = c(0,0.20,0.40,0.60,0.80,1),labels = scales::percent)+
  theme(legend.title = element_blank())
P1g=ggplotGrob(P1)
P2=ggplot(CellInfo, aes(Cluster , fill=Cluster))+geom_bar(stat="count",colour = "black",width = 0.7)+  scale_fill_manual(values = ClusterColors)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour = "black"), axis.text.y = element_blank(),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=10),
        axis.text=element_text(size=6.5, face="bold"), axis.title=element_text(size=15,face="bold"),
        plot.margin = unit(c(-0.1, 0.5,2.5, -0.5),"cm"))+labs(y ="number of cells", x= NULL)+ 
  scale_y_continuous(expand=c(0,0),trans ="log2",limits=c(256,131072),oob =  scales::squish)+coord_flip()+
  theme(legend.position = "none")
P2g=ggplotGrob(P2)

pdf(paste0(OutputDirectory,"barplots ",ObjName,Subset, "fig 1e.pdf"),width=12,height=5.4)

grid.arrange(grobs=list(P1g,P2g), widths = c(1,0.35),heights=c(0.19,1),layout_matrix = rbind(c(1, NA),c(1,2))) 
dev.off()

##Figure 1f
SeuratObj@meta.data$Type2=factor(SeuratObj@meta.data$Type, levels=c("Astrocytoma", "Oligodendroglioma", "GBM", "Recurrent GBM"))

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

mylegend<-get_legend(x3+theme_min()+ NoAxes()+guides(colour = guide_colourbar(title=Y,title.position = "top",title.theme = element_text(size=10))))

pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"piecharts type,sexandgrade and assignment fig1f.pdf")), height = 4, width =6,onefile = T)
ggpubr::ggarrange(x1+NoLegend(),x2+NoLegend(),x3+NoLegend(),ncol =1,nrow = 3,
                  legend.grob = mylegend,legend = "right") 
dev.off()


