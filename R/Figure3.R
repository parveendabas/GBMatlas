library(easypackages)
MyPackages<-c("dplyr","ggplot2","ggpubr","gridExtra","viridis","egg",
              "grid","lattice","gtools","Biobase","RColorBrewer","tibble",
              "Seurat","cowplot","patchwork","stringr","ComplexHeatmap","SingleCellExperiment", "ggmin")
libraries(MyPackages)

#Set working directory
Directory="~/Box/Yun lab projects/scRNAseq Project/Human Glioma/ALLSAMPLESJuly2021/"
#dir.create(Directory)

setwd(Directory)

#Change variables

ObjName= "HumanGlioma-AllSamples-"
Subset="TCells"
resolution= 0.5
dimz=1:30
abb="TC"

RobjDirectory=paste0(Directory,"R_Objects/")
OutputDirectory=paste0(Directory,"Output/",Subset,"/")
#dir.create(OutputDirectory)
featuremapDirectory=paste0(OutputDirectory,"Featuremaps/")
#dir.create(featuremapDirectory)

#Load File
SeuratObj=readRDS(paste0(RobjDirectory,"TcellClusters-7-22-21-goodwithoutcluster6 new patient names.rds"))

#Extract Meta Data
CellInfo <- SeuratObj@meta.data

#Assign colors
SampleColors=c(`LGG-01` = "#003F5C", `LGG-02` = "#008796", `LGG-03` = "#2F4B7C", 
               `LGG-04` = "#227ED4", `ndGBM-01` = "#665191", `ndGBM-02` = "#9E71C7", 
               `ndGBM-03` = "#A05195", `ndGBM-04` = "#E082C1", `ndGBM-05` = "#D45087", 
               `ndGBM-06` = "#F57689", `ndGBM-07` = "#BD3E3E", `ndGBM-08` = "#F95D6A", 
               `ndGBM-09` = "#E07B18", `rGBM-01` = "#FF7C43", `rGBM-02` = "#FFA600", 
               `rGBM-03` = "#F5D256", `rGBM-04` = "#42A665", `rGBM-05` = "#8DB032"
)


ClusterColors=c(TC01 = "#8DB032", TC02 = "#42A665", TC03 = "#9E71C7", TC04 = "#665191", 
                TC05 = "#227ED4", TC06 = "#2F4B7C", TC07 = "#008796", TC08 = "#003F5C"
)

FragmentColors=c(`LGG-01-A` = "#003F5C", `LGG-01-B` = "#005B72", `LGG-01-C` = "#007789", 
                 `LGG-01-D` = "#087B91", `LGG-02` = "#1B6486", `LGG-03` = "#2D4C7C", 
                 `LGG-04-1` = "#2A5D9C", `LGG-04-2` = "#2572BF", `LGG-04-3` = "#2D76C9", 
                 `ndGBM-01-A` = "#4764AE", `ndGBM-01-C` = "#625394", `ndGBM-01-D` = "#795CA3", 
                 `ndGBM-01-F` = "#8F68B9", `ndGBM-02-1` = "#9E6CC0", `ndGBM-02-2` = "#9F5FAC", 
                 `ndGBM-02-4` = "#9F5398", `ndGBM-02-5` = "#B460A3", `ndGBM-03-1` = "#CE74B4", 
                 `ndGBM-03-2` = "#DE7CBA", `ndGBM-03-3` = "#D968A3", `ndGBM-04` = "#D5548C", 
                 `ndGBM-05` = "#DD5B87", `ndGBM-06` = "#EB6A88", `ndGBM-07` = "#EF7082", 
                 `ndGBM-08` = "#D95A64", `ndGBM-09` = "#C34446", `rGBM-01-A` = "#CD464A", 
                 `rGBM-01-B` = "#E5525B", `rGBM-01-C` = "#F75F64", `rGBM-01-D` = "#ED6A43", 
                 `rGBM-02-2` = "#E37623", `rGBM-02-3` = "#E77B23", `rGBM-02-4` = "#F47B34", 
                 `rGBM-02-5` = "#FF7D3F", `rGBM-03-1` = "#FF8E25", `rGBM-03-2` = "#FF9F0A", 
                 `rGBM-03-3` = "#FCB013", `rGBM-04-1` = "#F8C136", `rGBM-04-2` = "#F0D056", 
                 `rGBM-04-3` = "#AABF5C", `rGBM-04-4` = "#63AE62", `rGBM-05-1` = "#51A85A", 
                 `rGBM-05-2` = "#6FAC46", `rGBM-05-3` = "#8DB032")
AssignmentColors=c(`CD4 TCells` = "#003F5C", `CD8 TCells` = "#74599E", `Naiive TCells` = "#E46388", 
                   NKCells = "#F77B38", Tregs = "#8DB032")
TypeColors=c(Astrocytoma = "#003F5C", GBM = "#795192", Oligodendroglioma = "#EC5873", 
             `Recurrent GBM` = "#FFA600")
GradeColors=c(II = "#FFA600", III = "#BA508E", IV = "#003F5C")

#Figure3a
Idents(SeuratObj)=CellInfo$Cluster
markers = FindAllMarkers(SeuratObj,logfc.threshold = 0.5,test.use = "wilcox",only.pos = T )
write.csv(markers,paste0(OutputDirectory,ObjName,Subset," cluster markers ","res",resolution,".csv"))
top20.markers <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top20.markers$cluster=as.character(top20.markers$cluster)
top20.markers=top20.markers[order(top20.markers$cluster),]
SeuratObj.sce=as.SingleCellExperiment(SeuratObj)
plot.data<-as.data.frame(assay(SeuratObj.sce, "logcounts"))
plot.data<-plot.data[top20.markers$gene,]
plot.data <- plot.data - rowMeans(plot.data)
plot.data=na.omit(plot.data)
#CellInfoS=subset@meta.data
column_annot <-CellInfo[,c("Cluster","Patient","Type"),drop=F]
column_annot$Patient = as.factor(as.character(column_annot$Patient))
column_annot=with(column_annot, column_annot[order(Patient), , drop=F])
column_annot=with(column_annot, column_annot[order(Cluster), , drop=F])
plot.data<-plot.data[,row.names(column_annot)]
column.col= SampleColors


column.colors=list()
column.colors[["Patient"]]<-column.col
column.colors[["Cluster"]]<-ClusterColors
column.colors[["Type"]]<-TypeColors

Patient=as.matrix(column_annot[,c("Patient"),drop=F])
Cluster=as.matrix(column_annot[,c("Cluster"),drop=F])
Type=as.matrix(column_annot[,c("Type"),drop=F])

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
lgd2=Legend(labels = levels(as.factor(column_annot$Patient)),title="Patient",legend_gp = gpar(fill=column.col,fontsize=5))
lgd4=Legend(labels = levels(as.factor(column_annot$Type)),title="Type",legend_gp = gpar(fill=TypeColors,fontsize=5))

pdf(paste0(OutputDirectory,"heatmap",ObjName,Subset,"res",resolution,"top20 genes per cluster figure3a.pdf"),width=8,height=4.5)
draw(HM,heatmap_legend_list = list( lgd,lgd1,lgd4, lgd2), heatmap_legend_side = "right")
dev.off()

## Figure 3b
U1=DimPlot(SeuratObj, reduction = "umap",group.by = "Cluster" ,cols =  ClusterColors,label=T)+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "right")+labs(title="Clusters")+
  guides(colour = guide_legend(override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=15),ncol =1))

pdf(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"Clusters UMAP fig 3b.pdf"),width = 7.5,height =5,family = "ArialMT")
U1
dev.off()


#Figure 3c

Genes=c('CD3E','CD4','CD8A','TCF7','CCR7','LAG3','PDCD1','FOXP3','IL2RA','CTLA4','TIGIT','KLRB1','HOPX','S100A4','JUNB','GZMB','IFNG','NKG7','GNLY','MKI67','GZMH','KLRD1')
pdf(paste0(OutputDirectory,"dotplot ",ObjName,Subset, " Fig 3c.pdf"),width=7,height=3.5)
DotPlot(SeuratObj,group.by = "Cluster",dot.scale = 5 ,features = Genes ,scale = T)+scale_colour_viridis_c(option = "plasma")+
  ggmin::theme_min()+ RotatedAxis()+theme(axis.title = element_blank())

dev.off()

#Figure 3d

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

pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"Figure 3d.pdf")), height = 4, width =6,onefile = T)
ggpubr::ggarrange(x1+NoLegend(),x2+NoLegend(),x3+NoLegend(),ncol =1,nrow = 3,
                  legend.grob = mylegend,legend = "right") 
dev.off()




#Figure 3e

for(i in 1:length(levels(as.factor(SeuratObj@meta.data$Patient)))){
  x=ggplot(SeuratObj@meta.data%>%filter(Patient==levels(as.factor(SeuratObj@meta.data$Patient))[i]), aes(Patient, fill=Assignment))+geom_bar(stat="count",position="fill",colour = "white",width =0.7)+ 
    scale_fill_manual(values = AssignmentColors)+coord_polar("y",start=0) + #facet_grid(~ Fragment) + #
    labs(title=levels(as.factor(SeuratObj@meta.data$Patient))[i])+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_blank(),axis.ticks = element_blank(),
          panel.background = element_blank(), axis.title = element_blank(),strip.background=element_rect(fill="white"),
          strip.text=element_text(face = "bold",color = "black",size=8),plot.margin = unit(c(0,-0.2,-0.2,0),'cm'))
  Final[[i]]=x
}
mylegend<-get_legend(Final[[1]]+theme_min()+ NoAxes()+guides(colour = guide_colourbar(title=Y,title.position = "top",title.theme = element_text(size=10))))

pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"piecharts Patient and Assignment figure 3e.pdf")), height = 6, width =10,onefile = T)

ggpubr::ggarrange(Final[[1]]+NoLegend(),Final[[2]]+NoLegend(),Final[[3]]+NoLegend(),Final[[4]]+NoLegend(),NA,
                  Final[[5]]+NoLegend(),Final[[6]]+NoLegend(),Final[[7]]+NoLegend(),Final[[8]]+NoLegend(),NA,
                  Final[[9]]+NoLegend(),Final[[10]]+NoLegend(),Final[[11]]+NoLegend(),Final[[12]]+NoLegend(),Final[[13]]+NoLegend(),
                  Final[[14]]+NoLegend(),Final[[15]]+NoLegend(),Final[[16]]+NoLegend(),Final[[17]]+NoLegend(),Final[[18]]+NoLegend(),
                  ncol = 5,nrow = 4,
                  legend.grob = mylegend,legend = "right",widths = c(1,1,1,1,1)) 
dev.off()




#Figure 3f
pdf(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"violinplot fig 3f.pdf"),width = 5.7,height =4,family = "ArialMT")
VlnPlot(SeuratObj,features="PDCD1", group.by = "Patient", cols = SampleColors,log = T)+theme_min()+RotatedAxis()
dev.off()