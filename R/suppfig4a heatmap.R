#Figure 4a
Idents(SeuratObj)=CellInfo$Assignment2
markers = FindAllMarkers(SeuratObj,logfc.threshold = 0.25,test.use = "wilcox",only.pos = T )
write.csv(markers,paste0(OutputDirectory,ObjName,Subset," Assignment2 markers ","res",resolution,".csv"))
markers=markers[order(markers$cluster),]
top20.markers <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top20.markers$cluster=as.character(top20.markers$cluster)
top20.markers=top20.markers[order(top20.markers$cluster),]
set.seed(12)
Idents(SeuratObj)=CellInfo$Sample
subset=subset(SeuratObj, downsample=1000)

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
top2 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
genes= top2$gene

rows=rowAnnotation(sel = anno_mark(at = match(genes,row.names(plot.data)), labels = genes,labels_gp =gpar(col = "black", fontsize = 9)))

col<- circlize::colorRamp2(breaks = c(-2, 0, 2), colors = c("#007dd1", "white", "#ab3000"))
HM=Heatmap(name="logcounts",as.matrix(plot.data),cluster_rows =F,cluster_columns = F,top_annotation = colanno, right_annotation = rows,row_names_gp = gpar(fontsize=5),
           col = col,show_column_names= F,show_row_names=F,border = F,show_heatmap_legend = F,use_raster = T)        
lgd=Legend(title = "logcounts", at=  c(-2,0, 2),col_fun = col)
lgd1=Legend(labels = levels(as.factor(column_annot$Cluster)),title="Cluster",legend_gp = gpar(fill=ClusterColors,fontsize=5))
lgd2=Legend(labels = levels(as.factor(column_annot$Patient)),title="Patient",legend_gp = gpar(fill=SampleColors,fontsize=5))
lgd4=Legend(labels = levels(as.factor(column_annot$Type)),title="Type",legend_gp = gpar(fill=TypeColors,fontsize=5))

pdf(paste0(OutputDirectory,"heatmap",ObjName,Subset,"res",resolution,"top20 genes per cluster Fig4a.pdf"),width=6,height=5)
draw(HM,heatmap_legend_list = list( lgd,lgd1,lgd4, lgd2), heatmap_legend_side = "right")
dev.off()

Final=list()

for(i in 1:length(levels(as.factor(SeuratObj@meta.data$Patient)))){
  x=ggplot(SeuratObj@meta.data%>%filter(Patient==levels(as.factor(SeuratObj@meta.data$Patient))[i]), aes(Patient, fill=Assignment))+geom_bar(stat="count",position="fill",colour = "white",width =0.7)+ 
    scale_fill_manual(values = AssignmentColors)+coord_polar("y",start=0) + #facet_grid(~ Fragment) + #
    labs(title=levels(as.factor(SeuratObj@meta.data$Patient))[i])+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_blank(),axis.ticks = element_blank(),
          panel.background = element_blank(), axis.title = element_blank(),strip.background=element_rect(fill="white"),
          strip.text=element_text(face = "bold",color = "black",size=8),plot.margin = unit(c(0,-0.2,-0.2,0),'cm'))
  Final[[i]]=x
}
mylegend<-get_legend(Final[[2]]+theme_min()+ NoAxes()+guides(colour = guide_colourbar(title=Y,title.position = "top",title.theme = element_text(size=10))))


pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"piecharts patient and Assignment.pdf")), height = 6, width =10,onefile = T)
ggpubr::ggarrange(Final[[1]]+NoLegend(),Final[[2]]+NoLegend(),Final[[3]]+NoLegend(),NA,NA,Final[[4]]+NoLegend(),
                  Final[[5]]+NoLegend(),Final[[6]]+NoLegend(),Final[[7]]+NoLegend(),Final[[8]]+NoLegend(),
                  Final[[9]]+NoLegend(),Final[[10]]+NoLegend(),Final[[11]]+NoLegend(),Final[[12]]+NoLegend(),Final[[13]]+NoLegend(),
                  Final[[14]]+NoLegend(),Final[[15]]+NoLegend(),Final[[16]]+NoLegend(),Final[[17]]+NoLegend(),Final[[18]]+NoLegend(),
                  ncol = 5,nrow = 4,
                  legend.grob = mylegend,legend = "right",widths = c(1,1,1,1,1)) 
dev.off()
