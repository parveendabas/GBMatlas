

SeuratObj=readRDS(paste0(RobjDirectory,"Allhuman-11-3-21 patients renamed.rds"))

ClusterColors=c(C1 = "#003F5C", C2 = "#196687", C3 = "#2879CD", C4 = "#8965B3",
C5 = "#AB599D", C6 = "#D75D96", C7 = "#E56674", C8 = "#EE5761",
C9 = "#EB7B27", C10 = "#FFA206", C11 = "#A3BD5C", C12 = "#8DB032"
)

pdf(paste0(OutputDirectory,ObjName,Subset,"Cluster UMAP Iteration by sample suppfig2a.pdf"),width = 20,height =10,family = "ArialMT")
DimPlot(SeuratObj, reduction = "umap",group.by = "Cluster",cols =ClusterColors,split.by = "Patient",ncol = 6, pt.size = 0.8)+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("") + ylab("")+labs(title="Clusters per patient")+
  FontSize(x.title = 16, y.title = 16,main = 18)+ theme(legend.position = "right",text =element_text(size=15) )+
  guides(colour = guide_legend(title="Clusters", override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=15)))
dev.off()



#suppfig2d
AssignmentColors=c(TCells = "#8DB032", BCells = "#F9BF31", Other = "#E47B1E", 
                   Myeloid = "#CD4D53", Endo = "#DC73B0", Oligo = "#966CBF", Pericytes = "#2960A1", 
                   Glioma = "#003F5C")
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
mylegend<-get_legend(Final[[1]]+theme_min()+ NoAxes()+guides(colour = guide_colourbar(title="Assignment",title.position = "top",title.theme = element_text(size=10))))

pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"piecharts Patient and Assignment suppfigure 2d.pdf")), height = 6, width =10,onefile = T)

ggpubr::ggarrange(Final[[1]]+NoLegend(),Final[[2]]+NoLegend(),NA,NA,NA,NA,Final[[3]]+NoLegend(),Final[[4]]+NoLegend(),
                  Final[[5]]+NoLegend(),Final[[6]]+NoLegend(),Final[[7]]+NoLegend(),Final[[8]]+NoLegend(),
                  Final[[9]]+NoLegend(),Final[[10]]+NoLegend(),Final[[11]]+NoLegend(),Final[[12]]+NoLegend(),Final[[13]]+NoLegend(),NA,
                  Final[[14]]+NoLegend(),Final[[15]]+NoLegend(),Final[[16]]+NoLegend(),Final[[17]]+NoLegend(),Final[[18]]+NoLegend(),NA,
                  ncol = 6,nrow = 4,
                  legend.grob = mylegend,legend = "right",widths = c(1,1,1,1,1,1)) 
dev.off()

###
SeuratObj=readRDS(paste0(RobjDirectory,"GliomaClusters-11-3-21 patients renamed.rds"))

SampleColors=c(`LGG-03` = "#003F5C", `LGG-04` = "#008796", `ndGBM-01` = "#2F4B7C", 
               `ndGBM-02` = "#227ED4", `ndGBM-03` = "#665191", `ndGBM-04` = "#9E71C7", 
               `ndGBM-05` = "#A05195", `ndGBM-06` = "#E082C1", `ndGBM-07` = "#D45087", 
               `ndGBM-08` = "#F57689", `ndGBM-09` = "#BD3E3E", `ndGBM-10` = "#F95D6A", 
               `ndGBM-11` = "#E07B18", `rGBM-01` = "#FF7C43", `rGBM-02` = "#FFA600", 
               `rGBM-03` = "#F5D256", `rGBM-04` = "#42A665", `rGBM-05` = "#8DB032"
)

CellInfo=SeuratObj@meta.data
set.seed(12)
Idents(SeuratObj)=CellInfo$Patient
sampled.cells <- sample(x = SeuratObj@active.ident, size = 10000, replace = F)
subset=subset(SeuratObj,cells=names(sampled.cells))

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
sc=read.csv("~/Box/Yun lab projects/scRNAseq Project/Human Glioma/ALLSAMPLESJuly2021/Output/Glioma/HumanGliomacellstates20000cells.csv",row.names = 1)
#get coordinates
h=hierarchy (sc, quadrants = c("AC","OPC","Mesenchymal","Proneural"), log.scale = T)
# make plots
##Perpatient
CellInfo=SeuratObj@meta.data
for( i in 1: length(names(SampleColors))){
  name=names(SampleColors)[i]
  color=SampleColors[i]
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
for (i in 1:length(names(SampleColors))) {
  PatientMD=SeuratObj@meta.data[SeuratObj@meta.data$Patient==paste0(names(SampleColors)[i]),]
  groups=PatientMD[,c("Patientcolor","Patient")]
  title=paste0(names(SampleColors)[i])
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
pdf(file =paste0(OutputDirectory,ObjName,Subset," subtype butterfly by Patient ordered by type3-suppfig5e.pdf"), height = 10, width =22,onefile = T)
grid.arrange(grobs=Final, widths = c(3,1,1,1,1,1,1),layout_matrix = rbind(c(1, 2,NA, NA,NA,NA,3),c(1,4,5,6,7,8,9),
                                                                        c(1,10,11,12,13,14,NA),
                                                                        c(1,15,16,17,18,19,NA))) 
dev.off()

#Suppfig4a heatmap
SeuratObj=readRDS(paste0(RobjDirectory,"MyeloidClusters-8-10-21-patients renamed.rds"))
ClusterColors=c(MC01 = "#003F5C", MC02 = "#E25159", MC03 = "#73589E", MC04 = "#8DB032", 
                MC05 = "#F6CC4B", MC06 = "#B863A5", MC07 = "#F77B38", MC08 = "#2D5186", 
                MC09 = "#E46388")
AssignmentColors=c(`a-microglia` = "#003F5C", `s-mac 1` = "#2D5187", `AP-microglia` = "#74599E", 
                   DCs = "#B863A5", `h-microglia` = "#E46388", `i-microglia` = "#E25159", 
                   `MDSC` = "#F77B38", Proliferating = "#F6CC4B", `s-mac 2` = "#8DB032"
)

#Figure 4a
Idents(SeuratObj)=CellInfo$Cluster
markers = FindAllMarkers(SeuratObj,logfc.threshold = 0.25,test.use = "wilcox",only.pos = T )
write.csv(markers,paste0(OutputDirectory,ObjName,Subset," cluster markers ","res",resolution,".csv"))

#markers =read.csv("~/Box/Yun lab projects/scRNAseq Project/Human Glioma/ALLSAMPLESJuly2021/Output/Myeloid-no Tcells or Astrocytes/HumanGlioma-AllSamples-Myeloid-no Tcells or Astrocytes cluster markers res0.25.csv")

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
top2 <- markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
genes= top2$gene

rows=rowAnnotation(sel = anno_mark(at = match(genes,row.names(plot.data)), labels = genes,labels_gp =gpar(col = "black", fontsize = 9)))

col<- circlize::colorRamp2(breaks = c(-2, 0, 2), colors = c("#007dd1", "white", "#ab3000"))
HM=Heatmap(name="logcounts",as.matrix(plot.data),cluster_rows =F,cluster_columns = F,top_annotation = colanno, right_annotation = rows,row_names_gp = gpar(fontsize=5),
           col = col,show_column_names= F,show_row_names=F,border = F,show_heatmap_legend = F,use_raster = T)        
lgd=Legend(title = "logcounts", at=  c(-2,0, 2),col_fun = col)
lgd1=Legend(labels = levels(as.factor(column_annot$Cluster)),title="Cluster",legend_gp = gpar(fill=ClusterColors,fontsize=5))
lgd2=Legend(labels = levels(as.factor(column_annot$Patient)),title="Patient",legend_gp = gpar(fill=SampleColors,fontsize=5))
lgd4=Legend(labels = levels(as.factor(column_annot$Type)),title="Type",legend_gp = gpar(fill=TypeColors,fontsize=5))

pdf(paste0(OutputDirectory,"heatmap",ObjName,Subset,"res",resolution,"top20 genes per cluster Fig4a.pdf"),width=8,height=5)
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
##suppfig 4b
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
mylegend<-get_legend(Final[[1]]+theme_min()+ NoAxes()+guides(colour = guide_colourbar(title="Assignment",title.position = "top",title.theme = element_text(size=10))))

pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"piecharts Patient and Assignment suppfigure 4b.pdf")), height = 6, width =10,onefile = T)

ggpubr::ggarrange(Final[[1]]+NoLegend(),Final[[2]]+NoLegend(),NA,NA,NA,NA,Final[[3]]+NoLegend(),Final[[4]]+NoLegend(),
                  Final[[5]]+NoLegend(),Final[[6]]+NoLegend(),Final[[7]]+NoLegend(),Final[[8]]+NoLegend(),
                  Final[[9]]+NoLegend(),Final[[10]]+NoLegend(),Final[[11]]+NoLegend(),Final[[12]]+NoLegend(),Final[[13]]+NoLegend(),NA,
                  Final[[14]]+NoLegend(),Final[[15]]+NoLegend(),Final[[16]]+NoLegend(),Final[[17]]+NoLegend(),Final[[18]]+NoLegend(),NA,
                  ncol = 6,nrow = 4,
                  legend.grob = mylegend,legend = "right",widths = c(1,1,1,1,1,1)) 
dev.off()

##suppfig4c
U1=DimPlot(SeuratObj, reduction = "umap",group.by = "Cluster" ,cols =  ClusterColors,label=T)+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="Clusters")+
  guides(colour = guide_legend(override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=15),ncol = 4))

U2=DimPlot(SeuratObj, reduction = "umap",group.by = "Patient",cols = SampleColors,shuffle = T )+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="Patients")+
  guides(colour = guide_legend(override.aes = list(size=4),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=10),ncol = 4))
U1+U2
 ##supp fig5a big heatmap
m_df_H<- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
Filtered.genesets=filter(m_df_H,gs_name=="HALLMARK_OXIDATIVE_PHOSPHORYLATION"|gs_name=="HALLMARK_TNFA_SIGNALING_VIA_NFKB"|gs_name=="HALLMARK_INTERFERON_GAMMA_RESPONSE"|gs_name=="HALLMARK_HYPOXIA")
geneset.all <-   Filtered.genesets[Filtered.genesets$gene_symbol %in% rownames(Counts),]
genesetlist=geneset.all%>% split(x = .$gene_symbol, f = .$gs_name)

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
pdf(paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"large Pathway combined heatmap fixed.pdf")),width=25,height=10)
draw(HM,heatmap_legend_list = list(lgd,lgd4, lgd2,lgd1), heatmap_legend_side = "right",legend_labels_gp =gpar(col = "black", fontsize = 20,fontface="bold"))
dev.off()


#Suppfig7b
pat2=c("LGG-04", "ndGBM-01", "ndGBM-02", "ndGBM-03","ndGBM-11",  "rGBM-01", "rGBM-02", "rGBM-03", "rGBM-04", "rGBM-05")

for(i in 1:length(pat2)){
  x=ggplot(SeuratObj@meta.data%>%filter(Patient==pat2[i]), aes(Patient, fill=Assignment))+geom_bar(stat="count",position="fill",colour = "white",width =0.7)+ 
    scale_fill_manual(values = AssignmentColors)+coord_polar("y",start=0) + facet_grid(~ Fragment) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_blank(),axis.ticks = element_blank(),
          panel.background = element_blank(), axis.title = element_blank(),strip.background=element_rect(fill="white"),
          strip.text=element_text(face = "bold",color = "black",size=8),plot.margin = unit(c(0,-0.2,-0.2,0),'cm'))
  Final[[i]]=x
}
mylegend<-get_legend(Final[[2]]+theme_min()+ NoAxes()+guides(colour = guide_colourbar(title="Assignment",title.position = "top",title.theme = element_text(size=10))))

pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"piecharts fragment and assignment suppfigure 5b.pdf")), height = 6, width =10,onefile = T)
ggpubr::ggarrange(Final[[1]]+NoLegend(),Final[[6]]+NoLegend(),Final[[2]]+NoLegend(),Final[[7]]+NoLegend(),Final[[3]]+NoLegend(),
                  Final[[8]]+NoLegend(),Final[[4]]+NoLegend(),Final[[9]]+NoLegend(),Final[[5]]+NoLegend(),Final[[10]],ncol = 2,nrow = 5,
                  legend.grob = mylegend,legend = "right",widths = c(1,1)) 
dev.off()


#suppfig5a
SeuratObj=readRDS(paste0(RobjDirectory,"GliomaClusters-11-3-21 patients renamed.rds"))
ClusterColors=c(GC01 = "#4BA75E", GC02 = "#003F5C", GC03 = "#235A82", GC04 = "#877EAE", 
                GC05 = "#4367B2", GC06 = "#8DB032", GC07 = "#007D8E", GC08 = "#266AB3", 
                GC09 = "#7B5DA5")
pat2=c("LGG-04", "ndGBM-01", "ndGBM-02", "ndGBM-03","ndGBM-11",  "rGBM-01", "rGBM-02", "rGBM-03", "rGBM-04", "rGBM-05")

for(i in 1:length(pat2)){
  x=ggplot(SeuratObj@meta.data%>%filter(Patient==pat2[i]), aes(Patient, fill=Cluster))+geom_bar(stat="count",position="fill",colour = "white",width =0.7)+ 
    scale_fill_manual(values = ClusterColors)+coord_polar("y",start=0) + facet_grid(~ Fragment) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_blank(),axis.ticks = element_blank(),
          panel.background = element_blank(), axis.title = element_blank(),strip.background=element_rect(fill="white"),
          strip.text=element_text(face = "bold",color = "black",size=8),plot.margin = unit(c(0,-0.2,-0.2,0),'cm'))
  Final[[i]]=x
}
mylegend<-get_legend(Final[[2]]+theme_min()+ NoAxes()+guides(colour = guide_colourbar(title="Cluster",title.position = "top",title.theme = element_text(size=10))))

pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"piecharts fragment and assignment suppfigure 5a.pdf")), height = 6, width =10,onefile = T)
ggpubr::ggarrange(Final[[1]]+NoLegend(),Final[[6]]+NoLegend(),Final[[2]]+NoLegend(),Final[[7]]+NoLegend(),Final[[3]]+NoLegend(),
                  Final[[8]]+NoLegend(),Final[[4]]+NoLegend(),Final[[9]]+NoLegend(),Final[[5]]+NoLegend(),Final[[10]],ncol = 2,nrow = 5,
                  legend.grob = mylegend,legend = "right",widths = c(1,1)) 
dev.off()


#suppfig7c
CellInfo=SeuratObj@meta.data
set.seed(12)
Idents(SeuratObj)=CellInfo$Patient
sampled.cells <- sample(x = SeuratObj@active.ident, size = 10000, replace = F)
subset=subset(SeuratObj,cells=names(sampled.cells))

Counts <- as.matrix(GetAssayData(object = subset, assay = "RNA", slot = "counts"))

#get scores
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
  annotate("text",x = xmin+2.5, y = ymax-1,label = "EMT",color="white",fontface="bold",size=9)+ 
  annotate("rect", xmin = xmax-5, xmax = xmax, ymin = ymax-2, ymax = ymax, fill= "black")  +
  annotate("text",x = xmax-2.5, y = ymax-1,label = "MYC-Tar",color="white",fontface="bold",size=9)+
  annotate("rect", xmin = xmin, xmax = xmin+5, ymin = ymin+2, ymax = ymin, fill= "black")  + 
  annotate("text",x = xmin+2.5, y = ymin+1,label = "Hypoxia",color="white",fontface="bold",size=9)+ 
  annotate("rect", xmin = xmax-5, xmax =xmax, ymin = ymin+2, ymax = ymin, fill= "black")  +
  annotate("text",x = xmax-2.5, y = ymin+1,label = "IFNg-res",color="white",fontface="bold",size=9)  
P0=ggplotGrob(p0)
##make the small plots--> one per Patient
Final <- list()
Final[[1]]=ggplotGrob(p0)
for (i in 1:length(names(SampleColors))) {
  PatientMD=SeuratObj@meta.data[SeuratObj@meta.data$Patient==paste0(names(SampleColors)[i]),]
  groups=PatientMD[,c("Patientcolor","Patient")]
  title=paste0(names(SampleColors)[i])
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
    annotate("text",x = xmin+2.5, y = ymax-1,label = "EMT",color="white",fontface="bold",size=4)+ 
    annotate("rect", xmin = xmax-5, xmax = xmax, ymin = ymax-2, ymax = ymax, fill= "black")  +
    annotate("text",x = xmax-2.5, y = ymax-1,label = "MYC-Tar",color="white",fontface="bold",size=4)+
    annotate("rect", xmin = xmin, xmax = xmin+5, ymin = ymin+2, ymax = ymin, fill= "black")  + 
    annotate("text",x = xmin+2.5, y = ymin+1,label = "Hypoxia",color="white",fontface="bold",size=4)+ 
    annotate("rect", xmin = xmax-5, xmax =xmax, ymin = ymin+2, ymax = ymin, fill= "black")  +
    annotate("text",x = xmax-2.5, y = ymin+1,label = "IFNg-res",color="white",fontface="bold",size=4)  
  Final[[i+1]] = ggplotGrob(P)
}
Final[[1]]=ggplotGrob(p0)
numofplots= length(Final)
pdf(file =paste0(OutputDirectory,ObjName,Subset," subtype butterfly by Patient ordered by type3-suppfig7c.pdf"), height = 10, width =22,onefile = T)
grid.arrange(grobs=Final, widths = c(3,1,1,1,1,1,1),layout_matrix = rbind(c(1, 2,NA, NA,NA,NA,3),c(1,4,5,6,7,8,9),
                                                                          c(1,10,11,12,13,14,NA),
                                                                          c(1,15,16,17,18,19,NA))) 
dev.off()

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
title="All Fragments"
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
  annotate("text",x = xmin+2, y = ymax-0.5,label = "EMT",color="white",fontface="bold",size=9)+ 
  annotate("rect", xmin = xmax-4, xmax = xmax, ymin = ymax-1, ymax = ymax, fill= "black")  +
  annotate("text",x = xmax-2, y = ymax-0.5,label = "MYC-Tar",color="white",fontface="bold",size=9)+
  annotate("rect", xmin = xmin, xmax = xmin+4, ymin = ymin+1, ymax = ymin, fill= "black")  + 
  annotate("text",x = xmin+2, y = ymin+0.5,label = "Hypoxia",color="white",fontface="bold",size=9)+ 
  annotate("rect", xmin = xmax-4, xmax =xmax, ymin = ymin+1, ymax = ymin, fill= "black")  +
  annotate("text",x = xmax-2, y = ymin+0.5,label = "IFNg-res",color="white",fontface="bold",size=9) 
P0=ggplotGrob(p0)
##make the small plots--> one per Fragment

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
pdf(file =paste0(OutputDirectory,ObjName,Subset," pathway butterfly by Fragment suppfigure 7c.pdf"), height = 15, width =37,onefile = T)
grid.arrange(grobs=Final, widths = c(4,1,1,1,1,1,1,1,1,1,1,1),layout_matrix = rbind(c(1, 2, 3,4,5,NA,NA,NA,NA,NA,NA,NA),c(1,6,7,8,9,10, 11, 12,13,14,15,16),
                                                                                    c(1,17,18,19, 20, 21,22,23,24,25,26,27),
                                                                                    c(1, 28,29, 30,31,32,33,34,35,36,37,38),
                                                                                    c(1, 39,40,41,42,43,44,45,NA,NA,NA,NA))) 
dev.off()