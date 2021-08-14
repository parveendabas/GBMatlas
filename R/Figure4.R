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
Subset="Myeloid-no Tcells or Astrocytes"
resolution= 0.25
dimz=1:30
abb="MC"

RobjDirectory=paste0(Directory,"R_Objects/")

OutputDirectory=paste0(Directory,"Output/",Subset,"/")
#dir.create(OutputDirectory)
featuremapDirectory=paste0(OutputDirectory,"Featuremaps/")
#dir.create(featuremapDirectory)

#Load File
SeuratObj=readRDS(paste0(RobjDirectory,"MyeloidClusters-8-10-21-patients renamed.rds"))

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


ClusterColors=c(MC01 = "#003F5C", MC02 = "#E25159", MC03 = "#73589E", MC04 = "#8DB032", 
                MC05 = "#F6CC4B", MC06 = "#B863A5", MC07 = "#F77B38", MC08 = "#2D5186", 
                MC09 = "#E46388")

AssignmentColors=c(`a-microglia` = "#003F5C", `s-mac1` = "#2D5187", `AP-microglia` = "#74599E", 
                   DCs = "#B863A5", `h-microglia` = "#E46388", `i-microglia` = "#E25159", 
                   `MDSC` = "#F77B38", Proliferating = "#F6CC4B", `s-mac2` = "#8DB032"
)


#Figure 4a
U1=DimPlot(SeuratObj, reduction = "umap",group.by = "Cluster" ,cols =  ClusterColors,label=T)+ggmin::theme_min()+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=NULL)+ xlab("UMAP1") + ylab("UMAP2")+
  FontSize(x.title = 20, y.title = 20,main = 20)+ theme(legend.position = "bottom")+labs(title="Clusters")+
  guides(colour = guide_legend(override.aes = list(size=5),title.theme = element_text(size=15, face="bold"),title.position = "top",label.theme = element_text(size=15),ncol = 5))

pdf(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"Clusters UMAP Figure 4b.pdf"),width = 5.5,height =6,family = "ArialMT")
U1
dev.off()


#Figure 4b
pdf(paste0(OutputDirectory,"dotplot ",ObjName,Subset, "by Cluster figure4c.pdf"),width=10,height=3.5)

DotPlot(SeuratObj,group.by = "Cluster" ,features = c('ITGAM','ITGAX','CD14','CD33','CD68',"MSR1"
                                                     ,'HLA-DRB1',"CD74","B2M",'IFNGR1','CD1C','S100A4','MIF','LYZ',
                                                     'CX3CR1','TMEM119','P2RY12','BHLHE41','SPRY1','CCL4',
                                                     'CCL3','ITGA4','IL1B','IL6','IL12A','IL23A','CXCL9',"NFKBIZ"
                                                     ,'NOS2','CD80','CD86','TNF','MRC1','IL1R1',
                                                     'CCL17','TGFB1','IGF1','FN1','CCL1','IL10','TNFSF14',
                                                     'MERTK','CXCL13','CD163','VEGFA',"MKI67"),
        scale = T)+RotatedAxis()+RotatedAxis()+scale_colour_viridis_c(option = "plasma")

dev.off()
#Figure 4c
pdf(paste0(OutputDirectory,ObjName,Subset, "bargraph fig4c.pdf"),width=8.5,height=5)

CellInfo %>%   ggplot(aes(x=Cluster, fill=fct_rev(Patient))) +scale_fill_manual(values = SampleColors)+
  geom_bar(color="black",position = "fill",width = 0.7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.text = element_text(color="black",size=13,face = "bold"),legend.title = element_text(color="black",size=13,face = "bold"),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15, face="bold"), axis.title=element_text(size=15,face="bold"), 
        plot.margin = unit(c(0.2, 0.5,0.2, 0.5),"cm"))+labs(y ="Composition(percentage of cells)", x= NULL)+ 
  scale_y_continuous(expand = c(0,0),limits = c(0,1),breaks = c(0,0.20,0.40,0.60,0.80,1),labels = scales::percent)+
  theme(legend.title = element_blank())

dev.off()
#Figure 4d-> made in prism
write.csv(as.matrix(table(SeuratObj@meta.data$Assignment,SeuratObj@meta.data$Patient)),
          file = paste0(OutputDirectory,ObjName,Subset," number of cells per Assignment and sample to make fig4d.csv"))
#Figure 4e
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
         main = paste("                                          ",name,objname," (NES)"),filename = paste("GSEA",name," (NES) top and bottom 2each-Heatmap Fig4e.pdf"))

#Figure 4f
library(gprofiler2)
library(stringr)
library(dplyr)
library(GOplot)
Myeloid.Markers=Clusters.genes
names(Myeloid.Markers)[1] <- "ID"
names(Myeloid.Markers)[8] <- "adj.P.Val"
names(Myeloid.Markers)[2] <- "cluster"

Clusters=levels(as.factor(CellInfo$Cluster))
for (i in 1: length(Clusters)){
  X=  Clusters[i]
  Genes<-Clusters.genes %>%
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
Final=list()
for (i in 1:9){
  MC=paste0("MC0",i)
  print(MC)
  DEgenes=Myeloid.Markers%>%dplyr::filter(cluster==MC)
  df=eval(parse(text=MC))
  circ<-circle_dat(df,DEgenes)
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
Genes=c("ITGAM","ITGAX","TGFBI","CD68","CD163","CD1C","BATF3","P2RY12")
FeaturePlot(SeuratObj,features=Genes,order = T,ncol = 4,label = T,label.size = 2,cols = Nour_cols(c("darkpurple","lightorange")))
for (i in 1:length(c(Genes))){
  Y=Genes[i]
  X=FeaturePlot(SeuratObj,features=Y,order = T,label = T,label.size = 2,cols = Nour_cols(c("darkpurple","lightorange")))+theme_min()
  pdf(paste0(featuremapDirectory,Y," featuremap ",ObjName," ",Subset,".pdf"),width = 4.5,height = 4)
  print(X)
  dev.off()
  
}


#Figure 4g

Idents(SeuratObj)=CellInfo$Cluster
sampled.cells <- sample(x = SeuratObj@active.ident, size = 20000, replace = F)
subset=subset(SeuratObj,cells=names(sampled.cells))

Counts <- as.matrix(GetAssayData(object = subset, assay = "RNA", slot = "counts"))



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
  annotate("text",x = xmin+2, y = ymax-0.5,label = "IFNG-res",color="white",fontface="bold",size=8)+ 
  annotate("rect", xmin = xmax-4, xmax = xmax, ymin = ymax-1, ymax = ymax, fill= "black")  +
  annotate("text",x = xmax-2, y = ymax-0.5,label = "Hypoxia",color="white",fontface="bold",size=8)+
  annotate("rect", xmin = xmin, xmax = xmin+4, ymin = ymin+1, ymax = ymin, fill= "black")  + 
  annotate("text",x = xmin+2, y = ymin+0.5,label = "Ox-Phos",color="white",fontface="bold",size=8)+ 
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
    annotate("text",x = xmin+2, y = ymax-0.5,label = "IFNG-res",color="white",fontface="bold",size=4)+ 
    annotate("rect", xmin = xmax-4, xmax = xmax, ymin = ymax-1, ymax = ymax, fill= "black")  +
    annotate("text",x = xmax-2, y = ymax-0.5,label = "Hypoxia",color="white",fontface="bold",size=4)+
    annotate("rect", xmin = xmin, xmax = xmin+4, ymin = ymin+1, ymax = ymin, fill= "black")  + 
    annotate("text",x = xmin+2, y = ymin+0.5,label = "Ox-Phos",color="white",fontface="bold",size=4)+ 
    annotate("rect", xmin = xmax-4, xmax =xmax, ymin = ymin+1, ymax = ymin, fill= "black")  +
    annotate("text",x = xmax-2, y = ymin+0.5,label = "TNFA-sig",color="white",fontface="bold",size=4)  
  Final[[i+1]] = ggplotGrob(P)
}
Final[[1]]=ggplotGrob(p0)
numofplots= length(Final)
pdf(file =paste0(OutputDirectory,ObjName,Subset," pathway butterfly by Cluster fig4g.pdf"), height = 6, width =19.5,onefile = T)
grid.arrange(grobs=Final, widths = c(2,1,1,1,1,1),layout_matrix = rbind(c(1, 2, 3,4,5,6),c(1,7,8,9,10,NA))) 
dev.off()




#Figure 4h

#Load Dataset
CGGA.EXP <- read.table("~/Downloads/CGGA.mRNAseq_325.RSEM-genes.20200506.txt",as.is=T,header = T,row.names=1)
AnnotationCGGA <- read.csv("~/Downloads/CGGA.mRNAseq_All_clinical.csv",as.is=T,header = T,row.names=1)
CGGA.GBM=filter(AnnotationCGGA,Histology=="GBM")
#got markers
AllMyeloid.Markers.surv= read.csv("~/Box/Yun lab projects/scRNAseq Project/Human Glioma/ALLSAMPLESJuly2021/Output/Myeloid-no Tcells or Astrocytes/HumanGlioma-AllSamples-Myeloid-no Tcells or Astrocytes cluster markers res0.25.csv")
Myeloid.Markers.surv=AllMyeloid.Markers.surv %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
#write.csv(Myeloid.Markers.surv,"~/Box/Yun lab manuscripts/GBM Single Cells/Nature Communications revision/RCode/files needed for figure 2e and f/Myeloidmarkers.csv")
#Myeloid.Markers.surv=top20.markers

#remove genes from signature gene lists that do not exist in expression matrix ##otherwise the code won't work
Counts=as.data.frame(CGGA.EXP)
geneset.all <-   Myeloid.Markers.surv[ Myeloid.Markers.surv$gene %in% rownames(Counts),]
geneset <- geneset.all[!duplicated(geneset.all$gene),]
genesetlist=geneset.all%>% split(x =.$gene, f = .$cluster,drop = F)
Scounts=Counts[geneset$gene,]
#get scores 
library(scrabble)

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
#sc=read.csv("~/Box/Yun lab projects/scRNAseq Project/Human Glioma/ALLSAMPLESJuly2021//Output/Myeloid-no Tcells or Astrocytes/HumanGlioma-AllSamples-Myeloid-no Tcells or Astrocytesres0.25cluster CGGAscores.csv",row.names = 1)
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
rnames=c("ExpLevelPositive", "GenderMale", "HistologyAnaplastic Oligodendrolgioma", 
         "HistologyAstrocytoma", "HistologyGBM", "HistologyOligodendroglioma", 
         "RecurrenceRecurrent", "RecurrenceSecondary", 
         "MGMT_statusun-methylated", "IDH_mutation_statusWildtype")
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
  mv.fname <- paste(OutputDirectory,ObjName,Subset,"res",resolution,data.subset[i],"_multivariat_results.csv",sep="")
  write.csv(mv.res,mv.fname)
  mv.res=as.data.frame(mv.res)
  Summtable[i]=mv.res$`Pr(>|z|)`
  colnames(Summtable)[i]=YY
  XX <- ggsurvplot(fit1, data = CGGA.sc2, pval = TRUE,pval.cExpression.Levelrd = c(50, 1),pval.size=6,legend="top",palette = c("Blue","Red"),risk.table = F,font.main = c(12, "bold"),legend.title = "status",
                   title=YY, legend.labs = c("Not-Enriched", "Enriched"),font.x = c(14, "bold"),font.y = c(14, "bold"),font.tickslab = c(12, "plain"),xlab="Time (Months)")
  Final[[i]] = XX
}
write.csv(Summtable,paste(OutputDirectory,ObjName,Subset,"res",resolution,"summary_multivariat_results supp table 10.csv",sep=""))

pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"survival plots CGGA allFigure4E.pdf")), height = 3.5, width =16)
arrange_ggsurvplots(Final,print = TRUE, ncol = 5,nrow = 1,surv.plot.height = NULL,risk.table.height = NULL,ncensor.plot.height = NULL)
dev.off()
### to get figure 4f- I subset the matrix and only get the scores for GBM patients
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
#sc=read.csv(paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"cluster GBMonlyCGGAscores.csv")),row.names = 1)
##Survival Analysis
library(survival)
library(survminer)
CGGA.sc=merge(sc,AnnotationCGGA,by.x=0, by.y=0)
CGGA.sc$OS<-as.numeric(CGGA.sc$survival)
CGGA.sc$status<-as.numeric(CGGA.sc$status)
surv_object <- Surv(time = CGGA.sc$OS, event = CGGA.sc$status)
#for lExpression.Levelp to assign "positive" or "negative" to cell scores then draw survival plots
data.subset=colnames(sc)
data.subset=sort(data.subset, decreasing = FALSE)
CGGA.sc2 =CGGA.sc
Final <- list()
rnames=c("ExpLevelPositive", "GenderMale", 
         "RecurrenceRecurrent", "RecurrenceSecondary", 
         "MGMT_statusun-methylated", "IDH_mutation_statusWildtype")
Summtable=as.data.frame(row.names =rnames,x=rnames  )

for( i in 1: length(data.subset))
{
  YY= data.subset[i]
  Expression.Level=paste(YY,"Levels",sep=".")
  CGGA.sc2 <- CGGA.sc2%>% mutate(Expression.Level = ifelse(CGGA.sc2[YY]>=0, "Positive", "Negative"))
  CGGA.sc2$Expression.Level <-factor(CGGA.sc2$Expression.Level)
  fit1 <- survfit(surv_object ~Expression.Level, data = CGGA.sc2)
  #next three lines sdded by Joshy
  fit1.mv <- coxph(surv_object  ~Expression.Level + Gender + Recurrence  + MGMT_status +IDH_mutation_status , data = CGGA.sc2)
  mv.res <- summary(fit1.mv)$coefficients
  mv.fname <- paste(OutputDirectory,ObjName,Subset,"res",resolution,data.subset[i],"_GBM only multivariat_results.csv",sep="")
  write.csv(mv.res,mv.fname)
  mv.res=as.data.frame(mv.res)
  Summtable[i]=mv.res$`Pr(>|z|)`
  colnames(Summtable)[i]=YY
  
  XX <- ggsurvplot(fit1, data = CGGA.sc2, pval = TRUE,pval.cExpression.Levelrd = c(50, 1),pval.size=6,legend="top",palette = c("Blue","Red"),risk.table = F,font.main = c(12, "bold"),legend.title = "status",
                   title=YY, legend.labs = c("Not-Enriched", "Enriched"),font.x = c(14, "bold"),font.y = c(14, "bold"),font.tickslab = c(12, "plain"),xlab="Time (Months)")
  Final[[i]] = XX
}
write.csv(Summtable,paste(OutputDirectory,ObjName,Subset,"res",resolution,"GBM onlysummary_multivariat_results supp table 10.csv",sep=""))

pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"GBM only survivalFigure4F.pdf")), height = 3.5, width =16)
arrange_ggsurvplots(Final,print = TRUE, ncol = 5,nrow = 1,surv.plot.height = NULL,risk.table.height = NULL,ncensor.plot.height = NULL)
dev.off()

