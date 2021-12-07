library(easypackages)
MyPackages<-c("dplyr","ggplot2","ggpubr","gridExtra","viridis","egg",
              "grid","lattice","gtools","Biobase","RColorBrewer","tibble",
              "Seurat","cowplot","patchwork","stringr","ComplexHeatmap","SingleCellExperiment", "ggmin")
libraries(MyPackages)
#Figure 6A
Directory="~/Box/Yun lab projects/scRNAseq Project/Human Glioma/ALLSAMPLESJuly2021/"
RobjDirectory=paste0(Directory,"R_Objects/")

#Figure5a
#load all cluster r object

SeuratObj=readRDS(paste0(RobjDirectory,"Allhuman-11-3-21 patients renamed.rds"))
AssignmentColors=c(TCells = "#8DB032", BCells = "#F9BF31", Other = "#E47B1E", 
                   Myeloid = "#CD4D53", Endo = "#DC73B0", Oligo = "#966CBF", Pericytes = "#2960A1", 
                   Glioma = "#003F5C")

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

pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"piecharts fragment and assignment figure 5a.pdf")), height = 6, width =10,onefile = T)
ggpubr::ggarrange(Final[[1]]+NoLegend(),Final[[6]]+NoLegend(),Final[[2]]+NoLegend(),Final[[7]]+NoLegend(),Final[[3]]+NoLegend(),
                  Final[[8]]+NoLegend(),Final[[4]]+NoLegend(),Final[[9]]+NoLegend(),Final[[5]]+NoLegend(),Final[[10]],ncol = 2,nrow = 5,
                  legend.grob = mylegend,legend = "right",widths = c(1,1)) 
dev.off()





#Figure5c
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

#load myeloid clusters R object
SeuratObj=readRDS(paste0(RobjDirectory,"MyeloidClusters-7-20-21-patients renamed.rds"))
set.seed(12)
Idents(SeuratObj)=CellInfo$Cluster
sampled.cells <- sample(x = SeuratObj@active.ident, size = 20000, replace = F)
subset=subset(SeuratObj,cells=names(sampled.cells))

Counts <- as.matrix(GetAssayData(object = subset, assay = "RNA", slot = "counts"))

library(msigdbr)
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
title="All Myeloid Fragments"
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
  annotate("text",x = xmin+2, y = ymax-0.5,label = "IFNG-res",color="white",fontface="bold",size=9)+ 
  annotate("rect", xmin = xmax-4, xmax = xmax, ymin = ymax-1, ymax = ymax, fill= "black")  +
  annotate("text",x = xmax-2, y = ymax-0.5,label = "Hypoxia",color="white",fontface="bold",size=9)+
  annotate("rect", xmin = xmin, xmax = xmin+4, ymin = ymin+1, ymax = ymin, fill= "black")  + 
  annotate("text",x = xmin+2, y = ymin+0.5,label = "Ox-Phos",color="white",fontface="bold",size=9)+ 
  annotate("rect", xmin = xmax-4, xmax =xmax, ymin = ymin+1, ymax = ymin, fill= "black")  +
  annotate("text",x = xmax-2, y = ymin+0.5,label = "TNFA-sig",color="white",fontface="bold",size=9) 
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
pdf(file =paste0(OutputDirectory,ObjName,Subset," pathway butterfly by Fragment figure 5b.pdf"), height = 15, width =35,onefile = T)
grid.arrange(grobs=Final, widths = c(4,1,1,1,1,1,1,1,1,1,1,1),layout_matrix = rbind(c(1, 2, 3,4,5,NA,NA,NA,NA,NA,NA,NA),c(1,6,7,8,9,10, 11, 12,13,14,15,16),
                                                                                c(1,17,18,19, 20, 21,22,23,24,25,26,27),
                                                                                c(1, 28,29, 30,31,32,33,34,35,36,37,38),
                                                                                c(1, 39,40,41,42,43,44,45,NA,NA,NA,NA))) 
dev.off()


##please insert code for cellphonedb