library(presto)
setwd("~/Box/Yun lab projects/scRNAseq Project/Human Glioma/ALLSAMPLESJuly2021/Output/Glioma/pathway")

Idents(SeuratObj)=SeuratObj@meta.data$Assignment
CellInfo=SeuratObj@meta.data
objname="Glioma"
Cluster= levels(as.factor(CellInfo$Cluster))
Iterate= levels(as.factor(CellInfo$Assignment))
#Iterate=c("Myeloid","Glioma","Bcells",  "Microglia",  "NKCells", "Tcells")
load("~/Box/Yun lab projects/ESC Analysis October 2020/r objects/C8mousegenesconverted.RData")
library(presto)
library(msigdbr)
library(dplyr)
library(fgsea)
library(tibble)
library(pheatmap)
Z="AllClusters"
pathways = c("H")
names= c("Hallmark")
for (j in 1:length(Iterate)){
  Z=Iterate[j]
  print(Z)
  Annot.subset <-subset(SeuratObj, cells=WhichCells(SeuratObj,idents =Z ))
  Annot=Annot.subset@meta.data
  #2 Run Test
  Clusters.genes <- wilcoxauc(Annot.subset , 'Group')
  write.csv(Clusters.genes,paste(Z,"DE genes between Samples in", Z,"-presto output.csv"))
  #check if test was successful
  dplyr::count(Clusters.genes, group)
  pathways = c("H","C7","CP:KEGG","CP:BIOCARTA")
  names= c("Hallmark","C7","KEGG","Biocarta")
  category=c("H","C7","C2","C2")
  for (i in 3: length(pathways)) {
    cat=pathways[i]
    print(cat)
    name=names[i]
    
    if (category[i]=="C8") {
      fgsea_sets<- mousemygo
      Annot.pathway2=as.data.frame(names(mousemygo))
      names(Annot.pathway2)="pathway"
    } else {
      if (category[i]=="C2") {
        m_df_H<- msigdbr(species = "Homo sapiens", category = "C2",subcategory = cat)
      } else {
        m_df_H<- msigdbr(species = "Homo sapiens", category = cat)
      }
      Annot.pathway2=as.data.frame(levels(as.factor(m_df_H$gs_name)))
      names(Annot.pathway2)="pathway"
      fgsea_sets<- m_df_H %>% split(x = .$gene_symbol, f = .$gs_name)
    }
    Clusters=levels(as.factor(Annot$Cluster))
    
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
      write.csv(hmc,paste(name, X,Z,".csv"))
      names(Y)[names(Y)=="NES"]=paste(X)
      Annot.pathway<-Y[,c("pathway",paste(X))]
      Annot.pathway2<-merge(Annot.pathway, Annot.pathway2, by.x="pathway", by.y="pathway")
      
      plot=ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
        geom_segment( aes(reorder(pathway, NES), xend=pathway, y=0, yend=NES)) +
        geom_point( size=3, aes( fill = Enrichment), shape=21, stroke=1) +
        scale_fill_manual(values = c("Down-regulated" = "dodgerblue", "Up-regulated" = "firebrick") ) +
        scale_y_continuous(limits = NULL,expand = expansion(mult = 0.2, add = 0))+
        theme(axis.text= element_text(size=8, face="bold"))+theme_min()+
        coord_flip() +labs(x="Pathway", y="Normalized Enrichment Score",title=paste(name,"-",X,Z))
      ggsave(plot, file=paste(name,"-",objname, X,Z,".pdf"), scale=1,width=10)
      
    }
    ##make a heatmap
    rownames(Annot.pathway2)=Annot.pathway2$pathway
    Annot.pathway2=Annot.pathway2[,-1]
    positions=Clusters

    Annot.pathway2=Annot.pathway2[,positions]
    Annot.pathway2[is.na(Annot.pathway2)]=0
    Annot.pathway3=rbind(Annot.pathway2%>%filter_all(any_vars(.<= 2)),Annot.pathway2%>%filter_all(any_vars(.<= -2)))
    redblackblu=CustomPalette(low = "blue", high = "red", mid = "white", k = 1000)
    pheatmap(Annot.pathway3,scale = "none",cluster_rows = T,cluster_cols = F,border_color = "black",treeheight_row = 0,legend_breaks = c(-5,-2,0,2,5),
             fontsize =8, cellwidth = 15,cellheight =10,color= plasma(10),angle_col = 90,
             main = paste("                                          ",name, Z,objname," (NES)"),filename = paste("GSEA",name,Z," (NES)-Heatmap.pdf"))
    
    
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
             main = paste("                                          ",name,Z,objname," (NES)"),filename = paste("GSEA",name,Z," (NES) top and bottom 2each-Heatmap.pdf"))
    pathnames=c(row.names(Annot.pathway2 %>% top_n(n = 5, wt = eval(parse(text=names(Annot.pathway2)[1])))),row.names(Annot.pathway2 %>% top_n(n = -5, wt = eval(parse(text=names(Annot.pathway2)[1])))))
    for (i in 2: length(names(Annot.pathway2))){
      pathnames1=c(row.names(Annot.pathway2 %>% top_n(n = 5, wt = eval(parse(text=names(Annot.pathway2)[i])))),row.names(Annot.pathway2 %>% top_n(n = -5, wt = eval(parse(text=names(Annot.pathway2)[i])))))
      pathnames=c(pathnames,pathnames1)
    }
    pathnames =pathnames[!duplicated(pathnames)]
    Annot.pathway4=Annot.pathway2[pathnames,]
    Annot.pathway4[is.na(Annot.pathway4)]=0
    redblackblu=CustomPalette(low = "blue", high = "red", mid = "white", k = 1000)
    pheatmap(Annot.pathway4,scale = "none",cluster_rows = T,cluster_cols = F,border_color = "black",treeheight_row = 0,legend_breaks = c(-5,-2,0,2,5),
             fontsize =8, cellwidth = 15,cellheight =10,color= plasma(10),angle_col = 90,
             main = paste("                                          ",name,Z,objname," (NES)"),filename = paste("GSEA",name,Z," (NES) top 5 and bottom 5each-Heatmap.pdf"))
    
  } 
  
  
}

pathwaynames=read.csv("../pathways.csv",header = F,as.is = T)
Annot.pathway2=as.data.frame(read.csv("../pathways.csv",header = F))
names(Annot.pathway2)="pathway"
Clusters=levels(as.factor(MB@meta.data$Group))
library("tidyverse") 
objname="X2808"
for (i in 1:length(Clusters)){
  ZZ=Clusters[i]
  filess <- list.files(pattern=paste(ZZ,"AllClusters .csv"))
  filessL <- as.list(filess)
  names(filessL) <- str_replace(filess,paste(ZZ,"AllClusters .csv"),"")
  filessL <- lapply(filessL, read.csv, row.names=1)
  filessDF=filessL %>% reduce(full_join) 
  filessDFS=filessDF %>% filter(filessDF$pathway %in%Annot.pathway2$pathway, )
  assign(paste0(ZZ,"pathways"), filessDFS)
  names(filessDFS)[names(filessDFS)=="NES"]=paste0(ZZ)
  Annot.pathway<-filessDFS[,c("pathway",paste0(ZZ))]
  Annot.pathway2<-merge(Annot.pathway, Annot.pathway2, by.x="pathway", by.y="pathway")  
}
Annot.pathway2=Annot.pathway2[!duplicated(Annot.pathway2$pathway),]
rownames(Annot.pathway2)=Annot.pathway2$pathway
Annot.pathway2=Annot.pathway2[,-1]
positions=Clusters
Annot.pathway2=Annot.pathway2[,positions]
Annot.pathway2[is.na(Annot.pathway2)]=0
write.csv(Annot.pathway2,"selectpathways.csv")

#Annot.pathway2=read.csv("selectpathways.csv",row.names = 1)
#Annot.pathway3=rbind(Annot.pathway2%>%filter_all(any_vars(.<= 2)),Annot.pathway2%>%filter_all(any_vars(.<= -2)))
redblackblu=CustomPalette(low = "blue", high = "red", mid = "white", k = 1000)
pheatmap(Annot.pathway2,scale = "none",cluster_rows = T,cluster_cols = F,border_color = "black",
         fontsize =8, cellwidth = 15,cellheight =10,color= plasma(10),angle_col = 90,
         main = paste("select pathways All Clusters (NES)"),filename = paste(objname, "GSEA select pathways allclusters (NES)-Heatmap.pdf"))
n=4
Annot.pathway3=filter(Annot.pathway2, across(, ~ .x >= n|.x <= -n))

pheatmap(Annot.pathway3,scale = "none",cluster_rows = T,cluster_cols = F,border_color = "black",
         fontsize =8, cellwidth = 15,cellheight =10,color= plasma(10),angle_col = 90,
         main = paste("select pathways All Clusters (NES)"),filename = paste0(objname, " GSEA select pathways allclusters (NES>",n,")-Heatmap.pdf"))

n=2
pathnames=c(row.names(Annot.pathway2 %>% top_n(n = n, wt = eval(parse(text=names(Annot.pathway2)[1])))),row.names(Annot.pathway2 %>% top_n(n = -n, wt = eval(parse(text=names(Annot.pathway2)[1])))))
for (i in 2: length(names(Annot.pathway2))){
  pathnames1=c(row.names(Annot.pathway2 %>% top_n(n = n, wt = eval(parse(text=names(Annot.pathway2)[i])))),row.names(Annot.pathway2 %>% top_n(n = -n, wt = eval(parse(text=names(Annot.pathway2)[i])))))
  pathnames=c(pathnames,pathnames1)
}
pathnames =pathnames[!duplicated(pathnames)]
Annot.pathway4=Annot.pathway2[pathnames,]
Annot.pathway4[is.na(Annot.pathway4)]=0
redblackblu=CustomPalette(low = "blue", high = "red", mid = "white", k = 1000)
pheatmap(Annot.pathway4,scale = "none",cluster_rows = T,cluster_cols = F,border_color = "black",
         fontsize =8, cellwidth = 15,cellheight =10,color= plasma(10),angle_col = 90,
         main = paste("select pathways All Clusters (NES)"),filename = paste(objname, "GSEA select top and bottom", n , "pathways all clusters (NES)-Heatmap.pdf"))


####
#Load Libraries
library(Seurat)
#install.packages("remotes")
#remotes::install_github("immunogenomics/presto")
library(presto)
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
library(tibble)
library(pheatmap)

# Load Your Seurat Object--> this is the same one from Parveen's Seurat objects in box
#SeuratObj <-readRDS("~/Downloads/Seurat_Formatted_Normalized_USE_GBM_Cluster_08_09_10_11.rds")
#SeuratObj=readRDS(paste0(RobjDirectory,"GliomaClusters-7-15-21.rds"))
#CellInfo <- SeuratObj@meta.data

#2 Run Test
SeuratObjClusters.genes <- wilcoxauc(SeuratObj , 'Cluster')
write.csv(SeuratObjClusters.genes,paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"wilcoxauc.csv")))
dplyr::count(SeuratObjClusters.genes, group)

# Load Hallmark Pathways
m_df_H<- msigdbr(species = "Homo sapiens", category = "H")

# Generate Supp Fig 4A- I made each one individually and combined them in Powepoint
Clusters=levels(as.factor(CellInfo$Cluster))
Glioma.Hallmark2=as.data.frame(levels(as.factor(m_df_H$gs_name)))
names(Glioma.Hallmark2)="pathway"
for (i in 1: length(Clusters)){
  X=  Clusters[i]
  Genes<-SeuratObjClusters.genes %>%
    dplyr::filter(group == X) %>%
    arrange(desc(auc)) %>% 
    dplyr::select(feature, auc)
  ranks<- deframe(Genes)
  head(ranks)
  fgsea_sets<- m_df_H %>% split(x = .$gene_symbol, f = .$gs_name)
  fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  fgseaResTidy %>% 
    dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
    arrange(padj) %>% 
    head()
  fgseaResTidy$Enrichment = ifelse(fgseaResTidy$NES > 0, "Up-regulated", "Down-regulated")
  Filtidy<-fgseaResTidy %>% filter(padj < 0.05) 
  filtRes = rbind(Filtidy %>% filter(Enrichment == "Up-regulated") %>% head(n= 10),
                  Filtidy %>% filter(Enrichment == "Down-regulated") %>% head(n= 10))
  
  plot=ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill= Enrichment)) +
    scale_fill_manual(values = c("Down-regulated" = "blue", "Up-regulated" = "red") ) +
    coord_flip() +
    labs(x=NULL, y="Normalized Enrichment Score",
         title=paste("GSEA-Hallmark-Human Glioma", X))+ 
    theme_minimal()
  ggsave(plot, file=paste("Hallmark-Human Glioma", X,".pdf"), scale=1)
  
  Y<-fgseaResTidy ## this part is to make the heatmap later
  hmc=as.data.frame(Y)
  hmc=apply(Y,2,as.character) 
  #write.csv(hmc,paste("Hallmark-Human Glioma", X,".csv"))
  names(Y)[names(Y)=="NES"]=paste(X)
  Glioma.Hallmark<-Y[,c("pathway",paste(X))]
  Glioma.Hallmark2<-merge(Glioma.Hallmark, Glioma.Hallmark2, by.x="pathway", by.y="pathway")
  
}