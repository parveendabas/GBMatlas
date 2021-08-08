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
CellInfo2=CellInfo
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

#change fragment names
CellInfo$IDfrag=CellInfo$Fragment
id= levels(as.factor(CellInfo$ID))
for (j in 1:length(id)){
  CellInfo$Frag[CellInfo$ID==id[j]]=str_replace(CellInfo$Fragment[CellInfo$ID==id[j]],id[j],CellInfo$Patient[CellInfo$ID==id[j]])
}
CellInfo$Fragment=CellInfo$Frag
SeuratObj@meta.data=CellInfo
CellInfo=SeuratObj@meta.data
CellInfo$IDfrag[which(str_detect(row.names(CellInfo), "57"))] <- "CNSTM-394-3"

CellInfo$Fragment[which(str_detect(row.names(CellInfo), "57"))] <- "ndGBM-03-3"
SeuratObj@meta.data=CellInfo

#set.seed(001)
AssignmentColors= (Nour_pal("all",reverse = T)(length(levels(as.factor(SeuratObj@meta.data$Assignment)))))
#names(AssignmentColors)= levels(as.factor(SeuratObj@meta.data$Assignment))
names(AssignmentColors)=c("TCells", "BCells","Other" ,  "Myeloid",  "Endo", "Oligo",  "Pericytes", "Glioma")
                          
Final=list()
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
#pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"piecharts fragment and assignment.pdf")), height = 6, width =10,onefile = T)

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

pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"piecharts type,sexandgrade and assignment.pdf")), height = 4, width =6,onefile = T)
ggpubr::ggarrange(x1+NoLegend(),x2+NoLegend(),x3+NoLegend(),ncol =1,nrow = 3,
                  legend.grob = mylegend,legend = "right") 
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
for(i in 1:length(pat2)){
  x=ggplot(SeuratObj@meta.data%>%filter(Patient==pat2[i]), aes(Patient, fill=Assignment))+geom_bar(stat="count",position="fill",colour = "white",width =0.7)+ 
    scale_fill_manual(values = AssignmentColors)+coord_polar("y",start=0) + facet_grid(~ Fragment) + #labs(title=pat2[i])+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_blank(),axis.ticks = element_blank(),
          panel.background = element_blank(), axis.title = element_blank(),strip.background=element_rect(fill="white"),
          strip.text=element_text(face = "bold",color = "black",size=8),plot.margin = unit(c(0,-0.2,-0.2,0),'cm'))
  Final[[i]]=x
}
mylegend<-get_legend(Final[[2]]+theme_min()+ NoAxes()+guides(colour = guide_colourbar(title=Y,title.position = "top",title.theme = element_text(size=10))))

pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"piecharts fragment and Assignment.pdf")), height = 6, width =10,onefile = T)
ggpubr::ggarrange(Final[[1]]+NoLegend(),Final[[6]]+NoLegend(),Final[[2]]+NoLegend(),Final[[7]]+NoLegend(),Final[[3]]+NoLegend(),
                  Final[[8]]+NoLegend(),Final[[4]]+NoLegend(),Final[[9]]+NoLegend(),Final[[5]]+NoLegend(),Final[[10]],ncol = 2,nrow = 5,
                  legend.grob = mylegend,legend = "right",widths = c(1,1)) 
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

pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"piecharts fragment and Assignment.pdf")), height = 6, width =10,onefile = T)
grid.arrange(grobs=Final, widths = c(1,1,1,1,1),layout_matrix = rbind(c(1, 2, 3,4,NA),c(5, 6, 7,8,NA),c(9,10,11,12,13),c(14,15,16,17,18))) 


ggpubr::ggarrange(Final[[1]]+NoLegend(),Final[[2]]+NoLegend(),Final[[3]]+NoLegend(),Final[[4]]+NoLegend(),NA,
                  Final[[5]]+NoLegend(),Final[[6]]+NoLegend(),Final[[7]]+NoLegend(),Final[[8]]+NoLegend(),NA,
                  Final[[9]]+NoLegend(),Final[[10]]+NoLegend(),Final[[11]]+NoLegend(),Final[[12]]+NoLegend(),Final[[13]]+NoLegend(),
                  Final[[14]]+NoLegend(),Final[[15]]+NoLegend(),Final[[16]]+NoLegend(),Final[[17]]+NoLegend(),Final[[18]]+NoLegend(),
                  ncol = 5,nrow = 4,
                  legend.grob = mylegend,legend = "right",widths = c(1,1,1,1,1)) 
dev.off()