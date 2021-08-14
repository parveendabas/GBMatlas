library(easypackages)
MyPackages<-c("dplyr","ggplot2","ggpubr","gridExtra","viridis","egg",
              "grid","lattice","gtools","Biobase","RColorBrewer","tibble",
              "Seurat","cowplot","patchwork","stringr","ComplexHeatmap","SingleCellExperiment", "ggmin")
libraries(MyPackages)

#Figure 6A
Directory="~/Box/Yun lab projects/scRNAseq Project/Human Glioma/ALLSAMPLESJuly2021/"
RobjDirectory=paste0(Directory,"R_Objects/")
#load myeloid clusters R object
SeuratObj=readRDS(paste0(RobjDirectory,"MyeloidClusters-7-20-21-patients renamed.rds"))
VlnPlot(SeuratObj,group.by = "Cluster",features = "S100A4",cols=ClusterColors,pt.size = 0)+theme_min2()+NoLegend()+theme(title = element_blank())


#Figure 6B
#load T cell clusters R object
SeuratObj=readRDS(paste0(RobjDirectory,"TcellClusters-7-22-21-goodwithoutcluster6 new patient names.rds"))
VlnPlot(SeuratObj,group.by = "Cluster",features = "S100A4",cols=ClusterColors,pt.size = 0)+theme_min2()+NoLegend()+theme(title = element_blank())

#Figure 6D

library(survival)
library(survminer)
#load data
CGGA.EXP <- read.table("~/Downloads/CGGA.mRNAseq_325.RSEM-genes.20200506.txt",as.is=T,header = T,row.names=1)
AnnotationCGGA <- read.csv("~/Downloads/CGGA.mRNAseq_All_clinical.csv",as.is=T,header = T,row.names=1)
CGGA.GBM=filter(AnnotationCGGA,Histology=="GBM")
#extract S100A4 expression and merge with annotation
data <- CGGA.EXP["S100A4",]
t.data=t(data)
CGGA.dat=merge(t.data,AnnotationCGGA,by.x=0, by.y=0)
#get KM plots
CGGA.dat$OS<-as.numeric(CGGA.dat$survival)
CGGA.dat$status<-as.numeric(CGGA.dat$status)
surv_object <- Surv(time = CGGA.dat$OS, event = CGGA.dat$status)
hist(CGGA.dat$S100A4)
median(CGGA.dat$S100A4)
CGGA.dat <-CGGA.dat %>% mutate(S100A4.Levels = ifelse(S100A4 >=14.91, "High", "Low"))
table(CGGA.dat$S100A4.Levels)
CGGA.dat$S100A4.Levels <-factor(CGGA.dat$S100A4.Levels)
fit1 <- survfit(surv_object ~ S100A4.Levels, data = CGGA.dat)
fit1.mv <- coxph(surv_object  ~S100A4.Levels + Gender + Histology + Recurrence  + MGMT_status +IDH_mutation_status , data = CGGA.dat)
mv.res <- summary(fit1.mv)$coefficients
mv.fname <- paste("S100a4 All glioma_multivariat_results.csv",sep="")
write.csv(mv.res,mv.fname)

summary(fit1)
P1=ggsurvplot(fit1, data = CGGA.dat, pval = TRUE,pval.coord = c(50, 1),pval.size=6,legend="top",palette = c("Red","Blue"),risk.table = F,font.main = c(12, "bold"),legend.title = "Expression",
              title="S100A4-AllGliomas-CGGA", legend.labs = c("High", "Low"),font.x = c(14, "bold"),font.y = c(14, "bold"),font.tickslab = c(12, "plain"))
#get GBM patients survival only
Counts=as.data.frame(CGGA.EXP)
data <- Counts["S100A4",]
t.data=t(data)
CGGA.dat=merge(t.data,CGGA.GBM,by.x=0, by.y=0)
CGGA.dat$OS<-as.numeric(CGGA.dat$survival)
CGGA.dat$status<-as.numeric(CGGA.dat$status)
surv_object <- Surv(time = CGGA.dat$OS, event = CGGA.dat$status)
hist(CGGA.dat$S100A4)
median(CGGA.dat$S100A4)
CGGA.dat <-CGGA.dat %>% mutate(S100A4.Levels = ifelse(S100A4 >=38.39, "High", "Low"))
table(CGGA.dat$S100A4.Levels)
CGGA.dat$S100A4.Levels <-factor(CGGA.dat$S100A4.Levels)
fit1 <- survfit(surv_object ~ S100A4.Levels, data = CGGA.dat)
fit1.mv <- coxph(surv_object  ~S100A4.Levels + Gender +  Recurrence  + MGMT_status +IDH_mutation_status , data = CGGA.dat)
mv.res <- summary(fit1.mv)$coefficients
mv.fname <- paste("GBM only_multivariat_results.csv",sep="")
write.csv(mv.res,mv.fname)
summary(fit1)
P2=ggsurvplot(fit1, data = CGGA.dat, pval = TRUE,pval.coord = c(50, 1),pval.size=6,legend="top",palette = c("Red","Blue"),risk.table = F,font.main = c(12, "bold"),legend.title = "Expression",
              title="S100A4-GBM only-CGGA", legend.labs = c("High", "Low"),font.x = c(14, "bold"),font.y = c(14, "bold"),font.tickslab = c(12, "plain"))
P=list(P1,P2)
pdf(file ="Figure3e.pdf", height = 3.5, width =14,onefile = T)
arrange_ggsurvplots(P,print = TRUE, ncol = 4,nrow = 1,surv.plot.height = NULL,risk.table.height = NULL,ncensor.plot.height = NULL)
dev.off()
