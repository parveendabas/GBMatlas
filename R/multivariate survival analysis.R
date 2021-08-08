library(Seurat)
library(scrabble)
library(dplyr)
##to get figure 2e
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
#sc=read.csv(paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"cluster CGGAscores.csv")),row.names = 1)
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
write.csv(Summtable,paste(OutputDirectory,ObjName,Subset,"res",resolution,"summary_multivariat_results.csv",sep=""))

pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"new cluster survival plots CGGA allFigure2E.pdf")), height = 3.5, width =16)
arrange_ggsurvplots(Final,print = TRUE, ncol = 5,nrow = 1,surv.plot.height = NULL,risk.table.height = NULL,ncensor.plot.height = NULL)
dev.off()
### to get figure 2f- I subset the matrix and only get the scores for GBM patients
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
write.csv(Summtable,paste(OutputDirectory,ObjName,Subset,"res",resolution,"GBM onlysummary_multivariat_results.csv",sep=""))

pdf(file =paste(paste0(OutputDirectory,ObjName,Subset,"res",resolution,"new cluster GBM only survivalFigure2F.pdf")), height = 3.5, width =16)
arrange_ggsurvplots(Final,print = TRUE, ncol = 5,nrow = 1,surv.plot.height = NULL,risk.table.height = NULL,ncensor.plot.height = NULL)
dev.off()



#################333
CGGA.EXP <- read.table("./CGGA.mRNAseq_325.RSEM-genes.20200506.txt",as.is=T,header = T,row.names=1)
AnnotationCGGA <- read.csv("./2020-08-20_CGGA_pheno.csv",as.is=T,header = T,row.names=1)
CGGA.GBM=filter(AnnotationCGGA,Histology=="GBM")
#got markers
Myeloid.Markers.surv= read.csv("Myeloidmarkers.csv")

#remove genes from signature gene lists that do not exist in expression matrix ##otherwise the code won't work
Counts=as.data.frame(CGGA.EXP)
geneset.all <-   Myeloid.Markers.surv[ Myeloid.Markers.surv$gene %in% rownames(Counts),]
geneset <- geneset.all[!duplicated(geneset.all$gene),]
genesetlist=geneset.all%>% split(x =.$gene, f = .$cluster,drop = F)
Scounts=Counts[geneset$gene,]
Scounts = log2(Scounts+1)
#get scores 
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
#write.csv(sc,"CGGAscores.csv")
#scores will slightly chage everytime you run them because of random sampling of control genes
sc=read.csv("CGGAscores.csv",row.names = 1)
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
for( i in 1: length(data.subset))
{
  YY= data.subset[i]
  Expression.Level=paste(YY,"Levels",sep=".")
  CGGA.sc2 <- CGGA.sc2%>% mutate(Expression.Level = ifelse(CGGA.sc2[YY]>=0, "Positive", "Negative"))
  CGGA.sc2$Expression.Level <-factor(CGGA.sc2$Expression.Level)
  fit1 <- survfit(surv_object ~Expression.Level, data = CGGA.sc2)
  #next three lines sdded by Joshy
  fit1.mv <- coxph(surv_object  ~Expression.Level + Gender + Histology +Subtype , data = CGGA.sc2)
  mv.res <- summary(fit1.mv)$coefficients
  mv.fname <- paste(data.subset[i],"_multivariat_results.csv",sep="")
  write.csv(mv.res,mv.fname)
  XX <- ggsurvplot(fit1, data = CGGA.sc2, pval = TRUE,pval.cExpression.Levelrd = c(50, 1),pval.size=6,legend="top",palette = c("Blue","Red"),risk.table = F,font.main = c(12, "bold"),legend.title = "status",
                   title=YY, legend.labs = c("Not-Enriched", "Enriched"),font.x = c(14, "bold"),font.y = c(14, "bold"),font.tickslab = c(12, "plain"),xlab="Time (Months)")
  Final[[i]] = XX
}
pdf(file ="Figure2E.pdf", height = 3.5, width =14,onefile = T)
arrange_ggsurvplots(Final,print = TRUE, ncol = 4,nrow = 1,surv.plot.height = NULL,risk.table.height = NULL,ncensor.plot.height = NULL)
dev.off()
### to get figure 2f- I subset the matrix and only get the scores for GBM patients
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
#write.csv(sc,"GBMonlyCGGAscores.csv")
#Since everytime you run the scoring the results are slightly different I provided the scoring file that reproduces the exact figures in the paper
sc=read.csv("GBMonlyCGGAscores.csv",row.names = 1)
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
for( i in 1: length(data.subset))
{
  YY= data.subset[i]
  Expression.Level=paste(YY,"Levels",sep=".")
  CGGA.sc2 <- CGGA.sc2%>% mutate(Expression.Level = ifelse(CGGA.sc2[YY]>=0, "Positive", "Negative"))
  CGGA.sc2$Expression.Level <-factor(CGGA.sc2$Expression.Level)
  fit1 <- survfit(surv_object ~Expression.Level, data = CGGA.sc2)
  XX <- ggsurvplot(fit1, data = CGGA.sc2, pval = TRUE,pval.cExpression.Levelrd = c(50, 1),pval.size=6,legend="top",palette = c("Blue","Red"),risk.table = F,font.main = c(12, "bold"),legend.title = "status",
                   title=YY, legend.labs = c("Not-Enriched", "Enriched"),font.x = c(14, "bold"),font.y = c(14, "bold"),font.tickslab = c(12, "plain"),xlab="Time (Months)")
  Final[[i]] = XX
}
pdf(file ="Figure2F.pdf", height = 3.5, width =14,onefile = T)
arrange_ggsurvplots(Final,print = TRUE, ncol = 4,nrow = 1,surv.plot.height = NULL,risk.table.height = NULL,ncensor.plot.height = NULL)
dev.off()


library(Seurat)
library(scrabble)
library(dplyr)
##to get figure 2e
#Load Dataset
CGGA.EXP <- read.table("./CGGA.mRNAseq_325.RSEM-genes.20200506.txt",as.is=T,header = T,row.names=1)
AnnotationCGGA <- read.csv("./2020-08-20_CGGA_pheno.csv",as.is=T,header = T,row.names=1)
CGGA.GBM=filter(AnnotationCGGA,Histology=="GBM")
#got markers
Myeloid.Markers.surv= read.csv("Myeloidmarkers.csv")

#remove genes from signature gene lists that do not exist in expression matrix ##otherwise the code won't work
Counts=as.data.frame(CGGA.EXP)
geneset.all <-   Myeloid.Markers.surv[ Myeloid.Markers.surv$gene %in% rownames(Counts),]
geneset <- geneset.all[!duplicated(geneset.all$gene),]
genesetlist=geneset.all%>% split(x =.$gene, f = .$cluster,drop = F)
Scounts=Counts[geneset$gene,]
Scounts = log2(Scounts+1)
#get scores 
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
#write.csv(sc,"CGGAscores.csv")
#scores will slightly chage everytime you run them because of random sampling of control genes
sc=read.csv("CGGAscores.csv",row.names = 1)
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
for( i in 1: length(data.subset))
{
  YY= data.subset[i]
  Expression.Level=paste(YY,"Levels",sep=".")
  CGGA.sc2 <- CGGA.sc2%>% mutate(Expression.Level = ifelse(CGGA.sc2[YY]>=0, "Positive", "Negative"))
  CGGA.sc2$Expression.Level <-factor(CGGA.sc2$Expression.Level)
  fit1 <- survfit(surv_object ~Expression.Level, data = CGGA.sc2)
  #next three lines sdded by Joshy
  fit1.mv <- coxph(surv_object  ~Expression.Level + Gender + Histology +Subtype , data = CGGA.sc2)
  mv.res <- summary(fit1.mv)$coefficients
  mv.fname <- paste(data.subset[i],"_multivariat_results.csv",sep="")
  write.csv(mv.res,mv.fname)
  XX <- ggsurvplot(fit1, data = CGGA.sc2, pval = TRUE,pval.cExpression.Levelrd = c(50, 1),pval.size=6,legend="top",palette = c("Blue","Red"),risk.table = F,font.main = c(12, "bold"),legend.title = "status",
                   title=YY, legend.labs = c("Not-Enriched", "Enriched"),font.x = c(14, "bold"),font.y = c(14, "bold"),font.tickslab = c(12, "plain"),xlab="Time (Months)")
  Final[[i]] = XX
}
pdf(file ="Figure2E.pdf", height = 3.5, width =14,onefile = T)
arrange_ggsurvplots(Final,print = TRUE, ncol = 4,nrow = 1,surv.plot.height = NULL,risk.table.height = NULL,ncensor.plot.height = NULL)
dev.off()
### to get figure 2f- I subset the matrix and only get the scores for GBM patients
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
#write.csv(sc,"GBMonlyCGGAscores.csv")
#Since everytime you run the scoring the results are slightly different I provided the scoring file that reproduces the exact figures in the paper
sc=read.csv("GBMonlyCGGAscores.csv",row.names = 1)
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
for( i in 1: length(data.subset))
{
  YY= data.subset[i]
  Expression.Level=paste(YY,"Levels",sep=".")
  CGGA.sc2 <- CGGA.sc2%>% mutate(Expression.Level = ifelse(CGGA.sc2[YY]>=0, "Positive", "Negative"))
  CGGA.sc2$Expression.Level <-factor(CGGA.sc2$Expression.Level)
  fit1 <- survfit(surv_object ~Expression.Level, data = CGGA.sc2)
  XX <- ggsurvplot(fit1, data = CGGA.sc2, pval = TRUE,pval.cExpression.Levelrd = c(50, 1),pval.size=6,legend="top",palette = c("Blue","Red"),risk.table = F,font.main = c(12, "bold"),legend.title = "status",
                   title=YY, legend.labs = c("Not-Enriched", "Enriched"),font.x = c(14, "bold"),font.y = c(14, "bold"),font.tickslab = c(12, "plain"),xlab="Time (Months)")
  Final[[i]] = XX
}
pdf(file ="Figure2F.pdf", height = 3.5, width =14,onefile = T)
arrange_ggsurvplots(Final,print = TRUE, ncol = 4,nrow = 1,surv.plot.height = NULL,risk.table.height = NULL,ncensor.plot.height = NULL)
dev.off()


