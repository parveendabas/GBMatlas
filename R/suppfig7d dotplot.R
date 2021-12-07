rm(list=ls()) # clear workspace

### cellphonedb plot dot_plot --rows in/rows.txt --columns in/columns.txt

library(ggplot2)
  #setwd()
  

pkWD <- "/Users/kumarpa/Desktop/Work/Jax/Kyuson/Upload_To_Github"
PatientDir <- paste0("/Users/kumarpa/Desktop/Work/Jax/Kyuson/CellCellComm/cpdb/Output_from_cpdb_ALL_Patients"); PatientDir

setwd(pkWD)
plotWD.main <- paste(getwd(),paste0("DotPlot_CellphoneDB"),sep="/"); print(plotWD.main)
dir.create(file.path(getwd(),paste0("DotPlot_CellphoneDB")), showWarnings = FALSE)


CTlist=c("Glioma|Myeloid", "Glioma|Tcells",  "Myeloid|Glioma", "Myeloid|Tcells",  "Tcells|Glioma", "Tcells|Myeloid")
Log2MeanCutoff=2
Patient <- c("Output_cpdb_")
CellCount="Analysis_cpdb_ALLcells_per_CT"
pvalcutoff=0.01

for(CTuse in CTlist){
  #CTuse="Glioma|Myeloid"
  
  sampleNames=dir(path = PatientDir, pattern = paste0(Patient,".*$")); sampleNames
  print(paste0("Processig for pval-cutoff:", pvalcutoff))
      
  pval.list=list(); mean.list=list(); mean.list.ORI=list(); pval.list.ORI=list(); 
  for(sample in sampleNames){
  #sample="Output_cpdb_ndGBM-09"
  
    print(paste0("Processing sample: ",sample))
    pkWD <- paste0(PatientDir,"/",sample,"/",CellCount); pkWD
      
  setwd(pkWD)
  means_path = './means.txt'
  pvalues_path = './pvalues.txt'
  significant_means_path='./significant_means.txt'
  means_separator = '\t'
  pvalues_separator = '\t'
  
  ## Process the files
  all_pval.Full = read.table(pvalues_path, header=T, stringsAsFactors = F, sep=means_separator, comment.char = '', check.names=F); print(dim(all_pval.Full))
  all_pval.Full <- all_pval.Full[!duplicated(all_pval.Full$interacting_pair),]
  rownames(all_pval.Full) <- all_pval.Full$interacting_pair; head(all_pval.Full)
  all_pval <- all_pval.Full
  intr_pairs = all_pval$interacting_pair; length(intr_pairs)
  all_pval = all_pval[,-c(1:11)]
  selected_columns = colnames(all_pval)
  selected_columns = selected_columns[selected_columns %in% paste0(CTuse)]; selected_columns
  all_pval <- all_pval[,selected_columns, drop=FALSE]; head(all_pval)
  sig_pval_intr_pairs <- rownames(all_pval[rowSums(all_pval<=pvalcutoff)>0,,drop=FALSE]); length(sig_pval_intr_pairs)
  
  all_means = read.table(means_path, header=T, stringsAsFactors = F, sep=pvalues_separator, comment.char = '', check.names=F); print(dim(all_means))
  all_means <- all_means[!duplicated(all_means$interacting_pair),]
  rownames(all_means) <- all_means$interacting_pair; head(all_means)
  all_means = all_means[,-c(1:11)]
  all_means <- all_means[,selected_columns,drop=FALSE]; head(all_means)
  all_means.Full <- all_means
  
  
  if(length(sig_pval_intr_pairs) > 0){
    
  sig_means = read.table(means_path, header=T, stringsAsFactors = F, sep=pvalues_separator, comment.char = '', check.names=F); print(dim(sig_means))
  sig_means <- sig_means[!duplicated(sig_means$interacting_pair),]
  rownames(sig_means) <- sig_means$interacting_pair; head(sig_means)
  sig_means = sig_means[,-c(1:11)]; head(sig_means)
  sig_means<- sig_means[,selected_columns,drop=FALSE]; head(sig_means)
  sig_means$Log2Mean <- log2(sig_means[,CTuse]); 
  sig_means <- sig_means[sig_means$Log2Mean > Log2MeanCutoff,]; head(sig_means); dim(sig_means)
  sig_means$Log2Mean <- NULL
  sig_means$sum <- rowSums(sig_means); head(sig_means)
  
  discarded_rows <- rownames(sig_means[sig_means$sum<=0,]); head(discarded_rows); length(discarded_rows)
  Temp_selected_rows <- rownames(sig_means[sig_means$sum>0,]); head(Temp_selected_rows); length(Temp_selected_rows)
  selected_rows_pval <- intersect(Temp_selected_rows, sig_pval_intr_pairs); length(selected_rows_pval)
  
  colnames(all_pval) <- paste0(colnames(all_pval),"|", gsub("Output_cpdb_","", sample))
  all_pval <- all_pval[selected_rows_pval,,drop=FALSE]
  all_pval$LR <- rownames(all_pval)
  pval.list[[sample]] <- all_pval
  
  colnames(all_means) <- paste0(colnames(all_means),"|", gsub("Output_cpdb_","", sample))
  all_means$LR <- rownames(all_means)
  mean.list[[sample]] <- all_means
  
  
  
  
  
  } else {
    print(paste0("No interactions found in sample: ",sample, " for interaction: ",CTuse))
    print(paste0("*************************"))
    Sys.sleep(0.5)
    print("Using empty data.frame")
    all_pval <- data.frame(A=as.numeric(),B=as.character()); colnames(all_pval) <- c(paste0(CTuse,"|",gsub("Output_cpdb_","", sample)), "LR")
    pval.list[[sample]] <- all_pval
    all_means <- all_pval
    mean.list[[sample]] <- all_means
    
    
  }
  
  
  if(length(selected_columns) > 0){
  print("Filling Full data.frames")
  all_pval.Full <- all_pval.Full[,selected_columns,drop=FALSE]
  colnames(all_pval.Full) <- paste0(colnames(all_pval.Full),"|", gsub("Output_cpdb_","", sample))
  all_pval.Full$LR <- rownames(all_pval.Full); head(all_pval.Full)
  pval.list.ORI[[sample]] <- all_pval.Full
  
  all_means.Full <- all_means.Full[,selected_columns,drop=FALSE]
  colnames(all_means.Full) <- paste0(colnames(all_means.Full),"|", gsub("Output_cpdb_","", sample))
  all_means.Full$LR <- rownames(all_means.Full)
  mean.list.ORI[[sample]] <- all_means.Full; head(all_means.Full)
  print("Done Filling Full data.frames")
  
  } else {
    
    print(paste0("No Columns itself found in sample: ",sample, " for interaction: ",CTuse))
    Sys.sleep(0.3)
    print("Using empty data.frame for Full")
    all_pval.Full <- data.frame(A=as.numeric(),B=as.character()); colnames(all_pval.Full) <- c(paste0(CTuse,"|",gsub("Output_cpdb_","", sample)), "LR")
    pval.list.ORI[[sample]] <- all_pval.Full
    all_means.Full <- all_pval.Full
    mean.list.ORI[[sample]] <- all_means.Full
  }
    
  }
  
  ## Merge data frames
  print(paste0("Merging ",length(sampleNames), " data frames"))
  CombinedLRs <- Reduce(function(...) merge(..., by = "LR", all=TRUE), pval.list)
  all_pval_patient <- Reduce(function(...) merge(..., by = "LR", all=TRUE), pval.list.ORI)
  rownames(all_pval_patient) <- all_pval_patient$LR
  all_pval_patient <- all_pval_patient[CombinedLRs$LR,]
  all_mean_patient <- Reduce(function(...) merge(..., by = "LR", all=TRUE), mean.list.ORI)
  rownames(all_mean_patient) <- all_mean_patient$LR
  all_mean_patient <- all_mean_patient[CombinedLRs$LR,]
  
  rownames(all_pval_patient) <- all_pval_patient$LR; all_pval_patient$LR <- NULL; head(all_pval_patient,1)
  selected_rows_pval <- rownames(all_pval_patient)
  head(all_pval_patient,1); dim(all_pval_patient)
  
  rownames(all_mean_patient) <- all_mean_patient$LR; all_mean_patient$LR <- NULL; head(all_mean_patient,1)
  selected_rows_mean <- rownames(all_mean_patient)
  head(all_mean_patient,1); dim(all_mean_patient)
  
  ColsOrder <- sort(colnames(all_pval_patient))
  all_pval_patient <- all_pval_patient[,ColsOrder]
  all_mean_patient <- all_mean_patient[,ColsOrder]
  
  
  PlotPval="YES"
  if(PlotPval=="YES"){
    
    setwd(plotWD.main)
    pdf(file=paste0("Plotting_",Patient,"ALL_",CellCount,"_pval_based_",pvalcutoff,"_",CTuse,"_Log2mean_above_",Log2MeanCutoff,"_Total_",length(selected_rows_pval),"_Pairs.pdf"),height = 13,width = 17)
    
    
    set <- split(selected_rows_pval, ceiling(seq_along(selected_rows_pval)/60)); length(set)
    selected_columns=colnames(all_pval_patient)
    intr_pairs=rownames(all_pval_patient)
    
    for(i in 1:length(set)){
      #i=1
      print(paste0("Set no. ",i," of total ",length(set)," sets"))
      
      Plot_rows <- c(set[[i]]); head(Plot_rows); length(Plot_rows)
      
      
      all_mean_patient.temp=all_mean_patient[intr_pairs,]
      all(rownames(all_pval_patient)==rownames(all_mean_patient.temp))
      
      sel_pval = all_pval_patient[match(Plot_rows, intr_pairs), selected_columns]
      sel_means = all_mean_patient.temp[match(Plot_rows, intr_pairs), selected_columns]
      
      df_names = expand.grid(Plot_rows, selected_columns)
      pval = unlist(sel_pval)
      pval[pval==0] = 0.0009
      plot.data = cbind(df_names,pval)
      pr = unlist(as.data.frame(sel_means))
      pr[pr==0] = 1
      plot.data = cbind(plot.data,log2(pr))
      colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')
      print(dim(plot.data))
      head(plot.data)
      
      my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)
      
      p1 <- ggplot(plot.data,aes(x=clusters,y=pair)) +
              geom_point(aes(size=-log10(pvalue),color=mean)) +
              scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
              scale_x_discrete(position = "top") +
              theme_bw() +
              theme(panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    axis.text=element_text(size=15, colour = "black"),
                    axis.text.x = element_text(size=15, angle = 90, hjust = 1),
                    axis.text.y = element_text(size=15, colour = "black"),
                    axis.title=element_blank(),
                    panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
      
      print(p1)
      
    }
    
    dev.off()
    
     selected_rows_pval
      head(all_pval_patient)
      write.table(all_pval_patient[selected_rows_pval,selected_columns],file=paste0("Details_",Patient,"ALL_",CellCount,"_pval_based_",pvalcutoff,"_",CTuse,"_Log2mean_above_",Log2MeanCutoff,"_Total_",length(selected_rows_pval),"_Pairs.txt"),row.names=T,quote=F,sep="\t")
                                                                            
    
  }
  
  
## CTlist
}
