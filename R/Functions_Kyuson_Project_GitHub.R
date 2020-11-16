
Gene.Labels.pheatmap <- function(pheatmap, kept.labels, repel.degree) {
  
  # repel.degree = number within [0, 1], which controls how much 
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis
  
  heatmap <- pheatmap$gtable
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  
  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels, new.label$label, "")
  
  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant
    
    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), "npc"))
  }
  new.y.positions <- repelled.y(new.label$y, d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)
  
  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions
  
  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap, grobs = new.flag, t = 4, l = 4
  )
  
  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  
  # plot result
  #grid.newpage()
  #grid.draw(heatmap)
  
  # return a copy of the heatmap invisibly
  invisible(heatmap)
}


Voilin_Plot_OneGenePerLine <- function(SeuratObject, mat, GeneCol, GroupCol, colPalette, TotalPlotsPerRow){
  
  print(paste0("Processing Violin plots"))
  #SeuratObject=SCdata
  #mat=PlotGenes
  #GeneCol="Genes"
  #GroupCol="CT"
  #colPalette=cbPalette.Cluster
  #print(paste0("SeuratObject:",SeuratObject))
  print(SeuratObject)
  #print(corner(mat))
  print(paste0("GeneCol:",GeneCol))
  print(paste0("GroupCol:",GroupCol))
  print(paste0("colPalette:",colPalette))
  print(paste0("TotalPlotsPerRow:",TotalPlotsPerRow))
  
  mat.object=list(); i = 1;
  for(Ctype in unique(mat[,GeneCol])){
    #Ctype="Epithelial"
    print(paste0("Processing ViolinPlot ",GeneCol," ---> ",Ctype))
    GeneToPlot <- as.character(mat[mat[,GeneCol] ==Ctype,GeneCol]); GeneToPlot; length(GeneToPlot)
    GPlist=list(); j = 1;
    for(GP in GeneToPlot){
      print(paste0(GeneCol,": ",GP))
      GPlist[[j]] <- VlnPlot(SeuratObject, features = GP, pt.size = 0.00, combine = TRUE, cols=colPalette) + 
        theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), 
              legend.position = "none", axis.text.y = element_blank()) + 
        stat_summary(fun=median, geom="point", size=1, color="black")  + labs(title = NULL)
      j=j+1
    }
    mat.object[[i]] <- plot_grid(plotlist=GPlist, nrow = 1, ncol = TotalPlotsPerRow)
    i = i +1;
  }
  return(mat.object)
  
}



