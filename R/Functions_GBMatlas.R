
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




.prepare_mat <- function(mat, FUN = rowMeans) {
  # If x is a matrix (rows by cols), transform x to row or column means
  if (!is.null(dim(mat))) x = FUN(mat)
  # If x is numeric, it will first be sorted.
  if (is.numeric(x)) x = sort(x)
  return(x)
}


.match_bins <- function(Group, bins) {
  .check_arg(arg = Group)
  .check_missing(Group = Group, ref = bins)

  binID = bins[Group]
  binmatches = sapply(binID, function(id) names(bins)[bins == id], simplify = F)
  binmatches
}


.sample_bins <- function(bins,
                         n = 100,
                         replace = FALSE) {

  .check_sample_size(bins = bins, n = n, replace = replace)
  sapply(bins, function(bin) {
    sample(bin, size = n, replace = replace) },
    simplify = F) %>%
    unlist(use.names = F)
}



binmatch <- function(Group,
                     mat = NULL,
                     x = NULL,
                     bins = NULL,
                     nbin = 30,
                     n = 100,
                     replace = FALSE,
                     ...) {
  .check_args_exist(x = x, data = mat, bins = bins)

  if (!is.null(x) | !is.null(mat)) {
    bins = bin(x = x, mat = mat, breaks = nbin, ...)
  }

  binmatches = .match_bins(Group = Group, bins = bins)
  binsamples = .sample_bins(bins = binmatches, n = n, replace = replace)
  binsamples
}



bin <- function(mat = NULL,
                x = NULL,
                breaks = 30,
                FUN = rowMeans) {
  if (!is.null(mat)) x = .prepare_mat(mat, FUN = FUN)
  .check_arg(arg = x)
  binIDs = cut(seq_along(x), breaks = breaks, labels = F, include.lowest = T)
  .name_binIDs(binIDs = binIDs, x = x)
}



hierarchy = function(m, quadrants = NULL, log.scale = T) {

  if (!is.null(quadrants)) {
    stopifnot(all(unlist(quadrants) %in% colnames(m)))
    dat = as.data.frame(sapply(quadrants, function(col) do.call(pmax, list(as.data.frame(m[, col])))))
  } else {
    stopifnot(ncol(m) == 4)
    dat = as.data.frame(m)
  }

  rows = rownames(m)
  colnames(dat) = c('bl', 'br', 'tl', 'tr')

  dat = dat %>%
    dplyr::mutate(bottom = pmax(bl, br),
                  top = pmax(tl, tr),
                  b.center = br - bl,
                  t.center = tr - tl,
                  x = ifelse(bottom > top, b.center, t.center), # dependent var
                  x.scaled = (sign(x) * log2(abs(x) + 1)),
                  y = top - bottom, # independent var
                  y.scaled = (sign(y) * log2(abs(y) + 1)))

  if (!log.scale) dat = dplyr::transmute(dat, X = x, Y = y)
  else dat = dplyr::transmute(dat, X = x.scaled, Y = y.scaled)
  rownames(dat) = rows
  class(dat) = append(class(dat), 'hierarchy')
  dat
}


plot_hierarchy = function(X,
                          quadrant.names = c('bl', 'br', 'tl', 'tr'),
                          main = NULL,
                          xlab = 'Relative meta-module score [log2(|SC1-SC2|+1)]',
                          ylab = 'Relative meta-module score [log2(|SC1-SC2|+1)]',
                          groups = NULL,
                          group.cols = NULL,
                          legend = T,
                          legend.pos = 'bottom',
                          legend.horiz = T) {

  if (is.null(groups)) col = 'darkred'
  else col = 'grey85'

  plot(X[,1], X[,2], pch = 20, col = col, main = main, xlab = xlab, ylab = ylab)

  if (is.null(groups)) legend = F
  else {
    stopifnot(!is.null(names(groups)))
    stopifnot(all(groups %in% rownames(X)))
    groups = split(groups, names(groups))
    Xgrp = sapply(groups, function(rows) X[rows,,drop = F], simplify = F)
    if (!is.null(group.cols)) colgrp = group.cols[names(groups)]
    else colgrp = rainbow(n = length(Xgrp))
    Map(points,
        x = sapply(Xgrp, `[[`, 1, simplify = F),
        y = sapply(Xgrp, `[[`, 2, simplify = F),
        col = colgrp,
        MoreArgs = list(pch = 20))
  }

  abline(v = 0, lty = 2)
  abline(h = 0, lty = 2)

  if (legend) {
    legend(legend.pos,
           fill = colgrp,
           legend = names(groups),
           horiz = legend.horiz,
           cex = 0.8,
           box.col = 'white',
           bg = 'white',
           box.lwd = 0)
  }

  Names = quadrant.names
  cex = 1.2
  mtext(side = 1, adj = 0, text = Names[1], cex = cex, line = cex - 1)
  mtext(side = 1, adj = 1, text = Names[2], cex = cex, line = cex - 1)
  mtext(side = 3, adj = 0, text = Names[3], cex = cex)
  mtext(side = 3, adj = 1, text = Names[4], cex = cex)
}
.check_missing <- function(Group, ref) {
  # all values in x should be in bins.
  if (is.numeric(ref)) ref = names(ref)
  are_missing = !Group %in% ref
  err = 'Error: Returning missing name(s)...'
  if (any(are_missing)) {
    missing = Group[are_missing]
    stop(cat(c(err, missing), sep = '\n'))
  }
}

.score <- function(mat, groups, controls = NULL, center = FALSE) {

  # Score matrix by groups
  # Scores are the column averages for the subset of rows specified in groups
  # Optional: Control group(s) scores to subtract from the real scores

  # Input check
  .check_missing(Group = unique(unlist(groups)), ref = rownames(mat))

  # Score without controls
  if (is.null(controls) | isFALSE(controls)) {
    s.mat = sapply(groups, function(p) colMeans(mat[p, ]))
  }

  # Score with controls
  else {
    if (isTRUE(controls)) {
      # Default control (all rows)
      control.vals = colMeans(mat)
      control.call = quote(control.vals[i])
    }
    else {
      # User-provided controls
      .check_missing(Group = unique(unlist(controls)), ref = rownames(mat))
      control.call = quote(colMeans(mat[controls[[i]], , drop = F]))
    }
    # Score
    groupids = setNames(1:length(groups), names(groups))
    score.call = quote(colMeans(mat[groups[[i]], , drop = F]))
    s.mat = sapply(groupids, function(i) eval(score.call) - eval(control.call))
  }

  # Center?
  if (center) {
    s.mat = colCenter(s.mat)
  }

  return(s.mat)
}



score <- function(mat,
                  groups,
                  binmat = NULL,
                  bins = NULL,
                  controls = NULL,
                  bin.control = F,
                  center = F,
                  nbin = 30,
                  n = 100,
                  replace = F) {

  # Wrapper for .score() with option to first generate control groups
  # Controls are generated with binmatch()

  # Generate controls if bins provided
  if (!is.null(bins|binmat)) {
    bin.control = TRUE
  }

  # Except if controls already provided
  if (!is.null(controls)) {
    bin.control = FALSE
  }

  if (bin.control && is.null(binmat) && all(unique(round(rowMeans(mat), 3)) == 0)) {
    message('Warning: if <mat> is row-centered and <binmat> was not provided, your bins will not be meaningful.')
  }
  # Get bin controls
  if (bin.control) {
    # Make bins if not provided
    if (is.null(bins)) {
      if (!is.null(binmat)) {
        bins = bin(mat = binmat, breaks = nbin)
      } else {
        bins = bin(mat = mat, breaks = nbin)
      }
    }

    # Make controls for groups using bins
    if (.arg_is_args(arg = groups)) {
      # If <groups> is many:
      controls = sapply(groups,
                        binmatch,
                        bins = bins,
                        n = n,
                        replace = replace,
                        simplify = F)
    } else {
      # If <groups> is one:
      controls = binmatch(Group = groups,
                          bins = bins,
                          n = n,
                          replace = replace) }
  }

  # Score (with controls or without if controls = NULL)
  .score(mat = mat, groups = groups, controls = controls, center = center)
}

.check_arg <- function(arg) {
  # x should contain one character vector only.
  if (is.list(arg) & all(lengths(arg) >= 1)) {
    stop('Please provide a single character vector.')
  }
}

.arg_is_args <- function(arg) {
  # TRUE if there is more than one char. vector in arg
  arg = ifelse (is.list(arg) & all(lengths(arg) >= 1), TRUE, FALSE)
  return(arg)
}

.check_sample_size <- function(bins, n, replace) {
  err1 = c('Error: Group is smaller than sample size <n> and <replace> == F.')
  err2 = c('Set <replace> = TRUE to proceed (or make <n> smaller).')
  errors <- c(err1, err2)
  if (any(lengths(bins) < n) & replace == FALSE) {
    stop(cat(errors, sep = '\n'))
  }
}

.name_binIDs <- function(binIDs, x) {
  if (!is.null(names(x))) {
    names(binIDs) = names(x)
  } else {
    names(binIDs) = x
  }
  return(binIDs)
}

.check_missing <- function(Group, ref) {
  # all values in x should be in bins.
  if (is.numeric(ref)) ref = names(ref)
  are_missing = !Group %in% ref
  err = 'Error: Returning missing name(s)...'
  if (any(are_missing)) {
    missing = Group[are_missing]
    stop(cat(c(err, missing), sep = '\n'))
  }
}

.check_args_exist <- function(x, data, bins) {
  args = list(x, data, bins)
  err1 = 'Not enough arguments. Provide one of:'
  err2 = '\t1. <data> : data to be transformed to <x>'
  err3 = '\t2. <x> : vector to be binned'
  err4 = '\t3. <bins>: bin vector for control <Group>'
  errors = c(err1, err2, err3, err4)
  if (all(sapply(args, is.null))) {
    stop(cat(errors, sep = '\n'))
  }
}

#' @title Center a matrix column-wise
#' @description Center a matrix column-wise
#' @param m a matrix or Matrix
#' @param by either "mean", "median" or a numeric vector of length equal to the number of columns of ‘m’. Default: "mean"
#' @return column-centered matrix
#' @rdname colcenter
#' @export
colCenter = function(m, by = 'mean') {
  m = as.matrix(m)
  if (by == 'mean')  by = T
  else if (by == 'median') by = matrixStats::colMedians(m)
  else stopifnot(is.numeric(by) & length(by) == ncol(m))
  scale(m, center = by, scale = F)
}
