#!/usr/bin/env Rscript

# rcircos_modified_functions.R
# Modified from RCircos package, v. 1.2.2, Henry Zhang for Yini Li lab use
# Jul. 2025


# Parameter reassignment function -----------------------------------------

# can modify more params with this
RCircos.Reset.Plot.Parameters.Custom <- function (new.params = NULL) {
  
  if (is.null(new.params)) 
    stop("Missing function argument.\n")
  
  old.params <- RCircos.Get.Plot.Parameters()
  
  # check for invalid inputs (can modify this later to allow more params to change)
  if (new.params$radius.len != old.params$radius.len || new.params$plot.radius != 
      old.params$plot.radius || new.params$tracks.inside != 
      old.params$tracks.inside || new.params$tracks.outside != 
      old.params$tracks.outside) { # consider redoing all these stupid if statements w/ logical vecs
    stop("Please use RCircos.Set.Core.Components() instead.\n")
  }
  
  if (new.params$highlight.pos != old.params$highlight.pos) {
    stop("Please use customized ideogram plot methods instead.\n")
  }
  
  RCircos.Validate.Plot.Parameters(new.params)
  
  if (new.params$chrom.width != old.params$chrom.width) {
    differ <- new.params$chrom.width - old.params$chrom.width
    new.params$highlight.pos <- old.params$highlight.pos + 
      differ
    new.name.pos <- old.params$chr.name.pos + differ
    
    if (new.params$chr.name.pos < new.name.pos) {
      new.params$chr.name.pos <- new.name.pos
    }
  }
  
  if (new.params$track.in.start >= new.params$chr.ideo.pos) {
    new.params$track.in.start <- new.params$chr.ideo.pos - 
    0.05
  }
  
  new.name.end <- new.params$chr.name.pos + 0.3
  
  if (new.params$track.out.start < new.name.end) {
    # new.params$track.out.start <- new.name.end
    new.params$track.out.start <- new.params$track.out.start
  }
  
  if (new.params$track.padding != old.params$track.padding || 
      new.params$track.height != old.params$track.height) {
    message(paste0("Track height and/or track padding have been ", 
                   "reset\n. Actual total data track plotted may differ.\n"))
  }
  
  if (old.params$base.per.unit != new.params$base.per.unit && 
      old.params$chrom.paddings == new.params$chrom.paddings) {
    
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
    
    band.len <- RCircos.Cyto$ChromEnd - RCircos.Cyto$ChromStart
    genome.len <- sum(as.numeric(band.len))
    padding.const <- RCircos.Get.Padding.Constant()
    total.units <- genome.len/new.params$base.per.unit
    new.padding <- round(padding.const * total.units, digits = 0)
    
    if (new.padding != new.params$base.per.unit) {
      message(paste("\nNote: chrom.padding", new.params$chrom.paddings, 
                    " was reset to", new.padding, "\n"))
      new.params$chrom.paddings <- new.padding
    }
  }
  
  RCircosEnvironment <- NULL
  RCircosEnvironment <- get("RCircos.Env", envir = globalenv())
  
  RCircosEnvironment[["RCircos.PlotPar"]] <- NULL
  RCircosEnvironment[["RCircos.PlotPar"]] <- new.params
  
  if (old.params$base.per.unit != new.params$base.per.unit || 
      old.params$chrom.paddings != new.params$chrom.paddings) {
    
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
    RCircos.Cyto <- RCircos.Cyto[, 1:5]
    
    RCircosEnvironment[["RCircos.Cytoband"]] <- NULL
    RCircos.Set.Cytoband.Data(RCircos.Cyto)
    
    RCircosEnvironment[["RCircos.Base.Position"]] <- NULL
    RCircos.Set.Base.Plot.Positions()
  }
}


# Ideogram plotting function ----------------------------------------------

# removes ticks and chr labels
RCircos.Chromosome.Ideogram.Plot.Custom <- function (tick.interval = 0) {
  RCircos.Draw.Chromosome.Ideogram.Custom()
}

# redefine ideogram draw function
# removes ticks
RCircos.Draw.Chromosome.Ideogram.Custom <- function (ideo.pos = NULL, ideo.width = NULL) 
{
  RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
  RCircos.Pos <- RCircos.Get.Plot.Positions()
  RCircos.Par <- RCircos.Get.Plot.Parameters()
  
  if (is.null(ideo.pos)) {
    ideo.pos <- RCircos.Par$chr.ideo.pos
  }
  
  if (is.null(ideo.width)) {
    ideo.width <- RCircos.Par$chrom.width
  }
  
  outerPos <- ideo.pos + ideo.width
  innerPos <- ideo.pos
  chromosomes <- unique(RCircos.Cyto$Chromosome)
  
  RCircos.Track.Outline.Custom(outerPos, innerPos, num.layers = 1, 
                        chromosomes, track.colors = rep(NA, length(chromosomes))) # clear not white
}


# Histogram function ------------------------------------------------------

# removes chr section outlines
RCircos.Histogram.Plot.Custom <- function (hist.data = NULL, data.col = 4, 
                                           track.num = NULL, side = c("in", "out"), 
                                           min.value = NULL, max.value = NULL, 
                                           inside.pos = NULL, outside.pos = NULL, 
                                           genomic.columns = 3, is.sorted = TRUE) {
  if (is.null(hist.data)) {
    stop("Genomic data missing in RCircos.Histogram.Plot().\n")
  }
  
  boundary <- RCircos.Get.Plot.Boundary(track.num, side, inside.pos, outside.pos, FALSE)
  
  outerPos <- boundary[1]
  innerPos <- boundary[2]
  
  if (is.null(genomic.columns) || genomic.columns < 2 || genomic.columns > 3) {
    stop("Incorrect number of columns for genomic position.\n")
  }
  
  if (is.null(data.col) || data.col <= genomic.columns) {
    stop("Hist data column must be greater than", genomic.columns, ".\n")
  }
  
  RCircos.Pos <- RCircos.Get.Plot.Positions()
  RCircos.Par <- RCircos.Get.Plot.Parameters()
  RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
  
  hist.data <- RCircos.Get.Single.Point.Positions(hist.data, genomic.columns)
  locations <- RCircos.Get.Start.End.Locations(hist.data, RCircos.Par$hist.width)
  histColors <- RCircos.Get.Plot.Colors(hist.data, RCircos.Par$hist.color)
  histValues <- as.numeric(hist.data[, data.col])
  
  if (is.null(max.value) || is.null(min.value)) {
    max.value <- max(histValues)
    min.value <- min(histValues)
    
  } else {
    
    if (min.value > max.value) {
      stop("min.value > max.value.")
    }
  }
  
  histHeight <- RCircos.Get.Data.Point.Height(
    histValues, 
    min.value, 
    max.value, 
    plot.type = "points", 
    outerPos - innerPos
    )
  
  # RCircos.Track.Outline(outerPos, innerPos, RCircos.Par$sub.tracks)
  for (aPoint in seq_len(nrow(hist.data))) {
    height <- innerPos + histHeight[aPoint]
    theStart <- locations[aPoint, 1]
    theEnd <- locations[aPoint, 2]
    
    polygonX <- c(RCircos.Pos[theStart:theEnd, 1] * height, 
                  RCircos.Pos[theEnd:theStart, 1] * innerPos)
    polygonY <- c(RCircos.Pos[theStart:theEnd, 2] * height, 
                  RCircos.Pos[theEnd:theStart, 2] * innerPos)
    
    polygon(polygonX, polygonY, col = histColors[aPoint], 
            border = NA)
  }
}


# Gene labels -------------------------------------------------------------

# allows tick marks on chromosome box things or whatever
RCircos.Gene.Connector.Plot.Custom <- function (genomic.data = NULL, track.num = NULL, side = "in", 
          inside.pos = NULL, outside.pos = NULL, genomic.columns = 3, 
          is.sorted = FALSE) 
{
  if (is.null(genomic.data)) 
    stop("Genomic data missing for RCircos.Gene.Connector.Plot().\n")
  boundary <- RCircos.Get.Plot.Boundary(track.num, side, inside.pos, 
                                        outside.pos, erase.area = FALSE)
  outerPos <- boundary[1] - 0.452  # MODIFICATION HERE
  innerPos <- boundary[2] - 0.198  # MODIFICATION HERE
  RCircos.Pos <- RCircos.Get.Plot.Positions()
  RCircos.Par <- RCircos.Get.Plot.Parameters()
  geneData <- RCircos.Get.Single.Point.Positions(genomic.data, 
                                                 genomic.columns)
  labelData <- RCircos.Get.Gene.Label.Locations(geneData, genomic.columns, 
                                                is.sorted)
  connectData <- data.frame(labelData$Location, labelData$LabelPosition)
  if (outerPos < RCircos.Par$chr.ideo.pos) {
    genomicCol <- ncol(connectData) - 1
    labelCol <- ncol(connectData)
  }
  else {
    genomicCol <- ncol(connectData)
    labelCol <- ncol(connectData) - 1
  }
  vHeight <- round((outerPos - innerPos)/10, digits = 4)
  hRange <- outerPos - innerPos - 2 * vHeight
  topLoc <- outerPos - vHeight
  botLoc <- innerPos + vHeight
  lineColors <- RCircos.Get.Plot.Colors(labelData, RCircos.Par$text.color)
  chroms <- unique(connectData[, 1])
  for (aChr in seq_along(chroms)) {
    chrRows <- which(connectData[, 1] == chroms[aChr])
    total <- length(chrRows)
    for (aPoint in seq_len(total)) {
      p1 <- connectData[chrRows[aPoint], genomicCol]
      p2 <- connectData[chrRows[aPoint], labelCol]
      lines(c(RCircos.Pos[p1, 1] * outerPos, RCircos.Pos[p1, 
                                                         1] * topLoc), c(RCircos.Pos[p1, 2] * outerPos, 
                                                                         RCircos.Pos[p1, 2] * topLoc), col = "red") # lineColors[chrRows[aPoint]])
      lines(c(RCircos.Pos[p2, 1] * botLoc, RCircos.Pos[p2, 
                                                       1] * innerPos), c(RCircos.Pos[p2, 2] * botLoc, 
                                                                         RCircos.Pos[p2, 2] * innerPos), col = "red") # lineColors[chrRows[aPoint]])
      lines(c(RCircos.Pos[p1, 1] * topLoc, RCircos.Pos[p2, 
                                                       1] * botLoc), c(RCircos.Pos[p1, 2] * topLoc, 
                                                                       RCircos.Pos[p2, 2] * botLoc), col = "red") # lineColors[chrRows[aPoint]])
    }
  }
}


# Chromosome outline function ---------------------------------------------

# get rid of chrY
RCircos.Track.Outline.Custom <- function (inside.pos = NULL, outside.pos = NULL, num.layers = 1, 
          chrom.list = NULL, track.colors = NULL) 
{
  if (is.null(outside.pos) || is.null(inside.pos)) 
    stop("Missing outside.pos/inside.pos in RCircos.Track.Outline().\n")
  RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
  RCircos.Pos <- RCircos.Get.Plot.Positions()
  RCircos.Par <- RCircos.Get.Plot.Parameters()
  subtrack.height <- (outside.pos - inside.pos)/num.layers
  chromosomes <- unique(as.character(RCircos.Cyto$Chromosome))
  if (!is.null(chrom.list)) {
    if (sum(chrom.list %in% chromosomes) != length(chrom.list)) {
      stop(paste("One or more chromosome is not", "in chromosome ideogram data.\n"))
    }
    chromosomes <- chrom.list
  }
  if (is.null(track.colors)) {
    track.colors <- rep(RCircos.Par$track.background, length(chromosomes))
  }
  else {
    if (length(track.colors) != length(chromosomes)) 
      track.colors <- rep(track.colors, length(chromosomes))
  }
  for (aChr in seq_len(length(chromosomes))) {
      if (!chromosomes[[aChr]] == "chrY") {
      chr.rows <- which(RCircos.Cyto$Chromosome == chromosomes[aChr])
      the.chr <- RCircos.Cyto[chr.rows, ]
      plot.start <- min(RCircos.Cyto$StartPoint[chr.rows])
      plot.end <- max(RCircos.Cyto$EndPoint[chr.rows])
      polygon.x <- c(RCircos.Pos[plot.start:plot.end, 1] * 
                       outside.pos, RCircos.Pos[plot.end:plot.start, 1] * 
                       inside.pos)
      polygon.y <- c(RCircos.Pos[plot.start:plot.end, 2] * 
                       outside.pos, RCircos.Pos[plot.end:plot.start, 2] * 
                       inside.pos)
      polygon(polygon.x, polygon.y, col = track.colors[aChr])
      if (num.layers > 1) {
        for (a.line in seq_len(num.layers - 1)) {
          height <- outside.pos - a.line * subtrack.height
          lines(RCircos.Pos[plot.start:plot.end, 1] * height, 
                RCircos.Pos[plot.start:plot.end, 2] * height, 
                col = RCircos.Par$grid.line.color)
        }
      }
    }
  }
}