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
  
  RCircos.Track.Outline(outerPos, innerPos, num.layers = 1, 
                        chromosomes, track.colors = rep("white", length(chromosomes)))
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
