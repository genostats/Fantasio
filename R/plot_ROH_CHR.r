plot_ROH_CHR <-
  function (ROH, unit, chr, list_id, regions, color2 = "green4") {
    if (unit == "cM") {
      pos1 <- which(colnames(ROH) == "POS1_cM")
      pos2 <- which(colnames(ROH) == "POS2_cM")
      myxlab <- "Position (cM)"
      coeff = 1
      if (length(pos1) == 0)
        stop("no genetic distances in ROH file")
    }
    else if (unit == "bases") {
      pos1 <- which(colnames(ROH) == "POS1")
      pos2 <- which(colnames(ROH) == "POS2")
      myxlab <- "Position (Mb)"
      coeff <- 1e6
    }
    else{
      stop("units parameter accepts 'cM' and 'bases' only")
    }
    
    start <- 0
    
    if (unit == "cM") {
      if ((chr == 13) || (chr == 14) || (chr == 15)) {
        start <- -20
      } else if ((chr == 21) || (chr == 22)) {
        start <- -15
      }
    }
    
    
    end <- quantsmooth::lengthChromosome(chr, unit) / coeff
    
    #plot vide
    y_max <- length(list_id) + 1
    
    plot(
      x = c(start, end),
      y <- c(0, y_max),
      type = "n",
      yaxt = "n",
      ylab = "",
      xlab = myxlab,
      main = paste("ROHs segments on chromosome ", chr, sep = "")
    )
    
    axis(
      2,
      at = (1:length(list_id)) + 0.25,
      list_id,
      col.ticks = 0,
      las = 2
    )
    
    #option regions
    
    #traitement de l'option regions
    if (!is.null(regions)) {
      if (nrow(regions) > 0) {
        #dessiner les regions
        for (i in 1:nrow(regions)) {
          polygon(
            x = regions[i, c("start", "end", "end", "start")] / coeff,
            y = c(0.75, 0.75, y_max, y_max),
            col = color2,
            border = color2,
            lwd = 2
          )
        }
      }
    }
    
    #dessiner le chromosome
    quantsmooth::paintCytobands(
      chr,
      units = unit,
      pos = c(0, 0.5),
      orientation = "h",
      legend = FALSE,
      length.out = end
    )
    
    for (j in 1:length(list_id)) {
      toplot <- ROH[ROH$IID == list_id[j],]
      
      if (nrow(toplot) > 0) {
        for (k in 1:nrow(toplot)) {
          polygon(
            x = c(toplot[k, pos1], toplot[k, pos2], toplot[k, pos2], toplot[k, pos1]) /
              coeff,
            y = c(j, j, j + 0.5, j + 0.5),
            col = color(toplot$PHE[k]),
            ### CHANEG !!!!
            lwd = 1
          )
        }
      }
    }
  }
