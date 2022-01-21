# =============================================================================
#' @title Run UMAP.
#' @description This function runs UMAP.
#' @param obj The \code{espresso} object.
#' @param gset Gene set.
#' @param umap_param List of UMAP parameters (Refer to \code{uwot} package).
#' @param seed Random seed.
#' @importFrom uwot umap
#' @export
runUMAP <- function(obj, umap_param = NULL, gset = NULL, seed = 0) {
  set.seed(as.numeric(seed))
  if (!is.null(gset)) {
    gset <- gset[!is.na(match(gset, colnames(obj@exprs)))]
    genes <- gset
  } else {
    genes <- obj@gset
  }
  umap.defaults <- list(n_neighbors = 15,
                        n_components = 2,
                        metric = "euclidean",
                        n_epochs = NULL,
                        learning_rate = 1,
                        scale = FALSE,
                        init = "spectral",
                        init_sdev = NULL,
                        spread = 1,
                        min_dist = 0.01,
                        set_op_mix_ratio = 1,
                        local_connectivity = 1,
                        bandwidth = 1,
                        repulsion_strength = 1,
                        negative_sample_rate = 5,
                        a = NULL,
                        b = NULL,
                        nn_method = NULL,
                        n_trees = 50,
                        search_k = 2 * 15 * 50,
                        approx_pow = FALSE,
                        y = NULL,
                        target_n_neighbors = 15,
                        target_metric = "euclidean",
                        target_weight = 0.5,
                        pca = NULL,
                        pca_center = TRUE,
                        pcg_rand = TRUE,
                        fast_sgd = FALSE,
                        ret_model = FALSE,
                        ret_nn = FALSE,
                        ret_extra = c(),
                        n_threads = NULL,
                        n_sgd_threads = 0,
                        grain_size = 1,
                        tmpdir = tempdir(),
                        verbose = getOption("verbose", TRUE)
  )
  custom_settings <- umap.defaults
  if (!is.null(umap_param)) {
    for (i in 1:length(umap_param)) {
      param <- names(umap_param)[i]
      if (is.na(match(param, names(custom_settings)))) {
        warning(paste0("'", param, "' is invalid for the UMAP setting."))
        next
      }
      custom_settings[[param]] <- umap_param[[param]]
    }
  }
  if(is.null(umap_param[["search_k"]])) {
    custom_settings[["search_k"]] <-  2 * custom_settings[["n_neighbors"]] * custom_settings[["n_trees"]]
  }
  if(is.null(umap_param[["target_n_neighbors"]])) {
    custom_settings[["target_n_neighbors"]] <-  custom_settings[["n_neighbors"]]
  }
  exprs <- obj@exprs[, match(genes, colnames(obj@exprs)), drop = FALSE]
  obj@umap <- umap(exprs, 
                   n_neighbors = custom_settings[["n_neighbors"]],
                   n_components = custom_settings[["n_components"]],
                   metric = custom_settings[["metric"]],
                   n_epochs = custom_settings[["n_epochs"]],
                   learning_rate = custom_settings[["learning_rate"]],
                   scale = custom_settings[["scale"]],
                   init = custom_settings[["init"]],
                   init_sdev = custom_settings[["init_sdev"]],
                   spread = custom_settings[["spread"]],
                   min_dist = custom_settings[["min_dist"]],
                   set_op_mix_ratio = custom_settings[["set_op_mix_ratio"]],
                   local_connectivity = custom_settings[["local_connectivity"]],
                   bandwidth = custom_settings[["bandwidth"]],
                   repulsion_strength = custom_settings[["repulsion_strength"]],
                   negative_sample_rate = custom_settings[["negative_sample_rate"]],
                   a = custom_settings[["a"]],
                   b = custom_settings[["b"]],
                   nn_method = custom_settings[["nn_method"]],
                   n_trees = custom_settings[["n_trees"]],
                   search_k = custom_settings[["search_k"]],
                   approx_pow = custom_settings[["approx_pow"]],
                   y = custom_settings[["y"]],
                   target_n_neighbors = custom_settings[["target_n_neighbors"]],
                   target_metric = custom_settings[["target_metric"]],
                   target_weight = custom_settings[["target_weight"]],
                   pca = custom_settings[["pca"]],
                   pca_center = custom_settings[["pca_center"]],
                   pcg_rand = custom_settings[["pcg_rand"]],
                   fast_sgd = custom_settings[["fast_sgd"]],
                   ret_model = custom_settings[["ret_model"]],
                   ret_nn = custom_settings[["ret_nn"]],
                   ret_extra = custom_settings[["ret_extra"]],
                   n_threads = custom_settings[["n_threads"]],
                   n_sgd_threads = custom_settings[["n_sgd_threads"]],
                   grain_size = custom_settings[["grain_size"]],
                   tmpdir = custom_settings[["tmpdir"]],
                   verbose = custom_settings[["verbose"]]
                   )
  return(obj)
}

# =============================================================================
#' @title Plot UMAP.
#' @description This function plots UMAP.
#' @param obj The \code{espresso} object.
#' @param file Path to the output file (default: "umap.png").
#' @param movie Whether to make movie (default: FALSE).
#' @param mname Name of movie file (default: "movie").
#' @param dir Path to the output directory for the movie ("default: ".").
#' @param windowRect Vector of four values indicating the 
#'                   left, top, right and bottom of 
#'                   the displayed window (in pixels).
#' @param lit Lpecifying if lighting calculation should take place on geometry (default: TRUE).
#' @param pch Symbols to plot (default: 16).
#' @param size Size of the plotted points (default: 0.5).
#' @param cex Character expansion factor (default: 2).
#' @param inset Inset distance(s) from the margins (default: c(0.02)).
#' @param view_mat Transformation matrix.
#' @param zoom Zoom factor (default: 0.9).
#' @param axis Axis of rotation (default: c(0, 0, 1)).
#' @param rpm Rotation speed in rotations per minute (default:  5).
#' @param duration Duration of the animation (default: 10).
#' @param domcol List indicating domains and colors.
#' @param rgl Logical value indicating that 'rgl' is used for 3D plotting.
#' @param save Logical value indicating that the plot is saved or not. (default: TRUE).
#' @param rglclose Logical value indicating that the  RGL window is closed after plotting (default: FALSE).
#' @param color_by_gene Gene name given as characters (default: NULL).
#' @param lcex Character expansion factor for the legend (default: 1.2) 
#' @param color_by_group Data.frame indicating samples and their group.
#' @param grpcol List indicating group and colors.
#' @param magnify Multiplicative factor to apply to size of window when producing legend (default: 1).
#' @importFrom grDevices rainbow
#' @importFrom rgl open3d plot3d legend3d rgl.viewpoint rgl.snapshot rgl.close
#'             rgl.postscript movie3d spin3d bgplot3d
#' @importFrom scatterplot3d scatterplot3d
#' @importFrom tagcloud smoothPalette
#' @export
plotUMAP <- function(obj, file = "umap.pdf", 
                     movie = FALSE, dir = ".",
                     mname = "movie",
                     seed = 0,
                     windowRect = c(0, 200, 1000, 1200), 
                     lit = TRUE,
                     pch = 16, size = 0.5,
                     cex = 2, inset = c(0.02),
                     view_mat = c(0.8662382, -0.4996285, -0.001548417, 0, 
                                  0.2102456, 0.3617015, 0.908277810, 0, 
                                  -0.4532415, -0.7871108, 0.418364525, 0, 
                                  0, 0, 0, 1),
                     zoom = 0.9, axis = c(0,0,1), rpm = 5, duration = 10,
                     domcol = NULL, rgl = TRUE, save = TRUE, rglclose = FALSE, 
                     color_by_gene = NULL, lcex = 1.2,
                     color_by_group = NULL, grpcol = NULL, magnify = 1) {
  if (length(grep("\\.png$", file)) == 1){
    fmt <- 'png'
  } else if (length(grep("\\.pdf$", file)) == 1) {
    fmt <- 'pdf'
  } else if (length(grep("\\.eps$", file)) == 1) {
    fmt <- 'eps'
  } else {
    fmt <- NA
  }
  sample2domain <- obj@asgmt[match(rownames(obj@exprs), obj@asgmt$sample), ]
  if (is.null(color_by_gene) && is.null(color_by_group)) {
    if (is.null(domcol)) {
      domains <- sort(unique(obj@asgmt$domain))
      colors <- rainbow(length(domains))[match(sample2domain$domain, domains)]
      led_colors <- rainbow(length(domains))
    } else if (is.list(domcol)) {
      domains <- names(domcol)
      colors <- unlist(domcol)[match(sample2domain$domain, domains)]
      led_colors <- unlist(domcol)
    } else {
      stop("'color' must be given as a list.")
    }
  } else if (!is.null(color_by_gene) && is.null(color_by_group)) {
    if (length(color_by_gene) != 1) {
      stop("The length of 'color_by_gene' must be 1.")
      return(1)
    }
    if (!is.character(color_by_gene)) {
      stop("'color_by_gene' must be given as 'character'.")
      return(1)
    }
    if (any(color_by_gene == obj@gset)) {
      exprs_level <- scale(obj@exprs[, match(color_by_gene, colnames(obj@exprs))])
      colpal <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))
      colors <- smoothPalette(exprs_level, palfunc = colpal, n = 1000)
    } else {
      stop("No such gene.")
      return(1)
    }
  } else if (is.null(color_by_gene) && !is.null(color_by_group)) {
    sample2group <- color_by_group[match(rownames(obj@exprs), color_by_group$sample), ]
    if (is.null(grpcol)) {
      groups <- sort(unique(color_by_group$group))
      colors <- rainbow(length(groups))[match(sample2group$group, groups)]
      led_colors <- rainbow(length(groups))
    } else if (is.list(grpcol)) {
      groups <- names(grpcol)
      colors <- unlist(grpcol)[match(sample2group$group, groups)]
      led_colors <- unlist(grpcol)
    } else {
      stop("'color' must be given as a list.")
    }
  } else {
    stop("Coloring must be done either domain, gene, or group.")
    return(1)
  }
  if (ncol(obj@umap) == 3 && rgl == TRUE) {
    open3d(windowRect = windowRect)
    plot3d(obj@umap, pch = pch, col = colors, type = "s", size = size, 
           xlab = "Dim1", ylab = "Dim2", zlab = "Dim3", lit = lit, box = FALSE)
    if (is.null(color_by_gene) && is.null(color_by_group)) {
      legend3d("topright", legend = domains, pch = pch,
               col = led_colors, cex = cex, inset = inset, magnify = magnify)
    } else if (!is.null(color_by_gene) && is.null(color_by_group)) {
      bgplot3d({
        plot.new()
        .color.legend(0.975, 0.2, 1, 0.8,
                     rect_col = smoothPalette(seq(1,1000,1), palfunc = colpal, n = 1000),
                     level = exprs_level, lcex = lcex
                     )
      })
    } else if (is.null(color_by_gene) && !is.null(color_by_group)) {
      legend3d("topright", legend = groups, pch = pch,
               col = led_colors, cex = cex, inset = inset, magnify = magnify)
    } else {
      stop("Coloring must be done either domain, gene, or group.")
      return(1)
    }

    view_mat <- matrix(view_mat, ncol = 4, byrow = TRUE)
    rgl.viewpoint(userMatrix = view_mat, zoom = zoom)
    if (save == TRUE) {
      if (is.na(fmt)) {
        warning("Your plot was not saved. Supported export format: png, pdf, eps.")
      } else if (fmt == 'png') {
        rgl.snapshot(filename = file, fmt = "png")
      } else if (fmt == 'pdf') {
        rgl.postscript(filename = file, fmt = "pdf")
      } else if (fmt == 'eps') {
        rgl.postscript(filename = file, fmt = "eps")
      }
    }
    if (movie == TRUE) {
      if (is.null(dir)) {
        stop("The output directory should be specified.")
      }
      movie3d(spin3d(axis = axis, rpm = rpm), 
              duration = duration, dir = dir, movie = mname)
    }
    if (rglclose == TRUE) {
      rgl.close()
    }
  } else if (ncol(obj@umap) == 3 && rgl == FALSE) {
    if (save == TRUE) {
      if (is.na(fmt)) {
        warning("Your plot was not saved. Supported export format: png, pdf, eps.")
      } else if (fmt == 'png') {
        png(file)
      } else if (fmt == 'pdf') {
        pdf(file)
      } else if (fmt == 'eps') {
        postscript(file, horizontal = F, width = 7, height = 7)
      }
    }
    s3d <- scatterplot3d(obj@umap, color = colors, pch = pch,
                         xlab = "Dim1", ylab = "Dim2", zlab = "Dim3",
                         mar = c(8, 4, 7, 8)
                         )
    legend(par()$usr[2], par()$usr[4], legend=domains, col = led_colors, 
           inset = -0.28, xpd = TRUE, horiz = F, pch = pch, bty = "n")
    if (!is.na(fmt) && save == TRUE) {
      dev.off()
    }
  } else if (ncol(obj@umap) == 2) {
    if (save == TRUE) {
      if (is.na(fmt)) {
        warning("Your plot was not saved. Supported export format: png, pdf, eps.")
      } else if (fmt == 'png') {
        png(file)
      } else if (fmt == 'pdf') {
        pdf(file)
      } else if (fmt == 'eps') {
        postscript(file, horizontal = F, width = 7, height = 7)
      }
    }
    cex <- cex / 2
    oldmar <- par()$mar
    par(mar = c(7, 4, 5 ,8))
    plot(obj@umap, xlab = "Dim1", ylab = "Dim2",
         col = colors, pch = pch, cex = cex)
    par(xpd = TRUE)
    legend(par()$usr[2] + 0.2, par()$usr[4], 
           pch = pch, legend = domains, col = led_colors, cex = cex, bty = "n")
    par(mar = oldmar)
    if (!is.na(fmt) && save == TRUE) {
      dev.off()
    }
  } else {
    stop("Specify 'n_component' = 2 or 3 to draw the UMAP plot.")
  }
}

# =============================================================================
#' @title Draw a color legend.
#' @description This function draws a color legend.
#' @param xl,yb,xr,yt The lower left and upper right coordinates of the rectangle of colors.
#' @param rect_col The colors that will fill the rectangle.
#' @param lcex Character expansion factor for the labels.
#' @param level The values corresponds to the colors.
#' @importFrom plotrix gradient.rect
.color.legend <- function (xl, yb, xr, yt, rect_col, lcex, level){
  oldcex <- par("cex")
  par(xpd = TRUE, cex = lcex)
  gradient.rect(xl, yb, xr, yt, col = rect_col, nslices = length(rect_col), gradient = "y")
  texty0 <- yb + (0 - min(level)) / (max(level) - min(level)) * (yt - yb)
  unit_length <- 1 / (max(level) - min(level)) * (yt - yb)
  textx <- xr + 0.2 * strwidth("O") + 0.02
  textadj <- c(0, 0.5)
  if (unit_length > (yt - yb) / 3) {
    step_width <- 0.5
  } else {
    step_width <- 1.0
  }
  i <- 0
  while ((texty0 + unit_length * i) <= yt) {
    texty <- texty0 + unit_length * i
    text(textx, texty, labels = as.character(i), adj = textadj)
    segments(xr, texty, xr + 0.01, texty)
    i <- i + step_width
  }
  i <- step_width
  while ((texty0 - unit_length * i) >= yb) {
    texty <- texty0 - unit_length * i
    text(textx, texty, labels = as.character(-i), adj = textadj)
    segments(xr, texty, xr + 0.01, texty)
    i <- i + step_width
  }
  par(xpd = FALSE, cex = oldcex)
}
