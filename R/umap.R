# =============================================================================
#' @title Plot 3D UMAP
#' @description This function plots UMAP
#' @param obj \code{espresso} object.
#' @param file path to the output directory.
#' @param gset vector of genes.
#' @param movie whether to make movie.
#' @param dir path to the output directory for the movie.
#' @param seed random seed
#' @param windowRect vector of four values indicating the 
#'                   left, top, right and bottom of 
#'                   the displayed window (in pixels).
#' @param lit specifying if lighting calculation should take place on geometry
#' @param pch symbols to plot
#' @param size size of the plotted points.
#' @param lpch symbols in the legend.
#' @param cex character expansion factor.
#' @param inset inset distance(s) from the margins.
#' @param view_mat the transformation matrix
#' @param zoom zoom factor.
#' @param axis axis of rotation.
#' @param rpm rotation speed in rotations per minute.
#' @param duration the duration of the animation.
#' @param domcol list indicating domains and colors.
#' @importFrom grDevices rainbow
#' @importFrom umap umap.defaults umap
#' @importFrom rgl open3d plot3d legend3d rgl.viewpoint rgl.snapshot 
#'             rgl.postscript movie3d spin3d
#' @export
#' 
plotUMAP <- function(obj, file = "umap.png", gset = NULL, 
                     movie = FALSE, dir = ".",
                     umap_param = NULL,
                     seed = 0,
                     windowRect = c(0, 200, 1000, 1200), 
                     lit =TRUE,
                     pch = 21, size = 0.5, lpch = 16,
                     cex = 2, inset = c(0.02),
                     view_mat = c(0.8662382, -0.4996285, -0.001548417, 0, 
                                  0.2102456, 0.3617015, 0.908277810, 0, 
                                  -0.4532415, -0.7871108, 0.418364525, 0, 
                                  0, 0, 0, 1),
                     zoom = 0.9, axis = c(0,0,1), rpm = 5, duration = 10,
                     domcol = NULL) {
  set.seed(as.numeric(seed))
  if (!is.null(gset)) {
    gset <- gset[!is.na(match(gset, colnames(obj@exprs)))]
    genes <- gset
  } else {
    genes <- obj@gset
  }
  custom_settings = umap.defaults
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
  sample2domain <- obj@asgmt[match(rownames(obj@exprs), obj@asgmt$sample), ]
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
  exprs <- obj@exprs[, match(genes, colnames(obj@exprs)), drop = FALSE]
  umap_out <- umap(exprs, config = custom_settings)
  open3d(windowRect = windowRect)
  plot3d(umap_out$layout, pch = pch, col = colors, type = "s", size = size, 
         xlab = "X", ylab = "Y", zlab = "Z", lit = lit, box = FALSE)
  legend3d("topright", legend = domains, pch = lpch, 
           col = led_colors, cex = cex, inset = inset)
  view_mat <- matrix(view_mat, ncol = 4, byrow = TRUE)
  rgl.viewpoint(userMatrix = view_mat, zoom = zoom)
  rgl.snapshot(filename = file, fmt = "png")
  if (movie == TRUE) {
    if (is.null(dir)) {
      stop("The output directory should be specified.")
    }
    movie3d(spin3d(axis = axis, rpm = rpm), 
            duration = duration, dir = dir, movie="movie")
  }
}