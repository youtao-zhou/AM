# BetterVisDotPlot Heatmap Seurat
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
BetterVis_DotPlot_Heatmap_Seurat = function(seurat_object,
                                            assay = "RNA",
                                            layer = "data",
                                            data.are.log = TRUE,
                                            features,
                                            split.by = NULL,
                                            idents = NULL,
                                            split.idents = NULL,
                                            scale = TRUE,
                                            rescale = FALSE,
                                            rescale.range = c(0, 3),
                                            rotate.axis = FALSE,
                                            dotplot = TRUE,
                                            dots.type = "square root",
                                            dots.size = 4,
                                            show.noexpr.dots = FALSE,
                                            col.min = ifelse(isTRUE(scale), -2, 0),
                                            col.max = ifelse(isTRUE(scale), 2, "q100"),
                                            data.colors = if (isTRUE(scale)) c("#35A5FF", "white", "red") else "Viridis",
                                            palette.reverse = FALSE,
                                            na.color = "grey40",
                                            background.color = "white",
                                            idents.colors = NULL,
                                            show.idents.names.colors = FALSE,
                                            show.idents.oppo.colors = TRUE,
                                            split.colors = NULL,
                                            show.split.names.colors = FALSE,
                                            show.split.oppo.colors = TRUE,
                                            order.idents = NULL,
                                            order.split = NULL,
                                            order.colors = TRUE,
                                            kmeans.repeats = 100,
                                            cluster.idents = TRUE,
                                            idents.kmeans = 1,
                                            idents.kmeans.numbers.size = 11,
                                            cluster.features = TRUE,
                                            features.kmeans = 1,
                                            features.kmeans.numbers.size = 11,
                                            idents.gap = 1,
                                            features.gap = 1,
                                            idents.names.size = 9,
                                            features.names.size = 9,
                                            features.names.style = "italic",
                                            row.names.side = "left",
                                            row.names.width = 15,
                                            column.names.angle = 45,
                                            column.names.side = "bottom",
                                            column.names.height = 15,
                                            inner.border = TRUE,
                                            outer.border = TRUE,
                                            data.legend.name = ifelse(isTRUE(scale), "Z-Score", "Average Expression"),
                                            data.legend.side = "bottom",
                                            data.legend.direction = "horizontal",
                                            data.legend.position = "topcenter",
                                            data.legend.width = 5,
                                            show.data.legend = TRUE,
                                            idents.legend.name = "Clusters",
                                            show.idents.legend = TRUE,
                                            split.legend.name = split.by,
                                            show.split.legend = TRUE,
                                            legend.title.size = 10,
                                            legend.text.size = 10,
                                            legend.gap = 10,
                                            output.data = FALSE,
                                            gene_markers = NULL,
                                            marker_colors = NULL,
                                            annotation_marker_legend = NULL,
                                            row.split = NULL,
                                            col.split = NULL,
                                            border_gp_control = FALSE,
                                            right.annotation = FALSE,
                                            bottom.annotation = FALSE,
                                            top.annotation = FALSE,
                                            keep_marker_order = FALSE,
                                            ...) {
  if (!(assay %in% names(seurat_object@assays))) {
    message("Assay '", assay, "' was not found in the Seurat object, using 'RNA' instead")
    assay <- "RNA"
  }
  if (!(layer %in% slotNames(seurat_object@assays[[assay]]))) {
    if ("data" %in% slotNames(seurat_object@assays[[assay]])) {
      message("Layer '", layer, "' was not found in the Seurat object's '", assay, "' assay, using 'data' instead")
      layer <- "data"
    } else {
      message("Layer '", layer, "' was not found in the Seurat object's '", assay, "' assay, using 'counts' instead")
      layer <- "counts"
    }
  }

  if (!is.character(split.by)) {
    ident <- "ident"
  } else {
    ident <- c("ident", split.by)
  }
  ident.1 <- levels(Idents(seurat_object))

  if (is.character(idents)) {
    ident.1 <- ident.1[ident.1 %in% idents]
    if (length(ident.1) == 0) stop("None of the identities supplied to idents were found")
    if (length(ident.1) < length(idents)) {
      message("The following identities supplied to idents were not found:\n", paste0(setdiff(idents, ident.1), collapse = ", "))
    }
  }
  if (is.character(order.idents) | is.numeric(order.idents)) {
    if (length(order.idents) == length(ident.1)) {
      if (is.character(order.idents)) ident.1 <- ident.1[order(match(ident.1, order.idents))] else ident.1 <- ident.1[order.idents]
    } else {
      if (length(order.idents) == 1 & any(order.idents == "reverse")) ident.1 <- rev(ident.1) else
        stop("order.idents needs to be 'reverse' or a character/numeric vector of same length as the number of identities")
    }
  }
  if (length(ident.1) == 1 & is.null(split.by)) stop("Only one identity in the active.ident identity, cannot compare features expression")

  if (is.character(split.by)) vars <- c("ident", split.by, features) else vars <- c("ident", features)
  DefaultAssay(seurat_object) <- assay
  data <- suppressWarnings(FetchData(object = seurat_object, vars = vars, layer = layer))
  features.removed <- setdiff(features, colnames(data))

  if (is.character(split.by)) {
    if (!split.by %in% colnames(data)) {
      message("The split.by identity was not found in the meta.data slot, data will not be split")
      split.by <- NULL
      ident.3 <- ident.1
    } else {
      Idents(seurat_object) <- split.by
      ident.2 <- levels(Idents(seurat_object))
      if (is.character(split.idents)) {
        ident.2 <- ident.2[ident.2 %in% split.idents]
        if (length(ident.2) == 0) stop("None of the identities supplied to split.idents were found")
        if (length(ident.2) < length(split.idents)) {
          message("The following identities supplied to split.idents were not found:\n", paste0(setdiff(split.idents, ident.2), collapse = ", "))
        }
      }
      if (is.character(order.split) | is.numeric(order.split)) {
        if (length(order.split) == length(ident.2)) {
          if (is.character(order.split)) ident.2 <- ident.2[order(match(ident.2, order.split))] else ident.2 <- ident.2[order.split]
        } else {
          if (length(order.split) == 1 & any(order.split == "reverse")) ident.2 <- rev(ident.2) else
            stop("order.split needs to be 'reverse' or a character/numeric vector of same length as the number of identities")
        }
      }
      ident.3 <- NULL
      for (i in ident.1) for (j in ident.2) ident.3 <- c(ident.3, paste(i, "in", j))
    }
    if (length(ident.3) == 1) stop("Only one identity in the active.ident and split.by identities, cannot compare features expression")
  } else {
    ident.3 <- ident.1
  }

  if (is.character(split.by)) {
    data$ident <- paste(data$ident, "in", data[, 2])
    data[, 2] <- NULL
  }
  if (ncol(data) == 1) stop("None of the features were found")

  mat <- dot <- list()
  for (i in 2:ncol(data)) {
    mat[[i - 1]] <- data[, c(1, i), drop = FALSE]
    avg <- pct <- list()
    for (j in ident.3) {
      avg[[j]] <- mat[[i - 1]][mat[[i - 1]]$ident == j, 2, drop = FALSE]
      pct[[j]] <- as.data.frame(length(avg[[j]][avg[[j]] > 0]) / length(avg[[j]][, 1]))
      if (isTRUE(data.are.log) & layer == "data") avg[[j]] <- as.data.frame(mean(expm1(avg[[j]][, 1]))) else
        avg[[j]] <- as.data.frame(mean(avg[[j]][, 1]))
      colnames(avg[[j]]) <- colnames(data)[i]; rownames(avg[[j]]) <- j
      colnames(pct[[j]]) <- colnames(data)[i]; rownames(pct[[j]]) <- j
    }
    mat[[i - 1]] <- do.call(rbind, avg)
    dot[[i - 1]] <- do.call(rbind, pct)
  }
  mat <- as.matrix(as.data.frame(mat))[ident.3, , drop = FALSE]; colnames(mat) <- colnames(data)[2:ncol(data)]
  dot <- as.matrix(as.data.frame(dot))[ident.3, , drop = FALSE]; colnames(dot) <- colnames(data)[2:ncol(data)]
  if (isTRUE(data.are.log) & layer == "data") mat <- log1p(mat)
  dot <- dot[rowSums(is.na(mat)) < ncol(mat), , drop = FALSE]
  mat <- mat[rowSums(is.na(mat)) < ncol(mat), , drop = FALSE]
  if (nrow(mat) < 2) stop("Less than two identities left, cannot compare features expression")
  idents.removed <- setdiff(ident.3, rownames(mat))
  dot <- dot[, colSums(mat) > 0, drop = FALSE]
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  if (ncol(mat) == 0) stop("None of the features were expressed in any cells")
  features.removed <- c(features.removed, setdiff(colnames(data), colnames(mat)))
  features.removed <- features.removed[!features.removed == "ident"]
  if (length(features.removed) > 0) {
    message("The following features were removed as they were not found or were not expressed in any cells:\n", paste0(features.removed, collapse = ", "))
  }

  if (isTRUE(scale)) {
    if (nrow(mat) < 5 & isFALSE(output.data)) {
      message("Less than 5 identities will be displayed, scaling may produce misleading visualization.\nUsing Cell_Heatmap() might be more appropriate.")
    }
    mat <- scale(mat)
  }
  if (isTRUE(rescale) & isFALSE(scale)) {
    if (nrow(mat) < 5 & isFALSE(output.data)) {
      message("Less than 5 identities will be displayed, rescaling may produce misleading visualization.\nUsing Cell_Heatmap() might be more appropriate.")
    }
    mat <- apply(mat, 2, function(x) rescale(x, to = rescale.range))
  }

  if (is.character(col.min)) {
    if (isTRUE(grepl("^q[0-9]{1,3}$", col.min, perl = TRUE))) {
      q1 <- as.numeric(sub("q", "", col.min)) / 100
      q1 <- quantile(mat, probs = q1, na.rm = TRUE)
      col.min <- q1
    }
  }
  if (is.character(col.max)) {
    if (isTRUE(grepl("^q[0-9]{1,3}$", col.max, perl = TRUE))) {
      q2 <- as.numeric(sub("q", "", col.max)) / 100
      q2 <- quantile(mat, probs = q2, na.rm = TRUE)
      col.max <- q2
    }
  }

  if (isTRUE(rotate.axis)) {
    mat <- t(mat); dot <- t(dot)
  }

  if (isTRUE(output.data)) return(list(Expression = round(mat, 4), Percent = round(dot, 4)))

  if (isTRUE(scale)) {
    col.mid <- 0
  } else {
    col.mid <- (col.min + col.max) / 2
  }
  if (!is.function(data.colors)) {
    if (length(data.colors) > 3) stop("data.colors must have 2 or 3 colors or must be a palette")
    if (length(data.colors) == 3) {
      data.colors <- colorRamp2(breaks = c(col.min, col.mid, col.max), colors = c(data.colors[1], data.colors[2], data.colors[3]))
    } else if (length(data.colors) == 2) {
      data.colors <- colorRamp2(breaks = c(col.min, col.max), colors = c(data.colors[1], data.colors[2]))
    } else {
      data.colors <- colorRamp2(breaks = c(col.min, col.max), hcl_palette = data.colors, reverse = palette.reverse)
    }
  }

  if (!is.character(split.colors) & is.character(split.by)) {
    current.seed <- .Random.seed; set.seed(1234)
    CustomColors <- colors()[grep('gr(a|e)y|white|black', colors(), invert = TRUE)]
    split.colors <- sample(CustomColors, length(ident.2))
    .Random.seed <- current.seed
  }

  if (!is.character(idents.colors)) idents.colors <- hue_pal()(n = length(ident.1))
  if (isTRUE(order.colors)) {
    if (!is.null(names(idents.colors))) idents.colors <- idents.colors[ident.1]
    if (is.character(order.idents)) {
      if (length(order.idents) > 1) {
        names(idents.colors) <- ident.1; idents.colors <- idents.colors[order.idents]
      } else {
        if (order.idents == "reverse") idents.colors <- rev(idents.colors)
      }
    }
    if (is.character(split.by)) {
      if (!is.null(names(split.colors))) split.colors <- split.colors[ident.2]
      if (is.character(order.split)) {
        if (length(order.split) > 1) {
          names(split.colors) <- ident.2; split.colors <- split.colors[order.split]
        } else {
          if (order.split == "reverse") split.colors <- rev(split.colors)
        }
      }
    }
  }
  idents.colors2 <- idents.colors
  if (is.character(split.by)) {
    split.colors2 <- split.colors
    dup.colors <- dup.colors2 <- list()
    for (i in 1:length(idents.colors)) {
      dup.colors[[i]] <- rep(idents.colors[i], length(ident.3) / length(ident.1))
      dup.colors2[[i]] <- split.colors
    }
    idents.colors <- unlist(dup.colors)
    split.colors <- unlist(dup.colors2)
    names(idents.colors) <- names(split.colors) <- ident.3
    idents.colors <- idents.colors[!names(idents.colors) %in% idents.removed]
    split.colors <- split.colors[!names(split.colors) %in% idents.removed]
  } else {
    names(idents.colors) <- ident.3
    idents.colors <- idents.colors[!names(idents.colors) %in% idents.removed]
  }

  idents.legend <- split.legend <- NULL
  if (isTRUE(show.idents.legend) & (isTRUE(show.idents.names.colors) | isTRUE(show.idents.oppo.colors))) {
    idents.legend <- Legend(
      at = ident.1,
      legend_gp = gpar(fill = idents.colors2),
      title = idents.legend.name,
      gap = unit(0.5, "cm"),
      border = TRUE,
      title_gp = gpar(fontsize = legend.title.size, fontface = "bold"),
      labels_gp = gpar(fontsize = legend.text.size)
    )
  }
  if (is.character(split.by) & isTRUE(show.split.legend) & (isTRUE(show.split.names.colors) | isTRUE(show.split.oppo.colors))) {
    split.legend <- Legend(
      at = ident.2,
      legend_gp = gpar(fill = split.colors2),
      title = split.legend.name,
      gap = unit(0.5, "cm"),
      border = TRUE,
      title_gp = gpar(fontsize = legend.title.size, fontface = "bold"),
      labels_gp = gpar(fontsize = legend.text.size)
    )
  }
  anno.legend <- list()
  if (!is.null(idents.legend)) anno.legend <- c(anno.legend, list(idents.legend))
  if (!is.null(split.legend)) anno.legend <- c(anno.legend, list(split.legend))

  if (isFALSE(rotate.axis)) annotation.side <- "row" else annotation.side <- "column"
  idents.labels3 <- idents.labels4 <- NULL
  row.dend.side <- row.title.side <- "right"
  column.dend.side <- column.title.side <- "top"
  idents.kmeans.color <- ifelse(idents.kmeans.numbers.size > 0, "black", "white")
  features.kmeans.color <- ifelse(features.kmeans.numbers.size > 0, "black", "white")

  if (is.character(split.by)) {
    if (isTRUE(show.split.names.colors) & isTRUE(show.idents.names.colors)) {
      if ((isFALSE(rotate.axis) & row.names.side == "left") | (isTRUE(rotate.axis) & column.names.side == "top")) {
        idents.labels <- HeatmapAnnotation(
          Identities = names(idents.colors),
          Split = names(split.colors),
          col = list(Identities = idents.colors, Split = split.colors),
          na_col = na.color, which = annotation.side,
          show_annotation_name = FALSE, show_legend = FALSE
        )
      } else {
        idents.labels <- HeatmapAnnotation(
          Split = names(split.colors),
          Identities = names(idents.colors),
          col = list(Split = split.colors, Identities = idents.colors),
          na_col = na.color, which = annotation.side,
          show_annotation_name = FALSE, show_legend = FALSE
        )
      }
    } else if (isTRUE(show.idents.names.colors)) {
      idents.labels <- HeatmapAnnotation(
        Identities = names(idents.colors),
        col = list(Identities = idents.colors),
        na_col = na.color, which = annotation.side,
        show_annotation_name = FALSE, show_legend = FALSE
      )
    } else if (isTRUE(show.split.names.colors)) {
      idents.labels <- HeatmapAnnotation(
        Split = names(split.colors),
        col = list(Split = split.colors),
        na_col = na.color, which = annotation.side,
        show_annotation_name = FALSE, show_legend = FALSE
      )
    } else {
      idents.labels <- NULL
    }
    if (isTRUE(show.split.oppo.colors) & isTRUE(show.idents.oppo.colors)) {
      if ((isFALSE(rotate.axis) & row.names.side == "left") | (isTRUE(rotate.axis) & column.names.side == "top")) {
        idents.labels2 <- HeatmapAnnotation(
          Split = names(split.colors),
          Identities = names(idents.colors),
          col = list(Split = split.colors, Identities = idents.colors),
          na_col = na.color, which = annotation.side,
          show_annotation_name = FALSE, show_legend = FALSE
        )
      } else {
        idents.labels2 <- HeatmapAnnotation(
          Identities = names(idents.colors),
          Split = names(split.colors),
          col = list(Identities = idents.colors, Split = split.colors),
          na_col = na.color, which = annotation.side,
          show_annotation_name = FALSE, show_legend = FALSE
        )
      }
    } else if (isTRUE(show.idents.oppo.colors)) {
      idents.labels2 <- HeatmapAnnotation(
        Identities = names(idents.colors),
        col = list(Identities = idents.colors),
        na_col = na.color, which = annotation.side,
        show_annotation_name = FALSE, show_legend = FALSE
      )
    } else if (isTRUE(show.split.oppo.colors)) {
      idents.labels2 <- HeatmapAnnotation(
        Split = names(split.colors),
        col = list(Split = split.colors),
        na_col = na.color, which = annotation.side,
        show_annotation_name = FALSE, show_legend = FALSE
      )
    } else {
      idents.labels2 <- NULL
    }
  } else {
    if (isTRUE(show.idents.names.colors)) {
      idents.labels <- HeatmapAnnotation(
        Identities = names(idents.colors),
        col = list(Identities = idents.colors),
        na_col = na.color, which = annotation.side,
        show_annotation_name = FALSE, show_legend = FALSE
      )
    } else {
      idents.labels <- NULL
    }
    if (isTRUE(show.idents.oppo.colors)) {
      idents.labels2 <- HeatmapAnnotation(
        Identities = names(idents.colors),
        col = list(Identities = idents.colors),
        na_col = na.color, which = annotation.side,
        show_annotation_name = FALSE, show_legend = FALSE
      )
    } else {
      idents.labels2 <- NULL
    }
  }
  if (row.names.side == "right") {
    if (isFALSE(rotate.axis)) {
      tmp <- idents.labels; idents.labels <- idents.labels2; idents.labels2 <- tmp
    }
    row.dend.side <- row.title.side <- "left"
  }
  idents.names.gp <- gpar(fontsize = idents.names.size)
  features.names.gp <- gpar(fontsize = features.names.size, fontface = features.names.style)
  idents.title.gp <- gpar(fontsize = idents.kmeans.numbers.size, col = idents.kmeans.color)
  features.title.gp <- gpar(fontsize = features.kmeans.numbers.size, col = features.kmeans.color)
  if (isTRUE(rotate.axis)) {
    idents.labels3 <- idents.labels; idents.labels <- NULL
    idents.labels4 <- idents.labels2; idents.labels2 <- NULL
    tmp <- idents.names.gp; idents.names.gp <- features.names.gp; features.names.gp <- tmp
    tmp <- idents.names.size; idents.names.size <- features.names.size; features.names.size <- tmp
    tmp <- cluster.idents; cluster.idents <- cluster.features; cluster.features <- tmp
    tmp <- idents.kmeans; idents.kmeans <- features.kmeans; features.kmeans <- tmp
    tmp <- idents.title.gp; idents.title.gp <- features.title.gp; features.title.gp <- tmp
    tmp <- idents.gap; idents.gap <- features.gap; features.gap <- tmp
  }
  if (column.names.side == "top") {
    if (isTRUE(rotate.axis)) {
      tmp <- idents.labels3; idents.labels3 <- idents.labels4; idents.labels4 <- tmp
    }
    column.dend.side <- column.title.side <- "bottom"
  }

  if (isFALSE(dotplot)) {
    if (sum(is.na(mat)) > 0) {
      na.label <- list(Legend(
        at = "NA", legend_gp = gpar(fill = na.color), title = NULL,
        gap = unit(0.5, "cm"), border = TRUE
      ))
    } else na.label <- list()
    anno.legend <- c(na.label, anno.legend)
  }
  ht_opt$ANNOTATION_LEGEND_PADDING <- unit(legend.gap, "mm")
  if (data.legend.side == "right" | data.legend.side == "left" | column.names.side == "top") {
    ht_opt$HEATMAP_LEGEND_PADDING <- unit(legend.gap, "mm")
  } else {
    ht_opt$HEATMAP_LEGEND_PADDING <- unit(2, "mm")
  }

  gene_annotation <- NULL
  if (!is.null(gene_markers)) {
    if (is.null(names(gene_markers)) && length(gene_markers) == nrow(mat)) names(gene_markers) <- rownames(mat)

    if (!is.null(names(gene_markers))) gene_markers <- gene_markers[rownames(mat)]

    if (isTRUE(keep_marker_order)) {
      genes_in_order <- names(gene_markers)[!is.na(gene_markers) & names(gene_markers) %in% rownames(mat)]
      if (length(genes_in_order) > 0) {
        mat <- mat[genes_in_order, , drop = FALSE]
        dot <- dot[genes_in_order, , drop = FALSE]
        gene_markers <- gene_markers[genes_in_order]
      }
    } else {
      if (!all(is.na(gene_markers))) {
        gm <- as.factor(gene_markers)
        ord <- order(gm, na.last = TRUE)
        mat <- mat[ord, , drop = FALSE]; dot <- dot[ord, , drop = FALSE]; gene_markers <- droplevels(gm[ord])
      }
    }

    if (is.null(marker_colors)) {
      lv <- levels(factor(gene_markers))
      auto_cols <- hue_pal()(length(lv)); names(auto_cols) <- lv
      marker_colors <- auto_cols
    }

    if (length(gene_markers) == nrow(mat) && !all(is.na(gene_markers))) {
      if (isTRUE(annotation_marker_legend)) {
        gene_annotation <- rowAnnotation(
          Marker = gene_markers,
          col = list(Marker = marker_colors),
          gap = unit(3, "mm"),
          show_annotation_name = FALSE
        )
      } else {
        gene_annotation <- rowAnnotation(
          Marker = gene_markers,
          col = list(Marker = marker_colors),
          gap = unit(3, "mm"),
          show_annotation_name = FALSE,
          annotation_legend_param = list(Marker = list(title = "", at = character(0), labels = character(0))),
          show_legend = c("Marker" = FALSE)
        )
      }
    } else {
      warning("`gene_markers` could not be aligned to rows; left annotation will be disabled.")
      gene_annotation <- NULL
    }
  }

  if (isTRUE(dotplot)) {
    zero.dot <- 0
    if (isTRUE(show.noexpr.dots)) {
      zero.dot <- min(dot[dot > 0], na.rm = TRUE); if (zero.dot > 0.001) zero.dot <- 0.001
      dot[dot == 0] <- zero.dot
    }
    if (dots.type == "radius") {
      dots.creation <- function(j, i, x, y, w, h, fill) {
        grid.rect(x = x, y = y, width = w, height = h, gp = gpar(col = NA, fill = background.color))
        ifelse(isTRUE(inner.border),
          grid.circle(
            x = x, y = y, r = pindex(dot, i, j) * unit(dots.size, "mm"),
            gp = gpar(fill = data.colors(pindex(mat, i, j)), col = "black", lty = 1)
          ),
          grid.circle(
            x = x, y = y, r = pindex(dot, i, j) * unit(dots.size, "mm"),
            gp = gpar(fill = data.colors(pindex(mat, i, j)), col = NA)
          )
        )
      }
      dots.legend <- Legend(
        labels = c(0, 10, 25, 50, 75, 100), title = "Percent\nof cells\nexpressing",
        row_gap = unit(legend.text.size / 15, "mm"),
        title_gp = gpar(fontsize = legend.title.size, fontface = "bold"),
        labels_gp = gpar(fontsize = legend.text.size),
        graphics = list(
          function(x, y, w, h) grid.circle(x = x, y = y, r = zero.dot * unit(legend.text.size / 6, "mm"), gp = gpar(fill = "black")),
          function(x, y, w, h) grid.circle(x = x, y = y, r = 0.1 * unit(legend.text.size / 6, "mm"), gp = gpar(fill = "black")),
          function(x, y, w, h) grid.circle(x = x, y = y, r = 0.25 * unit(legend.text.size / 6, "mm"), gp = gpar(fill = "black")),
          function(x, y, w, h) grid.circle(x = x, y = y, r = 0.5 * unit(legend.text.size / 6, "mm"), gp = gpar(fill = "black")),
          function(x, y, w, h) grid.circle(x = x, y = y, r = 0.75 * unit(legend.text.size / 6, "mm"), gp = gpar(fill = "black")),
          function(x, y, w, h) grid.circle(x = x, y = y, r = 1 * unit(legend.text.size / 6, "mm"), gp = gpar(fill = "black"))
        )
      )
    }
    if (dots.type == "square root") {
      dots.creation <- function(j, i, x, y, w, h, fill) {
        grid.rect(x = x, y = y, width = w, height = h, gp = gpar(col = "gray", fill = "white", lwd = 0.5))
        ifelse(isTRUE(inner.border),
          grid.circle(
            x = x, y = y, r = sqrt(pindex(dot, i, j)) * unit(dots.size, "mm"),
            gp = gpar(fill = data.colors(pindex(mat, i, j)), col = "black", lty = 1)
          ),
          grid.circle(
            x = x, y = y, r = sqrt(pindex(dot, i, j)) * unit(dots.size, "mm"),
            gp = gpar(fill = data.colors(pindex(mat, i, j)), col = NA)
          )
        )
      }
    }
    if (data.legend.side == "right") {
      show.data.legend <- FALSE
      if (isTRUE(scale)) data.legend.breaks <- c(col.min, col.min / 2, 0, col.max / 2, col.max) else
        data.legend.breaks <- c(col.min, col.min, col.max / 2, col.max)
      anno.legend <- c(list(Legend(
        col_fun = data.colors, title = gsub("\\ ", "\n", data.legend.name),
        legend_height = unit(data.legend.width, "cm"),
        title_gp = gpar(fontsize = legend.title.size, fontface = "bold"),
        labels_gp = gpar(fontsize = legend.text.size),
        grid_width = unit(0.5, "cm"),
        title_gap = unit(legend.title.size / 4, "mm")
      )), anno.legend)
    }

    ht <- Heatmap(
      mat,
      col = data.colors,
      na_col = na.color,
      name = data.legend.name,
      rect_gp = gpar(type = "none"),
      border = outer.border,
      layer_fun = dots.creation,
      left_annotation = if (!is.null(gene_annotation)) gene_annotation else idents.labels,
      border_gp = if (border_gp_control) gpar() else gpar(col = NA),
      row_split = row.split,
      row_title = NULL,
      column_split = col.split,
      column_title = NULL,
      right_annotation = if (right.annotation) idents.labels2 else NULL,
      bottom_annotation = if (bottom.annotation) idents.labels3 else NULL,
      top_annotation = if (top.annotation) idents.labels4 else NULL,
      cluster_rows = cluster.idents,
      cluster_columns = cluster.features,
      row_title_rot = 0,
      row_title_side = row.title.side,
      row_title_gp = idents.title.gp,
      row_km = idents.kmeans,
      row_km_repeats = kmeans.repeats,
      column_km = features.kmeans,
      column_km_repeats = kmeans.repeats,
      show_column_names = ifelse(features.names.size > 0, TRUE, FALSE),
      show_row_names = ifelse(idents.names.size > 0, TRUE, FALSE),
      column_names_gp = features.names.gp,
      row_names_gp = idents.names.gp,
      column_names_rot = column.names.angle,
      column_title_side = column.title.side,
      column_title_gp = features.title.gp,
      row_names_side = row.names.side,
      row_names_max_width = unit(row.names.width, "cm"),
      column_names_side = column.names.side,
      column_names_max_height = unit(column.names.height, "cm"),
      row_dend_side = row.dend.side,
      column_dend_side = column.dend.side,
      show_heatmap_legend = show.data.legend,
      heatmap_legend_param = list(
        legend_direction = data.legend.direction,
        title_position = data.legend.position,
        legend_width = unit(data.legend.width, "cm"),
        title_gp = gpar(fontsize = legend.title.size, fontface = "bold"),
        labels_gp = gpar(fontsize = legend.text.size)
      ),
      ...
    )
  } else {
    if (isTRUE(inner.border)) inner.border <- gpar(col = "black") else inner.border <- gpar(col = NA)
    ht <- Heatmap(
      mat,
      col = data.colors,
      na_col = na.color,
      name = data.legend.name,
      rect_gp = inner.border,
      border = outer.border,
      left_annotation = if (!is.null(gene_annotation)) gene_annotation else idents.labels,
      right_annotation = idents.labels2,
      bottom_annotation = idents.labels3,
      top_annotation = idents.labels4,
      cluster_rows = cluster.idents,
      cluster_columns = cluster.features,
      row_title_rot = 0,
      row_title_side = row.title.side,
      row_title_gp = gpar(fontsize = idents.kmeans.numbers.size),
      row_km = idents.kmeans,
      row_km_repeats = kmeans.repeats,
      column_km = features.kmeans,
      column_km_repeats = kmeans.repeats,
      column_names_gp = features.names.gp,
      row_names_gp = idents.names.gp,
      column_names_rot = column.names.angle,
      column_title_side = column.title.side,
      column_title_gp = gpar(fontsize = features.kmeans.numbers.size),
      row_names_side = row.names.side,
      row_names_max_width = unit(row.names.width, "cm"),
      column_names_side = column.names.side,
      column_names_max_height = unit(column.names.height, "cm"),
      row_dend_side = row.dend.side,
      column_dend_side = column.dend.side,
      row_gap = unit(idents.gap, "mm"),
      column_gap = unit(3, "mm"),
      show_heatmap_legend = show.data.legend,
      heatmap_legend_param = list(
        legend_direction = data.legend.direction,
        title_position = data.legend.position,
        legend_width = unit(data.legend.width, "cm"),
        title_gp = gpar(fontsize = legend.title.size, fontface = "bold"),
        labels_gp = gpar(fontsize = legend.text.size)
      ),
      ...
    )
  }

  return(draw(ht,
    heatmap_legend_side = data.legend.side,
    align_heatmap_legend = "heatmap_center",
    align_annotation_legend = "heatmap_center",
    annotation_legend_list = anno.legend,
    legend_grouping = "original"
  ))
}
