### Visualization helper functions

#' Plot a 2D Gamma Warp With Scalar Input
#'
#' Draw arrows from each `(u, v)` grid point to `gamma_fun(u, v)`. This is useful
#' for inspecting a registration warp in parameter space.
#'
#' @param gamma_fun Gamma function with scalar input `gamma_fun(u, v)` returning
#' a length-2 vector.
#' @param n Grid resolution used when `uv_grid` is not supplied.
#' @param arrow_color Arrow color.
#' @param point_color Original grid-point color.
#' @param warp_scale Multiplier applied to the displacement vectors for display.
#' @param uv_grid Optional data frame or matrix with first two columns `u`, `v`.
#'
#' @return A ggplot object.
#' @import ggplot2
#' @export
visualize_gamma_scalar_input <- function(gamma_fun,
                                         n = 25,
                                         arrow_color = "tomato",
                                         point_color = "grey40",
                                         warp_scale = 1,
                                         uv_grid = NULL) {
  if (is.null(uv_grid)) {
    u_seq <- seq(0, 1, length.out = n)
    v_seq <- seq(0, 1, length.out = n)
    uv_grid <- expand.grid(u = u_seq, v = v_seq)
  } else {
    uv_grid <- as.data.frame(uv_grid)
    if (ncol(uv_grid) < 2L) {
      stop("uv_grid must have at least two columns.")
    }
    uv_grid <- uv_grid[, 1:2, drop = FALSE]
    names(uv_grid) <- c("u", "v")
  }

  UV2 <- t(apply(uv_grid, 1, function(uv) gamma_fun(uv[1], uv[2])))
  UV2 <- as.data.frame(UV2)
  names(UV2) <- c("u2", "v2")

  df <- cbind(uv_grid, UV2)
  df$u2 <- df$u + (df$u2 - df$u) * warp_scale
  df$v2 <- df$v + (df$v2 - df$v) * warp_scale

  ggplot2::ggplot(df) +
    ggplot2::geom_segment(
      ggplot2::aes(x = u, y = v, xend = u2, yend = v2),
      alpha = 0.6,
      color = arrow_color,
      arrow = grid::arrow(length = grid::unit(0.08, "cm"))
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = u, y = v),
      size = 0.5,
      color = point_color
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Gamma warp (single-point input)",
      x = "u",
      y = "v"
    )
}

#' Plot a List of Simulated Surfaces
#'
#' @param data_list List of `surface_data` objects.
#' @param uv_scale Either `"global"` or `"each"`.
#' @param opacities Per-surface marker opacity.
#' @param legend_labels Per-surface legend labels.
#' @param marker_size Marker size.
#' @param title Plot title.
#' @param show_legend Logical; show plot legend.
#' @param legend_bg Legend background.
#' @param fix_axis Logical; use a shared cube-like axis range.
#'
#' @return A plotly figure.
#' @importFrom plotly plot_ly add_markers layout
#' @importFrom magrittr %>%
#' @export
plot_n_surface <- function(data_list,
                           uv_scale = c("global", "each"),
                           opacities = NULL,
                           legend_labels = NULL,
                           marker_size = 3,
                           title = "Simulated Surfaces",
                           show_legend = TRUE,
                           legend_bg = "rgba(255,255,255,0.6)",
                           fix_axis = TRUE) {
  uv_scale <- match.arg(uv_scale)

  if (!is.list(data_list) || length(data_list) < 1) {
    stop("data_list must be a non-empty list of surface objects.")
  }

  k <- length(data_list)
  if (is.null(opacities)) {
    opacities <- rep(0.85, k)
  }
  if (length(opacities) != k) {
    stop("Length of opacities must match length(data_list).")
  }

  if (is.null(legend_labels)) {
    legend_labels <- paste0("Surface ", seq_len(k))
  }
  if (length(legend_labels) != k) {
    stop("Length of legend_labels must match length(data_list).")
  }

  to_01 <- function(x) {
    x <- as.numeric(x)
    r <- range(x, na.rm = TRUE)
    if (!is.finite(r[1]) || !is.finite(r[2]) || r[1] == r[2]) {
      return(rep(0.5, length(x)))
    }
    (x - r[1]) / (r[2] - r[1])
  }

  extract_df <- function(data, idx = NA_integer_) {
    if (is.null(data$XYZ) || is.null(data$UV_grid)) {
      stop(sprintf("data_list[[%s]] must contain $XYZ and $UV_grid.", idx))
    }
    df <- cbind(data$XYZ, data$UV_grid)

    if (!all(c("u", "v") %in% colnames(df))) {
      stop(sprintf("data_list[[%s]]$UV_grid must have columns 'u' and 'v'.", idx))
    }
    if (!all(c("x", "y", "z") %in% colnames(df))) {
      stop(sprintf("data_list[[%s]]$XYZ must have columns 'x', 'y', 'z'.", idx))
    }

    df
  }

  dfs <- lapply(seq_len(k), function(i) extract_df(data_list[[i]], i))

  if (uv_scale == "each") {
    dfs <- lapply(dfs, function(df) {
      df$u01 <- to_01(df$u)
      df$v01 <- to_01(df$v)
      df
    })
  } else {
    u_all <- unlist(lapply(dfs, `[[`, "u"), use.names = FALSE)
    v_all <- unlist(lapply(dfs, `[[`, "v"), use.names = FALSE)
    u01_all <- to_01(u_all)
    v01_all <- to_01(v_all)

    lens <- vapply(dfs, nrow, integer(1))
    starts <- c(1L, 1L + cumsum(lens))[seq_len(k)]
    ends <- cumsum(lens)

    for (i in seq_len(k)) {
      dfs[[i]]$u01 <- u01_all[starts[i]:ends[i]]
      dfs[[i]]$v01 <- v01_all[starts[i]:ends[i]]
    }
  }

  colors <- lapply(dfs, function(df) {
    rgb(df$u01, df$v01, 1 - df$u01, maxColorValue = 1)
  })

  if (fix_axis) {
    x_all <- unlist(lapply(dfs, `[[`, "x"), use.names = FALSE)
    y_all <- unlist(lapply(dfs, `[[`, "y"), use.names = FALSE)
    z_all <- unlist(lapply(dfs, `[[`, "z"), use.names = FALSE)

    x_range_raw <- range(x_all, na.rm = TRUE)
    y_range_raw <- range(y_all, na.rm = TRUE)
    z_range_raw <- range(z_all, na.rm = TRUE)

    x_mid <- mean(x_range_raw)
    y_mid <- mean(y_range_raw)
    z_mid <- mean(z_range_raw)
    max_span <- max(diff(x_range_raw), diff(y_range_raw), diff(z_range_raw))

    x_range <- c(x_mid - max_span / 2, x_mid + max_span / 2)
    y_range <- c(y_mid - max_span / 2, y_mid + max_span / 2)
    z_range <- c(z_mid - max_span / 2, z_mid + max_span / 2)
  }

  p <- plotly::plot_ly()

  for (i in seq_len(k)) {
    p <- p %>%
      plotly::add_markers(
        data = dfs[[i]],
        x = ~x,
        y = ~y,
        z = ~z,
        color = I(colors[[i]]),
        marker = list(size = marker_size, opacity = opacities[i]),
        name = legend_labels[i],
        showlegend = show_legend
      )
  }

  p %>%
    plotly::layout(
      title = title,
      scene = list(
        xaxis = list(
          title = "x",
          range = if (fix_axis) x_range else NULL,
          autorange = !fix_axis
        ),
        yaxis = list(
          title = "y",
          range = if (fix_axis) y_range else NULL,
          autorange = !fix_axis
        ),
        zaxis = list(
          title = "z",
          range = if (fix_axis) z_range else NULL,
          autorange = !fix_axis
        ),
        aspectmode = "cube"
      ),
      legend = list(bgcolor = legend_bg)
    )
}
