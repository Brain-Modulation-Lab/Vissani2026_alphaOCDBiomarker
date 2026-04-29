library(svglite)
# =============================================================================
# 1. I/O FUNCTIONS
# =============================================================================

FONT_FAMILY    <- "Helvetica"
load_file <- function(path) {
  if (!file.exists(path)) stop(paste("File not found:", path), call. = FALSE)
  ext <- tolower(tools::file_ext(path))
  df <- switch(ext,
               "csv"  = readr::read_csv(path, show_col_types = FALSE),
               "xlsx" = readxl::read_excel(path),
               "xls"  = readxl::read_excel(path),
               stop("Unsupported extension"))
  return(df)
}


# 1. SVG: Using a simplified call that is version-independent
# We wrap svglite to ensure ggsave's extra arguments (like bg) don't crash it
save_fig <- function(plot, name, path, w = 10, h = 8) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE)
    
    # 1. SVG: Using svglite directly to avoid X11 module calls
    ggsave(
      filename = file.path(path, paste0(name, ".svg")),
      plot = plot,
      device = function(file, width, height, ...) {
        svglite::svglite(file, width = width, height = height, ...)
      },
      width = w, height = h, units = "in"
    )
    
    # 2. PDF: Using cairo_pdf (since capabilities confirmed it's TRUE)
    # This is much better for Illustrator than standard pdf()
    ggsave(
      filename = file.path(path, paste0(name, ".pdf")), 
      plot = plot, 
      device = grDevices::cairo_pdf, 
      width = w, height = h, units = "in"
    )
    
    # 3. PNG: Standard high-res
    ggsave(
      filename = file.path(path, paste0(name, ".png")), 
      plot = plot, 
      dpi = 300, width = w, height = h, units = "in"
    )
  }
# =============================================================================
# 2. THEMES
# =============================================================================

COW_THEME <- theme_cowplot(font_size = 12, font_family = FONT_FAMILY)

MATLAB_THEME <- theme_cowplot(font_size = 12, font_family = FONT_FAMILY) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.8), 
    axis.ticks.length = unit(0.3, "cm"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5)
  )

# =============================================================================
# 2. COMPTUATION METHODS
# =============================================================================

safe_max  <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)


get_best_session <- function(session_ids, days, scores) {
  eligible_indices <- which(days > 0)
  if (length(eligible_indices) == 0) return(NA_real_)
  best_idx <- eligible_indices[which.min(scores[eligible_indices])]
  return(session_ids[best_idx])
}

spearman_boot_fn <- function(data, indices) {
  d <- data[indices, ]
  if(nrow(d) < 3) return(NA)
  cor(d$x, d$y, method = "spearman", use = "complete.obs")
}