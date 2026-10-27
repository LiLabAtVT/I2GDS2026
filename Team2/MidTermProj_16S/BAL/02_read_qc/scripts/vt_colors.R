# Virginia Tech brand colors, shared by all plotting scripts.
# Source: https://www.assets.cms.vt.edu/brand/brand_color_sheet.txt (brand.vt.edu color tool)
# Other scripts load this file with source(), so every figure uses the same
# colours and look. To use your own colours, change the hex codes below.
suppressPackageStartupMessages(library(ggplot2))

# Named colour palette (hex codes)
vt <- c(
  maroon    = "#861F41",  # Chicago Maroon (primary)
  orange    = "#E5751F",  # Burnt Orange (primary)
  stone     = "#75787B",  # Hokie Stone (primary)
  purple    = "#642667",  # Pylon Purple
  pink      = "#CE0058",  # Boundless Pink
  yellow    = "#F7EA48",  # Triumphant Yellow
  teal      = "#508590",  # Sustainable Teal
  turquoise = "#2CD5C4",  # Vibrant Turquoise
  grey      = "#D7D2CB",  # Land Grant Grey
  smoke     = "#E5E1E6"   # Skipper Smoke
)

# Fixed colour for each group, so canker is always maroon and healthy orange
status_cols <- c(Canker = vt[["maroon"]], Healthy = vt[["orange"]])

# Fixed colour for each cultivar; the order here is also the order used in plots
cultivar_cols <- c(
  "Aztec Fuji"       = vt[["maroon"]],
  "Cripps Pink"      = vt[["orange"]],
  "Gala"             = vt[["teal"]],
  "Granniwinkle"     = vt[["purple"]],
  "Hewe's Crab"      = vt[["turquoise"]],
  "Unknown cultivar" = vt[["stone"]]
)

# Sequential (white -> maroon) and diverging (orange -> white -> maroon) fills
# (used for heatmaps in later steps)
scale_fill_vt_seq <- function(...) scale_fill_gradient(low = "#FFFFFF", high = vt[["maroon"]], ...)
scale_fill_vt_div <- function(...) scale_fill_gradient2(low = vt[["orange"]], mid = "#FFFFFF",
                                                        high = vt[["maroon"]], midpoint = 0, ...)

# Shared plot style: clean white background, no minor grid lines,
# bold maroon titles, legend at the bottom
theme_vt <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = vt[["smoke"]], colour = NA),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold", colour = vt[["maroon"]]),
          legend.position = "bottom")
}

# Save a figure as vector PDF (for print/reports) and 300-dpi PNG (for slides/docs).
# PNGs go in results/png/ alongside results/pdf/.
save_fig <- function(p, pdf_path, w = 8, h = 5) {
  ggsave(pdf_path, p, width = w, height = h)
  png_dir <- file.path(dirname(dirname(pdf_path)), "png")
  dir.create(png_dir, showWarnings = FALSE)
  ggsave(file.path(png_dir, sub(".pdf", ".png", basename(pdf_path), fixed = TRUE)), p,
         width = w, height = h, dpi = 300, bg = "white")
}
