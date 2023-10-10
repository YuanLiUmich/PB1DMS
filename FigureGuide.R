# Recommend to save figures as pdf (vectorized) then transform to tiff (not-vectorized)
# Publication parameters

fig_res <- 300
fig_wid_one_third <- 5.50
fig_wid_half <- 8.30
fig_wid_full <- 17.35
fig_wid_max <- 17.4
fig_hgt_max <- 23.0

line_wid <- 1.0
point_size <- 2

font <- "ArialMT"
fz_general <- 11
fz_panel_title <- 12
fz_axis_title <- 11
fz_axis_tick <- 10
fz_annotation <- 10
 
palette <- c("#5C4B51", "#1A3A48", "#8CBEB2", "#048789", "#F2EBBF", 
             "#F3B562", "#F06060", "#96505B", "#7E8AA2")

palette2 <- c("A" = "#8CFF8C", "G" = "#FFFFFF", "L" = "#455E45", "S" = "#FF7042",
          "V" = "#D695C2", "T" = "#B84C00", "K" = "#4747B8", "D" = "#A00042",
          "I" = "#004C00", "N" = "#FF7C70", "E" = "#660000", "P" = "#525252",
          "R" = "#00007C", "F" = "#543C42", "Q" = "#FF4C4C", "Y" = "#8C704C",
          "H" = "#7070FF", "C" = "#FFFF70", "M" = "#B8A042", "W" = "#4F4600",
          "X" = "#B8B8B8")

palette_691 <- c("F" = "#D695C2", "K" = "#7E8AA2", "N" = "#F06060", 
                 "R" = "#96505B", "X" = "#525252")

compression <- "lzw"


# # For pdf
# # The default settings can be checked with:
# pdf.options()
# # provide all device specifications and changes directly in the device setup call:
# # The default resolution of R and RStudio images are exported at 72 ppi,
# # insufficient for publication.
# pdf("some_filename.png", width = 10, height = 10, units = "cm", res = 300)
# 
# print(p_ggplot)
# invisible(dev.off())
# 
# 
# # For tiff
# tiff("test.tif", family = "Arial", units = "cm",
#     width = 17.35, height = 23, pointsize = 18, 
#     res = 300, compression = "lzw")
# # query base graphics graphical parameter text size (pointsize)
# par()$ps
# # query ggplot2 text size
# theme_bw()$text$size
# # ggplot2 ignores any parameter passed to the device via 'pointsize', to update
# tiff("test_gg.tif", width = 17.35, height = 23.35, units = "cm", res = 300,
#      compression = "lzw")
# theme_set(theme_bw(base_size = 10))
# theme_update(axis.text = element_text(size = 17.5, face = "italic"))
# 
# print(p_ggplot)
# invisible(dev.off())