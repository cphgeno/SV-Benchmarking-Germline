library(tidyverse)
library(here)
library(patchwork)
library(ggh4x)

# capital letters for recall and F1

here::i_am("r-scripts/plot_combine.R")
main_folder <- here("pedigree-results", "truvari", "type-ignored")
fig1a_file <- here(main_folder, "fig1_combined_barplot_median.Rds")
fig1b_file <- here(main_folder, "plots_gt_sample.Rds")
fig1c1_file <- here(main_folder, "fig1_f1_all.Rds")
fig1c2_file <- here(main_folder, "fig1_recall_all.Rds")
fig1c3_file <- here(main_folder, "fig1_recall.Rds")
fig1c4_file <- here(main_folder, "fig1_count_split_TP.Rds")
fig1d_file <- here(main_folder, "NA12878", "chr-f1_all.Rds")


source(here("r-scripts", "01-truvari-summary.R"))
data %>%
  group_by(caller) %>%
  summarize(
    f1 = median(f1), f1_sd = sd(f1, na.rm = TRUE),
    recall = median(recall), recall_sd = sd(recall, na.rm = TRUE),
    precision = median(precision), precision_sd = sd(precision, na.rm = TRUE)
  )
gglayer_med_bar <- c(
  gglayer_bar_theme,
  stat_summary(aes(fill = caller), fun = "median", geom = "col"),
  stat_summary(fun.data = mean_se, geom = "linerange", linewidth = 1.5, alpha = 0.5)
)
f1 <- ggplot(data, aes(x = fct_relevel(caller, f1_order), y = f1)) +
  gglayer_med_bar +
  labs(y = "F1") +
  scale_y_continuous(limits = c(0, max(data$precision)))
f1
recall <- (ggplot(data, aes(x = fct_relevel(caller, f1_order), y = recall)) +
  gglayer_med_bar +
  scale_y_continuous(limits = c(0, max(data$precision)))) + labs(y = "Recall")
precision <- (ggplot(data, aes(x = fct_relevel(caller, f1_order), y = precision)) +
  gglayer_med_bar +
  labs(y = "Precision") +
  scale_y_continuous(limits = c(0, max(data$precision))))


fig1a <- f1 / recall / precision +
  plot_layout(guides = "collect", axes = "collect") &
  theme(
    legend.position = "right",
    legend.title = element_text(size = 7, margin = margin(b = 4, l = -4), face = "bold"),
    legend.text = element_text(size = 6, hjust = -1, margin = margin(l = -4, b = -2, t = 2)),
    legend.key.spacing.x = unit(0, units = "mm"),
    legend.text.position = "top",
    legend.margin = margin(l = 5, r = -25),
    # legend.box.margin=margin(0),
    legend.box.spacing = unit(1, units = "mm")
  )
fig1a
gt_table <- readRDS(file = fig1b_file) %>% opt_vertical_padding(scale = .8)
gc()

if (!file.exists(fig1c1_file) | !file.exists(fig1c2_file) | !file.exists(fig1c3_file) | !file.exists(fig1c4_file)) { # make sure these files exist
  source(here("r-scripts", "svlen.R"), local = TRUE)
}
fig1c1 <- readRDS(file = fig1c1_file) + labs(y = "F1", fill = "Caller")
fig1c2 <- readRDS(file = fig1c2_file) + labs(y = "Recall", fill = "Caller")
fig1c3 <- readRDS(file = fig1c3_file) + labs(y = " F1 ", fill = "Caller")
fig1c4 <- readRDS(file = fig1c4_file) + labs(fill = "Caller")
gc()

fig1d <- readRDS(file = fig1d_file) + guides(color = "none") + labs(title = NULL, y = "F1") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


if (!file.exists(fig1d_file)) { # make sure these files exist
  source(here("r-scripts", "02-vcf-investigation.R"), local = TRUE)
}
gc()

layout_fig1c <- c(
  area(1, 1, 1, 3), area(1, 4, 1, 6),
  area(2, 1, 2, 6),
  area(3, 1, 3, 6)
)
height_stretch <- 2
width_stretch <- 1
layout <- c(
  area(1, 1, 3 * height_stretch, 2), area(1, 3, 3 * height_stretch, 6), # row 1
  area(3 * height_stretch + 1, 1, 4 * height_stretch, 3), area(3 * height_stretch + 1, 4, 4 * height_stretch, 6),
  area(4 * height_stretch + 1, 1, 5 * height_stretch, 6),
  area(5 * height_stretch + 1, 1, 6 * height_stretch, 6),
  area(6 * height_stretch + 1, 1, 7 * height_stretch, 6),
  area(7 * height_stretch + 1, 1, height_stretch * 7 + 1, 6)
)
width <- 210
height <- 297
width_unit <- width / 6
height_unit <- height / 7

plot(layout)


fig1b <- gt_table %>% opt_horizontal_padding(scale = .9)
ggsave(fig1a & theme(legend.margin = margin(l = 5, r = 0)), filename = "fig1a.pdf", height = height_unit * 3, width = (width_unit * 2) + 25, units = "mm", dpi = 300, )
ggsave(wrap_table(fig1b, panel = "body"), filename = "fig1b.pdf")
ggsave(fig1c1 + fig1c2 + fig1c3 + fig1c4 +
  plot_layout(guides = "collect", axes = "collect", design = layout_fig1c) +
  plot_annotation(tag_levels = "A"), filename = "fig1c.pdf", height = height_unit * 3, width = width, units = "mm", dpi = 300)
ggsave(fig1d + guides(color = NULL), filename = "fig1d.pdf", height = height_unit, width = width, unit = "mm", dpi = 300)



print("Making overall figures....!")

patchwork <-
  # row 1
  wrap_elements(plot = fig1a, clip = FALSE) + wrap_table(fig1b, panel = "body") +
    fig1c1 + fig1c2 + # row 2 (f1)
    fig1c3 + # row 3 (recall)
    fig1c4 + # row 4 (count)
    fig1d + # row 5 (chromosome)
    guide_area() +
    plot_layout(guides = "collect", axes = "collect", design = layout) +
    plot_annotation(tag_levels = "A") &
    guides(fill = guide_legend(nrow = 2, ncol = 7, byrow = TRUE)) & # legend
    theme(legend.position = "bottom", legend.byrow = TRUE, legend.direction = "horizontal")


patchwork[[1]] <- patchwork[[1]] + theme(legend.position = "right") + plot_layout(guides = "keep")
# patchwork
ggsave(patchwork,
  filename = "fig1_overall.png",
  width = width, height = height,
  units = "mm", dpi = 300
)

fig1b <- gt_table %>% opt_horizontal_padding(scale = .8)

patchwork <-
  # row 1
  wrap_elements(plot = fig1a, clip = FALSE) + wrap_table(fig1b, panel = "body") +
    fig1c1 + fig1c2 + # row 2 (f1)
    fig1c3 + # row 3 (recall)
    fig1c4 + # row 4 (count)
    fig1d + # row 5 (chromosome)
    guide_area() +
    plot_layout(guides = "collect", axes = "collect", design = layout) +
    plot_annotation(tag_levels = "A") &
    guides(fill = guide_legend(nrow = 2, ncol = 7, byrow = TRUE)) & # legend
    theme(legend.position = "bottom", legend.byrow = TRUE, legend.direction = "horizontal")
patchwork[[1]] <- patchwork[[1]] + theme(legend.position = "right") + plot_layout(guides = "keep")
ggsave(patchwork,
  filename = "fig1_overall.pdf",
  width = width, height = height,
  units = "mm", dpi = 300
)
patchwork
