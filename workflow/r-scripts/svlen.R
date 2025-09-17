library(here)
library(tidyverse)

here::i_am("r-scripts/svlen.R")
binned_data <- read.delim(here("summary_counts.tsv"))

plot_width <- 16
plot_height <- 9
options(repr.plot.width = plot_width, repr.plot.height = plot_height, scipen = 999)

theme_set(theme_bw())
truvari_folder <- here::here("pedigree-results", "truvari")
data_folder <- "type-ignored"

truvari_name <- str_replace(data_folder, "-", " ")
folders <-
  list.dirs(
    path = here::here(truvari_folder, data_folder),
    full.names = FALSE,
    recursive = TRUE
  )
print(folders)
callers <- unique(str_split_i(str_remove(folders, data_folder), "/", 1))[-1]
samples <- unique(str_split_i(str_remove(folders, data_folder), "/", 2))[-1]
print(samples)
print(callers)

save_folder <- here::here(truvari_folder, data_folder)
# quit(save="ask")

summary_file <- here::here(truvari_folder, data_folder, "data.Rds")
color_file <- here::here(truvari_folder, "color_pal.Rds")
if (!file.exists(summary_file) | !file.exists(color_file)) { # make sure these files exist
  source(here("r-scripts", "01-truvari-report.R"))
}
# load the f1 summary data
summary_data <- readRDS(file = summary_file) %>% filter(sample == "NA12878")
# get the color and orders of the callers
color_pal <- readRDS(color_file)
color_pal <- c(color_pal, "#000" = "TRUTH")
print(color_pal)

# sanity check
sanity <- binned_data %>%
  dplyr::rename(caller = tool) %>%
  group_by(caller) %>%
  summarize(FN_sanity = sum(nFN), FP_sanity = sum(nFP), TP_sanity = sum(nTP))
check_tbl <- sanity %>%
  inner_join(summary_data, by = "caller") %>%
  mutate(
    FN_match = FN_sanity == FN,
    FP_match = FP_sanity == FP,
    TP_match = TP_sanity == TP.base
  ) %>%
  filter(!(FN_match & FP_match & TP_match)) %>%
  dplyr::select(caller, TP_match, TP.base, TP_sanity, FP_match, FP, FP_sanity, FN_match, FN, FN_sanity)

# 3.  Inspect
check_tbl
unique(binned_data$SVTYPE)
binned_data <- binned_data %>%
  mutate(caller = factor(str_to_lower(tool), levels = names(color_pal), )) %>%
  mutate(SVLEN = factor(SVLEN_range, levels = c("50-100", "100-500", "500-1000", "1000+"))) %>%
  filter(!is.na(caller)) # remove anything not mapped properly
f1_all <- ggplot(binned_data, aes(x = SVLEN, y = f1, fill = caller)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = color_pal)
saveRDS(f1_all, file = here(save_folder, "fig1_f1_all.Rds"))
recall_all <- ggplot(binned_data, aes(x = SVLEN, y = recall, fill = caller)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = color_pal)
saveRDS(recall_all, file = here(save_folder, "fig1_recall_all.Rds"))

binned_data_type <- binned_data %>%
  mutate(SVTYPE = factor(case_when(
    SVTYPE == "DUP:TANDEM" ~ "INS",
    SVTYPE == "DUP:INV" ~ "INS",
    SVTYPE == "DUP" ~ "INS",
    TRUE ~ SVTYPE
  ), levels = c("DEL", "INS", "INV"))) %>%
  filter(!is.na(SVTYPE))

recall_split <- ggplot(binned_data_type, aes(x = SVLEN, y = recall, fill = caller)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = color_pal) +
  facet_wrap(~SVTYPE, scales = "fixed")
recall_split
saveRDS(recall_split, file = here(save_folder, "fig1_recall.Rds"))

library(ggh4x)
count_split_TP <- ggplot(binned_data_type, aes(x = SVLEN, y = nTP, fill = caller)) +
  geom_col(position = "dodge") +
  labs(y = "# of True Positives") +
  scale_fill_manual(values = color_pal) +
  facet_wrap(~SVTYPE, scales = "free_y", axes = "all_y") +
  scale_y_continuous(limits = c(0, 2035)) +
  scale_y_facet(PANEL == 2, labels = NULL, limits = c(0, 2035)) +
  scale_y_facet(SVTYPE == "INV", limits = c(0, 20.35))
saveRDS(count_split_TP, file = here(save_folder, "fig1_count_split_TP.Rds"))
count_split_TP
count_split_FP <- ggplot(binned_data_type, aes(x = SVLEN, y = nFP, fill = caller)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = color_pal) +
  facet_row(~SVTYPE, scales = "free")
count_split_FP

library(patchwork)
recall_split / count_split_TP + plot_layout(guides = "collect", axes = "collect")
