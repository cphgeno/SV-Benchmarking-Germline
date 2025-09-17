#!/usr/bin/env Rscript
library(tidyverse) 
library(VariantAnnotation) 
library(ggrepel)
library(here)
library("optparse")
library(patchwork)
library(janitor)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(Gviz)
library(data.table)



# options(error=traceback)
option_list <- list(
  make_option(c("-i", "--input"),
    type = "character", default = "pedigree-results",
    help = "directory with results [default = %default", metavar = "character"
  ),
  make_option(c("-r", "--rds"),
    type = "character", default = "FALSE",
    help = "if you want to reuse or name [default = %default", metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

here::i_am("workflow/r-scripts/02-vcf-investigation.R")
truvari_folder <- here(opt$input, "truvari")

# load vcf processing functions
source(here("workflow","r-scripts", "00-vcf_process.R"))

args <- commandArgs(trailingOnly = TRUE)
# useRDS <- opt$rds
useRDS <- TRUE
options(repr.plot.width = 16, repr.plot.height = 9)

theme_set(theme_bw())


# Loading Data ------------------------------------------------------------



data_folder <-
  "type-ignored/"
truvari_name <- str_replace(str_remove(data_folder, "/"), "-", " ")
folders <-
  list.dirs(
    path = here(truvari_folder, data_folder),
    full.names = FALSE,
    recursive = TRUE
  )
folders <- folders[!grepl("phab_bench", folders)] # get a list of all the summary files
folders <- folders[!grepl("temp", folders)] # get a list of all the summary files
print(folders)
callers <- unique(str_split_i(str_remove(folders, data_folder), "/", 1))[-1]
samples <- unique(str_split_i(str_remove(folders, data_folder), "/", 2))[-1]
print(samples)
print(callers)
tp_type <- "tp-base"
what <- c(tp_type, "fp", "fn")
sample_name <- samples[1]

# quit(save="ask")

summary_file <- here(truvari_folder, data_folder, "data.Rds")
color_file <- here(truvari_folder, "color_pal.Rds")
if (!file.exists(summary_file) | !file.exists(color_file)) { # make sure these files exist
  source(here("r-scripts", "01-truvari-report.R"))
}
# load the f1 summary data
summary_data <- readRDS(file = summary_file)
# get the color and orders of the callers
color_pal <- readRDS(color_file)
print(color_pal)

samples <- c("NA12878")

caller_i <- "manta"
sample_i <- sample_name
overall_data <- data.frame()
rds_location <- here(truvari_folder, data_folder, tp_type, "vcf_data.Rds")
if (useRDS == "TRUE" & file.exists(rds_location)) { # use that hard won vcf data :D
  print("LOADING OLD RDS DATA.....")
  overall_data <- readRDS(file = rds_location)
  print(head(overall_data, 1))
} else {
  for (caller_i in callers) {
    for (sample_i in samples) {
      print(caller_i)
      vcfs <- vcf_read(sample_i, caller_i, include = str_c(what, collapse = "|"))
      print(str_c("VCFs loaded for sample ", sample_i, " and caller ", caller_i))
      if (length(vcfs) == length(what)) { # make sure we're loading all 3

        for (i in seq(1, length(vcfs))) {
          print(str_c("Processing ", names(vcfs)[i]))
          data <- vcf_process(vcfs[[i]], status = names(vcfs)[i]) |>
            mutate(sample = sample_i) |>
            mutate(caller = caller_i)

          if (length(overall_data) > 1) {
            data <- data |> dplyr::select(any_of(names(overall_data)))
            print(head(data, 1))
          }
          overall_data <- bind_rows(overall_data, data) 
          overall_data <- overall_data[rowSums(is.na(overall_data)) != ncol(overall_data),,drop=FALSE]
          overall_data <- overall_data[,colSums(is.na(overall_data)) != nrow(overall_data),drop=FALSE]
        }
      }
    }
  }

  # save that hard won vcf data :)
  print("Saving RDS data...")
  print(head(overall_data, 1))
  saveRDS(overall_data, file = rds_location)
}



# based on hue_pal()2, then grey
group.colors <- c(
  "True Positive" = "#00BFC4",
  "False Positive" = "#F8766D",
  "False Negative" = "grey"
)
names(overall_data)
overall_data <- dplyr::rename(overall_data, c("type" = "SVTYPE"))


if (length(unique(overall_data$caller) != length(names(color_pal)))) {
  print(str_c("We are missing ", unique(overall_data$caller)[!(unique(overall_data$caller) %in% names(color_pal))], " in the summary data"))
  print(str_c("We are missing ", names(color_pal)[!(names(color_pal) %in% unique(overall_data$caller))], " in the vcfs"))
}

overall_data <- overall_data %>%
  mutate(sample = as_factor(overall_data$sample)) %>%
  mutate(caller = factor(str_to_lower(caller), levels = names(color_pal), )) %>%
  filter(!is.na(caller)) # remove anything not mapped properly
print(levels(overall_data$sample))

# remove extra typing
overall_data <- overall_data |>
  mutate(type = case_when(caller == "tardis" ~ str_split_i(type, ":", 1), type == "DUP:TANDEM" ~ "DUP", TRUE ~ type))

### 0 Sanity Check
print("Running sanity checks for read in data....")
sanity <- overall_data %>%
  group_by(caller, sample, status) |>
  dplyr::count() |>
  pivot_wider(
    names_from = status,
    values_from = n,
    values_fill = 0,
  ) |>
  clean_names() |>
  full_join(
    summary_data |>
      dplyr::select(caller, sample,
        tp_comp = TP.comp,
        fn = FN, fp = FP, f1, precision, recall
      ),
    by = c("caller", "sample"), suffix = c(".vcf", ".json")
  ) |>
  mutate(across(fn.vcf:f1, ~ ifelse(is.na(.x), 0, .x)))

# print("Removing any caller sample set found here:")
# sanity[!(sanity$tp_comp.json == sanity$tp_comp.vcf), ]
# sanity[!(sanity$fn.json == sanity$fn.vcf), ]
# sanity[!(sanity$fp.json == sanity$fp.vcf), ]


overall_data <- overall_data |>
  #  anti_join(sanity[!(sanity$fn.json == sanity$fn.vcf), ], by = c("caller", "sample")) |>
  #  anti_join(sanity[!(sanity$fp.json == sanity$fp.vcf), ], by = c("caller", "sample")) |>
  #  anti_join(sanity[!(sanity$tp_comp.json == sanity$tp_comp.vcf), ], by = c("caller", "sample")) |>
  mutate(status = case_when(status == "fn" ~ "False Negative", status == "fp" ~ "False Positive", status == tp_type ~ "True Positive"))

head(overall_data)


# Type Counts -------------------------------------------------------------


type_counts <- overall_data %>%
  group_by(sample, caller, status, type) %>%
  dplyr::count() |>
  mutate(caller = factor(caller, levels = unique(overall_data$caller))) |>
  mutate(n_label = case_when(n > 1000 ~ NA, TRUE ~ n))

print(levels(type_counts$caller))

sample_names <- unique(overall_data$sample)
for (sample_name in sample_names) {
  type_count_graph <- ggplot(
    type_counts %>% filter(sample == sample_name),
    aes(x = reorder(type, -n), y = n, fill = status)
  ) +
    geom_col() +
    scale_fill_manual(values = group.colors) +
    scale_color_manual(values = group.colors) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(title = str_c("Types called and missed, with ", truvari_name, " on sample ", sample_name), x = "Type", y = "Counts")


  ggplot2::ggsave(create.dir = TRUE, width = 16, height = 9, here(truvari_folder, data_folder, tp_type, sample_name, "type_counts_raw.png"))

  type_count_graph + facet_wrap(vars(caller))
  ggplot2::ggsave(create.dir = TRUE, width = 16, height = 9, here(truvari_folder, data_folder, tp_type, sample_name, "type_counts.png"))

  type_count_graph + facet_wrap(vars(caller), scales = "free_x") +
    geom_label_repel(max.overlaps = 50, nudge_y = 100, aes(label = n_label, color = status), fill = "white", direction = "y")
  ggplot2::ggsave(create.dir = TRUE, width = 16, height = 9, here(truvari_folder, data_folder, tp_type, sample_name, "type_counts_freex_label.png"))


  type_count_graph + facet_wrap(vars(caller), scales = "free_x")
  ggplot2::ggsave(create.dir = TRUE, width = 16, height = 9, here(truvari_folder, data_folder, tp_type, sample_name, "type_counts_freexy.png"))

  type_count_graph + facet_wrap(vars(caller), scales = "free_x") + geom_label_repel(max.overlaps = 50, nudge_y = 100, aes(label = n_label, color = status), fill = "white", direction = "y")
  ggplot2::ggsave(create.dir = TRUE, width = 16, height = 9, here(truvari_folder, data_folder, tp_type, sample_name, "type_counts_freexy_label.png"))
  print(sample_names)
  type_count_graph <- ggplot(
    type_counts %>% filter(sample == sample_name) %>%
      filter(status != "False Negative"),
    aes(x = reorder(type, -n), y = n, fill = status)
  ) +
    geom_col() +
    scale_fill_manual(values = group.colors) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    facet_wrap(vars(caller), scales = "free_x") +
    labs(title = str_c("Types called with ", truvari_name, " on sample ", sample_name), x = "Type", y = "Counts")
  type_count_graph
  ggplot2::ggsave(create.dir = TRUE, width = 16, height = 9, here(truvari_folder, data_folder, tp_type, sample_name, "type_counts.png"))

  type_count_graph +
    scale_color_manual(values = group.colors) +
    geom_label_repel(max.overlaps = 50, nudge_y = 100, aes(label = n_label, color = status), fill = "white", direction = "y")
  ggplot2::ggsave(create.dir = TRUE, width = 16, height = 9, here(truvari_folder, data_folder, tp_type, sample_name, "type_counts_label.png"))
}


# F1 By Type --------------------------------------------------------------

print("== 2 - F1 BY TYPE")
f1_stats <- type_counts |>
  dplyr::select(-n_label) |>
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  pivot_wider(
    names_from = status,
    values_from = n,
    values_fill = 0
  ) %>%
  clean_names() %>%
  mutate(precision = true_positive / (true_positive + false_positive)) %>%
  mutate(recall = true_positive / (true_positive + false_negative)) %>%
  mutate(f1 = 2 * ((recall * precision)) / (recall + precision))


ggplot(f1_stats, aes(x = reorder(type, -f1), y = f1, fill = caller)) +
  scale_fill_manual(values = color_pal, guide = "none") +
  stat_summary(fun = "median", fun.args = list(na.rm = TRUE), geom = "col") +
  facet_wrap(~caller, scale = "free_x") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  geom_point(aes(shape = sample))

ggplot2::ggsave(create.dir = TRUE, width = 16, height = 9, here(truvari_folder, data_folder, tp_type, "type_f1.png"))

ggplot(f1_stats, aes(x = reorder(type, -precision), y = precision, fill = caller)) +
  scale_fill_manual(values = color_pal, guide = "none") +
  stat_summary(fun = "median", fun.args = list(na.rm = TRUE), geom = "col") +
  facet_wrap(~caller, scale = "free_x") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  geom_point(aes(shape = sample))
ggplot2::ggsave(create.dir = TRUE, width = 16, height = 9, here(truvari_folder, data_folder, tp_type, "type_precision.png"))

ggplot(f1_stats, aes(x = reorder(type, -recall), y = recall, fill = caller)) +
  scale_fill_manual(values = color_pal, guide = "none") +
  stat_summary(fun = "median", fun.args = list(na.rm = TRUE), geom = "col") +
  facet_wrap(~caller, scale = "free_x") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  geom_point(aes(shape = sample))

ggplot2::ggsave(create.dir = TRUE, width = 16, height = 9, here(truvari_folder, data_folder, tp_type, "type_recall.png"))




# Overall Counts ----------------------------------------------------------


counts <- overall_data %>%
  group_by(sample, caller, status) %>%
  dplyr::count()

ggplot(
  counts,
  aes(
    x = fct_reorder2(caller, status, n),
    y = n,
    fill = status
  )
) +
  geom_col() +
  facet_grid(. ~ sample) +
  scale_fill_manual(values = group.colors) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggplot2::ggsave(create.dir = TRUE, width = 16, height = 9, here(truvari_folder, data_folder, tp_type, "overall-positive.png"))


for (sample_name in sample_names) {
  ggplot(
    counts %>% filter(sample == sample_name),
    aes(
      x = fct_reorder2(caller, status, n),
      y = n,
      fill = status
    )
  ) +
    geom_col() +
    scale_fill_manual(values = group.colors) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(x = "Callers", y = "Count", title = str_c("Number of Variants called in Sample ", sample_name), subtitle = "Status based on Truvari run with type ignored")
  ggplot2::ggsave(create.dir = TRUE, width = 16, height = 9, here(truvari_folder, data_folder, tp_type, sample_name, "overall-positive.png"))
}


# Extra Info --------------------------------------------------------------

## Looking at Width
overall_data <- overall_data |>
  mutate(width = as.integer(width)) |>
  mutate(binned_width = cut_interval(width, length = 30))
overall_data <- overall_data |> mutate(type_cat = (case_when(
  str_detect(type, "INS") ~ "INS",
  str_detect(type, "DEL") ~ "DEL",
  str_detect(type, "DUP") ~ "INS",
  TRUE ~ "OTHER"
)))
ggplot(overall_data |> filter(type_cat != "INS"), aes(x = width, y = caller)) +
  geom_boxplot(outliers = FALSE) +
  facet_grid(status ~ type_cat)
ggplot2::ggsave(create.dir = TRUE, width = 16, height = 9, here(truvari_folder, data_folder, tp_type, "basic_width_noins.png"))
ggplot(overall_data, aes(x = width, y = caller)) +
  geom_boxplot(outliers = FALSE) +
  facet_grid(status ~ type_cat)
ggplot2::ggsave(create.dir = TRUE, width = 16, height = 9, here(truvari_folder, data_folder, tp_type, "basic_width_withins.png"))
## Looking at matching stats

### PctSizeSimilarity
ggplot(overall_data, aes(x = caller, y = PctSizeSimilarity)) +
  # geom_hline(yintercept = 0.7, cor = "red") +
  geom_jitter(aes(color = caller), size = 0.1, alpha = 0.4) +
  geom_violin(outliers = FALSE, draw_quantiles = c(.25, .5, .75), alpha = 0.7) +
  scale_color_manual(values = color_pal, guide = "none") +
  facet_wrap(status ~ .) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(title = "Pct size similarity between this variant and its closest match")
ggplot2::ggsave(create.dir = TRUE, width = 16, height = 9, here(truvari_folder, data_folder, tp_type, "pctsizesimilarity_violin.png"))



# By Chromosome -----------------------------------------------------------

top_5_callers <-
  (summary_data %>% group_by(caller) %>% summarise(f1 = median(f1)) %>% arrange(-f1))$caller[1:5]

## Binned Location Scores
for (sample_name in sample_names) {
  vcfs_smaller <-
    overall_data %>% filter(caller %in% top_5_callers, sample == sample_name)

  chr_stats <-
    vcfs_smaller %>%
    group_by(seqnames, sample, caller, status) %>%
    dplyr::count() %>%
    mutate(n = ifelse(is.na(n), 0, n)) %>%
    pivot_wider(names_from = status, values_from = n, values_fill = 0) %>%
    clean_names() %>%
    mutate(precision = true_positive / (true_positive + false_positive)) %>%
    mutate(recall = true_positive / (true_positive + false_negative)) %>%
    mutate(f1 = 2 * ((recall * precision)) / (recall + precision)) %>%
    na.omit()

  print(head(chr_stats))
  print(chr_stats)

  ggplot(chr_stats |> filter(sample == sample_name), aes(
    x = seqnames,
    y = f1,
    color = caller,
    group = caller
  )) +
    scale_color_manual(values = color_pal) +
    geom_line() +
    geom_point() +
    labs(
      title = str_c("F1 Scores Across Chromosomes in sample ", sample_name),
      subtitle = "using only the top 5 callers based on F1 score",
      x = "Chromosome", y = "F1"
    )
  ggplot2::ggsave(create.dir = TRUE, width = 16, height = 9, here(truvari_folder, data_folder, tp_type, sample_name, "chr-f1.png"))

  ggplot(chr_stats |> filter(sample == sample_name), aes(
    x = seqnames,
    y = precision,
    color = caller,
    group = caller
  )) +
    scale_color_manual(values = color_pal) +
    geom_line() +
    geom_point() +
    labs(
      title = str_c("Precision Across Chromosomes in sample ", sample_name),
      subtitle = "using only the top 5 callers based on F1 score",
      x = "Chromosome"
    )
  ggplot2::ggsave(create.dir = TRUE, width = 16, height = 9, here(truvari_folder, data_folder, tp_type, sample_name, "chr-pr.png"))

  ggplot(chr_stats |> filter(sample == sample_name), aes(
    x = seqnames,
    y = recall,
    color = caller,
    group = caller
  )) +
    scale_color_manual(values = color_pal) +
    geom_line() +
    geom_point() +
    labs(
      title = str_c("Recall Across Chromosomes in sample ", sample_name),
      subtitle = "using only the top 5 callers based on F1 score",
      x = "Chromosome"
    )
  ggplot2::ggsave(create.dir = TRUE, width = 16, height = 9, here(truvari_folder, data_folder, tp_type, sample_name, "chr-recall.png"))

  chr_stats_larger <-
    overall_data %>%
    group_by(seqnames, sample, caller, status) %>%
    dplyr::count() %>%
    mutate(n = ifelse(is.na(n), 0, n)) %>%
    pivot_wider(names_from = status, values_from = n, values_fill = 0) %>%
    clean_names() %>%
    mutate(precision = true_positive / (true_positive + false_positive)) %>%
    mutate(recall = true_positive / (true_positive + false_negative)) %>%
    mutate(f1 = 2 * ((recall * precision)) / (recall + precision)) %>%
    na.omit() %>%
    filter(sample == sample)

  chr_large <- ggplot(chr_stats_larger |> filter(sample == sample_name), aes(
    x = seqnames,
    y = f1,
    color = caller,
    group = caller
  )) +
    scale_color_manual(values = color_pal) +
    geom_line() +
    geom_point() +
    labs(title = str_c("F1 Scores Across Chromosomes in sample ", sample_name, " with ", truvari_name), x = "Chromosome")
  ggplot2::ggsave(chr_large, create.dir = TRUE, width = 16, height = 9, filename = here(truvari_folder, data_folder, tp_type, sample_name, "chr-f1_all.png"))
  saveRDS(chr_large, here(truvari_folder, data_folder, tp_type, sample_name, "chr-f1_all.Rds"))

  ggplot(chr_stats_larger |> filter(sample == sample_name), aes(
    x = seqnames,
    y = precision,
    color = caller,
    group = caller
  )) +
    scale_color_manual(values = color_pal) +
    geom_line() +
    geom_point() +
    labs(title = str_c("F1 Scores Across Chromosomes in sample ", sample_name), x = "Chromosome")
  ggplot2::ggsave(create.dir = TRUE, width = 16, height = 9, here(truvari_folder, data_folder, tp_type, sample_name, "chr-precision_all.png"))

  ggplot(chr_stats_larger |> filter(sample == sample_name), aes(
    x = seqnames,
    y = recall,
    color = caller,
    group = caller
  )) +
    scale_color_manual(values = color_pal) +
    geom_line() +
    geom_point() +
    labs(title = str_c("F1 Scores Across Chromosomes in sample ", sample_name), x = "Chromosome")
  ggplot2::ggsave(create.dir = TRUE, width = 16, height = 9, here(truvari_folder, data_folder, tp_type, sample_name, "chr-recall_all.png"))
}


# Location  ---------------------------------------------------------------
## Binning  ---------------------------------------------------------------

info <- scanVcfHeader(file = "pedigree-results/truvari/type-ignored/cnvpytor/NA12878/tp-comp.vcf.gz")

chr_len <- info@header$contig
chr_len <- data.frame(chrom = rownames(chr_len), size = as.numeric(chr_len$length)) %>%
  mutate(chrom = factor(chrom, levels = c(paste0("chr", 1:22), "chrX", "chrY", "chrM"))) %>%
  filter(!is.na(chrom)) %>%
  arrange(chrom)
chr_len$sums <- 0
for (i in 2:length(chr_len$sums)) {
  chr_len$sums[i] <- sum(chr_len$size[1:i - 1])
}
chr_sum_list <- chr_len$sums
chr_len_list <- chr_len$size
names(chr_sum_list) <- names(chr_len_list) <- chr_len$chrom
overall_data <- overall_data %>% mutate(adjusted_start = start + chr_sum_list[seqnames])
cut_width <- 5000000

find_start_buckets <- function(chr_len, cut_width = 30000000) {
  buckets <- c()
  for (i in 1:25) {
    print(chr_len[i, ])
    start <- chr_len[i, 3]
    end <- start + chr_len[i, 2]
    new_buckets <- c()
    print(paste("The end is:", end))
    for (j in 1:round(chr_len[i, 2] / cut_width)) {
      new_max <- start + (j * cut_width)
      new_max <- ifelse(new_max > end, end, new_max)
      new_buckets <- c(new_buckets, new_max)
    }
    print(paste("j was", j))
    new_buckets <- unique(new_buckets)
    print(new_buckets)
    buckets <- c(buckets, new_buckets)
  }
  buckets <- unique(buckets)
  return(buckets)
}
overall_data <- overall_data %>% filter(!is.na(adjusted_start))
thresholds <- c(1, find_start_buckets(chr_len, cut_width = cut_width)) %>% sort()
idx <- findInterval(x = overall_data$adjusted_start, vec = thresholds)
overall_data <- overall_data %>% mutate(
  adjusted_start_bin = thresholds[idx], # shift by 1; 0 → NA
  adjusted_end_bin = thresholds[idx + 1],
  start_bin = thresholds[idx] - chr_sum_list[seqnames],
  end_bin = thresholds[idx + 1] - chr_sum_list[seqnames]
)
unique(overall_data$adjusted_start_bin)
unique(overall_data$adjusted_end_bin)


chr_stats_larger <-
  overall_data %>%
  mutate(seqnames = factor(seqnames, levels = c(paste0("chr", 1:22), "chrX", "chrY", "chrM"))) %>%
  filter(!is.na(seqnames)) %>%
  group_by(seqnames, adjusted_start_bin, adjusted_end_bin, start_bin, end_bin, sample, caller, status) %>%
  dplyr::count() %>%
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  pivot_wider(names_from = status, values_from = n, values_fill = 0) %>%
  clean_names() %>%
  mutate(precision = true_positive / (true_positive + false_positive)) %>%
  mutate(recall = true_positive / (true_positive + false_negative)) %>%
  mutate(f1 = 2 * ((recall * precision)) / (recall + precision)) %>%
  na.omit() %>%
  filter(sample == "NA12878") %>%
  mutate(
    start_bin = as.integer(start_bin),
    end_bin   = as.integer(end_bin)
  ) %>%
  # Drop anything with missing or nonsensical coordinates
  filter(
    !is.na(seqnames),
    !is.na(start_bin),
    !is.na(end_bin),
    start_bin > 0,
    end_bin > 0,
    start_bin <= 2147483647, # stay inside int32
    end_bin <= 2147483647,
    start_bin <= end_bin # sanity check
  )

library(zoo)

chr_stats_larger <- chr_stats_larger %>%
  mutate(
    precision = ifelse(false_positive == 0, NA, precision),
    f1 = ifelse(false_positive == 0, NA, f1)
  )

ggplot(chr_stats_larger, aes(x = adjusted_start_bin, y = recall, fill = caller)) +
  geom_area(size = 0.8, alpha = 0.5) +
  scale_fill_manual(values = color_pal)


ggplot(chr_stats_larger, aes(x = adjusted_start_bin, y = recall, color = caller)) +
  geom_line(size = 0.8, alpha = 0.8) +
  scale_color_manual(values = color_pal)

## GVIS by Chromosome  ---------------------------------------------------------------

session <- browserSession("UCSC")
genome(session) <- "hg38"
black_tab <- "encBlacklist"

if (!exists("ENCODE_BLACKLIST")) {
  query <- ucscTableQuery(session, table = black_tab)
  ENCODE_BLACKLIST <- getTable(query)
}
ENCODE_BLACKLIST <- ENCODE_BLACKLIST %>% mutate(fill = ifelse(name == "Low Mappability", "#0000ff20", "#ff000020"))

plot_chromosome <- function(df,
                            value_var = "recall",
                            chr = "chr1",
                            genome = "hg38",
                            color_pal = NULL,
                            show_plot = TRUE,
                            plot_type = "s",
                            line_alpha = 0.7, add = FALSE,
                            show_title = TRUE, legend = TRUE,
                            xmin = NULL, xmax = NULL) {
  ## ---- Safety checks ----------------------------------------------------
  req_cols <- c("seqnames", "start_bin", "end_bin", "caller", value_var)
  stopifnot(all(req_cols %in% names(df)))


  ## viewport limits (same as ‘from’ and ‘to’ in plotTracks)
  if (is.null(xmin)) {
    xmin <- 0
  }
  if (is.null(xmax)) {
    xmax <- chr_len_list[chr]
  }

  if (!is.null(color_pal)) {
    stopifnot(
      is.character(color_pal),
      !is.null(names(color_pal)),
      all(nzchar(names(color_pal)))
    )
  }

  ## helper: add/replace alpha in a hex colour -----------------------------
  add_alpha <- function(hex, alpha = 1) {
    alpha_hex <- sprintf("%02X", round(alpha * 255))
    sub("^#([0-9A-Fa-f]{6})([0-9A-Fa-f]{2})?$", paste0("#\\1", alpha_hex), hex)
  }

  ## ---- Data wrangling ---------------------------------------------------


  wide <- df %>%
    filter(start_bin > xmin, end_bin < xmax) %>%
    filter(seqnames == chr) %>%
    ungroup() %>%
    dplyr::select(
      start = start_bin, end = end_bin,
      caller, value = all_of(value_var)
    ) %>%
    pivot_wider(names_from = caller, values_from = value) %>%
    arrange(start)

  ## Build GRanges (one row per bin)
  gr <- GRanges(chr, IRanges(wide$start, wide$end), genome = genome)
  genome(gr) <- genome

  ## Matrix: rows = callers, cols = bins
  data_mat <- t(as.matrix(wide %>% dplyr::select(-start, -end)))
  callers <- rownames(data_mat) # callers present in data

  ## ---- Colour handling & row re-ordering --------------------------------
  if (is.null(color_pal)) {
    # default palette
    col_vec <- RColorBrewer::brewer.pal(length(callers), "Set1")
    row_order <- callers # keep original order
  } else {
    # callers that have colours AND are in the data — keep their order as in the palette
    keep <- names(color_pal) %in% callers
    row_order <- names(color_pal)[keep]

    # callers present in data but missing in palette → append them
    missing_pal <- setdiff(callers, row_order)
    if (length(missing_pal)) {
      warning(
        "No colour for caller(s): ",
        paste(missing_pal, collapse = ", "),
        ". Using black."
      )
      row_order <- c(row_order, missing_pal)
      color_pal <- c(color_pal, setNames(
        rep("#000000", length(missing_pal)),
        missing_pal
      ))
    }

    # callers that had colours but are absent from data
    extra <- names(color_pal)[!names(color_pal) %in% callers]
    if (length(extra)) {
      warning(
        "Palette has colours for caller(s) not in data: ",
        paste(extra, collapse = ", ")
      )
    }

    # re-order data matrix rows and colour vector
    data_mat <- data_mat[row_order, , drop = FALSE]
    col_vec <- color_pal[row_order]
  }
  col_vec <- add_alpha(col_vec, line_alpha)




  ## ---- DataTrack --------------------------------------------------------
  dtrack <- DataTrack(
    range      = gr,
    data       = data_mat,
    groups     = factor(row_order, levels = row_order), # keep order
    type       = plot_type,
    name       = tools::toTitleCase(value_var),
    genome     = genome,
    chromosome = chr,
    legend     = legend
  )
  displayPars(dtrack) <- list(col = col_vec, lwd = 2)
  ## ---- Plot -------------------------------------------------------------
  itrack <- IdeogramTrack(genome = genome, chromosome = chr, name = chr)
  gtrack <- GenomeAxisTrack()


  bl <- ENCODE_BLACKLIST[ENCODE_BLACKLIST$chrom == chr, ] %>% filter(chromStart > xmin, chromEnd < xmax)


  highlight_track <- HighlightTrack(
    trackList  = list(gtrack, dtrack), # list *inside* the highlight
    start      = bl$chromStart,
    end        = bl$chromEnd,
    chromosome = chr,
    genome     = genome,
    col        = NA, # border
    fill       = bl$fill # 20(hex) ≈ 12% opaque black
  )


  plot_title <- sprintf(
    "%s per caller on %s",
    tools::toTitleCase(value_var), chr
  )
  plot_title <- ifelse(show_title, title, "")

  plot_tracks <- plotTracks(
    list(itrack, highlight_track),
    from = xmin,
    to = xmax,
    main = plot_title,
    showId = TRUE,
    showBandId = TRUE,
    add = add
  )
  invisible(dtrack)
  return(plot_tracks)
}
grob_chrom <- function(chr, value_var = "recall") {
  print(paste("Plotting chromosome... ", chr))
  if (!chr %in% chr_stats_larger$seqnames) {
    message(sprintf("Skipping %s – no rows in input data.\n", chr))
    return(grid::rectGrob(gp = grid::gpar(col = NA))) # empty placeholder
  }
  plot_chromosome(chr_stats_larger, chr,
    value_var = value_var,
    color_pal = color_pal, line_alpha = 0.6,
    show_title = FALSE, legend = FALSE
  ) %>% grid::grid.grabExpr()
}

make_legend_grob <- function(color_pal, line_alpha = 1) {
  callers <- names(color_pal)
  cols <- sprintf(
    "%s%02X", substr(color_pal, 1, 7),
    round(line_alpha * 255)
  ) # append alpha if wanted
  n <- length(callers)

  segs <- grid::segmentsGrob(
    x0 = rep(0.02, n),
    y0 = rev(seq_len(n)) / (n + 1),
    x1 = rep(0.12, n),
    y1 = rev(seq_len(n)) / (n + 1),
    gp = grid::gpar(col = cols, lwd = 2)
  )

  labs <- grid::textGrob(
    callers,
    x = 0.15,
    y = rev(seq_len(n)) / (n + 1),
    just = "left",
    gp = grid::gpar(cex = 0.75)
  )

  grid::grobTree(grid::rectGrob(gp = grid::gpar(col = NA)), segs, labs)
}

plot_chrs <- chr_len %>% filter(chrom != "chrX", chrom != "chrY", chrom != "chrM")
legend_grob <- make_legend_grob(color_pal, line_alpha = 0.6)
legend_patch <- wrap_elements(full = legend_grob) + theme_void()
grobs <- lapply(as.character(plot_chrs$chrom), grob_chrom, value_var = "recall")
panels <- wrap_plots(c(grobs, list(legend_patch)), ncol = 4)
full_recall <- panels + plot_annotation(title = "Per-caller Recall across All Chromosomes")
ggsave(full_recall, filename = here(truvari_folder, data_folder, tp_type, "full_recall.pdf"), width = 32, height = 16)
ggsave(full_recall, filename = here(truvari_folder, data_folder, tp_type, "full_recall.png"), width = 32, height = 16)

grobs <- lapply(as.character(plot_chrs$chrom), grob_chrom, value_var = "f1")
panels <- wrap_plots(c(grobs, list(legend_patch)), ncol = 4)
full_recall <- panels + plot_annotation(title = "Per-caller F1 across All Chromosomes")
ggsave(full_recall, filename = here(truvari_folder, data_folder, tp_type, "full_f1.pdf"), width = 32, height = 16)
ggsave(full_recall, filename = here(truvari_folder, data_folder, tp_type, "full_f1.png"), width = 32, height = 16)

## Difficult Region Zoom-In  ---------------------------------------------------------------

difficult_regions <- (chr_stats_larger %>% group_by(seqnames, adjusted_start_bin, adjusted_end_bin, start_bin, end_bin) %>% summarize(recall = mean(recall)) %>% arrange(recall))

for (i in 1:20) {
  dr <- difficult_regions[i, ]
  region_name <- (paste0("dr", i, "_", dr$seqnames, "_", dr$start_bin, "-", dr$end_bin))
  print(region_name)
  dr_data <- overall_data %>%
    filter(
      start > dr$start_bin,
      end < dr$end_bin,
      seqnames == as.character(dr$seqnames)
    ) %>%
    dplyr::select(-END) %>%
    mutate(type_color = group.colors[status])
  ## Build GRanges (one row per bin)
  head(dr_data, n = 1)
  gr <- makeGRangesFromDataFrame(as.data.frame(dr_data, stringsAsFactors = FALSE),
    keep.extra.columns = TRUE,
    start.field = "start", end.field = "end",
    seqnames.field = "seqnames",
    seqinfo = chr_len_list, na.rm = TRUE
  )
  genome(gr) <- "hg38"

  itrack <- IdeogramTrack(genome = "hg38", chromosome = dr$seqnames)
  gtrack <- GenomeAxisTrack()
  track_list <- list()
  for (caller_i in names(color_pal)) {
    gr_data <- makeGRangesFromDataFrame(dr_data %>% filter(caller == caller_i),
      seqinfo = chr_len_list,
      keep.extra.columns = TRUE
    )
    gr_data$type[is.na(gr_data$type)] <- "BND"
    atrack <- AnnotationTrack(gr_data, group = gr_data$type, feature = gr_data$status, name = caller_i)
    displayPars(atrack) <- list(
      fill = group.colors, col = "#00000011",
      background.panel = adjustcolor(color_pal[caller_i], alpha.f = 0.1),
      background.title = color_pal[caller_i],
      cex.id = 0.2, just.id = "above", fontsize.id = 8
    )
    track_list[[length(track_list) + 1L]] <- atrack
  }
  xmin <- min(dr_data$start)
  xmax <- max(dr_data$end)
  bl <- ENCODE_BLACKLIST[ENCODE_BLACKLIST$chrom == dr$seqnames, ] %>% filter(chromStart <= xmax, chromEnd >= xmin)

  if (dim(bl)[1] > 0) {
    highlight_track <- HighlightTrack(
      trackList  = c(gtrack, track_list), # list *inside* the highlight
      start      = bl$chromStart,
      end        = bl$chromEnd,
      chromosome = dr$seqnames,
      genome     = "hg38",
      col        = NA, # border
      fill       = bl$fill # 20(hex) ≈ 12% opaque black
    )
  } else {
    highlight_track <- HighlightTrack(
      trackList  = c(gtrack, track_list), # list *inside* the highlight
      start      = 1,
      end        = 1,
      chromosome = dr$seqnames,
      genome     = "hg38",
      col        = NA, # border
      fill       = "#00000000" # 20(hex) ≈ 12% opaque black
    )
  }

  pdf(here(truvari_folder, data_folder, tp_type, paste0(region_name, ".pdf")), width = 16, height = 8)
  plotTracks(list(itrack, highlight_track),
    from = xmin,
    to = xmax,
    stacking = "dense",
    groupAnnotation = "group", min.width = 1, showOverplotting = TRUE
  )
  dev.off()
  png(here(truvari_folder, data_folder, tp_type, paste0(region_name, ".png")), width = 8, height = 4, units = "in", res = 300)
  plotTracks(list(itrack, highlight_track),
    from = xmin,
    to = xmax,
    stacking = "dense",
    groupAnnotation = "group", min.width = 1, showOverplotting = TRUE
  )
  dev.off()
  pdf(here(truvari_folder, data_folder, tp_type, paste0(region_name, "_EXPAND.pdf")), width = 16, height = 16)
  plotTracks(list(itrack, highlight_track),
    from = xmin,
    to = xmax,
    groupAnnotation = "group", min.width = 1, cex.feature = 0.4, just.feature = "above", showOverplotting = TRUE
  )
  dev.off()
}



# Specific SVs ------------------------------------------------------------


sv_colors <- c(
  "DEL" = "#D55E00", # Muted red
  "INS" = "#0072B2", # Muted blue
  "INV" = "#CC79A7" # Bright magenta
)
truth_only <- overall_data %>%
  filter(status != "False Positive") %>%
  group_by(seqnames, start, end, type, status, caller) %>%
  mutate(status = (status == "True Positive"))
truth_only_wide <- truth_only %>%
  pivot_wider(
    id_cols = c(seqnames, start, end, type),
    names_from = caller,
    values_from = status,
    values_fn = ~ any(.x), unused_fn = ~ unique(.x)
  ) %>%
  mutate(popdel = ifelse(is.na(popdel), FALSE, popdel)) %>%
  arrange(seqnames, start)
truth_only_wide <- truth_only_wide %>%
  mutate(n_callers = rowSums(pick(where(is.logical)))) %>%
  arrange(n_callers)


truth_only_wide <- truth_only_wide %>%
  mutate(source_1 = str_split(SOURCES_1, "[-._:]+") |> map_chr(~ str_c(.x[1:2], collapse = "_"))) %>%
  mutate(source_2 = str_split(SOURCES_2, "[-._:]+") |> map_chr(~ str_c(.x[1:2], collapse = "_"))) %>%
  mutate(source_3 = str_split(SOURCES_3, "[-._:]+") |> map_chr(~ str_c(.x[1:2], collapse = "_"))) %>%
  mutate(source_4 = str_split(SOURCES_4, "[-._:]+") |> map_chr(~ str_c(.x[1:2], collapse = "_"))) %>%
  rowwise() %>% # work row-by-row
  mutate(source = { # new combined column
    v <- c_across(starts_with("source_")) |> na.omit() # grab & drop NA
    if (length(v)) str_c(sort(v), collapse = ", ") else NA_character_
  }) %>%
  ungroup()

unique(truth_only_wide$source)

ratios <- truth_only_wide %>%
  group_by(n_callers) %>%
  summarise(n = n())
print(paste0(round((ratios[1, ]$n / sum(ratios$n)) * 100, 2), "% of calls are missed"))
ggplot(truth_only_wide, aes(x = start, y = n_callers)) +
  facet_wrap(vars(seqnames), scale = "free_x") +
  stat_bin_hex(bins = 50)
ggsave("wow.pdf", width = 64, height = 32, limitsize = FALSE)

ggplot(truth_only_wide %>% mutate(source = forcats::fct_reorder(source, n_callers, .fun = mean, .desc = TRUE)), aes(x = source, y = n_callers)) +
  geom_jitter(aes(color = type), size = 0.8, alpha = 0.3) +
  scale_color_manual(values = sv_colors) +
  geom_violin() +
  labs(y = "Number of Callers that Found It", y = "Source from Truth Set") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave("source.pdf", width = 16, height = 9, limitsize = FALSE)
