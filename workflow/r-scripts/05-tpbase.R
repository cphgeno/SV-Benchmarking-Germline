#!/usrbin/env Rscript
library(janitor)
library(ggrepel)
library(ggforce)
library(here)
library(optparse)

option_list <- list(
  make_option(c("-R", "--rds"),
    action = "store_true", default = FALSE,
    help = "Reuse an existing RDS [default]"
  ),
  make_option(c("-i", "--input"),
    type = "character", default = "pedigree-results",
    help = "The folder containing the results [default %default]"
  ),
  make_option(c("-t", "--truvari"),
    type = "character", default = "type-ignored",
    help = "The truvari run to use [default %default]"
  ),
  make_option(c("-T", "--truth"),
    type = "character", default = "merged_hg38.svs.sort.oa.vcf.gz",
    help = "The VCF with the truth set [default %default]"
  ),
  make_option(c("-s", "--sample"),
    type = "character", default = "NA12878",
    help = "The sample to focus on [default %default]"
  )
)

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)

args <- parse_args(parser)
useRDS <- args$rds

## THIS IS ERRORING
here::i_am("r-scripts/05-tpbase.R")


# load vcf processing functions
source(here::here("r-scripts", "00-vcf_process.R"))

plot_width <- 16
plot_height <- 9
options(repr.plot.width = plot_width, repr.plot.height = plot_height, scipen = 999)

theme_set(theme_bw())
truvari_folder <- here::here(args$input, "truvari")
data_folder <- args$truvari

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
focus_sample <- args$sample

save_folder <- here::here(truvari_folder, data_folder, focus_sample)
# quit(save="ask")

summary_file <- here::here(truvari_folder, data_folder, "data.Rds")
color_file <- here::here(truvari_folder, "color_pal.Rds")
if (!file.exists(summary_file) | !file.exists(color_file)) { # make sure these files exist
  source(here("r-scripts", "01-truvari-report.R"))
}
# load the f1 summary data
summary_data <- readRDS(file = summary_file)
# get the color and orders of the callers
color_pal <- readRDS(color_file)
color_pal <- c(color_pal, "#000" = "TRUTH")
print(color_pal)

overall_data <- vcf_process_all_raw(
  truvari_folder = truvari_folder, data_folder = data_folder,
  sample = focus_sample,
  what = "tp-base", restrict = FALSE, debug = TRUE
)
overall_data <- overall_data %>% rename_with(~ gsub("_1", "", .), ends_with("_1"))
f1_order <- c(levels((summary_data |> mutate(caller = fct_reorder(caller, f1)))$caller), "TRUTH")
overall_data <- overall_data |>
  mutate(caller = fct_relevel(caller, f1_order)) |>
  arrange(caller) # reorder by f1
head(overall_data)
summary_all <- overall_data %>%
  group_by(caller) %>%
  dplyr::summarise(mean_SVLEN = mean(SVLEN), median_SVLEN = median(SVLEN), max_SVLEN = max(SVLEN), min_SVLEN = min(SVLEN), count = n())
summary_all_type <- overall_data %>%
  group_by(caller, SVTYPE) %>%
  dplyr::summarise(mean_SVLEN = mean(SVLEN), median_SVLEN = median(SVLEN), max_SVLEN = max(SVLEN), min_SVLEN = min(SVLEN), count = n())
summary_order <- c(levels((summary_all %>% mutate(caller = fct_reorder(caller, median_SVLEN)))$caller))


summary_del <- overall_data %>%
  filter(SVTYPE == "DEL") %>%
  group_by(caller) %>%
  dplyr::summarise(mean_SVLEN = mean(SVLEN), median_SVLEN = median(SVLEN), max_SVLEN = max(SVLEN), min_SVLEN = min(SVLEN), count = n())
del_levels <- levels((summary_del %>% mutate(caller = fct_reorder(caller, median_SVLEN)))$caller)
del_order <- c(del_levels[del_levels != "TRUTH"], "TRUTH")

standard_theme <- theme(legend.position = "none", panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())

overall_data <- overall_data %>% mutate(SVLEN_range = case_when(
  SVLEN <= 100 ~ "50-100",
  SVLEN <= 500 ~ "100-500",
  SVLEN <= 1000 ~ "500-1000",
  SVLEN > 1000 ~ "1000+",
))
binned_svlen <- overall_data %>%
  group_by(caller, SVLEN_range, SVTYPE) %>%
  summarise(nFN = sum())
