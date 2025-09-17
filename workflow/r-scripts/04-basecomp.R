library(tidyverse) 
# need openssl 1.1.1, make sure to load it
library(StructuralVariantAnnotation) 
library(janitor)

options(repr.plot.width = 16, repr.plot.height = 9)

theme_set(theme_bw())

data_folder <-
  "typeignore/"
folders <-
  list.dirs(
    path = data_folder,
    full.names = FALSE,
    recursive = TRUE
  )
print(folders)
callers <- unique(str_split_i(str_remove(folders, data_folder), "/", 1))[-1]
samples <- unique(str_split_i(str_remove(folders, data_folder), "/", 2))[-1]
print(samples)
print(callers)
what <- c("tp-comp", "tp-base")
sample <- samples[1]
caller <- "manta"
# quit(save="ask")

summary_file <- str_c(data_folder, "data.Rds")
if (!file.exists(summary_file)) {
  source("../../r-scripts/01-truvari-report.R")
}
data <- readRDS(file = summary_file)

vcf_read <- function(sample, caller) {
  vcfs <- NULL
  vcf_files <-
    list.files(
      path = data_folder,
      recursive = TRUE,
      pattern = str_c("*.vcf.gz$"),
      full.names = TRUE
    ) # get a list of all the summary files
  vcf_files <- lapply(what, function(x) str_c(data_folder, "/", caller, "/", sample, "/", x, ".vcf.gz"))
  print(paste("vcf files for", sample, caller))
  print(vcf_files)
  vcfs <- lapply(vcf_files, readVcf,
    genome = "hg38",
    ScanVcfParam(fixed = NA, geno = NA, info = c("SVTYPE", "MatchId", "TruScore"))
  )

  names(vcfs) <- what

  for (vcf in vcfs) {
    seqlevels(vcf) <-
      str_remove(seqlevels(vcf), "chr") # change the chr names for plotting
    # try(end(vcf) <- vcf@info$END)
    # properly add the ends
  }
  return(vcfs)
}

# caller <- "manta"
overall_data <- data.frame()
for (sample in samples) {
  for (caller in callers) {
    tryCatch({
      print(caller)
      vcfs <- vcf_read(sample, caller)
      print("VCFs loaded")

      for (i in seq(1, length(vcfs))) {
        print(str_c("Processing ", names(vcfs)[i]))
        vcf <- vcfs[[i]]
        test <- as.data.frame(rowRanges(vcf), optional = TRUE, row.names = NULL)
        type <- as.data.frame(vcf@info) %>% dplyr::select(!where(is.list))
        data <- bind_cols(test, type)
        data <- data %>%
          rownames_to_column() %>%
          mutate(sample = sample) %>%
          mutate(caller = caller) %>%
          mutate(status = names(vcfs)[i])
        overall_data <- bind_rows(overall_data, data, )
      }
    })
  }
}

# do cleanup
rm(data)
rm(vcf)
rm(vcfs)

# save that hard won vcf data :)
# saveRDS(overall_data, file = str_c(data_folder, "/vcf_data.Rds"))

group.colors <- c(
  "tp-comp" = "blue",
  "fp" = "red",
  "fn" = "grey"
)

overall_data <- rename(overall_data, c("SVTYPE" = "type"))
# overall_data <- overall_data |>  filter(!is.na(type))

overall_data <- overall_data %>%
  mutate(sample = as_factor(sample)) %>%
  mutate(caller = as_factor(str_to_lower(caller)))
print(levels(overall_data$sample))
sample_names <- c("200915", "200921", "NA24385")
sample <- sample_names[1]
levels(overall_data$sample) <- sample_names
overall_data <- overall_data |> filter(sample == sample)

overall_data <- overall_data |> mutate(type = as_factor(type))

matches <- overall_data |>
  dplyr::select(type, MatchId, sample, caller, status, TruScore) |>
  pivot_wider(names_from = status, values_from = type)

ggplot(matches, aes(x = `tp-comp`, fill = `tp-base`)) +
  geom_bar() +
  facet_wrap(~caller, scales = "free") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ggrepel::geom_label_repel(
    data = matches |> filter(`tp-comp` != `tp-base` & !(`tp-comp` == "DUP" & `tp-base` == "INS")),
    mapping = aes(x = `tp-comp`, fill = `tp-base`, label = ..count..), stat = "count",
    direction = "y"
  ) +
  labs(
    title = "Type Comparison in True Positive variants called with type ignore in sample 200915",
    subtitle = "Count of matches other than INS to INS or DUP and DEL to DEL shown",
    fill = "Type in Truth Set", x = "Type in Caller"
  )

matches <- matches |>
  mutate(is_match = case_when(
    `tp-comp` == `tp-base` ~ "FULL",
    `tp-comp` == "DUP" & `tp-base` == "INS" ~ "DUPTOINS",
    (caller == "tardis" | (caller == "tiddit" & `tp-comp` != "DEL")) & `tp-base` == "INS" ~ "PARTIAL",
    `tp-comp` == "INS" | `tp-comp` == "DEL" ~ "MISMATCH", TRUE ~ "OTHER"
  )) |>
  mutate(is_match = factor(is_match, levels = c("FULL", "DUPTOINS", "PARTIAL", "OTHER", "MISMATCH")))
match_counts <- matches |>
  group_by(is_match) |>
  summarise(n = n(), Med_TruScore = median(TruScore))

ggplot(matches, aes(x = is_match, y = TruScore)) +
  geom_jitter(aes(color = caller), size = 0.5, alpha = 0.8) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.8, adjust = 1.5) +
  labs(title = "True Positive Matches and their TruVari score, based on type similarity")

ggplot(matches, aes(x = is_match, y = TruScore)) +
  geom_jitter(size = 0.3, alpha = 0.5) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.8, adjust = 1.5) +
  facet_wrap(~caller, scale = "free_x") +
  labs(title = "True Positive Matches and their TruVari score, based on type similarity")
