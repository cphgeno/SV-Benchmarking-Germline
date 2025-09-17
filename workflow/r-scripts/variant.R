library(tidyverse)
library(StructuralVariantAnnotation)
library(ggforce)
library(here)
library(janitor)

here::i_am("r-scripts/variant.R")

truvari <- "type-ignored"

regionChr <- "chrX"
regionStart <- 11022525
regionEnd <- 11025465
regionName <- str_c(regionChr, ":", regionStart, "-", regionEnd)
param <- ScanVcfParam(which = GRanges(regionChr, IRanges(regionStart, regionEnd)))

group.colors <- c(
  "True Positive" = "#00BFC4",
  "False Positive" = "#F8766D",
  "False Negative" = "grey"
)

overall_data <- vcf_process_all(
  data_folder = truvari, param = param,
  what = c("tp-comp", "fp", "fn")
)
overall_data <- overall_data %>%
  # mutate(TruScore=ifelse(is.na(TruScore), 0, TruScore)) %>%
  filter(start >= regionStart, end <= regionEnd) |>
  remove_constant() |>
  remove_empty()

ggplot(overall_data, aes(x = caller, fill = status)) +
  geom_bar(na.rm = FALSE) +
  scale_x_discrete(drop = FALSE) +
  labs(
    title = str_c("Variants in ", regionName),
    subtitle = str_c("Range of ", round((regionEnd - regionStart) / 1000, digits = 2), "kb")
  ) +
  scale_fill_manual(values = group.colors)
ggplot2::ggsave(here("results", "truvari", truvari, "igv", str_c(regionChr, "_", regionStart, "-", regionEnd, "_tp-fp.png")))

ggplot(overall_data, aes(x = caller, fill = type)) +
  facet_grid(sample ~ status, scales = "free_x") +
  geom_bar(na.rm = FALSE) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(
    title = str_c("Variants in ", regionName),
    subtitle = str_c("Range of ", round((regionEnd - regionStart) / 1000, digits = 2), "kb")
  )
ggplot2::ggsave(here("results", "truvari", truvari, "igv", str_c(regionChr, "_", regionStart, "-", regionEnd, "_type.png")))

View(overall_data |>
  mutate(caller = case_when(status == "False Negative" ~ "TRUTH", TRUE ~ caller)) |>
  remove_constant() |>
  remove_empty() |>
  unique() |>
  dplyr::select(-starts_with("Cluster")) |>
  dplyr::select(-starts_with("Dist")) |>
  dplyr::select(-starts_with("HG")) |>
  dplyr::select(-starts_with("Num")) |>
  dplyr::select(-starts_with("Multi")) |>
  dplyr::select(-ends_with("calls")) |>
  dplyr::select(-starts_with("nator")) |>
  # dplyr::select(-(ExactMatchIDs)) |>
  # dplyr::select(-REFWIDENED) |>
  arrange(sample, status, caller))

ggplot(overall_data)
