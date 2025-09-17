library(here)
library(tidyverse)
library(tidyr)
library(purrr)
library(lubridate)
library(ggpubr)
library(patchwork)
library(ggforce)

theme_set(theme_bw(base_size = 9) + theme(
					   axis.title.x=element_text(size=9),
					   axis.title.y=element_text(size=9),
					   axis.text.x = element_text(size=6)))

options(scipen = 999)
# read in all the data
here::i_am("workflow/r-scripts/03-benchmarks.R")
data_folder <- here("pedigree-results", "benchmarks")
files <- list.files(path = data_folder, full.names = TRUE, pattern = "*.tsv", recursive = TRUE)
print(files)
data <- files |>
  map(read_delim, col_types = "d_ddddddd", delim = "\t") |>
  reduce(bind_rows)
data$file <- str_remove(str_remove(files, data_folder), ".tsv")

# get the callers and samples
data <- data |>
  mutate(file = str_remove(file, data_folder)) |>
  separate(file, c(NA, "caller", "sample"), sep = "/")
data <- data |> mutate(sample = as.factor(sample))

# some callers are spread across multiple benchmarks
data <- data |> separate(caller, c("caller", "subprocess"), sep = "-")

print(data)
# sum up the subprocesses
data <- data |>
  group_by(caller, sample) |>
  summarise(s = sum(s), cpu_time = sum(cpu_time), max_rss = max(max_rss), max_vms = max(max_vms), max_uss = max(max_uss), max_pss = max(max_pss), io_in = sum(io_in), io_out = sum(io_out), mean_load = mean(mean_load))

# pretty format seconds
data <- data |> mutate(time = lubridate::seconds_to_period(s), cpu_time_pretty = lubridate::seconds_to_period(cpu_time))

write.csv(data, file = here(data_folder, "caller-summary.csv"))

print(data)
print(levels(data$sample))




# grab the colors
color_pal <- readRDS(here("pedigree-results", "truvari", "color_pal.Rds"))
# and subset only to valid
print(unique(data$caller))
print(names(color_pal))
data <- data %>%
  filter(caller %in% names(color_pal)) %>%
  mutate(caller = factor(caller, levels = names(color_pal)))
# and plot!
sec_to_hours <- function(x) {
  sprintf("%.1f", x / 3600)
}

sec_breaks <- function(lims) {
  range_sec <- lims[2] - lims[1]
  if (range_sec <= 12 * 3600) {
    step <- 3600     # 1 hour 
  } else if (range_sec <= 48 * 3600) {
    step <- 4 * 3600   # 4 hours
  } else {
    step <- 8 * 3600  # 8 hours
  }
  seq(0, ceiling(lims[2] / step) * step, by = step)
}
median_sdl <- function(x, mult = 1) {
  m <- median(x, na.rm = TRUE)
  sd <- sd(x, na.rm = TRUE)
  data.frame(y = m, ymin = m - mult * sd, ymax = m + mult * sd)
}
benchmark_plot <- list(
  stat_summary(fun = "median", geom = "col"), 
  stat_summary(fun.data = median_sdl, geom="linerange", linewidth=0.8,alpha=0.5), 
  scale_fill_manual(values = color_pal, guide = "none"),
  #geom_point(aes(shape = sample),position=position_jitter(width=0.05,height=0)),
  scale_y_continuous(expand=c(0,0)),
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)))

# Running time (s)
runtime <- ggplot(data %>% filter(s > 0), aes(fill = caller, x = caller, y = s)) +
    benchmark_plot + 
  scale_y_continuous(expand=c(0,0),  breaks=sec_breaks, labels = function(x) sprintf("%.1f", x / 3600)) + 
    labs(x = "Callers", y = "Running Time (h)") 
ggplot2::ggsave(runtime,width = 16, height = 9, file = here(data_folder, "max_sec.png"))
# CPU Time (s)
cputime <- ggplot(data %>% filter(cpu_time > 30), aes(fill = caller, x = caller, y = cpu_time)) +
    benchmark_plot + 
  scale_y_continuous(expand=c(0,0), breaks=sec_breaks, labels = function(x) as.integer(round(x / 3600))) + 
  labs(x = "Callers", y = "CPU Time (h)") 
ggplot2::ggsave(cputime, width = 16, height = 9, file = here(data_folder, "cpu_time.png"))

run_to_cpu <- ggplot(data, aes(color = caller, x = s, y = cpu_time, shape = sample)) +
  geom_point() +
  scale_color_manual(values = color_pal) +
  labs(x = " Run time (s)", y = "CPU time (s)")
ggplot2::ggsave(run_to_cpu, width = 16, height = 9, file = here(data_folder, "run_to_cpu.png"))

cpu_to_mem <- ggplot(data, aes(color = caller, x = cpu_time, y = max_rss / 1024, shape = sample)) +
  geom_point() +
  scale_color_manual(values = color_pal) +
  labs(x = "CPU time (s)", y = "RAM allocated (Gb)")
ggplot2::ggsave(cpu_to_mem, width = 16, height = 9, file = here(data_folder, "cpu_to_rss.png"))

# RSS - RAM allocated during allocation, including preloaded libraries
rss <- ggplot(data %>% filter(max_rss > 0), aes(fill = caller, x = caller, y = max_rss / 1024)) +
  labs(x = "Callers", y = "Max RSS Memory Usage (Gb)") +
    benchmark_plot  
ggplot2::ggsave(rss,width = 16, height = 9, file = here(data_folder, "max_rss.png"))

vms <- ggplot(data %>% filter(max_vms > 0), aes(fill = caller, x = caller, y = max_vms / 1024)) +
  labs(x = "Callers", y = "Max Virtual Memory Usage (Gb)") +
    benchmark_plot  
ggplot2::ggsave(vms, width = 16, height = 9, file = here(data_folder, "max_vms.png"))

# Unique Set Size libraries and pages allocated only to this process
uss <- ggplot(data, aes(fill = caller, x = caller, max_uss)) +
  labs(x = "Callers", y = "USS (Mb)") +
    benchmark_plot  
ggplot2::ggsave(uss, width = 16, height = 9, file = here(data_folder, "max_uss.png"))

# Mean load
load <- ggplot(data, aes(fill = caller, x = caller, y = mean_load)) +
  labs(x = "Callers", y = "Mean Load") +
    benchmark_plot  
ggplot2::ggsave(load, width = 16, height = 9, file = here(data_folder, "mead_load.png"))

# IO in
ioin <- ggplot(data, aes(fill = caller, x = caller, y = io_in)) +
  labs(x = "Callers", y = "Max IO In") +
    benchmark_plot  
ggplot2::ggsave(ioin, width = 16, height = 9, file = here(data_folder, "max_ioin.png"))

summary <- runtime + cputime + rss + vms + plot_annotation() + plot_layout(guides='collect',axes='collect')
ggplot2::ggsave(summary, width = 16, height = 9, file = here(data_folder, "summary.png"))
    ggplot2::ggsave(summary, width = 170, height = 100, units="mm", dpi=300, file = here(data_folder, "summary.pdf"))

padded_num <- function(x, pad_width = 3) {
  formatted <- sprintf("%d", as.integer(x)) # integer only
  stringr::str_pad(formatted, width = pad_width, side = "left")
}

cputime_zoom <- ggplot(data %>% filter(cpu_time > 30), aes(fill = caller, x = caller, y = cpu_time)) +
    benchmark_plot +
    scale_y_continuous(
        expand = c(0,0),
	breaks = sec_breaks, 
	labels = function(x) padded_num(x / 3600)
        ) +
    facet_zoom(y = caller == 'manta') +
    labs(x = "Callers", y = "CPU Time (h)") + 
    theme(panel.spacing.x = unit(0.1, "lines"), 
	  axis.text.x = element_blank(), 
	  axis.ticks.x = element_blank(), 
	  axis.ticks.x.top = element_blank(), 
	  axis.text.x.top = element_blank())

cputime_manta <- max((data %>% filter(caller == 'manta'))$cpu_time)
# Base plot without default x-axis text
vms_zoom <- vms +
  facet_zoom(ylim = c(0, 80)) +
    theme(panel.spacing.x = unit(0.1, "lines")) 


summary_zoom <- (runtime + coord_cartesian(ylim=c(0,cputime_manta)))  + cputime_zoom + (rss + coord_cartesian(ylim=c(0,80))) + vms_zoom + plot_annotation() + plot_layout(guides='collect',axes='collect_x')
ggplot2::ggsave(summary_zoom, width = 16, height = 9, file = here(data_folder, "summary_zoom.png"))
ggplot2::ggsave(summary_zoom, width = 170, height = 100, units="mm", dpi=300, file = here(data_folder, "summary_zoom.pdf"))
