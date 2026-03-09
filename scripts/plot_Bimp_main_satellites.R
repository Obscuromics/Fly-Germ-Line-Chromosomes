# Plot results of blasting satellites from short arm of the chrX to GRCs
################################################################################
require(valr)
require(ggplot2)
require(dplyr)
################################################################################
# load chromosome sizes
load_chrom_sizes <- function(file){
  df <- read.csv(file, header = FALSE, sep = "\t")
  colnames(df) <- c("chrom", "size")
  return(df)
}

# read raw TRASH file
read_trash_raw_file <- function(file){
  satellites <- read.csv(file, header = TRUE)[,c(1,2,3,8)]
  colnames(satellites) <- c("start", "end", "rep", "chrom")
  satellites <- satellites %>% relocate(chrom, .before = start)
  return(satellites)
}
################################################################################
setwd("/Users/ab66/Documents/sanger_work/diptera/analysis_on_curated_genomes")
home <- getwd()
data <- "data"
figures <- "figures"

# main satellites in Bimp
satellite_rep_to_plot <- c(156, 222, 234)
################################################################################
# read chromosome sizes
sizes <- load_chrom_sizes(
  file = file.path(home, data, "chrom_sizes/idBraImpa2.1.primary.chrom_sizes.tsv"))

sizes <- sizes %>%
  dplyr::mutate(start = rep(0, nrow(sizes))) %>% 
  dplyr::rename("end" = "size") %>%
  dplyr::relocate(end, .after = start)

# read satellite annotation
satellites <- read_trash_raw_file(
  file = file.path(home, data, "trash_satellites/all.repeats.from.bimp.fa.csv"))
satellites <- satellites %>%
  filter(chrom %in% sizes$chrom) %>%
  filter(rep %in% satellite_rep_to_plot)

# calculate satellites in windows
windows <- bed_makewindows(sizes, win_size = 50000)
intersect <- bed_intersect(windows, satellites)
agg_intersect <- aggregate(
  intersect$chrom,
  by = list(intersect$chrom, intersect$start.x, intersect$end.x, intersect$rep.y),
  FUN = length)
colnames(agg_intersect) <- c("chrom", "start", "end", "rep", "count")
agg_intersect$rep <- as.character(agg_intersect$rep)
agg_intersect <- left_join(windows, agg_intersect) # force all windows to be present 
agg_intersect$chrom <- factor(
  agg_intersect$chrom, 
  levels = c("SUPER_1", "SUPER_2", "SUPER_3", "SUPER_X", "SUPER_GRC"))

# plot
cp_repeats <- c(
  "156" = "#8cc49f",
  "222" = "#efa714",
  "234" = "#cd8a00"
)

p <- agg_intersect %>%
  ggplot(aes(x = start/1000000, y = count, colour = rep)) +
  geom_segment(aes(x = start/1000000, xend = start/1000000, y = 0, yend = count), linewidth = 0.8) +
  geom_point(size = 1) +
  scale_color_manual(values = cp_repeats, name = "Satellite\nrepeat") +
  scale_fill_manual(values = cp_repeats, name = "Satellite\nrepeat") +
  labs(x = "Genome position (Mbp; window size = 50 kbp)", y = "Number of repeat copies") +
  facet_wrap(~chrom, nrow = 1, scales = "free_x") +
  theme_bw()

ggsave(
  file = file.path(home, figures, "Bimp_main_satellites.svg"), 
  plot = p, width = 10, height = 2.5)
