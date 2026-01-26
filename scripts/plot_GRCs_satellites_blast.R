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
################################################################################
setwd("/Users/ab66/Documents/sanger_work/diptera/analysis_on_curated_genomes")
home <- getwd()
data <- "data"
figures <- "figures"
blast_res <- file.path(home, "/satellites/filtered_blast")
################################################################################
# read chromosome sizes
sizes <- load_chrom_sizes(
  file = file.path(home, data, "chrom_sizes/idBraCopr2.1.chrom_sizes.tsv"))
sizes <- sizes %>%
  dplyr::mutate(start = rep(0, nrow(sizes))) %>% 
  dplyr::rename("end" = "size") %>%
  dplyr::relocate(end, .after = start) %>%
  dplyr::filter(chrom %in% c("SUPER_GRC1", "SUPER_GRC2"))

# load blast results
toGRCs_185 <- read.csv(
  file = file.path(blast_res, "Bcop_GRC_185_bp_to_GRCs.blastn_match.fmt6.out.tsv"),
  header = FALSE, sep = "\t")[,c(2,3,4,10,11)]
colnames(toGRCs_185) <- c("rep", "chrom", "identity", "start", "end")
toGRCs_185 <- bed_merge(toGRCs_185)
toGRCs_185$rep <- rep("185", nrow(toGRCs_185))

toGRCs_others <- read.csv(
  file = file.path(blast_res, "Bcop_SUPER_X_sat_consensus_to_GRCs.blastn_match.fmt6.out.tsv"),
  header = FALSE, sep = "\t")[,c(2,3,4,10,11)]
colnames(toGRCs_others) <- c("rep", "chrom", "identity", "start", "end")
toGRCs_others <- toGRCs_others %>%
  group_by(rep) %>% bed_merge()

blast_res <- as.data.frame(rbind(toGRCs_185, toGRCs_others))
################################################################################
windows <- bed_makewindows(sizes, win_size = 50000)
intersect <- bed_intersect(windows, blast_res)
agg_intersect <- aggregate(intersect$chrom,
                           by = list(intersect$chrom, intersect$start.x, intersect$end.x, intersect$rep.y),
                           FUN = length)
colnames(agg_intersect) <- c("chrom", "start", "end", "rep", "count")
agg_intersect$rep <- as.character(agg_intersect$rep)
agg_intersect <- left_join(windows, agg_intersect)

cp_repeats <- c(
  "62"  = "#5f8accff",
  "118" = "#a775b7ff",
  "185" = "#f3dc77ff",
  "237" = "#c0618cff",
  "94"  = "#8cc49fff"
)

p <- agg_intersect %>%
  ggplot(aes(x = start/1000000, y = count, colour = rep)) +
  annotate("rect", xmin = 68424739/1000000, xmax = 68433161/1000000,
           ymin = 0, ymax = Inf, fill = "grey", alpha = 0.4) +
  geom_segment(aes(x = start/1000000, xend = start/1000000, y = 0, yend = count), linewidth = 0.8) +
  geom_point(size = 1) +
  scale_color_manual(values = cp_repeats, name = "Satellite\nrepeat") +
  scale_fill_manual(values = cp_repeats, name = "Satellite\nrepeat") +
  labs(x = "Genome position (Mbp; window size = 50 kbp)", y = "Number of repeat copies") +
  facet_wrap(~chrom, nrow = 1, scales = "free_x") +
  theme_bw()

ggsave(
  file = file.path(home, figures, "Bcop_satellite_blast_res.svg"), 
  plot = p)#, width = 10, height = 2.5)
