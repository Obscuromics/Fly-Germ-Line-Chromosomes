# Plot results of blasting satellites from short arm of the chrX to GRCs
################################################################################
require(valr)
require(ggplot2)
require(ggnewscale)
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

# read paf alignment file
load_paf_file <- function(file, t_chroms = NULL, t_ranges = NULL, t_sp = NULL,
                          q_chroms = NULL, q_ranges = NULL, q_sp = NULL){
  
  # read file
  ali <- read.table(file = file, header = FALSE)[,c(1:12)]
  colnames(ali) <- c("target", "t_length", "t_start", "t_end", "t_string",
                     "query", "q_length", "q_start", "q_end",
                     "n_matches", "n_bases", "score")
  
  # filter target chromosomes
  if(!is.null(t_chroms)){
    ali <- ali %>%
      filter(target %in% t_chroms)
  }
  
  # filter query chromosomes
  if(!is.null(q_chroms)){
    ali <- ali %>%
      filter(query %in% q_chroms)
  }
  
  # filter target ranges
  if(!is.null(t_ranges)){
    ali <- ali %>%
      filter((t_start > t_ranges[1]) & (t_end < t_ranges[2]))
  }
  
  # filter q_ranges
  if(!is.null(q_ranges)){
    ali <- ali %>%
      filter((q_start > q_ranges[1]) & (q_end < q_ranges[2]))
  }
  
  ali <- ali %>% arrange(target, t_start)
  ali$t_species <- rep(t_sp, nrow(ali))
  ali$q_species <- rep(q_sp, nrow(ali))
  return(ali)
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
  dplyr::relocate(end, .after = start)

# read centromeres file
centromeres <- read.csv(
  file = file.path(home, data, "centromeres/Bcop.centromeres.bed"), 
  sep = "\t", header = FALSE)
colnames(centromeres) <- c("chrom", "start", "end")

# read satellite annotation
satellites <- read_trash_raw_file(
  file = file.path(home, data, "trash_satellites/all.repeats.from.bcop.fa.csv"))
sat_chrX_short_arm <- satellites %>% 
  filter(chrom == "SUPER_X") %>%
  filter(start >= centromeres$start[centromeres$chrom == "SUPER_X"]) %>%
  filter(rep %in% c(155, 162, 73, 94, 118, 237))
sat_chrX_short_arm$rep <- as.character(sat_chrX_short_arm$rep)
sat_chrX_short_arm$rep[sat_chrX_short_arm$rep == "73"] <- "62"

# load synteny blocks with GRCs
sbs <- read.csv(
  file = file.path(home, data, "controlling_element/Bcop_GRCs2chrX.bed.tsv"),
  header = TRUE, sep = "\t")
sbs_in_chrX <- sbs %>% filter(chrom == "SUPER_X")
sbs_in_chrX$length <- as.character(sbs_in_chrX$length)

# load blast results of 185 bp sat to the core
toCore_185 <- read.csv(
  file = file.path(blast_res, "Bcop_GRC_185_bp_to_core.blastn_match.fmt6.out.tsv"),
  header = FALSE, sep = "\t")[,c(2,3,4,10,11)]
colnames(toCore_185) <- c("rep", "chrom", "identity", "start", "end")
toCore_185$rep <- as.character(toCore_185$rep)

################################################################################
chrX_sizeToPlot <- sizes[sizes$chrom == "SUPER_X",]
chrX_sizeToPlot$start <- centromeres$start[centromeres$chrom == "SUPER_X"]

cp_repeats <- c(
  "155" = "#44b5cf",
  "162" = "#5badb9",
  "185" = "#f3dc77",
  "118" = "#ffcc6a",
  "237" = "#ffa28c",
  "62"  = "#de86e4",
  "94"  = "#ff88bb"
  )

cp_sbs <- c(
  "10072" = "grey1",
  "355"   = "grey30",
  "8422"  = "grey50"
  )

p <- ggplot() +
  geom_rect(data = chrX_sizeToPlot, 
            aes(xmin = start/1000, xmax = end/1000, ymin = 0, ymax = 1), 
            colour = "darkgrey", fill = "white", linewidth = 0.1) +
  geom_rect(data = sat_chrX_short_arm,
            aes(xmin = start/1000, xmax = end/1000, ymin = 0, ymax = 1,
                colour = rep, fill = rep)) +
  geom_rect(data = toCore_185,
            aes(xmin = start/1000, xmax = end/1000, ymin = 0, ymax = 1,
                colour = rep, fill = rep)) +
  scale_color_manual(values = cp_repeats, name = "repeats") +
  scale_fill_manual(values = cp_repeats, name = "repeats") +
  ggnewscale::new_scale_fill() +
  ggnewscale::new_scale_color() +
  geom_rect(data = sbs_in_chrX, 
            aes(xmin = start/1000, xmax = end/1000, ymin = 0, ymax = 1,
                colour = length, fill = length)) +
  scale_color_manual(values = cp_sbs, name = "synteny block") +
  scale_fill_manual(values = cp_sbs, name = "synteny block") +
  geom_rect(data = centromeres %>% filter(chrom == "SUPER_X"), 
            aes(xmin = start/1000, xmax = end/1000, ymin = 0, ymax = 1), 
            colour = "red", fill = NA, alpha = 0.1) +
  theme_bw()

p <- p + theme(legend.position = "none")
  
ggsave(
  file = file.path(home, figures, "Bcop_SUPER_X_short_arm.svg"), 
  plot = p, width = 10, height = 2)
################################################################################
windows <- bed_makewindows(sizes, win_size = 1000)

sat_chrX_short_arm <- rbind(
  sat_chrX_short_arm, 
  data.frame("chrom" = toCore_185$chrom, "start" = toCore_185$start,
             "end" = toCore_185$end, "rep" = toCore_185$rep))

intersect <- bed_intersect(windows, sat_chrX_short_arm)
agg_intersect <- aggregate(intersect$chrom,
                           by = list(intersect$chrom, intersect$start.x, intersect$end.x, intersect$rep.y),
                           FUN = length)
colnames(agg_intersect) <- c("chrom", "start", "end", "rep", "count")
agg_intersect$rep <- as.character(agg_intersect$rep)
agg_intersect <- agg_intersect %>% filter(start > 68000000)

p2 <- agg_intersect %>%
  ggplot(aes(x = start/1000000, y = count, colour = rep)) +
  annotate("rect", xmin = 68424739/1000000, xmax = 68433161/1000000,
           ymin = 0, ymax = Inf, fill = "grey", alpha = 0.4) +
  geom_segment(aes(x = start/1000000, xend = start/1000000, y = 0, yend = count), linewidth = 0.8) +
  geom_point(size = 1) +
  scale_color_manual(values = cp_repeats, name = "repeats") +
  scale_fill_manual(values = cp_repeats, name = "repeats") +
  labs(x = "SUPER_X (Mbp; window size = 1 kbp)", y = "Number of repeat copies") +
  
  theme_bw()

ggsave(
  file = file.path(home, figures, "Bcop_SUPER_X_short_arm.zoomed.svg"), 
  plot = p2, width = 10, height = 3.5)
################################################################################
# synteny blocks between the short arm of the chrX and both GRCs
alignment_file <- file.path(home, data, "alignments/Bcop_Bcop.1aln.paf")

# synteny blocks to plot
sbs <- list(
  
  "X_GRC1" = list("target" = "SUPER_X", 
                   "query" = "SUPER_GRC1", 
                   t_coord = c(66331842-1000, 66345636+1000), 
                   q_coord = c(2260413-1000, 2275104+1000)),
  "X_GRC2" = list("target" = "SUPER_X", 
                   "query" = "SUPER_GRC2", 
                   t_coord = c(66331837-1000, 66345636+1000),
                   q_coord = c(43356409-1000, 43366882+1000))
)

for(i in names(sbs)){
  sb <- sbs[[i]]
  ali <- load_paf_file(file = alignment_file,
                       t_chroms = sb$target, q_chroms = sb$query,
                       t_ranges = sb$t_coord, q_ranges = sb$q_coord)
  
  ### --- plotting --- ###
  p <- ggplot(data.frame(
    start = ali$t_start,
    end = ali$t_end,
    strand = ali$t_string,
    query.start = ifelse(ali$t_string == "+", ali$q_start, ali$q_end),
    query.end = ifelse(ali$t_string == "+", ali$q_end, ali$q_start)
  )) +
    
    aes(x = start/1000, y = query.start/1000, xend = end/1000, yend = query.end/1000)
  
  p <- p + geom_segment(lineend = "round", linewidth = 0.4) +
    labs(x = paste0(unique(ali$target), " (Kbp)"), 
         y = paste0(unique(ali$query), " (Kbp)")) +
    theme_bw() +
    coord_fixed()
  
  plt_file_name <- paste0("Bcop", i, ".synteny_block")
  
  ggsave(plot = p , filename = file.path(home, figures, paste0(plt_file_name, ".svg")), 
         device = "svg", width = 4.5, height = 3.5)
}