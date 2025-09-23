# Plot GC content and satellite distribution along chromosomes with annotated
# centromeres
################################################################################
require(valr)
require(ggplot2)
require(dplyr)
require(scales)
require(ggpubr)
require(pheatmap)
################################################################################
# Functions
# read GC content distribution file
read_gc_content <- function(file){
  gc <- read.csv(file, header = TRUE, sep = "\t")
  gc <- gc[,c(1:3,5)]
  colnames(gc) <- c("chrom", "start", "end", "frac")
  gc["feature"] <- rep("GC", nrow(gc))
  return(gc)
}

# read raw TRASH file
read_raw_trash_file <- function(file){
  satellites <- read.csv(file, header = TRUE)
  satellites <- satellites %>% select(seq.name, start, end, width, seq)
  colnames(satellites) <- c("chrom", "start", "end", "width", "seq")
  return(satellites)
}

# read summary TRASH file
read_sum_trash_file <- function(file){
  trash_summary <- read.csv(file, header = TRUE)[, c(2, 3, 4, 6, 7, 8)]
  colnames(trash_summary) <- c("chrom", "start", "end", "rep", "consensus", "count")
  trash_summary$rep <- as.character(trash_summary$rep)
  return(trash_summary)
}

# load chromosome sizes
load_chrom_sizes <- function(file){
  df <- read.csv(file, header = FALSE, sep = "\t")
  colnames(df) <- c("chrom", "size")
  return(df)
}

# read file with candidate centromeric regions
read_centromeres <- function(file, chrom_sizes){
  centromeres <- read.csv(file, header = TRUE, sep = "\t")[,c(1:3)]
  centromeres$startToPlot <- NA
  centromeres$endToPlot <- NA
  
  for(i in 1:nrow(centromeres)){
    df <- centromeres[i,]
    chrom_length <- sizes[sizes$chrom == df$chrom,]$size
    if(df$start - 10000 > 0){
      centromeres[i,]$startToPlot <- df$start - 10000
    }else{
      centromeres[i,]$startToPlot <- 0
    }
    
    if(df$end + 10000 < chrom_length){
      centromeres[i,]$endToPlot <- df$end + 10000
    }else{
      centromeres[i,]$endToPlot <- chrom_length
    }
  }
  return(centromeres)
}

# Generate list of kmers for a satellite repeat
generate_kmers <- function(sequence, k = 3){
  count = 1
  kmers <- list()
  while((count + k - 1) <= nchar(sequence)){
    kmer <- substr(sequence, count, count + k - 1)
    kmers <- append(kmers, kmer)
    count <- count + 1
  }
  return(unlist(kmers))
}

# Generate jaccard similarity index for a pair of sequences
jaccard_similarity <- function(kmers_a, kmers_b){
  intersection <- length(intersect(kmers_a, kmers_b))
  union <- length(union(kmers_a, kmers_b))
  similarity = intersection / union
  return(similarity)
}

# generate reverse function
rev_sequence <- function(sequence) {
  alphabets <- strsplit(as.character(sequence), split = "")[[1]]
  return(rev(alphabets))
}

# generate reverse and complementary sequence
compl_sequence <- function(sequence) {
  sequence <- rev_sequence(sequence)
  cmplvec <- sapply(sequence, function(base) switch(base, "A" = "T", "C" = "G", "G" = "C", "T" = "A"))
  return(paste(cmplvec, collapse = ""))
}
################################################################################
setwd("/Users/ab66/Documents/sanger_work/diptera/analysis_on_curated_genomes")
home <- getwd()
data <- "data"
figures <- "figures"

# current species to plot
species <- 
  "Bcop"
  #"Bimp"
  #"Ling"
################################################################################
gc_files <- c(
  "Bcop" = "gc/idBraCopr2.1.primary.w50kb.GC.SUPERonly.bed",
  "Bimp" = "gc/idBraImpa2.1.primary.w50kb.GC.SUPERonly.bed",
  "Ling" = "gc/idLycInge5.1.primary.w50kb.GC.SUPERonly.bed")

trash_raw_files <- c(
  "Bcop" = "trash_satellites/all.repeats.from.bcop.fa.csv",
  "Bimp" = "trash_satellites/all.repeats.from.bimp.fa.csv",
  "Ling" = "trash_satellites/all.repeats.from.lyco.fa.csv")

chrom_sizes_files <- c(
  "Bcop" = "chrom_sizes/idBraCopr2.1.chrom_sizes.tsv",
  "Bimp" = "chrom_sizes/idBraImpa2.1.primary.chrom_sizes.tsv",
  "Ling" = "chrom_sizes/idLycInge5.1.primary.chrom_sizes.tsv")

centromeres_files <- c(
  "Bcop" = "centromeres/Bcop.candidate_centromeres.tsv",
  "Bimp" = "centromeres/Bimp.candidate_centromeres.tsv",
  "Ling" = "centromeres/Ling.candidate_centromeres.tsv")

trash_summary_files <- c(
  "Bcop" = "trash_satellites/Summary.of.repetitive.regions.bcop.fa.csv",
  "Bimp" = "trash_satellites/Summary.of.repetitive.regions.bimp.fa.csv",
  "Ling" = "trash_satellites/Summary.of.repetitive.regions.lyco.fa.csv")

# repeats to plot
repeats_cp <- list(
  "Bcop" = c("46"  = "#6565ab", 
             "176" = "#a775b7", 
             "56"  = "#b34074", 
             "156" = "#db7ac8", 
             "162" = "#5badb9", # core 
             "155" = "#44b5cf", # core 
             "189" = "#5f8acc", # core 
             "185" = "#f3dc77"  # GRCs
             ),
  
  "Bimp" = c("11"  = "#cc4968",
             "14"  = "#a356a3",
             "22"  = "#0079a7", 
             "122" = "#e5e578", 
             "156" = "#8cc49f"
             ),
  
  "Ling" = c("39"  = "#95430b", # core
             "40"  = "#d36c54", # core
             "148" = "#3897c7", # core
             "296" = "#8882c3",
             "149" = "#3f669d", # core
             "126" = "#f9e195", # GRC2
             "178" = "#1c6759",
             "94"  = "#8bbdd5", # GRC1
             "109" = "#9bc589", # GRC1
             "111" = "#bfdb98", # GRC1
             "142" = "#eb985d",
             "228" = "#af5d86"  # GRC1
             )
)
################################################################################
# Load and process files
### --- read GC distribution file --- ###
gc <- read_gc_content(file = file.path(home, data, gc_files[[species]]))
windows <- gc[,c("chrom", "start", "end")]

### --- read raw TRASH file --- ###
satellites <- read_raw_trash_file(
  file = file.path(home, data, trash_raw_files[[species]]))

# overlap raw TRASH file with the windows to calculate tandem repeat coverage
intersect <- bed_intersect(windows, bed_merge(satellites))
satellites_in_windows <- aggregate(
  intersect$.overlap, by=list(intersect$chrom, intersect$start.x, intersect$end.x), FUN=sum)
colnames(satellites_in_windows) <- c("chrom", "start", "end", "cov")
satellites_in_windows["frac"] <- satellites_in_windows$cov / 50000
satellites_in_windows["feature"] <- rep("TR", length(satellites_in_windows$chrom))
satellites_in_windows <- satellites_in_windows %>%
  select("chrom", "start", "end", "frac", "feature")

### --- load chromosome sizes --- ###
sizes <- load_chrom_sizes(file = file.path(home, data, chrom_sizes_files[[species]]))

### --- load centromere coordinates --- ###
centromeres <- read_centromeres(
  file = file.path(home, data, centromeres_files[[species]]),
  chrom_sizes = sizes)
centromeres$annot <- paste(centromeres$chrom, c(1:nrow(centromeres)), sep = ".")
################################################################################
# Plotting
features_all <- rbind(gc, satellites_in_windows)
features_all$chrom <- factor(
  features_all$chrom, levels = c("SUPER_1", "SUPER_2", "SUPER_3", "SUPER_X",
                                 "SUPER_GRC1", "SUPER_GRC2", "SUPER_GRC"))
centromeres$chrom <- factor(
  centromeres$chrom, levels = c("SUPER_1", "SUPER_2", "SUPER_3", "SUPER_X",
                                "SUPER_GRC1", "SUPER_GRC2", "SUPER_GRC"))

### --- plot GC-content distribution --- ### 
gc_plt <- features_all %>%
  filter(feature == "GC") %>%
  ggplot(aes(x = start/1000000, y = frac)) +
  geom_rect(
    data = centromeres, 
    mapping = aes(xmin = start/1000000, xmax = end/1000000,ymin = -Inf, ymax = Inf, 
                  x = NULL, y = NULL), fill = "red", alpha = 0.6) +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) + 
  geom_line(size = 0.3, col = "black") +
  labs(x = "Genome position (Mbp)", y = "GC (%)") +
  facet_wrap(.~chrom, scales = "free_x", nrow = 1) +
  theme_bw()

### --- plot satellite distribution --- ###
tr_plt <- features_all %>%
  filter(feature == "TR") %>%
  ggplot(aes(x = start/1000000, y = frac)) +
  geom_rect(
    data = centromeres, 
    mapping = aes(xmin = start/1000000, xmax = end/1000000,ymin = -Inf, ymax = Inf, 
                  x = NULL, y = NULL), fill = "red", alpha = 0.6) +
  geom_point(colour = "black", size = 0.3) +
  geom_segment(aes(x = start/1000000, xend = start/1000000, y = 0, yend = frac), 
               col = "black", size = 0.3) +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) + 
  labs(x = "Genome position (Mbp)", y = "Tandem Repeats (%)") +
  facet_wrap(.~chrom, scales = "free_x", nrow = 1) +
  theme_bw()

# merge together
plt <- ggarrange(gc_plt + theme(axis.title.x = element_blank()), 
                 tr_plt, nrow = 2, common.legend = TRUE, legend = "none")

ggsave(
  plot = plt, 
  filename = file.path(home, figures, paste0(species, "_chrAll_GC_vs_satellites.svg")),
  device = "svg", units = "cm", width = 30, height = 10)
################################################################################
# Plot dominant satellites in the centromeres
centromeres$annot <- factor(centromeres$annot, levels = centromeres$annot)

# read TRASH summary file
trash_summary <- read_sum_trash_file(
  file = file.path(home, data, trash_summary_files[[species]]))

# find overlap between TRASH satellites and centromeric regions
centromere_repeats <- bed_intersect(
  centromeres[,c("chrom", "start", "end", "annot")], trash_summary)
colnames(centromere_repeats) <- c(
  "chrom", "chr_start", "chr_end", "annot", 
  "start", "end", "rep", "consensus", "count", "overlap")

# plotting
sat_plt <- ggplot() +
  geom_rect(data = centromeres, 
            aes(xmin = start/1000000, xmax = end/1000000, ymin = 0, ymax = 1), 
            colour = "darkgrey", fill = "white", size = 0.1) +
  geom_rect(data = centromere_repeats, 
            aes(xmin = start/1000000, xmax = end/1000000, ymin = 0, ymax = 1), 
            fill = "darkgrey") +
  geom_rect(data = centromere_repeats %>% filter(rep %in% names(repeats_cp[[species]])), 
            aes(xmin = start/1000000, xmax = end/1000000, ymin = 0, ymax = 1, 
                fill = factor(rep, levels = names(repeats_cp[[species]])))) +
  scale_fill_manual(values = repeats_cp[[species]], name = "Repeat") +
  labs(x = "Genome position (Mbp)", y = "satellite") +
  facet_wrap(.~annot, scales = "free_x", nrow = 1) +
  theme_bw()

ggsave(
  plot = sat_plt, 
  filename = file.path(home, figures, paste0(species, "_chrAll_centr_satellites.svg")),
  device = "svg", units = "cm", width = 30, height = 10)
################################################################################
# Calculate and plot kmer-based jaccard similarity matrix

# extract all consensus sequences for satellite repeats of a certain length 
df <- centromere_repeats[centromere_repeats$rep %in% names(repeats_cp[[species]]),]
df <- df %>% filter(!consensus == "none_identified")
consensuses_list <- df$consensus
names(consensuses_list) <- paste(
  species, df$annot, df$start, df$end, paste0("len", df$rep), sep = "_")

# generate kmers
kmers <- sapply(consensuses_list, generate_kmers, k = 6)
kmers_compl <- sapply(sapply(consensuses_list, compl_sequence), generate_kmers, k = 6)

# calculate jaccard similarity matrix
m_size <- length(names(kmers))
jaccard_m <- matrix(nrow = m_size, ncol = m_size)
colnames(jaccard_m) <- names(kmers)
rownames(jaccard_m) <- names(kmers)

for(seq_a in colnames(jaccard_m)){
  for(seq_b in rownames(jaccard_m)){
    a_index <- match(seq_a, colnames(jaccard_m))
    b_index <- match(seq_b, rownames(jaccard_m))
    ab_sim <- jaccard_similarity(kmers_a = kmers[[seq_a]], kmers_b = kmers[[seq_b]])
    ab_compl_sim <- jaccard_similarity(kmers_a = kmers[[seq_a]], kmers_b = kmers_compl[[seq_b]])
    jaccard_m[a_index, b_index] <- max(ab_sim, ab_compl_sim)
  }
}

# plot jaccard similarity matrix
jaccard_plt <- pheatmap(
  jaccard_m, fontsize = 5, method = "complete", 
  color = colorRampPalette(c("#440154FF", "#238A8DFF", "#FDE725FF"))(1000))

ggsave(
  plot = jaccard_plt, 
  filename = file.path(home, figures, paste0(species, "_jaccard_similarity.centr_all.viridis.pdf")),
  device = "pdf", units = "cm", width = 40, height = 40)