# Plot genome statistics for three germ-line genomes (per each chromosome): 
# chromosome size, repeat and satellite fraction, BUSCO statistics and genes
################################################################################
# load packages
require(dplyr)
require(ggplot2)
require(stringr)
require(valr)
require(viridis)
require(ggpubr)
################################################################################
# Functions
# read chromosome size files
load_chrom_sizes <- function(file, species){
  df <- read.csv(file, header = FALSE, sep = "\t")
  colnames(df) <- c("chrom", "size")
  df["species"] <- rep(species, length(df$chrom))
  df["chrom_annot"] <- paste(df$species, df$chrom, sep = "_")
  return(df)
}

# generate busco stats
generate_busco_stats <- function(file, sizes, species, total_busco_number){
  buscos <- read.csv(file, sep='\t', comment.char = '#', header = FALSE)[,c(0:5)]
  colnames(buscos) <- c("busco", "status", "chrom", "start", "end")
  agg <- aggregate(buscos$busco, by=list(buscos$chrom, buscos$status), FUN=length)
  colnames(agg) <- c("chrom", "status", "count")
  agg <- agg %>% filter(!status == "Missing")
  agg["species"] <- rep(species, length(agg$chrom))
  agg["chrom_annot"] <- paste(agg$species, agg$chrom, sep = "_")
  agg <- agg %>% filter(chrom %in% sizes$chrom) %>% arrange(chrom)
  agg["busco_frac"] <- round(agg$count / total_busco_number * 100, 2)
  return(agg)
}

# read and process TE files
generate_tes_stats <- function(file, split = FALSE, split_by = NULL, chrom_sizes, species){
  bed_tes <- read.table(file, sep = "\t", header = FALSE)
  colnames(bed_tes) <- c("chrom", "start", "end", "family", "score", "strand")
  
  if(isTRUE(split)){
    bed_tes[c('family', 'subfamily')] <- str_split_fixed(bed_tes$family, split_by, 2)
    bed_tes[bed_tes$subfamily == "Penelope",]$family <- "PLE"
  }
  
  # filter by chromosomes
  bed_tes <- bed_tes[bed_tes$chrom %in% chrom_sizes$chrom,]
  
  df <- NULL
  for(family in unique(bed_tes$family)){
    rep_df <- bed_tes[bed_tes$family == family,]
    rep_df <- bed_merge(rep_df)
    rep_df["family"] <- rep(family, length(rep_df$chrom))
    df <- rbind(df, rep_df)
  }
  
  df["length"] <- df$end - df$start + 1
  agg <- aggregate(df$length, by = list(df$chrom, df$family), FUN=sum)
  colnames(agg) <- c("chrom", "family", "sum")
  
  frac <- c()
  
  for (row in 1:nrow(agg)){
    chr <- agg[row, "chrom"]
    size <- chrom_sizes[chrom_sizes$chrom == chr, "size"]
    frac <- append(frac, round(agg[row, "sum"] / size * 100, 2))
  }
  
  agg$chrFrac <- frac
  agg <- agg %>% arrange(chrom)
  agg["species"] <- rep(species, length(agg$chrom))
  agg["chrom_annot"] <- paste(agg$species, agg$chrom, sep = "_")
  return(agg)
}

# read and process satellite files (TRASH results)
generate_satellites_stats <- function(file, chrom_sizes){
  satellites <- read.csv(file, header = TRUE)
  satellites <- satellites[,c("name", "start", "end", 
                              "most.freq.value.N", "consensus.primary", "consensus.count")]
  colnames(satellites) <- c("chrom", "start", "end", "rep", "consensus", "count")
  satellites <- satellites[!duplicated(satellites),]
  satellites <- bed_merge(satellites)
  satellites["width"] <- satellites$end - satellites$start + 1
  agg <- aggregate(satellites$width, by = list(satellites$chrom), FUN = sum)
  colnames(agg) <- c("chrom", "satellites_width")
  agg <- merge(chrom_sizes, agg, by = "chrom")
  agg["satellites_frac"] <- round(agg$satellites_width / agg$size * 100, 2)
  return(agg)
}

# read and process gene annotations
generate_gene_stats <- function(file, species, chrom_list){
  genes <- read.csv(file, header = FALSE, sep = "\t")[,c(0:3)]
  colnames(genes) <- c("chrom", "start", "end")
  genes <- genes %>% filter(chrom %in% chrom_list)
  agg <- aggregate(genes$chrom, by = list(genes$chrom), FUN = length)
  colnames(agg) <- c("chrom", "n_genes")
  agg["species"] <- rep(species, length(agg$chrom))
  agg["chrom_annot"] <- paste(agg$species, agg$chrom, sep = "_")  
  return(agg)
}
################################################################################
home <- getwd()
data <- "data"
figures <- "figures"

chrom_list <- c("SUPER_1", "SUPER_2", "SUPER_3", "SUPER_X", 
                "SUPER_GRC1", "SUPER_GRC2", "SUPER_GRC")

chrom_annot_list <- c(
  "Bcop_SUPER_1", "Bcop_SUPER_2", "Bcop_SUPER_3", "Bcop_SUPER_X", 
  "Bcop_SUPER_GRC1", "Bcop_SUPER_GRC2",
  "Bimp_SUPER_1", "Bimp_SUPER_2", "Bimp_SUPER_3", "Bimp_SUPER_X", 
  "Bimp_SUPER_GRC",
  "Ling_SUPER_1", "Ling_SUPER_2", "Ling_SUPER_3", "Ling_SUPER_X", 
  "Ling_SUPER_GRC1", "Ling_SUPER_GRC2")

diptera_total <- 4867
################################################################################
# Prepare tables
### --- read chromosome sizes --- ###
bcop_sizes <- load_chrom_sizes(
  file = file.path(home, data, "chrom_sizes/idBraCopr2.1.chrom_sizes.tsv"), 
  species = "Bcop")

bimp_sizes <- load_chrom_sizes(
  file = file.path(home, data, "chrom_sizes/idBraImpa2.1.primary.chrom_sizes.tsv"), 
  species = "Bimp")

ling_sizes <- load_chrom_sizes(
  file = file.path(home, data, "chrom_sizes/idLycInge5.1.primary.chrom_sizes.tsv"), 
  species = "Ling")

all_sizes <- rbind(bcop_sizes, bimp_sizes, ling_sizes)

### --- read and process busco files --- ### 
bcop_buscos <- generate_busco_stats(
  file = file.path(home, data, "buscos/BraCopr_buscos.diptera_odb12.tsv"), 
  sizes = bcop_sizes, species = "Bcop", total_busco_number = diptera_total)
bimp_buscos <- generate_busco_stats(
  file = file.path(home, data, "buscos/BraImpa_buscos.diptera_odb12.tsv"), 
  sizes = bimp_sizes, species = "Bimp", total_busco_number = diptera_total)
ling_buscos <- generate_busco_stats(
  file = file.path(home, data, "buscos/LycInge_buscos.diptera_odb12.tsv"), 
  sizes = ling_sizes, species = "Ling", total_busco_number = diptera_total)

all_buscos <- rbind(bcop_buscos, bimp_buscos, ling_buscos)

### --- read and process transposable elements --- ### 
bcop_tes <- generate_tes_stats(
  file = file.path(home, data, "TEs_bed_files/idBraCopr2.filteredRepeats.bed"), 
  split = TRUE, split_by = "/", chrom_sizes = bcop_sizes, species = "Bcop")
bimp_tes <- generate_tes_stats(
  file = file.path(home, data, "TEs_bed_files/idBraImpa2.filteredRepeats.bed"), 
  split = TRUE, split_by = "/", chrom_sizes = bimp_sizes, species = "Bimp")
ling_tes <- generate_tes_stats(
  file = file.path(home, data, "TEs_bed_files/idLyncInge5.filteredRepeats.bed"), 
  split = TRUE, split_by = "/", chrom_sizes = ling_sizes, species = "Ling")

all_tes <- rbind(bcop_tes, bimp_tes, ling_tes)
#all_tes <- all_tes %>% filter(!family %in% c("Low_complexity", "Satellite", "Simple_repeat"))
all_tes["family_annot"] <- all_tes$family
all_tes$family_annot[
  all_tes$family_annot %in% c("Low_complexity", "Satellite", "Simple_repeat")] <- "Other"

### --- read and process satellite annotations --- ###
bcop_satellites <- generate_satellites_stats(
  file = file.path(home, data, "trash_satellites/Summary.of.repetitive.regions.bcop.fa.csv"), 
  chrom_sizes = bcop_sizes)
bimp_satellites <- generate_satellites_stats(
  file = file.path(home, data, "trash_satellites/Summary.of.repetitive.regions.bimp.fa.csv"), 
  chrom_sizes = bimp_sizes)
ling_satellites <- generate_satellites_stats(
  file = file.path(home, data, "trash_satellites/Summary.of.repetitive.regions.lyco.fa.csv"),
  chrom_sizes = ling_sizes)
all_satellites <- rbind(bcop_satellites, bimp_satellites, ling_satellites)

### --- read and process gene annotations --- ###
bcop_core <- generate_gene_stats(
  file = file.path(home, data, "gene_annotations/bcop_core.gene.bed"), 
  species = "Bcop", chrom_list = bcop_sizes$chrom)
bcop_grc <- generate_gene_stats(
  file = file.path(home, data, "gene_annotations/bcop_grc.gene.bed"), 
  species = "Bcop", chrom_list = bcop_sizes$chrom)

bimp_core <- generate_gene_stats(
  file = file.path(home, data, "gene_annotations/bimp_core.gene.bed"), 
  species = "Bimp", chrom_list = bimp_sizes$chrom)
bimp_grc <- generate_gene_stats(
  file = file.path(home, data, "gene_annotations/bimp_grc.gene.bed"), 
  species = "Bimp", chrom_list = bimp_sizes$chrom)

ling_core <- generate_gene_stats(
  file = file.path(home, data, "gene_annotations/ling_core.gene.bed"), 
  species = "Ling", chrom_list = ling_sizes$chrom)
ling_grc <- generate_gene_stats(
  file = file.path(home, data, "gene_annotations/ling_grc.gene.bed"), 
  species = "Ling", chrom_list = ling_sizes$chrom)

all_genes <- rbind(bcop_core, bcop_grc,
                   bimp_core, bimp_grc, 
                   ling_core, ling_grc)
################################################################################
# plotting

### --- sizes --- ###
plt_sizes <- ggplot(all_sizes) +
  geom_col(aes(x = size/1000000, y = factor(chrom_annot, rev(chrom_annot_list)))) + 
  labs(x = "size (Mbp)", y = NULL) +
  theme(plot.title = element_text(size=10)) + 
  facet_wrap(~species, scale = "free_y", nrow = 3) +
  theme_bw()

### --- buscos -- ###

BUSCO_colors <-c("Fragmented" = "#ffdd55ff", "Duplicated" = "#56b4e9", "Complete" = "#cc79a7")

plt_buscos <- ggplot(all_buscos) +
  geom_col(aes(x = count, y = factor(chrom_annot, rev(chrom_annot_list)), fill = forcats::fct_rev(status))) + 
  labs(x = "# of BUSCOs", y = NULL) +
  theme(plot.title = element_text(size=10)) + 
  scale_fill_manual(values = BUSCO_colors, name = "Status") +
  facet_wrap(~species, scale = "free_y", nrow = 3) +
  theme_bw()

### --- transposable elements --- ###
TE_colors <- c("Unknown" = "lightgrey", "Other" = "#332288", "SINE" = "#88CCEE", "Satellite" = "#44AA99", 
               "RC" = "#999933", "PLE" = "#DDCC77", "LTR" = "#CC6677", 
               "LINE" = "#882255", "DNA" = "#AA4499")
TE_colors_order <- c("Unknown", "Other", "SINE", "RC", "PLE", "LTR", "LINE", "DNA")

all_tes$family_annot <- factor(all_tes$family_annot, levels = rev(TE_colors_order))

plt_tes <- ggplot(all_tes) +
  geom_col(aes(x = chrFrac, y = factor(chrom_annot, rev(chrom_annot_list)), fill = forcats::fct_rev(family_annot))) + 
  labs(x = "repeat fraction (%)", y = NULL) +
  theme(plot.title = element_text(size=10)) + 
  scale_fill_manual(values = TE_colors, name = "TE family") + 
  facet_wrap(~species, scale = "free_y", nrow = 3) +
  theme_bw()

### --- satellites --- ###
plt_satellites <- ggplot(all_satellites) +
  geom_col(aes(x = satellites_frac, y = factor(chrom_annot, rev(chrom_annot_list))), fill = "#44AA99") + 
  labs(x = "satellites (%)", y = NULL) +
  theme(plot.title = element_text(size=10)) + 
  facet_wrap(~species, scale = "free_y", nrow = 3) +
  theme_bw()

### --- genes --- ###
plt_genes <- ggplot(all_genes) +
  geom_col(aes(x = n_genes, y = factor(chrom_annot, rev(chrom_annot_list))), fill = "#4c5a79") + 
  labs(x = "# of genes", y = NULL) +
  theme(plot.title = element_text(size=10)) + 
  facet_wrap(~species, scale = "free_y", nrow = 3) +
  theme_bw()

### --- merge plots together --- ###
plt_all <- ggarrange(
  plt_sizes + theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
                    axis.title.x = element_text(size = 22), axis.text.x = element_text(size = 18)),
  plt_tes + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), 
                  axis.title.x = element_text(size = 22), axis.text.x = element_text(size = 18)),
  plt_satellites + theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
                         axis.title.x = element_text(size = 22), axis.text.x = element_text(size = 18)),
  plt_buscos + theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
                     axis.title.x = element_text(size = 22), axis.text.x = element_text(size = 18)),
  plt_genes + theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
                    axis.title.x = element_text(size = 22), axis.text.x = element_text(size = 18)),
  nrow = 1, legend = "none")

ggsave(plot = plt_all , filename = file.path(home, figures, "GRCs_chrom_stats.svg"), 
       device = "svg", units = "cm", width = 100, height = 40)

ggsave(plot = plt_all , filename = file.path(home, figures, "GRCs_chrom_stats.png"), 
       device = "png", units = "cm", width = 100, height = 40)
