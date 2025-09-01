# Calculate distribution of different features in windows relative to 
# the centromere position: transposable elements, genes
################################################################################
# Load packages
require("plyranges")
require("dplyr")
require("valr")
require("ggplot2")
require("scales")
require("ggpubr")
################################################################################
# Functions
# read chromosome sizes
read_chrom_sizes <- function(file){
  chromosomes <- read.table(file = file, sep = "\t", header = FALSE, col.names = c("chrom", "end"))
  chromosomes$start <- rep(0, length(chromosomes$chrom))
  chromosomes <- chromosomes[,c("chrom", "start", "end")]
  return(chromosomes)
}

# read and process TEs divergence files
read_and_process_te_gffs <- function(file, div_cutoff = 0.5){
  # Read in data
  divergence_eg_gff <- read_gff(file)
  
  # Breakdown classification of repeats
  divergence_eg_gff <- divergence_eg_gff %>%
    mutate(subclass = sub("/.", "", type),
           superfamily = sub("-.*", "", sub(".*/", "", type)))
  
  # fix Penelopes
  divergence_eg_gff <- divergence_eg_gff %>%
    dplyr::mutate(subclass = ifelse(superfamily == "Penelope", "PLE", subclass)) %>%
    dplyr::mutate(subclass = ifelse(
      subclass %in% c("DNA", "LINE", "LTR", "PLE", "RC", "SINE", "Satellite", "Unknown"), 
      subclass, "Other")) %>%
    dplyr::mutate(named_subclass = case_when(
      subclass == "DNA" ~ "DNA Transposon", subclass == "LTR" ~ "LTR Retrotransposon",
      subclass == "PLE" ~ "Penelope", subclass == "RC" ~ "Rolling Circle",
      .default = subclass))
  
  divergence_eg_gff_plt  <- divergence_eg_gff %>%
    filter(!is.na(KIMURA80)) %>%
    mutate(KIMURA80 = as.numeric(KIMURA80)) %>%
    dplyr::filter(KIMURA80 <= div_cutoff) %>%
    as_tibble() %>%
    dplyr::mutate(KIMURA80 = round(x = KIMURA80, digits = 2)) %>%
    group_by(named_subclass, KIMURA80) %>%
    mutate(KIMURA_SUM = sum(width)) %>%
    ungroup() %>%
    dplyr::select(seqnames, start, end, subclass, named_subclass, KIMURA80, KIMURA_SUM) %>%
    base::unique()
  
  colnames(divergence_eg_gff_plt) <- c("chrom", "start", "end", "subclass", "named_subclass",
                                       "KIMURA80", "KIMURA_SUM")
  return(divergence_eg_gff_plt)
}

# read and process genes gffs
read_and_process_gene_gffs <- function(file, chroms){
  genes <- read_gff(file)
  genes <- genes %>%
    as_tibble() %>%
    filter(seqnames %in% chroms) %>%
    filter(type == "gene") %>%
    select(seqnames, start, end, strand, ID)
  colnames(genes) <- c("chrom", "start", "end", "strand", "gene_id")
  return(genes)
}

# split chromosomes into arms based on centromere position 
chroms_to_arms <- function(sizes, centromeres, frac_threshold = 35){
  # join chromosome infor with centromeres
  chromosomes <- left_join(sizes, centromeres, by = "chrom")
  
  # left arm
  #chrom_arms_I <- chromosomes[,c("chrom", "start", "centrStart", "end")]
  chrom_arms_I <- chromosomes[,c("chrom", "start", "centrMid", "end")]
  colnames(chrom_arms_I) <- c("chrom", "start", "end", "chrom_size")
  chrom_arms_I["arm"] <- rep("I", length(chrom_arms_I$chrom))
  
  # right arm
  #chrom_arms_II <- chromosomes[,c("chrom", "centrEnd", "end", "end")]
  chrom_arms_II <- chromosomes[,c("chrom", "centrMid", "end", "end")]
  colnames(chrom_arms_II) <- c("chrom", "start", "end", "chrom_size")
  chrom_arms_II["arm"] <- rep("II", length(chrom_arms_II$chrom))
  
  chrom_arms <- rbind(chrom_arms_I, chrom_arms_II) %>% arrange(chrom, start) 
  chrom_arms["length"] <- chrom_arms$end - chrom_arms$start + 1
  chrom_arms["frac"] <- chrom_arms$length / chrom_arms$chrom_size * 100
  
  # remove arms that occupy less than 35% of the chromosome
  chrom_arms <- chrom_arms %>% filter(frac >= frac_threshold)
  
  return(chrom_arms)
}

# calculate coverage in windows
calculate_coverage <- function(regions, features_df, type){
  cov <- bed_coverage(regions, features_df)
  cov["feature"] <- rep(type, length(cov$chrom))
  return(cov)
}

read_gc_content <- function(file, intervals){
  gc <- read.csv(file, header = TRUE, sep = "\t")
  gc <- gc[,c(1:5)]
  colnames(gc) <- c("chrom", "start", "end", "at", ".frac")
  gc <- gc %>% filter(chrom %in% chrom_arms$chrom) %>% select(chrom, start, end, .frac)
  gc["feature"] <- rep("GC", length(gc$chrom))
  gc <- left_join(intervals, gc, by = c("chrom", "start", "end"))
  return(gc)
}

################################################################################
# Main directories
#home <- "/Users/ab66/Documents/sanger_work/diptera/analysis_on_curated_genomes"
home <- getwd()
data <- "data"
figures <- "figures"
################################################################################
# Load files

### --- chromosome sizes --- ###
sizes <- read_chrom_sizes(
  #file = file.path(home, data, "chrom_sizes/idBraCopr2.1.chrom_sizes.tsv")) # bcop
  #file = file.path(home, data, "chrom_sizes/idBraImpa2.1.primary.chrom_sizes.tsv")) # bimp
  file = file.path(home, data, "chrom_sizes/idLycInge5.1.primary.chrom_sizes.tsv")) # bimp


### --- centromeres --- ###
centromeres <- read.table(
  #file = file.path(home, data, "centromeres/Bcop.centromeres.bed"),
  #file = file.path(home, data, "centromeres/Bimp.centromeres.bed"),
  file = file.path(home, data, "centromeres/Ling.centromeres.bed"),
  sep = "\t", header = FALSE, col.names = c("chrom", "centrStart", "centrEnd"))

centromeres["centrMid"] <- 
  ((centromeres$centrEnd - centromeres$centrStart) / 2) + centromeres$centrStart

### --- transposable elements --- ###
tes_dir <- file.path(home, data, "TEs_kimura_gff")
TE_files <- c(
  
  # bcop
  #"idBraCopr2.filteredRepeats.autosomes.out.gff", "idBraCopr2.filteredRepeats.chrX.out.gff",
  #"idBraCopr2.filteredRepeats.GRC1.out.gff", "idBraCopr2.filteredRepeats.GRC2.out.gff")
  
  # bimp
  #"idBraImpa2.filteredRepeats.autosomes.out.gff", "idBraImpa2.filteredRepeats.chrX.out.gff",
  #"idBraImpa2.filteredRepeats.GRC.out.gff")
  
  # ling
  "idLyncInge5.filteredRepeats.autosomes.out.gff", "idLyncInge5.filteredRepeats.chrX.out.gff",
  "idLyncInge5.filteredRepeats.GRC1.out.gff", "idLyncInge5.filteredRepeats.GRC2.out.gff")

div_tes_all <-  do.call(
  "rbind", lapply(file.path(tes_dir, TE_files), read_and_process_te_gffs))

### --- genes --- ###
gene_dir <- file.path(home, data, "gene_annotations")
gene_files <- c(
  #"bcop_core.gff3", "bcop_grc.gff3") # bcop
  #"bimp_core.gff3", "bimp_grc.gff3") # bimp
  "ling_core.gff3", "ling_grc.gff3") # ling

genes <- do.call(
  "rbind", lapply(file.path(gene_dir, gene_files), read_and_process_gene_gffs, chroms = sizes$chrom))

################################################################################
# Process files
# split chromosomes into arms based on centromere position
chrom_arms <- chroms_to_arms(sizes = sizes, centromeres = centromeres)

# split chromosome arms into equal intervals
### micro intervals ###
chr_regions_small <- bed_makewindows(chrom_arms, num_win = 100)

# change order in arm I - distance from the centromere
chr_regions_small <- chr_regions_small %>%
  filter(arm == "I") %>%
  group_by(chrom) %>%
  mutate(.win_id = c(-1:-100)) %>%
  ungroup() %>%
  rbind(chr_regions_small[chr_regions_small$arm == "II",]) %>%
  arrange(chrom, start) %>%
  select(chrom, start, end, arm, .win_id)

# calculate features cov distribution in the intervals
feature_cov <- list()
feature_cov$tes <- calculate_coverage(
  regions = chr_regions_small, features_df = div_tes_all, type = "TEs")
feature_cov$genes <- calculate_coverage(
  regions = chr_regions_small, features_df = genes, type = "genes")

# bind TEs and genes dfs
feature_cov_small <- rbind(feature_cov$tes, feature_cov$genes) %>% 
  select(chrom, start, end, arm, .win_id, .frac, feature)

# infer the order of plotting
feature_cov_small$chrom <- factor(
  feature_cov_small$chrom, 
  levels = c("SUPER_1",
             "SUPER_2",
             "SUPER_3",
             "SUPER_X",
             "SUPER_GRC1",
             "SUPER_GRC2"))
             #"SUPER_GRC"))

features_colour <- c("genes" = "#4c5a79", "TEs" = "grey25")

# plot as lines
plt_lines <- feature_cov_small %>%
  #filter(chrom %in% c("SUPER_1", "SUPER_GRC1", "SUPER_GRC2")) %>%
  ggplot(aes(x = .win_id / 100, y = .frac * 100, col = feature)) +
  geom_point(alpha = 0.2, size = 0.3) +
  geom_smooth(span = 0.5, se = TRUE, aes(group=feature, col = feature, fill = feature)) +
  labs(x = "Distance", y = "Fraction (%)") +
  scale_color_manual(values = features_colour, name = "Feature") +
  scale_fill_manual(values = features_colour, name = "Feature") +
  geom_vline(xintercept = 0, linetype="dashed", color = "darkred", size = 1) +
  #geom_line(data = feature_cov$gc, aes(x = .win_id, y = .frac * 100)) +
  scale_x_continuous(labels = label_number(accuracy = 0.1)) + 
  ylim(0,60) +
  #facet_wrap(.~chrom, strip.position = "left", ncol = 1) + # chrom-scale
  facet_wrap(.~chrom, scales = "free_x", strip.position = "left", ncol = 1) + # arm-scale
  #ggtitle("Bradysia comprophila") +
  #ggtitle("Bradysia impatiens") +
  ggtitle("Lycoriella ingenua") +
  theme_bw()

ggsave(file.path(home, figures,
                 #"Bcop_chrom_organization_lines_vertical.png"),
                 #"Bimp_chrom_organization_lines_vertical.png"),
                 "Ling_chrom_organization_lines_vertical.png"),
       plt_lines, width = 5, height = 10)
ggsave(file.path(home, figures, 
                 #"Bcop_chrom_organization_lines_vertical.svg"),
                 #"Bimp_chrom_organization_lines_vertical.svg"),
                 "Ling_chrom_organization_lines_vertical.svg"),
       plt_lines, width = 5, height = 10)