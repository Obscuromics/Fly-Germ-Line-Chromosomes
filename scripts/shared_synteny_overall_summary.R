# Plot stats for shared synteny blocks across chromosomes and species
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

# load and process alignments within a species
generate_sb_stats <- function(file, chrom_sizes, species){
  ali <- read.table(file, header = FALSE)
  colnames(ali) <- c("query", "q_length", "q_start", "q_end", "q_string",
                     "target", "t_length", "t_start", "t_end", "n_matches", 
                     "aln_length", "map_quality", "dv", "df")
  ali <- ali %>% filter(query %in% chrom_sizes$chrom & target %in% chrom_sizes$chrom)
  
  sb_stats = NULL
  
  for(chrom in unique(chrom_sizes$chrom)){
    ali_flt <- ali %>% filter(query == chrom & !(target == chrom)) %>%
      select(query, q_start, q_end, target)
    colnames(ali_flt) <- c("chrom", "start", "end", "target")
    for(chrom_target in unique(ali_flt$target)){
      df <- ali_flt[ali_flt$target == chrom_target,]
      df <- bed_merge(df)
      df["width"] <- df$end - df$start + 1
      df <- df %>% filter(width >= 1000)
      target_sum <- sum(df$width)
      if(chrom_target %in% c("SUPER_GRC1", "SUPER_GRC2", "SUPER_GRC")){
        target_annot <- paste(species, chrom_target, sep = "_") 
      }else{target_annot <- chrom_target}
      sb_stats <- rbind(sb_stats, data.frame(chrom = chrom, target = chrom_target, 
                                             target_annot = target_annot, width_sum = target_sum))
    }
    no_sb_length <- chrom_sizes[chrom_sizes$chrom == chrom,]$size - 
      sum(sb_stats[sb_stats$chrom == chrom,]$width_sum)
    sb_stats <- rbind(sb_stats, data.frame(chrom = chrom, target = "Other", target_annot = "Other", 
                                           width_sum = no_sb_length))
  }
  sb_stats["species"] <- rep(species, length(sb_stats$chrom))
  sb_stats["chrom_annot"] <- paste(sb_stats$species, sb_stats$chrom, sep = "_")
  return(sb_stats)
}

# load and process alignments across sciarid species
generate_sb_diff_sp_stats <- function(file, sp1_chrom_list, sp2_chrom_list, sp1, sp2, rev = FALSE){
  ali <- read.table(file, header = FALSE)
  colnames(ali) <- c("query", "q_length", "q_start", "q_end", "q_string",
                     "target", "t_length", "t_start", "t_end", "n_matches", 
                     "aln_length", "map_quality", "dv", "df")
  if(rev){
    ali <- ali %>% filter(query %in% sp2_chrom_list & target %in% sp1_chrom_list) %>%
      select(target, t_start, t_end, query)
    colnames(ali) <- c("query", "q_start", "q_end", "target")
    #print(ali)
  }else{
    ali <- ali %>% filter(query %in% sp1_chrom_list & target %in% sp2_chrom_list) %>%
      select(query, q_start, q_end, target)
    #print(ali)
  }
  
  sb_stats <- NULL
  
  for(chrom in unique(ali$query)){
    ali_flt <- ali %>% filter(query == chrom)
    colnames(ali_flt) <- c("chrom", "start", "end", "target")
    for(chrom_target in unique(ali_flt$target)){
      df <- ali_flt[ali_flt$target == chrom_target,]
      df <- bed_merge(df)
      df["width"] <- df$end - df$start + 1
      df <- df %>% filter(width >= 1000)
      target_sum <- sum(df$width)
      if(chrom_target %in% c("SUPER_GRC1", "SUPER_GRC2", "SUPER_GRC")){
        target_annot <- paste(sp2, chrom_target, sep = "_") 
      }else{target_annot <- chrom_target}
      sb_stats <- rbind(sb_stats, data.frame(chrom = chrom, target = chrom_target, 
                                             target_annot = target_annot, 
                                             width_sum = target_sum))
    }
  }
  
  sb_stats["species"] <- rep(sp1, length(sb_stats$chrom))
  sb_stats["chrom_annot"] <- paste(sb_stats$species, sb_stats$chrom, sep = "_")
  return(sb_stats)
}

# load and process alignments with aaph
generate_sb_with_aaph_stats <- function(file, chrom_sizes, species){
  ali <- read.table(file, header = FALSE)
  colnames(ali) <- c("query", "q_length", "q_start", "q_end", "q_string",
                     "target", "t_length", "t_start", "t_end", "n_matches", 
                     "aln_length", "map_quality", "dv", "df")
  
  ali <- ali %>% filter(target %in% chrom_sizes$chrom) %>%
    select(target, t_start, t_end, query)
  colnames(ali) <- c("query", "q_start", "q_end", "target")
  ali[which(ali$target %in% c("CM059996.1", "CM059999.1")), "target_annot"] <- "Aaph_chrA"
  ali[which(ali$target %in% c("CM059997.1", "CM059998.1")), "target_annot"] <- "Aaph_chrX"
  
  sb_stats = NULL
  for(chrom in unique(chrom_sizes$chrom)){
    ali_flt <- ali %>% filter(query == chrom)
    colnames(ali_flt) <- c("chrom", "start", "end", "target", "target_annot")
    for(chrom_target in unique(ali_flt$target)){
      df <- ali_flt[ali_flt$target == chrom_target,]
      target_annot <- unique(df$target_annot)
      df <- bed_merge(df)
      df["width"] <- df$end - df$start + 1
      df <- df %>% filter(width >= 1000)
      target_sum <- sum(df$width)
      sb_stats <- rbind(sb_stats, data.frame(chrom = chrom, target = chrom_target, 
                                             target_annot = target_annot, width_sum = target_sum))
    }
  }
  
  chrom_pairs <- data.frame(
    chrom = c(rep("SUPER_1", 4), rep("SUPER_2", 4), rep("SUPER_3", 4), rep("SUPER_X", 4)),
    target = c(rep(c("CM059996.1", "CM059997.1", "CM059998.1", "CM059999.1"), 4)))
  
  if(species %in% c("Bcop", "Ling")){
    chrom_pairs <- rbind(chrom_pairs, data.frame(
      chrom = c(rep("SUPER_GRC1", 4), rep("SUPER_GRC2", 4)),
      target = c(rep(c("CM059996.1", "CM059997.1", "CM059998.1", "CM059999.1"), 2))))
  }else if(species == "Bimp"){
    chrom_pairs <- rbind(chrom_pairs, data.frame(
      chrom = c(rep("SUPER_GRC", 4)),
      target = c("CM059996.1", "CM059997.1", "CM059998.1", "CM059999.1")))
  }
  
  sb_stats <- left_join(chrom_pairs, sb_stats)
  
  sb_stats[which(sb_stats$target %in% c("CM059996.1", "CM059999.1")), "target_annot"] <- "Aaph_chrA"
  sb_stats[which(sb_stats$target %in% c("CM059997.1", "CM059998.1")), "target_annot"] <- "Aaph_chrX"
    
  sb_stats["species"] <- rep(species, length(sb_stats$chrom))
  sb_stats["chrom_annot"] <- paste(sb_stats$species, sb_stats$chrom, sep = "_")
  return(sb_stats)
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

### --- GRCs to core within the same genome --- ###
bcop_sbs <- generate_sb_stats(
  file = file.path(home, data, "alignments/Bcop_Bcop.m.1aln.paf"), 
  chrom_sizes = bcop_sizes, species = "Bcop")
bimp_sbs <- generate_sb_stats(
  file = file.path(home, data, "alignments/Bimp_Bimp.m.1aln.paf"), 
  chrom_sizes = bimp_sizes, species = "Bimp")
ling_sbs <- generate_sb_stats(
  file = file.path(home, data, "alignments/Ling_Ling.m.1aln.paf"), 
  chrom_sizes = ling_sizes, species = "Ling")

sb_self <- rbind(bcop_sbs, bimp_sbs, ling_sbs)

### --- genome alignments across species --- ###
bcop_bimp <- generate_sb_diff_sp_stats(
  file = file.path(home, data, "alignments/Bcop_Bimp.m.1aln.paf"), 
  sp1_chrom_list = bcop_sizes$chrom, 
  sp2_chrom_list = c("SUPER_GRC"), sp1 = "Bcop", sp2 = "Bimp", rev = FALSE)

bimp_bcop <- generate_sb_diff_sp_stats(
  file = file.path(home, data, "alignments/Bcop_Bimp.m.1aln.paf"), 
  sp1_chrom_list = bimp_sizes$chrom, 
  sp2_chrom_list = c("SUPER_GRC1", "SUPER_GRC2"), sp1 = "Bimp", sp2 = "Bcop", rev = TRUE)

bcop_ling <- generate_sb_diff_sp_stats(
  file = file.path(home, data, "alignments/Bcop_Ling.m.1aln.paf"), 
  sp1_chrom_list = bcop_sizes$chrom, 
  sp2_chrom_list = c("SUPER_GRC1", "SUPER_GRC2"), sp1 = "Bcop", sp2 = "Ling", rev = FALSE)

ling_bcop <- generate_sb_diff_sp_stats(
  file = file.path(home, data, "alignments/Bcop_Ling.m.1aln.paf"), 
  sp1_chrom_list = ling_sizes$chrom, 
  sp2_chrom_list = c("SUPER_GRC1", "SUPER_GRC2"), sp1 = "Ling", sp2 = "Bcop", rev = TRUE)

bimp_ling <- generate_sb_diff_sp_stats(
  file = file.path(home, data, "alignments/Bimp_Ling.m.1aln.paf"), 
  sp1_chrom_list = bimp_sizes$chrom, 
  sp2_chrom_list = c("SUPER_GRC1", "SUPER_GRC2"), sp1 = "Bimp", sp2 = "Ling", rev = FALSE)

ling_bimp <- generate_sb_diff_sp_stats(
  file = file.path(home, data, "alignments/Bimp_Ling.m.1aln.paf"), 
  sp1_chrom_list = ling_sizes$chrom, 
  sp2_chrom_list = c("SUPER_GRC"), sp1 = "Ling", sp2 = "Bimp", rev = TRUE)

sbs_diff_sp <- rbind(bcop_bimp, bimp_bcop, 
                     bcop_ling, ling_bcop, 
                     bimp_ling, ling_bimp)

### --- sciarid chromosomes to aaph genome --- ###
aaph_bcop_sbs <- generate_sb_with_aaph_stats(
  file = file.path(home, data, "alignments/Aaph_Bcop.m.1aln.paf"), 
  chrom_sizes = bcop_sizes, species = "Bcop")

aaph_bimp_sbs <- generate_sb_with_aaph_stats(
  file = file.path(home, data, "alignments/Aaph_Bimp.m.1aln.paf"), 
  chrom_sizes = bimp_sizes, species = "Bimp")

aaph_ling_sbs <- generate_sb_with_aaph_stats(
  file = file.path(home, data, "alignments/Aaph_Ling.m.1aln.paf"), 
  chrom_sizes = ling_sizes, species = "Ling")

sbs_with_aaph <- rbind(aaph_bcop_sbs, aaph_bimp_sbs, aaph_ling_sbs)


################################################################################
# Plotting
sb_colours <- c(
  "SUPER_1" = "#4cceaf", "SUPER_2" = "#03a487", 
  "SUPER_3" = "#007c61", "SUPER_X" = "#00563e",
  "Aaph_chrA" = "#ce8eda", "Aaph_chrX" = "#a86bb4",
  "Bcop_SUPER_GRC1" = "#ff8ec5", "Bcop_SUPER_GRC2" = "#d4679e", 
  "Bimp_SUPER_GRC" = "#ffdb6b", 
  "Ling_SUPER_GRC1" = "#00a7d1", "Ling_SUPER_GRC2" = "#0071cc",
  "Other" = "lightgrey")

# plotting order
chr_order <- c("Ling_SUPER_GRC1", "Ling_SUPER_GRC2", 
               "Bimp_SUPER_GRC",
               "Bcop_SUPER_GRC1", "Bcop_SUPER_GRC2", 
               "Aaph_chrA", "Aaph_chrX", 
               "SUPER_1", "SUPER_2", "SUPER_3", "SUPER_X", "Other")

### --- GRCs to core within the same genome --- ###
plt_sbs_self <- sb_self %>% 
  filter(!target == "Other") %>%
  ggplot() +
  geom_col(aes(x = width_sum/1000000, 
               y = factor(chrom_annot, rev(chrom_annot_list)), 
               fill = factor(target_annot, names(sb_colours)))) + 
  labs(x = "sum length of SBs (Mbp)", y = NULL) +
  theme(plot.title = element_text(size=10)) +
  scale_fill_manual(values = sb_colours, name = "Chromosome") +
  xlim(0,20) + 
  facet_wrap(~species, scale = "free_y", nrow = 3) +
  theme_bw()

### --- genome alignments across species --- ###
plt_sbs_diff_sp <- sbs_diff_sp %>% 
  filter(!target == "Other") %>%
  ggplot() +
  geom_col(aes(x = width_sum/1000000, y = factor(chrom_annot, rev(chrom_annot_list)), 
               fill = factor(target_annot, names(sb_colours)))) + 
  labs(x = "sum length of SBs (Mbp)", y = NULL) +
  theme(plot.title = element_text(size=10)) +
  scale_fill_manual(values = sb_colours, name = "Chromosome") +
  xlim(0,20) +
  facet_wrap(~species, scale = "free_y", nrow = 3) +
  theme_bw()

### --- sciarid chromosomes to aaph genome --- ###
plt_sbs_with_aaph <- sbs_with_aaph %>% 
  filter(!target == "Other") %>%
  ggplot() +
  geom_col(aes(x = width_sum/1000000, 
               y = factor(chrom_annot, rev(chrom_annot_list)), 
               fill = factor(target_annot, names(sb_colours)))) + 
  labs(x = "sum length of SBs (Mbp)", y = NULL) +
  theme(plot.title = element_text(size=10)) +
  scale_fill_manual(values = sb_colours, name = "Chromosome") +
  xlim(0,20) +
  facet_wrap(~species, scale = "free_y", nrow = 3) +
  theme_bw()

### --- merge plots together --- ###
plt_all <- ggarrange(
  plt_sbs_self + theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
                       axis.title.x = element_text(size = 22), axis.text.x = element_text(size = 18)),
  plt_sbs_diff_sp + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), 
                          axis.title.x = element_text(size = 22), axis.text.x = element_text(size = 18)),
  plt_sbs_with_aaph + theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
                            axis.title.x = element_text(size = 22), axis.text.x = element_text(size = 18)),
  nrow = 1, legend = "none")

ggsave(plot = plt_all , filename = file.path(home, figures, "synteny_block_summary.svg"), 
       device = "svg", units = "cm", width = 80, height = 40)
ggsave(plot = plt_all , filename = file.path(home, figures, "synteny_block_summary.png"), 
       device = "png", units = "cm", width = 80, height = 40)