# Plot results of blasting satellites from short arm of the chrX to GRCs
################################################################################
require(valr)
require(ggplot2)
require(dplyr)
require(plyranges)
library(karyoploteR)
################################################################################
# Functions

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

# read chromosome files
load_chrom_sizes <- function(file, species){
  df <- read.csv(file, header = FALSE, sep = "\t")
  colnames(df) <- c("chrom", "size")
  df["species"] <- rep(species, length(df$chrom))
  return(df)
}

# convert genomic coordinates to linear format
prepare_linear_coordinates <- function(ali, t_sp, q_sp, chrom_list){
  # read chromosome sizes
  chrom_sizes <- list()
  for(sp in c(t_sp, q_sp)){
    chrom_sizes[[sp]] <- load_chrom_sizes(file = chrom_list[[sp]], species = sp)
  }
  
  # merge chromosomes sizes with the alignment table
  ali <- left_join(ali, chrom_sizes[[t_sp]], 
                   by = c("target" = "chrom", "t_species" = "species"))
  
  ali <- left_join(ali, chrom_sizes[[q_sp]],
                   by = c("query" = "chrom", "q_species" = "species"))
  
  names(ali)[names(ali) == "size.x"] <- "t_chrom_size"
  names(ali)[names(ali) == "size.y"] <- "q_chrom_size"
  
  # transform coordinates of the target to linear format
  ali <- ali %>% arrange(target, t_start)
  
  chr_offset <- 0
  ali_lin <- NULL
  
  for(i in unique(ali$target)){
    ali_t <- ali[ali$target == i,]
    size_t <- unique(ali_t$t_chrom_size)
    
    ali_t$t_linear_start <- chr_offset + ali_t$t_start
    ali_t$t_linear_end <- chr_offset + ali_t$t_end
    ali_lin <- rbind(ali_lin, ali_t)
    
    chr_offset <- chr_offset + size_t
  }
  
  # transform coordinates of the query to linear format
  ali_lin <- ali_lin %>% arrange(query, q_start)
  
  chr_offset <- 0
  ali_lin_both <- NULL
  
  for(i in unique(ali_lin$query)){
    ali_q <- ali_lin[ali_lin$query == i,]
    size_q <- unique(ali_q$q_chrom_size)
    
    ali_q$q_linear_start <- chr_offset + ali_q$q_start
    ali_q$q_linear_end <- chr_offset + ali_q$q_end
    ali_lin_both <- rbind(ali_lin_both, ali_q)
    
    chr_offset <- chr_offset + size_q
  }
  
  ali_lin_both <- ali_lin_both %>% arrange(target, t_linear_start) %>%
    select(target, t_linear_start, t_linear_end, t_species, t_chrom_size, t_string,
           query, q_linear_start, q_linear_end, q_species, q_chrom_size)
  
  return(ali_lin_both)
}

# calculate break positions for plotting
zipWithNext <- function(x, step = 1) {
  Pairs(
    x,
    c(tail(x, step * -1), rep(NA, step))
  )
}

calcLabelPosition <- function(breakPos) {
  # Return the position of the midpoint of each chr
  # in the context of the merged object
  breakPos |> zipWithNext() |> as.data.frame() |>
    rowMeans() |> head(-1)
}

# read BUSCO files
read_busco_file <- function(file_name, species, buscos_to_origin, chrom_list){
  df <- read.csv(file_name, sep = '\t', comment.char = '#', header = FALSE,
                 na.strings = c("", "NA"))[,c(0:6)]
  colnames(df) <- c("busco", "status", "chr", "start", "end", "strand")
  
  # swap start and end for buscos on "-" strand
  df_new <- df %>% filter(strand == "+")
  df_new <- rbind(
    df_new, df %>% dplyr::filter(strand == "-") %>% 
      dplyr::rename(start = end, end = start)) %>%
    arrange(chr, start)
  df <- df_new
  # merge buscos with origin
  buscos_to_origin <- buscos_to_origin %>% 
    filter(sp == species) %>% select(busco, chr, start, end, origin, sp)
  df <- left_join(df, buscos_to_origin)
  df[which(is.na(df$origin)), "origin"] <- "Other"
  df[which(is.na(df$sp)), "sp"] <- species
  
  # read chromosome sizes and merge with busco df
  chrom_sizes <- load_chrom_sizes(file = chrom_list[[species]], species = species)
  df <- left_join(df, chrom_sizes, by = c("chr" = "chrom", "sp" = "species"))
  
  # add colour
  df[which(df$origin == "Sciaridae"), "colour"] <- "#44AA99"
  df[which(df$origin == "Cecidomyiidae"), "colour"] <- "#AA4499"
  df[which(df$origin == "Other"), "colour"] <- "#AA4499" <- "grey90"
  
  colnames(df) <- c("busco", "status", "chrom", "start", 
                    "end", "strand", "origin", "sp", "chrom_size", "colour")
  return(df)
}

# read and process TEs divergence files
read_and_process_te_gffs <- function(file, div_cutoff = 0.5){
  # Read in data
  divergence_eg_gff <- read_gff(file)
  
  # Breakdown classification of repeats
  divergence_eg_gff <- divergence_eg_gff %>%
    mutate(subclass = sub("/.*", "", type),
           superfamily = sub("-.*", "", sub(".*/", "", type)))
  
  # fix Penelopes
  divergence_eg_gff <- divergence_eg_gff %>%
    dplyr::mutate(subclass = ifelse(superfamily == "Penelope", "PLE", subclass))%>%
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
    dplyr::select(seqnames, start, end, subclass, superfamily,
                  named_subclass, KIMURA80, KIMURA_SUM) %>%
    base::unique()
  
  colnames(divergence_eg_gff_plt) <- c("chrom", "start", "end", "subclass",
                                       "superfamily",
                                       "named_subclass",
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

read_chrom_sizes <- function(file){
  chromosomes <- read.table(file = file, sep = "\t", header = FALSE, col.names = c("chrom", "end"))
  chromosomes$start <- rep(0, length(chromosomes$chrom))
  chromosomes <- chromosomes[,c("chrom", "start", "end")]
  return(chromosomes)
}

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

target_species <- "Bcop"
target_chrom2plot <- c("SUPER_GRC2")

query_species <- "Bcop"
query_chrom2plot <- c("SUPER_GRC1")

chrom_files <- c("idBraCopr2.1.chrom_sizes.tsv",
                 "idBraImpa2.1.primary.chrom_sizes.tsv",
                 "idLycInge5.1.primary.chrom_sizes.tsv")

chrom_list <- file.path(home, data, "chrom_sizes", chrom_files)
names(chrom_list) <- c("Bcop", "Bimp", "Ling")

busco_files <- c("BraCopr_buscos.diptera_odb12.tsv",
                 "BraImpa_buscos.diptera_odb12.tsv",
                 "LycInge_buscos.diptera_odb12.tsv")

busco_list <- file.path(home, data, "buscos", busco_files)
names(busco_list) <- c("Bcop", "Bimp", "Ling")

TE_files <- c(
    "TEs_kimura_gff/idBraCopr2.filteredRepeats.GRC1.out.gff", 
    "TEs_kimura_gff/idBraCopr2.filteredRepeats.GRC2.out.gff")

gene_annotation_files <- c(
  "gene_annotations/bcop_core.gff3",
  "gene_annotations/bcop_grc.gff3")

alignment_file <- "alignments/Bcop_Bcop.m.1aln.paf"
################################################################################
### --- load paf alignment --- ###
ali <- load_paf_file(
  file = file.path(home, data, alignment_file),
  t_chroms = target_chrom2plot, t_sp = target_species, 
  q_chroms = query_chrom2plot, q_sp = query_species)

# zoom in to the synteny block
sb <- ali %>% 
  filter(t_start >= 39192040 & t_end <= 42263740) %>%
  filter(q_start >= 62068860 & q_end <= 65694342)

print(paste0("Number of alignments: ", nrow(sb)))
print(paste0("Alignmnt length: ", sum(sb$n_bases)/1000000, " (Mbp)"))
print(paste0("Identity: ", round(sum(sb$n_matches)/sum(sb$n_bases) * 100, 2), "%"))

### --- plot the synteny block --- ###
p_zoom <- ggplot(data.frame(
  start = sb$t_start,
  end = sb$t_end,
  strand = sb$t_string,
  query.start = ifelse(sb$t_string == "+", sb$q_start, sb$q_end),
  query.end = ifelse(sb$t_string == "+", sb$q_end, sb$q_start)
)) +
  
  aes(x = start/1000, y = query.start/1000, xend = end/1000, yend = query.end/1000)

p_zoom <- p_zoom + geom_segment(lineend = "round", linewidth = 0.4) +
  labs(x = paste0(unique(sb$target), " (Kbp)"), 
       y = paste0(unique(sb$query), " (Kbp)")) +
  theme_bw() +
  coord_fixed()
################################################################################
### --- satellite repeats --- ###
satellites <- read_trash_raw_file(
  file = file.path(home, data, "trash_satellites/all.repeats.from.bcop.fa.csv"))
satellites_GRC1 <- satellites %>%
  filter(chrom == "SUPER_GRC1" & start >= 62068860 & end <= 65694342)
satellites_GRC2 <- satellites %>%
  filter(chrom == "SUPER_GRC2" & start >= 39192040 & end <= 42263740)

satellites_to_plot <- rbind(satellites_GRC1, satellites_GRC2)
satellites_to_plot <- left_join(
  satellites_to_plot, data.frame(rep = c(94, 30, 27),
                                 col = c("#ff88bb", "#7489d9", "#007c98")))
#satellites_to_plot[which(is.na(satellites_to_plot$col)),]$col <- "grey85"

sat_grc1 <- satellites_to_plot %>% filter(chrom == "SUPER_GRC1")
sat_grc2 <- satellites_to_plot %>% filter(chrom == "SUPER_GRC2")

### --- transposable elements --- ###
te_colours <- data.frame(
  "subclass" = c("Unknown", "Other", "SINE", "Satellite", 
                 "RC", "PLE", "LTR", "LINE", "DNA"),
  "col" = c("lightgrey", "#332288", "#88CCEE", "#44AA99",
            "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499"))

div_tes_all <-  do.call(
  "rbind", lapply(file.path(home, data, TE_files), read_and_process_te_gffs))

div_tes_all <- left_join(div_tes_all, te_colours)

tes_grc1 <- div_tes_all %>% 
  filter(chrom == "SUPER_GRC1" &  start >= 62068860 & end <= 65694342)

tes_grc2 <- div_tes_all %>% 
  filter(chrom == "SUPER_GRC2" & start >= 39192040 & end <= 42263740)

### --- genes --- ###
sizes <- read_chrom_sizes(
  file = file.path(home, data, "chrom_sizes/idBraCopr2.1.chrom_sizes.tsv")) %>%
  filter(chrom %in% c("SUPER_GRC1", "SUPER_GRC2"))

genes <- do.call(
  "rbind", lapply(file.path(home, data, gene_annotation_files), 
                  read_and_process_gene_gffs, chroms = sizes$chrom))

genes_grc1 <- genes %>% 
  filter(chrom == "SUPER_GRC1" &  start >= 62068860 & end <= 65694342)
genes_grc2 <- genes %>% 
  filter(chrom == "SUPER_GRC2" & start >= 39192040 & end <= 42263740)
################################################################################
# plotting
df <- data.frame(chr = c("GRC1", "GRC2"),
                 start = c(62068860, 39192040),
                 end = c(65694342, 42263740))

sb[which(sb$t_string == "+"), "t_col"] <- "#f79943ff"
sb[which(sb$t_string == "-"), "t_col"] <- "#8f83b9ff"

custome.genome <- toGRanges(df)

pdf(file.path(home, figures, "Bcop_GRC1_to_GRC2_PAR_content.pdf"))

kp <- plotKaryotype(genome = custome.genome)
kpAddBaseNumbers(kp, tick.dist = 500000)

# satellites
kpRect(kp, chr = "GRC1", x0 = sat_grc1$start, x1 = sat_grc1$end,
       y0 = 0, y1 = 0.1, border = sat_grc1$col, col = sat_grc1$col, lwd = 0.3)
kpRect(kp, chr = "GRC1", x0 = df$start[df$chr == "GRC1"], x1 = df$end[df$chr == "GRC1"],
       y0 = 0, y1 = 0.1, border = "black", lwd = 0.3)

kpRect(kp, chr = "GRC2", x0 = grc2$start, x1 = grc2$end,
       y0 = 0, y1 = 0.1, border = grc2$col, col = grc2$col, lwd = 0.3)
kpRect(kp, chr = "GRC2", x0 = df$start[df$chr == "GRC2"], x1 = df$end[df$chr == "GRC2"],
       y0 = 0, y1 = 0.1, border = "black", lwd = 0.3)

kpAddLabels(kp, labels="satellites", r0=0, r1=0.1, data.panel = 1)

# tes
kpRect(kp, chr = "GRC1", x0 = tes_grc1$start, x1 = tes_grc1$end,
       y0 = 0.2, y1 = 0.3, border = tes_grc1$col, col = tes_grc1$col, lwd = 0.3)
kpRect(kp, chr = "GRC1", x0 = df$start[df$chr == "GRC1"], x1 = df$end[df$chr == "GRC1"],
       y0 = 0.2, y1 = 0.3, border = "black", lwd = 0.3)

kpRect(kp, chr = "GRC2", x0 = tes_grc2$start, x1 = tes_grc2$end,
       y0 = 0.2, y1 = 0.3, border = tes_grc2$col, col = tes_grc2$col, lwd = 0.3)
kpRect(kp, chr = "GRC2", x0 = df$start[df$chr == "GRC2"], x1 = df$end[df$chr == "GRC2"],
       y0 = 0.2, y1 = 0.3, border = "black", lwd = 0.3)

kpAddLabels(kp, labels="TEs", r0=0.2, r1=0.3, data.panel = 1)

# genes
kpRect(kp, chr = "GRC1", x0 = genes_grc1$start, x1 = genes_grc1$end,
       y0 = 0.4, y1 = 0.5, border = "#4c5a79ff", col = "#4c5a79ff", lwd = 0.3)
kpRect(kp, chr = "GRC1", x0 = df$start[df$chr == "GRC1"], x1 = df$end[df$chr == "GRC1"],
       y0 = 0.4, y1 = 0.5, border = "black", lwd = 0.3)

kpRect(kp, chr = "GRC2", x0 = genes_grc2$start, x1 = genes_grc2$end,
       y0 = 0.4, y1 = 0.5, border = "#4c5a79ff", col = "#4c5a79ff", lwd = 0.3)
kpRect(kp, chr = "GRC2", x0 = df$start[df$chr == "GRC2"], x1 = df$end[df$chr == "GRC2"],
       y0 = 0.4, y1 = 0.5, border = "black", lwd = 0.3)
kpAddLabels(kp, labels="genes", r0=0.4, r1=0.5, data.panel = 1)

# synteny blocks
kpRect(kp, chr = "GRC1", x0 = sb$q_start, x1 = sb$q_end,
       y0 = 0.6, y1 = 0.7, border = "black", col = sb$t_col, lwd = 0.3)
kpRect(kp, chr = "GRC2", x0 = sb$t_start, x1 = sb$t_end,
       y0 = 0.6, y1 = 0.7, border = "black", col = "#f79943ff", lwd = 0.3)
kpAddLabels(kp, labels="SBs", r0=0.6, r1=0.7, data.panel = 1)

dev.off()
################################################################################

