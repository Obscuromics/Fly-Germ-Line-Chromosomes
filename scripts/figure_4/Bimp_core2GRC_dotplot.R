# Plot chromosome alignments as dot plots
################################################################################
require(dplyr)
require(ggplot2)
require(S4Vectors)
require(stringr)
################################################################################
# Functions
# read paf alignment file
load_paf_file <- function(file, t_chroms = NULL, t_ranges = NULL, t_sp = NULL,
                          q_chroms = NULL, q_ranges = NULL, q_sp = NULL){
  
  # read file
  ali <- read.table(file = file, header = FALSE)[,c(1:9)]
  colnames(ali) <- c("target", "t_length", "t_start", "t_end", "t_string",
                     "query", "q_length", "q_start", "q_end")
  
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
    select(target, t_linear_start, t_linear_end, t_species, t_chrom_size,
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
  df[which(df$origin == "Sciaridae"), "colour"] <- "#4CCEAF"
  df[which(df$origin == "Cecidomyiidae"), "colour"] <- "#CE8EDA"
  df[which(df$origin == "Other"), "colour"] <- "#CE8EDA" <- "grey90"
  
  colnames(df) <- c("busco", "status", "chrom", "start", 
                    "end", "strand", "origin", "sp", "chrom_size", "colour")
  return(df)
}
################################################################################
setwd("/Users/ab66/Documents/sanger_work/diptera/analysis_on_curated_genomes")
home <- getwd()
data <- "data"
figures <- "figures"

teal <- rgb(76, 206, 175, maxColorValue = 255)
magenta <- rgb(206, 142, 218, maxColorValue = 255)

target_species <- "Bimp"
target_chrom2plot <- c("SUPER_1", "SUPER_2", "SUPER_3", "SUPER_X")

query_species <- "Bimp"
query_chrom2plot <- c("SUPER_GRC")

chrom_files <- c("idBraCopr2.1.chrom_sizes.tsv",
                 "idBraImpa2.1.primary.chrom_sizes.tsv",
                 "idLycInge5.1.primary.chrom_sizes.tsv")

chrom_list <- file.path(home, data, "chrom_sizes", chrom_files)
names(chrom_list) <- c("Bcop", "Bimp", "Ling")

busco_files <- c("BraCopr_buscos.diptera_odb10.tsv",
                 "BraImpa_buscos.diptera_odb10.tsv",
                 "LycInge_buscos.diptera_odb10.tsv")

busco_list <- file.path(home, data, "buscos", busco_files)
names(busco_list) <- c("Bcop", "Bimp", "Ling")

alignment_file <- "alignments/Bimp_Bimp.m.1aln.paf"
################################################################################
### --- load paf alignment --- ###
ali <- load_paf_file(
  file = file.path(home, data, alignment_file),
  t_chroms = target_chrom2plot, t_sp = target_species, 
  q_chroms = query_chrom2plot, q_sp = query_species)

### --- convert coordinates to linear format --- ###
ali_lin <- prepare_linear_coordinates(
  ali = ali, t_sp = target_species, q_sp = query_species, chrom_list = chrom_list)

### --- plotting --- ###
p <- ggplot(ali_lin)  +
  aes(x = t_linear_start, y = q_linear_start, xend = t_linear_end, yend = q_linear_end)

p <- p + geom_segment(lineend = "round", linewidth = 0.4) +
  labs(x = NULL, y = NULL) +
  theme_bw()

# calculate breaks and labels positions
breaks <- list()
ticks  <- list()
labels <- list()

chr_info_t <- ali_lin %>%
  arrange(target) %>%
  distinct(target, t_chrom_size) %>%
  mutate(cum_end = cumsum(t_chrom_size))

chr_info_q <- ali_lin %>%
  arrange(query) %>%
  distinct(query, q_chrom_size) %>%
  mutate(cum_end = cumsum(q_chrom_size))

breaks$target <- c(0, chr_info_t$cum_end)
ticks$target <- calcLabelPosition(breaks$target)
labels$target <- chr_info_t$target

# add breaks and labels
p <- p +
  scale_x_continuous(expand = c(0, 0), minor_breaks = NULL,
                     breaks = breaks$target, labels = NULL, position = 'top',
                     sec.axis = dup_axis(breaks=ticks$target, labels=labels$target))

breaks$query <- c(0, chr_info_q$cum_end)
ticks$query <- calcLabelPosition(breaks$query)
labels$query <- chr_info_q$query

p <- p +
  scale_y_continuous(guide = guide_axis(angle = 90), 
                     expand = c(0, 0), minor_breaks = NULL,
                     breaks = breaks$query, labels = NULL, position = 'right',
                     sec.axis = dup_axis(breaks=ticks$query, labels=labels$query)) +
  coord_fixed()

p <- p + expand_limits(x = 0, y = 0) # force to plot from 0
################################################################################
# add chromosome coloured based on origin
origin <- read.table(
  file.path(home, data, "phylogeny/busco_grc_classification_odb10_diptera.tsv"),
  sep = "\t", header = TRUE)[,c(0:5)]
origin[c("sp", "chr")] <- str_split_fixed(origin$spchr, "_", 2)
colnames(origin) <- c("busco", "spchr", "start", "end", "origin", "sp", "chr")

buscos <- read_busco_file(
  file_name = busco_list[[query_species]], species = query_species, 
  buscos_to_origin = origin, chrom_list = chrom_list)

# transform busco coordinates of the query to linear format
q_buscos <- buscos %>%
  filter(chrom %in% query_chrom2plot) %>%
  arrange(chrom, start)

chr_offset <- 0
q_buscos_lin <- NULL

for(i in unique(q_buscos$chrom)){
  df <- q_buscos[q_buscos$chrom == i,]
  size <- unique(df$chrom_size)
  
  # calculate linear coordinates
  df$linear_start <- chr_offset + df$start
  df$linear_end <- chr_offset + df$end
  
  # bind with the main df
  q_buscos_lin <- rbind(q_buscos_lin, df)
  
  # change chr offset
  chr_offset <- chr_offset + size
}

### -- plotting buscos with origin --- ###
if(nrow(chr_info_q) == 1){
  chr_info_q$cum_start <- 0
}else{
  chr_info_q$cum_start <- c(0, chr_info_q$cum_end[1:nrow(chr_info_q)-1])
}

for(i in chr_info_q$query){
  df <- chr_info_q[chr_info_q$query == i,]
  df_buscos <- q_buscos_lin[q_buscos_lin$chrom == i,]
  
  # plot buscos
  p <- p + annotate(
    "rect", ymin = df_buscos$linear_start, ymax = df_buscos$linear_end,
    xmin = -7000000, xmax = -1000, 
    colour = df_buscos$colour, fill = df_buscos$colour, linewidth = 0.3)
  
  p <- p + annotate(
    "rect", ymin = df$cum_start, ymax = df$cum_end,
    xmin = -7000000, xmax = -1000, 
    colour = "black", fill = NA, linewidth = 0.3)
}

plt_file_name <- paste0(target_species, "2", 
                        query_species, ".dotplot_with_buscos")

ggsave(plot = p , filename = file.path(home, figures, paste0(plt_file_name, ".svg")), 
       device = "svg")#, units = "cm", width = 7, height = 4)

ggsave(plot = p , filename = file.path(home, figures, paste0(plt_file_name, ".png")), 
       device = "png")#, units = "cm", width = 100, height = 40)
