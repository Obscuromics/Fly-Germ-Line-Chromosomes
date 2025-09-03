# Plot synteny between GRCs using buscos and colour based on origin of the genes
################################################################################
require(stringr)
require(dplyr)
################################################################################
# Functions
# read BUSCO files
read_busco_file <- function(file_name, prefix, species, buscos_to_origin, chrom){
  chr_label <- paste0('chr', prefix)
  df <- read.csv(file_name, sep = '\t', comment.char = '#', header = FALSE,
                 na.strings = c("", "NA"))[,c(0:6)]
  colnames(df) <- c("busco", "status", "chr", "start", "end", "strand")
  df <- df %>% filter(chr %in% chrom$chr)
  buscos_to_origin <- buscos_to_origin %>% 
    filter(sp == species) %>% select(busco, chr, start, end, origin)#, colour)
  df <- left_join(df, buscos_to_origin)
  #df[which(is.na(df$colour)), "colour"] <- "grey90"
  #df <- df %>% na.omit() %>% arrange(chr) 
  colnames(df) <- c('busco', 'status', chr_label, paste0(prefix, 'start'),
                    paste0(prefix, 'end'), 'strand', 'origin')#, 'colour')
  return(df)
}
################################################################################
home <- getwd()
data <- "data"
figures <- "figures"

busco_files <- c("LycInge_buscos.diptera_odb10.tsv",
                 "BraCopr_buscos.diptera_odb10.tsv",
                 "BraImpa_buscos.diptera_odb10.tsv")

chrom_files <- c("LycInge.chrom_sizes.tsv",
                 "BraCopr.chrom_sizes.tsv",
                 "BraImpa.chrom_sizes.tsv")
################################################################################
# read classification file
origin <- read.table(
  file.path(home, data, "phylogeny/busco_grc_classification_odb10_diptera.tsv"),
  sep = "\t", header = TRUE)[,c(0:5)]

origin[c("sp", "chr")] <- str_split_fixed(origin$spchr, "_", 2)
#origin[which(origin$origin == "Sciaridae"), "col"] <- "#4CCEAF"
#origin[which(origin$origin == "Cecidomyiidae"), "col"] <- "#CE8EDA"
colnames(origin) <- c("busco", "spchr", "start", "end", "origin",
                      "sp", "chr")#, "colour")

# remove busco genes that have different origin
#origin_flt <- NULL
#diff_origin <- NULL

#for(busco in unique(origin$busco)){
#  df <- origin[origin$busco == busco,]
#  if(length(unique(df$origin)) == 1){
#    origin_flt <- rbind(origin_flt, df)
#  }else{
#    diff_origin <- rbind(diff_origin, df)
#  }
#}

# read busco and chromosome files and make alignments
busco_list <- file.path(home, data, "buscos", busco_files)
names(busco_list) <- c("Ling", "Bcop", "Bimp")

chrom_list <- file.path(home, data, "buscos", chrom_files)
names(chrom_list) <- c("Ling", "Bcop", "Bimp")

minimum_buscos = 1

# load synteny_plotter functions
devtools::source_url("https://github.com/Obscuromics/synteny_plotter/blob/dev/scripts/helper_functions.R?raw=TRUE")

# initiate reference
ref_chroms <- read.table(chrom_list[1], sep = '\t', header = TRUE)
ref_chroms <- ref_chroms %>% arrange(order)
ref_df <- read_busco_file(
  file_name = busco_list[1], prefix = "R", species = "Ling", 
  buscos_to_origin = origin,
  chrom = ref_chroms)

chr_offset <- max(ref_chroms$length) / 2
processed_Q_list <- list()
max_ends <- list()

# save ref as temp_ref
temp_ref_chroms <- ref_chroms
temp_ref_df <- ref_df


for (file in busco_list[-1]){
  i <- match(file, busco_list)
  sp <- names(busco_list)[i]
  query_chroms <- read.table(chrom_list[i], sep = '\t', header = TRUE)
  query_df <- read_busco_file(
    file = file, prefix = "Q", species = sp, buscos_to_origin = origin,
    chrom = query_chroms)
  
  processed_Q <- make_alignment_table(
    temp_ref_df, temp_ref_chroms, query_df, query_chroms, chr_offset)
  
  alignments <- processed_Q$alignments
  processed_Q_list <- append(processed_Q_list, processed_Q)
  max_ends <- append(max_ends, max(alignments$Rend))
  max_ends <- append(max_ends, max(alignments$Qend))
  temp_ref_df <- query_df
  colnames(temp_ref_df) <- c(
    'busco', 'status', 'chrR', 'Rstart',
    'Rend', 'strand', 'origin')#, 'colour')
  temp_ref_chroms <- query_chroms
}
################################################################################
max_end <- max(unlist(max_ends))
plot_length <- max_end
gap <- 5
alpha = 0.6
show_outline = TRUE

pdf(file.path(home, figures, "GRCs_ribbon_plot.pdf"))
print('[+] Generating plot')
plot(0,cex = 0, xlim = c(1, plot_length), 
     #ylim = c(((gap+1)*-1*length(busco_list)*2),((gap+1)*length(busco_list)*2)),
     ylim = c(((gap+1)*-1*2*2*2),((gap+1)*2*2*2)),
     #ylim = c(-40, 40),
     xlab = "", ylab = "", bty = "n", yaxt="n", xaxt="n")

main_counter <- 1
y_offset <- 0
y_increment <- 11

for (file in busco_list[-1]){
  print(file)
  j <- match(file, busco_list)
  query_chroms <- read.table(chrom_list[j], sep = '\t', header = TRUE)
  
  alignments <- processed_Q_list[[main_counter]]
  chr_order_R <- processed_Q_list[[main_counter+1]]
  chr_order_Q <- processed_Q_list[[main_counter+2]]
  offset_list_R <- processed_Q_list[[main_counter+3]]
  offset_list_Q <- processed_Q_list[[main_counter+4]]
  
  if(main_counter != (length(processed_Q_list) - 4)){ # plot only reference and alignments
    
    if (max(alignments$Qend) != max_end){ # i.e. if this is the longest chr_set
      if (max(alignments$Rend) != max_end){
        adjustment_length_R <- (max_end - max(alignments$Rend)) / 2 
        adjustment_length_Q <- (max_end - max(alignments$Qend)) / 2 
      }
      else{
        adjustment_length_R <- 0
        adjustment_length_Q <- (max_end - max(alignments$Qend)) / 2 
      }
    }
    else{
      adjustment_length_Q <- 0
      adjustment_length_R <- (max_end - max(alignments$Rend)) / 2 
    }
    
    ### --- plot alignments --- ###
    alignments[which(alignments$origin.x == alignments$origin.y &  
                       alignments$origin.x == "Sciaridae"), "colour"] <- "#4CCEAF"
    
    alignments[which(alignments$origin.x == alignments$origin.y &  
                       alignments$origin.x == "Cecidomyiidae"), "colour"] <- "#CE8EDA"
    
    busco_to_origin <- alignments[, c("busco", "colour")]
    busco_to_origin[which(is.na(busco_to_origin$colour)), "colour"] <- "grey90"
    #busco_to_origin$colour <- adjustcolor(busco_to_origin$colour, alpha.f = 0.8)
    
    # plot lines with different/unknown origin
    aln_grey <- alignments[is.na(alignments$colour),]
    plot_one_ref_chr(aln_grey, adjustment_length_R, adjustment_length_Q, y_offset, 
                     busco_to_origin, lwd = 0.2)
    
    # plot lines with the same origin
    aln_col <- alignments[!is.na(alignments$colour),]
    plot_one_ref_chr(aln_col, adjustment_length_R, adjustment_length_Q, y_offset, 
                     busco_to_origin, lwd = 0.4)
    
    ### --- plotting reference chromosomes --- ###
    counter <- 0
    offset <- 0
    
    ref_buscos <- read_busco_file(
      file_name = busco_list[j-1], prefix = 'R', species = names(busco_list[j-1]), 
      buscos_to_origin = origin, chrom = chr_order_R)
    ref_buscos[which(ref_buscos$origin == "Sciaridae"), "colour"] <- "#4CCEAF"
    ref_buscos[which(ref_buscos$origin == "Cecidomyiidae"), "colour"] <- "#CE8EDA"
    
    for (i in chr_order_R$chr){
      chr_length <- chr_order_R[chr_order_R$chr == i,]$length
      chr_buscos <- ref_buscos[ref_buscos$chrR == i,]
      
      Rfirst <- offset
      Rlast <- chr_order_R[chr_order_R$chr == i,]$length + offset
      
      if(counter != 0){
        Rfirst <- offset  # allows for accumulative chr positions
        Rlast <- chr_length + offset # allows for accumulative chr positions
        
        chr_buscos$Rstart <- chr_buscos$Rstart + offset
        chr_buscos$Rend <- chr_buscos$Rend + offset
      }
      
      offset <- offset + chr_length + chr_offset # accumulative offset
      counter <- counter + 1
      
      # plot chromosome outline
      rect(Rfirst+adjustment_length_R, gap-y_offset,
           Rlast+adjustment_length_R, gap-y_offset+2, col = "white", lwd = 0.5)
      
      # plot buscos with NA origin first
      chr_buscos_NA <- chr_buscos[is.na(chr_buscos$origin),]
      rect(chr_buscos_NA$Rstart+adjustment_length_R, gap-y_offset,
           chr_buscos_NA$Rend+adjustment_length_R, gap-y_offset+2, 
           col = "grey90", border = "grey90", lwd = 0.2)
      
      # plot buscos with known origin
      chr_buscos_org <- chr_buscos[!is.na(chr_buscos$origin),]
      rect(chr_buscos_org$Rstart+adjustment_length_R, gap-y_offset,
           chr_buscos_org$Rend+adjustment_length_R, gap-y_offset+2, 
           col = chr_buscos_org$colour, border = chr_buscos_org$colour, lwd = 0.4)
      
      text(x = ((Rlast+Rfirst+1)/2)+adjustment_length_R, y = gap-y_offset, 
           label = ref_chroms[ref_chroms$chr == i,]$annot,
           srt = 0, cex = 0.5, col = "black")
    }
  }
  
  if(main_counter == (length(processed_Q_list) - 4)){
    
    adjustment_length_R <- adjustment_length_Q
    adjustment_length_Q <- (max_end - max(alignments$Qend)) / 2 
    
    ### --- plot alignments --- ###
    alignments[which(alignments$origin.x == alignments$origin.y &  
                       alignments$origin.x == "Sciaridae"), "colour"] <- "#4CCEAF"
    
    alignments[which(alignments$origin.x == alignments$origin.y &  
                       alignments$origin.x == "Cecidomyiidae"), "colour"] <- "#CE8EDA"
    
    busco_to_origin <- alignments[, c("busco", "colour")]
    busco_to_origin[which(is.na(busco_to_origin$colour)), "colour"] <- "grey90"
    #busco_to_origin$colour <- adjustcolor(busco_to_origin$colour, alpha.f = 0.8)
    
    # plot lines with different/unknown origin
    aln_grey <- alignments[is.na(alignments$colour),]
    plot_one_ref_chr(aln_grey, adjustment_length_R, adjustment_length_Q, y_offset, 
                     busco_to_origin, lwd = 0.2)
    
    # plot lines with the same origin
    aln_col <- alignments[!is.na(alignments$colour),]
    plot_one_ref_chr(aln_col, adjustment_length_R, adjustment_length_Q, y_offset, 
                     busco_to_origin, lwd = 0.4)
    
    ### --- plotting reference chromosomes --- ###
    counter <- 0
    offset <- 0
    
    ref_buscos <- read_busco_file(
      file_name = busco_list[j-1], prefix = 'R', species = names(busco_list[j-1]), 
      buscos_to_origin = origin, chrom = chr_order_R)
    ref_buscos[which(ref_buscos$origin == "Sciaridae"), "colour"] <- "#4CCEAF"
    ref_buscos[which(ref_buscos$origin == "Cecidomyiidae"), "colour"] <- "#CE8EDA"
    
    for (i in chr_order_R$chr){
      chr_length <- chr_order_R[chr_order_R$chr == i,]$length
      chr_buscos <- ref_buscos[ref_buscos$chrR == i,]
      
      Rfirst <- offset
      Rlast <- chr_order_R[chr_order_R$chr == i,]$length + offset
      
      if(counter != 0){
        Rfirst <- offset  # allows for accumulative chr positions
        Rlast <- chr_length + offset # allows for accumulative chr positions
        
        chr_buscos$Rstart <- chr_buscos$Rstart + offset
        chr_buscos$Rend <- chr_buscos$Rend + offset
      }
      
      offset <- offset + chr_length + chr_offset # accumulative offset
      counter <- counter + 1
      
      # plot chromosome outline
      rect(Rfirst+adjustment_length_R, gap-y_offset,
           Rlast+adjustment_length_R, gap-y_offset+2, col = "white", lwd = 0.5)
      
      # plot buscos with NA origin first
      chr_buscos_NA <- chr_buscos[is.na(chr_buscos$origin),]
      rect(chr_buscos_NA$Rstart+adjustment_length_R, gap-y_offset,
           chr_buscos_NA$Rend+adjustment_length_R, gap-y_offset+2, 
           col = "grey90", border = "grey90", lwd = 0.2)
      
      # plot buscos with known origin
      chr_buscos_org <- chr_buscos[!is.na(chr_buscos$origin),]
      rect(chr_buscos_org$Rstart+adjustment_length_R, gap-y_offset,
           chr_buscos_org$Rend+adjustment_length_R, gap-y_offset+2,
           col = chr_buscos_org$colour, border = chr_buscos_org$colour, lwd = 0.4)
      
      text(x = ((Rlast+Rfirst+1)/2)+adjustment_length_R, y = gap-y_offset, 
           label = ref_chroms[ref_chroms$chr == i,]$annot,
           srt = 0, cex = 0.5, col = "black")
    }
    
    ### --- plotting query chromosomes --- ###
    counter <- 0
    offset <- 0
    
    query_buscos <- read_busco_file(
      file_name = busco_list[j], prefix = 'Q', species = names(busco_list[j]), 
      buscos_to_origin = origin, chrom = chr_order_Q)
    query_buscos[which(query_buscos$origin == "Sciaridae"), "colour"] <- "#4CCEAF"
    query_buscos[which(query_buscos$origin == "Cecidomyiidae"), "colour"] <- "#CE8EDA"
    
    for (i in chr_order_Q$chr){
      chr_length <- chr_order_Q[chr_order_Q$chr == i,]$length
      chr_buscos <- query_buscos[query_buscos$chrQ == i,]
      
      Qfirst <- offset
      Qlast <- chr_order_Q[chr_order_Q$chr == i,]$length + offset
      
      if (counter != 0){ # only need to offset start/end if this is not the first chr
        Qfirst <- offset
        Qlast <- chr_length + offset
        
        chr_buscos$Qstart <- chr_buscos$Qstart + offset
        chr_buscos$Qend <- chr_buscos$Qend + offset
      }
      
      offset <- offset + chr_length + chr_offset # accumulative offset
      counter <- counter + 1
      
      # draw outline of the chromosome
      rect(Qfirst+adjustment_length_Q, 1-gap-y_offset-2,
           Qlast+adjustment_length_Q, 1-gap-y_offset-2+2, col = "white", lwd = 0.5)
      
      # plot buscos with NA origin first
      chr_buscos_NA <- chr_buscos[is.na(chr_buscos$origin),]
      rect(chr_buscos_NA$Qstart+adjustment_length_Q, 1-gap-y_offset-2,
           chr_buscos_NA$Qend+adjustment_length_Q, 1-gap-y_offset-2+2, 
           col = "grey90", border = "grey90", lwd = 0.2)
      
      # plot buscos with known origin
      chr_buscos_org <- chr_buscos[!is.na(chr_buscos$origin),]
      rect(chr_buscos_org$Qstart+adjustment_length_Q, 1-gap-y_offset-2,
           chr_buscos_org$Qend+adjustment_length_Q, 1-gap-y_offset-2+2, 
           col = chr_buscos_org$colour, border = chr_buscos_org$colour, lwd = 0.4)
      
      text(x = ((Qlast+Qfirst+1)/2)+adjustment_length_Q, y = 1-gap-y_offset, 
           label = query_chroms[query_chroms$chr == i,]$annot,
           srt = 0, cex = 0.5, col = "black")
    }
  }
  
  main_counter <- main_counter + 5
  y_offset <- y_offset + y_increment
  ref_chroms <- query_chroms
}

dev.off()
################################################################################