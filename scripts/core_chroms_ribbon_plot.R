# Plot synteny between GRCs using buscos and colour based on origin of the genes
################################################################################
require(stringr)
require(dplyr)
################################################################################
# Functions
# read BUSCO files
read_busco_file <- function(file_name, prefix, species, chrom){
  chr_label <- paste0('chr', prefix)
  df <- read.csv(file_name, sep = '\t', comment.char = '#', header = FALSE,
                 na.strings = c("", "NA"))[,c(0:6)]
  colnames(df) <- c("busco", "status", "chr", "start", "end", "strand")
  
  # keep only complete genes
  df <- df[df$status %in% c("Complete"),]
  
  # swap start and end for buscos on "-" strand
  df_new <- df %>% filter(strand == "+")
  df_new <- rbind(
    df_new, df %>% dplyr::filter(strand == "-") %>% 
      dplyr::rename(start = end, end = start)) %>%
    arrange(chr, start)
  df <- df_new
  
  # filter for duplicated genes
  #duplicated_buscos <- df$busco[duplicated(df$busco)]
  #df[which(df$busco %in% duplicated_buscos), "origin"] <- NA
  
  colnames(df) <- c('busco', 'status', chr_label, paste0(prefix, 'start'),
                    paste0(prefix, 'end'), 'strand')
  
  return(df)
}
################################################################################
home <- getwd()
data <- "data"
figures <- "figures"

busco_files <- c("LycInge_buscos.diptera_odb12.tsv",
                 "BraCopr_buscos.diptera_odb12.tsv",
                 "BraImpa_buscos.diptera_odb12.tsv")

chrom_files <- c("LycInge.core.chrom_sizes.tsv",
                 "BraCopr.core.chrom_sizes.tsv",
                 "BraImpa.core.chrom_sizes.tsv")
################################################################################
# read busco and chromosome files and make alignments
busco_list <- file.path(home, data, "buscos", busco_files)
names(busco_list) <- c("Ling", "Bcop", "Bimp")

chrom_list <- file.path(home, data, "buscos", chrom_files)
names(chrom_list) <- c("Ling", "Bcop", "Bimp")

minimum_buscos = 1

# load synteny_plotter functions
#devtools::source_url("https://github.com/Obscuromics/synteny_plotter/blob/dev/scripts/helper_functions.R?raw=TRUE")

# initiate reference
ref_chroms <- read.table(chrom_list[1], sep = '\t', header = TRUE)
ref_chroms <- ref_chroms %>% arrange(order)
ref_df <- read_busco_file(
  file_name = busco_list[1], prefix = "R", species = "Ling",
  chrom = ref_chroms)

# add colours
col_list <- c("#4cceaf","#03a487", "#007c61", "#00563e")
ref_chroms$colour <- col_list
busco2colour <- ref_chroms[,c(1,6)]
temp <- ref_df[,c(1,3)]
colnames(temp) <- c('busco', 'chr')
busco2colour <- merge(busco2colour, temp, by = 'chr')

chr_offset <- max(ref_chroms$length) / 2
processed_Q_list <- list()
max_ends <- list()

# save ref as temp_ref
temp_ref_chroms <- ref_chroms
temp_ref_df <- ref_df

source('/Users/ab66/Documents/sanger_work/Tools/synteny_plotter/scripts/helper_functions.R')

for (file in busco_list[-1]){
  i <- match(file, busco_list)
  sp <- names(busco_list)[i]
  query_chroms <- read.table(chrom_list[i], sep = '\t', header = TRUE)
  query_df <- read_busco_file(
    file = file, prefix = "Q", species = sp,
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
    'Rend', 'strand')
  temp_ref_chroms <- query_chroms
}
################################################################################
max_end <- max(unlist(max_ends))
plot_length <- max_end
gap <- 5
alpha = 0.6
show_outline = TRUE

pdf(file.path(home, figures, "core_chroms_ribbon_plot.pdf"))
print('[+] Generating plot')
plot(0,cex = 0, xlim = c(1, plot_length), 
     ylim = c(((gap+1)*-1*2*2*2),((gap+1)*2*2*2)),
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
    y1 <- gap-y_offset-y_increment
    y2 <- gap-y_offset
    
    plot_one_ref_chr(alignments, adjustment_length_R, adjustment_length_Q,
                     y1, y2, busco2colour, lwd = 0.5)
    
    ### --- plotting reference chromosomes --- ###
    counter <- 0
    offset <- 0
    
    ref_buscos <- read_busco_file(
      file_name = busco_list[j-1], prefix = 'R', species = names(busco_list[j-1]), 
      chrom = chr_order_R)
    
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
      #rect(Rfirst+adjustment_length_R, gap-y_offset,
       #    Rlast+adjustment_length_R, gap-y_offset+2, col = "white", lwd = 0.5)
      
      segments(Rfirst+adjustment_length_R, gap-y_offset, 
               Rlast+adjustment_length_R, gap-y_offset, lwd = 10)
      
      text(x = ((Rlast+Rfirst+1)/2)+adjustment_length_R, y = gap-y_offset, 
           label = ref_chroms[ref_chroms$chr == i,]$annot,
           srt = 0, cex = 0.5, col = "grey")
    }
  }
  
  if(main_counter == (length(processed_Q_list) - 4)){
    
    adjustment_length_R <- adjustment_length_Q
    adjustment_length_Q <- (max_end - max(alignments$Qend)) / 2 
    
    ### --- plot alignments --- ###
    y1 <- gap-y_offset-y_increment
    y2 <- gap-y_offset
    
    plot_one_ref_chr(alignments, adjustment_length_R, adjustment_length_Q,
                     y1, y2, busco2colour, lwd = 0.5)
    
    ### --- plotting reference chromosomes --- ###
    counter <- 0
    offset <- 0
    
    ref_buscos <- read_busco_file(
      file_name = busco_list[j-1], prefix = 'R', species = names(busco_list[j-1]), 
      chrom = chr_order_R)
    
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
      
      # plot chromosomes
      #rect(Rfirst+adjustment_length_R, gap-y_offset,
       #    Rlast+adjustment_length_R, gap-y_offset+2, col = "white", lwd = 0.5)
      
      segments(Rfirst+adjustment_length_R, gap-y_offset, 
               Rlast+adjustment_length_R, gap-y_offset, lwd = 10)
      
      text(x = ((Rlast+Rfirst+1)/2)+adjustment_length_R, y = gap-y_offset, 
           label = ref_chroms[ref_chroms$chr == i,]$annot,
           srt = 0, cex = 0.5, col = "grey")
    }
    
    ### --- plotting query chromosomes --- ###
    counter <- 0
    offset <- 0
    
    query_buscos <- read_busco_file(
      file_name = busco_list[j], prefix = 'Q', species = names(busco_list[j]), 
      chrom = chr_order_Q)
    
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
      #rect(Qfirst+adjustment_length_Q, 1-gap-y_offset-2,
      #     Qlast+adjustment_length_Q, 1-gap-y_offset-2+2, col = "white", lwd = 0.5)
      
      segments(Qfirst+adjustment_length_Q, gap-y_offset-y_increment, 
               Qlast+adjustment_length_Q, gap-y_offset-y_increment, lwd = 10)
      
      text(x = ((Qlast+Qfirst+1)/2)+adjustment_length_Q, y = 1-gap-y_offset-2, 
           label = query_chroms[query_chroms$chr == i,]$annot,
           srt = 0, cex = 0.5, col = "grey")
    }
  }
  
  main_counter <- main_counter + 5
  y_offset <- y_offset + y_increment
  ref_chroms <- query_chroms
}

dev.off()
################################################################################