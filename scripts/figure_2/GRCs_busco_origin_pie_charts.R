# Plot synteny between GRCs using buscos and colour based on origin of the genes
################################################################################
require(stringr)
require(dplyr)
################################################################################
# Functions
# read BUSCO files
read_busco_file <- function(file_name, species, buscos_to_origin, chrom){
  df <- read.csv(file_name, sep = '\t', comment.char = '#', header = FALSE,
                 na.strings = c("", "NA"))[,c(0:6)]
  colnames(df) <- c("busco", "status", "chr", "start", "end", "strand")
  df_new <- df %>% filter(strand == "+")
  df_new <- rbind(
    df_new, df %>% filter(strand == "-") %>% rename(start = end, end = start)) %>%
    arrange(chr, start)
  df_new <- df_new %>% filter(chr %in% chrom)
  buscos_to_origin <- buscos_to_origin %>% 
    filter(sp == species) %>% select(busco, chr, start, end, origin)
  df_new <- left_join(df_new, buscos_to_origin)
  return(df_new)
}
################################################################################
home <- getwd()
data <- "data"
figures <- "figures"

busco_files <- c("LycInge_buscos.diptera_odb10.tsv",
                 "BraCopr_buscos.diptera_odb10.tsv",
                 "BraImpa_buscos.diptera_odb10.tsv")
################################################################################
# read classification file
origin <- read.table(
  file.path(home, data, "phylogeny/busco_grc_classification_odb10_diptera.tsv"),
  sep = "\t", header = TRUE)[,c(0:5)]

origin[c("sp", "chr")] <- str_split_fixed(origin$spchr, "_", 2)
colnames(origin) <- c("busco", "spchr", "start", "end", "origin", "sp", "chr")

# read busco files
busco_list <- file.path(home, data, "buscos", busco_files)
names(busco_list) <- c("Ling", "Bcop", "Bimp")
stats <- list()

for(file in busco_list){
  i <- match(file, busco_list)
  sp <- names(busco_list)[i]
  df <- read_busco_file(
    file_name = file, species = sp, buscos_to_origin = origin, 
    chrom = c("SUPER_GRC1", "SUPER_GRC2", "SUPER_GRC"))
  df[which(is.na(df$origin)), "origin"] <- "Other"
  #duplicated_buscos <- df$busco[duplicated(df$busco)]
  #df <- df %>% filter(!busco %in% duplicated_buscos)
  agg <- aggregate(df$busco, by = list(df$origin), FUN = length)
  colnames(agg) <- c("origin", "count")
  agg$colour <- c("#CE8EDA", "grey90", "#4CCEAF")
  stats[[sp]] <- agg
}

### --- plotting --- ###
for(sp in names(stats)){
  pdf(file = file.path(home, figures, paste0(sp, "_origin_pie_chart.pdf")))
  pie(stats[[sp]]$count, labels = stats[[sp]]$origin, col = stats[[sp]]$colour,
      main = sp)
  dev.off()
}