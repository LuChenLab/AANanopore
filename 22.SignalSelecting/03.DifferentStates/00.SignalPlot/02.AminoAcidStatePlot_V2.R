setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(shiny)
library(plotly)
library(ggpubr)
get_density <- function(x, y, n = 1000, ...) {
  dens <- MASS::kde2d(x, y, n = n, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii] / max(dens$z[ii]))
}

AABlockade <- data.table(AA = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"), 
                         Blockade = c(0.14718, 0.17148, 0.16516, 0.21409, 0.18622, 0.24528, 0.1882, 0.12072, 0.24652, 0.20722, 0.1995, 0.16875, 0.19772, 0.22018, 0.2183, 0.13131, 0.16101, 0.22744, 0.21276, 0.19044))
AABlockade$amino_acid <- plyr::mapvalues(AABlockade$AA, names(Biostrings::AMINO_ACID_CODE), Biostrings::AMINO_ACID_CODE)


meta1 <- data.table(openxlsx::read.xlsx("./data/单氨基酸实验信息表20210517.xlsx", sheet = 1, cols = 1:7))
meta2 <- data.table(openxlsx::read.xlsx("./data/单氨基酸实验信息表20210517.xlsx", sheet = 2, cols = 1:7))
meta3 <- data.table(openxlsx::read.xlsx("./data/单氨基酸实验信息表20210517.xlsx", sheet = 3, cols = 1:7))
meta4 <- data.table(openxlsx::read.xlsx("./data/单氨基酸实验信息表20210517.xlsx", sheet = 4, cols = 1:7))
meta <- rbind(meta1, meta2, meta3, meta4, use.names = FALSE)
colnames(meta) <- c("file_name", "date", "amino_acid", "concentration", "start_time", "end_time", "type")
meta[amino_acid == "cys", amino_acid := "Cys"]
meta <- meta[concentration == 0 | amino_acid == "blank"]
meta$file_path <- mapply(meta$file_name, FUN = function(x) list.files("./data", recursive = TRUE, full.names = TRUE, pattern = as.character(x))[1])
meta <- meta[!is.na(file_path)]
meta[, start_time := as.numeric(start_time)]
meta[, end_time := as.numeric(end_time)]
# meta <- meta[!grepl("abf", file_name)]
meta[, file_name := gsub(".abf$", "", file_name)]
meta <- meta[file_name %in% gsub("RawSignal_", "", gsub(".MainL0.txt", "", list.files("./analysis/22.SignalSelecting/03.DifferentStates/Background/")))]

for (aat in AABlockade$amino_acid) {
  bgsignal <- unique(paste0("./analysis/22.SignalSelecting/03.DifferentStates/Background/RawSignal_", meta[amino_acid == aat, file_name], ".MainL0.txt"))
  bgsignal <- do.call(rbind, lapply(bgsignal, fread))
  bgsignal$D <- bgsignal[, get_density(x = Blockade, y = log10(DwellTime))]
  
  dataset <- do.call(rbind, lapply(list.files(paste0("./analysis/22.SignalSelecting/03.DifferentStates/", aat), "MainL0.txt", full.names = T), fread))
  dataset$File <- paste0("File", stringr::str_remove_all(dataset$ID, "_([[:digit:]]+)$"))
  State1 <- tryCatch(fread(paste0("./analysis/22.SignalSelecting/03.DifferentStates/", aat, "/", aat, "_State1.txt")), error = function(e) NULL)
  State2 <- tryCatch(fread(paste0("./analysis/22.SignalSelecting/03.DifferentStates/", aat, "/", aat, "_State2.txt")), error = function(e) NULL)
  dataset[ID %in% State1$ID, Group := "State1"]
  dataset[ID %in% State2$ID, Group := "State2"]
  dataset[is.na(Group), Group := "Noise"]
  dataset$D <- dataset[, get_density(x = Blockade, y = log10(DwellTime))]
  dataset[, Group := factor(Group, levels = c("Noise", "State1", "State2"))]
  setkey(dataset, Group)
  data1 <- dataset[DwellTime <= max(c(State1$DwellTime, State2$DwellTime)) & Blockade < 0.6]
  
  ggplot(data = bgsignal, mapping = aes(x = Blockade, y = DwellTime)) + 
    geom_point(aes(colour = D), size = 1) + 
    scale_y_log10(limits = data1[, range(DwellTime)]) + 
    scale_x_continuous(limits = data1[, range(Blockade)]) + 
    scale_colour_gradient2(high = "red", low = "white") + 
    theme_light(base_size = 15) + 
    theme(panel.grid = element_blank(), 
          axis.title.x = element_blank(), 
          axis.text.x = element_blank(),
          legend.position = "none") -> p1
  
  ggplot(data = data1, mapping = aes(x = Blockade, y = DwellTime)) + 
    geom_point(aes(colour = D), size = 1) + 
    scale_y_log10() + 
    scale_colour_gradient2(high = "red", low = "white") + 
    theme_light(base_size = 15) + 
    theme(panel.grid = element_blank(), 
          axis.title.x = element_blank(), 
          axis.text.x = element_blank(),
          legend.position = "none") -> p2
  
  ggplot(data = data1, mapping = aes(x = Blockade, y = DwellTime)) +
    geom_density_2d_filled(aes(fill = ..level..), contour_var = "ndensity", breaks = seq(0, 1, length.out = 100)) + 
    scale_y_log10() + 
    scale_fill_manual(values = c(grey(99:0/99), "red", "red"), aesthetics = c("fill", "color")) + 
    geom_point(aes(colour = Group, alpha = D), size = 1) + 
    theme_light(base_size = 15) + 
    theme(panel.grid = element_blank(), legend.position = "none") -> p3
  p <- cowplot::plot_grid(p1, p2, p3, ncol = 1)
  ggsave(paste0("analysis/22.SignalSelecting/03.DifferentStates/00.SignalPlot/", aat, "_V2.pdf"), p, width = 6, height = 9)
}




aat <- "CbC"
bgsignal <- list.files('./analysis/22.SignalSelecting/03.DifferentStates/CbC', "Background", full.names = T)
bgsignal <- do.call(rbind, lapply(bgsignal, fread))
bgsignal$D <- bgsignal[, get_density(x = Blockade, y = log10(DwellTime))]

dataset <- do.call(rbind, lapply(list.files('./analysis/22.SignalSelecting/03.DifferentStates/CbC', "CbC_RawSignal", full.names = T), fread))
dataset$File <- paste0("File", stringr::str_remove_all(dataset$ID, "_([[:digit:]]+)$"))
State1 <- tryCatch(fread(paste0("./analysis/22.SignalSelecting/03.DifferentStates/", aat, "/", aat, "_State1.txt")), error = function(e) NULL)
State2 <- tryCatch(fread(paste0("./analysis/22.SignalSelecting/03.DifferentStates/", aat, "/", aat, "_State2.txt")), error = function(e) NULL)
dataset[ID %in% State1$ID, Group := "State1"]
dataset[ID %in% State2$ID, Group := "State2"]
dataset[is.na(Group), Group := "Noise"]
dataset$D <- dataset[, get_density(x = Blockade, y = log10(DwellTime))]
dataset[, Group := factor(Group, levels = c("Noise", "State1", "State2"))]
setkey(dataset, Group)
data1 <- dataset[DwellTime <= max(c(State1$DwellTime, State2$DwellTime)) & Blockade < 0.6]

ggplot(data = bgsignal, mapping = aes(x = Blockade, y = DwellTime)) + 
  geom_point(aes(colour = D), size = 1) + 
  scale_y_log10(limits = data1[, range(DwellTime)]) + 
  scale_x_continuous(limits = data1[, range(Blockade)]) + 
  scale_colour_gradient2(high = "red", low = "white") + 
  theme_light(base_size = 15) + 
  theme(panel.grid = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        legend.position = "none") -> p1

ggplot(data = data1, mapping = aes(x = Blockade, y = DwellTime)) + 
  geom_point(aes(colour = D), size = 1) + 
  scale_y_log10() + 
  scale_colour_gradient2(high = "red", low = "white") + 
  theme_light(base_size = 15) + 
  theme(panel.grid = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        legend.position = "none") -> p2

ggplot(data = data1, mapping = aes(x = Blockade, y = DwellTime)) +
  geom_density_2d_filled(aes(fill = ..level..), contour_var = "ndensity", breaks = seq(0, 1, length.out = 100)) + 
  scale_y_log10() + 
  scale_fill_manual(values = c(grey(99:0/99), "red", "red"), aesthetics = c("fill", "color")) + 
  geom_point(aes(colour = Group, alpha = D), size = 1) + 
  theme_light(base_size = 15) + 
  theme(panel.grid = element_blank(), legend.position = "none") -> p3

p <- cowplot::plot_grid(p1, p2, p3, ncol = 1)
ggsave(paste0("analysis/22.SignalSelecting/03.DifferentStates/00.SignalPlot/", aat, "_V2.pdf"), p, width = 6, height = 9)




