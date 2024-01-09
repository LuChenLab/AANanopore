setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(ggplot2)
get_density <- function(x, y, n = 1000, ...) {
  dens <- MASS::kde2d(x, y, n = n, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii] / max(dens$z[ii]))
}

set.seed(19)
AA_Cols <- sample(c(RColorBrewer::brewer.pal(n = 12, "Paired"), RColorBrewer::brewer.pal(n = 8, "Dark2")), 20)
names(AA_Cols) <- Biostrings::AA_STANDARD
AA_Cols[AA_Cols == "#FFFF99"] <- "#1A1A1A"
names(AA_Cols) <- AMINO_ACID_CODE[Biostrings::AA_STANDARD]
names(AA_Cols)[names(AA_Cols) == "Cys"] <- "CbC"

files1 <- list.files("./analysis/22.SignalSelecting/03.DifferentStates", "_State1.txt", recursive = T, full.names = T)
files2 <- list.files("./analysis/22.SignalSelecting/03.DifferentStates", "_State2.txt", recursive = T, full.names = T)

State1 <- do.call(rbind, lapply(files1, fread))
State1[, State := "State1"]

State2 <- do.call(rbind, lapply(files2, fread))
State2[, State := "State2"]




  

for (aat in State1[, unique(A)][6:20]) {
  file_name <- State1[A == aat, gsub("File", "", unique(File))]
  if(aat == "CbC") {
    SignalCurrent <- mapply(function(x) list.files("./analysis/21.ABFProcessing/03.MixtureAA/SignalCurrent", x, recursive = T, full.names = T), file_name)
  } else {
    SignalCurrent <- mapply(function(x) list.files("./analysis/21.ABFProcessing/01.StandardAA/SignalCurrent", x, recursive = T, full.names = T), file_name)
  }
  SignalCurrent <- lapply(SignalCurrent, readRDS)
  SignalCurrent <- do.call(rbind, SignalCurrent)
  SignalCurrent <- SignalCurrent[, .(Current, x = seq_len(nrow(.SD))), ID]
  
  SignalCurrent_State1 <- SignalCurrent[ID %in% State1[A == aat, ID]]
  SignalCurrent_State1$x <- SignalCurrent_State1$x/max(SignalCurrent_State1$x)*State1[A == aat, max(DwellTime)]
  ggplot(SignalCurrent_State1, aes(x, y = abs(Current - 1), colour = ID)) + 
    geom_line(alpha = .5) + 
    scale_colour_manual(values = rep("black", SignalCurrent_State1[, length(unique(ID))])) + 
    scale_y_reverse() + 
    guides(colour = "none") + 
    theme_light(base_size = 15) + 
    labs(y = "Blockade", x = "Dwell time (ms)", title = aat) -> p
  ggsave(paste0("./analysis/22.SignalSelecting/03.DifferentStates/00.SignalPlot/", aat, "_SignalEvent.pdf"), p, width = 4, height = 4)
}



fill_cols <- colorRampPalette(c(rep("white", 1), RColorBrewer::brewer.pal(n = 9, name = "Reds")[4:9]))(100)

for (aat in State1[, unique(A)]) {
  file_name <- State1[A == aat, gsub("File", "", unique(File))]
  if(aat == "CbC") {
    SignalCurrent <- mapply(function(x) list.files("./analysis/21.ABFProcessing/03.MixtureAA/SignalCurrent", x, recursive = T, full.names = T), file_name)
  } else {
    SignalCurrent <- mapply(function(x) list.files("./analysis/21.ABFProcessing/01.StandardAA/SignalCurrent", x, recursive = T, full.names = T), file_name)
  }
  SignalCurrent <- lapply(SignalCurrent, readRDS)
  SignalCurrent <- do.call(rbind, SignalCurrent)
  SignalCurrent <- SignalCurrent[, .(Current, x = seq_len(nrow(.SD))), ID]
  
  SignalCurrent_State1 <- SignalCurrent[ID %in% State1[A == aat, ID]]
  SignalCurrent_State1$x <- SignalCurrent_State1$x/max(SignalCurrent_State1$x)*State1[A == aat, max(DwellTime)]
  
  SignalCurrent_State1$D <- SignalCurrent_State1[, get_density(x = x, y = Current)]
  
  ggplot(data = SignalCurrent_State1, mapping = aes(x = x, y = Current)) + 
    geom_line(aes(colour = ID)) + 
    geom_density_2d_filled(aes(fill = ..level..), contour_var = "ndensity", breaks = seq(0, 1, length.out = 100), alpha = .5) + 
    scale_colour_manual(values = rep("black", SignalCurrent_State1[, length(unique(ID))])) +
    scale_fill_manual(values = fill_cols, aesthetics = c("fill")) + 
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = rev(c("0.00", "0.25", "0.50", "0.75", "1.00"))) + 
    guides(colour = "none", fill = "none") + 
    labs(y = "Blockade", x = "Dwell time (ms)", title = aat) + 
    theme_light(base_size = 15) -> p
  # ggsave(paste0("./analysis/22.SignalSelecting/03.DifferentStates/00.SignalPlot/", aat, "_SignalEvent.pdf"), p, width = 4, height = 4)
  png(paste0("./analysis/22.SignalSelecting/03.DifferentStates/00.SignalPlot/", aat, "_SignalEvent.png"), width = 480, height = 480)
  print(p)
  dev.off()
}





