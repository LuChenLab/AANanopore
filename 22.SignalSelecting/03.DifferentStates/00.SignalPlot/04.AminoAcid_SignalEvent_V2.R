setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(patchwork)
library(ggplot2)
library(ggrepel)

get_density <- function(x, y, n = 1000, ...) {
  dens <- MASS::kde2d(x, y, n = n, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii] / max(dens$z[ii]))
}

Peaks <- function(x, y) {
  stopifnot(length(x) == length(y))
  uR <- IRanges(diff(y) >= 0)
  dR <- IRanges(diff(y) <= 0)
  res <- do.call(c, lapply(seq_along(uR), function(i) {
    reduce(c(uR[i], dR[i]))
  }))
  
  Ps <- mapply(seq_along(res), FUN = function(i) {
    x[start(res[i]):end(res[i])][which.max(y[start(res[i]):end(res[i])])]
  })
  res <- IRanges(start = x[start(res)], end = x[end(res)])
  end(res[length(res)]) <- max(x)
  attr(res, "peaks") <- Ps
  return(res)
}

densityPeak <- function(x, bw = 0.01, n = 1024, from = 0, to = 1, plot = T, ...) {
  den <- density(x, bw = bw, n = n, ...)
  ps <- Peaks(den$x, den$y)
  dx <- attr(ps, "peaks")
  dy <- den$y[den$x %in% attr(ps, "peaks")]
  res <- c(dx[which.max(dy)], dx[dy == max(dy[dx < dx[which.max(dy)]])])
  if(plot) {
    plot(den)
    abline(v = res, lty = 2, col = 2)
  } else {
    return(res)
  }
}


set.seed(19)
AA_Cols <- sample(c(RColorBrewer::brewer.pal(n = 12, "Paired"), RColorBrewer::brewer.pal(n = 8, "Dark2")), 20)
names(AA_Cols) <- Biostrings::AA_STANDARD
AA_Cols[AA_Cols == "#FFFF99"] <- "#1A1A1A"
names(AA_Cols) <- AMINO_ACID_CODE[Biostrings::AA_STANDARD]
names(AA_Cols)[names(AA_Cols) == "Cys"] <- "CbC"

files1 <- list.files("./analysis/22.SignalSelecting/03.DifferentStates", "_State1.txt", recursive = T, full.names = T)
files1 <- grep("His2", files1, value = T, invert = T)
files2 <- list.files("./analysis/22.SignalSelecting/03.DifferentStates", "_State2.txt", recursive = T, full.names = T)

State1 <- do.call(rbind, lapply(files1, fread))
State1[, State := "State1"]

ggplot(State1[A %in% c("Ser", "Ala", "Thr", "Lys", "Gln", "Leu", "Ile", "Ile", "Pro", "Trp", "Glu", "His", "CbC")], 
       aes(x = Blockade, y = DwellTime)) + 
  geom_point(aes(colour = A)) + 
  scale_y_log10() + 
  scale_colour_manual(values = AA_Cols) + 
  theme_light(base_size = 15)













fill_cols <- colorRampPalette(c(rep("white", 1), RColorBrewer::brewer.pal(n = 9, name = "Reds")[4:9]))(100)

for (aat in State1[, unique(A)]) {
  print(aat)
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
    geom_line(aes(colour = ID), size = .1) + 
    geom_density_2d_filled(aes(fill = ..level..), contour_var = "ndensity", breaks = seq(0, 1, length.out = 100), alpha = .5) + 
    scale_colour_manual(values = rep("black", SignalCurrent_State1[, length(unique(ID))])) +
    scale_fill_manual(values = fill_cols, aesthetics = c("fill")) + 
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1), labels = rev(c("0.00", "0.25", "0.50", "0.75", "1.00"))) + 
    guides(colour = "none", fill = "none") + 
    labs(y = "Blockade", x = "Dwell time (ms)", title = aat) + 
    theme_light(base_size = 15) -> p1
  
  den <- density(SignalCurrent_State1[, Current], from = 0, to = 1, bw = 0.01, n = 1000)
  den <- with(den, data.table(x, y))
  hist <- hist(SignalCurrent_State1[, Current], breaks = 100, plot = F)
  den$y <- den$y/max(den$y)*max(hist$counts)
  # hist(SignalCurrent_State1[, Current], breaks = 100)
  # lines(den$x, den$y)
  # densityPeak(x = SignalCurrent_State1[, Current], plot = T)
  
  ps <- densityPeak(x = SignalCurrent_State1[, Current], plot = F)
  ps <- data.table(x = ps, y = mapply(ps, FUN = function(i) den[which.min(abs(x - i)), y]))
  ps[which.max(y), L := "State1"]
  ps[which.min(y), L := "State2"]
  ps[, B := round(1 - x, 3)]
  
  ggplot(SignalCurrent_State1, aes(x = Current)) + 
    geom_histogram(binwidth = 0.01, fill = "grey") +
    scale_x_continuous(limits = c(0, 1), breaks = ps[, x], labels = ps[, L]) + 
    scale_y_continuous(limits = c(0, ps[, max(y)] * 1.25)) + 
    geom_line(data = den, mapping = aes(x = x, y = y), colour = "#A50F15") + 
    geom_text_repel(data = ps, mapping = aes(x = x, y = y, label = B), 
                    min.segment.length = 0, nudge_y = 1, 
                    direction = "x") + 
    coord_flip() + 
    theme_minimal() + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text.x = element_blank()) -> p2
  
  p12 <- p1 + p2 + plot_layout(widths = c(2, 1), nrow = 1)
  ggsave(paste0("./analysis/22.SignalSelecting/03.DifferentStates/00.SignalPlot/", aat, "_SignalEvent.pdf"), p12, width = 6, height = 4)
  # ggsave(paste0("./analysis/22.SignalSelecting/03.DifferentStates/00.SignalPlot/", aat, "_SignalEvent.pdf"), p, width = 4, height = 4)
  # png(paste0("./analysis/22.SignalSelecting/03.DifferentStates/00.SignalPlot/", aat, "_SignalEvent.png"), width = 480, height = 480)
  # print(p)
  # dev.off()
}



