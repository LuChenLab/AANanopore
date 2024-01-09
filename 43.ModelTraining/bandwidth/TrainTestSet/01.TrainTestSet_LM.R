setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(shiny)
library(plotly)
library(Biostrings)

set.seed(19)
AA_Cols <- sample(c(RColorBrewer::brewer.pal(n = 12, "Paired"), RColorBrewer::brewer.pal(n = 8, "Dark2")), 20)
names(AA_Cols) <- Biostrings::AA_STANDARD
AA_Cols[AA_Cols == "#FFFF99"] <- "#1A1A1A"

get_density <- function(x, y, n = 1000, ...) {
  dens <- MASS::kde2d(x, y, n = n, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii] / max(dens$z[ii]))
}


AABlockade <- data.table(AA = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"), 
                         Blockade = c(0.14718, 0.17148, 0.16516, 0.21409, 0.16, 0.24528, 0.1882, 0.12072, 0.24652, 0.20722, 0.1995, 0.16875, 0.19772, 0.22018, 0.2183, 0.13131, 0.16101, 0.22744, 0.21276, 0.19044))
AABlockade$amino_acid <- plyr::mapvalues(AABlockade$AA, names(Biostrings::AMINO_ACID_CODE), Biostrings::AMINO_ACID_CODE)

AAGroup <- list(GSA = c('G', 'S', 'A'), 
                TNRK = c("T", "N", "R", "K"), 
                QVML = c("Q", "V", "M", "L"), 
                IYDPFW = c("I", "Y", "D", "P", "F", "W"), 
                EH = c("E", "H"))
AAGroup <- data.table(Group = rep(names(AAGroup), mapply(length, AAGroup)), AA = unlist(AAGroup))
AAGroup[, amino_acid := plyr::mapvalues(AA, names(AMINO_ACID_CODE), AMINO_ACID_CODE)]


Sigs <- lapply(AAGroup[, amino_acid], function(aat) {
  sigs <- do.call(rbind, lapply(list.files(paste0("./analysis/42.SignalSelecting/01.StandardAA/", aat), full.names = T), fread))
  sigs$File <- stringr::str_remove_all(sigs$ID, "_([[:digit:]]+)$")
  sigs$D <- sigs[, get_density(x = Blockade, y = log10(DwellTime))]
  sigs
})
Sigs <- do.call(rbind, Sigs)
Sigs <- merge(Sigs, AAGroup, by.x = "A", by.y = "amino_acid")

ggplot(Sigs[Group == "GSA" & D > .1], aes(x = Blockade, y = DwellTime, colour = AA)) + 
  geom_point() + 
  scale_y_log10()
ggplot(Sigs[Group == "GSA" & D > .1], aes(x = Blockade, fill = AA)) + 
  geom_histogram(binwidth = 0.0003)


ggplot(Sigs[Group == "TNRK" & D > .1], aes(x = Blockade, y = DwellTime, colour = AA)) + 
  geom_point() + 
  scale_y_log10()
ggplot(Sigs[Group == "TNRK" & D > .1], aes(x = Blockade, fill = AA)) + 
  geom_histogram(binwidth = 0.0003)


ggplot(Sigs[Group == "QVML" & D > .1], aes(x = Blockade, y = DwellTime, colour = AA)) + 
  geom_point() + 
  scale_y_log10()
ggplot(Sigs[Group == "QVML" & D > .1], aes(x = Blockade, fill = AA)) + 
  geom_histogram(binwidth = 0.0003)



ggplot(Sigs[Group == "IYDPFW" & D > .1], aes(x = Blockade, y = DwellTime, colour = AA)) + 
  geom_point() + 
  scale_y_log10()
ggplot(Sigs[Group == "IYDPFW" & D > .1], aes(x = Blockade, fill = AA)) + 
  geom_histogram(binwidth = 0.0003)


ggplot(Sigs[Group == "EH" & D > .1], aes(x = Blockade, y = DwellTime, colour = AA)) + 
  geom_point() + 
  scale_y_log10()
ggplot(Sigs[Group == "EH" & D > .1], aes(x = Blockade, fill = AA)) + 
  geom_histogram(binwidth = 0.0003)


ggplot(Sigs[D > .15], aes(x = Blockade, y = DwellTime, colour = AA)) + 
  geom_point(alpha = .5) + 
  scale_y_log10() + 
  facet_wrap(~ Group, scales = 'free_x') + 
  scale_colour_manual(breaks = names(AA_Cols), values = AA_Cols)

ggplot(Sigs[D > .15], aes(x = Blockade, fill = AA, colour = AA)) + 
  geom_density(alpha = .5, adjust = 2) + 
  facet_wrap(~ Group, scales = 'free_x') + 
  scale_fill_manual(breaks = names(AA_Cols), values = AA_Cols) + 
  scale_colour_manual(breaks = names(AA_Cols), values = AA_Cols)

SigsID <- Sigs[D > .1, .(ID, AA)]
SigsID_TU <- SigsID[AA %in% c('L', 'M')]

set.seed(123)
Test_Set <- SigsID_TU[, .SD[sample(.N, .N - 2000)], AA]
Train_Set <- SigsID_TU[!ID %in% Test_Set[, ID]]

save(Test_Set, Train_Set, file = "./analysis/43.ModelTraining/bandwidth/TrainTestSet/TrainTestSet.RData")





SigsID_TU <- SigsID[AA %in% c('N', 'R')]

set.seed(123)
Test_Set <- SigsID_TU[, .SD[sample(.N, .N - 2000)], AA]
Train_Set <- SigsID_TU[!ID %in% Test_Set[, ID]]

save(Test_Set, Train_Set, file = "./analysis/43.ModelTraining/bandwidth/TrainTestSet/TrainTestSet_NR.RData")


