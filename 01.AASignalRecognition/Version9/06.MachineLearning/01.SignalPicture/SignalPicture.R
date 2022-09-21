setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(ggplot2)
library(parallel)

mean2 <- function(x, q1 = 0.25, q2 = 0.75) {
  q <- quantile(x, c(q1, q2))
  mean(x[x >= min(q) & x <= max(q)])
}

plotSig <- function(x, title = NULL) {
  ggplot() +
    geom_line(data = x, mapping = aes(x = Time, y = pA)) +
    theme_minimal(base_size = 15) +
    labs(title = title)
}

# c("21205012", "21205021", "21201008", "21303003")

meta <- do.call(rbind, lapply(list.files("./analysis/01.AASignalRecognition/Version9/02.RangesL0L1/", ".txt", full.names = TRUE, recursive = T), fread))
meta <- meta[!file_name %in% c("21205012", "21205021", "21201008", "21303003")]
meta <- meta[, .(file_name, amino_acid, concentration, start_time, end_time)]
meta[, file_name := as.character(file_name)]


files <- list.files("./analysis/01.AASignalRecognition/Version9/04.FinalSignal", pattern = "_Sigs.txt", recursive = TRUE, full.names = T)
AAs <- lapply(files, fread)
names(AAs) <- gsub("_Sigs.txt", "", basename(files))
AAs <- AAs[!names(AAs) %in% c("21205012", "21205021", "21201008", "21303003")]

for(i in seq_along(AAs)) AAs[[i]] <- data.table(file_name = names(AAs)[i], ID = paste0(names(AAs)[i], "_", sprintf("%05d", seq_len(nrow(AAs[[i]])))), AAs[[i]])
AA <- do.call(rbind, AAs)

AA <- merge(AA, meta[, .(file_name, amino_acid)], by = "file_name")
setkey(AA, ID)


files <- list.files("./analysis/01.AASignalRecognition/Version9/03.BUBSignal", pattern = "_BUB.Rds", recursive = TRUE, full.names = T)
BUBs <- mclapply(files, readRDS, mc.cores = 10)
names(BUBs) <- gsub("_BUB.Rds", "", basename(files))
BUBs <- BUBs[!names(BUBs) %in% c("21205012", "21205021", "21201008", "21303003")]
for(i in seq_along(BUBs)) names(BUBs[[i]]) <- paste0(names(BUBs)[i], "_", sprintf("%05d", seq_along(BUBs[[i]])))
names(BUBs) <- NULL
BUBs <- do.call(c, BUBs)
BUBs <- BUBs[AA$ID]

BUBs <- mclapply(BUBs, function(x) {
  basemean <- x[L == "B", mean2(pA)]
  x[, pA := pA / basemean]
  x[, Sm := NULL]
  x
}, mc.cores = 10)

for(i in seq_along(BUBs)) BUBs[[i]] <- data.table(ID = names(BUBs)[i], BUBs[[i]])


set.seed(123)
AA_tu <- AA[, .SD[sample(.N, 660), ], amino_acid]
BUBs_tu <- BUBs[AA_tu$ID]

mclapply(seq_along(BUBs_tu), function(i) {
  p <- plotSig(BUBs_tu[[i]], title = names(BUBs_tu)[i]) + lims(y = c(0, 1.2))
  output <- paste0("./analysis/01.AASignalRecognition/Version9/06.MachineLearning/01.SignalPicture2/", AA_tu[i, amino_acid], "_", names(BUBs_tu)[i], ".pdf")
  ggsave(output, p)
}, mc.cores = 10)

BUBs_tu_tab <- do.call(rbind, BUBs_tu)
BUBs_tu_tab
fwrite(BUBs_tu_tab, "./analysis/01.AASignalRecognition/Version9/06.MachineLearning/01.SignalPicture2/RawSignal.txt", sep = "\t", row.names=F, quote=F)

fwrite(AA_tu, "./analysis/01.AASignalRecognition/Version9/06.MachineLearning/01.SignalPicture2/AASignal.txt", sep = "\t", row.names=F, quote=F)

BUBs_tu_tab <- fread("./analysis/01.AASignalRecognition/Version9/06.MachineLearning/01.SignalPicture2/RawSignal.txt")

mclapply(BUBs_tu_tab[, unique(ID)], function(x) {
  output <- paste0("./analysis/01.AASignalRecognition/Version9/06.MachineLearning/01.SignalPicture2/RawSignal/Signal_", x, ".txt")
  fwrite(BUBs_tu_tab[ID == x, 2:5], output, sep = "\t", row.names=F, quote=F)
}, mc.cores = 10)











