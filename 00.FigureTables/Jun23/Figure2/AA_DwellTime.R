setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(openxlsx)
library(ggplot2)
library(data.table)
library(Biostrings)
library(ggpubr)

meta <- do.call(rbind, lapply(list.files("./analysis/01.AASignalRecognition/Version9/02.RangesL0L1/", ".txt", full.names = TRUE, recursive = T), fread))
# meta <- meta[, .(file_name, amino_acid, concentration, start_time, end_time, TimeRatio)]
meta[, file_name := as.character(file_name)]
meta <- meta[!file_name %in% c("21205012", "21205021", "21201008", "21303003")]

files <- list.files("./analysis/01.AASignalRecognition/Version9/04.FinalSignal", pattern = "_Sigs.txt", recursive = TRUE, full.names = T)
AAs <- lapply(files, fread)
names(AAs) <- gsub("_Sigs.txt", "", basename(files))
for(i in seq_along(AAs)) AAs[[i]] <- data.table(file_name = names(AAs)[i], AAs[[i]])
AA <- do.call(rbind, AAs)
AA <- AA[!file_name %in% c("21205012", "21205021", "21201008", "21303003")]

AA <- merge(meta, AA, by = "file_name")
AA <- AA[!is.na(Blockade) & Outer == 0 & AreaRatio_L1 > 0.1]
AA[, .N, amino_acid]
AA[, aa := plyr::mapvalues(amino_acid, AMINO_ACID_CODE, names(AMINO_ACID_CODE))]

AA[aa %in% c("E", "D", "H", "R", "K"), Class := "charged"]
AA[aa %in% c("L", "I", "M", "V", "A", "F", "G", "W", "P"), Class := "nonpolar"]
AA[aa %in% c("S", "N", "Q", "T", "Y", "C"), Class := "polar"]


load(file = "./analysis/03.MachineLearning/01.data/Version6/AAInfo_RawSig.RData")
# AAInfo <- AAInfo[!amino_acid %in% c("Cys", "Pro")]
AAInfo$AllTime <- mapply(AAInfo$ID, FUN = function(x) {
  RawSig[ID == x, diff(range(Time))] * 1000
})
AAInfo[, aa := plyr::mapvalues(amino_acid, AMINO_ACID_CODE, names(AMINO_ACID_CODE))]
AAInfo[aa %in% c("E", "D", "H", "R", "K"), Class := "charged"]
AAInfo[aa %in% c("L", "I", "M", "V", "A", "F", "G", "W", "P"), Class := "nonpolar"]
AAInfo[aa %in% c("S", "N", "Q", "T", "Y", "C"), Class := "polar"]

od <- AAInfo[, mean(DwellTime), c("aa", "Class")][order(Class, V1), as.character(aa)]
AAInfo[, aa := factor(aa, levels = od)]

saveRDS(AAInfo, file = "./analysis/00.FigureTables/Jun23/Figure2/AA_DwellTime.Rds")

