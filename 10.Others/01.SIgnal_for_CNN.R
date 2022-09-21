setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")

mean2 <- function(x, q1 = 0.25, q2 = 0.75) {
  q <- quantile(x, c(q1, q2))
  mean(x[x >= min(q) & x <= max(q)])
}


library(data.table)
library(openxlsx)
library(parallel)


xlsxfiles <- list.files("./analysis/05.MixedAA/Version1", ".xlsx", full.names = TRUE)
Signalsfiles <- list.files("./analysis/05.MixedAA/Version1", ".Rds", full.names = TRUE)


# Mixed AA


lapply(seq_along(Signalsfiles), function(i) {
  Signal <- as.data.table(read.xlsx(xlsxfiles[i]))
  Signal$File <- gsub(".xlsx", "", basename(xlsxfiles[i]))
  Signal$ID <- paste0(gsub(".xlsx", "", basename(xlsxfiles[i])), "_", sprintf("%04d", seq_len(Signal[, .N])))
  
  BUBs <- readRDS(Signalsfiles[i])
  names(BUBs) <- Signal$ID
  BUBs <- mclapply(BUBs, function(x) {
    basemean <- x[L == "B", mean2(pA)]
    x[, pA := pA / basemean]
    x[, Sm := NULL]
    x[L == "U"]
  }, mc.cores = 10)
  
  BUBs <- data.table(ID = rep(names(BUBs), mapply(nrow, BUBs)), do.call(rbind, BUBs))
  return(BUBs)
}) -> Sig_List
Sigs <- do.call(rbind, Sig_List)
Sigs[, L := NULL]


lapply(seq_along(Signalsfiles), function(i) {
  Signal <- as.data.table(read.xlsx(xlsxfiles[i]))
  Signal$File <- gsub(".xlsx", "", basename(xlsxfiles[i]))
  Signal$ID <- paste0(gsub(".xlsx", "", basename(xlsxfiles[i])), "_", sprintf("%04d", seq_len(Signal[, .N])))
  return(Signal)
}) -> Info_List
Infos <- do.call(rbind, Info_List)

fwrite(Sigs, "/mnt/raid61/Personal_data/tangchao/Temp/MixedAA_Current_for_MZ.csv", row.names = F, sep = ",")
fwrite(Infos, "/mnt/raid61/Personal_data/tangchao/Temp/MixedAA_Infos_for_MZ.csv", row.names = F, sep = ",")






# hydrolysate


# Polypeptide1 

Signal <- as.data.table(read.xlsx("./analysis/02.PolypeptideSequencing/20211025/Version2/Polypeptide1_APRLRFYSL.xlsx"))
Signal$ID <- paste0(Signal[, Sample], "_", sprintf("%04d", do.call(c, lapply(Signal[, .N, Sample][, N], function(x) seq_len(x)))))

BUBs <- paste0("./analysis/02.PolypeptideSequencing/20211025/Version2/BUB", Signal[, unique(Sample)], ".Rds")
BUBs <- lapply(BUBs, readRDS)
BUBs <- do.call(c, BUBs)
names(BUBs) <- Signal$ID

BUBs <- mclapply(BUBs, function(x) {
  basemean <- x[L == "B", mean2(pA)]
  x[, pA := pA / basemean]
  x[, Sm := NULL]
  x[L == "U"]
}, mc.cores = 10)

BUBs <- data.table(ID = rep(names(BUBs), mapply(nrow, BUBs)), do.call(rbind, BUBs))
BUBs[, L := NULL]

fwrite(BUBs, "/mnt/raid61/Personal_data/tangchao/Temp/Polypeptide1_Current_for_MZ.csv", row.names = F, sep = ",")
fwrite(Signal, "/mnt/raid61/Personal_data/tangchao/Temp/Polypeptide1_Infos_for_MZ.csv", row.names = F, sep = ",")


# Polypeptide2

Signal <- as.data.table(read.xlsx("./analysis/02.PolypeptideSequencing/20211025/Version2/Polypeptide2_RPVKVYPNGAEDESAEAFPLEF.xlsx"))
Signal$ID <- paste0(Signal[, Sample], "_", sprintf("%04d", do.call(c, lapply(Signal[, .N, Sample][, N], function(x) seq_len(x)))))

BUBs <- paste0("./analysis/02.PolypeptideSequencing/20211025/Version2/BUB", Signal[, unique(Sample)], ".Rds")
BUBs <- lapply(BUBs, readRDS)
BUBs <- do.call(c, BUBs)
names(BUBs) <- Signal$ID

BUBs <- mclapply(BUBs, function(x) {
  basemean <- x[L == "B", mean2(pA)]
  x[, pA := pA / basemean]
  x[, Sm := NULL]
  x[L == "U"]
}, mc.cores = 10)

BUBs <- data.table(ID = rep(names(BUBs), mapply(nrow, BUBs)), do.call(rbind, BUBs))
BUBs[, L := NULL]

fwrite(BUBs, "/mnt/raid61/Personal_data/tangchao/Temp/Polypeptide2_Current_for_MZ.csv", row.names = F, sep = ",")
fwrite(Signal, "/mnt/raid61/Personal_data/tangchao/Temp/Polypeptide2_Infos_for_MZ.csv", row.names = F, sep = ",")



# Polypeptide3

Signal <- as.data.table(read.xlsx("./analysis/02.PolypeptideSequencing/20211025/Version2/Polypeptide3_DRVYIHPFHL.xlsx"))
Signal$ID <- paste0(Signal[, Sample], "_", sprintf("%04d", do.call(c, lapply(Signal[, .N, Sample][, N], function(x) seq_len(x)))))

BUBs <- paste0("./analysis/02.PolypeptideSequencing/20211025/Version2/BUB", Signal[, unique(Sample)], ".Rds")
BUBs <- lapply(BUBs, readRDS)
BUBs <- do.call(c, BUBs)
names(BUBs) <- Signal$ID

BUBs <- mclapply(BUBs, function(x) {
  basemean <- x[L == "B", mean2(pA)]
  x[, pA := pA / basemean]
  x[, Sm := NULL]
  x[L == "U"]
}, mc.cores = 10)

BUBs <- data.table(ID = rep(names(BUBs), mapply(nrow, BUBs)), do.call(rbind, BUBs))
BUBs[, L := NULL]

fwrite(BUBs, "/mnt/raid61/Personal_data/tangchao/Temp/Polypeptide3_Current_for_MZ.csv", row.names = F, sep = ",")
fwrite(Signal, "/mnt/raid61/Personal_data/tangchao/Temp/Polypeptide3_Infos_for_MZ.csv", row.names = F, sep = ",")


