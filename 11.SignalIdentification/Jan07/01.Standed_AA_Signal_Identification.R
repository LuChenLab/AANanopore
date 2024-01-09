setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")

library(data.table)
library(Biostrings)
library(IRanges)
library(ggplot2)
library(parallel)
library(changepoint)

KBins2 <- function(sig, minseglen1 = 100, minseglen2 = 10, pen.value = 1e-3) {
  ansmean1 <- suppressWarnings(changepoint::cpt.mean(sig, penalty = "MBIC", method = "PELT", minseglen = minseglen1))
  Tab1 <- data.table::data.table(P = sig, B1 = rep(seq_len(changepoint::nseg(ansmean1)), changepoint::seg.len(ansmean1)))
  Tab2 <- Tab1[, .(P = median(P)), by = "B1"]
  
  ansmean2 <- suppressWarnings(changepoint::cpt.meanvar(Tab2[, P], penalty = "Asymptotic", pen.value = pen.value, method = "PELT", minseglen = minseglen2))
  Tab2[, B := rep(seq_len(changepoint::nseg(ansmean2)), changepoint::seg.len(ansmean2))]
  res <- merge(Tab1, Tab2[, .(B1, B)], by = "B1")
  res[, B1 := NULL]
  res[, .(P = median(P), N = .N), "B"][, rep(P, N)]
}

mean2 <- function(x, q1 = 0.25, q2 = 0.75) {
  q <- quantile(x, c(q1, q2))
  mean(x[x > min(q) & x < max(q)])
}

WhiskerRange <- function(x) {
  iqr <- 1.5 * IQR(x, na.rm = TRUE)
  c(min(x[x >= quantile(x, 1/4, na.rm = T) - iqr], na.rm = TRUE), max(x[x <= quantile(x, 3/4, na.rm = T) + iqr], na.rm = TRUE))
}

GetSignal <- function(x, minblockade = 0.1) {
  abf2 <- x[rep(runLength(Rle(x$Sm)) >= 20, runLength(Rle(x$Sm)))]
  gr <- IRanges(end = cumsum(runLength(Rle(abf2[, Sm]))), width = runLength(Rle(abf2[, Sm])))
  Sm <- runValue(Rle(abf2[, Sm]))
  S0 <- list()
  i <- 1
  while (i < length(Sm)) {
    if(Sm[i + 1] < Sm[i] * (1 - minblockade)) {
      l0 <- which(Sm > Sm[i] * (1 - minblockade))
      l0 <- l0[l0 > i]
      if(length(l0) != 0) {
        Ei <- min(l0)
        S0 <- append(S0, list(c(end(gr[i]), start(gr[Ei]))))
        i <- Ei
      } else {
        i <- i + 1
      }
    } else {
      i <- i + 1
    }
  }
  
  Sigs <- lapply(S0, function(x) {
    L <- gr[end(gr) == x[1]]
    R <- gr[start(gr) == x[2]]
    gri <- gr[queryHits(findOverlaps(gr, IRanges(x[1], x[2]), type = "within"))]
    SignalCurrent <- with(density(abf2[(x[1] + 1):(x[2] - 1), pA], n = 10000, adjust = 1), x[which.max(y)])
    SignalCurrent.5 <- with(density(abf2[(x[1] + 1):(x[2] - 1), pA], n = 10000, adjust = .75), x[which.max(y)])
    L2Ratio = mean(abf2[(x[1] + 1):(x[2] - 1), Sm] < SignalCurrent * (1 - minblockade))
    StageSD = abf2[(x[1] + 1):(x[2] - 1), sd(pA), Sm][, mean(V1)]
    SignalCurrentPercent <- mean(abs(abf2[(x[1] + 1):(x[2] - 1), Sm] - SignalCurrent) < StageSD) * 100
    data.table(StartTime = abf2[x[1], Time], EndTime = abf2[x[2], Time], 
               LeftLength = width(L), RightLength = width(R), 
               BaseMean = mean2(c(abf2[start(L):end(L), pA], abf2[start(R):end(R), pA])),
               DeltaMean = abs(mean2(abf2[start(L):end(L), pA]) - mean2(abf2[start(R):end(R), pA])),
               StageSD = StageSD, 
               CurrentSD = abf2[(x[1] + 1):(x[2] - 1), sd(pA)], 
               Segments = sum(prop.table(width(gri)) > 0.1), Valid = all(abf2[start(L):end(R), Valid]), 
               SignalCurrent = SignalCurrent, SignalCurrent.5 = SignalCurrent.5, L2Ratio = L2Ratio, 
               SignalCurrentPercent = SignalCurrentPercent)
  })
  Sigs <- do.call(rbind, Sigs)
  return(Sigs)
}

SignalCurrent <- function(x, abf) {
  target <- abf[, inrange(Time, x$StartTime, x$EndTime, incbounds = FALSE)]
  target <- IRanges(target)
  target <- as.data.table(target)
  target$ID <- seq_len(nrow(target))
  target <- target[, .(Loci = list(start:end)), ID]
  target <- target[, Loci]
  names(target) <- x$ID
  abf2 <- abf[do.call(c, target)]
  abf2$ID <- rep(names(target), mapply(length, target))
  abf2$BaseMean <- rep(x[, BaseMean], mapply(length, target))
  abf2[, Current := pA / BaseMean]
  abf2[, .(ID, Current)]
}




meta0 <- data.table(openxlsx::read.xlsx("./data/单氨基酸实验信息表20210517.xlsx", sheet = 1, cols = 8:14))
meta1 <- data.table(openxlsx::read.xlsx("./data/单氨基酸实验信息表20210517.xlsx", sheet = 1, cols = 15:21))
meta2 <- data.table(openxlsx::read.xlsx("./data/单氨基酸实验信息表20210517.xlsx", sheet = 2, cols = 1:7))
meta3 <- data.table(openxlsx::read.xlsx("./data/单氨基酸实验信息表20210517.xlsx", sheet = 3, cols = 1:7))
meta4 <- data.table(openxlsx::read.xlsx("./data/单氨基酸实验信息表20210517.xlsx", sheet = 4, cols = 1:7))
meta <- rbind(meta0, meta1, meta2, meta3, meta4, use.names = FALSE)
colnames(meta) <- c("file_name", "date", "amino_acid", "concentration", "start_time", "end_time", "type")
meta <- meta[!is.na(file_name)]
meta[amino_acid == "cys", amino_acid := "Cys"]
meta <- meta[amino_acid %in% Biostrings::AMINO_ACID_CODE]
meta <- meta[concentration != 0]
meta$file_path <- mapply(meta$file_name, FUN = function(x) list.files("./data", recursive = TRUE, full.names = TRUE, pattern = as.character(x))[1])
meta <- meta[!is.na(file_path)]
meta[, start_time := as.numeric(start_time)]
meta[, end_time := as.numeric(end_time)]

meta_Val <- rbind(data.table(openxlsx::read.xlsx("./data/meta_info_and_base_line_20210525.xlsx")), data.table(openxlsx::read.xlsx("./data/meta_info_and_base_line_20210525additional.xlsx")))[amino_acid == "Val"]
meta_Val$file_path <- mapply(meta_Val$file_name, FUN = function(x) list.files("./data", recursive = TRUE, full.names = TRUE, pattern = as.character(x))[1])
meta_Val[, L1min := NULL]
meta_Val[, L1max := NULL]
meta <- rbind(meta, meta_Val)
meta <- meta[!grepl("abf", file_name)]

mclapply(1:nrow(meta), function(i) {
  print(i)
  if(!file.exists(paste0("./analysis/11.SignalIdentification/Dec27/ABF_", meta[i, file_name], ".Rds"))) return()
  abf <- readRDS(paste0("./analysis/11.SignalIdentification/Dec27/ABF_", meta[i, file_name], ".Rds"))
  
  Sigs <- GetSignal(x = abf)
  Sigs[, DwellTime := EndTime - StartTime]
  Sigs <- data.table(ID = paste(meta[i, file_name], Sigs[, seq_len(.N)], sep = "_"), Sigs)
  fwrite(Sigs, file = paste0("./analysis/11.SignalIdentification/Jan07/RawSignal_", meta[i, file_name], ".txt"), sep = "\t", row.names = F, quote = F)
}, mc.cores = 10)

meta[!file_name %in% gsub(".txt", "", gsub("RawSignal_", "", list.files("./analysis/11.SignalIdentification/Jan07")))]
meta <- meta[file_name %in% gsub(".txt", "", gsub("RawSignal_", "", list.files("./analysis/11.SignalIdentification/Jan07")))]

meta[, File := paste0("./analysis/11.SignalIdentification/Jan07/RawSignal_", file_name, ".txt")]

lapply(seq_len(nrow(meta)), function(i) {
  sig <- fread(meta[i, File])
  sig <- sig[Valid == TRUE]
  sig[, Blockade := 1 - SignalCurrent / BaseMean]
  sig[, DwellTime := DwellTime * 1000]
  sig[, A := meta[i, amino_acid]]
  fwrite(sig, gsub("/Jan07/", "/Jan07/AASignal/", meta[i, File]), sep = "\t", row.names = F, quote = F)
})




# mclapply(1:nrow(meta), function(i) {
#   print(i)
#   if(file.exists(paste0("./analysis/11.SignalIdentification/Dec27/SignalCurrent_", meta[i, file_name], ".Rds"))) return(NULL)
#   if(!file.exists(paste0("./analysis/11.SignalIdentification/Dec27/ABF_", meta[i, file_name], ".Rds"))) return(NULL)
#   if(!file.exists(paste0("./analysis/11.SignalIdentification/Dec27/RawSignal_", meta[i, file_name], ".txt"))) return(NULL)
#   abf <- readRDS(paste0("./analysis/11.SignalIdentification/Dec27/ABF_", meta[i, file_name], ".Rds"))
#   Sigs <- fread(paste0("./analysis/11.SignalIdentification/Dec27/RawSignal_", meta[i, file_name], ".txt"))
#   Current <- SignalCurrent(Sigs, abf = abf, cores = 20)
#   saveRDS(Current, file = paste0("./analysis/11.SignalIdentification/Dec27/SignalCurrent_", meta[i, file_name], ".Rds"))
# }, mc.cores = 1)




meta1 <- data.table(openxlsx::read.xlsx("./data/单氨基酸实验信息表20210517.xlsx", sheet = 1, cols = 1:7))
meta2 <- data.table(openxlsx::read.xlsx("./data/单氨基酸实验信息表20210517.xlsx", sheet = 2, cols = 1:7))
meta3 <- data.table(openxlsx::read.xlsx("./data/单氨基酸实验信息表20210517.xlsx", sheet = 3, cols = 1:7))
meta4 <- data.table(openxlsx::read.xlsx("./data/单氨基酸实验信息表20210517.xlsx", sheet = 4, cols = 1:7))
meta <- rbind(meta1, meta2, meta3, meta4, use.names = FALSE)
colnames(meta) <- c("file_name", "date", "amino_acid", "concentration", "start_time", "end_time", "type")
meta[amino_acid == "cys", amino_acid := "Cys"]
meta <- meta[amino_acid %in% Biostrings::AMINO_ACID_CODE]
meta <- meta[concentration == 0]
meta$file_path <- mapply(meta$file_name, FUN = function(x) list.files("./data", recursive = TRUE, full.names = TRUE, pattern = as.character(x))[1])
meta <- meta[!is.na(file_path)]
meta[, start_time := as.numeric(start_time)]
meta[, end_time := as.numeric(end_time)]
meta <- meta[!grepl("abf", file_name)]


mclapply(1:nrow(meta), function(i) {
  print(i)
  if(!file.exists(paste0("./analysis/11.SignalIdentification/Dec27/ABF_", meta[i, file_name], ".Rds"))) return()
  abf <- readRDS(paste0("./analysis/11.SignalIdentification/Dec27/ABF_", meta[i, file_name], ".Rds"))
  
  Sigs <- GetSignal(x = abf)
  Sigs[, DwellTime := EndTime - StartTime]
  Sigs <- data.table(ID = paste(meta[i, file_name], Sigs[, seq_len(.N)], sep = "_"), Sigs)
  fwrite(Sigs, file = paste0("./analysis/11.SignalIdentification/Jan07/RawSignal_", meta[i, file_name], ".txt"), sep = "\t", row.names = F, quote = F)
}, mc.cores = 10)


meta[!file_name %in% gsub(".txt", "", gsub("RawSignal_", "", list.files("./analysis/11.SignalIdentification/Jan07")))]
meta <- meta[file_name %in% gsub(".txt", "", gsub("RawSignal_", "", list.files("./analysis/11.SignalIdentification/Jan07")))]

meta[, File := paste0("./analysis/11.SignalIdentification/Jan07/RawSignal_", file_name, ".txt")]

lapply(seq_len(nrow(meta)), function(i) {
  sig <- fread(meta[i, File])
  sig <- sig[Valid == TRUE]
  sig[, Blockade := 1 - SignalCurrent / BaseMean]
  sig[, DwellTime := DwellTime * 1000]
  sig[, A := meta[i, amino_acid]]
  fwrite(sig, gsub("/Jan07/", "/Jan07/BackgroundSignal/", meta[i, File]), sep = "\t", row.names = F, quote = F)
})

