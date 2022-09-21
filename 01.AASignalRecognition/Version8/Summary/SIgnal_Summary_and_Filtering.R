library(data.table)
library(IRanges)
library(openxlsx)

mean2 <- function(x, q1 = 0.25, q2 = 0.75) {
  q <- quantile(x, c(q1, q2))
  mean(x[x > min(q) & x < max(q)])
}

mean3 <- function(x, HT = 10) {
  if(length(x) >= 5*HT) {
    mean(x[(HT+1):(length(x) - HT)])
  } else {
    HT <- floor(length(x)/5)
    mean(x[(HT+1):(length(x) - HT)])
  }
}

UpPeak <- function(x) {
  SigMean <- x[L == "U", mean3(pA)]
  rg <- range(which(x[, L == "U" & pA < SigMean]))
  x[rg[1]:rg[2], max(pA)] - SigMean
}


aas <- sort(unique(openxlsx::read.xlsx("/mnt/raid61/Personal_data/tangchao/AANanopore/data/meta_info_and_base_line_20210525.xlsx")[[3]]))

lapply(aas, function(AA) {
  print(paste0(AA, ": ", which(aas == AA), " of ", length(aas)))
  InDir <- file.path("/mnt/raid61/Personal_data/tangchao/AANanopore/analysis/01.AASignalRecognition/Version6", AA)
  OutDir <- file.path("/mnt/raid61/Personal_data/tangchao/AANanopore/analysis/01.AASignalRecognition/Version8/Summary/")
  if(!dir.exists(OutDir)) dir.create(OutDir, recursive = TRUE)
  BUBs <- readRDS(file.path(InDir, "BUB.Rds"))
  
  summa <- data.table(File = rep(names(BUBs), mapply(length, BUBs)), 
                      No = do.call(c, lapply(BUBs, seq_along)))
  
  BUBs <- do.call(c, BUBs)
  ty <- mapply(function(x) class(x)[1], BUBs)
  summa <- summa[ty == "data.table"]
  BUBs <- BUBs[ty == "data.table"]
  
  summa$StartTime <- mapply(BUBs, FUN = function(x) x[, min(Time)])
  summa$EndTime <- mapply(BUBs, FUN = function(x) x[, max(Time)])
  summa$BaseLineLength <- mapply(BUBs, FUN = function(x) min(width(IRanges(x[, L == "B"]))))
  summa$BaseLineMean <- mapply(BUBs, FUN = function(x) x[L == "B", mean2(pA)])
  summa$DwellTime <- mapply(BUBs, FUN = function(x) x[L == "U", diff(range(Time)) * 1000])
  summa$SignalMean <- mapply(BUBs, FUN = function(x) x[L == "U", mean3(pA)])
  summa$SignalSD <- mapply(BUBs, FUN = function(x) x[L == "U", sd(pA)])
  summa$SignalMAD <- mapply(BUBs, FUN = function(x) x[L == "U", mad(pA)])
  
  summa$BlockadeMaximum <- mapply(BUBs, FUN = function(x) 1 - x[L == "U", min(pA)]/x[L == "B", mean2(pA)])
  summa$SignalRangeDown <- mapply(BUBs, FUN = function(x) abs(x[L == "U", min(pA)] - x[L == "U", mean3(pA)]))
  summa$SignalRangeUp <- mapply(BUBs, FUN = function(x) UpPeak(x))
  summa$SignalRange <- apply(summa[, .(SignalRangeDown, SignalRangeUp)], 1, max)
  summa$BaseLineRange <- mapply(BUBs, FUN = function(x) mean(abs(x[L == "B", pA] - x[L == "B", mean2(pA)])))
  
  summa <- summa[SignalRangeUp > 0, ]
  summa[, SignalRangeDown := NULL]
  summa[, SignalRangeUp := NULL]
  summa[, Blockade := 1 - SignalMean/BaseLineMean]
  
  openxlsx::write.xlsx(summa, paste0(OutDir, AA, "_0.xlsx"))
  
  summa1 <- summa[SignalSD < 3.5, ]
  openxlsx::write.xlsx(summa1, paste0(OutDir, AA, "_1.xlsx"))
  
  summa2 <- summa[SignalRange / BaseLineMean / Blockade < 0.35, ]
  openxlsx::write.xlsx(summa2, paste0(OutDir, AA, "_2.xlsx"))
  
  summa3 <- summa[SignalRange / BaseLineMean / Blockade < 0.35 & SignalSD < 3.5, ]
  openxlsx::write.xlsx(summa3, paste0(OutDir, AA, "_3.xlsx"))
})


ind <- "/mnt/raid61/Personal_data/tangchao/AANanopore/analysis/01.AASignalRecognition/Version8/Summary"
alls <- list.files(ind, ".xlsx$", full.names = T)

lapply(alls, function(x) {
  summa <- openxlsx::read.xlsx(x)
  Lis <- split.data.frame(summa, summa$File)
  for(i in seq_along(Lis)) {
    out <- paste0(ind, "/", gsub("_", paste0("_", names(Lis)[i], "_"), basename(x)))
    openxlsx::write.xlsx(Lis[[i]], out)
  }
})


