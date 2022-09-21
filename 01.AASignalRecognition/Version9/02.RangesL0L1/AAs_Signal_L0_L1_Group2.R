setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(parallel)
library(openxlsx)
library(plotly)
library(crosstalk)


meta1 <- data.table(openxlsx::read.xlsx("./data/meta_info_and_base_line_20210525.xlsx"))
meta2 <- data.table(openxlsx::read.xlsx("./data/meta_info_and_base_line_20210525additional.xlsx"))
meta <- rbind(meta1, meta2)
meta[amino_acid == "cys", amino_acid := "Cys"]
meta1 <- meta[amino_acid %in% c("Lys", "Cys", "Gly", "Thr")]
meta1$file_path <- mapply(meta1$file_name, FUN = function(x) list.files("./data", recursive = TRUE, full.names = TRUE, pattern = as.character(x)))

meta1 <- meta1[file_name %in% gsub(".Rds", "", list.files("./analysis/01.AASignalRecognition/Version9/01.SignalPolish/Group2/"))]

for(i in seq_len(nrow(meta1))) {
  print(paste0(i, " of ", nrow(meta1)))
  print(meta1[i, ])
  
  if(file.exists(paste0("./analysis/01.AASignalRecognition/Version9/02.RangesL0L1/Group2/", meta1[i, file_name], ".txt"))) next()
  
  L1min <- meta1[i, L1min]
  L1max <- meta1[i, L1max]
  
  MinL0 <- tryCatch(meta1[i, L0min], error = function(e) NA)
  MaxL0 <- tryCatch(meta1[i, L0max], error = function(e) NA)
  
  abf <- readRDS(paste0("./analysis/01.AASignalRecognition/Version9/01.SignalPolish/Group2/", meta1[i, file_name], ".Rds"))
  
  if(anyNA(c(MinL0, MaxL0))) {
    print("We need identify the baseline from data.")
    pA <- abf[, pA]
    pADen <- do.call(rbind, lapply(1:10/5, function(i) {
      d <- density(pA, adjust = i, n = 1024)
      data.table(x = d$x, y = d$y, adjust = i)
    }))
    
    tx <- highlight_key(pADen)
    widgets <- bscols(
      widths = c(12),
      filter_select("adjust", "Adjust", tx, ~ adjust, multiple = F)
    )
    print(bscols(widths = c(2, 8), widgets, plot_ly(tx, x = ~ x, y = ~ y, showlegend = FALSE) %>% 
                   add_lines(color = ~ adjust, colors = "black", showlegend = FALSE)))
    
    # askYesNo(msg = "Record the range of base lines and continue?")
    
    adjust <- readline("The adjust: ")
    adjust <- as.numeric(adjust)
    
    MinL0 <- readline("The minimum of baseline(pA): ")
    MinL0 <- as.numeric(MinL0)
    
    MaxL0 <- readline("The maximum of baseline(pA): ")
    MaxL0 <- as.numeric(MaxL0)
    
    DenSm <- density(abf[, pA], adjust = adjust)
    
    # L0 <- CISelect(abf[, pA], adjust = adjust, plot = F)
    # MinL0 <- min(L0)
    # MaxL0 <- max(L0)
    L0 <- c(MinL0, MaxL0)
    
    L0Min <- abf[pA > MinL0 & pA < MaxL0, median(Sm)] * 0.9
    
    plot(DenSm, xlab = "Current (pA)", main = "density of current")
    abline(v = L0, lty = 2, col = 2)
    abline(v = L0Min, lty = 2, col = 3)
    askYesNo(msg = "Continue?")
  }
  
  MinL1 <- tryCatch(meta1[i, MinL1], error = function(e) NA)
  MaxL1 <- tryCatch(meta1[i, MaxL1], error = function(e) NA)
  
  if(anyNA(c(MinL1, MaxL1))) {
    print("Now, we need identify the L1 from data.")
    
    pA <- abf[Sm < max(L0) * 0.9 & Sm > 50, pA]
    
    pADen <- do.call(rbind, lapply(1:10/10, function(i) {
      d <- density(pA, adjust = i, n = 1024)
      data.table(x = d$x, y = d$y, adjust = i)
    }))
    
    if(meta1[i, amino_acid] == "Cys") {
      pADen$x2 <- MinL0 * (1 - L1max)
      pADen$x3 <- MaxL0 * (1 - L1min)
    } else {
      pADen$x2 <- abf[pA > MinL0 & pA < MaxL0, median(Sm)] * c(1 - L1min)
      pADen$x3 <- abf[pA > MinL0 & pA < MaxL0, median(Sm)] * c(1 - L1max)
    }
    
    tx <- highlight_key(pADen)
    widgets <- bscols(
      widths = c(12),
      filter_select("adjust", "Adjust", tx, ~ adjust, multiple = F)
    )
    
    print(bscols(widths = c(2, 8), widgets, plot_ly(tx, x = ~ x, y = ~ y, showlegend = FALSE) %>% 
                   add_lines(color = ~ adjust, colors = "black", showlegend = FALSE) %>% 
                   add_lines(x = ~ x2) %>%
                   add_lines(x = ~ x3)))
    
    pd <- select.list(choices = c("YES", "NO"), title = "The theoretical value is OK?")
    
    if(pd == "YES") {
      MinL1 <- MinL0 * (1 - L1max)
      MaxL1 <- MaxL0 * (1 - L1min)
    } else {
      MinL1 <- readline("The minimum of L1(pA): ")
      MinL1 <- as.numeric(MinL1)
      
      MaxL1 <- readline("The maximum of L1(pA): ")
      MaxL1 <- as.numeric(MaxL1)
    }
  }
  
  res <- data.table(meta1[i, ], MinL0 = MinL0, MaxL0 = MaxL0, MinL1 = MinL1, MaxL1 = MaxL1)
  fwrite(res, paste0("./analysis/01.AASignalRecognition/Version9/02.RangesL0L1/Group2/", meta1[i, file_name], ".txt"), row.names = F, sep = "\t", quote = F)
  rm(list = setdiff(ls(), "meta1")); gc()
}







