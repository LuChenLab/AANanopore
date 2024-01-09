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

SignalCurrent <- function(x, abf, cores = 10) {
  Sigs_Multiple_Current <- mclapply(seq_len(nrow(x)), FUN = function(i) {
    data.table(ID = x[i, ID], Current = abf[Time > x[i, StartTime] & Time < x[i, EndTime], pA] / x[i, BaseMean])
  }, mc.cores = cores)
  do.call(rbind, Sigs_Multiple_Current)
}

meta <- data.table(openxlsx::read.xlsx("./data/ChenShanchuan/20230521/副本副本数据记录20230521新.xlsx", sheet = 2))
meta <- meta[实验目的 == "阿尔兹海默症突变位点"]
colnames(meta)[1:4] <- c("file_name", "start_time", "end_time", "amino_acid")
meta <- meta[!grepl("blank", amino_acid)]
files <- list.files("./analysis/11.SignalIdentification", "RawSignal", full.names = T, recursive = T)
meta$file_path <- as.character(mapply(meta[, paste0(file_name, "_", start_time, "_", end_time)], FUN = function(x) grep(pattern = x, x = files, value = T)))
meta[file_path == "character(0)", file_path := NA]
meta <- meta[file_path != "character(0)"]
meta[, start_time := as.numeric(start_time)]
meta[, end_time := as.numeric(end_time)]
meta <- meta[amino_acid == "LVFAG"]

mclapply(1:nrow(meta), function(i) {
  print(i)
  abf <- readRDS(paste0("./analysis/11.SignalIdentification/Apr02/ABF_", meta[amino_acid == "LVFAG", paste(amino_acid, file_name, start_time, end_time, sep = "_")], ".Rds"))
  Sigs <- GetSignal(x = abf)
  Sigs[, DwellTime := EndTime - StartTime]
  Sigs <- data.table(ID = paste(meta[i, file_name], Sigs[, seq_len(.N)], sep = "_"), Sigs)
  fwrite(Sigs, file = paste0("./analysis/11.SignalIdentification/Apr02/RawSignal_", meta[amino_acid == "LVFAG", paste(amino_acid, file_name, start_time, end_time, sep = "_")], "_V2.txt"), sep = "\t", row.names = F, quote = F)
}, mc.cores = 1)



library(shiny)
library(plotly)

ui <- fluidPage(
  plotlyOutput('myPlot1'),
  verbatimTextOutput("sp"),
  plotlyOutput('myPlot2')
)

server <- function(input, output){
  output$myPlot1 = renderPlotly({
    plot_ly(data = Sigs, x = ~ BaseMean, type = "histogram") %>%
      layout(dragmode = "select")
  })
  
  output$sp <- renderPrint({
    d <- (event_data("plotly_selected"))
    if(!is.null(d)) {
      cat("Range of selected:", with(d, paste(round(range(x), 3), collapse = " - ")))
    }
  })
  
  output$myPlot2 = renderPlotly({
    d <- (event_data("plotly_selected"))
    if(!is.null(d)) {
      fwrite(d, "/mnt/raid61/Personal_data/tangchao/Temp/subset.points.txt", sep = "\t")
      ggplot(data = Sigs[BaseMean > with(d, min(x)) & BaseMean < with(d, max(x))], aes(x = Blockade, colour = cut_interval(SignalCurrentPercent, n = 10))) + 
        geom_line(stat = "density", bw = .003, n = 10000, na.rm = T) + 
        scale_colour_manual(values = colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "Reds"))(10), guide = "none") + 
        geom_vline(xintercept = AABlockade[AA %in% unlist(strsplit("LVFGKA", "")), Blockade]) + 
        theme(legend.position = "none") -> p
      ggplotly(p)
    }
  })
}


Sigs <- fread(file = paste0("./analysis/11.SignalIdentification/Apr02/RawSignal_LVFAG_20230330 _0006_0_23_V2.txt"))
Sigs[, Blockade := 1 - SignalCurrent / BaseMean]
Sigs[, Blockade2 := 1 - SignalCurrent.5 / BaseMean]
Sigs <- Sigs[Valid == TRUE & BaseMean > 90 & BaseMean < 130]

shinyApp(ui, server)
L0 <- range(fread("/mnt/raid61/Personal_data/tangchao/Temp/subset.points.txt")[[2]])
Sigs1 <- Sigs[BaseMean > min(L0) & BaseMean < max(L0)]
Sigs1 <- Sigs1[ID %in% unique(Current$ID)]


SignalCurrent <- function(x, abf, cores = 10) {
  Sigs_Multiple_Current <- mclapply(seq_len(nrow(x)), FUN = function(i) {
    data.table(ID = x[i, ID], Current = abf[Time > x[i, StartTime] & Time < x[i, EndTime], pA] / x[i, BaseMean])
  }, mc.cores = cores)
  do.call(rbind, Sigs_Multiple_Current)
}


abf <- readRDS(paste0("./analysis/11.SignalIdentification/Apr02/ABF_LVFAG_20230330 _0006_0_23.Rds"))
Current <- SignalCurrent(Sigs1, abf = abf, cores = 10)

Feature_Mat <- lapply(Sigs1[, ID], function(j) {
  print(j)
  if(Current[ID == j, .N] == 0) return(NULL)
  D <- density(Current[ID == j, Current], from = 0, to = 1, n = 500, adjust = 0.5)$y
  pA <- sort(Current[ID == j, Current])
  pA <- pA[round(seq(1, length(pA), length.out = 100))]
  D <- c(D, pA)
  cbind(t(data.frame(round(D, 4), row.names = c(paste0("X", sprintf("%03d", 1:500)), paste0("P", sprintf("%03d", 1:100))))), 
        Sigs1[ID == j, .(DeltaMean, StageSD, CurrentSD, Segments, SignalCurrentPercent, L2Ratio, DwellTime, Blockade, Blockade2, ID)])
})

Feature_Mat <- do.call(rbind, Feature_Mat)
save(Current, Feature_Mat, file = "./analysis/13.MachineLearning/01.DataPreparation/Apr02/ABF_LVFAG_20230330 _0006_0_23.RData")

Fit_L1_X <- readRDS("./analysis/13.MachineLearning/01.DataPreparation/Jan10/Model0/RF/RF_Fit_L1_X.Rds")

pred <- predict(Feature_Mat, object = Fit_L1_X, type = "prob")
Sigs1$Pred <- predict(Feature_Mat, object = Fit_L1_X)
Sigs1[, Pred := plyr::mapvalues(Pred, AABlockade[, amino_acid], AABlockade[, AA])]
Sigs1$Prob <- as.numeric(apply(pred, 1, max))

ggplot(Sigs1[Blockade < 0.3 & Pred != "Noise"], aes(x = Blockade, y = DwellTime * 1000)) + 
  geom_point() + 
  scale_y_log10() + 
  geom_vline(xintercept = AABlockade[AA %in% unlist(strsplit("LVFGKA", "")), Blockade])


ggplot(Sigs1[Blockade < 0.3 & Pred != "Noise"], aes(x = Blockade, y = DwellTime * 1000, label = Pred, colour = Prob)) + 
  geom_text() + 
  scale_y_log10()

M2 <- Sigs1[Blockade < 0.3 & Pred != "Noise"][, .N, Pred]
ggplot(M2, aes(x = Pred, y = N, fill = !Pred %in% unlist(strsplit("LVFGKA", "")))) + 
  geom_col() + 
  labs(title = "LVFGKA") + 
  theme(legend.position = "none")

Sigs1[Blockade < 0.3 & Pred != "Noise" & Prob > .5, mean(Pred %in% unlist(strsplit("LVFGKA", "")))]

