setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(shiny)
library(plotly)
library(data.table)
set.seed(19)
AA_Cols <- sample(c(RColorBrewer::brewer.pal(n = 12, "Paired"), RColorBrewer::brewer.pal(n = 8, "Dark2")), 20)
names(AA_Cols) <- Biostrings::AA_STANDARD
AA_Cols[AA_Cols == "#FFFF99"] <- "#1A1A1A"
# AA_Cols <- c(AA_Cols, "black")

AABlockade <- data.table(AA = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"), 
                         Blockade = c(0.14718, 0.17148, 0.16516, 0.21409, 0.18622, 0.24528, 0.1882, 0.12072, 0.24652, 0.20722, 0.1995, 0.16875, 0.19772, 0.22018, 0.2183, 0.13131, 0.16101, 0.22744, 0.21276, 0.19044))
AABlockade$amino_acid <- plyr::mapvalues(AABlockade$AA, names(Biostrings::AMINO_ACID_CODE), Biostrings::AMINO_ACID_CODE)

PointStat <- function(dataset, target = NULL, outfile = NULL) {
  if(is.null(target)) {
    As <- NULL
  } else {
    As <- intersect(Biostrings::AA_STANDARD, unlist(strsplit(target, "")))
  }
  
  shinyApp(
    ui = fluidPage(
      fluidRow(style = "padding-bottom: 20px;",
               column(1, numericInput('TimeStart', 'TimeStart', value = 0, min = 0, max = dataset[, ceiling(max(StartTime)/60)], step = 1)),
               column(1, numericInput('TimeEnd', 'TimeEnd', value = dataset[, ceiling(max(EndTime)/60)], min = 0, max = dataset[, ceiling(max(EndTime)/60)], step = 1)),
               column(1, numericInput('binwidth', 'binwidth', value = .1, min = 0, max = 1, step = .1)),
               column(4, sliderInput('BaseMean', 'BaseMean region', value = c(90, 130), min = 80, max = max(dataset[, round(BaseMean)])))
      ),
      plotlyOutput('BaseMeanSelect'),
      fluidRow(style = "padding-bottom: 20px;",
               column(4, sliderInput("DeltaMean", "Max DeltaMean", value = 1, min = 0, max = 2, step = 0.1)),
               column(4, sliderInput("CurrentSD", "CurrentSD", value = 15, min = 0, max = 40, step = 1)), 
               column(4, sliderInput("StageSD", "StageSD", value = 3.5, min = 0, max = 20, step = 0.1))
      ),
      fluidRow(style = "padding-bottom: 20px;",
               column(4, sliderInput("SignalCurrentWidth", "SignalCurrentWidth", value = c(5, 15), min = 0, max = ceiling(dataset[, max(SignalCurrentWidth)]), step = 0.1)),
               column(4, sliderInput("SignalCurrentPercent", "SignalCurrentPercent", min = 0, max = 100, value = c(80, 100), step = 5)),
               column(4, sliderInput("DwellTime", "DwellTime (ms)", value = c(0, 20), min = 0, max = 50, step = .1))
      ),
      fluidRow(style = "padding-bottom: 20px;",
               column(4, sliderInput("Blockade", "Blockade", value = c(.1, .3), min = 0, max = 1, step = .01)), 
               column(4, selectInput('Model', 'Model', choices = c('M1', 'M2', 'M3'))), 
               column(4, sliderInput("prob", "Probability", value = c(0, 1), min = 0, max = 1, step = .01))
      ),
      fluidRow(style = "padding-bottom: 20px;",
               column(width = 4, numericInput('size', 'size', value = 3, min = 0, max = 10, step = 1)),
               column(width = 4, numericInput('adjust', 'adjust', value = 0.5, min = 0.1, max = Inf, step = 0.1)),
               column(width = 4, sliderInput('binwidth2', 'binwidth', value = 0.001, min = 1e-03, max = .1, step = 1e-03))
      ),
      fluidRow(column(6, plotlyOutput('PointSelect')), column(6, plotlyOutput('hist'))),
      fluidRow(column(6, plotlyOutput('density')), column(6, plotlyOutput('barchart'))),
      actionButton("goButton", 'Save', class = "btn-success"), 
      verbatimTextOutput("value")
    ), 
    server = function(input, output, session) {
      selectedData <- reactive({
        dataset[StartTime >= input$TimeStart * 60 & EndTime <= input$TimeEnd * 60 & BaseMean > input$BaseMean[1] & BaseMean < input$BaseMean[2]]
      })
      
      output$BaseMeanSelect = renderPlotly({
        ggplot(selectedData(), aes(x = BaseMean)) + 
          geom_histogram(binwidth = input$binwidth) + 
          theme_minimal() -> p1
        ggplotly(p1) %>% layout(dragmode = "select")
      })
      
      selectedData2 <- reactive({
        d <- (event_data("plotly_selected"))
        if(is.null(d)) {
          dataset2 <- dataset[StartTime >= input$TimeStart * 60 & EndTime <= input$TimeEnd * 60 & 
                                DeltaMean < input$DeltaMean & 
                                CurrentSD < input$CurrentSD & 
                                StageSD < input$StageSD & 
                                Blockade > input$Blockade[1] & Blockade < input$Blockade[2] &
                                SignalCurrentWidth >= input$SignalCurrentWidth[1] & SignalCurrentWidth <= input$SignalCurrentWidth[2] & 
                                SignalCurrentPercent >= input$SignalCurrentPercent[1] & SignalCurrentPercent <= input$SignalCurrentPercent[2] & 
                                DwellTime > input$DwellTime[1] / 1000 & DwellTime < input$DwellTime[2] / 1000 & 
                                BaseMean > min(input$BaseMean) & BaseMean < max(input$BaseMean)]
        } else {
          dataset2 <- dataset[StartTime >= input$TimeStart * 60 & EndTime <= input$TimeEnd * 60 & 
                                DeltaMean < input$DeltaMean & 
                                CurrentSD < input$CurrentSD & 
                                StageSD < input$StageSD & 
                                Blockade > input$Blockade[1] & Blockade < input$Blockade[2] &
                                SignalCurrentWidth >= input$SignalCurrentWidth[1] & SignalCurrentWidth <= input$SignalCurrentWidth[2] & 
                                SignalCurrentPercent >= input$SignalCurrentPercent[1] & SignalCurrentPercent <= input$SignalCurrentPercent[2] & 
                                DwellTime > input$DwellTime[1] / 1000 & DwellTime < input$DwellTime[2] / 1000 & 
                                BaseMean > with(d, min(x)) & BaseMean < with(d, max(x))]
        }
        
        setnames(dataset2, paste0(input$Model, "Pred"), "Pred")
        setnames(dataset2, paste0(input$Model, "Prob"), "Prob")
        dataset2[Prob >= min(input$prob) & Prob <= max(input$prob) & Pred != "Noise"]
      })
      
      output$PointSelect = renderPlotly({
        ggplot(selectedData2(), aes(x = Blockade, y = DwellTime * 1000, colour = Pred, label = Pred)) + 
          geom_text(size = input$size) + 
          scale_color_manual(values = AA_Cols) +
          scale_y_log10() + 
          theme_minimal() + 
          # theme(legend.position = "none") + 
          labs(y = "DwellTime (ms)") -> p2
        if(!is.null(As)) {
          p2 <- p2 + geom_vline(xintercept = AABlockade[AA %in% As, Blockade], colour = "grey")
        }
        ggplotly(p2)
      })
      
      output$hist = renderPlotly({
        ggplot(selectedData2(), aes(x = Blockade, fill = Pred)) + 
          geom_histogram(binwidth = input$binwidth2) + 
          scale_fill_manual(values = AA_Cols) + 
          theme_minimal() -> p3
        ggplotly(p3)
      })
      
      output$density = renderPlotly({
        ggplot(selectedData2(), aes(x = Blockade)) + 
          geom_line(stat = "density", adjust = input$adjust) + 
          theme_minimal() -> p4
        if(!is.null(As)) {
          p4 <- p4 + geom_vline(xintercept = AABlockade[AA %in% As, Blockade], colour = "grey")
        }
        ggplotly(p4)
      })
      
      output$barchart = renderPlotly({
        tab <- selectedData2()[, .N, Pred]
        tab <- rbind(tab, data.table(Pred = Biostrings::AA_STANDARD, N = 0))[, .SD[which.max(N), ], Pred]
        tab[, Pred := factor(Pred, levels = tab[order(N, decreasing = T), Pred])]
        
        ggplot(tab, aes(x = Pred, y = N)) + 
          geom_col(aes(fill = !Pred %in% As)) + 
          geom_text(aes(label = N, colour = !Pred %in% As), position = position_nudge(y = tab[, max(N)*.02])) + 
          theme_classic() + 
          theme(legend.position = "none", 
                axis.title = element_blank(), 
                axis.line.y = element_blank(), 
                axis.text.y = element_blank(), axis.ticks.y = element_blank()) -> p5
        ggplotly(p5)
      })
      
      output$value <- renderPrint({
        if(is.null(outfile)) {
          cat(paste0("Target AAs (%): ", selectedData2()[, round(mean(Pred %in% As) * 100, 2)]))
        } else {
          if(input$goButton == 0) {
            cat(paste0("Target AAs (%): ", selectedData2()[, round(mean(Pred %in% As) * 100, 2)]))
          } else {
            fwrite(selectedData2(), file = outfile, sep = "\t", quote = F)
            cat("Finished")
          }
        }
      })
    }
  )
}


RFC <- readRDS("./analysis/43.ModelTraining/Model1/RF4.Rds")

fm <- fread("analysis/21.ABFProcessing/03.MixtureAA/FeatureMatrix/FeatureMatrix_20230417_0005.txt")
fm <- data.table(ID = fm[, ID], M1Pred = predict(RFC, fm), M1Prob = apply(predict(RFC, fm, type = "prob"), 1, max))
dataset <- fread("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_20230417_0005.txt")
dataset <- merge(dataset, fm, by = "ID")
PointStat(dataset, target = "LWQSTIFN")
setkey(dataset, StartTime)

abf <- readRDS("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_20230417_0005.Rds")
sort(table(dataset[M1Pred != "Noise", round(StartTime)]))
ggplot(abf[Time > 505 - 5 & Time < 505 + 5], aes(x = Time, y = pA)) + 
  geom_step(colour = 'grey80') + 
  geom_text(data = dataset[StartTime > 505 - 5 & EndTime < 505 + 5 & M1Prob > .2 & M1Pred != "Noise"], 
            aes(x = StartTime + DwellTime / 2, y = SignalCurrent - 6, label = M1Pred), colour = 'red') + 
  theme_light(base_size = 15) + 
  theme(panel.grid = element_blank())




SigTime <- dataset[StartTime > 505 - 5 & EndTime < 505 + 5 & M1Prob > .2 & M1Pred != "Noise", ]
abf1 <- lapply(seq_len(nrow(SigTime)), function(i) {
  print(i)
  Ci <- abf[Time >= SigTime$StartTime[i] - min(0.001, 1e-05 * SigTime$LeftLength[i]) & Time <= SigTime$EndTime[i] + min(0.001, 1e-05 * SigTime$RightLength[i])]
  Ci$ID <- SigTime$ID[i]
  Ci$M1Pred <- SigTime$M1Pred[i]
  Ci
})
abf1 <- do.call(rbind, abf1)

ggplot(abf1, aes(x = Time, y = pA, colour = M1Pred)) + 
  geom_step() + 
  facet_wrap(~ID, nrow = 1, scales = "free_x") + 
  theme_minimal() + 
  theme(strip.text = element_blank(),
        panel.grid = element_blank(), 
        axis.text.x = element_blank(),
        axis.line.x = element_blank()) + 
  scale_colour_manual(values = AA_Cols) + 
  geom_text(data = SigTime, aes(x = StartTime + DwellTime / 2, y = SignalCurrent - 6, label = M1Pred), colour = 'red')
  





fm <- fread("analysis/21.ABFProcessing/03.MixtureAA/FeatureMatrix/FeatureMatrix_20230419_0002.txt")
fm <- data.table(ID = fm[, ID], M1Pred = predict(RFC, fm), M1Prob = apply(predict(RFC, fm, type = "prob"), 1, max))
dataset <- fread("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_20230419_0002.txt")
dataset <- merge(dataset, fm, by = "ID")
PointStat(dataset, target = "DE")
setkey(dataset, StartTime)

ggplot(dataset, aes(x = M1Pred, y = M1Prob)) + 
  geom_boxplot() + 
  geom_jitter(height = 0)

abf <- readRDS("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_20230419_0002.Rds")
sort(table(dataset[M1Pred != "Noise", round(StartTime)]))
ggplot(abf[Time > 729 - .5 & Time < 729 + .5], aes(x = Time, y = pA)) + 
  geom_step(colour = 'grey80') + 
  geom_text(data = dataset[StartTime > 729 - .5 & EndTime < 729 + .5 & M1Prob > .2], 
            aes(x = StartTime + DwellTime / 2, y = SignalCurrent - 6, label = M1Pred), colour = 'red') + 
  theme_light(base_size = 15) + 
  theme(panel.grid = element_blank())





fm <- fread("analysis/21.ABFProcessing/03.MixtureAA/FeatureMatrix/FeatureMatrix_20230512_0004.txt")
fm <- data.table(ID = fm[, ID], M1Pred = predict(RFC, fm), M1Prob = apply(predict(RFC, fm, type = "prob"), 1, max))
dataset <- fread("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_20230512_0004.txt")
dataset <- merge(dataset, fm, by = "ID")
PointStat(dataset, target = "GAVLIFWYNQMSTPDE")
setkey(dataset, StartTime)

ggplot(dataset, aes(x = M1Pred, y = M1Prob)) + 
  geom_boxplot() + 
  geom_jitter(height = 0)


abf <- readRDS("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_20230512_0004.Rds")
sort(table(dataset[M1Pred != "Noise", round(StartTime)]))
ggplot(abf[Time > 701 - .5 & Time < 701 + .5], aes(x = Time, y = pA)) + 
  geom_step(colour = 'grey80') + 
  geom_text(data = dataset[StartTime > 701 - .5 & EndTime < 701 + .5], 
            aes(x = StartTime + DwellTime / 2, y = SignalCurrent - 6, label = M1Pred), colour = 'red') + 
  theme_light(base_size = 15) + 
  theme(panel.grid = element_blank())

dataset2 <- dataset[CurrentSD < 5 & M1Pred != 'Noise']
ggplot(abf[Time > 701 - .5 & Time < 701 + .5], aes(x = Time, y = pA)) + 
  geom_step(colour = 'grey80') + 
  geom_text(data = dataset2[StartTime > 701 - .5 & EndTime < 701 + .5], 
            aes(x = StartTime + DwellTime / 2, y = SignalCurrent - 6, label = M1Pred), colour = 'red') + 
  theme_light(base_size = 15) + 
  theme(panel.grid = element_blank())








PredFun <- function(mat, model, p = .9) {
  res <- data.table(ID = mat[, ID], Pred = predict(model, mat), Prob = apply(predict(model, mat, type = "prob"), 1, max))
  res0 <- res[!Pred %in% c('T', 'N', 'R', 'M', 'L', 'F', 'D', 'P')]
  Id1 <- which(res[, Pred %in% c('T', 'N', 'R') & Prob > p])
  Id2 <- which(res[, Pred %in% c('M', 'L') & Prob > p])
  Id3 <- which(res[, Pred %in% c('F', 'D', 'P') & Prob > p])
  if(length(Id1) > 0) {
    res1 <- data.table(ID = mat[Id1, ID], Pred = predict(RFC3.1, mat[Id1, ]), Prob = apply(predict(RFC3.1, mat[Id1, ], type = "prob"), 1, max))
    res0 <- rbind(res0, res1)
  }
  
  if(length(Id2) > 0) {
    res2 <- data.table(ID = mat[Id2, ID], Pred = predict(RFC3.2, mat[Id2, ]), Prob = apply(predict(RFC3.2, mat[Id2, ], type = "prob"), 1, max))
    res0 <- rbind(res0, res2)
  }
  
  if(length(Id3) > 0) {
    res3 <- data.table(ID = mat[Id3, ID], Pred = predict(RFC3.3, mat[Id3, ]), Prob = apply(predict(RFC3.3, mat[Id3, ], type = "prob"), 1, max))
    res0 <- rbind(res0, res3)
  }
  return(res0)
}
RFC3.1 <- readRDS("./analysis/43.ModelTraining/Model1/RF3.1.Rds")
RFC3.2 <- readRDS("./analysis/43.ModelTraining/Model1/RF3.2.Rds")
RFC3.3 <- readRDS("./analysis/43.ModelTraining/Model1/RF3.3.Rds")


fm <- fread("analysis/21.ABFProcessing/03.MixtureAA/FeatureMatrix/FeatureMatrix_21o19001.txt")
fm <- data.table(ID = fm[, ID], M1Pred = predict(RFC, fm), M1Prob = apply(predict(RFC, fm, type = "prob"), 1, max))

fm <- PredFun(mat = fm, model = RFC, p = .2); colnames(fm) <- c("ID", "M1Pred", "M1Prob")

dataset <- fread("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_21o19001.txt")
dataset <- merge(dataset, fm, by = "ID")
PointStat(dataset, target = "FYSL")
setkey(dataset, StartTime)

ggplot(dataset, aes(x = M1Pred, y = M1Prob)) + 
  geom_boxplot() + 
  geom_jitter(height = 0)


abf <- readRDS("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_21o19001.Rds")
sort(table(dataset[M1Pred != "Noise", round(StartTime)]))
ggplot(abf[Time > 701 - .5 & Time < 701 + .5], aes(x = Time, y = pA)) + 
  geom_step(colour = 'grey80') + 
  geom_text(data = dataset[StartTime > 701 - .5 & EndTime < 701 + .5], 
            aes(x = StartTime + DwellTime / 2, y = SignalCurrent - 6, label = M1Pred), colour = 'red') + 
  theme_light(base_size = 15) + 
  theme(panel.grid = element_blank())

dataset2 <- dataset[CurrentSD < 5 & M1Pred != 'Noise']
ggplot(abf[Time > 701 - .5 & Time < 701 + .5], aes(x = Time, y = pA)) + 
  geom_step(colour = 'grey80') + 
  geom_text(data = dataset2[StartTime > 701 - .5 & EndTime < 701 + .5], 
            aes(x = StartTime + DwellTime / 2, y = SignalCurrent - 6, label = M1Pred), colour = 'red') + 
  theme_light(base_size = 15) + 
  theme(panel.grid = element_blank())


