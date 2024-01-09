setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(parallel)
library(Biostrings)

files <- list.files("./analysis/11.SignalIdentification/Jan07/AASignal", pattern = "_Feature_Mat.txt", full.names = TRUE)
Feature_Mat <- lapply(files, fread)
Feature_Mat <- do.call(rbind, Feature_Mat)

files <- list.files("./analysis/11.SignalIdentification/Jan07/AASignal", pattern = "Selected_V2", full.names = TRUE)
Selected_Sig <- lapply(files, fread)
Selected_Sig <- do.call(rbind, Selected_Sig)
mean(Selected_Sig[, ID] %in% Feature_Mat[, ID])

files2 <- list.files("./analysis/11.SignalIdentification/Jan07/BackgroundSignal", pattern = "_Feature_Mat.txt", full.names = TRUE)
Feature_Mat0 <- lapply(files2, fread)
Feature_Mat0 <- do.call(rbind, Feature_Mat0)

files <- list.files("./analysis/11.SignalIdentification/Jan07/BackgroundSignal", pattern = "Selected", full.names = TRUE)
Selected_Sig0 <- lapply(files, fread)
Selected_Sig0 <- do.call(rbind, Selected_Sig0)
Selected_Sig0[, A := "Noise"]

Selected_Sig0 <- Selected_Sig0[Blockade <= Selected_Sig[, max(Blockade)]]

Feature_Mat1 <- Feature_Mat[ID %in% Selected_Sig[, ID]]
Feature_Mat2 <- Feature_Mat0[ID %in% Selected_Sig0[, ID]]
Feature_Mat2[, A := "Noise"]


Mat <- rbind(Feature_Mat1, Feature_Mat2)

set.seed(3456)
Train <- rbind(Mat[A %in% Mat[, .N, A][N >= 2000, A], .SD[sample(.N, 2000), ], A], 
               Mat[A %in% Mat[, .N, A][N < 2000, A], .SD[sample(.N, 2000, replace = T), ], A])
Train <- Train[, c("ID", grep("^X", colnames(Train), value = T), "A"), with = F]
Train[, A := plyr::mapvalues(A, Biostrings::AMINO_ACID_CODE, names(Biostrings::AMINO_ACID_CODE))]

Train <- data.frame(Train[, 2:501], Class = factor(Train[[502]]))
Train <- Train[, 201:501]

saveRDS(Train, "./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model2/Train.Rds")



library(data.table)
library(Biostrings)
library(ggplot2)
library(parallel)
library(caret)

fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10)

library(doParallel)
cl <- makePSOCKcluster(10)
registerDoParallel(cl)


Train <- readRDS("./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model2/Train.Rds")
set.seed(825)
RFC <- train(Class ~ ., data = Train, 
             # preProc = c("center", "scale", "YeoJohnson", "nzv"), 
             method = "rf", 
             trControl = fitControl,
             verbose = FALSE,
             ## to evaluate:
             tuneGrid = expand.grid(mtry = 30),
             # tuneLength = 50,
             metric = "Accuracy", 
             allowParallel = TRUE)
saveRDS(RFC, file = "./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model2/RF_Fit.Rds")


RFC <- readRDS("./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model2/RF_Fit.Rds")


files <- list.files("./analysis/11.SignalIdentification/Jan07/AASignal", pattern = "_Feature_Mat.txt", full.names = TRUE)
Feature_Mat <- lapply(files, fread)
Feature_Mat <- do.call(rbind, Feature_Mat)
Feature_Mat$File <- mapply(function(x) x[1], strsplit(Feature_Mat$ID, "_"))
files <- list.files("./analysis/11.SignalIdentification/Jan07/AASignal", pattern = "Selected_V2", full.names = TRUE)

Feature_Mat <- Feature_Mat[File %in% gsub("RawSignal_", "", gsub("_Selected_V2.txt", "", basename(files)))]

Feature_Mat$Pred <- predict(RFC, Feature_Mat)
Feature_Mat$Prob <- apply(predict(RFC, Feature_Mat, type = "prob"), 1, max)
Feature_Mat <- Feature_Mat[, .(ID, DeltaMean, StageSD, CurrentSD, Segments, SignalCurrentPercent, L2Ratio, DwellTime, Blockade, A, Pred, Prob)]
Feature_Mat[, A := plyr::mapvalues(A, AMINO_ACID_CODE, names(AMINO_ACID_CODE))]

hist(Feature_Mat[, Prob], breaks = 100)
Feature_Mat[, mean(Prob == 1)]
Feature_Mat[Prob == 1, .N, A]
Feature_Mat[Prob == 1, .N, Pred]

Feature_Mat[Prob > .9, .N, Pred]
Feature_Mat[Blockade < 0.3 & Pred != "Noise" & Prob > .9, .N, Pred]





ggplot(Feature_Mat[Blockade < 0.3, ], aes(x = Blockade, y = DwellTime, label = Pred, colour = Pred)) + 
  geom_text() + 
  scale_y_log10()

ggplot(Feature_Mat[Blockade < 0.3 & Pred != "Noise" & Prob > .95, ], aes(x = Blockade, y = DwellTime, label = Pred, colour = Pred)) + 
  geom_text() + 
  scale_y_log10()


ggplot(Feature_Mat1[Blockade < 0.3 & Pred != "Noise" & Prob > .5, ], aes(x = Blockade, y = DwellTime, label = A, colour = A)) + 
  geom_text() + 
  scale_y_log10()


merge(Feature_Mat1[, median(DwellTime), A], by.x = "A", by.y = "Pred", 
Feature_Mat[Blockade < 0.3 & Pred != "Noise" & Prob > .8, median(DwellTime), Pred])



fm <- fread("./analysis/11.SignalIdentification/FeatureMatrix/FeatureMatrix_20230330 _0006.txt")
fm <- data.table(ID = fm[, ID], M1Pred = predict(RFC, fm), M1Prob = apply(predict(RFC, fm, type = "prob"), 1, max))
dataset <- fread("./analysis/11.SignalIdentification/Signal/RawSignal_20230330 _0006.txt")
dataset <- merge(dataset, fm, by = "ID")
PointStat(dataset, target = "LVFAG")

fm <- fread("./analysis/11.SignalIdentification/FeatureMatrix/FeatureMatrix_20230330 _0007.txt")
fm <- data.table(ID = fm[, ID], M1Pred = predict(RFC, fm), M1Prob = apply(predict(RFC, fm, type = "prob"), 1, max))
dataset <- fread("./analysis/11.SignalIdentification/Signal/RawSignal_20230330 _0007.txt")
dataset <- merge(dataset, fm, by = "ID")
PointStat(dataset, target = "LVFGKA")

fm <- fread("./analysis/11.SignalIdentification/FeatureMatrix/FeatureMatrix_20230417_0007.txt")
fm <- data.table(ID = fm[, ID], M1Pred = predict(RFC, fm), M1Prob = apply(predict(RFC, fm, type = "prob"), 1, max))
dataset <- fread("./analysis/11.SignalIdentification/Signal/RawSignal_20230417_0007.txt")
dataset <- merge(dataset, fm, by = "ID")
PointStat(dataset, target = "LWQSTIFNDE")

fm <- fread("./analysis/11.SignalIdentification/FeatureMatrix/FeatureMatrix_20230508_0014.txt")
fm <- data.table(ID = fm[, ID], M1Pred = predict(RFC, fm), M1Prob = apply(predict(RFC, fm, type = "prob"), 1, max))
dataset <- fread("./analysis/11.SignalIdentification/Signal/RawSignal_20230508_0014.txt")
dataset <- merge(dataset, fm, by = "ID")
PointStat(dataset, target = "GAVLIFWYNQMSTPDE")

fm <- fread("./analysis/11.SignalIdentification/FeatureMatrix/FeatureMatrix_20230512_0004.txt")
fm <- data.table(ID = fm[, ID], M1Pred = predict(RFC, fm), M1Prob = apply(predict(RFC, fm, type = "prob"), 1, max))
dataset <- fread("./analysis/11.SignalIdentification/Signal/RawSignal_20230512_0004.txt")
dataset <- merge(dataset, fm, by = "ID")
PointStat(dataset, target = "GAVLIFWYNQMSTPDE")
#
library(shiny)
library(plotly)

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
               column(4, sliderInput("StageSD", "StageSD", value = 5, min = 0, max = 20, step = 0.1))
      ),
      fluidRow(style = "padding-bottom: 20px;",
               column(4, sliderInput("L2Ratio", "L2Ratio", value = c(0, 1), min = 0, max = 1, step = 0.01)),
               column(4, sliderInput("SignalCurrentPercent", "SignalCurrentPercent", min = 0, max = 100, value = c(50, 100), step = 5)),
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
               column(width = 4, sliderInput('binwidth2', 'binwidth', value = 0.01, min = 1e-03, max = .1, step = 1e-03))
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
                                L2Ratio >= input$L2Ratio[1] & L2Ratio <= input$L2Ratio[2] & 
                                SignalCurrentPercent >= input$SignalCurrentPercent[1] & SignalCurrentPercent <= input$SignalCurrentPercent[2] & 
                                DwellTime > input$DwellTime[1] / 1000 & DwellTime < input$DwellTime[2] / 1000 & 
                                BaseMean > min(input$BaseMean) & BaseMean < max(input$BaseMean)]
        } else {
          dataset2 <- dataset[StartTime >= input$TimeStart * 60 & EndTime <= input$TimeEnd * 60 & 
                                DeltaMean < input$DeltaMean & 
                                CurrentSD < input$CurrentSD & 
                                StageSD < input$StageSD & 
                                Blockade > input$Blockade[1] & Blockade < input$Blockade[2] &
                                L2Ratio >= input$L2Ratio[1] & L2Ratio <= input$L2Ratio[2] & 
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
        tab[, Pred := factor(Pred, levels = tab[order(N, decreasing = T), Pred])]
        
        ggplot(tab, aes(x = Pred, y = N)) + 
          geom_col(aes(fill = !Pred %in% As)) + 
          geom_text(aes(label = N), position = position_nudge(y = tab[, max(N)*.01])) + 
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
