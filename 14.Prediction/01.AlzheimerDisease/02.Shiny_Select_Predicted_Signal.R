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


dataset <- fread("./analysis/14.Prediction/01.ShinySelected/Model1/SignalPredicted_20230417_0007.txt")
PointStat(dataset, target = "LWQSTIFNDE")









