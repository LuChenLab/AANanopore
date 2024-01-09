library(shiny)
library(plotly)

AABlockade <- data.table(AA = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"), 
                         Blockade = c(0.14718, 0.17148, 0.16516, 0.21409, 0.18622, 0.24528, 0.1882, 0.12072, 0.24652, 0.20722, 0.1995, 0.16875, 0.19772, 0.22018, 0.2183, 0.13131, 0.16101, 0.22744, 0.21276, 0.19044))
AABlockade$amino_acid <- plyr::mapvalues(AABlockade$AA, names(Biostrings::AMINO_ACID_CODE), Biostrings::AMINO_ACID_CODE)

StandardAA <- function(SigFile) {
  dataset <- fread(SigFile)
  shinyApp(
    ui = fluidPage(
      fluidRow(style = "padding-bottom: 20px;",
               column(2, numericInput('size', 'Point size', value = .3, min = .1, max = 1, step = .1)),
               column(3, sliderInput('nbins', 'Number of bins', value = 100, min = 100, max = 1000, step = 50)),
               column(3, sliderInput('BaseMean', 'BaseMean region', value = c(90, 130), min = 80, max = max(dataset[, round(BaseMean)])))
      ),
      fluidRow(column(6, plotlyOutput('BaseMeanSelect1')), column(6, plotlyOutput('BaseMeanSelect2'))),
      fluidRow(style = "padding-bottom: 20px;",
               column(4, sliderInput("DeltaMean", "Max DeltaMean", value = 1, min = 0, max = 2, step = 0.1)),
               column(4, sliderInput("CurrentSD", "CurrentSD", value = 15, min = 0, max = 40, step = 1)), 
               column(4, sliderInput("StageSD", "StageSD", value = 5, min = 0, max = 20, step = 0.1))
      ),
      fluidRow(style = "padding-bottom: 20px;",
               column(4, sliderInput("L2Ratio", "L2Ratio", value = c(0, 1), min = 0, max = 1, step = 0.01)),
               column(4, sliderInput("SignalCurrentPercent", "SignalCurrentPercent", min = 0, max = 100, value = c(50, 100), step = 5)),
               column(4, sliderInput("DwellTime", "DwellTime (ms)", value = c(0, 20), min = 0, max = 50, step = .1)),
               column(3, sliderInput('Blockade', 'Blockade region', value = c(.1, .5), min = 0, max = 1))
      ),
      plotlyOutput('PointSelect'),
      verbatimTextOutput("sp"),
      actionButton("goButton", 'Save', class = "btn-success")
    ), 
    server = function(input, output, session) {
      selectedData <- reactive({
        dataset[BaseMean > input$BaseMean[1] & BaseMean < input$BaseMean[2]]
      })
      
      output$BaseMeanSelect1 = renderPlotly({
        ggplot(selectedData(), aes(x = BaseMean, y = Blockade)) + 
          geom_point(size = input$size) + 
          geom_hline(yintercept = AABlockade[amino_acid %in% selectedData()[, unique(A)], Blockade], color = "red") + 
          theme_minimal() -> p
        # plot_ly(data = selectedData(), x = ~BaseMean, y = ~Blockade) %>% plot_ly(y = AABlockade[, Blockade])
        ggplotly(p)
      })
      
      output$BaseMeanSelect2 = renderPlotly({
        ggplot(selectedData(), aes(x = BaseMean)) + 
          geom_histogram(bins = input$nbins) + 
          theme_minimal() -> p2
        ggplotly(p2) %>% layout(dragmode = "select")
      })
      
      selectedData2 <- reactive({
        d <- (event_data("plotly_selected"))
        if(is.null(d)) {
          dataset[DeltaMean < input$DeltaMean & 
                    CurrentSD < input$CurrentSD & 
                    StageSD < input$StageSD & 
                    Blockade > 0 & Blockade < 0.5 &
                    L2Ratio >= input$L2Ratio[1] & L2Ratio <= input$L2Ratio[2] & 
                    SignalCurrentPercent >= input$SignalCurrentPercent[1] & SignalCurrentPercent <= input$SignalCurrentPercent[2] & 
                    DwellTime > input$DwellTime[1] & DwellTime < input$DwellTime[2]]
        } else {
          dataset[DeltaMean < input$DeltaMean & 
                    CurrentSD < input$CurrentSD & 
                    StageSD < input$StageSD & 
                    Blockade > min(input$Blockade) & Blockade < max(input$Blockade) &
                    L2Ratio >= input$L2Ratio[1] & L2Ratio <= input$L2Ratio[2] & 
                    SignalCurrentPercent >= input$SignalCurrentPercent[1] & SignalCurrentPercent <= input$SignalCurrentPercent[2] & 
                    DwellTime > input$DwellTime[1] & DwellTime < input$DwellTime[2] & 
                    BaseMean > with(d, min(x)) & BaseMean < with(d, max(x))]
        }
      })
      output$PointSelect <- renderPlotly({
        ggplot(selectedData2(), aes(x = Blockade, y = DwellTime)) + 
          geom_point(size = input$size) + 
          geom_vline(xintercept = AABlockade[amino_acid %in% selectedData()[, unique(A)], Blockade], color = "red") + 
          theme_minimal() + 
          labs(title = selectedData()[, unique(A)]) -> p3
        ggplotly(p3)
      })
      
      output$sp <- renderPrint({
        cat("Number of selected points:", nrow(selectedData2()), "\n")
        if(input$goButton != 0) {
          subdata <- selectedData2()
          fwrite(subdata, gsub(".txt$", "_Selected.txt", SigFile), sep = "\t", quote = F)
          cat("Finished")
        }
      })
    }
  )
}

files <- list.files("./analysis/11.SignalIdentification/Jan07/BackgroundSignal", full.names = TRUE)

i = 74
StandardAA(files[i])




