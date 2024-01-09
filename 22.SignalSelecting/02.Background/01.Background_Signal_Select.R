setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")

library(data.table)
meta1 <- data.table(openxlsx::read.xlsx("./data/单氨基酸实验信息表20210517.xlsx", sheet = 1, cols = 1:7))
meta2 <- data.table(openxlsx::read.xlsx("./data/单氨基酸实验信息表20210517.xlsx", sheet = 2, cols = 1:7))
meta3 <- data.table(openxlsx::read.xlsx("./data/单氨基酸实验信息表20210517.xlsx", sheet = 3, cols = 1:7))
meta4 <- data.table(openxlsx::read.xlsx("./data/单氨基酸实验信息表20210517.xlsx", sheet = 4, cols = 1:7))
meta <- rbind(meta1, meta2, meta3, meta4, use.names = FALSE)
colnames(meta) <- c("file_name", "date", "amino_acid", "concentration", "start_time", "end_time", "type")
meta[amino_acid == "cys", amino_acid := "Cys"]
meta <- meta[concentration == 0 | amino_acid == "blank"]
meta$file_path <- mapply(meta$file_name, FUN = function(x) list.files("./data", recursive = TRUE, full.names = TRUE, pattern = as.character(x))[1])
meta <- meta[!is.na(file_path)]
meta[, start_time := as.numeric(start_time)]
meta[, end_time := as.numeric(end_time)]
# meta <- meta[!grepl("abf", file_name)]
meta[, file_name := gsub(".abf$", "", file_name)]


StandardAA <- function(SigFile, outdir = "analysis/22.SignalSelecting/02.Background", target = NULL) {
  if(is.null(outdir)) stop()
  if(is.null(target)) stop()
  outdir <- gsub("//", "/", paste0(outdir, "/", target))
  if(!dir.exists(outdir)) dir.create(outdir)
  dataset <- fread(SigFile)
  dataset[, A := target]
  dataset[, DwellTime := DwellTime * 1000]
  shinyApp(
    ui = fluidPage(
      fluidRow(style = "padding-bottom: 20px;",
               column(2, numericInput('size', 'Point size', value = .3, min = .1, max = 1, step = .1)),
               column(3, sliderInput('nbins', 'Number of bins', value = 100, min = 100, max = 1000, step = 50)),
               column(3, sliderInput('BaseMean', 'BaseMean region', value = c(90, 130), min = 80, max = max(dataset[, round(BaseMean)]))), 
               column(2, numericInput('BaseMeanTU1', 'BaseMean (Min)', value = dataset[, .N, round(BaseMean)][which.max(N), round] - 1, min = 80, max = max(dataset[, round(BaseMean)]), step = .25)),
               column(2, numericInput('BaseMeanTU2', 'BaseMean (Max)', value = dataset[, .N, round(BaseMean)][which.max(N), round] + 1, min = 80, max = max(dataset[, round(BaseMean)]), step = .25))
      ),
      fluidRow(column(6, plotlyOutput('BaseMeanSelect1')), column(6, plotlyOutput('BaseMeanSelect2'))),
      fluidRow(style = "padding-bottom: 20px;",
               column(4, sliderInput("DeltaMean", "Max DeltaMean", value = 1, min = 0, max = 2, step = 0.1)),
               column(4, sliderInput("CurrentSD", "CurrentSD", value = 15, min = 0, max = 40, step = 1)), 
               column(4, sliderInput("StageSD", "StageSD", value = 3.5, min = 0, max = 20, step = 0.1))
      ),
      fluidRow(style = "padding-bottom: 20px;",
               column(4, sliderInput("SignalCurrentWidth", "SignalCurrentWidth", value = c(5, 15), min = 0, max = ceiling(dataset[, max(SignalCurrentWidth)]), step = 0.1)),
               column(4, sliderInput("SignalCurrentPercent", "SignalCurrentPercent", min = 0, max = 100, value = c(80, 100), step = 5)),
               column(4, sliderInput("DwellTime", "DwellTime (ms)", value = c(0, 30), min = 0, max = 100, step = .1))
      ),
      plotlyOutput('PointSelect'),
      actionButton("goButton", 'Save', class = "btn-success"),
      verbatimTextOutput("sp")
    ), 
    server = function(input, output, session) {
      selectedData <- reactive({
        dataset[BaseMean > input$BaseMean[1] & BaseMean < input$BaseMean[2]]
      })
      
      output$BaseMeanSelect1 = renderPlotly({
        ggplot(selectedData(), aes(x = BaseMean, y = Blockade)) + 
          geom_point(size = input$size) + 
          # geom_hline(yintercept = AABlockade[amino_acid %in% selectedData()[, unique(A)], Blockade], color = "red") + 
          theme_minimal() -> p
        # plot_ly(data = selectedData(), x = ~BaseMean, y = ~Blockade) %>% plot_ly(y = AABlockade[, Blockade])
        ggplotly(p)
      })
      
      output$BaseMeanSelect2 = renderPlotly({
        ggplot(selectedData(), aes(x = BaseMean)) + 
          geom_histogram(bins = input$nbins) + 
          geom_vline(xintercept = c(input$BaseMeanTU1, input$BaseMeanTU2), color = "red") + 
          theme_minimal() -> p2
        ggplotly(p2)
      })
      
      selectedData2 <- reactive({
        dataset[DeltaMean < input$DeltaMean & 
                  CurrentSD < input$CurrentSD & 
                  StageSD < input$StageSD & 
                  Blockade > 0 & Blockade < 0.5 &
                  SignalCurrentWidth > min(input$SignalCurrentWidth) & SignalCurrentWidth <= max(input$SignalCurrentWidth) & 
                  SignalCurrentPercent >= min(input$SignalCurrentPercent) & SignalCurrentPercent <= max(input$SignalCurrentPercent) & 
                  DwellTime > min(input$DwellTime) & DwellTime < max(input$DwellTime) & 
                  BaseMean > input$BaseMeanTU1 & BaseMean < input$BaseMeanTU2]
      })
      
      output$PointSelect <- renderPlotly({
        ggplot(selectedData2(), aes(x = Blockade, y = DwellTime)) + 
          geom_point(size = input$size) + 
          # geom_vline(xintercept = AABlockade[amino_acid %in% selectedData()[, unique(A)], Blockade], color = "red") + 
          theme_minimal() + 
          labs(title = selectedData()[, unique(A)]) -> p3
        ggplotly(p3) %>% layout(dragmode = "select")
      })
      
      output$sp <- renderPrint({
        d <- (event_data("plotly_selected"))
        if(input$goButton != 0) {
          if(!is.null(d)) {
            subdata <- selectedData2()[d$pointNumber + 1]
            fwrite(subdata, paste0(outdir, "/", basename(SigFile)), sep = "\t", quote = F)
            cat("Finished")
          } else {
            cat("Number of selected points:", nrow(d))
          }
        } else {
          cat("Number of selected points:", nrow(d))
        }
      })
    }
  )
}



meta[, file := paste0("./analysis/21.ABFProcessing/02.Background/RawSignal/RawSignal_", file_name, ".txt")]
meta[, file.exists(file)]
meta <- meta[file.exists(file)]
meta <- meta[, .SD[sample(.N, 1)], amino_acid]

StandardAA(SigFile = meta[21, file], target = "noise")



noise <- do.call(rbind, lapply(list.files("analysis/22.SignalSelecting/02.Background/noise", full.names = T), fread))

ggplot(noise, aes(x = Blockade, y = DwellTime * 1000)) + 
  geom_point() + 
  scale_y_log10()






