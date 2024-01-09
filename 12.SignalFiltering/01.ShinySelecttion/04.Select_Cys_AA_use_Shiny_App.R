setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(shiny)
library(plotly)

AABlockade <- data.table(AA = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"), 
                         Blockade = c(0.14718, 0.17148, 0.16516, 0.21409, 0.18622, 0.24528, 0.1882, 0.12072, 0.24652, 0.20722, 0.1995, 0.16875, 0.19772, 0.22018, 0.2183, 0.13131, 0.16101, 0.22744, 0.21276, 0.19044))
AABlockade$amino_acid <- plyr::mapvalues(AABlockade$AA, names(Biostrings::AMINO_ACID_CODE), Biostrings::AMINO_ACID_CODE)

StandardAA <- function(SigFile, outdir = "./analysis/12.SignalFiltering/AASignal/Shiny3/Cys") {
  dataset <- fread(SigFile)
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
               column(4, sliderInput("L2Ratio", "L2Ratio", value = c(0, 1), min = 0, max = 1, step = 0.01)),
               column(4, sliderInput("SignalCurrentPercent", "SignalCurrentPercent", min = 0, max = 100, value = c(50, 100), step = 5)),
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
          geom_hline(yintercept = AABlockade[amino_acid %in% selectedData()[, unique(A)], Blockade], color = "red") + 
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
                  L2Ratio >= input$L2Ratio[1] & L2Ratio <= input$L2Ratio[2] & 
                  SignalCurrentPercent >= input$SignalCurrentPercent[1] & SignalCurrentPercent <= input$SignalCurrentPercent[2] & 
                  DwellTime > input$DwellTime[1] & DwellTime < input$DwellTime[2] & 
                  BaseMean > input$BaseMeanTU1 & BaseMean < input$BaseMeanTU2]
      })
      
      output$PointSelect <- renderPlotly({
        ggplot(selectedData2(), aes(x = Blockade, y = DwellTime)) + 
          geom_point(size = input$size) + 
          geom_vline(xintercept = AABlockade[amino_acid %in% selectedData()[, unique(A)], Blockade], color = "red") + 
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

files <- list.files("./analysis/11.SignalIdentification/Jan07/AASignal", full.names = TRUE)
files <- grep(pattern = "Selected", files, invert = T, value = T)
files <- grep(pattern = "Feature_Mat", files, invert = T, value = T)

files_info <- data.table(file = files, AA = mapply(files, FUN = function(x) fread(x)[, unique(A)]))
files_info <- files_info[AA == "Cys"]

StandardAA(SigFile = files_info[7, file])


files <- list.files("./analysis/12.SignalFiltering/AASignal/Shiny3/Cys", full.names = T)
f <- files[1]
RawSignal <- fread(f)
Current <- readRDS(paste0("./analysis/11.SignalIdentification/Dec27/ABF_", gsub("RawSignal_", "", gsub(".txt", "", basename(f))), ".Rds"))

i = 23
ggplot(Current[Time >= RawSignal[i, StartTime] - 0.005 & Time <= RawSignal[i, EndTime] + 0.005], aes(x = Time, y = pA)) + 
  geom_step() + 
  geom_hline(yintercept = RawSignal[i, BaseMean] * (1 - AABlockade[AA == "C", Blockade]), col = 2) -> p


for (j in seq_along(files)) {
  f <- files[j]
  RawSignal <- fread(f)
  Current <- readRDS(paste0("./analysis/11.SignalIdentification/Dec27/ABF_", gsub("RawSignal_", "", gsub(".txt", "", basename(f))), ".Rds"))
  for (i in seq_len(nrow(RawSignal))) {
    print(paste(j, i))
    ggplot(Current[Time >= RawSignal[i, StartTime] - 0.005 & Time <= RawSignal[i, EndTime] + 0.005], aes(x = Time, y = pA)) + 
      geom_step() + 
      geom_hline(yintercept = RawSignal[i, BaseMean] * (1 - AABlockade[AA == "C", Blockade]), col = 2) -> p
    ggsave(paste0("./analysis/12.SignalFiltering/AASignal/Shiny3/Cys/Plot/", RawSignal[i, ID], ".pdf"), p, width = 4, height = 3)
  }
}


files <- list.files("./analysis/12.SignalFiltering/AASignal/Shiny3/Cys", "RawSignal_", full.names = T)
RawSignals <- do.call(rbind, lapply(files, fread))
Cys_ID <- fread("./analysis/12.SignalFiltering/AASignal/Shiny3/Cys/Cys_ID.txt")

ggplot(RawSignals, aes(x = Blockade, y = DwellTime, colour = ID %in% Cys_ID[N > 1, ID])) + 
  geom_point() + 
  scale_y_log10() + 
  scale_x_continuous(limits = c(0.1, 0.3))

Cys_ID[N > 1, ]

ggplot(RawSignals[ID %in% Cys_ID[N >= 1, ID]], aes(x = Blockade, y = DwellTime)) + 
  geom_point() + 
  scale_y_log10() + 
  scale_x_continuous(limits = c(0.1, 0.3))

ggplot(RawSignals[ID %in% Cys_ID[N >= 1, ID]], aes(x = Blockade)) + 
  geom_histogram(binwidth = .005) + 
  scale_x_continuous(limits = c(0.1, 0.3))


RawSignals2 <- RawSignals[ID %in% Cys_ID[N >= 1, ID]]

fwrite(RawSignals2, "./analysis/12.SignalFiltering/AASignal/Shiny3/Cys/RawSignal.txt", sep = "\t", quote = FALSE)

