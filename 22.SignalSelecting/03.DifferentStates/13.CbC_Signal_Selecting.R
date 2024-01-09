setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(shiny)
library(plotly)
library(ggpubr)
get_density <- function(x, y, n = 1000, ...) {
  dens <- MASS::kde2d(x, y, n = n, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii] / max(dens$z[ii]))
}

BaseMeanFiltering <- function(SigFile, outdir = "analysis/22.SignalSelecting/03.DifferentStates", target = NULL) {
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
               column(3, radioButtons(inputId = "LogScale", label = "log scale:", choices = c("TRUE", "FALSE"), inline = T))
      ),
      fluidRow(column(6, plotlyOutput('BaseMeanBlockade')), column(6, plotlyOutput('BlockadeDwellTime'))),
      plotlyOutput('BaseMeanHistogram'),
      actionButton("goButton", 'Save', class = "btn-success"),
      verbatimTextOutput("sp")
    ), 
    server = function(input, output, session) {
      selectedData <- reactive({
        dataset[BaseMean > input$BaseMean[1] & BaseMean < input$BaseMean[2]]
      })
      
      output$BaseMeanBlockade = renderPlotly({
        ggplot(selectedData(), aes(x = BaseMean, y = Blockade)) + 
          geom_point(size = input$size) + 
          # geom_hline(yintercept = AABlockade[amino_acid %in% selectedData()[, unique(A)], Blockade], color = "red") + 
          theme_minimal() -> p
        ggplotly(p)
      })
      
      output$BaseMeanHistogram = renderPlotly({
        ggplot(selectedData(), aes(x = BaseMean)) + 
          geom_histogram(bins = input$nbins) + 
          theme_minimal() -> p2
        ggplotly(p2) %>% layout(dragmode = "select")
      })
      
      output$BlockadeDwellTime = renderPlotly({
        d <- (event_data("plotly_selected"))
        if(is.null(d)) {
          selectedData2 <- reactive({
            dataset
          })
        } else {
          selectedData2 <- reactive({
            dataset[BaseMean > min(d$x) & BaseMean < max(d$x)]
          })
        }
        ggplot(selectedData2(), aes(x = Blockade, y = DwellTime)) + 
          geom_point(size = input$size) + 
          # geom_vline(xintercept = AABlockade[amino_acid %in% selectedData2()[, unique(A)], Blockade], color = "red") + 
          theme_minimal() + 
          labs(title = selectedData2()[, unique(A)]) -> p3
        if(input$LogScale == "TRUE") {
          p3 <- p3 + scale_y_log10()
        }
        ggplotly(p3)
      })
      
      output$sp <- renderPrint({
        d <- (event_data("plotly_selected"))
        if(is.null(d)) {
          cat("Number of selected points:", nrow(d))
        } else {
          selectedData2 <- reactive({
            dataset[BaseMean > min(d$x) & BaseMean < max(d$x)]
          })
          if(input$goButton != 0) {
            if(!is.null(d)) {
              fwrite(selectedData2(), paste0(outdir, "/", gsub(".txt", ".MainL0.txt", basename(SigFile))), sep = "\t", quote = F)
              cat("Saved", nrow(selectedData2()), "points.")
            } else {
              cat("Number of selected points:", nrow(selectedData2()))
            }
          } else {
            cat("Number of selected points:", nrow(selectedData2()))
          }
        }
      })
    }
  )
}

StateSignalSelecting <- function(workdir = "analysis/22.SignalSelecting/03.DifferentStates", dataset = NULL, target = NULL) {
  if(is.null(workdir)) stop()
  if(is.null(target)) stop()
  stopifnot(is.element(target, list.dirs(workdir, full.names = F)))
  outdir <- gsub("//", "/", paste0(workdir, "/", target))
  if(nrow(dataset) == 0) stop()
  dataset$File <- paste0("File", stringr::str_remove_all(dataset$ID, "_([[:digit:]]+)$"))
  
  shinyApp(
    ui = fluidPage(
      fluidRow(style = "padding-bottom: 20px;",
               column(3, selectInput("Files", "File:", unique(dataset$File), multiple = TRUE)),
               column(3, sliderInput('Blockade', 'Blockade region', value = c(0, 1), min = 0, max = 1, step = 0.05)),
               column(3, sliderInput('bandwidth', 'Smoothing', value = 10, min = 1, max = 100, step = 1)), 
               column(3, numericInput('size', 'Point size', value = .3, min = .1, max = 1, step = .1))
      ),
      fluidRow(column(6, plotlyOutput('FileSelect')), column(6, plotlyOutput('PointSelect'))),
      fluidRow(column(6, plotlyOutput('SelectedPoint')), column(6, sliderInput('density', 'confidence', value = .1, min = 0, max = 1, step = 0.01)), 
               column(6, radioButtons(inputId = "State", label = "Saving", choices = c("Selecting", "State1", "State2"), inline = T)),
               column(6, verbatimTextOutput("sp"))
      )
    ),
    server = function(input, output, session) {
      selectedData <- reactive({
        if (length(input$Files) != 0) {
          dataset <- dataset[File %in% input$Files]
        }
        dataset[Blockade > min(input$Blockade) & Blockade < max(input$Blockade)]
      })
      
      output$FileSelect <- renderPlotly({
        ggplot(selectedData(), aes(x = Blockade, colour = File)) + 
          geom_line(stat = "density", bw = input$bandwidth/1000) + 
          # geom_vline(xintercept = AABlockade[amino_acid %in% target, Blockade], color = "red") + 
          theme_minimal() -> p1
        ggplotly(p1)
      })
      
      output$PointSelect <- renderPlotly({
        ggplot(selectedData(), aes(x = Blockade, y = log10(DwellTime), colour = get_density(x = Blockade, y = log10(DwellTime)))) + 
          geom_point(size = input$size) + 
          guides(colour = "none") + 
          # scale_y_log10() + 
          theme_minimal() -> p2
        ggplotly(p2) %>% layout(dragmode = "select")
      })
      
      output$SelectedPoint <- renderPlotly({
        d <- (event_data("plotly_selected"))
        if(!is.null(d)) {
          subdata <- selectedData()[d$pointNumber + 1]
          subdata$D <- subdata[, get_density(x = Blockade, y = log10(DwellTime))]
          b1 <- as.numeric(ggpubr::mean_sd(subdata[, Blockade])[2])
          b2 <- as.numeric(ggpubr::mean_sd(subdata[, Blockade])[3])
          d2 <- as.numeric(ggpubr::mean_sd(subdata[, log10(DwellTime)])[1])
          selectedData()[, Selected := (ID %in% subdata[D >= input$density, ID]) | (between(Blockade, b1, b2) & log10(DwellTime) > d2)]
        } else {
          subdata <- data.table(selectedData(), D = selectedData()[, get_density(x = Blockade, y = log10(DwellTime))])
          selectedData()[, Selected := ID %in% subdata[D >= input$density, ID]]
        }
        
        ggplot(selectedData(), aes(x = Blockade, y = DwellTime, colour = Selected)) + 
          geom_point(size = input$size) + 
          scale_y_log10() + 
          theme_minimal() -> p3
        ggplotly(p3)
      })
      
      output$sp <- renderPrint({
        if(input$State %in% c("Selecting")) {
          cat("Number of selected points:", selectedData()[Selected == TRUE, .N])
        } else {
          fwrite(selectedData()[Selected == TRUE, ], paste0(outdir, "/", target, "_", input$State, ".txt"), sep = "\t", quote = F)
          cat("Selected", selectedData()[Selected == TRUE, .N], "points.")
        }
      })
    }
  )
}


BG1 <- "./analysis/21.ABFProcessing/05.CarboxymethylCys/RawSignal/RawSignal_20230404_0004.txt"

BaseMeanFiltering(SigFile = BG1, target = aat)
file.rename(from = "./analysis/22.SignalSelecting/03.DifferentStates/CbC/RawSignal_20230404_0004.MainL0.txt", 
            to = "./analysis/22.SignalSelecting/03.DifferentStates/CbC/Background_RawSignal_20230404_0004.MainL0.txt")

BG2 <- "./analysis/21.ABFProcessing/05.CarboxymethylCys/RawSignal/RawSignal_20230407_0000.txt"
BaseMeanFiltering(SigFile = BG2, target = aat)
file.rename(from = "./analysis/22.SignalSelecting/03.DifferentStates/CbC/RawSignal_20230407_0000.MainL0.txt", 
            to = "./analysis/22.SignalSelecting/03.DifferentStates/CbC/Background_RawSignal_20230407_0000.MainL0.txt")

CbC1 <- "./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_20230404_0005.txt"
BaseMeanFiltering(SigFile = CbC1, target = aat)
file.rename(from = "./analysis/22.SignalSelecting/03.DifferentStates/CbC/RawSignal_20230404_0005.MainL0.txt", 
            to = "./analysis/22.SignalSelecting/03.DifferentStates/CbC/CbC_RawSignal_20230404_0005.MainL0.txt")

CbC2 <- "./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_20230407_0000.txt"
BaseMeanFiltering(SigFile = CbC2, target = aat)
file.rename(from = "./analysis/22.SignalSelecting/03.DifferentStates/CbC/RawSignal_20230407_0000.MainL0.txt", 
            to = "./analysis/22.SignalSelecting/03.DifferentStates/CbC/CbC_RawSignal_20230407_0000.MainL0.txt")

CbC3 <- "./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_20230407_0001.txt"
BaseMeanFiltering(SigFile = CbC3, target = aat)
file.rename(from = "./analysis/22.SignalSelecting/03.DifferentStates/CbC/RawSignal_20230407_0001.MainL0.txt", 
            to = "./analysis/22.SignalSelecting/03.DifferentStates/CbC/CbC_RawSignal_20230407_0001.MainL0.txt")



CbC_dataset <- do.call(rbind, lapply(list.files(outdir, "CbC", full.names = T), fread))
StateSignalSelecting(dataset = CbC_dataset, target = aat)


