setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(patchwork)
library(plotly)
library(ggpubr)
library(shiny)

BaseMeanFiltering <- function(SigFile, starttime = NULL, endtime = NULL, outdir = NULL, target = NULL, product = NULL, label = NULL) {
  if(is.null(outdir)) stop()
  if(is.null(target)) stop()
  dataset <- fread(SigFile)
  if(!is.null(starttime) & !is.null(endtime)) {
    dataset <- dataset[StartTime > starttime*60 & EndTime < endtime*60]
  }
  
  dataset[, A := target]
  dataset[, DwellTime := DwellTime * 1000]
  shinyApp(
    ui = fluidPage(
      fluidRow(style = "padding-bottom: 20px;",
               column(2, numericInput('size', 'Point size', value = .3, min = .1, max = 1, step = .1)),
               column(3, sliderInput('nbins', 'Number of bins', value = 500, min = 100, max = 1000, step = 50)),
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
          # geom_vline(xintercept = product, color = "red") +
          theme_minimal() + 
          labs(title = selectedData2()[, unique(A)]) -> p3
        if(!is.null(product)) {
          if(length(product) > 0) {
            if(!is.null(label) & length(label) == length(product)) {
              p3 <- p3 + geom_vline(xintercept = product, color = "red") + 
                scale_x_continuous(breaks = product, labels = label)
            } else {
              p3 <- p3 + geom_vline(xintercept = product, color = "red")
            }
          }
        }
        
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
              if(!dir.exists(outdir)) dir.create(outdir)
              fwrite(selectedData2(), paste0(outdir, "/", target, ".MainL0.txt"), sep = "\t", quote = F)
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



AABlockade <- lapply(list.files("./analysis/61.SignalSelecting/01.StandardAA/02.SelectedSignals", "_V2.Rds", full.names = T), function(x) readRDS(x)$Summary)
AABlockade <- do.call(rbind, AABlockade)
AABlockade <- AABlockade[State == "State1"]

meta <- as.data.table(openxlsx::read.xlsx("data/ChenShanchuan/20231024/20231024_data.xlsx", sheet = 1))
meta <- meta[实验目的 == "阿尔兹海默症突变位点" & grepl("20231024", 文件名)]
colnames(meta) <- c("file_name", "start_time", "end_time", "amino_acid", "file_id", "purpose", "sample", "baseline", "note")
meta[, file_name := gsub(".abf", "", file_name)]
setkey(meta, file_name)
meta <- meta[, .SD[, .(amino_acid, start_time, end_time, file_id = paste0(file_name, ".", seq_len(.N)))], file_name]
meta[, sg_files := paste0("./analysis/81.ABFProcessing/RawSignal/RawSignal_", file_name, ".txt")]
meta[, file.exists(sg_files)]
meta <- meta[file.exists(sg_files), ]

outdir <- file.path("analysis/84.ADpolypeptide/21.SelectedL0")

meta <- meta[!file_id %in% gsub(".MainL0.txt", "", list.files("analysis/84.ADpolypeptide/21.SelectedL0"))]
# i = 1
BaseMeanFiltering(SigFile = meta[i, sg_files], target = meta[i, file_id], outdir = outdir, starttime = meta[i, start_time], endtime = 10, 
                  product = as.numeric(AABlockade[AA %in% AMINO_ACID_CODE[unlist(strsplit(meta[i, amino_acid], ""))], Blockade]), 
                  label = AABlockade[AA %in% AMINO_ACID_CODE[unlist(strsplit(meta[i, amino_acid], ""))], AA])
(i <- i + 1); meta[i]

