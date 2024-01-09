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

AABlockade <- data.table(AA = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"), 
                         Blockade = c(0.14718, 0.17148, 0.16516, 0.21409, 0.18622, 0.24528, 0.1882, 0.12072, 0.24652, 0.20722, 0.1995, 0.16875, 0.19772, 0.22018, 0.2183, 0.13131, 0.16101, 0.22744, 0.21276, 0.19044))
AABlockade$amino_acid <- plyr::mapvalues(AABlockade$AA, names(Biostrings::AMINO_ACID_CODE), Biostrings::AMINO_ACID_CODE)

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

meta[, file := paste0("./analysis/21.ABFProcessing/02.Background/RawSignal/RawSignal_", file_name, ".txt")]
meta[, file.exists(file)]
meta <- meta[file.exists(file)]

BaseMeanFiltering(SigFile = meta[128, file], target = "Background")


