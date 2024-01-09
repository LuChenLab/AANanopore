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
meta <- meta[amino_acid %in% Biostrings::AMINO_ACID_CODE]
meta <- meta[concentration != 0]
meta$file_path <- mapply(meta$file_name, FUN = function(x) list.files("./data", recursive = TRUE, full.names = TRUE, pattern = as.character(x))[1])
meta <- meta[!is.na(file_path)]
meta[, start_time := as.numeric(start_time)]
meta[, end_time := as.numeric(end_time)]
# meta <- meta[!grepl("abf", file_name)]
meta[, file_name := gsub(".abf$", "", file_name)]

meta_Val <- rbind(data.table(openxlsx::read.xlsx("./data/meta_info_and_base_line_20210525.xlsx")), data.table(openxlsx::read.xlsx("./data/meta_info_and_base_line_20210525additional.xlsx")))[amino_acid == "Val"]
meta_Val$file_path <- mapply(meta_Val$file_name, FUN = function(x) list.files("./data", recursive = TRUE, full.names = TRUE, pattern = as.character(x))[1])
meta_Val[, L1min := NULL]
meta_Val[, L1max := NULL]
meta <- rbind(meta, meta_Val)
meta <- meta[!grepl("abf", file_name)]


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
          geom_hline(yintercept = AABlockade[amino_acid %in% selectedData()[, unique(A)], Blockade], color = "red") + 
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
          geom_vline(xintercept = AABlockade[amino_acid %in% selectedData2()[, unique(A)], Blockade], color = "red") + 
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

StateSignalSelecting <- function(workdir = "analysis/22.SignalSelecting/03.DifferentStates", target = NULL) {
  if(is.null(workdir)) stop()
  if(is.null(target)) stop()
  stopifnot(is.element(target, list.dirs(workdir, full.names = F)))
  outdir <- gsub("//", "/", paste0(workdir, "/", target))
  dataset <- do.call(rbind, lapply(list.files(outdir, "MainL0.txt", full.names = T), fread))
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
          geom_vline(xintercept = AABlockade[amino_acid %in% target, Blockade], color = "red") + 
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
          selectedData()[, Selected := (ID %in% subdata[D >= input$density, ID]) | (between(Blockade, b1, b2) & log10(DwellTime) > d2 & ID %in% subdata[, ID])]
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

# /mnt/raid61/Personal_data/tangchao/AANanopore/analysis/22.SignalSelecting/01.StandardAA

files <- list.files("./analysis/21.ABFProcessing/01.StandardAA/RawSignal", full.names = TRUE)
files <- grep(pattern = "Selected", files, invert = T, value = T)
files <- grep(pattern = "Feature_Mat", files, invert = T, value = T)


aat <- "Ala"
files_info <- meta[amino_acid == aat]
files_info[, file := paste0("./analysis/21.ABFProcessing/01.StandardAA/RawSignal/RawSignal_", file_name, ".txt")]
files_info[, file.exists(file)]
BaseMeanFiltering(SigFile = files_info[1, file], target = aat)
BaseMeanFiltering(SigFile = files_info[4, file], target = aat)

StateSignalSelecting(target = aat)





meta[, unique(amino_acid)]
aat <- "Arg"
files_info <- meta[amino_acid == aat]
files_info[, file := paste0("./analysis/21.ABFProcessing/01.StandardAA/RawSignal/RawSignal_", file_name, ".txt")]
files_info[, file.exists(file)]
BaseMeanFiltering(SigFile = files_info[1, file], target = aat)
BaseMeanFiltering(SigFile = files_info[2, file], target = aat)
BaseMeanFiltering(SigFile = files_info[3, file], target = aat)
BaseMeanFiltering(SigFile = files_info[5, file], target = aat)
BaseMeanFiltering(SigFile = files_info[6, file], target = aat)

StateSignalSelecting(target = aat)


meta[, unique(amino_acid)]
aat <- "Asn"
files_info <- meta[amino_acid == aat]
files_info[, file := paste0("./analysis/21.ABFProcessing/01.StandardAA/RawSignal/RawSignal_", file_name, ".txt")]
files_info[, file.exists(file)]
BaseMeanFiltering(SigFile = files_info[1, file], target = aat)
# BaseMeanFiltering(SigFile = files_info[2, file], target = aat)
BaseMeanFiltering(SigFile = files_info[3, file], target = aat)
BaseMeanFiltering(SigFile = files_info[4, file], target = aat)
BaseMeanFiltering(SigFile = files_info[5, file], target = aat)
BaseMeanFiltering(SigFile = files_info[6, file], target = aat)

StateSignalSelecting(target = aat)


meta[, unique(amino_acid)]
aat <- "Asp"
files_info <- meta[amino_acid == aat]
files_info[, file := paste0("./analysis/21.ABFProcessing/01.StandardAA/RawSignal/RawSignal_", file_name, ".txt")]
files_info[, file.exists(file)]
BaseMeanFiltering(SigFile = files_info[1, file], target = aat)
BaseMeanFiltering(SigFile = files_info[2, file], target = aat)
BaseMeanFiltering(SigFile = files_info[3, file], target = aat)
# BaseMeanFiltering(SigFile = files_info[4, file], target = aat)
# BaseMeanFiltering(SigFile = files_info[5, file], target = aat)
# BaseMeanFiltering(SigFile = files_info[6, file], target = aat)

StateSignalSelecting(target = aat)




meta[, unique(amino_acid)]
aat <- "Gln" # 1 & 4 很奇怪
files_info <- meta[amino_acid == aat]
files_info[, file := paste0("./analysis/21.ABFProcessing/01.StandardAA/RawSignal/RawSignal_", file_name, ".txt")]
files_info[, file.exists(file)]
# BaseMeanFiltering(SigFile = files_info[1, file], target = aat)
BaseMeanFiltering(SigFile = files_info[2, file], target = aat)
BaseMeanFiltering(SigFile = files_info[3, file], target = aat)
# BaseMeanFiltering(SigFile = files_info[4, file], target = aat)
BaseMeanFiltering(SigFile = files_info[5, file], target = aat)
BaseMeanFiltering(SigFile = files_info[6, file], target = aat)

StateSignalSelecting(target = aat)



meta[, unique(amino_acid)]
aat <- "Glu"
files_info <- meta[amino_acid == aat]
files_info[, file := paste0("./analysis/21.ABFProcessing/01.StandardAA/RawSignal/RawSignal_", file_name, ".txt")]
files_info[, file.exists(file)]
BaseMeanFiltering(SigFile = files_info[1, file], target = aat)
BaseMeanFiltering(SigFile = files_info[2, file], target = aat)
BaseMeanFiltering(SigFile = files_info[3, file], target = aat)
BaseMeanFiltering(SigFile = files_info[4, file], target = aat)
BaseMeanFiltering(SigFile = files_info[5, file], target = aat)
BaseMeanFiltering(SigFile = files_info[6, file], target = aat)

StateSignalSelecting(target = aat)


meta[, unique(amino_acid)]
aat <- "Gly" # 2 & 5 很奇怪
files_info <- meta[amino_acid == aat]
files_info[, file := paste0("./analysis/21.ABFProcessing/01.StandardAA/RawSignal/RawSignal_", file_name, ".txt")]
files_info[, file.exists(file)]
BaseMeanFiltering(SigFile = files_info[1, file], target = aat)
# BaseMeanFiltering(SigFile = files_info[2, file], target = aat)
BaseMeanFiltering(SigFile = files_info[3, file], target = aat)
BaseMeanFiltering(SigFile = files_info[4, file], target = aat)
# BaseMeanFiltering(SigFile = files_info[5, file], target = aat)
BaseMeanFiltering(SigFile = files_info[6, file], target = aat)

StateSignalSelecting(target = aat)


meta[, unique(amino_acid)]
aat <- "His" # SignalCurrentPercent >= 50
files_info <- meta[amino_acid == aat]
files_info[, file := paste0("./analysis/21.ABFProcessing/01.StandardAA/RawSignal/RawSignal_", file_name, ".txt")]
files_info[, file.exists(file)]
BaseMeanFiltering(SigFile = files_info[1, file], target = aat)
BaseMeanFiltering(SigFile = files_info[2, file], target = aat)
BaseMeanFiltering(SigFile = files_info[3, file], target = aat)
BaseMeanFiltering(SigFile = files_info[4, file], target = aat)
BaseMeanFiltering(SigFile = files_info[5, file], target = aat)
BaseMeanFiltering(SigFile = files_info[6, file], target = aat)

StateSignalSelecting(target = aat)
file.rename("./analysis/22.SignalSelecting/03.DifferentStates/His/His_State1.txt", 
            "./analysis/22.SignalSelecting/03.DifferentStates/His/His2_State1.txt")
StateSignalSelecting(target = aat)


meta[, unique(amino_acid)]
aat <- "Ile" 
files_info <- meta[amino_acid == aat]
files_info[, file := paste0("./analysis/21.ABFProcessing/01.StandardAA/RawSignal/RawSignal_", file_name, ".txt")]
files_info[, file.exists(file)]
BaseMeanFiltering(SigFile = files_info[1, file], target = aat)
BaseMeanFiltering(SigFile = files_info[2, file], target = aat)
BaseMeanFiltering(SigFile = files_info[3, file], target = aat)
BaseMeanFiltering(SigFile = files_info[4, file], target = aat)
BaseMeanFiltering(SigFile = files_info[5, file], target = aat)
BaseMeanFiltering(SigFile = files_info[6, file], target = aat)

StateSignalSelecting(target = aat)


meta[, unique(amino_acid)]
aat <- "Leu" 
files_info <- meta[amino_acid == aat]
files_info[, file := paste0("./analysis/21.ABFProcessing/01.StandardAA/RawSignal/RawSignal_", file_name, ".txt")]
files_info[, file.exists(file)]
BaseMeanFiltering(SigFile = files_info[1, file], target = aat)
BaseMeanFiltering(SigFile = files_info[2, file], target = aat)
BaseMeanFiltering(SigFile = files_info[3, file], target = aat)
BaseMeanFiltering(SigFile = files_info[4, file], target = aat)
BaseMeanFiltering(SigFile = files_info[5, file], target = aat)
BaseMeanFiltering(SigFile = files_info[6, file], target = aat)

StateSignalSelecting(target = aat)


meta[, unique(amino_acid)]
aat <- "Lys" 
files_info <- meta[amino_acid == aat]
files_info[, file := paste0("./analysis/21.ABFProcessing/01.StandardAA/RawSignal/RawSignal_", file_name, ".txt")]
files_info[, file.exists(file)]
BaseMeanFiltering(SigFile = files_info[1, file], target = aat)
BaseMeanFiltering(SigFile = files_info[2, file], target = aat)
BaseMeanFiltering(SigFile = files_info[3, file], target = aat)
# BaseMeanFiltering(SigFile = files_info[4, file], target = aat)
BaseMeanFiltering(SigFile = files_info[5, file], target = aat)
BaseMeanFiltering(SigFile = files_info[6, file], target = aat)

StateSignalSelecting(target = aat)


meta[, unique(amino_acid)]
aat <- "Met" # 2 & 5 很奇怪
files_info <- meta[amino_acid == aat]
files_info[, file := paste0("./analysis/21.ABFProcessing/01.StandardAA/RawSignal/RawSignal_", file_name, ".txt")]
files_info[, file.exists(file)]
BaseMeanFiltering(SigFile = files_info[1, file], target = aat)
# BaseMeanFiltering(SigFile = files_info[2, file], target = aat)
BaseMeanFiltering(SigFile = files_info[3, file], target = aat)
BaseMeanFiltering(SigFile = files_info[4, file], target = aat)
# BaseMeanFiltering(SigFile = files_info[5, file], target = aat)
BaseMeanFiltering(SigFile = files_info[6, file], target = aat)

StateSignalSelecting(target = aat)


meta[, unique(amino_acid)]
aat <- "Phe" # 1 很奇怪
files_info <- meta[amino_acid == aat]
files_info[, file := paste0("./analysis/21.ABFProcessing/01.StandardAA/RawSignal/RawSignal_", file_name, ".txt")]
files_info[, file.exists(file)]
# BaseMeanFiltering(SigFile = files_info[1, file], target = aat)
BaseMeanFiltering(SigFile = files_info[2, file], target = aat)
BaseMeanFiltering(SigFile = files_info[3, file], target = aat)
BaseMeanFiltering(SigFile = files_info[4, file], target = aat)
BaseMeanFiltering(SigFile = files_info[5, file], target = aat)
BaseMeanFiltering(SigFile = files_info[6, file], target = aat)

StateSignalSelecting(target = aat)


meta[, unique(amino_acid)]
aat <- "Pro"
files_info <- meta[amino_acid == aat]
files_info[, file := paste0("./analysis/21.ABFProcessing/01.StandardAA/RawSignal/RawSignal_", file_name, ".txt")]
files_info[, file.exists(file)]
BaseMeanFiltering(SigFile = files_info[1, file], target = aat)
BaseMeanFiltering(SigFile = files_info[2, file], target = aat)
# BaseMeanFiltering(SigFile = files_info[3, file], target = aat)
# BaseMeanFiltering(SigFile = files_info[4, file], target = aat)
# BaseMeanFiltering(SigFile = files_info[5, file], target = aat)
BaseMeanFiltering(SigFile = files_info[6, file], target = aat)

StateSignalSelecting(target = aat)



meta[, unique(amino_acid)]
aat <- "Ser" # 1 很奇怪
files_info <- meta[amino_acid == aat]
files_info[, file := paste0("./analysis/21.ABFProcessing/01.StandardAA/RawSignal/RawSignal_", file_name, ".txt")]
files_info[, file.exists(file)]
# BaseMeanFiltering(SigFile = files_info[1, file], target = aat)
BaseMeanFiltering(SigFile = files_info[2, file], target = aat)
BaseMeanFiltering(SigFile = files_info[3, file], target = aat)
BaseMeanFiltering(SigFile = files_info[4, file], target = aat)
BaseMeanFiltering(SigFile = files_info[5, file], target = aat)
BaseMeanFiltering(SigFile = files_info[6, file], target = aat)

StateSignalSelecting(target = aat)


meta[, unique(amino_acid)]
aat <- "Thr" # 1 & 2很奇怪
files_info <- meta[amino_acid == aat]
files_info[, file := paste0("./analysis/21.ABFProcessing/01.StandardAA/RawSignal/RawSignal_", file_name, ".txt")]
files_info[, file.exists(file)]
# BaseMeanFiltering(SigFile = files_info[1, file], target = aat)
# BaseMeanFiltering(SigFile = files_info[2, file], target = aat)
BaseMeanFiltering(SigFile = files_info[3, file], target = aat)
BaseMeanFiltering(SigFile = files_info[4, file], target = aat)
BaseMeanFiltering(SigFile = files_info[5, file], target = aat)


StateSignalSelecting(target = aat)


meta[, unique(amino_acid)]
aat <- "Trp" 
files_info <- meta[amino_acid == aat]
files_info[, file := paste0("./analysis/21.ABFProcessing/01.StandardAA/RawSignal/RawSignal_", file_name, ".txt")]
files_info[, file.exists(file)]
BaseMeanFiltering(SigFile = files_info[1, file], target = aat)
BaseMeanFiltering(SigFile = files_info[2, file], target = aat)
BaseMeanFiltering(SigFile = files_info[3, file], target = aat)
BaseMeanFiltering(SigFile = files_info[4, file], target = aat)
BaseMeanFiltering(SigFile = files_info[5, file], target = aat)
BaseMeanFiltering(SigFile = files_info[6, file], target = aat)

StateSignalSelecting(target = aat)


meta[, unique(amino_acid)]
aat <- "Tyr" # SignalCurrentPercent >= 60
files_info <- meta[amino_acid == aat]
files_info[, file := paste0("./analysis/21.ABFProcessing/01.StandardAA/RawSignal/RawSignal_", file_name, ".txt")]
files_info[, file.exists(file)]
BaseMeanFiltering(SigFile = files_info[1, file], target = aat)
BaseMeanFiltering(SigFile = files_info[2, file], target = aat)
BaseMeanFiltering(SigFile = files_info[3, file], target = aat)
BaseMeanFiltering(SigFile = files_info[4, file], target = aat)
BaseMeanFiltering(SigFile = files_info[5, file], target = aat)
BaseMeanFiltering(SigFile = files_info[6, file], target = aat)

StateSignalSelecting(target = aat)


meta[, unique(amino_acid)]
aat <- "Val"
files_info <- meta[amino_acid == aat]
files_info[, file := paste0("./analysis/21.ABFProcessing/01.StandardAA/RawSignal/RawSignal_", file_name, ".txt")]
files_info[, file.exists(file)]
# BaseMeanFiltering(SigFile = files_info[1, file], target = aat)
BaseMeanFiltering(SigFile = files_info[2, file], target = aat)
BaseMeanFiltering(SigFile = files_info[3, file], target = aat)
# BaseMeanFiltering(SigFile = files_info[4, file], target = aat)
BaseMeanFiltering(SigFile = files_info[5, file], target = aat)
BaseMeanFiltering(SigFile = files_info[6, file], target = aat)

StateSignalSelecting(target = aat)



meta[, unique(amino_acid)]
aat <- "Cys"
files_info <- meta[amino_acid == aat]
files_info[, file := paste0("./analysis/21.ABFProcessing/01.StandardAA/RawSignal/RawSignal_", file_name, ".txt")]
files_info[, file.exists(file)]
files_info <- files_info[file.exists(file)]
BaseMeanFiltering(SigFile = files_info[1, file], target = aat)
BaseMeanFiltering(SigFile = files_info[2, file], target = aat)
BaseMeanFiltering(SigFile = files_info[3, file], target = aat)
BaseMeanFiltering(SigFile = files_info[4, file], target = aat)
# BaseMeanFiltering(SigFile = files_info[5, file], target = aat)
BaseMeanFiltering(SigFile = files_info[6, file], target = aat)

StateSignalSelecting(target = aat)






