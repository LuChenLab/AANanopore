setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(IRanges)
library(ggplot2)
library(parallel)
library(changepoint)

meta0 <- data.table(openxlsx::read.xlsx("./data/单氨基酸实验信息表20210517.xlsx", sheet = 1, cols = 8:14))
meta1 <- data.table(openxlsx::read.xlsx("./data/单氨基酸实验信息表20210517.xlsx", sheet = 1, cols = 15:21))
meta2 <- data.table(openxlsx::read.xlsx("./data/单氨基酸实验信息表20210517.xlsx", sheet = 2, cols = 1:7))
meta3 <- data.table(openxlsx::read.xlsx("./data/单氨基酸实验信息表20210517.xlsx", sheet = 3, cols = 1:7))
meta4 <- data.table(openxlsx::read.xlsx("./data/单氨基酸实验信息表20210517.xlsx", sheet = 4, cols = 1:7))
meta <- rbind(meta0, meta1, meta2, meta3, meta4, use.names = FALSE)
colnames(meta) <- c("file_name", "date", "amino_acid", "concentration", "start_time", "end_time", "type")
meta <- meta[!is.na(file_name)]
meta[amino_acid == "cys", amino_acid := "Cys"]
meta <- meta[amino_acid %in% Biostrings::AMINO_ACID_CODE]
meta <- meta[concentration != 0]
meta$file_path <- mapply(meta$file_name, FUN = function(x) list.files("./data", recursive = TRUE, full.names = TRUE, pattern = as.character(x))[1])
meta <- meta[!is.na(file_path)]
meta[, start_time := as.numeric(start_time)]
meta[, end_time := as.numeric(end_time)]

meta_Val <- rbind(data.table(openxlsx::read.xlsx("./data/meta_info_and_base_line_20210525.xlsx")), data.table(openxlsx::read.xlsx("./data/meta_info_and_base_line_20210525additional.xlsx")))[amino_acid == "Val"]
meta_Val$file_path <- mapply(meta_Val$file_name, FUN = function(x) list.files("./data", recursive = TRUE, full.names = TRUE, pattern = as.character(x))[1])
meta_Val[, L1min := NULL]
meta_Val[, L1max := NULL]
meta <- rbind(meta, meta_Val)
meta <- meta[!grepl("abf", file_name)]

meta[!file_name %in% gsub(".txt", "", gsub("RawSignal_", "", list.files("./analysis/11.SignalIdentification/Jan07")))]
meta <- meta[file_name %in% gsub(".txt", "", gsub("RawSignal_", "", list.files("./analysis/11.SignalIdentification/Jan07")))]


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
            fwrite(subdata, gsub(".txt$", "_Selected_V2.txt", SigFile), sep = "\t", quote = F)
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
StandardAA(files[119])

file.exists(gsub(".txt$", "_Selected_V2.txt", files[146]))

files_info <- data.table(file = files, AA = mapply(files, FUN = function(x) fread(x)[, unique(A)]))
files_info[AA == "Cys"]
StandardAA(files_info[AA == "Cys"][7, file])



files2 <- list.files("./analysis/11.SignalIdentification/Jan07/AASignal", "_Selected_V2.txt", full.names = TRUE)
file.remove(grep("21201008", files2, value = T))


Selec2 <- lapply(files2, fread)
Selec2 <- do.call(rbind, Selec2)
Selec2[, AA := plyr::mapvalues(A, Biostrings::AMINO_ACID_CODE, names(Biostrings::AMINO_ACID_CODE))]

ggplot(Selec2, aes(x = Blockade, y = DwellTime, label = AA, colour = AA)) + 
  geom_text() + 
  scale_y_log10()

Selec2[, AA := factor(AA, levels = Selec2[, mean(Blockade), AA][order(V1), AA])]
ggplot(Selec2, aes(x = AA, y = Blockade)) + 
  stat_summary(fun.data = "mean_sd")

Selec2$File <- mapply(function(x) x[1], strsplit(Selec2[, ID], "_"))

ggplot(Selec2[AA == AA_STANDARD[20]], aes(x = Blockade, y = DwellTime, label = AA, colour = File)) + 
  geom_text() + 
  scale_y_log10()


