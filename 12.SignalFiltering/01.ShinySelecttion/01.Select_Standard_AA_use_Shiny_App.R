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
               column(4, sliderInput("StageSD", "StageSD", value = 5, min = 0, max = 20, step = 0.1))
      ),
      fluidRow(style = "padding-bottom: 20px;",
               column(4, sliderInput("L2Ratio", "L2Ratio", value = c(0, 1), min = 0, max = 1, step = 0.01)),
               column(4, sliderInput("SignalCurrentPercent", "SignalCurrentPercent", min = 0, max = 100, value = c(50, 100), step = 5)),
               column(4, sliderInput("DwellTime", "DwellTime (ms)", value = c(0, 20), min = 0, max = 50, step = .1))
      ),
      plotlyOutput('PointSelect'),
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
        cat("Number of selected points:", nrow(d))
        if(!is.null(d)) {
          subdata <- selectedData2()[d$pointNumber + 1]
          fwrite(subdata, gsub(".txt$", "_Selected.txt", SigFile), sep = "\t", quote = F)
        }
      })
    }
  )
}

files <- list.files("./analysis/11.SignalIdentification/Jan07/AASignal", full.names = TRUE)

i = 119
StandardAA(files[100])
ggplot(fread(files[i]), aes(x = Blockade, y = log10(DwellTime))) + 
  geom_point() + 
  geom_vline(xintercept = AABlockade[amino_acid %in% fread(files[i])[, unique(A)], Blockade], color = "red")
fread(files[i])


files2 <- list.files("./analysis/11.SignalIdentification/Jan07/AASignal", pattern = "_Selected.txt", full.names = TRUE)
files <- setdiff(files, gsub("_Selected", "", files2))

i = 17
StandardAA(files[i])

files2 <- list.files("./analysis/11.SignalIdentification/Jan07/AASignal", pattern = "_Selected.txt", full.names = TRUE)
files2 <- do.call(rbind, lapply(files2, fread))
files2[, .N, A]

tab <- files2[, mean_sd(Blockade), A]
tab <- merge(AABlockade, tab, by.y = "A", by.x = "amino_acid")
tab[, cor.test(Blockade, y)]
files2[, A := factor(A, levels = tab[order(y), amino_acid])]

ggplot(files2, aes(x = A, y = Blockade)) + 
  geom_violin() + 
  stat_summary(fun.data = "median_mad") + 
  geom_point(data = AABlockade, aes(x = amino_acid, y = Blockade), col = "red")

ggplot(files2, aes(x = Blockade, y = DwellTime)) + 
  geom_point(size = .1) + 
  scale_y_log10() + 
  geom_vline(data = AABlockade, aes(xintercept = Blockade), col = "red") + 
  facet_wrap(facets = ~ A, scales = "free_y", ncol = 1, strip.position = "right")


# 时间的上限可以调一下
# Gly

ggplot(files2, aes(x = Blockade)) + 
  geom_histogram(binwidth = .0025) + 
  facet_wrap(facets = ~ A, scales = "free_y")

ggplot(files2, aes(x = Blockade)) + 
  geom_line(stat = "density", adjust = 4) + 
  facet_wrap(facets = ~ A, scales = "free_y") + 
  geom_hline(yintercept = 0, lwd = 2, color = "white") + 
  theme_pubclean() + 
  scale_x_continuous(limits = c(0.08, .28))


