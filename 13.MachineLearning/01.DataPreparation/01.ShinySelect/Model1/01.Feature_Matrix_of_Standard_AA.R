setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(parallel)

files2 <- list.files("./analysis/11.SignalIdentification/Jan07/AASignal", pattern = "_Selected.txt", full.names = TRUE)
List <- lapply(files2, fread)
names(List) <- gsub("RawSignal_", "", gsub("_Selected.txt", "", basename(files2)))

Mat_List <- mclapply(seq_along(List), function(i) {
  print(i)
  Current <- readRDS(paste0("./analysis/11.SignalIdentification/Dec27/SignalCurrent_", names(List)[i], ".Rds"))
  Mat <- List[[i]]
  Feature_Mat <- lapply(Mat$ID, function(j) {
    D <- density(Current[ID == j, Current], from = 0, to = 1, n = 500, adjust = 0.5)$y
    pA <- sort(Current[ID == j, Current])
    pA <- pA[round(seq(1, length(pA), length.out = 100))]
    D <- c(D, pA)
    cbind(t(data.frame(round(D, 4), row.names = c(paste0("X", sprintf("%03d", 1:500)), paste0("P", sprintf("%03d", 1:100))))), 
          Mat[ID == j, .(DeltaMean, StageSD, CurrentSD, Segments, SignalCurrentPercent, L2Ratio, DwellTime, Blockade, A, ID)])
  })
  do.call(rbind, Feature_Mat)
}, mc.cores = 1)
Mat_List1 <- do.call(rbind, Mat_List)




files2 <- list.files("./analysis/11.SignalIdentification/Jan07/BackgroundSignal", pattern = "_Selected.txt", full.names = TRUE)
List <- lapply(files2, fread)
names(List) <- gsub("RawSignal_", "", gsub("_Selected.txt", "", basename(files2)))

Mat_List <- mclapply(seq_along(List), function(i) {
  print(i)
  Current <- readRDS(paste0("./analysis/11.SignalIdentification/Dec27/SignalCurrent_", names(List)[i], ".Rds"))
  Mat <- List[[i]]
  Feature_Mat <- lapply(Mat$ID, function(j) {
    D <- density(Current[ID == j, Current], from = 0, to = 1, n = 500, adjust = 0.5)$y
    pA <- sort(Current[ID == j, Current])
    pA <- pA[round(seq(1, length(pA), length.out = 100))]
    D <- c(D, pA)
    cbind(t(data.frame(round(D, 4), row.names = c(paste0("X", sprintf("%03d", 1:500)), paste0("P", sprintf("%03d", 1:100))))), 
          Mat[ID == j, .(DeltaMean, StageSD, CurrentSD, Segments, SignalCurrentPercent, L2Ratio, DwellTime, Blockade, A, ID)])
  })
  do.call(rbind, Feature_Mat)
}, mc.cores = 1)
Mat_List2 <- do.call(rbind, Mat_List)
Mat_List2[, A := "Noise"]

ggplot(Mat_List2, aes(x = Blockade, y = DwellTime)) + 
  geom_point() + 
  scale_y_log10()

ggplot(Mat_List1, aes(x = Blockade, y = DwellTime)) + 
  geom_point() + 
  scale_y_log10()

FeatueMat <- rbind(Mat_List1, Mat_List2)
FeatueMat[, A := plyr::mapvalues(A, Biostrings::AMINO_ACID_CODE, names(Biostrings::AMINO_ACID_CODE))]

saveRDS(FeatueMat, file = "./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model1/FeatueMatrix.Rds")
