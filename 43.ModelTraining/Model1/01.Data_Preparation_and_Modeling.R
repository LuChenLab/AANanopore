setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(shiny)
library(plotly)
library(Biostrings)

get_density <- function(x, y, n = 1000, ...) {
  dens <- MASS::kde2d(x, y, n = n, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii] / max(dens$z[ii]))
}


AABlockade <- data.table(AA = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"), 
                         Blockade = c(0.14718, 0.17148, 0.16516, 0.21409, 0.16, 0.24528, 0.1882, 0.12072, 0.24652, 0.20722, 0.1995, 0.16875, 0.19772, 0.22018, 0.2183, 0.13131, 0.16101, 0.22744, 0.21276, 0.19044))
AABlockade$amino_acid <- plyr::mapvalues(AABlockade$AA, names(Biostrings::AMINO_ACID_CODE), Biostrings::AMINO_ACID_CODE)

AAGroup <- list(GSA = c('G', 'S', 'A'), 
                TNRK = c("T", "N", "R", "K"), 
                QVML = c("Q", "V", "M", "L"), 
                IYDPFW = c("I", "Y", "D", "P", "F", "W"), 
                EH = c("E", "H"))
AAGroup <- data.table(Group = rep(names(AAGroup), mapply(length, AAGroup)), AA = unlist(AAGroup))
AAGroup[, amino_acid := plyr::mapvalues(AA, names(AMINO_ACID_CODE), AMINO_ACID_CODE)]


Sigs <- lapply(AAGroup[, amino_acid], function(aat) {
  sigs <- do.call(rbind, lapply(list.files(paste0("./analysis/42.SignalSelecting/01.StandardAA/", aat), full.names = T), fread))
  sigs$File <- stringr::str_remove_all(sigs$ID, "_([[:digit:]]+)$")
  sigs$D <- sigs[, get_density(x = Blockade, y = log10(DwellTime))]
  sigs
})
Sigs <- do.call(rbind, Sigs)
Sigs <- merge(Sigs, AAGroup, by.x = "A", by.y = "amino_acid")

ggplot(Sigs[Group == "GSA" & D > .9], aes(x = Blockade, y = DwellTime, colour = AA)) + 
  geom_point() + 
  scale_y_log10()
ggplot(Sigs[Group == "GSA" & D > .9], aes(x = Blockade, fill = AA)) + 
  geom_histogram(binwidth = 0.0003)


ggplot(Sigs[Group == "TNRK" & D > .95], aes(x = Blockade, y = DwellTime, colour = AA)) + 
  geom_point() + 
  scale_y_log10()
ggplot(Sigs[Group == "TNRK" & D > .95], aes(x = Blockade, fill = AA)) + 
  geom_histogram(binwidth = 0.0001)


ggplot(Sigs[Group == "QVML" & D > .95], aes(x = Blockade, y = DwellTime, colour = AA)) + 
  geom_point() + 
  scale_y_log10()
ggplot(Sigs[Group == "QVML" & D > .95], aes(x = Blockade, fill = AA)) + 
  geom_histogram(binwidth = 0.0001)



ggplot(Sigs[Group == "IYDPFW" & D > .95], aes(x = Blockade, y = DwellTime, colour = AA)) + 
  geom_point() + 
  scale_y_log10()
ggplot(Sigs[Group == "IYDPFW" & D > .95], aes(x = Blockade, fill = AA)) + 
  geom_histogram(binwidth = 0.0001)


ggplot(Sigs[Group == "EH" & D > .9], aes(x = Blockade, y = DwellTime, colour = AA)) + 
  geom_point() + 
  scale_y_log10()
ggplot(Sigs[Group == "EH" & D > .9], aes(x = Blockade, fill = AA)) + 
  geom_histogram(binwidth = 0.0001)



Sigs999 <- rbind(Sigs[Group == "GSA" & D > .9], 
                 Sigs[(Group == "TNRK" & D > .95 & AA != "K") | (Group == "TNRK" & D > .95 & AA == "K" & Blockade > 0.17) ], 
                 Sigs[Group == "QVML" & D > .95], 
                 Sigs[Group == "IYDPFW" & D > .95], 
                 Sigs[Group == "EH" & D > .9])

ggplot(Sigs999, aes(x = Blockade, y = DwellTime)) + 
  geom_text(aes(label = AA, colour = AA)) + 
  scale_y_log10()












# Modeling


setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(parallel)
library(Biostrings)

files <- list.files("./analysis/21.ABFProcessing/02.Background/FeatureMatrix", pattern = "", full.names = TRUE, recursive = TRUE)
Feature_Mat0 <- lapply(files, fread)
Feature_Mat0 <- do.call(rbind, Feature_Mat0)
Feature_Mat0 <- na.omit(Feature_Mat0)

files <- list.files("./analysis/21.ABFProcessing/01.StandardAA/FeatureMatrix", pattern = "", full.names = TRUE, recursive = TRUE)
Feature_Mat <- lapply(files, fread)
Feature_Mat <- do.call(rbind, Feature_Mat)

setkey(Feature_Mat0, ID)
Feature_Mat <- rbind(Feature_Mat, Feature_Mat0)
setkey(Feature_Mat, ID)







set.seed(123)
Selected_Sig_Upsample <- Sigs999[, .SD[sample(.N, 2000, replace = T), ], AA]


files <- list.files("./analysis/22.SignalSelecting/02.Background/noise", full.names = TRUE)
Selected_Sig0 <- lapply(files, fread)
Selected_Sig0 <- do.call(rbind, Selected_Sig0)
Selected_Sig0[, A := "Noise"]
Selected_Sig0 <- Selected_Sig0[ID %in% Feature_Mat0$ID]




Selected_Sig_TU <- rbind(Selected_Sig_Upsample, Selected_Sig0, fill = TRUE)
Selected_Sig_TU[, A := plyr::mapvalues(A, Biostrings::AMINO_ACID_CODE, names(Biostrings::AMINO_ACID_CODE))]
Selected_Sig_TU[, A := factor(A, levels = c(setdiff(AA_STANDARD, "C"), "Noise"))]
setkey(Selected_Sig_TU, A)

Train <- setDF(Feature_Mat[Selected_Sig_TU[, ID], ])
Train <- Train[, c(paste0("X", sprintf("%04d", 601:1000)), "StageSD", "CurrentSD", "Segments", "SignalCurrentPercent", "SignalCurrentWidth", "Blockade")]
Train$Class <- Selected_Sig_TU[, A]




saveRDS(Train, "./analysis/43.ModelTraining/Model1/Train.Rds")

setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(ggplot2)
library(parallel)
library(caret)

fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10)

library(doParallel)
cl <- makePSOCKcluster(10)
registerDoParallel(cl)

Train <- readRDS("./analysis/43.ModelTraining/Model1/Train.Rds")
set.seed(825)
RFC <- train(Class ~ ., data = Train, 
             # preProc = c("center", "scale", "YeoJohnson", "nzv"), 
             method = "rf", 
             trControl = fitControl,
             verbose = FALSE,
             ## to evaluate:
             tuneGrid = expand.grid(mtry = 30),
             # tuneLength = 50,
             metric = "Accuracy", 
             allowParallel = TRUE)
saveRDS(RFC, file = "./analysis/43.ModelTraining/Model1/RF.Rds")









Train <- setDF(Feature_Mat[Sigs999[, ID], ])
Train <- Train[, c(paste0("X", sprintf("%04d", 601:1000)), paste0("P", sprintf("%03d", 1:100)), "StageSD", "CurrentSD", "Segments", "SignalCurrentPercent", "SignalCurrentWidth", "Blockade")]
row.names(Train) <- Sigs999[, ID]

pca.res <- FactoMineR::PCA(Train)
pca.coord <- merge(Sigs999, as.data.table(pca.res$ind$coord, keep.rownames = "ID"), by = "ID")

ggplot(pca.coord, aes(x = Dim.3, y = Dim.4)) + 
  geom_text(aes(colour = AA, label = AA))

pca.coord[, plotly::plot_ly(x = Blockade, y = log10(DwellTime), z = Dim.1, mode = "markers", text = AA)]

