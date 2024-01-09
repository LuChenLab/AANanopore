setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(ggplot2)
library(parallel)
library(caret)
library(dplyr)

M1 <- readRDS(file = "./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model0/RF/RF_Fit_L1.Rds")
M2 <- readRDS(file = "./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model0/RF/RF_Fit_L1_P.Rds")
M3 <- readRDS(file = "./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model0/RF/RF_Fit_L1_X.Rds")

files <- list.files("./analysis/11.SignalIdentification/FeatureMatrix", full.names = T)

lapply(files, FUN = function(f) {
  print(which(files == f))
  Tab <- na.omit(fread(f))
  Tab <- data.table(ID = Tab[, ID], 
                    M1Pred = plyr::mapvalues(predict(M1, Tab), Biostrings::AMINO_ACID_CODE, names(Biostrings::AMINO_ACID_CODE)), 
                    M1Prob = apply(predict(M1, Tab, type = "prob"), 1, max), 
                    M2Pred = plyr::mapvalues(predict(M2, Tab), Biostrings::AMINO_ACID_CODE, names(Biostrings::AMINO_ACID_CODE)), 
                    M2Prob = apply(predict(M2, Tab, type = "prob"), 1, max), 
                    M3Pred = plyr::mapvalues(predict(M3, Tab), Biostrings::AMINO_ACID_CODE, names(Biostrings::AMINO_ACID_CODE)), 
                    M3Prob = apply(predict(M3, Tab, type = "prob"), 1, max))
  Sig <- fread(gsub("/FeatureMatrix/FeatureMatrix", "/Signal/RawSignal", f))
  Sig <- na.omit(Sig)
  if(identical(Tab[, ID], Sig[, ID])) {
    res <- cbind(Sig, Tab[, 2:7])
    fwrite(res, paste0("./analysis/14.Prediction/01.ShinySelected/Model0", gsub("FeatureMatrix", "SignalPredicted", basename(f))), sep = "\t", quote = F)
  }
})


M1 <- readRDS(file = "./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model1/RF/RF_Fit_L1.Rds")
M2 <- readRDS(file = "./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model1/RF/RF_Fit_L1_P.Rds")
M3 <- readRDS(file = "./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model1/RF/RF_Fit_L1_X.Rds")

lapply(files, FUN = function(f) {
  print(which(files == f))
  Tab <- na.omit(fread(f))
  Tab <- data.table(ID = Tab[, ID], 
                    M1Pred = predict(M1, Tab), 
                    M1Prob = apply(predict(M1, Tab, type = "prob"), 1, max), 
                    M2Pred = predict(M2, Tab), 
                    M2Prob = apply(predict(M2, Tab, type = "prob"), 1, max), 
                    M3Pred = predict(M3, Tab), 
                    M3Prob = apply(predict(M3, Tab, type = "prob"), 1, max))
  Sig <- fread(gsub("/FeatureMatrix/FeatureMatrix", "/Signal/RawSignal", f))
  Sig <- na.omit(Sig)
  if(identical(Tab[, ID], Sig[, ID])) {
    res <- cbind(Sig, Tab[, 2:7])
    fwrite(res, paste0("./analysis/14.Prediction/01.ShinySelected/Model1/", gsub("FeatureMatrix", "SignalPredicted", basename(f))), sep = "\t", quote = F)
  }
})




