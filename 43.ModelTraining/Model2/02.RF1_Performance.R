setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(ggplot2)
library(parallel)
library(caret)

Test <- readRDS("./analysis/43.ModelTraining/Model2/Test1.Rds")
RFC <- readRDS("./analysis/43.ModelTraining/Model2/RF1.Rds")
RFC1.1 <- readRDS("./analysis/43.ModelTraining/Model2/RF1_1.Rds")
RFC1.2 <- readRDS("./analysis/43.ModelTraining/Model2/RF1_2.Rds")
RFC1.3 <- readRDS("./analysis/43.ModelTraining/Model2/RF1_3.Rds")
RFC1.4 <- readRDS("./analysis/43.ModelTraining/Model2/RF1_4.Rds")
RFC1.5 <- readRDS("./analysis/43.ModelTraining/Model2/RF1_5.Rds")

Res <- data.table(Obse = Test$Class, Pred = predict(RFC, Test), Prob = apply(predict(RFC, Test, type = "prob"), 1, max))
hist(Res[, Prob], breaks = 100)

Res[, Obse := factor(Obse, levels = c('G', 'S', 'A', 'T', 'N', 'R', 'K', 'V', 'Q', 'L', 'M', 'I', 'Y', 'D', 'P', 'F', 'W', 'E', 'H', 'Noise'))]
Res[, Pred := factor(Pred, levels = c('G', 'S', 'A', 'T', 'N', 'R', 'K', 'V', 'Q', 'L', 'M', 'I', 'Y', 'D', 'P', 'F', 'W', 'E', 'H', 'Noise'))]
cfM <- as.matrix(Res[, confusionMatrix(Obse, Pred)])
cfM <- t(t(cfM)/colSums(cfM))
pheatmap::pheatmap(cfM, cluster_rows = F, cluster_cols = F, display_numbers = T)

PredFun <- function(mat, model, p = 0) {
  res <- data.table(ID = rownames(Test), Pred = predict(model, mat), Prob = apply(predict(model, mat, type = "prob"), 1, max))
  res0 <- res[!Pred %in% c('G', 'Noise', 'T', 'N', 'R', 'K', 'Q', 'V', 'L', 'M', 'I', 'Y', 'D', 'P', 'F', 'W')]
  Id1 <- which(res[, Pred %in% c('G', 'Noise') & Prob >= p])
  Id2 <- which(res[, Pred %in% c('T', 'N', 'R', 'K') & Prob >= p])
  Id3 <- which(res[, Pred %in% c('Q', 'V') & Prob >= p])
  Id4 <- which(res[, Pred %in% c('L', 'M') & Prob >= p])
  Id5 <- which(res[, Pred %in% c('I', 'Y', 'D', 'P', 'F', 'W') & Prob >= p])
  if(length(Id1) > 0) {
    res1 <- data.table(ID = rownames(mat)[Id1], Pred = predict(RFC1.1, mat[Id1, ]), Prob = apply(predict(RFC1.1, mat[Id1, ], type = "prob"), 1, max))
    res0 <- rbind(res0, res1)
  }
  
  if(length(Id2) > 0) {
    res2 <- data.table(ID = rownames(mat)[Id2], Pred = predict(RFC1.2, mat[Id2, ]), Prob = apply(predict(RFC1.2, mat[Id2, ], type = "prob"), 1, max))
    res0 <- rbind(res0, res2)
  }
  
  if(length(Id3) > 0) {
    res3 <- data.table(ID = rownames(mat)[Id3], Pred = predict(RFC1.3, mat[Id3, ]), Prob = apply(predict(RFC1.3, mat[Id3, ], type = "prob"), 1, max))
    res0 <- rbind(res0, res3)
  }
  
  if(length(Id4) > 0) {
    res4 <- data.table(ID = rownames(mat)[Id4], Pred = predict(RFC1.4, mat[Id4, ]), Prob = apply(predict(RFC1.4, mat[Id4, ], type = "prob"), 1, max))
    res0 <- rbind(res0, res4)
  }
  
  if(length(Id5) > 0) {
    res5 <- data.table(ID = rownames(mat)[Id5], Pred = predict(RFC1.5, mat[Id5, ]), Prob = apply(predict(RFC1.5, mat[Id5, ], type = "prob"), 1, max))
    res0 <- rbind(res0, res5)
  }
  setkey(res0, ID)
  return(res0[res$ID])
}


Res2 <- PredFun(mat = Test, model = RFC)
Res2$Obse <- Test$Class
Res2[, Obse := factor(Obse, levels = c('G', 'S', 'A', 'T', 'N', 'R', 'K', 'V', 'Q', 'L', 'M', 'I', 'Y', 'D', 'P', 'F', 'W', 'E', 'H', 'Noise'))]
Res2[, Pred := factor(Pred, levels = c('G', 'S', 'A', 'T', 'N', 'R', 'K', 'V', 'Q', 'L', 'M', 'I', 'Y', 'D', 'P', 'F', 'W', 'E', 'H', 'Noise'))]

Res[, mean(Pred == Obse)]
Res2[, mean(Pred == Obse)]

