setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(parallel)
library(Biostrings)

get_density <- function(x, y, n = 1000, ...) {
  dens <- MASS::kde2d(x, y, n = n, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


files <- list.files("./analysis/32.SignalSelecting/01.StandardAA", full.names = TRUE, recursive = TRUE)
Selected_Sig <- lapply(files, fread)
Selected_Sig <- do.call(rbind, Selected_Sig)
Selected_Sig$File <- stringr::str_remove_all(Selected_Sig$ID, "_([[:digit:]]+)$")

ggplot(Selected_Sig, aes(x = Blockade, y = DwellTime, colour = File)) + 
  geom_point(size = .3) + 
  scale_y_log10() + 
  facet_wrap(~ A, scales = "free") + 
  theme(legend.position = "none")

ggplot(Selected_Sig, aes(x = Blockade, colour = File)) + 
  geom_line(stat = "density") + 
  facet_wrap(~ A, scales = "free") + 
  theme(legend.position = "none")

Selected_Sig[, AA:= plyr::mapvalues(A, AMINO_ACID_CODE, names(AMINO_ACID_CODE))]

Selected_Sig[, AA := factor(AA, levels = Selected_Sig[, mean(Blockade), AA][order(V1), AA])]
ggplot(Selected_Sig, aes(x = AA, y = Blockade)) + 
  geom_violin() + 
  stat_summary(fun.data = "mean_sd")

AAGroup <- list(GSA = c('G', 'S', 'A'), 
                TNRK = c("T", "N", "R", "K"), 
                QVML = c("Q", "V", "M", "L"), 
                IYDPFW = c("I", "Y", "D", "P", "F", "W"), 
                EH = c("E", "H"))
AAGroup <- data.table(Group = rep(names(AAGroup), mapply(length, AAGroup)), AA = unlist(AAGroup))

Selected_Sig <- merge(Selected_Sig, AAGroup, by = "AA")
ggplot(Selected_Sig, aes(x = Blockade, y = DwellTime, colour = AA)) + 
  geom_text(aes(label = AA)) + 
  scale_y_log10() + 
  facet_wrap(~ Group, scales = "free")






Selected_Sig_List <- split(Selected_Sig, Selected_Sig[, AA])
Selected_Sig_List <- lapply(Selected_Sig_List, function(x) {
  x$Den <- x[, get_density(Blockade, log10(DwellTime))]
  x$Den <- x$Den / max(x$Den)
  x
})
Selected_Sig_List <- do.call(rbind, Selected_Sig_List)

Selected_Sig_List_tu <- Selected_Sig_List[Den >= .99]
Selected_Sig_List_tu[, .N, AA]
ggplot(Selected_Sig_List_tu, aes(x = Blockade, y = DwellTime, colour = AA)) + 
  geom_text(aes(label = AA)) + 
  scale_y_log10() + 
  facet_wrap(~ Group, scales = "free")





Dis2Center <- lapply(Selected_Sig_List_tu[, unique(AA)], function(aat) {
  print(aat)
  Fl <- Selected_Sig[AA == aat, unique(File)]
  FM <- do.call(rbind, lapply(paste0("./analysis/21.ABFProcessing/01.StandardAA/FeatureMatrix/FeatureMatrix_", Fl, ".txt"), fread))
  FM <- FM[ID %in% Selected_Sig[AA == aat, ID]]
  FM999 <- FM[ID %in% Selected_Sig_List_tu[AA == aat, ID]]
  FM999 <- colMeans(FM999[, paste0("X", sprintf("%04d", 601:1000)), with = F])
  
  FM <- data.frame(FM[, paste0("X", sprintf("%04d", 601:1000)), with = F], row.names = FM[, ID])
  FM100 <- rbind(FM, t(data.frame(ID999 = FM999)))
  CorFM100 <- cor(t(FM100), method = "spearman")
  CorFM100 <- data.table(ID = row.names(CorFM100), Cor = CorFM100[, "ID999"])[ID != "ID999"]
  
  CosFM100 <- lsa::cosine(as.matrix(t(FM100)))
  CosFM100 <- data.table(ID = row.names(CosFM100), Cos = CosFM100[, "ID999"])[ID != "ID999"]
  
  # DisFM100 <- as.matrix(dist(as.matrix(FM100), method = "euclidean"))
  # DisFM100 <- data.table(ID = row.names(DisFM100), Dis = DisFM100[, "ID999"])[ID != "ID999"]
  # DisFM100$Dis <- DisFM100$Dis - min(DisFM100$Dis)
  # DisFM100$Dis <- DisFM100$Dis / max(DisFM100$Dis)
  
  # Simi <- merge(merge(CorFM100, CosFM100), DisFM100)
  merge(CorFM100, CosFM100)
})
Dis2Center <- do.call(rbind, Dis2Center)

Selected_Sig_Dis <- merge(Selected_Sig, Dis2Center, by = "ID")
Selected_Sig_Dis[, .N, AA]
Selected_Sig_Dis_tu <- Selected_Sig_Dis[Cos >= .995]

Selected_Sig_Dis_tu <- lapply(Selected_Sig_Dis[, unique(AA)], function(aat) {
  sub <- Selected_Sig_Dis[AA == aat][order(Cos)]
  if(sub[, sum(Cos >= .995)] > 100) {
    tail(sub[Cos >= .995], 100)
  } else {
    if(sub[, sum(Cos >= .995)] < 20) {
      tail(sub, 20)
    } else {
      sub[Cos >= .995]
    }
  }
})
Selected_Sig_Dis_tu <- do.call(rbind, Selected_Sig_Dis_tu)
Selected_Sig_Dis_tu <- unique(rbind(Selected_Sig_Dis_tu, Selected_Sig_Dis_tu[ID %in% Selected_Sig_List_tu[, ID]]))


Selected_Sig_Dis_tu[, .N, AA]
ggplot(Selected_Sig_Dis_tu, aes(x = Blockade, y = DwellTime, colour = AA)) + 
  geom_text(aes(label = AA)) + 
  scale_y_log10() + 
  facet_wrap(~ Group, scales = "free")

ggplot(Selected_Sig_Dis_tu, aes(x = Blockade, colour = AA)) + 
  geom_line(stat = "density", adjust = 1) + 
  facet_wrap(~ Group, scales = "free")











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
Selected_Sig_Upsample <- Selected_Sig_Dis_tu[, .SD[sample(.N, 2000, replace = T), ], AA]

ggplot(Selected_Sig_Upsample, aes(x = Blockade, colour = AA)) + 
  geom_line(stat = "density", adjust = 3)


files <- list.files("./analysis/22.SignalSelecting/02.Background/noise", full.names = TRUE)
Selected_Sig0 <- lapply(files, fread)
Selected_Sig0 <- do.call(rbind, Selected_Sig0)
Selected_Sig0[, A := "Noise"]
Selected_Sig0 <- Selected_Sig0[ID %in% Feature_Mat0$ID]

ggplot(Selected_Sig0, aes(x = Blockade, y = DwellTime, label = A, colour = A)) + 
  geom_text() + 
  scale_y_log10()

ggplot(Feature_Mat0[Blockade > 0 & Blockade < .3], aes(x = Blockade, y = DwellTime * 1000)) + 
  geom_point() + 
  scale_y_log10()



# Selected_Sig0[, Blockade2 := round(Blockade * 100)]
# set.seed(123)
# Selected_Sig0_2 <- rbind(Selected_Sig0[, .SD[sample(.N, 100, replace = T), ], Blockade2], 
#                          Selected_Sig0[Blockade < 0.3, .SD[sample(.N, 200, replace = T), ], Blockade2])
# Selected_Sig0_2[, Blockade2 := NULL]
# ggplot(Selected_Sig0_2, aes(x = Blockade)) + 
#   geom_histogram()
# ggplot(Selected_Sig0_2, aes(x = Blockade, y = DwellTime, label = A, colour = A)) + 
#   geom_text() + 
#   scale_y_log10()
# 

Selected_Sig_TU <- rbind(Selected_Sig_Upsample, Selected_Sig0, fill = TRUE)
Selected_Sig_TU[, A := plyr::mapvalues(A, Biostrings::AMINO_ACID_CODE, names(Biostrings::AMINO_ACID_CODE))]
Selected_Sig_TU[, A := factor(A, levels = c(setdiff(AA_STANDARD, "C"), "Noise"))]
setkey(Selected_Sig_TU, A)

Train <- setDF(Feature_Mat[Selected_Sig_TU[, ID], 2:1001])
Train <- Train[, c(paste0("X", sprintf("%04d", 601:1000)), "StageSD", "CurrentSD", "Segments", "SignalCurrentPercent", "SignalCurrentWidth", "DwellTime", "Blockade")]
Train$Class <- Selected_Sig_TU[, A]




saveRDS(Train, "./analysis/33.ModelTraining/Model1/Train.Rds")

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

Train <- readRDS("./analysis/33.ModelTraining/Model1/Train.Rds")
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
saveRDS(RFC, file = "./analysis/33.ModelTraining/Model1/RF.Rds")









