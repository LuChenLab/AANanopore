setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)

outdir <- paste0("./analysis/81.ABFProcessing/SelectedSignals/01.StandardAA/06.Modeling/Version2/L1/SubAAs/")
if(!file.exists(outdir)) dir.create(outdir)

files <- list.files("./analysis/81.ABFProcessing/SelectedSignals/01.StandardAA/03.GroupedSignals", ".signal.txt", full.names = T)
Sigs <- do.call(rbind, lapply(files, fread))[State == "Sington", ]
Sigs[, ID0 := ID]
for (a in Sigs[, unique(AA)]) {
  Sigs[, ID0 := gsub(paste0(a, "_"), "", ID0)]
}
Sigs[, AA := gsub("CbC", "CMC", AA)]
Sigs[, ID := gsub("CbC", "CMC", ID)]

Sigs <- Sigs[AA %in% c(AMINO_ACID_CODE[c("A", "F", "G", "K", "L", "S", "V")])]


meta1 <- as.data.table(openxlsx::read.xlsx("./analysis/91.RealTimeHydrolysis/meta_ZhangMing.xlsx"))
meta2 <- as.data.table(openxlsx::read.xlsx("./analysis/91.RealTimeHydrolysis/meta_ZhangMing2.xlsx"))

meta <- rbind(meta1[is.na(amino_acid)], meta2[is.na(amino_acid)])[experiment != "EFG Rep1"]
meta[, sg_files := paste0("./analysis/91.RealTimeHydrolysis/01.SelectedL0/", file_id, ".MainL0.txt")]
meta[, file.exists(sg_files)]
Sig0 <- lapply(meta$sg_files, fread)
Sig0 <- do.call(rbind, Sig0)[Blockade > 0 & Blockade < 0.3]

ggplot(Sigs, aes(x = Blockade, y = DwellTime)) + 
  geom_point() + 
  scale_y_log10()
ggplot(Sig0, aes(x = Blockade, y = DwellTime)) + 
  geom_point() + 
  scale_y_log10()


FMfs <- list.files("./analysis/81.ABFProcessing/FeatureMatrix", "FeatureMatrix_", full.names = T, recursive = T)
FMfs <- unlist(lapply(Sigs[, unique(File)], FUN = function(x) grep(x, FMfs, value = T)))
FMs <- lapply(FMfs, fread)
FMs <- do.call(rbind, FMs)
if(!all(Sigs$ID0 %in% FMs$ID)) return(NULL)
FMs <- FMs[ID %in% Sigs$ID0]
setkey(FMs, ID)


FM0 <- lapply(unique(meta[, paste0("./analysis/81.ABFProcessing/FeatureMatrix/FeatureMatrix_", file_name, ".txt")]), fread)
FM0 <- do.call(rbind, FM0)
FM0 <- FM0[ID %in% Sig0$ID]
setkey(FM0, ID)

set.seed(123)
Sig_tu <- rbind(Sigs[AA == "Ala", .(AA, ID = ID0)], 
                Sigs[AA == "Ala", .SD[sample(.N, 2000 - 857, replace = T), .(ID = ID0)], AA], 
                Sigs[AA == "Lys", .(AA, ID = ID0)], 
                Sigs[AA == "Lys", .SD[sample(.N, 2000 - 650, replace = T), .(ID = ID0)], AA], 
                Sigs[AA == "Gly", .(AA, ID = ID0)], 
                Sigs[AA == "Gly", .SD[sample(.N, 2000 - 1149, replace = T), .(ID = ID0)], AA], 
                Sigs[!AA %in% c("Ala", "Lys", "Gly"), .SD[sample(.N, 2000, replace = F), .(ID = ID0)], AA], 
                Sig0[, .SD[sample(.N, 2000, replace = F), .(AA = "Noise", ID)]])
setnames(Sig_tu, "AA", "Class")
FM_tu <- rbind(FMs, FM0)[ID %in% Sig_tu$ID]
FM_tu <- merge(FM_tu, Sig_tu, by = "ID")

Sigs_Train_Upsample_FM <- as.data.frame(FM_tu[, grepl("X", colnames(FM_tu)) | colnames(FM_tu) %in% c("DeltaMean", "StageSD", "CurrentSD", "Segments", "SignalCurrentPercent", "DwellTime", "Blockade", "Class"), with = F])
Sigs_Train_Upsample_FM$Class <- factor(Sigs_Train_Upsample_FM$Class)
save(Sigs_Train_Upsample_FM, file = "./analysis/91.RealTimeHydrolysis/00.Models/01.AFGKLSV/Sigs_Train_Upsample_FM.RData")

library(caret)
library(doParallel)
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10)

cl <- makePSOCKcluster(30)
registerDoParallel(cl)


set.seed(825)
Fit1 <- train(Class ~ ., 
              data = Sigs_Train_Upsample_FM, 
              # preProc = c("center", "scale", "YeoJohnson", "nzv"), 
              method = "rf", 
              trControl = fitControl,
              verbose = FALSE,
              ## to evaluate:
              tuneGrid = expand.grid(mtry = 60),
              # tuneLength = 6,
              metric = "Accuracy", 
              allowParallel = TRUE)
saveRDS(Fit1, "./analysis/91.RealTimeHydrolysis/00.Models/01.AFGKLSV/Model_AFGKLSV_Noise.Rds")

