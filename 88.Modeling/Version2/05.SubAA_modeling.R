setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(caret)
library(doParallel)
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10)

cl <- makePSOCKcluster(30)
registerDoParallel(cl)

ms <- c("EFL", "FHL", "FLSY", "EFKS", "AFGLV", "AEFLN", "LQSTW", "DELQSTW", "FILNQSTW", "DEFILNQSTW", "DE", "AFGILMNPQSTVWY", "ADEFGILMNPQSTVWY", "ELS", "ELQSTW", "ABEILPQSTW", "ABEHIKLPQSTW", 
        "QTW", "BLQSTW", "ABILPQSTW", "ABHIKLPQSTW", "BLS", "AIPQTW", "ABQT", "ABFILQT", "ABFILQSTW", "ABEFILQSTW", "BQT", "BGIQST", "ABGIMQSTW", "ABEGIMQSTW", "ABEGIMQRSTW", "ACDEFGHIKLMNPQRSTVWY", 
        "LSGV", "LFGV", "LFSGV", "LVFAG", "LVFAK", "LVFAKG", "FMQSN")
ms <- unique(ms)
ms <- strsplit(ms, "")
ms <- mapply(function(x) paste(sort(x), collapse = ""), ms)
ms <- unique(ms)

lapply(ms, function(i) {
  print(i)
  AAs <- strsplit(i, "")[[1]]
  AAs <- plyr::mapvalues(AAs, names(AMINO_ACID_CODE), AMINO_ACID_CODE, warn_missing = F)
  if(is.element("Asx", AAs)) AAs <- gsub("Asx", "CMC", AAs)
  
  outdir <- paste0("./analysis/81.ABFProcessing/SelectedSignals/01.StandardAA/06.Modeling/Version2/L1/SubAAs/Model", i)
  if(!dir.exists(outdir)) dir.create(outdir)
  if(file.exists(paste0(outdir, "/RFmodel.Rds"))) return(NULL)
  
  load("./analysis/81.ABFProcessing/SelectedSignals/01.StandardAA/06.Modeling/Version2/L1/SubAAs/Modeling_Data.RData")
  rm(list = c("Sigs_Test_FM", "Sigs_Test", "Sigs_Train_Upsample"))
  Sigs_Train_Upsample_FM <- subset(Sigs_Train_Upsample_FM, Class %in% AAs); gc()
  Sigs_Train_Upsample_FM <- droplevels(Sigs_Train_Upsample_FM)
  
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
  saveRDS(Fit1, paste0(outdir, "/RFmodel.Rds"))
})

stopCluster(cl)
