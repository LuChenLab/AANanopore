setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(caret)
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10)

library(doParallel)
cl <- makePSOCKcluster(20)
registerDoParallel(cl)


outdir <- paste0("./analysis/81.ABFProcessing/SelectedSignals/01.StandardAA/06.Modeling/Version3/L1")
load(paste0(outdir, "/Modeling_Data.RData"))
rm(list = c("Sigs_Test_FM", "Sigs_Test", "Sigs_Train_Upsample", "Sigs_Valid", "Sigs_Valid_FM")); gc()


model1 <- train(Class ~ ., data = Sigs_Train_Upsample_FM, 
                # preProc = c("center", "scale", "YeoJohnson", "nzv"), 
                method = "rf", 
                trControl = fitControl,
                verbose = FALSE,
                ## to evaluate:
                tuneGrid = expand.grid(mtry = 60),
                # tuneLength = 6,
                metric = "Accuracy", 
                allowParallel = TRUE)
saveRDS(model1, file = paste0(outdir, "/RF_model.Rds"))
rm("model1"); gc()
