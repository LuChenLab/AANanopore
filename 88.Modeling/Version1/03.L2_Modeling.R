setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(caret)
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10)

library(doParallel)
cl <- makePSOCKcluster(30)
registerDoParallel(cl)

outdir <- paste0("./analysis/81.ABFProcessing/SelectedSignals/01.StandardAA/06.Modeling/Version1/L2")
load(paste0(outdir, "/Modeling_Data.RData"))
rm(list = c("Sigs_Test_FM", "Sigs_Test", "Sigs_Train_Upsample")); gc()

# RF ----

model1 <- train(Class ~ ., data = Sigs_Train_Upsample_FM, 
                # preProc = c("center", "scale", "YeoJohnson", "nzv"), 
                method = "rf", 
                trControl = fitControl,
                verbose = FALSE,
                ## to evaluate:
                # tuneGrid = expand.grid(mtry = 30),
                tuneLength = 6,
                metric = "Accuracy", 
                allowParallel = TRUE)
saveRDS(model1, file = paste0(outdir, "/RF_model.Rds"))
rm("model1"); gc()

# NB ----

model2 <- train(Class ~ ., 
                data = Sigs_Train_Upsample_FM, 
                method = "naive_bayes", 
                # tuneGrid = expand.grid(data.frame(size = 12, decay = 0.1)),
                tuneLength = 6, 
                metric = "Accuracy", 
                trControl = fitControl)
saveRDS(model2, file = paste0(outdir, "/NB_model.Rds"))
rm("model2"); gc()

# NNet ----

model3 <- train(Class ~ ., 
                data = Sigs_Train_Upsample_FM, 
                method = "pcaNNet", 
                # tuneGrid = expand.grid(data.frame(size = 12, decay = 0.1)),
                tuneLength = 6, 
                metric = "Accuracy", 
                trControl = fitControl)
saveRDS(model3, file = paste0(outdir, "/NNet_model.Rds"))
rm("model3"); gc()

# KNN ----

model4 <- train(Class ~ ., 
                data = Sigs_Train_Upsample_FM, 
                method = "knn", 
                # tuneGrid = expand.grid(data.frame(size = 12, decay = 0.1)),
                tuneLength = 6, 
                metric = "Accuracy", 
                trControl = fitControl)
saveRDS(model4, file = paste0(outdir, "/KNN_model.Rds"))
rm("model4"); gc()

# CART ----

model5 <- train(Class ~ ., 
                data = Sigs_Train_Upsample_FM, 
                method = "treebag", 
                # tuneGrid = expand.grid(data.frame(size = 12, decay = 0.1)),
                tuneLength = 6, 
                metric = "Accuracy", 
                trControl = fitControl)
saveRDS(model5, file = paste0(outdir, "/CART_model.Rds"))
rm("model5"); gc()

# AdaBoost ----

model6 <- train(Class ~ ., 
                data = Sigs_Train_Upsample_FM, 
                method = "AdaBoost.M1", 
                # tuneGrid = expand.grid(data.frame(size = 12, decay = 0.1)),
                tuneLength = 6, 
                metric = "Accuracy", 
                trControl = fitControl)
saveRDS(model6, file = paste0(outdir, "/AdaBoost_model.Rds"))
rm("model6"); gc()

stopCluster(cl)