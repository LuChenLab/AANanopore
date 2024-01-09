setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(patchwork)
library(parallel)
library(ggplot2)
library(ggrepel)
library(IRanges)
library(dbscan)

KnnSelecting <- function(signals, dists, dist_method = "euclidean", k = 20, noise = 0) {
  stopifnot(noise < k)
  if(is(dists, "dist")) {
    knn1 <- dbscan::kNN(dists, k = k)
    knn2 <- apply(knn1$id, 2, function(x) x %in% which(attr(dists, "Labels") %in% grep("^Noise", rownames(knn1$id), value = TRUE)))
  } else {
    knn1 <- dbscan::kNN(dists[[dist_method]], k = k)
    knn2 <- apply(knn1$id, 2, function(x) x %in% which(attr(dists[[dist_method]], "Labels") %in% grep("^Noise", rownames(knn1$id), value = TRUE)))
  }
  
  tu <- row.names(knn1$id)[which(rowSums(knn2) <= noise)]
  # tu <- union(tu, knn2[which(rowSums(knn2) <= noise), 1])
  grep("^Noise", tu, invert = TRUE, value = TRUE)
}

RemoveOutlier <- function(tree, n = 10, Keep = 0.95) {
  h <- length(tree$height)
  g <- cutree(tree, h = tree$height[h])
  g <- data.table(tip = names(g), cluster = g)
  
  while (g[, .N, cluster][N >= n, sum(N)/length(tree$labels)] > Keep) {
    h <- h - 1
    g <- cutree(tree, h = tree$height[h])
    g <- data.table(tip = names(g), cluster = g)
  }
  cat(paste0("Keep: ", round(g[, .N, cluster][N >= n, sum(N)/length(tree$labels)] * 100, 2), "%"))
  g[cluster %in% g[, .N, cluster][N >= n, cluster], tip]
}

AAs <- c(AMINO_ACID_CODE[1:20], "CbC")


# KNN = 10 ----

mclapply(AAs, function(aat) {
  if(file.exists(paste0("analysis/82.ABFProcessing/SelectedSignals/01.StandardAA/02.SelectedSignals/", aat, ".signal.dist.after.knn10.Rds"))) return(NULL)
  bgsignal <- fread(paste0("./analysis/82.ABFProcessing/SelectedSignals/", aat, "_background.txt"))
  bgsignal[, file_id := as.character(file_id)]
  bgsignal[, File := gsub("\\.[0-9]", "", file_id)]
  bgsignal[, ID := paste0("Noise_", ID)]
  
  dataset <- fread(paste0("./analysis/82.ABFProcessing/SelectedSignals/", aat, "_signal.txt"))
  dataset[, file_id := as.character(file_id)]
  dataset[, File := gsub("\\.[0-9]", "", file_id)]
  dataset[, ID := paste0(aat, "_", ID)]
  
  FeatureMatrix <- paste0("./analysis/82.ABFProcessing/FeatureMatrix/FeatureMatrix_", unique(dataset$File), ".txt")
  FeatureMatrix <- do.call(rbind, lapply(FeatureMatrix, fread))
  FeatureMatrix[, ID := paste0(aat, "_", ID)]
  FeatureMatrix <- FeatureMatrix[ID %in% dataset[, ID]]
  FeatureMatrix[, AA := aat]
  setkey(FeatureMatrix, ID)
  FeatureMatrix <- FeatureMatrix[dataset$ID, ]
  FeatureMatrix <- data.frame(FeatureMatrix[, grepl("^X", colnames(FeatureMatrix)) | grepl("DwellTime", colnames(FeatureMatrix)), with = F], row.names = FeatureMatrix[[1]])
  
  LineMat <- melt.data.table(as.data.table(FeatureMatrix, keep.rownames = "ID"), id.vars = "ID")[variable != "DwellTime"]
  LineMat[, x := as.numeric(gsub("^X", "", variable))]
  LineMat[, Blockade := 1 - x/1000]
  LineMat <- LineMat[, .SD[, .(alpha = value/sum(value), Blockade)], ID]
  
  alldists <- readRDS(file.path("./analysis/82.ABFProcessing/SelectedSignals/01.StandardAA/01.SignalDistance", paste0(aat, ".signal.dist.Rds")))
  sigs <- KnnSelecting(signals = dataset, dists = alldists, k = 10, noise = 0)
  
  subtab <- subset.data.frame(FeatureMatrix, row.names(FeatureMatrix) %in% sigs)
  res.ds <- stats::dist(subtab)
  res.hc <- stats::hclust(res.ds, method = "single")
  saveRDS(res.hc, file = paste0("analysis/82.ABFProcessing/SelectedSignals/01.StandardAA/02.SelectedSignals/", aat, ".signal.dist.after.knn10.Rds"))
}, mc.cores = 21)


# KNN = 15 ----


mclapply(AAs, function(aat) {
  if(file.exists(paste0("analysis/82.ABFProcessing/SelectedSignals/01.StandardAA/02.SelectedSignals/", aat, ".signal.dist.after.knn15.Rds"))) return(NULL)
  bgsignal <- fread(paste0("./analysis/82.ABFProcessing/SelectedSignals/", aat, "_background.txt"))
  bgsignal[, file_id := as.character(file_id)]
  bgsignal[, File := gsub("\\.[0-9]", "", file_id)]
  bgsignal[, ID := paste0("Noise_", ID)]
  
  dataset <- fread(paste0("./analysis/82.ABFProcessing/SelectedSignals/", aat, "_signal.txt"))
  dataset[, file_id := as.character(file_id)]
  dataset[, File := gsub("\\.[0-9]", "", file_id)]
  dataset[, ID := paste0(aat, "_", ID)]
  
  FeatureMatrix <- paste0("./analysis/82.ABFProcessing/FeatureMatrix/FeatureMatrix_", unique(dataset$File), ".txt")
  FeatureMatrix <- do.call(rbind, lapply(FeatureMatrix, fread))
  FeatureMatrix[, ID := paste0(aat, "_", ID)]
  FeatureMatrix <- FeatureMatrix[ID %in% dataset[, ID]]
  FeatureMatrix[, AA := aat]
  setkey(FeatureMatrix, ID)
  FeatureMatrix <- FeatureMatrix[dataset$ID, ]
  FeatureMatrix <- data.frame(FeatureMatrix[, grepl("^X", colnames(FeatureMatrix)) | grepl("DwellTime", colnames(FeatureMatrix)), with = F], row.names = FeatureMatrix[[1]])
  
  LineMat <- melt.data.table(as.data.table(FeatureMatrix, keep.rownames = "ID"), id.vars = "ID")[variable != "DwellTime"]
  LineMat[, x := as.numeric(gsub("^X", "", variable))]
  LineMat[, Blockade := 1 - x/1000]
  LineMat <- LineMat[, .SD[, .(alpha = value/sum(value), Blockade)], ID]
  
  alldists <- readRDS(file.path("./analysis/82.ABFProcessing/SelectedSignals/01.StandardAA/01.SignalDistance", paste0(aat, ".signal.dist.Rds")))
  sigs <- KnnSelecting(signals = dataset, dists = alldists, k = 15, noise = 0)
  
  subtab <- subset.data.frame(FeatureMatrix, row.names(FeatureMatrix) %in% sigs)
  res.ds <- stats::dist(subtab)
  res.hc <- stats::hclust(res.ds, method = "single")
  saveRDS(res.hc, file = paste0("analysis/82.ABFProcessing/SelectedSignals/01.StandardAA/02.SelectedSignals/", aat, ".signal.dist.after.knn15.Rds"))
}, mc.cores = 21)


# KNN = 18 ----


mclapply(AAs, function(aat) {
  if(file.exists(paste0("analysis/82.ABFProcessing/SelectedSignals/01.StandardAA/02.SelectedSignals/", aat, ".signal.dist.after.knn18_1.Rds"))) return(NULL)
  dataset <- fread(paste0("./analysis/82.ABFProcessing/SelectedSignals/", aat, "_signal.txt"))
  dataset[, file_id := as.character(file_id)]
  dataset[, File := gsub("\\.[0-9]", "", file_id)]
  dataset[, ID := paste0(aat, "_", ID)]
  
  FeatureMatrix <- paste0("./analysis/82.ABFProcessing/FeatureMatrix/FeatureMatrix_", unique(dataset$File), ".txt")
  FeatureMatrix <- do.call(rbind, lapply(FeatureMatrix, fread))
  FeatureMatrix[, ID := paste0(aat, "_", ID)]
  FeatureMatrix <- FeatureMatrix[ID %in% dataset[, ID]]
  FeatureMatrix[, AA := aat]
  setkey(FeatureMatrix, ID)
  FeatureMatrix <- FeatureMatrix[dataset$ID, ]
  FeatureMatrix <- data.frame(FeatureMatrix[, grepl("^X", colnames(FeatureMatrix)) | grepl("DwellTime", colnames(FeatureMatrix)), with = F], row.names = FeatureMatrix[[1]])
  
  alldists <- readRDS(file.path("./analysis/82.ABFProcessing/SelectedSignals/01.StandardAA/01.SignalDistance", paste0(aat, ".signal.dist.Rds")))
  sigs <- KnnSelecting(signals = dataset, dists = alldists, k = 18, noise = 1)
  
  subtab <- subset.data.frame(FeatureMatrix, row.names(FeatureMatrix) %in% sigs)
  res.ds <- stats::dist(subtab)
  res.hc <- stats::hclust(res.ds, method = "single")
  saveRDS(res.hc, file = paste0("analysis/82.ABFProcessing/SelectedSignals/01.StandardAA/02.SelectedSignals/", aat, ".signal.dist.after.knn18_1.Rds"))
}, mc.cores = 21)


# KNN = 20 ----


mclapply(AAs, function(aat) {
  if(file.exists(paste0("analysis/82.ABFProcessing/SelectedSignals/01.StandardAA/02.SelectedSignals/", aat, ".signal.dist.after.knn20.Rds"))) return(NULL)
  bgsignal <- fread(paste0("./analysis/82.ABFProcessing/SelectedSignals/", aat, "_background.txt"))
  bgsignal[, file_id := as.character(file_id)]
  bgsignal[, File := gsub("\\.[0-9]", "", file_id)]
  bgsignal[, ID := paste0("Noise_", ID)]
  
  dataset <- fread(paste0("./analysis/82.ABFProcessing/SelectedSignals/", aat, "_signal.txt"))
  dataset[, file_id := as.character(file_id)]
  dataset[, File := gsub("\\.[0-9]", "", file_id)]
  dataset[, ID := paste0(aat, "_", ID)]
  
  FeatureMatrix <- paste0("./analysis/82.ABFProcessing/FeatureMatrix/FeatureMatrix_", unique(dataset$File), ".txt")
  FeatureMatrix <- do.call(rbind, lapply(FeatureMatrix, fread))
  FeatureMatrix[, ID := paste0(aat, "_", ID)]
  FeatureMatrix <- FeatureMatrix[ID %in% dataset[, ID]]
  FeatureMatrix[, AA := aat]
  setkey(FeatureMatrix, ID)
  FeatureMatrix <- FeatureMatrix[dataset$ID, ]
  FeatureMatrix <- data.frame(FeatureMatrix[, grepl("^X", colnames(FeatureMatrix)) | grepl("DwellTime", colnames(FeatureMatrix)), with = F], row.names = FeatureMatrix[[1]])
  
  LineMat <- melt.data.table(as.data.table(FeatureMatrix, keep.rownames = "ID"), id.vars = "ID")[variable != "DwellTime"]
  LineMat[, x := as.numeric(gsub("^X", "", variable))]
  LineMat[, Blockade := 1 - x/1000]
  LineMat <- LineMat[, .SD[, .(alpha = value/sum(value), Blockade)], ID]
  
  alldists <- readRDS(file.path("./analysis/82.ABFProcessing/SelectedSignals/01.StandardAA/01.SignalDistance", paste0(aat, ".signal.dist.Rds")))
  sigs <- KnnSelecting(signals = dataset, dists = alldists, k = 20, noise = 0)
  
  subtab <- subset.data.frame(FeatureMatrix, row.names(FeatureMatrix) %in% sigs)
  res.ds <- stats::dist(subtab)
  res.hc <- stats::hclust(res.ds, method = "single")
  saveRDS(res.hc, file = paste0("analysis/82.ABFProcessing/SelectedSignals/01.StandardAA/02.SelectedSignals/", aat, ".signal.dist.after.knn20.Rds"))
}, mc.cores = 21)



















mclapply(AAs, function(aat) {
  if(file.exists(paste0("analysis/82.ABFProcessing/SelectedSignals/01.StandardAA/02.SelectedSignals/", aat, ".signal.remove.outliers18_1.Rds"))) return(NULL)
  dataset <- fread(paste0("./analysis/82.ABFProcessing/SelectedSignals/", aat, "_signal.txt"))
  dataset[, file_id := as.character(file_id)]
  dataset[, File := gsub("\\.[0-9]", "", file_id)]
  dataset[, ID := paste0(aat, "_", ID)]
  
  res.hc <- readRDS(paste0("analysis/82.ABFProcessing/SelectedSignals/01.StandardAA/02.SelectedSignals/", aat, ".signal.dist.after.knn18_1.Rds"))
  if(aat == "Gln") {
    sigs <- RemoveOutlier(tree = res.hc, n = 60, Keep = 0.95)
  } else {
    sigs <- RemoveOutlier(tree = res.hc, n = 50, Keep = 0.95)
  }
  saveRDS(dataset[ID %in% sigs], file = paste0("analysis/82.ABFProcessing/SelectedSignals/01.StandardAA/02.SelectedSignals/", aat, ".signal.remove.outliers18_1.Rds"))
}, mc.cores = 21)


