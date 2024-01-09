setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)

meta <- as.data.table(openxlsx::read.xlsx("data/MetaInfomation/整理数据(20230831).xlsx", sheet = 2, cols = 1:9))
colnames(meta) <- c("file_name", "start_time", "end_time", "amino_acid", "file_id", "purpose", "sample", "base_line", "note")
meta[file_name == "20230330 _0005", file_name := "20230330_0005"]
setkey(meta, file_name, start_time)
meta1 <- meta[, .SD[, .(amino_acid, start_time, end_time, file_id = paste0(file_name, ".", seq_len(.N)))], file_name]

meta <- as.data.table(openxlsx::read.xlsx("data/MetaInfomation/20231020 氨肽酶超滤.xlsx", sheet = 1))
colnames(meta) <- c("file_name", "start_time", "end_time", "amino_acid", "file_id", "purpose", "sample", "base_line", "note")
setkey(meta, file_name, start_time)
meta2 <- meta[, .SD[, .(amino_acid, start_time, end_time, file_id = paste0(file_name, ".", seq_len(.N)))], file_name]

meta <- as.data.table(openxlsx::read.xlsx("data/ChenShanchuan/20231024/20231024_data.xlsx", sheet = 1))
meta <- meta[实验目的 == "阿尔兹海默症突变位点" & grepl("20231024", 文件名)]
colnames(meta) <- c("file_name", "start_time", "end_time", "amino_acid", "file_id", "purpose", "sample", "baseline", "note")
meta[, file_name := gsub(".abf", "", file_name)]
setkey(meta, file_name)
meta3 <- meta[, .SD[, .(amino_acid, start_time, end_time, file_id = paste0(file_name, ".", seq_len(.N)))], file_name]

meta <- as.data.table(openxlsx::read.xlsx("data/WangZichun/20231025/20231020起实验数据记录.xlsx", sheet = 1))
meta <- meta[实验目的 == "阿尔兹海默症突变位点" & grepl("20231025", 文件名)]
colnames(meta) <- c("file_name", "start_time", "end_time", "amino_acid", "file_id", "purpose", "sample", "baseline", "note")
meta[, file_name := gsub(".abf", "", file_name)]
setkey(meta, file_name)
meta4 <- meta[, .SD[, .(amino_acid, start_time, end_time, file_id = paste0(file_name, ".", seq_len(.N)))], file_name]

meta <- as.data.table(openxlsx::read.xlsx("data/ChenShanchuan/20231026P1/20231026_P1_2.xlsx", sheet = 1))
meta <- meta[实验目的 == "阿尔兹海默症突变位点" & grepl("20231026", 文件名)]
colnames(meta) <- c("file_name", "start_time", "end_time", "amino_acid", "file_id", "purpose", "sample", "baseline", "note")
meta[, file_name := gsub(".abf", "", file_name)]
setkey(meta, file_name)
meta5 <- meta[, .SD[, .(amino_acid, start_time, end_time, file_id = paste0(file_name, ".", seq_len(.N)))], file_name]

metas <- rbind(meta1, meta2, meta3, meta4, meta5)
metas[amino_acid == "AP", amino_acid := NA]

metas[, sig_file := paste0("./analysis/81.ABFProcessing/ABF/ABF_", file_name, ".Rds")]
metas <- metas[file.exists(sig_file)]

metas[, sig_file := paste0("./analysis/81.ABFProcessing/RawSignal/RawSignal_", file_name, ".txt")]
metas <- metas[file.exists(sig_file)]

metas[, sig_file := NULL]
sig_file <- mapply(metas[, file_id], FUN = function(x) grep(".MainL0.txt", list.files("./analysis/84.ADpolypeptide", x, full.names = T, recursive = T), value = T))
metas <- merge(metas, data.table(file_id = names(unlist(sig_file)), sig_file = unlist(sig_file)), all.x = TRUE)

metas <- metas[!is.na(sig_file)]

openxlsx::write.xlsx(metas, "./analysis/84.ADpolypeptide/60.MergeAll/metainformation.xlsx")
