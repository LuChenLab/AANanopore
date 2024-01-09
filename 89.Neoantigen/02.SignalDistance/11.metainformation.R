setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)

meta <- as.data.table(openxlsx::read.xlsx("data/MetaInfomation/neoantigen水解20220314.xlsx", sheet = 1))
colnames(meta) <- c("file_name", "date", "amino_acid", "concentration", "start_time", "end_time", "base_line")
meta[, file_name := gsub(".abf", "", file_name)]
setkey(meta, file_name, start_time)
meta1 <- meta[, .SD[, .(amino_acid, start_time, end_time, file_id = paste0(file_name, ".", seq_len(.N)))], file_name]

meta <- as.data.table(openxlsx::read.xlsx("data/MetaInfomation/neoantigen水解20220406.xlsx", sheet = 1))[, 1:7]
colnames(meta) <- c("file_name", "date", "amino_acid", "concentration", "start_time", "end_time", "base_line")
meta[, file_name := gsub(".abf", "", file_name)]
meta <- meta[!grepl("mV", amino_acid)]
setkey(meta, file_name, start_time)
meta2 <- meta[, .SD[, .(amino_acid, start_time, end_time, file_id = paste0(file_name, ".", seq_len(.N)))], file_name]

meta <- as.data.table(openxlsx::read.xlsx("data/WangZichun/20231028/20231020起实验数据记录.xlsx", sheet = 1))[, 1:7]
colnames(meta) <- c("file_name", "start_time", "end_time", "amino_acid", "file_id", "purpose", "sample")
meta[, file_name := gsub(".abf", "", file_name)]
meta <- meta[grepl("肿瘤新抗原肽", purpose)]
setkey(meta, file_name, start_time)
meta3 <- meta[, .SD[, .(amino_acid, start_time, end_time, file_id = paste0(file_name, ".", seq_len(.N)))], file_name]

metas <- rbind(meta1, meta2, meta3)
metas[grepl("lank", amino_acid), amino_acid := NA]

metas[, sig_file := paste0("./analysis/81.ABFProcessing/ABF/ABF_", file_name, ".Rds")]
metas <- metas[file.exists(sig_file)]

metas[, sig_file := paste0("./analysis/81.ABFProcessing/RawSignal/RawSignal_", file_name, ".txt")]
metas <- metas[file.exists(sig_file)]

metas[, sig_file := paste0("./analysis/89.Neoantigen/01.SelectedL0/", file_id, ".MainL0.txt")]
metas <- metas[file.exists(sig_file)]
metas[, amino_acid := gsub("KLVVVGAGGV,", "", amino_acid)]
metas[, amino_acid := gsub("KLVVVGADGV,", "", amino_acid)]
metas[, amino_acid := gsub("VLLGVKLSGV,", "", amino_acid)]
metas[, amino_acid := gsub("VLLGVKLFGV,", "", amino_acid)]
metas[, amino_acid := gsub("VGDAL", "VGALD", amino_acid)]
metas[, amino_acid := gsub(" ", "", amino_acid)]

metas[, sig_file := NULL]

openxlsx::write.xlsx(metas, "./analysis/89.Neoantigen/02.SignalDistance/meta.xlsx")














