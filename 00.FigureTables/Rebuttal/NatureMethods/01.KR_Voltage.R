setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(ggplot2)

AABlockade <- data.table(AA = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"), 
                         Blockade = c(0.14718, 0.17148, 0.16516, 0.21409, 0.16, 0.24528, 0.1882, 0.12072, 0.24652, 0.20722, 0.1995, 0.16875, 0.19772, 0.22018, 0.2183, 0.13131, 0.16101, 0.22744, 0.21276, 0.19044))
AABlockade$amino_acid <- plyr::mapvalues(AABlockade$AA, names(Biostrings::AMINO_ACID_CODE), Biostrings::AMINO_ACID_CODE)


meta <- as.data.table(openxlsx::read.xlsx("./data/Voltage/标准品混合物1.xlsx"))
colnames(meta) <- c("file_name", "amino_acid", "voltage", "concentration", "start_time", "end_time", "baseline_mean")
meta[, start_time := as.numeric(start_time)]
meta[, end_time := as.numeric(end_time)]

files <- list.files("./analysis/21.ABFProcessing/04.Voltage/RawSignal", full.names = T)
Sigs <- lapply(files, function(x) {
  y <- data.table(AA = mapply(function(x) x[2], strsplit(x, "_")), Voltage = gsub("mV.txt", "", mapply(function(x) x[3], strsplit(x, "_"))), fread(x))
  starttime <- y[, min(StartTime)]
  endtime <- starttime + 120 # 2 minutes
  y[StartTime > starttime & EndTime < endtime]
})

Sigs <- do.call(rbind, Sigs)
Sigs[, Voltage := paste(Voltage, "(mV)")]
Sigs[, Voltage := factor(Voltage, levels = c("50 (mV)", "75 (mV)", "100 (mV)"))]

ggplot(Sigs[DwellTime < .1], aes(x = Blockade, y = DwellTime * 1000)) + 
  geom_point(size = .4) + 
  scale_y_log10() + 
  facet_grid(Voltage ~ AA) + 
  labs(y = "DwellTime (ms)") + 
  theme_light(base_size = 15)



