setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(ggplot2)

mean2 <- function(x, q1 = 0.25, q2 = 0.75) {
  q <- quantile(x, c(q1, q2))
  mean(x[x >= min(q) & x <= max(q)])
}

meta <- do.call(rbind, lapply(list.files("./analysis/01.AASignalRecognition/Version9/02.RangesL0L1/", ".txt", full.names = TRUE, recursive = T), fread))
meta <- meta[!file_name %in% c("21205012", "21205021", "21201008", "21303003")]

L0TimeRatio <- mcmapply(seq_len(nrow(meta)), FUN = function(i) {
  abf <- readRDS(grep(as.character(meta[i, file_name]), list.files("./analysis/01.AASignalRecognition/Version9/01.SignalPolish/", full.names = T, recursive = T), value = T))
  abf[, mean(Sm > meta[i, MinL0] & Sm < meta[i, MaxL0])]
}, mc.cores = 10)
L0TimeRatio <- data.table(file_name = meta$file_name, L0TimeRatio)
L0TimeRatio[, file_name := as.character(file_name)]

L1TimeRatio <- mcmapply(seq_len(nrow(meta)), FUN = function(i) {
  abf <- readRDS(grep(as.character(meta[i, file_name]), list.files("./analysis/01.AASignalRecognition/Version9/01.SignalPolish/", full.names = T, recursive = T), value = T))
  abf[, mean(Sm > meta[i, MinL1] & Sm < meta[i, MaxL1])]
}, mc.cores = 10)
L1TimeRatio <- data.table(file_name = meta$file_name, L1TimeRatio)
L1TimeRatio[, file_name := as.character(file_name)]

aa_L2 <- as.data.table(read.xlsx("./data/aa_L2.xlsx"))
meta <- merge(meta, aa_L2, by.x = "amino_acid", by.y = "name")

L2TimeRatio <- mcmapply(seq_len(nrow(meta)), FUN = function(i) {
  abf <- readRDS(grep(as.character(meta[i, file_name]), list.files("./analysis/01.AASignalRecognition/Version9/01.SignalPolish/", full.names = T, recursive = T), value = T))
  BaseLineMean <- abf[Sm > meta[i, MinL0] & Sm < meta[i, MaxL0], mean2(pA)]
  L2 <- BaseLineMean * (1 - c(meta[i, L2min], meta[i, L2max]))
  abf[, mean(Sm > min(L2) & Sm < max(L2))]
}, mc.cores = 10)
L2TimeRatio <- data.table(file_name = meta$file_name, L2TimeRatio)
L2TimeRatio[, file_name := as.character(file_name)]

meta$TimeRatio <- L0TimeRatio$L0TimeRatio + L1TimeRatio$L1TimeRatio + L2TimeRatio$L2TimeRatio


meta <- meta[, .(file_name, amino_acid, concentration, start_time, end_time, TimeRatio)]
meta[, file_name := as.character(file_name)]

files <- list.files("./analysis/01.AASignalRecognition/Version9/04.FinalSignal", pattern = "_Sigs.txt", recursive = TRUE, full.names = T)
AAs <- lapply(files, fread)
names(AAs) <- gsub("_Sigs.txt", "", basename(files))
for(i in seq_along(AAs)) AAs[[i]] <- data.table(file_name = names(AAs)[i], AAs[[i]])
AA <- do.call(rbind, AAs)

AA <- merge(meta, AA, by = "file_name")

AA <- AA[!is.na(Blockade) & Outer == 0 & AreaRatio_L1 > 0.1]

AA[, .N, amino_acid]
AA[, aa := plyr::mapvalues(amino_acid, AMINO_ACID_CODE, names(AMINO_ACID_CODE))]

AA[aa %in% c("E", "D", "H", "R", "K"), Class := "charged"]
AA[aa %in% c("L", "I", "M", "V", "A", "F", "G", "W", "P"), Class := "polar"]
AA[aa %in% c("S", "N", "Q", "T", "Y", "C"), Class := "nonpolar"]


AA[AreaRatio_L2 > 0, Count := 2]
AA[is.na(Count),  Count := 1]

PercentL2 <- AA[, .(PercentL2 = mean(Count == 2) * 100), aa]

AA_N <- merge(AA[, .(N = sum(Count)), file_name], meta, by = "file_name")
AA_N[, N2 := N/concentration]
AA_N[, N3 := N2 / ((end_time - start_time) * TimeRatio)]
AA_N[, aa := plyr::mapvalues(amino_acid, AMINO_ACID_CODE, names(AMINO_ACID_CODE))]

od <- AA_N[, median(N3), aa][order(V1), aa]
AA_N[, aa := factor(aa, levels = rev(od))]

AA_N[aa %in% c("E", "D", "H", "R", "K"), Class := "charged"]
AA_N[aa %in% c("L", "I", "M", "V", "A", "F", "G", "W", "P"), Class := "polar"]
AA_N[aa %in% c("S", "N", "Q", "T", "Y", "C"), Class := "nonpolar"]
AA_N[, Class := factor(Class, levels = c("charged", "polar", "nonpolar"))]

ggplot(AA_N[aa != "C"], aes(x = aa, y = N3, colour = Class)) + 
  geom_boxplot() + 
  geom_point() + 
  theme_bw(base_size = 15) + 
  labs(y = expression(Velocity~of~amono~acid~"("*Min^"-1"*")")) + 
  scale_colour_brewer(palette = "Dark2") + 
  theme(legend.position = "top", axis.title.x = element_blank())
ggsave("./analysis/01.AASignalRecognition/Version9/05.plot/AA_velocity.pdf", width = 5, height = 4)




ggplot(AA_N, aes(x = aa, y = N3, colour = Class)) + 
  geom_boxplot() + 
  geom_point() + 
  theme_bw(base_size = 15) + 
  labs(y = expression(Velocity~of~amono~acid~"("*Min^"-1"*")")) + 
  scale_colour_brewer(palette = "Dark2") + 
  theme(legend.position = "top", axis.title.x = element_blank())
ggsave("./analysis/01.AASignalRecognition/Version9/05.plot/AA_velocity_with_Cys.pdf", width = 5, height = 4)


# AA_N2 <- AA_N[aa != "C", mean(N3), by = c("Class", "aa")]
# od <- AA_N2[, median(V1), Class][order(V1), Class]
# AA_N2[, Class := factor(Class, levels = rev(od))]

library(ggpubr)
my_comp <- list(c("charged", "polar"), c("charged", "nonpolar"), c("polar", "nonpolar"))

AA_N2 <- AA_N[aa != "C", mean(N3), by = c("Class", "aa")]

ggplot(AA_N2, aes(x = Class, y = V1, colour = Class)) + 
  geom_boxplot() + 
  geom_point() + 
  theme_bw(base_size = 15) + 
  labs(y = expression(Velocity~of~amono~acid~"("*Min^"-1"*")")) + 
  scale_colour_brewer(palette = "Dark2") + 
  theme(legend.position = "none", axis.title.x = element_blank()) + 
  stat_compare_means(comparisons = my_comp)
ggsave("./analysis/01.AASignalRecognition/Version9/05.plot/AA_velocity_Compare.pdf", width = 4, height = 4)


AA_N2 <- AA_N[, mean(N3), by = c("Class", "aa")]

ggplot(AA_N2, aes(x = Class, y = V1, colour = Class)) + 
  geom_boxplot() + 
  geom_point() + 
  theme_bw(base_size = 15) + 
  labs(y = expression(Velocity~of~amono~acid~"("*Min^"-1"*")")) + 
  scale_colour_brewer(palette = "Dark2") + 
  theme(legend.position = "none", axis.title.x = element_blank()) + 
  stat_compare_means(comparisons = my_comp)
ggsave("./analysis/01.AASignalRecognition/Version9/05.plot/AA_velocity_Compare_with_Cys.pdf", width = 4, height = 4)


A_V <- data.table(V = c(87.8, 188.2, 120.1, 115.4, 105.4, 140.9, 145.1, 59.9, 156.3, 166.1, 168, 172.7, 165.2, 189.7, 123.3, 91.7, 118.3, 227.9, 191.2, 138.8), 
                  aa = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"))
A_N_V <- merge(AA_N2, A_V, by = "aa")




ggplot(A_N_V[aa != "C"], aes(x = V, y = V1, colour = Class)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = F) + 
  stat_cor() + 
  theme_bw(base_size = 15) + 
  labs(y = expression(Velocity~of~amono~acid~"("*Min^"-1"*")"), 
       x = expression("Volume of amino acid (" *10^"-3"*nm^3*")")) + 
  scale_colour_brewer(palette = "Dark2") + 
  theme(legend.position = "top")
ggsave("./analysis/01.AASignalRecognition/Version9/05.plot/AA_velocity_Vs_Volume.pdf", width = 5, height = 5)


ggplot(A_N_V, aes(x = V, y = V1, colour = Class)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = F) + 
  stat_cor() + 
  theme_bw(base_size = 15) + 
  labs(y = expression(Velocity~of~amono~acid~"("*Min^"-1"*")"), 
       x = expression("Volume of amino acid (" *10^"-3"*nm^3*")")) + 
  scale_colour_brewer(palette = "Dark2") + 
  theme(legend.position = "top")
ggsave("./analysis/01.AASignalRecognition/Version9/05.plot/AA_velocity_Vs_Volume_with_Cys.pdf", width = 5, height = 5)



