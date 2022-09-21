library(ggrepel)
setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(ggplot2)

mean2 <- function(x, q1 = 0.25, q2 = 0.75) {
  q <- quantile(x, c(q1, q2))
  mean(x[x >= min(q) & x <= max(q)])
}

# c("21205012", "21205021", "21201008", "21303003")

meta <- do.call(rbind, lapply(list.files("./analysis/01.AASignalRecognition/Version9/02.RangesL0L1/", ".txt", full.names = TRUE, recursive = T), fread))
meta <- meta[, .(file_name, amino_acid, concentration, start_time, end_time, TimeRatio)]
meta[, file_name := as.character(file_name)]
meta <- meta[!file_name %in% c("21205012", "21205021", "21201008", "21303003")]

files <- list.files("./analysis/01.AASignalRecognition/Version9/04.FinalSignal", pattern = "_Sigs.txt", recursive = TRUE, full.names = T)
AAs <- lapply(files, fread)
names(AAs) <- gsub("_Sigs.txt", "", basename(files))
for(i in seq_along(AAs)) AAs[[i]] <- data.table(file_name = names(AAs)[i], AAs[[i]])
AA <- do.call(rbind, AAs)
AA <- AA[!file_name %in% c("21205012", "21205021", "21201008", "21303003")]

AA <- merge(meta, AA, by = "file_name")
AA <- AA[!is.na(Blockade) & Outer == 0 & AreaRatio_L1 > 0.1]
AA[, .N, amino_acid]
AA[, aa := plyr::mapvalues(amino_acid, AMINO_ACID_CODE, names(AMINO_ACID_CODE))]

AA[aa %in% c("E", "D", "H", "R", "K"), Class := "charged"]
AA[aa %in% c("L", "I", "M", "V", "A", "F", "G", "W", "P"), Class := "polar"]
AA[aa %in% c("S", "N", "Q", "T", "Y", "C"), Class := "nonpolar"]

ggplot(data = AA[AreaRatio_L1 == 1 & DwellTime > 0.75 & SignalSD < 3], 
       mapping = aes(x = Blockade, y = DwellTime, colour = aa)) + 
  geom_point(size = 0.1) + 
  scale_x_continuous(limits = c(0.1, 0.26)) + 
  facet_wrap(~ Class, ncol = 1)

library(ggrepel)
aaBlockade <- lapply(AA[, unique(aa)], FUN = function(a) {
  den <- AA[DwellTime > 0.75 & SignalSD < 3.5 & aa == a, density(Blockade, adjust = 2)]
  SD <- AA[DwellTime > 0.75 & SignalSD < 3.5 & aa == a, sd(Blockade)]
  data.table(aa = a, x = den$x[which.max(den$y)], y = max(den$y), SD)
})
aaBlockade <- do.call(rbind, aaBlockade)
aaBlockade[aa %in% c("E", "D", "H", "R", "K"), Class := "charged"]
aaBlockade[aa %in% c("L", "I", "M", "V", "A", "F", "G", "W", "P"), Class := "polar"]
aaBlockade[aa %in% c("S", "N", "Q", "T", "Y", "C"), Class := "nonpolar"]

ggplot(data = AA[DwellTime > 0.75 & SignalSD < 3.5], 
       mapping = aes(x = Blockade, colour = aa)) + 
  geom_line(stat = "Density", adjust = 3) + 
  scale_x_continuous(limits = c(0.1, 0.26)) +
  facet_wrap(~ Class, ncol = 1, scales = "free_y", strip.position = "right") + 
  theme_bw(base_size = 15) + 
  theme(legend.position = "none") + 
  geom_text_repel(data = aaBlockade, aes(x = x, y = y, label = aa)) + 
  scale_colour_manual(values = ggsci::pal_igv()(20))
ggsave("./analysis/01.AASignalRecognition/Version9/05.plot/AA_Blockade_Density_with_Cys.pdf", width = 6, height = 6)

ggplot(data = AA[DwellTime > 0.75 & SignalSD < 3.5 & aa != "C"], 
       mapping = aes(x = Blockade, colour = aa)) + 
  geom_line(stat = "Density", adjust = 3) + 
  scale_x_continuous(limits = c(0.1, 0.26)) +
  facet_wrap(~ Class, ncol = 1, scales = "free_y", strip.position = "right") + 
  theme_bw(base_size = 15) + 
  theme(legend.position = "none") + 
  geom_text_repel(data = aaBlockade[aa != "C"], aes(x = x, y = y, label = aa)) + 
  scale_colour_manual(values = ggsci::pal_igv()(20))
ggsave("./analysis/01.AASignalRecognition/Version9/05.plot/AA_Blockade_Density.pdf", width = 6, height = 6)



ggplot(data = AA[DwellTime > 0.75 & SignalSD < 3.5], 
       mapping = aes(x = Blockade, colour = aa)) + 
  geom_line(stat = "Density", adjust = 4) + 
  scale_x_continuous(limits = c(0.1, 0.26)) +
  facet_wrap(~ Class, ncol = 1, scales = "free_y", strip.position = "right") + 
  theme_bw(base_size = 15) + 
  theme(legend.position = "none") + 
  geom_text_repel(data = aaBlockade, aes(x = x, y = y, label = aa)) + 
  scale_colour_manual(values = ggsci::pal_igv()(20))


AAs <- unique(AA$aa)
AAs <- setdiff(AAs, "C")

AreaDrop <- lapply(AAs, function(x) {
  A1 <- x
  A2 <- setdiff(AAs, x)
  area <- mapply(A2, FUN = function(y) {
    r1 <- AA[aa == x, Blockade]
    r2 <- AA[aa == y, Blockade]
    
    Fn1 <- ecdf(r1) # a *function*
    Fn2 <- ecdf(r2) # a *function*
    
    d1 <- density(r1, from = min(c(r1, r2)), to = max(c(r1, r2)), n = 2048, adjust = 1)
    d2 <- density(r2, from = min(c(r1, r2)), to = max(c(r1, r2)), n = 2048, adjust = 1)
    identical(d1$x, d2$x)
    ds <- data.table(p = seq_along(d1$x), x = d1$x, y1 = d1$y, y2 = d2$y)
    ds[, d := abs(y1 - y2)]
    p1 <- ds[, which.max(y1)]
    p2 <- ds[, which.max(y2)]
    solution <- ds[p1:p2, .SD[which.min(d), x]]
    min(Fn1(solution), 1 - Fn1(solution))
  })
  area <- c(1, area)
  names(area)[1] <- A1
  return(area)
})
names(AreaDrop) <- AAs

AreaDropTab <- do.call(rbind, lapply(AreaDrop, function(x) as.data.table(as.data.frame(x), keep.rownames = "A2")))
AreaDropTab$A1 <- rep(names(AreaDrop), mapply(length, AreaDrop))
AreaDropMat <- dcast.data.table(AreaDropTab, A1 ~ A2, value.var = "x")
AreaDropTab[x == 1, x := NA]
AreaDropTab[, Percent := x * 100]

ggplot(AreaDropTab, aes(x = A1, y = A2, fill = Percent)) + 
  geom_tile() + 
  geom_text(aes(label = round(Percent, 2))) +
  scale_fill_gradient2(low = "red", high = "red", midpoint = 0.0001) + 
  theme_classic(base_size = 15) + 
  theme(axis.title = element_blank())



library(RColorBrewer)
library(pheatmap)
AreaDropMat <- data.frame(AreaDropMat[, -1], row.names = AreaDropMat[[1]])
AreaDropMat <- AreaDropMat * 100

fold <- (100 - max(AreaDropMat[AreaDropMat < 100]))/max(AreaDropMat[AreaDropMat < 100])
n <- ceiling(max(AreaDropMat[AreaDropMat < 100]) * 100/1)

color1 = colorRampPalette(brewer.pal(n = 9, name = "OrRd")[1:7])(n)
color2 = colorRampPalette(brewer.pal(n = 9, name = "OrRd")[8])(n * fold)




pdf("./analysis/01.AASignalRecognition/Version9/05.plot/Similarity.pdf", width = 8, height = 8)
pheatmap(AreaDropMat, 
         color = c("grey98", color1, color2), 
         cluster_rows = F, 
         angle_col = 0, 
         display_numbers = round(AreaDropMat, 2), 
         cluster_cols = F, 
         legend = FALSE)
dev.off()




pdf("./analysis/01.AASignalRecognition/Version9/05.plot/Similarity_2.pdf", width = 8, height = 8)
pheatmap(AreaDropMat, 
         color = c("grey98", color1, color2), 
         cluster_rows = T, 
         angle_col = 0, 
         display_numbers = round(AreaDropMat, 2), 
         cluster_cols = T, 
         legend = FALSE)
dev.off()





A_V <- data.table(V = c(87.8, 188.2, 120.1, 115.4, 105.4, 140.9, 145.1, 59.9, 156.3, 166.1, 168, 172.7, 165.2, 189.7, 123.3, 91.7, 118.3, 227.9, 191.2, 138.8), 
                  aa = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"))

library(ggpubr)
aaBlockade <- merge(aaBlockade, A_V, by = "aa")
aaBlockade[, Class := factor(Class, levels = c("charged", "polar", "nonpolar"))]
ggplot(aaBlockade, aes(x = V, y = x)) + 
  geom_pointrange(mapping = aes(ymin = x - SD, ymax = x + SD, colour = Class)) + 
  geom_smooth(method = "lm") + 
  geom_text_repel(aes(label = aa)) + 
  stat_cor() + 
  theme_bw(base_size = 15) + 
  labs(y = "Blockade", x = expression("Volume of amino acid (" *10^"-3"*nm^3*")")) + 
  theme(legend.position = "top") + 
  scale_colour_brewer(palette = "Dark2")
ggsave("./analysis/01.AASignalRecognition/Version9/05.plot/Blockade_Vs_Volume.pdf", width = 4, height = 4)


ggplot(aaBlockade, aes(x = V, y = x)) + 
  geom_pointrange(mapping = aes(ymin = x - SD, ymax = x + SD, colour = Class)) + 
  geom_smooth(method = "lm") + 
  geom_text_repel(aes(label = aa)) + 
  stat_cor(label.y.npc = "bottom", label.x.npc = "centre") + 
  theme_bw(base_size = 15) + 
  labs(y = "Blockade", x = expression("Volume of amino acid (" *10^"-3"*nm^3*")")) + 
  theme(legend.position = "top") + 
  scale_colour_brewer(palette = "Dark2")
ggsave("./analysis/01.AASignalRecognition/Version9/05.plot/Blockade_Vs_Volume2.pdf", width = 4.2, height = 4.5)


ggplot(aaBlockade[!aa %in% c("K", "R", "D", "E", "P", "H")], aes(x = V, y = x)) + 
  geom_pointrange(mapping = aes(ymin = x - SD, ymax = x + SD, colour = Class)) + 
  geom_smooth(method = "lm", se = F, color = "grey") + 
  geom_text_repel(aes(label = aa)) + 
  stat_cor(label.y.npc = "bottom", label.x.npc = "centre") + 
  theme_bw(base_size = 15) + 
  labs(y = "Blockade", x = expression("Volume of amino acid (" *10^"-3"*nm^3*")")) + 
  theme(legend.position = "top") + 
  scale_colour_manual(values = RColorBrewer::brewer.pal(n = 3, name = "Dark2")[2:3])
ggsave("./analysis/01.AASignalRecognition/Version9/05.plot/Blockade_Vs_Volume4_noKRDEPH_for_lm.pdf", width = 4.2, height = 4.5)




ggplot(aaBlockade, aes(x = V, y = x)) + 
  geom_pointrange(mapping = aes(ymin = x - SD, ymax = x + SD, colour = Class)) + 
  geom_smooth(method = "lm", se = F, color = "grey") + 
  geom_text_repel(aes(label = aa)) + 
  stat_cor(label.y.npc = "bottom", label.x.npc = "centre") + 
  theme_bw(base_size = 15) + 
  labs(y = "Blockade", x = expression("Volume of amino acid (" *10^"-3"*nm^3*")")) + 
  theme(legend.position = "top") + 
  scale_colour_brewer(palette = "Dark2")
ggsave("./analysis/01.AASignalRecognition/Version9/05.plot/Blockade_Vs_Volume3.pdf", width = 4.2, height = 4.5)








ZM <- na.omit(as.data.table(read.xlsx("./data/blockade_vs_volume_gauss.xlsx")))
colnames(ZM) <- c("aa", "x", "SD", "N", "n", "V")
ZM[, V := as.numeric(V)]
ZM[aa %in% c("E", "D", "H", "R", "K"), Class := "charged"]
ZM[aa %in% c("L", "I", "M", "V", "A", "F", "G", "W", "P"), Class := "polar"]
ZM[aa %in% c("S", "N", "Q", "T", "Y", "C"), Class := "nonpolar"]

Mat <- rbind(data.table(Method = "TC", aaBlockade), data.table(Method = "ZM", ZM), fill = TRUE)

ggplot(Mat, aes(x = V, y = x, shape = Method)) + 
  geom_pointrange(mapping = aes(ymin = x - SD, ymax = x + SD, colour = Class)) + 
  geom_smooth(method = "lm", se = F) + 
  geom_text_repel(aes(label = aa)) + 
  stat_cor(label.y.npc = "bottom", label.x.npc = "centre") + 
  theme_bw(base_size = 15) + 
  labs(y = "Blockade", x = expression("Volume of amino acid (" *10^"-3"*nm^3*")")) + 
  theme(legend.position = "top") + 
  scale_colour_brewer(palette = "Dark2")










