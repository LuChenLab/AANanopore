setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(ggplot2)

files <- list.files("./analysis/04.Simulation/03.RepeatedSubsampling", full.names = T)

succMatch <- lapply(files, fread)
succMatch <- mapply(succMatch, FUN = function(x) x[, sum(Match == 1, na.rm = T)])

succMatch <- data.table(Kmers = mapply(function(x) x[1], strsplit(basename(files), "_")), Match = succMatch)

succMatch[, Kmers := factor(Kmers, levels = paste0("kmer", 5:11))]
ggplot(succMatch, aes(x = Kmers, y = Match/100, colour = 1)) + 
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.2) +
  theme_bw(base_size = 15) + 
  labs(y = "Theoretical recovery (%)") + 
  theme(legend.position = "none")
ggsave("./analysis/04.Simulation/03.RepeatedSubsampling/Success_Rate.pdf", width = 5, height = 3.5)


