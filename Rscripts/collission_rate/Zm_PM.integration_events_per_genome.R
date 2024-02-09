library(dplyr)
library(tidyverse)
library(ggridges)



cell_metrics <- read_delim("/Users/feilab/Desktop/Cell_metrics.zm_proso.tsv", delim='\t')

final <- cell_metrics %>% 
    dplyr::select(V4, plate,total_insertions, assigned_genome, Pm, Zm) %>%  
    pivot_longer(cols = c(Pm,Zm), names_to = "insertions_per_genome", values_to = "insertions")



insertions_per_cell_per_genome <- ggplot(final, aes(x = insertions_per_genome, y = log(insertions))) + 
  geom_violin() + 
  geom_boxplot(width=.05) + 
  geom_line(aes(group=V4), alpha = .1, color = "gray") + 
  facet_grid(plate~assigned_genome) + 
  xlab("Genome") + ylab("Log2(Number of Insertions)") + 
  theme_bw()

insertions_per_cell_per_genome

ggsave("/Users/feilab/Projects/05.comparative_single_cell/Mendieta_et_al_comparative_single_cell/imgs/collision_rate/Zm.proso.insertions_per_genome_per_cell.png", insertions_per_cell_per_genome, height = 10,  width = 8, units = "in")


?ggsave

