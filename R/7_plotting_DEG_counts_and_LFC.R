library(ggplot2)
library(tidyverse)
library(here)

mature.DE.list = readRDS(file = here("Rout/MA.clust.mature.DE.SAvsPBS.list.13.10.2024.RDS") )

names(mature.DE.list)[which(names(mature.DE.list) == "Granulocyte" )] = "Granulocyte 1" 
names(mature.DE.list)[which(names(mature.DE.list) == "Neutrophil" )] = "Granulocyte 2" 
names(mature.DE.list)[which(names(mature.DE.list) == "Cytotoxic T&NK cells")] = "Cytotoxic T and NK cells"
names(mature.DE.list)[which(names(mature.DE.list) == "Th2 CD4 T cell")] = "Th2 CD4+ T cell"
    
HSC.DE.list = readRDS(file = here("Rout/MA.clust.HSC.DE.SAvsPBS.list.14.11.2024.RDS") )

names(HSC.DE.list)[which(names(HSC.DE.list) == "erythroid megakaryocyte progenitor")] = "Erythroid megakaryocyte progenitor"


mature.all.DE = do.call(rbind, lapply(names(mature.DE.list), function(x) {
    dt = mature.DE.list[[x]]
    dt$celltype = x
    return(dt) } ))

HSC.all.DE = do.call(rbind, lapply(names(HSC.DE.list), function(x) {
    dt = HSC.DE.list[[x]]
    dt$celltype = x
    return(dt) } ))


manual.annot.all.DE = rbind(HSC.all.DE, mature.all.DE)

pp = manual.annot.all.DE %>% 
    mutate(celltype = factor(celltype,
                             levels = c("CLP",
                                        "Pro B cell",
                                        "B cell", 
                                        "Th2 CD4+ T cell", 
                                        "Cytotoxic T and NK cells",
                                        "Erythroid megakaryocyte progenitor",
                                        "Erythroblast",
                                        "GMP CMP",
                                        "GMP",
                                        "Proliferating GMP",
                                        "Granulocyte progenitor",
                                        "Early granulocyte", 
                                        "Granulocyte 1",
                                        "Granulocyte 2",
                                        "Activated monocyte",
                                        "Activated macrophage", 
                                        "Activated DC and macrophage", 
                                        "Dendritic cell",
                                        "MC basophil") %>% rev() ) ) %>%
    ggplot((aes(x = log2FoldChange, y = celltype, color = log2FoldChange, alpha = -log(padj) ) )) + 
    geom_jitter(height = 0.1) + 
    theme_classic() +
    scale_color_gradient2(low = "dodgerblue4", mid = "gray80", high = "darkred", name = "log2FC") +
    scale_alpha_continuous(range = c(0.1, 1), limits = c(0,300), breaks = c(0, 100, 200, 300)) +
    theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray80", linetype = 1, linewidth = 0.3),
        axis.title.y = element_blank()) # +
    guides(alpha = "none")


ggsave(plot = pp, filename = file.path("figures", "DEG_scatterplot_all_celltypes.pdf"),
       width = 4.5, height = 3)

ggsave(plot = pp, filename = file.path("figures", "DEG_scatterplot_all_celltypes.png"),
       width = 4.5, height = 3, dpi = 600)
