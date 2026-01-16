# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(ggh4x)
library(patchwork)
library(ragg)


# Set paths --------------------------------------------------------------------
data_dir <- file.path("results", "scRNAseq", "02_ORA")
res_dir <- file.path("results", "scRNAseq", "03_plots")


# Set ggplot2 theme ------------------------------------------------------------
theme_set(
  theme_bw() +
    theme(
      legend.position = "none",
      strip.text = element_text(size = rel(1.1)),
      axis.text = element_text(color = "black", size = rel(0.9))))


# Input data -------------------------------------------------------------------
## GO ORA results with gene log2FC values 
go_lipid_revigo_imp_2genetable_lfc <- read_tsv(file.path(data_dir, "go_lipid_revigo_imp_degs_lfc.txt"))

## Calculate average log2FC per dataset 
go_lipid_revigo_imp_2genetable_lfc_avg <- go_lipid_revigo_imp_2genetable_lfc %>% 
  dplyr::filter(Tissue %in% c("Lung", "Airway")) %>%
  mutate(celltype_tidy = str_remove_all(celltype, "^[A-Z0-9]*_|_[A-Z0-9]*$")) %>% 
  dplyr::select(celltype_tidy, avg_log2FC, status, dataset) %>% 
  distinct() %>% 
  group_by(celltype_tidy, status, dataset) %>% 
  summarise(
    n = n(),
    avg_avg_log2FC = mean(avg_log2FC)) %>% 
  ungroup() 

go_lipid_revigo_imp_2genetable_lfc_avg %>% 
  group_by(celltype_tidy, status) %>% 
  summarise(avg_avg_log2FC = mean(avg_avg_log2FC)) %>% 
  ungroup() 


# Plot -------------------------------------------------------------------------
go_lipid_revigo_imp_2genetable_lfc_avg %>% 
  mutate(
    abs_avg_avg_log2FC = abs(avg_avg_log2FC),
    celltype_tidy = fct_reorder(celltype_tidy, abs_avg_avg_log2FC, .fun = mean, .desc = FALSE),
    status = fct_rev(status)) %>% 
  ggplot(aes(x = avg_avg_log2FC, y = celltype_tidy, color = status)) +
  geom_point() +
  stat_summary(
    fun.y = mean,
    fun.ymin = function(x) mean(x) - sd(x), 
    fun.ymax = function(x) mean(x) + sd(x), 
    geom = "pointrange", alpha = 0.8) +
  facet_wrap2(
    vars(status), 
    scales = "free_x", 
    labeller = as_labeller(c("Up" = "Upregulated", "Down" = "Downregulated"))) +
  facetted_pos_scales(x = list(
    scale_x_continuous(breaks = c(1, 2.5, 4)),
    scale_x_continuous(limits = c(-3, -1), breaks = seq(-3, -1, 1)))) +
  scale_color_manual(values = c("#f4a261", "#0a9396")) + 
  labs(x = "Average log2FC", y = NULL) 

ggsave(file.path(res_dir, "log2FC_by_celltype.png"), device = agg_png, dpi = 300, width = 6, height = 8, units = "cm", scaling = 0.6)
ggsave(file.path(res_dir, "log2FC_by_celltype.svg"), width = 6, height = 8, units = "cm", scale = 1/0.6)
