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

go_lipid_revigo_imp_2genetable_lfc <- go_lipid_revigo_imp_2genetable_lfc %>% 
  dplyr::filter(Tissue %in% c("Lung", "Airway", "Nose")) %>%
  mutate(
    celltype_tidy = str_remove_all(celltype, "^[A-Z0-9]*_|_[A-Z0-9]*$"),
    celltype_tidy = str_replace(celltype_tidy, "_", " "),
    celltype_tidy = str_replace(celltype_tidy, "\\|", "/"), # T/NK cells
    celltype_tidy = str_remove(celltype_tidy, "s$"),
    celltype_tidy = str_remove(celltype_tidy, "^Sec "),
    celltype_tidy = ifelse(celltype_tidy == "Basal Suprabasal", "Suprabasal", celltype_tidy),
    celltype_tidy = ifelse(celltype_tidy == "Ciliated Ciliated diff", "Ciliated", celltype_tidy),
    celltype_tidy = ifelse(celltype_tidy == "Ciliated Multiciliated", "Multiciliated", celltype_tidy),
    celltype_tidy = ifelse(celltype_tidy == "Secretory diff", "Secretory", celltype_tidy),
    celltype_tidy = ifelse(celltype_tidy == "Epithelial Secretory", "Secretory", celltype_tidy))

## Calculate distinct GO terms by cell type
go_lipid_revigo_imp_2genetable_lfc %>% 
  distinct(ID, Description, status, celltype_tidy) %>% 
  count(celltype_tidy, sort = TRUE)

## Calculate average log2FC per dataset (Lung, Airway, and UTR samples)
go_lipid_revigo_imp_2genetable_lfc_avg <- go_lipid_revigo_imp_2genetable_lfc %>% 
  distinct(celltype_tidy, avg_log2FC, status, Dataset) %>% 
  group_by(celltype_tidy, status, Dataset) %>% 
  summarise(
    n = n(),
    avg_avg_log2FC = mean(avg_log2FC)) %>% 
  ungroup() 

go_lipid_revigo_imp_2genetable_lfc_avg %>% 
  group_by(celltype_tidy, status) %>% 
  summarise(avg_avg_log2FC = mean(avg_avg_log2FC)) %>% 
  ungroup() %>% 
  as.data.frame()


# Plot -------------------------------------------------------------------------
go_lipid_revigo_imp_2genetable_lfc_avg %>% 
  mutate(
    abs_avg_avg_log2FC = abs(avg_avg_log2FC),
    celltype_tidy = fct_reorder(celltype_tidy, abs_avg_avg_log2FC, .fun = mean, .desc = FALSE),
    status = fct_rev(status)) %>% 
  ggplot(aes(x = avg_avg_log2FC, y = celltype_tidy, color = status)) +
  stat_summary(
    fun.y = mean,
    fun.ymin = function(x) mean(x) - sd(x), 
    fun.ymax = function(x) mean(x) + sd(x), 
    geom = "pointrange", alpha = 0.6) +
  geom_point() +
  facet_wrap2(
    vars(status), 
    scales = "free_x", 
    labeller = as_labeller(c("Up" = "Upregulated", "Down" = "Downregulated"))) +
  facetted_pos_scales(x = list(
    scale_x_continuous(limits = c(0.6, 4.8), breaks = c(1, 2.5, 4)),
    scale_x_continuous(limits = c(-3, -1), breaks = seq(-2.5, -1.5, 1)))) +
  scale_color_manual(values = c("#f4a261", "#0a9396")) + 
  labs(x = "Average"~log[2]~"fold change", y = NULL) 

ggsave(file.path(res_dir, "log2FC_by_celltype.png"), device = agg_png, dpi = 300, width = 6, height = 6, units = "cm", scaling = 0.6)
ggsave(file.path(res_dir, "log2FC_by_celltype.svg"), width = 6, height = 6, units = "cm", scale = 1/0.6)
