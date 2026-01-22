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
      legend.justification = c("right"),
      strip.text = element_text(size = rel(1.04)),
      axis.text = element_text(color = "black", size = rel(0.9)),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 0.95, size = rel(0.85))))


# Input data -------------------------------------------------------------------
## Datasets meta ---------------------------------------------------------------
ds_meta <- read_tsv("sc_dataset_meta.txt")

### Annotate Lung/Airway datasets with tissue type info ------------------------
ds_meta <- ds_meta %>% 
  mutate(
    Tissue = if_else(Tissue == "Nose", "UTR", Tissue), # UTR = Upper Respiratory Tract
    Tissue_type = case_when(
      # Sample %in% c("Lung", "Trachea") ~ "Solid\ntissue",
      Sample == "Lung" ~ "Solid tissue",
      Sample == "Trachea" ~ "Solid\ntissue",
      Sample == "BALF" ~ "BALF",
      Sample == "Endotracheal aspirates" ~ "ETA",
      Dataset == "GSE168215" ~ "Bronchial\nbrushing",
      Sample == "Nasal swabs" ~ "Nasal swabs",
      .default = NA_character_),
    Tissue_type = factor(Tissue_type, levels = c("Solid tissue", "Solid\ntissue", "BALF", "ETA", "Bronchial\nbrushing", "Nasal swabs")))


## GO ORA results --------------------------------------------------------------
### Redundant terms are substituted with representative ones
go_lipid_revigo_imp_res <- read_tsv(file.path(data_dir, "go_lipid_revigo_imp_res.txt"))

### Combine GO ORA results with datasets metadata
go_lipid_revigo_imp_2plot <- go_lipid_revigo_imp_res %>% 
  left_join(ds_meta, by = "Dataset") %>% 
  # Add labels for Up/Down-regulated GO terms
  mutate(
    symbol = if_else(status == "Up", "\u25D7", "\u25D6"),
    color = if_else(status == "Up", "#f4a261", "#0a9396")) 


## KEGG ORA results ------------------------------------------------------------
kegg_lipid <- read_tsv(file.path(data_dir, "kegg_lipid.txt"))

### Combine GO ORA results with datasets metadata
kegg_lipid_2plot <- kegg_lipid %>% 
  left_join(ds_meta, by = "Dataset") %>% 
  # Add labels for Up/Down-regulated GO terms
  mutate(
    symbol = if_else(status == "Up", "\u25D7", "\u25D6"),
    color = if_else(status == "Up", "#f4a261", "#0a9396")) 


# Semicircle plots -------------------------------------------------------------
## By tissue -------------------------------------------------------------------
tissues <- c("Lung", "Airway", "UTR", "Liver", "Blood", "Lymph node", "Brain", "Heart", "Spleen", "Kidney", "Urine", "Placenta")


### GO BP ----------------------------------------------------------------------
go_lipid_revigo_imp_2plot_by_tissue <- go_lipid_revigo_imp_2plot %>% 
  distinct(Description_revigo, Dataset, status, Tissue, symbol, color) %>% 
  # At least one cell type in the dataset has Up/Down-regulated GO term
  count(Description_revigo, status, Tissue, symbol, color) %>% 
  mutate(
    Tissue = factor(Tissue, levels = tissues),
    Description_revigo = fct_rev(as.factor(Description_revigo)))

min_size_go_t <- min(go_lipid_revigo_imp_2plot_by_tissue$n)
max_size_go_t <- max(go_lipid_revigo_imp_2plot_by_tissue$n)


go_lipid_revigo_imp_2plot_by_tissue %>%
  ggplot(aes(x = Tissue, y = Description_revigo)) +
  geom_text(aes(color = color, label = symbol, size = log2(n+1)*3), family = "Arial Unicode MS", key_glyph = "point") +
  scale_x_discrete(expand = expansion(add = .8)) +
  scale_y_discrete(expand = expansion(add = .8)) +
  scale_color_identity(
    "",
    breaks = c("#f4a261", "#0a9396"),
    labels = c("Upregulated", "Downregulated"),
    guide = "legend") +
  scale_size_identity(
    "Number of\ndatasets",
    breaks = log2(c(min_size_go_t, c(min_size_go_t:max_size_go_t)[c(FALSE, TRUE)])+1)*3,
    labels = c(min_size_go_t, c(min_size_go_t:max_size_go_t)[c(FALSE, TRUE)]),
    guide = "legend") +
  guides(color = guide_legend(override.aes = list(size = 5), order = 1)) +
  guides(size = guide_legend(override.aes = list(color = "#0a9396"), order = 2)) +
  labs(x = "", y = "") 
  
ggsave(file.path(res_dir, "GOBP_by_tissue_lipid_keywords.HC.png"), device = agg_png, dpi = 300, width = 12, height = 10, units = "cm", scaling = 0.6)
ggsave(file.path(res_dir, "GOBP_by_tissue_lipid_keywords.HC.svg"), width = 12, height = 10, units = "cm", scale = 1/0.6)


### KEGG -----------------------------------------------------------------------
kegg_lipid_2plot_by_tissue <- kegg_lipid_2plot %>% 
  distinct(Description, Dataset, status, Tissue, symbol, color) %>% 
  # At least one cell type in the dataset has Up/Down-regulated GO term
  count(Description, status, Tissue, symbol, color) %>% 
  mutate(
    Tissue = factor(Tissue, levels = tissues),
    Description = fct_rev(as.factor(Description)))

min_size_k_t <- min(kegg_lipid_2plot_by_tissue$n)
max_size_k_t <- max(kegg_lipid_2plot_by_tissue$n)


kegg_lipid_2plot_by_tissue %>%
  ggplot(aes(x = Tissue, y = Description)) +
  geom_text(aes(color = color, label = symbol, size = log2(n+1)*3), family = "Arial Unicode MS", key_glyph = "point") +
  scale_x_discrete(expand = expansion(add = .8)) +
  scale_y_discrete(expand = expansion(add = .8)) +
  scale_color_identity(
    "",
    breaks = c("#f4a261", "#0a9396"),
    labels = c("Upregulated", "Downregulated"),
    guide = "legend") +
  scale_size_identity(
    "Number of\ndatasets",
    breaks = log2(c(min_size_k_t, c(min_size_k_t:max_size_k_t)[c(FALSE, TRUE)])+1)*3,
    labels = c(min_size_k_t, c(min_size_k_t:max_size_k_t)[c(FALSE, TRUE)]),
    guide = "legend") +
  guides(color = guide_legend(override.aes = list(size = 5), order = 1)) +
  guides(size = guide_legend(override.aes = list(color = "#0a9396"), order = 2)) +
  labs(x = "", y = "") 

ggsave(file.path(res_dir, "KEGG_by_tissue_lipid_keywords.HC.png"), device = agg_png, dpi = 300, width = 12, height = 9, units = "cm", scaling = 0.6)
ggsave(file.path(res_dir, "KEGG_by_tissue_lipid_keywords.HC.svg"), width = 12, height = 9, units = "cm", scale = 1/0.6)


## By dataset - only Lung and UTR ----------------------------------------------
datasets_order <- c("Melms et al., 2021", "Liao et al., 2020", "Ahn et al., 2021")


### GO BP ----------------------------------------------------------------------
go_lipid_revigo_imp_2plot_by_dataset <- go_lipid_revigo_imp_2plot %>% 
  dplyr::filter(Tissue %in% c("Lung", "Airway", "UTR")) %>% 
  distinct(Description_revigo, Dataset, `Dataset name for plot`, celltype, status, Tissue, Tissue_type, symbol, color) %>% 
  count(Description_revigo, Dataset, `Dataset name for plot`, status, Tissue, Tissue_type, symbol, color) %>% 
  mutate(
    Tissue = factor(Tissue, levels = tissues),
    Description_revigo = fct_reorder(Description_revigo, n),
    `Dataset name for plot` = factor(`Dataset name for plot`, levels = datasets_order))

min_size_go_ds <- min(go_lipid_revigo_imp_2plot_by_dataset$n)
max_size_go_ds <- max(go_lipid_revigo_imp_2plot_by_dataset$n)


go_lipid_revigo_imp_2plot_by_dataset %>% 
  ggplot(aes(x = `Dataset name for plot`, y = Description_revigo)) +
  geom_text(aes(color = color, label = symbol, size = log2(n+1)*3), family = "Arial Unicode MS", key_glyph = "point") +
  facet_nested(~ Tissue + Tissue_type, scales = "free_x") +
  force_panelsizes(cols = c(1, 1, 1)) +
  scale_x_discrete(expand = expansion(add = .8)) +
  scale_y_discrete(expand = expansion(add = .8)) +
  scale_color_identity(
    "",
    breaks = c("#f4a261", "#0a9396"),
    labels = c("Upregulated", "Downregulated"),
    guide = "legend") +
  scale_size_identity(
    "Number of\ncell types",
    breaks = log2(c(min_size_go_ds, c(min_size_go_ds:max_size_go_ds)[c(FALSE, TRUE)])+1)*3,
    labels = c(min_size_go_ds, c(min_size_go_ds:max_size_go_ds)[c(FALSE, TRUE)]),
    guide = "legend") +
  guides(color = guide_legend(override.aes = list(size = 5), order = 1)) +
  guides(size = guide_legend(override.aes = list(color = "#0a9396"), order = 2)) +
  labs(x = "", y = "") 

ggsave(file.path(res_dir, "GOBP_by_dataset_lipid_keywords.HC.png"), device = agg_png, dpi = 300, width = 14, height = 10, units = "cm", scaling = 0.6)
ggsave(file.path(res_dir, "GOBP_by_dataset_lipid_keywords.HC.svg"), width = 14, height = 10, units = "cm", scale = 1/0.6)


### KEGG ----------------------------------------------------------------------
kegg_lipid_2plot_by_dataset <- kegg_lipid_2plot %>% 
  dplyr::filter(Tissue %in% c("Lung", "Airway", "UTR")) %>% 
  distinct(Description, Dataset, `Dataset name for plot`, celltype, status, Tissue, Tissue_type, symbol, color) %>% 
  count(Description, Dataset, `Dataset name for plot`, status, Tissue, Tissue_type, symbol, color) %>% 
  mutate(
    Tissue = factor(Tissue, levels = tissues),
    Description = fct_reorder(Description, n),
    `Dataset name` = factor(`Dataset name for plot`, levels = datasets_order))

min_size_k_ds <- min(kegg_lipid_2plot_by_dataset$n)
max_size_k_ds <- max(kegg_lipid_2plot_by_dataset$n)


kegg_lipid_2plot_by_dataset %>% 
  ggplot(aes(x = `Dataset name for plot`, y = Description)) +
  geom_text(aes(color = color, label = symbol, size = log2(n+1)*3), family = "Arial Unicode MS", key_glyph = "point") +
  facet_nested(~ Tissue + Tissue_type, scales = "free_x") +
  force_panelsizes(cols = c(4, 2, 2, 2, 2)) +
  scale_x_discrete(expand = expansion(add = .8)) +
  scale_y_discrete(expand = expansion(add = .8)) +
  scale_color_identity(
    "",
    breaks = c("#f4a261", "#0a9396"),
    labels = c("Upregulated", "Downregulated"),
    guide = "legend") +
  scale_size_identity(
    "Number of\ncell types",
    breaks = log2(c(min_size_k_ds, c(min_size_k_ds:max_size_k_ds)[c(FALSE, FALSE, FALSE, FALSE, TRUE)])+1)*3,
    labels = c(min_size_k_ds, c(min_size_k_ds:max_size_k_ds)[c(FALSE, FALSE, FALSE, FALSE, TRUE)]),
    guide = "legend") +
  guides(color = guide_legend(override.aes = list(size = 5), order = 1)) +
  guides(size = guide_legend(override.aes = list(color = "#0a9396"), order = 2)) +
  labs(x = "", y = "") 

ggsave(file.path(res_dir, "KEGG_by_dataset_lipid_keywords.HC.png"), device = agg_png, dpi = 300, width = 14, height = 5, units = "cm", scaling = 0.6)
ggsave(file.path(res_dir, "KEGG_by_dataset_lipid_keywords.HC.svg"), width = 14, height = 5, units = "cm", scale = 1/0.6)


## By cell type - only Lung and Airway -----------------------------------------
### GO BP ----------------------------------------------------------------------
go_lipid_revigo_imp_2plot_by_celltype <- go_lipid_revigo_imp_2plot %>% 
  dplyr::filter(Tissue %in% c("Lung", "Airway", "UTR")) %>%
  mutate(
    Tissue = factor(Tissue, levels = tissues),
    celltype_tidy = str_remove_all(celltype, "^[A-Z0-9]*_|_[A-Z0-9]*$"),
    celltype_tidy = str_replace(celltype_tidy, "_", " "),
    celltype_tidy = str_replace(celltype_tidy, "\\|", "/"), # T/NK cells
    celltype_tidy = str_remove(celltype_tidy, "s$"),
    celltype_tidy = str_remove(celltype_tidy, "^Sec "),
    celltype_tidy = ifelse(celltype_tidy == "Basal Suprabasal", "Suprabasal", celltype_tidy),
    celltype_tidy = ifelse(celltype_tidy == "Ciliated Ciliated diff", "Ciliated", celltype_tidy),
    celltype_tidy = ifelse(celltype_tidy == "Ciliated Multiciliated", "Multiciliated", celltype_tidy),
    celltype_tidy = ifelse(celltype_tidy == "Secretory diff", "Secretory", celltype_tidy),
    celltype_tidy = ifelse(celltype_tidy == "Epithelial Secretory", "Secretory", celltype_tidy),
    `Dataset name for plot` = factor(`Dataset name for plot`, levels = datasets_order)) %>% 
  distinct(Description_revigo, Dataset, `Dataset name for plot`, celltype_tidy, status, Tissue, Tissue_type, symbol, color)
  
go_lipid_revigo_imp_2plot_by_celltype %>% 
  ggplot(aes(x = celltype_tidy, y = Description_revigo)) +
  geom_text(aes(color = color, label = symbol), size = 5, family = "Arial Unicode MS", key_glyph = "point") +
  facet_nested(~ Tissue + Tissue_type + `Dataset name for plot`, scales = "free_x") +
  force_panelsizes(cols = c(6, 7, 7)) +
  scale_x_discrete(expand = expansion(add = .8)) +
  scale_y_discrete(expand = expansion(add = .8)) +
  scale_color_identity(
    "",
    breaks = c("#f4a261", "#0a9396"),
    labels = c("Upregulated", "Downregulated"),
    guide = "legend") +
  guides(color = guide_legend(override.aes = list(size = 5), order = 1)) +
  labs(x = "", y = "") +
  theme(legend.position = "top")

ggsave(file.path(res_dir, "GOBP_by_celltype_lipid_keywords.HC.png"), device = agg_png, dpi = 300, width = 14, height = 9, units = "cm", scaling = 0.6)
ggsave(file.path(res_dir, "GOBP_by_celltype_lipid_keywords.HC.svg"), width = 14, height = 9, units = "cm", scale = 1/0.6)


### KEGG ----------------------------------------------------------------------
kegg_lipid_2plot_2plot_by_celltype <- kegg_lipid_2plot %>% 
  dplyr::filter(Tissue %in% c("Lung", "Airway", "UTR")) %>%
  mutate(
    Tissue = factor(Tissue, levels = tissues),
    celltype_tidy = str_remove_all(celltype, "^[A-Z0-9]*_|_[A-Z0-9]*$"),
    celltype_tidy = str_replace(celltype_tidy, "_", " "),
    celltype_tidy = str_replace(celltype_tidy, "\\|", "/"), # T/NK cells
    celltype_tidy = str_remove(celltype_tidy, "s$"),
    celltype_tidy = str_remove(celltype_tidy, "^Sec "),
    celltype_tidy = ifelse(celltype_tidy == "Basal Suprabasal", "Suprabasal", celltype_tidy),
    celltype_tidy = ifelse(celltype_tidy == "Ciliated Ciliated diff", "Ciliated", celltype_tidy),
    celltype_tidy = ifelse(celltype_tidy == "Ciliated Multiciliated", "Multiciliated", celltype_tidy),
    celltype_tidy = ifelse(celltype_tidy == "Secretory diff", "Secretory", celltype_tidy),
    celltype_tidy = ifelse(celltype_tidy == "Epithelial Secretory", "Secretory", celltype_tidy),
    `Dataset name` = factor(`Dataset name for plot`, levels = datasets_order)) %>% 
  distinct(Description, Dataset, `Dataset name for plot`, celltype_tidy, status, Tissue, Tissue_type, symbol, color)


kegg_lipid_2plot_2plot_by_celltype %>% 
  ggplot(aes(x = celltype_tidy, y = Description)) +
  geom_text(aes(color = color, label = symbol), size = 5, family = "Arial Unicode MS", key_glyph = "point") +
  facet_nested(~ Tissue + Tissue_type + `Dataset name for plot`, scales = "free_x") +
  force_panelsizes(cols = c(8, 7, 3)) +
  scale_x_discrete(expand = expansion(add = .8)) +
  scale_y_discrete(expand = expansion(add = .8)) +
  scale_color_identity(
    "",
    breaks = c("#f4a261", "#0a9396"),
    labels = c("Upregulated", "Downregulated"),
    guide = "legend") +
  guides(color = guide_legend(override.aes = list(size = 5), order = 1)) +
  labs(x = "", y = "") +
  theme(legend.position = "top")

ggsave(file.path(res_dir, "KEGG_by_celltype_lipid_keywords.HC.png"), device = agg_png, dpi = 300, width = 14, height = 6, units = "cm", scaling = 0.6)
ggsave(file.path(res_dir, "KEGG_by_celltype_lipid_keywords.HC.svg"), width = 16, height = 6, units = "cm", scale = 1/0.6)

