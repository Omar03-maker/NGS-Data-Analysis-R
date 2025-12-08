t# Package prérequis
if (!require("phyloseq")) install.packages("phyloseq")
if (!require("readxl")) install.packages("readxl")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tibble")) install.packages("tibble")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
if (!require("reshape2")) install.packages("reshape2")
if (!require("gridExtra")) install.packages("gridExtra")
if (!require("grid")) install.packages("grid")
if (!require("gtable")) install.packages("gtable")
if (!require("ggdendro")) install.packages("ggdendro")
if (!require("pheatmap")) install.packages("pheatmap")
# Chargement des packages 
library("phyloseq")
library("readxl")
library("ggplot2")
library("dplyr")
library("tibble")
library("RColorBrewer")
library("reshape2")
library("gridExtra")
library("grid")
library("gtable")
library("ggdendro")
library("pheatmap")

# CHARGEMENT DES DONNÉES (adaptez les chemins d'accès)
OTU_MAT <- read_excel("Chemin vers /Table_OTU.xlsx", sheet = "OTU_Counts")
TAX_MAT <- read_excel("Chemin vers /Table_taxonomy.xlsx")
SAMPLES_DF <- read_excel("Chemin vers /Table_sample.xlsx", sheet = "sample")

# PRÉPARATION DES DONNÉES POUR PHYLOSEQ
OTU_MAT <- OTU_MAT %>% tibble::column_to_rownames("OTU")
TAX_MAT <- TAX_MAT %>% tibble::column_to_rownames("OTU")
SAMPLES_DF <- SAMPLES_DF %>% tibble::column_to_rownames("echantillon")

OTU_MAT <- as.matrix(OTU_MAT)
TAX_MAT <- as.matrix(TAX_MAT)

# CRÉATION DE L'OBJET PHYLOSEQ
Otu <- otu_table(OTU_MAT, taxa_are_rows = TRUE)
Tax <- tax_table(TAX_MAT)
SAMPLES <- sample_data(SAMPLES_DF)
Phylo <- phyloseq(Otu, Tax, SAMPLES)

# FONCTION VOIR LES NIVEAUX TAXONOMIQUES PARTAGÉES
calculate_shared_species <- function(phylo_obj, tax_level = "family", min_abundance = 1) {
  
  # Filtrer les taxa avec abondance minimale
  phylo_filtered <- filter_taxa(phylo_obj, function(x) sum(x >= min_abundance) > 0, TRUE)
  
  # Extraire les données OTU et taxonomie
  otu_data <- as.data.frame(otu_table(phylo_filtered))
  tax_data <- as.data.frame(tax_table(phylo_filtered))
  
  # Vérifier l'orientation des données
  if (!taxa_are_rows(phylo_filtered)) { otu_data <- t(otu_data)}
  
  # Agréger par niveau taxonomique choisi
  if (tax_level %in% colnames(tax_data)) {
    # Combiner OTU et taxonomie
    combined_data <- cbind(tax_data[, tax_level, drop = FALSE], otu_data)
    
    # Agréger par niveau taxonomique
    aggregated <- combined_data %>%
      group_by(!!sym(tax_level)) %>%
      summarise(across(where(is.numeric), sum), .groups = 'drop') %>%
      column_to_rownames(tax_level)
  } else {
    aggregated <- otu_data}
  
  # Convertir en présence/absence (binaire)
  binary_data <- (aggregated > 0) * 1
  
  # Calculer la matrice de similarité (espèces partagées)
  samples <- colnames(binary_data)
  n_samples <- length(samples)
  
  shared_matrix <- matrix(0, nrow = n_samples, ncol = n_samples,
                          dimnames = list(samples, samples))
  
  for (i in 1:n_samples) {for (j in 1:n_samples) {
      # Nombre d'espèces partagées entre échantillons i et j
      shared_matrix[i, j] <- sum(binary_data[, i] == 1 & binary_data[, j] == 1)}}
  
  return(shared_matrix)}

# CRÉATION DE LA HEATMAP 
create_shared_species_heatmap <- function(phylo_obj, tax_level = "species", 
                                          min_abundance = 1, title = NULL) {
  
  # Calculer la matrice des espèces partagées
  shared_matrix <- calculate_shared_species(phylo_obj, tax_level, min_abundance)
  
  # Convertir en format long pour ggplot2
  melted_data <- melt(shared_matrix)
  colnames(melted_data) <- c("Sample1", "Sample2", "Shared_Species")
  
  # Titre par défaut
  if (is.null(title)) {title <- paste(" GROUPE TAXONOMIQUE PARTAGES ", 
                                      stringr::str_to_title(tax_level))}
  
  # Parametre de la heatmap
  p <- ggplot(melted_data, aes(x = Sample1, y = Sample2, fill = Shared_Species)) +
    geom_tile(color = "white", size = 0.1) +
    scale_fill_gradient2(low = "lightblue", mid = "steelblue", high = "darkblue",
                         midpoint = median(melted_data$Shared_Species),
                         name = " Nombre partagées") +
    labs(title = title,
         x = "Échantillons", 
         y = "Échantillons") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right",
      panel.grid = element_blank()
    ) +
    coord_fixed()
  
  return(p)}

# GÉNÉRATION DES GRAPHIQUES POUR DIFFÉRENTS NIVEAUX TAXONOMIQUES

# Heatmap basée sur les familles
heatmap_family <- create_shared_species_heatmap(Phylo, tax_level = "family", 
                                                title = "Familles partagées entre échantillons")
print(heatmap_family)

# Heatmap basée sur les ordres
heatmap_order <- create_shared_species_heatmap(Phylo, tax_level = "other", 
                                               title = "Ordres partagés entre échantillons")
print(heatmap_order)

# Heatmap basée sur les phylum
heatmap_class<- create_shared_species_heatmap(Phylo, tax_level = "class", 
                                                title = "Class partagés entre échantillons")
print(heatmap_class)

# Heatmap basée sur les phylum
heatmap_genius <- create_shared_species_heatmap(Phylo, tax_level = "genius", 
                                                title = "genius partagés entre échantillons")
print(heatmap_genius)

# Heatmap basée sur les phylum
heatmap_phylum <- create_shared_species_heatmap(Phylo, tax_level = "phylum", 
                                                title = "phylum partagés entre échantillons")
print(heatmap_phylum)


# VERSION AVEC CLUSTERING HIÉRARCHIQUE ET DENDROGRAMMES 
if (!require("pheatmap")) install.packages("pheatmap")
library("pheatmap")

create_phylogenetic_heatmap <- function(phylo_obj, tax_level = "family", min_abundance = 1) {
  
  shared_matrix <- calculate_shared_species(phylo_obj, tax_level, min_abundance)
  
  # Créer la heatmap avec dendrogrammes automatiques
  pheatmap(shared_matrix,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean", 
           clustering_method = "average",
           color = colorRampPalette(c("lightblue", "steelblue", "darkblue"))(100),
           show_rownames = TRUE,
           show_colnames = TRUE,
           fontsize = 10,
           fontsize_row = 8,
           fontsize_col = 8,
           main = paste("Espèces partagées - Clustering phylogénétique\n", 
                        stringr::str_to_title(tax_level)),
           legend_title = "Espèces\npartagées",
           border_color = "white",
           cellwidth = 20,
           cellheight = 20)
}

# Heatmap avec clustering phylogénétique 
create_phylogenetic_heatmap(Phylo, tax_level = "family")
create_phylogenetic_heatmap(Phylo, tax_level = "phylum")
create_phylogenetic_heatmap(Phylo, tax_level = "other")

# AFFICHER LES STATISTIQUES RÉSUMÉES
print("=== RÉSUMÉ DES ESPÈCES PARTAGÉES ===")
shared_stats <- calculate_shared_species(Phylo, "family")
cat("Nombre moyen partagés:", round(mean(shared_stats[upper.tri(shared_stats)]), 2), "family")
cat("Nombre maximum partagés:", max(shared_stats), "family")
cat("Nombre minimum partagés:", min(shared_stats[upper.tri(shared_stats)]), "family")
