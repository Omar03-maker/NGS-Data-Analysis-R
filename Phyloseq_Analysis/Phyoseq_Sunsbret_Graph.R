# Installation des packages et libraires : 
install.packages("plotly")     
install.packages("readxl")    
install.packages("ggplot2") 
install.packages("plotly")
install.packages("sunburstR")
library("phyloseq")
library("ggplot2")      
library("readxl")       
library("dplyr")        
library("tibble")
library("dplyr") 
library("tidyr")
library("RColorBrewer")
library("sunburstR")


# Charement des fichiers : 
OTU_MAT <- read_excel("Chemin vers /Table_OTU.xlsx", sheet = "OTU_Counts")
TAX_MAT <- read_excel("Chemin vers /Table_taxonomy.xlsx", sheet = "taxonenew")
SAMPLES_DF <- read_excel("Chemin vers /Table_sample.xlsx", sheet = "sample")

head(OTU_MAT)
head(TAX_MAT)
head(SAMPLES_DF)

OTU_MAT<- OTU_MAT %>%
  tibble::column_to_rownames("OTU") 
TAX_MAT <- TAX_MAT%>% 
  tibble::column_to_rownames("OTU")
SAMPLES_DF<- SAMPLES_DF%>% 
  tibble::column_to_rownames("echantillon") 

OTU_MAT <- as.matrix(OTU_MAT)
TAX_MAT <- as.matrix(TAX_MAT)

# Création de l'objet phyloseq 
Otu = otu_table(OTU_MAT, taxa_are_rows = TRUE)
Tax = tax_table(TAX_MAT)
SAMPLES = sample_data(SAMPLES_DF)
Phylo <- phyloseq(Otu, Tax, SAMPLES)
Phylo

# Visualisation des données : 
sample_names(Phylo)
rank_names(Phylo)
sample_variables(Phylo)

# Fonction pour créer un diagramme sunburst à partir d'un objet phyloseq
create_sunburst_fix <- function(Phylo, sample_name = NULL, top_n = 50) {
  
  # Extraire les données directement
  tax_data <- as.data.frame(tax_table(Phylo))
  otu_data <- as.data.frame(otu_table(Phylo))
  
  # Vérifier l'orientation et corriger si nécessaire
  if (taxa_are_rows(Phylo)) {
    # Taxa en lignes - normal
    if (!is.null(sample_name) && sample_name %in% colnames(otu_data)) {
      abundance <- otu_data[, sample_name]
    } else {
      abundance <- rowSums(otu_data)}
  } else {
    # Taxa en colonnes - transposer
    otu_data <- t(otu_data)
    if (!is.null(sample_name) && sample_name %in% colnames(otu_data)) {
      abundance <- otu_data[, sample_name]
    } else {
      abundance <- rowSums(otu_data)}}
  
  # Combiner taxonomie et abondance
  combined_data <- data.frame(tax_data, abundance = abundance)
  
  # Trier par abondance décroissante et prendre le top_n
  combined_data <- combined_data[order(combined_data$abundance, decreasing = TRUE), ]
  combined_data <- head(combined_data, top_n)
  
  # Remplacer les NA par "Unclassified"
  combined_data[is.na(combined_data)] <- "Unclassified"
  
  # Enlever les lignes avec abondance nulle
  combined_data <- combined_data[combined_data$abundance > 0, ]
  
  # Créer le chemin taxonomique pour sunburst
  # Utiliser les noms de colonnes réels de votre table taxonomique
  tax_cols <- colnames(tax_data)
  
  if (length(tax_cols) >= 6) {
    # Si vous avez au moins 6 niveaux taxonomiques
    combined_data$path <- paste(
      combined_data[, tax_cols[1]], # Kingdom/Domain
      combined_data[, tax_cols[2]], # Phylum
      combined_data[, tax_cols[3]], # Class
      combined_data[, tax_cols[4]], # Order
      combined_data[, tax_cols[5]], # Family
      combined_data[, tax_cols[6]], # Genus
      sep = "-")
  } else {
    # Utiliser tous les niveaux disponibles
    path_parts <- do.call(paste, c(combined_data[, tax_cols], sep = "-"))
    combined_data$path <- path_parts}
  
  # Préparer les données pour sunburst
  sunburst_data <- combined_data %>%
    select(path, abundance) %>%
    filter(abundance > 0)
  
  # Création du sunburst graphique
  sunburst(sunburst_data, width = "100%", height = 500)}

# Affichage de la figure 
fig_P01 <- create_sunburst_fix(CARBOM, sample_name = "P01", top_n = 50)
fig_P01
