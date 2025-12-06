# 🧬 Pipeline d'analyse et visualisation de données NGS avec R 
## Description
Il s'agit d'une pipeline complet pour l'analyse et la visualisation de données de séquençage nouvelle génération (NGS) avec le logiciel. Elle utilise deux jeux de données d'exemple pour vous permettre de tester le workflow :
- Données 16S rRNA d'échantillons microbiome intestinal de souris (Mothur MiSeq SOP)
- Données 18S rRNA de plancton océanique de la croisière CARBOM (Brésil, 2018)
 
Le workflow complet transforme les fichiers FASTQ bruts en visualisations et analyses statistiques de la diversité microbienne. Elle utilise notamment les packages Dada2 et phyloseq de R : 
- DADA2 pour le traitement des séquences brutes :
- Phyloseq pour l'analyse des données statistique et la visualisation :

## 📥 Prérequis et Téléchargement des données d'exemples :
Pour utiliser la pipeline il est essentiel d'installer R et R studio ainsi que les packages nécessaires Phyloqeq et Dada2.   
## 1 - Données pour Dada2 : 
### Pour l'analyse et le controle qualité télécharger les fichiers FASTQ du 16S rRNA inclus dans le fichier MiSeqSOPData : 
- wget http://www.mothur.org/w/images/d/d6/MiSeqSOPData.zip
- unzip MiSeqSOPData.zip -d data/MiSeq_SOP/
### Pour l'assignation taxonomique télécharger les fichiers via le lien suivant : https://zenodo.org/records/4587955
- assignTaxonomy (silva_nr_v128_train_set.fa.gz)
- addspecies (silva_species_assignment_v128.fa.gz) 
## 2 - Données pour Phyloseq : 
Le package phyloseq permet l'analyse et la visualisation de la diversité microbienne en utilisant 3 fichiers : 
- Table_OTU : Contient les échantillons codé en format OTU
- Taxonomy_Table : Classification taxonomique des OTUs
- Table_Sample : Contient les métadonnées des échantillons
## Utilsation de vos propres fichiers : Vous pouvez utiliser les fichiers correspondant a votre analyse NGS et remplacer dans les dossiers correspondant     

## 🚀 Utilisation de la pipeline
### Partie 1 : Traitement de donnnées avec DADA2
Étapes principales :
- Inspection de la qualité des séquences brutes
- Filtrage et trimming basé sur les scores de qualité
- Apprentissage du modèle d'erreur et débruitage 
- Détection et suppression des chimères
- Assignation taxonomique 
### Partie 2 : Analyses avec Phyloseq
#### 📊 Visualisations de la diverstité et composition des échantillons
- Bar plots et heatmpas : Composition taxonomique et abondance dans les échantillons (Personnalisables par niveau taxonomique (Phylum, Classe, Genre, etc.)
#### 📈 Analyses de diversité
- Alpha et beta diversité : Richesse en espèces de chaque échantillon et Différences entre échantillons
- Indices de diverstié : Shannon, Simpson, Chao1
 
### Les scripts vous guide à travers chaque type d'analyse
- Importation des données
- Filtrage et normalisation des données
- Visualisation et analyses statistique

## 🎯 En résumé
Cette pipeline vous permet de passer de fichiers FASTQ bruts à des analyses biologiques complètes en quelques étapes :
- DADA2 nettoie vos données et identifie les variants biologiques réels
- Phyloseq transforme ces variants en insights écologiques visuels
Que vous soyez débutant ou expert en bioinformatique, les scripts sont commentés et structurés pour faciliter la compréhension et l'adaptation à vos propres projets ! 🚀

## Auteur
- El Hadji Omar Dia
- GitHub: @Omar03-maker

## 📄 Citation : 
Si vous utilisez cette pipeline dans vos travaux, veuiller citer :
- Pour DADA2 : Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP (2016). DADA2: High-resolution sample inference from Illumina amplicon data. Nature Methods 13:581-583.
- Pour Phyloseq : McMurdie PJ, Holmes S (2013). phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. PLOS ONE 8(4):e61217. 

# ⭐ Si cela vous a été utile n'oubliez pas donner une étoile au repo ! 
