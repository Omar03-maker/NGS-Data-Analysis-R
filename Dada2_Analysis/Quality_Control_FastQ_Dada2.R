# INSTALLATION DES PACKAGES 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")
BiocManager::install("phyloseq")
install.packages("ggplot2")
install.packages("fastqcr")
library(phyloseq)
library(dada2)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(fastqcr)

# DÉFINIR LE CHEMIN DES DONNÉES
path <- "Chemin vers /miseqsopdata"
list.files(path)

# RÉCUPÉRER LES FICHIERS R1 ET R2 SÉPARÉMENT
# Les fichiers dans le dossiers miseqsopdata sont des exemples R1_001.fastq (Reads F) et R2_001.fastq (Reads R) : 
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Fonction pour créer le graphique de controle de la qualité 
add_quality_zones <- function(p) {
  plot_data <- ggplot_build(p)$data[[1]]
  x_min <- min(plot_data$x)
  x_max <- max(plot_data$x)
  
  p_new <- p + 
    annotate("rect", xmin = x_min, xmax = x_max, ymin = 0, ymax = 20,
             fill = "pink", alpha = 0.4) +
    annotate("rect", xmin = x_min, xmax = x_max, ymin = 20, ymax = 28,
             fill = "lightgreen", alpha = 0.4) +
    annotate("rect", xmin = x_min, xmax = x_max, ymin = 28, ymax = 42,
             fill = "grey", alpha = 0.4) +
    geom_hline(yintercept = 20, linetype = "dashed", color = "red", alpha = 0.7) +
    geom_hline(yintercept = 28, linetype = "dashed", color = "orange", alpha = 0.7) +
    ylim(0, 42) +
   
    geom_point(aes(x = -Inf, y = -Inf, color = "Bonne qualité"), alpha = 0) +
    geom_point(aes(x = -Inf, y = -Inf, color = "Qualité raisonnable"), alpha = 0) +
    geom_point(aes(x = -Inf, y = -Inf, color = "Mauvaise qualité"), alpha = 0) +
    scale_color_manual(name = "Qualité",
                       values = c("Bonne qualité" = "grey", 
                                  "Qualité raisonnable" = "lightgreen",
                                  "Mauvaise qualité" = "pink"),
      guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) +
    theme(legend.position = "right",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 13, face = "bold"),
          legend.key = element_blank(),
          legend.justification = "top")
  
  return(p_new)
}

# EXTRAIRE LES NOMS D'ÉCHANTILLONS ET VERIFIER LA QUALITÉ 
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

p <- plotQualityProfile(fnRs[1:2])
p_final <- add_quality_zones (p)
print(p_final)

# DÉFINIR LES CHEMINS DE SORTIE
## - Reads forwards : F_filt.fastq.gz (R1_001.fastq filtered) 
## - Reads Reverse : R_filt.fastq.gz (R1_001.fastq filtered)
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# FILTRER ET TRIMMER
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen=c(240, 160),  # 240pb pour Forward, 160pb pour Reverse
                     maxN=0, maxEE=c(2,2), truncQ=2, 
                     rm.phix=TRUE, compress=TRUE, multithread=FALSE)
head(out)

# VISUALISER LA QUALITÉ DES FICHIERS FILTRÉS
plotQualityProfile(filtFs[1:2])  
plotQualityProfile(filtRs[1:2])  

pf <- plotQualityProfile(filtRs[1:2])
pf_final <- add_quality_zones (pf)
print(pf_final)

# Afficher le résumé du filtrage
cat("Résumé du filtrage:")
head(out)
cat("Total des reads avant filtrage:", sum(out[,1]), "")
cat("Total des reads après filtrage:", sum(out[,2]), "")
cat("Pourcentage de reads conservés:", round(sum(out[,2])/sum(out[,1])*100,2),"%")

# APPRENDRE LES ERREURS
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# DÉRÉPLICATION
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# AFFICHER STATISTIQUES DE DÉRÉPLICATION
cat("Les échantillons forward suivants ont été dérepliqués:/n", paste(names(derepFs), collapse=", "), "/n")
cat("Les échantillons reverse suivants ont été dérepliqués:/n", paste(names(derepRs), collapse=", "), "/n")
cat("Statistiques de déréplication pour le premier échantillon forward:/n"); print(derepFs[[1]])
cat("Statistiques de déréplication pour le premier échantillon reverse:/n"); print(derepRs[[1]])

# NOMBRE DE SÉQUENCES UNIQUES
forward_unique_seqs <- sapply(derepFs, function(x) length(x$uniques))
reverse_unique_seqs <- sapply(derepRs, function(x) length(x$uniques))
cat("Nombre de séquences uniques pour les échantillons forward:/n"); print(forward_unique_seqs)
cat("Nombre de séquences uniques pour les échantillons reverse:/n"); print(reverse_unique_seqs)

# INFÉRENCE DES ÉCHANTILLONS 
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]
dadaRs[[1]]

# FUSION DES PAIRES
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])

# CONSTRUCTION DE LA TABLE DE SÉQUENCES
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

# SUPPRESSION DES CHIMÈRES
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
head(seqtab.nochim)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) # Pourcentage séquences non-chimériques

# SUIVI DU PIPELINE
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN),
               sapply(dadaRs, getN), 
               sapply(mergers, getN), 
               rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR",
                     "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# ATTRIBUTION DE TAXONOMIE
taxa <- assignTaxonomy(seqtab.nochim,
                       "C:/Users/elhad/Documents/El_Hadji_Omar/DIC_G2B/DIC2/ANALYSE_NGS/miseqsopdata/MiSeq_SOP/silva_nr_v128_train_set.fa.gz", 
                       multithread=TRUE)
taxa <- addSpecies(taxa,
                   "C:/Users/elhad/Documents/El_Hadji_Omar/DIC_G2B/DIC2/ANALYSE_NGS/miseqsopdata/MiSeq_SOP/silva_species_assignment_v128.fa.gz")

# ÉVALUATION DE LA PRÉCISION
unqs.mock <- seqtab.nochim["Mock",]
if(length(unqs.mock) > 0) {
  unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE)
  cat("DADA2 a inféré", length(unqs.mock), "séquences d'échantillons présentes dans la communauté Mock./n")
  mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
  match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
  cat("Parmi celles-ci,", sum(match.ref), "étaient des correspondances exactes avec les séquences de référence attendues./n")
}

# CRÉATION DE L'OBJET PHYLOSEQ
theme_set(theme_bw())
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject, 1, 1)
subject <- substr(subject, 2, 999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day, row.names=samples.out)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"

# CORRECTION: S'assurer que les noms de séquences sont cohérents
taxa.clean <- taxa
rownames(taxa.clean) <- colnames(seqtab.nochim)

# CRÉATION DE L'OBJET PHYLOSEQ
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa.clean))

# FILTRAGE DE L'ÉCHANTILLON MOCK
ps <- prune_samples(sample_names(ps) != "Mock", ps)

# VISUALISATION DE L'ALPHA-DIVERSITÉ
plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When") +
  geom_point(size=3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Jour", y = "Mesure de diversité alpha")

# TRANSFORMATION DES DONNÉES POUR L'ANALYSE BRAY-CURTIS
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS") +
  geom_point(size=3) +
  theme(legend.position="right")

# GRAPHIQUE À BARRES 
# Sélection des 20 taxa les plus abondants
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)

# Génération du graphique
plot_bar(ps.top20, x="Day", fill="Family") + 
  facet_wrap(~When, scales="free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  labs(x = "Jour", y = "Abondance relative", fill = "Famille")

# VISUALISATION AU NIVEAU DES PHYLUM
tax_table(ps)
phylum_counts <- tax_glom(ps, "Phylum")
phylum_counts_rel <- transform_sample_counts(phylum_counts, function(OTU) OTU/sum(OTU))
plot_bar(phylum_counts_rel, x="Day", fill="Phylum") + 
  facet_wrap(~When, scales="free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Jour", y = "Abondance relative", fill = "Phylum")
