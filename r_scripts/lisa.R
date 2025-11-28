################################################################################
# Script Name: lisa.R
# Description: Fait tourner une analyse des LISA avec les paramètres
#              fournis en entrée
# Author: Grégoire Le Campion, Antoine Le Doeuff
# Date created: 2025-11-28
################################################################################

# Suppression des warnings
options(warn = -1)

# Chargement des bibliothèques
suppressPackageStartupMessages({
  library(rgeoda)
  library(sf)
})

# //////////////////////////////////////////////////////////////////////////////
# Lire les arguments -----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

INPUT_PATH         <- args[1]
OUTPUT_PATH        <- args[2]
VARIABLE           <- args[3]
VARIABLE_2         <- args[4]
ANALYSIS_TYPE      <- args[5]
STANDARDIZE        <- as.logical(as.integer(args[6]))
CONTIGUITY_TYPE    <- args[7]
ORDER              <- as.integer(args[8])
STANDARDIZE_WEIGHT <- as.logical(as.integer(args[9]))
SIGNIFICANCE       <- as.numeric(args[10])

cat("=== Résumé des paramètres fournis ===\n")
cat("Chemin d'entrée              :", INPUT_PATH, "\n")
cat("Chemin de sortie             :", OUTPUT_PATH, "\n")
cat("Variable                     :", VARIABLE, "\n")
cat("Variable 2                   :", VARIABLE_2, "\n")
cat("Type d'analyse               :", ANALYSIS_TYPE, "\n")
cat("Standardisation              :", STANDARDIZE, "\n")
cat("Type de contiguïté           :", CONTIGUITY_TYPE, "\n")
cat("Ordre                        :", ORDER, "\n")
cat("Poids de standardisation     :", STANDARDIZE_WEIGHT, "\n")
cat("Alpha (p-value)              :", SIGNIFICANCE, "\n")
cat("=====================================\n")

# //////////////////////////////////////////////////////////////////////////////
# Lecture des données ----------------------------------------------------------
cat("=== DEBUT DE L'ANALYSE LISA ===\n\n")
# Lecture des données
cat("Lecture des données...\n")
sf_data <- st_read(INPUT_PATH, quiet=TRUE)
cat(sprintf("Nombre d'entités: %d\n", nrow(sf_data)))
cat(sprintf("Nombre de colonnes: %d\n", ncol(sf_data)))

# Vérification de la variable ..................................................
cat("\nVérification de la variable...\n")
if (!VARIABLE %in% names(sf_data)) {
  stop(sprintf("Variable %s non trouvée dans les données", VARIABLE))
}

if (ANALYSIS_TYPE == "bivariate") {
  if (!VARIABLE_2 %in% names(sf_data)) {
    stop(sprintf("Variable %s non trouvée dans les données", VARIABLE_2))
  }
}

cat("Variable(s) présente(s)\n")

# Standardisation ..............................................................
sf_data_work <- sf_data
if (STANDARDIZE) {
  sf_data_work[[VARIABLE]] <- as.numeric(scale(sf_data_work[[VARIABLE]]))

  if (ANALYSIS_TYPE == "bivariate") {
    sf_data_work[[VARIABLE]] <- as.numeric(scale(sf_data_work[[VARIABLE]]))
    sf_data_work[[VARIABLE_2]] <- as.numeric(scale(sf_data_work[[VARIABLE_2]]))
  }
}

# Vérification des NA ..........................................................
cat("\nVérification des valeurs manquantes...\n")
n_na_var1 <- sum(is.na(sf_data_work[[VARIABLE]]))
if (n_na_var1 > 0) {
  cat(sprintf("ATTENTION: %d valeurs manquantes pour %s\n", n_na_var1, VARIABLE))
}

if (ANALYSIS_TYPE == "bivariate") {
  n_na_var2 <- sum(is.na(sf_data_work[[VARIABLE_2]]))
  if (n_na_var2 > 0) {
    cat(sprintf("ATTENTION: %d valeurs manquantes pour %s\n", n_na_var2, VARIABLE_2))
  }
}

# Suppression des lignes avec NA
if (ANALYSIS_TYPE == "univariate") {
  complete_cases <- complete.cases(sf_data_work[[VARIABLE]])
} else {
  complete_cases <- complete.cases(
    sf_data_work[[VARIABLE]],
    sf_data_work[[VARIABLE_2]]
  )
}

if (sum(!complete_cases) > 0) {
  cat(sprintf("Suppression de %d lignes avec valeurs manquantes\n", sum(!complete_cases)))
  sf_data_work <- sf_data_work[complete_cases, ]
  sf_data <- sf_data[complete_cases, ]
}

# //////////////////////////////////////////////////////////////////////////////
# Prétraitement ----------------------------------------------------------------
# Création de la matrice de poids spatiaux .....................................
cat("\nCréation de la matrice de poids spatiaux...\n")
cat(sprintf("Type de contiguïté: %s\n", CONTIGUITY_TYPE))
cat(sprintf("Ordre: %s\n", ORDER))
cat("Standardisation des poids:", STANDARDIZE_WEIGHT,  "\n")

# Créer les poids avec rgeoda
if (CONTIGUITY_TYPE == "queen") {
  weights <- queen_weights(sf_data_work, order = ORDER,
                          include_lower_order = TRUE)
} else {
  weights <- rook_weights(sf_data_work, order = ORDER,
                         include_lower_order = TRUE)
}

cat(sprintf("Matrice de poids créée avec succès\n"))

# Préparation des données pour LISA - CRÉER UN DATA.FRAME EXPLICITE
cat("\nPréparation des données pour LISA...\n")

# Retirer la géométrie et créer un pur data.frame
data_no_geom <- st_drop_geometry(sf_data_work)

# Extraction des variables .....................................................
if (ANALYSIS_TYPE == "univariate") {
  var_data <- sf_data_work[[VARIABLE]]
  cat(sprintf("Variable %s extraite: %d valeurs\n", VARIABLE, length(var_data)))
} else {
  var_data1 <- sf_data_work[[VARIABLE]]
  var_data2 <- sf_data_work[[VARIABLE_2]]
  cat(sprintf(
    "Variables %s et %s extraites: %d valeurs chacune\n",
    VARIABLE, VARIABLE_2, length(var_data1)
  ))
}

# Créer le data.frame pour rgeoda
if (ANALYSIS_TYPE == "univariate") {
  var_df <- data.frame(var = data_no_geom[[VARIABLE]])
} else {
  var_df <- data.frame(
    var1 = data_no_geom[[VARIABLE]],
    var2 = data_no_geom[[VARIABLE_2]]
  )
}

# Calcul LISA
cat("\n=== CALCUL LISA ===\n")

if (ANALYSIS_TYPE == "univariate") {
  cat("Analyse univariée de la variable", VARIABLE, "\n")
  lisa_result <- local_moran(weights, var_df)
} else {
  cat("Analyse bivariée des variables", VARIABLE, "et", VARIABLE_2, "\n")
  lisa_result <- local_bimoran(weights, var_df)
}
cat("Calcul LISA terminé\n")

# Extraction des résultats
cat("\n=== EXTRACTION DES RESULTATS ===\n")
cat("Seuil de significativité utilisé:", SIGNIFICANCE, "\n")

# Extraire directement avec les fonctions rgeoda
lisa_i <- lisa_values(lisa_result)
pvals <- lisa_pvalues(lisa_result)
clusters <- lisa_clusters(lisa_result, cutoff = SIGNIFICANCE)
labels <- lisa_labels(lisa_result)

cat(sprintf("Statistiques extraites pour %d observations\n", length(lisa_i)))
cat(sprintf("Nombre de p-values: %d\n", length(pvals)))
cat(sprintf("Clusters identifiés\n"))

# Préparer les résultats - COPIE PROPRE de sf_data d'origine
result_sf <- sf_data

# Ajouter les résultats LISA
result_sf$LISA_I <- lisa_i
result_sf$LISA_pvalue <- pvals
result_sf$LISA_cluster <- clusters

# Créer une colonne catégorielle pour les clusters
cluster_labels <- c(
  "Not significant",
  "High-High",
  "Low-Low",
  "Low-High",
  "High-Low",
  "Undefined",
  "Isolated"
)
result_sf$LISA_category <- factor(clusters, levels = 0:6, labels = cluster_labels)

# Statistiques descriptives
cat("\n=== STATISTIQUES DESCRIPTIVES ===\n")
cat(sprintf("I de Moran moyen: %.4f\n", mean(lisa_i, na.rm=TRUE)))
cat(sprintf("I min: %.4f\n", min(lisa_i, na.rm=TRUE)))
cat(sprintf("I max: %.4f\n", max(lisa_i, na.rm=TRUE)))

cat("\n=== REPARTITION DES CLUSTERS ===\n")
cluster_table <- table(result_sf$LISA_category)
for (i in seq_along(cluster_table)) {
  cat(sprintf(
    "%s: %d (%.1f%%)\n",
    names(cluster_table)[i],
    cluster_table[i],
    100 * cluster_table[i] / nrow(result_sf)
  ))
}

# Statistiques par type de cluster
cat("\n=== STATISTIQUES PAR CLUSTER ===\n")
for (cat_name in names(cluster_table)) {
  if (cluster_table[cat_name] > 0) {
    subset_data <- result_sf[result_sf$LISA_category == cat_name, ]
    cat(sprintf("\n%s:\n", cat_name))
    cat(sprintf("  I moyen: %.4f\n", mean(subset_data$LISA_I, na.rm=TRUE)))
    cat(sprintf("  p-value moyenne: %.4f\n", mean(subset_data$LISA_pvalue, na.rm=TRUE)))
  }
}

# Sauvegarde
cat("\nSauvegarde des résultats...\n")
st_write(result_sf, OUTPUT_PATH, delete_dsn=TRUE, quiet=TRUE)

cat("\n=== ANALYSE LISA TERMINEE AVEC SUCCES ===\n")
