################################################################################
# Script Name: gwr.R
# Description: Fait tourner une GWR avec les paramètres fournis en entrée
# Author: Grégoire Le Campion, Antoine Le Doeuff
# Date created: 2025-11-26
################################################################################

# Suppression des warnings
options(warn = -1)

# Chargement des bibliothèques
suppressPackageStartupMessages({
  library(GWmodel)
  library(sp)
  library(sf)
})

# //////////////////////////////////////////////////////////////////////////////
# Read Args --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

INPUT_PATH       <- args[1]
OUTPUT_PATH      <- args[2]
DEPENDENT_VAR    <- args[3]
INDEPENDENT_VARS <- unlist(strsplit(args[4], split = ",", fixed = TRUE))
KERNEL           <- args[5]
APPROACH         <- args[6]
ADAPTATIVE       <- as.logical(as.integer(args[7]))
BANDWIDTH        <- as.numeric(args[8])
NEIGHBORS        <- as.integer(args[9])
STANDARDIZE      <- as.logical(as.integer(args[10]))
ROBUST           <- as.logical(as.integer(args[11]))

cat("=== Résumé des paramètres fournis ===\n")
cat("Chemin d'entrée (INPUT_PATH)        :", INPUT_PATH, "\n")
cat("Chemin de sortie (OUTPUT_PATH)      :", OUTPUT_PATH, "\n")
cat("Variable dépendante                 :", DEPENDENT_VAR, "\n")
cat("Variables indépendantes             :", paste(INDEPENDENT_VARS, collapse = ", "), "\n")
cat("Noyau choisi (KERNEL)               :", KERNEL, "\n")
cat("Approche (APPROACH)                 :", APPROACH, "\n")
cat("Bande passante adaptative           :", ADAPTATIVE, "\n")
cat("Bande passante (BANDWIDTH)          :", BANDWIDTH, "\n")
cat("Nombre de voisins (NEIGHBORS)       :", NEIGHBORS, "\n")
cat("Standardisation (STANDARDIZE)       :", STANDARDIZE, "\n")
cat("Méthode robuste (ROBUST)            :", ROBUST, "\n")
cat("=====================================\n")

# NOTE: La validation des paramètres est effectuée dans la méthode
# GWRAnalysisModule.write_args()

# //////////////////////////////////////////////////////////////////////////////
# Lecture des données ----------------------------------------------------------
cat("=== DEBUT DE L'ANALYSE GWR ===\n\n")

cat("Lecture des données...\n")
sf_data <- st_read(INPUT_PATH, quiet=TRUE)
cat(sprintf("Nombre d'entités: %d\n", nrow(sf_data)))
cat(sprintf("Nombre de colonnes: %d\n", ncol(sf_data)))

# //////////////////////////////////////////////////////////////////////////////
# Prétraitement ----------------------------------------------------------------
# Vérification des variables ...................................................
cat("\nVérification des variables...\n")
all_vars <- c(DEPENDENT_VAR, INDEPENDENT_VARS)
missing_vars <- all_vars[!all_vars %in% names(sf_data)]

if (length(missing_vars) > 0) {
  stop(paste("Variables manquantes:", paste(missing_vars, collapse=", ")))
}

cat("Toutes les variables sont présentes\n")

# Conversion en Spatial ........................................................
cat("\nConversion en objet Spatial...\n")
sp_data <- as(sf_data, "Spatial")

# Standardisation ..............................................................
if (STANDARDIZE) {
  # Standardisation des variables sur une COPIE
  sp_data_work <- sp_data
  for (var in all_vars) {
    if (var %in% names(sp_data_work@data)) {
      sp_data_work@data[[var]] <- as.numeric(scale(sp_data_work@data[[var]]))
    }
  }
} else {
  sp_data_work <- sp_data
}

# Vérification des NAs .........................................................
cat("\nVérification des valeurs manquantes...\n")
for (var in all_vars) {
  n_na <- sum(is.na(sp_data_work@data[[var]]))
  if (n_na > 0) {
    cat(sprintf("ATTENTION: %d valeurs manquantes pour %s\n", n_na, var))
  }
}

# Suppression des lignes avec NA
complete_cases <- complete.cases(sp_data_work@data[, all_vars])
if (sum(!complete_cases) > 0) {
  cat(sprintf(
    "Suppression de %d lignes avec valeurs manquantes\n",
    sum(!complete_cases)
  ))
  sp_data_work <- sp_data_work[complete_cases, ]
  sf_data <- sf_data[complete_cases, ]
}

# //////////////////////////////////////////////////////////////////////////////
# Modélisation -----------------------------------------------------------------
# Régression OLS globale .......................................................
# Formule
formula <-  as.formula(paste(DEPENDENT_VAR, "~", paste0(INDEPENDENT_VARS, collapse = "+")))
cat(sprintf("\nFormule: %s\n", deparse(formula)))

cat("\n=== REGRESSION OLS GLOBALE ===\n")
ols_model <- lm(formula, data=sp_data_work@data)
print(summary(ols_model))

# Calcul de la Bandwidth .......................................................
if (APPROACH != "None") {
  # Calcul de la bande passante optimale
  cat("Calcul de la bande passante optimale en cours...\n")
  bw <- bw.gwr(
    formula  = formula,
    data     = sp_data_work,
    approach = APPROACH,
    kernel   = KERNEL,
    adaptive = ADAPTATIVE,
    p        = 2,
    longlat  = FALSE
  )
  cat(sprintf("Bande passante optimale (%s): %.4f\n", APPROACH, bw))
} else {
  if (ADAPTATIVE) {
    bw <- NEIGHBORS
  } else {
    bw <- BANDWIDTH
    cat(sprintf("Bande passante utilisée: %.4f\n", bw))
  }
}

# Calcul de la GWR .............................................................
cat("\n=== REGRESSION GWR ===\n")
cat("Calcul en cours (peut prendre plusieurs minutes)...\n")

start_time <- Sys.time()
if (ROBUST) {
  gwr_model <- gwr.robust(
    formula  = formula,
    data     = sp_data_work,
    bw       = bw,
    kernel   = KERNEL,
    adaptive = ADAPTATIVE,
    p        = 2,
    longlat  = FALSE
  )
} else {
  gwr_model <- gwr.basic(
    formula  = formula,
    data     = sp_data_work,
    bw       = bw,
    kernel   = KERNEL,
    adaptive = ADAPTATIVE,
    p        = 2,
    longlat  = FALSE
  )
}
end_time <- Sys.time()
elapsed <- end_time - start_time
cat("GWR procedure completed in", round(elapsed, 3), "seconds\n")

print(gwr_model)

# //////////////////////////////////////////////////////////////////////////////
# Extraction des résultats -----------------------------------------------------
cat("\n=== EXTRACTION DES RESULTATS ===\n")
gwr_sdf <- gwr_model$SDF

# Préparation des données de sortie - COPIE PROPRE de sf_data d'origine
result_sf <- sf_data
result_sf$GWR_yhat <- gwr_sdf$yhat
result_sf$GWR_residual <- gwr_sdf$residual
result_sf$GWR_localR2 <- gwr_sdf$Local_R2

# Ajout des coefficients locaux
coef_cols <- grep(
  "^(?!yhat|residual|Local_R2|sum\\.w|gwr\\.e)",
  names(gwr_sdf@data),
  perl=TRUE,
  value=TRUE
)

for (col in coef_cols) {
  new_name <- paste0("GWR_", col)
  result_sf[[new_name]] <- gwr_sdf@data[[col]]
}

# Statistiques de diagnostic
cat("\n=== STATISTIQUES DE DIAGNOSTIC ===\n")
cat(sprintf("AIC Global OLS: %.2f\n", AIC(ols_model)))
cat(sprintf("AIC GWR: %.2f\n", gwr_model$GW.diagnostic$AICc))
cat(sprintf("R² Global OLS: %.4f\n", summary(ols_model)$r.squared))
cat(sprintf("R² moyen GWR: %.4f\n", mean(result_sf$GWR_localR2, na.rm=TRUE)))
cat(sprintf("R² médian GWR: %.4f\n", median(result_sf$GWR_localR2, na.rm=TRUE)))
cat(sprintf("R² min GWR: %.4f\n", min(result_sf$GWR_localR2, na.rm=TRUE)))
cat(sprintf("R² max GWR: %.4f\n", max(result_sf$GWR_localR2, na.rm=TRUE)))

# Sauvegarde
cat("\nSauvegarde des résultats...\n")
st_write(result_sf, OUTPUT_PATH, delete_dsn=TRUE, quiet=TRUE)

cat("\n=== ANALYSE TERMINEE AVEC SUCCES ===\n")
