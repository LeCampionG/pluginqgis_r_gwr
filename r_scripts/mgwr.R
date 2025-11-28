################################################################################
# Script Name: mgwr.R
# Description: Fait tourner une MGWR avec les paramètres fournis en entrée
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
ADAPTATIVE       <- as.logical(as.integer(args[6]))
STANDARDIZE      <- as.logical(as.integer(args[7]))
BW_APPROACH      <- args[8]
CRITERION_MGWR   <- args[9]
MAX_ITER         <- as.integer(args[10])
TOLERANCE        <- as.numeric(args[11])

cat("=== Résumé des paramètres fournis ===\n")
cat("Chemin d'entrée                     :", INPUT_PATH, "\n")
cat("Chemin de sortie                    :", OUTPUT_PATH, "\n")
cat("Variable dépendante                 :", DEPENDENT_VAR, "\n")
cat("Variables indépendantes             :", paste(INDEPENDENT_VARS, collapse = ", "), "\n")
cat("Noyau choisi                        :", KERNEL, "\n")
cat("Bande passante adaptative           :", ADAPTATIVE, "\n")
cat("Standardisation                     :", STANDARDIZE, "\n")
cat("Approche pour la Bande Passante     :", BW_APPROACH, "\n")
cat("Critère d'optimisation pour la MGWR :", CRITERION_MGWR, "\n")
cat("Nombre d'itérations max             :", MAX_ITER, "\n")
cat("Tolérance                           :", TOLERANCE, "\n")
cat("=====================================\n")

# //////////////////////////////////////////////////////////////////////////////
# Lecture des données ----------------------------------------------------------
cat("=== DEBUT DE L'ANALYSE MGWR ===\n\n")

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
# Formule
formula <-  as.formula(paste(DEPENDENT_VAR, "~", paste0(INDEPENDENT_VARS, collapse = "+")))
cat(sprintf("\nFormule: %s\n", deparse(formula)))

# Régression OLS globale .......................................................
cat("\n=== REGRESSION OLS GLOBALE ===\n")
ols_model <- lm(formula, data=sp_data_work@data)
print(summary(ols_model))

# Calcul de la Bandwidth .......................................................
cat("\n=== CALCUL DE LA BANDE PASSANTE GLOBALE INITIALE ===\n")
bw_global <- bw.gwr(
  formula  = formula,
  data     = sp_data_work,
  approach = BW_APPROACH,
  kernel   = KERNEL,
  adaptive = ADAPTATIVE,
  p        = 2,
  longlat  = FALSE
)

cat(sprintf("BW globale: %.2f\n", bw_global))

# Créer un vecteur de BW initiales (une par variable X + intercept)
n_vars <- length(INDEPENDENT_VARS) + 1
bw_init <- rep(bw_global, n_vars)
var_names <- c("Intercept", INDEPENDENT_VARS)

cat(sprintf("Nombre de variables (avec intercept): %d\n", n_vars))

# Calcul de la MGWR avec gwr.multiscale
cat("\n=== REGRESSION MGWR (gwr.multiscale) ===\n")
cat("Calcul en cours (peut prendre plusieurs minutes)...\n")
cat("Paramètres:\n")
cat("  - Criterion        : ", CRITERION_MGWR, "\n")
cat("  - Kernel           : ", KERNEL, "\n")
cat("  - Adaptatif        : ", ADAPTATIVE, "\n")
cat("  - Max itérations   : ", MAX_ITER, "\n")
cat("  - Tolérance        : ", TOLERANCE, "\n\n")

tryCatch({
  # Calculer la matrice de distances
  cat("Calcul de la matrice de distances...\n")
  dMat <- gw.dist(dp.locat = coordinates(sp_data_work))
  cat(sprintf("Matrice de distances: %d x %d\n", nrow(dMat), ncol(dMat)))

  # Appel à gwr.multiscale
  cat("Lancement de gwr.multiscale...\n")
  cat("Configuration:\n")
  cat(sprintf("  - Nombre de variables: %d\n", n_vars))
  cat(sprintf("  - Toutes les variables utilisent la même matrice de distances\n"))

  mgwr_model <- gwr.multiscale(
    formula            = formula,
    data               = sp_data_work,
    criterion          = CRITERION_MGWR,
    kernel             = KERNEL,
    adaptive           = ADAPTATIVE,
    bws0               = bw_init,
    bw.seled           = rep(TRUE, n_vars),
    dMats              = list(dMat),
    predictor.centered = rep(TRUE, n_vars),
    var.dMat.indx      = rep(1, n_vars),
    max.iterations     = MAX_ITER,
    threshold          = TOLERANCE,
    verbose            = FALSE
  )

  cat("\n=== RESULTATS MGWR ===\n")
  print(mgwr_model)

  # Extraction des résultats
  cat("\n=== EXTRACTION DES RESULTATS ===\n")

  # Prédictions et résidus - COPIE PROPRE de sf_data d'origine
  result_sf <- sf_data
  result_sf$MGWR_yhat <- mgwr_model$SDF$yhat
  result_sf$MGWR_residual <- mgwr_model$SDF$residual
  result_sf$MGWR_localR2 <- mgwr_model$SDF$Local_R2

  # Coefficients locaux
  # NOTE: il n'y a pas de R2 locaux comme avec gwr.basic ou gwr.robust
  cat("\nExtraction des coefficients locaux...\n")
  coef_names <- names(mgwr_model$SDF@data)
  coef_names <- coef_names[!coef_names %in% c(
    "yhat",
    "residual",
    "Local_R2",
    "sum.w",
    "gwr.e"
  )]

  for (coef_name in coef_names) {
    new_name <- paste0("MGWR_", coef_name)
    result_sf[[new_name]] <- mgwr_model$SDF@data[[coef_name]]
  }

  # Bandes passantes optimales
  cat("\n=== BANDES PASSANTES OPTIMALES ===\n")
  if (!is.null(mgwr_model$GW.arguments$bws)) {
    bws <- mgwr_model$GW.arguments$bws

    # Créer un data.frame avec les BW
    bw_df <- data.frame(Variable = var_names, Bandwidth = bws)
    print(bw_df)

    # Sauvegarder les bandes passantes dans le résultat
    for (i in seq_along(var_names)) {
      bw_col_name <- paste0("MGWR_BW_", gsub("[^A-Za-z0-9_]", "_", var_names[i]))
      result_sf[[bw_col_name]] <- rep(bws[i], nrow(result_sf))
    }
  }

  # Statistiques de diagnostic
  cat("\n=== STATISTIQUES DE DIAGNOSTIC ===\n")
  cat(sprintf("AIC Global OLS: %.2f\n", AIC(ols_model)))

  if (!is.null(mgwr_model$GW.diagnostic$AICc)) {
    cat(sprintf("AICc MGWR: %.2f\n", mgwr_model$GW.diagnostic$AICc))
  }

  cat(sprintf("R² Global OLS: %.4f\n", summary(ols_model)$r.squared))

  cat(sprintf("R² moyen MGWR: %.4f\n", mean(result_sf$MGWR_localR2, na.rm=TRUE)))
  cat(sprintf("R² médian MGWR: %.4f\n", median(result_sf$MGWR_localR2, na.rm=TRUE)))
  cat(sprintf("R² min MGWR: %.4f\n", min(result_sf$MGWR_localR2, na.rm=TRUE)))
  cat(sprintf("R² max MGWR: %.4f\n", max(result_sf$MGWR_localR2, na.rm=TRUE)))

  # Sauvegarde
  cat("\nSauvegarde des résultats...\n")
  st_write(result_sf, OUTPUT_PATH, delete_dsn=TRUE, quiet=TRUE)

  cat("\n=== ANALYSE MGWR TERMINEE AVEC SUCCES ===\n")
}, error = function(e) {
  cat("\n=== ERREUR LORS DU CALCUL MGWR ===\n")
  cat(sprintf("Message d'erreur: %s\n", e$message))
  cat(sprintf("\nTraceback:\n"))
  print(traceback())
  stop(e)
})
