# -*- coding: utf-8 -*-

"""
Module d'analyse MGWR (Multiscale Geographically Weighted Regression)
Utilise le package R GWmodel - fonction gwr.multiscale
"""

from qgis.core import (
    QgsVectorLayer,
    QgsVectorFileWriter,
    QgsFields,
    QgsField,
    QgsFeature
)

import subprocess
import os
import tempfile


class MGWRAnalysisModule:
    """Classe pour effectuer l'analyse MGWR avec R - GWmodel::gwr.multiscale"""

    @staticmethod
    def create_safe_field_mapping(layer_fields, selected_fields):
        """
        Crée un mapping sécurisé des noms de champs pour éviter les collisions
        lors de la troncature à 10 caractères (limite shapefile)
        
        Args:
            layer_fields: QgsFields de la couche
            selected_fields: Liste des noms de champs sélectionnés pour l'analyse
            
        Returns:
            dict: {nom_original: nom_court_unique}
        """
        field_mapping = {}
        used_names = set()
        
        for field_name in selected_fields:
            short_name = field_name[:10]
            
            if short_name in used_names:
                counter = 1
                while f"{short_name[:8]}_{counter}" in used_names and counter < 99:
                    counter += 1
                short_name = f"{short_name[:8]}_{counter}"
            
            field_mapping[field_name] = short_name
            used_names.add(short_name)
        
        for field in layer_fields:
            field_name = field.name()
            if field_name not in field_mapping:
                short_name = field_name[:10]
                
                if short_name in used_names:
                    counter = 1
                    while f"{short_name[:8]}_{counter}" in used_names and counter < 99:
                        counter += 1
                    short_name = f"{short_name[:8]}_{counter}"
                
                field_mapping[field_name] = short_name
                used_names.add(short_name)
        
        return field_mapping

    @staticmethod
    def export_layer_with_field_mapping(layer, output_path, field_mapping):
        """
        Exporte une couche en renommant les champs selon le mapping fourni
        
        Args:
            layer: QgsVectorLayer à exporter
            output_path: Chemin du fichier de sortie
            field_mapping: dict {nom_original: nom_court}
            
        Returns:
            tuple: (success, error_message)
        """
        new_fields = QgsFields()
        for field in layer.fields():
            new_field = QgsField(field)
            original_name = field.name()
            if original_name in field_mapping:
                new_field.setName(field_mapping[original_name])
            new_fields.append(new_field)
        
        save_options = QgsVectorFileWriter.SaveVectorOptions()
        save_options.driverName = "ESRI Shapefile"
        save_options.fileEncoding = "UTF-8"
        
        writer = QgsVectorFileWriter.create(
            output_path,
            new_fields,
            layer.wkbType(),
            layer.crs(),
            layer.transformContext(),
            save_options
        )
        
        if writer.hasError() != QgsVectorFileWriter.NoError:
            return False, f"Erreur lors de la création du shapefile: {writer.errorMessage()}"
        
        for feature in layer.getFeatures():
            new_feature = QgsFeature(new_fields)
            new_feature.setGeometry(feature.geometry())
            
            for i, field in enumerate(layer.fields()):
                new_feature.setAttribute(i, feature.attribute(field.name()))
            
            writer.addFeature(new_feature)
        
        del writer
        
        return True, None
    
    @staticmethod
    def write_r_script_to_file(script_path, shapefile_path, output_path, dependent_var, independent_vars,
                               kernel_name, adaptive, standardize, criterion, max_iter, tolerance):
        """
        Écrit le script R pour MGWR dans un fichier
        Utilise GWmodel::gwr.multiscale()
        
        Args:
            script_path: Où écrire le script R
            shapefile_path: Chemin du shapefile d'entrée
            output_path: Chemin du shapefile de sortie
            dependent_var: Nom de la variable dépendante
            independent_vars: Liste des noms des variables indépendantes
            kernel_name: Type de kernel ("gaussian" ou "bisquare")
            adaptive: True si bande passante adaptative
            standardize: True pour standardiser les variables
            criterion: Critère d'optimisation ("AICc", "AIC", "BIC", "CV")
            max_iter: Nombre maximal d'itérations
            tolerance: Tolérance de convergence
        """
        
        # Convertir les chemins Windows en format R
        shapefile_path = shapefile_path.replace("\\", "/")
        output_path = output_path.replace("\\", "/")
        
        formula_vars = ' + '.join(independent_vars)
        
        if standardize:
            vars_list = ', '.join([f'"{v}"' for v in [dependent_var] + independent_vars])
            standardize_code = f"""
# Standardisation des variables sur une COPIE
vars_to_scale <- c({vars_list})
sp_data_work <- sp_data
for (var in vars_to_scale) {{
  if (var %in% names(sp_data_work@data)) {{
    sp_data_work@data[[var]] <- as.numeric(scale(sp_data_work@data[[var]]))
  }}
}}
"""
        else:
            standardize_code = """
# Pas de standardisation - créer une copie de travail
sp_data_work <- sp_data
"""
        
        # Pour gwr.multiscale : criterion = "CVR" ou "dCVR", approach = "AICc" ou "CV"
        # On utilise CVR comme criterion et on adapte approach selon le choix utilisateur
        if criterion == "AICc":
            approach_gwmodel = "AICc"
        elif criterion == "AIC":
            approach_gwmodel = "AICc"  # Pas d'AIC seul, on utilise AICc
        elif criterion == "CV":
            approach_gwmodel = "CV"
        else:
            approach_gwmodel = "AICc"

        # criterion pour gwr.multiscale dCVR utilisé par défault
        criterion_gwmodel = "dCVR"
        
        script_content = f"""
# Suppression des warnings
options(warn = -1)

# Chargement des bibliothèques
suppressPackageStartupMessages({{
  library(GWmodel)
  library(sp)
  library(sf)
}})

cat("=== DEBUT DE L'ANALYSE MGWR ===\\n\\n")

# Lecture des données
cat("Lecture des données...\\n")
sf_data <- st_read("{shapefile_path}", quiet=TRUE)
cat(sprintf("Nombre d'entités: %d\\n", nrow(sf_data)))
cat(sprintf("Nombre de colonnes: %d\\n", ncol(sf_data)))

# Vérification des variables
cat("\\nVérification des variables...\\n")
required_vars <- c("{dependent_var}", {', '.join([f'"{v}"' for v in independent_vars])})
missing_vars <- required_vars[!required_vars %in% names(sf_data)]

if (length(missing_vars) > 0) {{
  stop(paste("Variables manquantes:", paste(missing_vars, collapse=", ")))
}}

cat("Toutes les variables sont présentes\\n")

# Conversion en Spatial
cat("\\nConversion en objet Spatial...\\n")
sp_data <- as(sf_data, "Spatial")

{standardize_code}

# Vérification des NA
cat("\\nVérification des valeurs manquantes...\\n")
for (var in required_vars) {{
  n_na <- sum(is.na(sp_data_work@data[[var]]))
  if (n_na > 0) {{
    cat(sprintf("ATTENTION: %d valeurs manquantes pour %s\\n", n_na, var))
  }}
}}

# Suppression des lignes avec NA
complete_cases <- complete.cases(sp_data_work@data[, required_vars])
if (sum(!complete_cases) > 0) {{
  cat(sprintf("Suppression de %d lignes avec valeurs manquantes\\n", sum(!complete_cases)))
  sp_data_work <- sp_data_work[complete_cases, ]
  sf_data <- sf_data[complete_cases, ]
}}

# Formule
formula <- {dependent_var} ~ {formula_vars}
cat(sprintf("\\nFormule: %s\\n", deparse(formula)))

# Régression OLS globale
cat("\\n=== REGRESSION OLS GLOBALE ===\\n")
ols_model <- lm(formula, data=sp_data_work@data)
print(summary(ols_model))

# Calcul de la bande passante globale initiale
cat("\\n=== CALCUL DE LA BANDE PASSANTE GLOBALE INITIALE ===\\n")
bw_global <- bw.gwr(formula = formula,
                    data = sp_data_work,
                    approach = "{approach_gwmodel}",
                    kernel = "{kernel_name}",
                    adaptive = {str(adaptive).upper()},
                    p = 2,
                    longlat = FALSE)

cat(sprintf("BW globale: %.2f\\n", bw_global))

# Créer un vecteur de BW initiales (une par variable X + intercept)
n_vars <- length(c({', '.join([f'"{v}"' for v in independent_vars])})) + 1  # +1 pour intercept
bw_init <- rep(bw_global, n_vars)
var_names <- c("Intercept", {', '.join([f'"{v}"' for v in independent_vars])})

cat(sprintf("Nombre de variables (avec intercept): %d\\n", n_vars))

# Calcul de la MGWR avec gwr.multiscale
cat("\\n=== REGRESSION MGWR (gwr.multiscale) ===\\n")
cat("Calcul en cours (peut prendre plusieurs minutes)...\\n")
cat(sprintf("Paramètres:\\n"))
cat(sprintf("  - Criterion: {criterion_gwmodel}\\n"))
cat(sprintf("  - Kernel: {kernel_name}\\n"))
cat(sprintf("  - Adaptatif: {str(adaptive).upper()}\\n"))
cat(sprintf("  - Max itérations: {max_iter}\\n"))
cat(sprintf("  - Tolérance: {tolerance}\\n\\n"))

tryCatch({{
  # Calculer la matrice de distances
  cat("Calcul de la matrice de distances...\\n")
  dMat <- gw.dist(dp.locat = coordinates(sp_data_work))
  cat(sprintf("Matrice de distances: %d x %d\\n", nrow(dMat), ncol(dMat)))
  
  # Appel à gwr.multiscale
  cat("Lancement de gwr.multiscale...\\n")
  cat("Configuration:\\n")
  cat(sprintf("  - Nombre de variables: %d\\n", n_vars))
  cat(sprintf("  - Toutes les variables utilisent la même matrice de distances\\n"))
  
  mgwr_model <- gwr.multiscale(formula = formula,
                                data = sp_data_work,
                                criterion = "{criterion_gwmodel}",
                                kernel = "{kernel_name}",
                                adaptive = {str(adaptive).upper()},
                                bws0 = bw_init,
                                bw.seled = rep(TRUE, n_vars),
                                dMats = list(dMat),
                                predictor.centered = rep(TRUE, n_vars),
                                var.dMat.indx = rep(1, n_vars),
                                max.iterations = {max_iter},
                                threshold = {tolerance},
                                verbose = FALSE)
  
  cat("\\n=== RESULTATS MGWR ===\\n")
  print(mgwr_model)
  
  # Extraction des résultats
  cat("\\n=== EXTRACTION DES RESULTATS ===\\n")
  
  # Prédictions et résidus - COPIE PROPRE de sf_data d'origine
  result_sf <- sf_data
  result_sf$MGWR_yhat <- mgwr_model$SDF$yhat
  result_sf$MGWR_residual <- mgwr_model$SDF$residual
  result_sf$MGWR_localR2 <- mgwr_model$SDF$Local_R2
  
  # Coefficients locaux
  cat("\\nExtraction des coefficients locaux...\\n")
  coef_names <- names(mgwr_model$SDF@data)
  coef_names <- coef_names[!coef_names %in% c("yhat", "residual", "Local_R2", "sum.w", "gwr.e")]
  
  for (coef_name in coef_names) {{
    new_name <- paste0("MGWR_", coef_name)
    result_sf[[new_name]] <- mgwr_model$SDF@data[[coef_name]]
  }}
  
  # Bandes passantes optimales
  cat("\\n=== BANDES PASSANTES OPTIMALES ===\\n")
  if (!is.null(mgwr_model$GW.arguments$bws)) {{
    bws <- mgwr_model$GW.arguments$bws
    
    # Créer un data.frame avec les BW
    bw_df <- data.frame(
      Variable = var_names,
      Bandwidth = bws
    )
    print(bw_df)
    
    # Sauvegarder les bandes passantes dans le résultat
    for (i in seq_along(var_names)) {{
      bw_col_name <- paste0("MGWR_BW_", gsub("[^A-Za-z0-9_]", "_", var_names[i]))
      result_sf[[bw_col_name]] <- rep(bws[i], nrow(result_sf))
    }}
  }}
  
  # Statistiques de diagnostic
  cat("\\n=== STATISTIQUES DE DIAGNOSTIC ===\\n")
  cat(sprintf("AIC Global OLS: %.2f\\n", AIC(ols_model)))
  
  if (!is.null(mgwr_model$GW.diagnostic$AICc)) {{
    cat(sprintf("AICc MGWR: %.2f\\n", mgwr_model$GW.diagnostic$AICc))
  }}
  
  cat(sprintf("R² Global OLS: %.4f\\n", summary(ols_model)$r.squared))
  cat(sprintf("R² moyen MGWR: %.4f\\n", mean(result_sf$MGWR_localR2, na.rm=TRUE)))
  cat(sprintf("R² médian MGWR: %.4f\\n", median(result_sf$MGWR_localR2, na.rm=TRUE)))
  cat(sprintf("R² min MGWR: %.4f\\n", min(result_sf$MGWR_localR2, na.rm=TRUE)))
  cat(sprintf("R² max MGWR: %.4f\\n", max(result_sf$MGWR_localR2, na.rm=TRUE)))
  
  # Sauvegarde
  cat("\\nSauvegarde des résultats...\\n")
  st_write(result_sf, "{output_path}", delete_dsn=TRUE, quiet=TRUE)
  
  cat("\\n=== ANALYSE MGWR TERMINEE AVEC SUCCES ===\\n")
  
}}, error = function(e) {{
  cat("\\n=== ERREUR LORS DU CALCUL MGWR ===\\n")
  cat(sprintf("Message d'erreur: %s\\n", e$message))
  cat(sprintf("\\nTraceback:\\n"))
  print(traceback())
  stop(e)
}})
"""
        
        with open(script_path, 'w', encoding='utf-8') as f:
            f.write(script_content)

    @staticmethod
    def run_analysis(r_path, layer, dependent_var, independent_vars, kernel_type,
                    adaptive, standardize, criterion, max_iter, tolerance):
        """
        Exécute l'analyse MGWR avec R - GWmodel::gwr.multiscale
        
        Args:
            r_path: Chemin vers Rscript
            layer: Couche QGIS à analyser
            dependent_var: Nom de la variable dépendante
            independent_vars: Liste des variables indépendantes
            kernel_type: Type de kernel ("gaussian" ou "bisquare")
            adaptive: True pour bande passante adaptative
            standardize: True pour standardiser
            criterion: Critère d'optimisation
            max_iter: Nombre maximal d'itérations
            tolerance: Tolérance de convergence
            
        Returns:
            tuple: (result_layer, message) ou (None, error_message)
        """
        temp_dir = tempfile.mkdtemp()
        input_shp = os.path.join(temp_dir, "input.shp")
        output_shp = os.path.join(temp_dir, "output.shp")
        script_path = os.path.join(temp_dir, "mgwr_analysis.R")

        try:
            # Créer un mapping sécurisé des noms de champs
            all_selected_fields = [dependent_var] + independent_vars
            field_mapping = MGWRAnalysisModule.create_safe_field_mapping(
                layer.fields(), 
                all_selected_fields
            )
            
            # Adapter les noms de variables pour le shapefile
            dependent_var_shp = field_mapping[dependent_var]
            independent_vars_shp = [field_mapping[v] for v in independent_vars]
            
            # Afficher le mapping
            print("\n" + "="*60)
            print("MAPPING DES NOMS DE COLONNES (SHAPEFILE) - MGWR")
            print("="*60)
            print(f"Variable dépendante:")
            print(f"  '{dependent_var}' → '{dependent_var_shp}'")
            print(f"\nVariables indépendantes:")
            for orig, short in zip(independent_vars, independent_vars_shp):
                print(f"  '{orig}' → '{short}'")
            print("="*60 + "\n")
            
            # Exporter la couche
            success, error_msg = MGWRAnalysisModule.export_layer_with_field_mapping(
                layer, input_shp, field_mapping
            )
            
            if not success:
                return None, f"Erreur lors de l'export: {error_msg}"

            # Nom du kernel
            kernel_names = {
                "gaussian": "gaussian",
                "bisquare": "bisquare",
                "tricube": "tricube",
                "exponential": "exponential",
                "boxcar": "boxcar"
            }
            kernel_name = kernel_names.get(kernel_type, "gaussian")

            # Créer le script R
            MGWRAnalysisModule.write_r_script_to_file(
                script_path, input_shp, output_shp, dependent_var_shp, independent_vars_shp,
                kernel_name, adaptive, standardize, criterion, max_iter, tolerance
            )

            # Exécuter le script R
            result = subprocess.run(
                [r_path, script_path],
                capture_output=True,
                text=True,
                timeout=3600  # 60 minutes pour MGWR
            )

            if result.returncode != 0:
                return None, f"Erreur R :\n{result.stderr}\n\nSortie stdout:\n{result.stdout}"

            # Charger le résultat
            result_layer = QgsVectorLayer(output_shp, "MGWR_Results", "ogr")

            if not result_layer.isValid():
                return None, "Erreur lors du chargement des résultats"

            return result_layer, result.stdout

        except subprocess.TimeoutExpired:
            return None, "L'analyse a dépassé le temps limite (60 minutes)"
        except Exception as e:
            return None, f"Erreur: {str(e)}"
