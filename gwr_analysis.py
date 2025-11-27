# -*- coding: utf-8 -*-

"""
Module d'analyse GWR (Geographically Weighted Regression)
Utilise le package R GWmodel
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
import shutil


class GWRAnalysisModule:
    """Classe pour effectuer l'analyse GWR avec R - GWmodel"""

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
        
        # D'abord, traiter les champs sélectionnés pour l'analyse
        for field_name in selected_fields:
            short_name = field_name[:10]
            
            if short_name in used_names:
                counter = 1
                while f"{short_name[:8]}_{counter}" in used_names and counter < 99:
                    counter += 1
                short_name = f"{short_name[:8]}_{counter}"
            
            field_mapping[field_name] = short_name
            used_names.add(short_name)
        
        # Ensuite, traiter les autres champs
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
                               kernel_name, approach, adaptive, bandwidth_value, neighbors, standardize, robust):
        """
        Écrit le script R dans un fichier
        
        Args:
            script_path: Où écrire le script R
            shapefile_path: Chemin du shapefile d'entrée
            output_path: Chemin du shapefile de sortie
            dependent_var: Nom de la variable dépendante
            independent_vars: Liste des noms des variables indépendantes
            kernel_name: Type de kernel
            approach: Méthode de calcul de bande passante ("CV", "AIC", ou None)
            adaptive: True si bande passante adaptative
            bandwidth_value: Valeur manuelle de la bande passante
            neighbors: Nombre de voisins (si adaptive=True et si Bandwith_value renseigner)
            standardize: True pour standardiser les variables
            robust: True pour utiliser gwr.robust
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
        
        if approach:
            bandwidth_code = f"""
# Calcul de la bande passante optimale par {approach}
cat("Calcul de la bande passante optimale en cours...\\n")
bw <- bw.gwr(formula = formula,
             data = sp_data_work,
             approach = "{approach}",
             kernel = "{kernel_name}",
             adaptive = {str(adaptive).upper()},
             p = 2,
             longlat = FALSE)
cat(sprintf("Bande passante optimale ({approach}): %.4f\\n", bw))
"""
        else:
            if adaptive:
                bw_value = neighbors
            else:
                bw_value = bandwidth_value
            bandwidth_code = f"""
# Bande passante manuelle
bw <- {bw_value}
cat(sprintf("Bande passante utilisée: %.4f\\n", bw))
"""
        
        gwr_function = "gwr.robust" if robust else "gwr.basic"
        
        script_content = f"""
# Suppression des warnings
options(warn = -1)

# Chargement des bibliothèques
suppressPackageStartupMessages({{
  library(GWmodel)
  library(sp)
  library(sf)
}})

cat("=== DEBUT DE L'ANALYSE GWR ===\\n\\n")

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

# Formule
formula <- {dependent_var} ~ {formula_vars}
cat(sprintf("\\nFormule: %s\\n", deparse(formula)))

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

{bandwidth_code}

# Régression OLS globale
cat("\\n=== REGRESSION OLS GLOBALE ===\\n")
ols_model <- lm(formula, data=sp_data_work@data)
print(summary(ols_model))

# Calcul de la GWR
cat("\\n=== REGRESSION GWR ===\\n")
cat("Calcul en cours (peut prendre plusieurs minutes)...\\n")

gwr_model <- {gwr_function}(formula = formula,
                             data = sp_data_work,
                             bw = bw,
                             kernel = "{kernel_name}",
                             adaptive = {str(adaptive).upper()},
                             p = 2,
                             longlat = FALSE)

print(gwr_model)

# Extraction des résultats
cat("\\n=== EXTRACTION DES RESULTATS ===\\n")
gwr_sdf <- gwr_model$SDF

# Préparation des données de sortie - COPIE PROPRE de sf_data d'origine
result_sf <- sf_data
result_sf$GWR_yhat <- gwr_sdf$yhat
result_sf$GWR_residual <- gwr_sdf$residual
result_sf$GWR_localR2 <- gwr_sdf$Local_R2

# Ajout des coefficients locaux
coef_cols <- grep("^(?!yhat|residual|Local_R2|sum\\\\.w|gwr\\\\.e)", 
                  names(gwr_sdf@data), 
                  perl=TRUE, 
                  value=TRUE)

for (col in coef_cols) {{
  new_name <- paste0("GWR_", col)
  result_sf[[new_name]] <- gwr_sdf@data[[col]]
}}

# Statistiques de diagnostic
cat("\\n=== STATISTIQUES DE DIAGNOSTIC ===\\n")
cat(sprintf("AIC Global OLS: %.2f\\n", AIC(ols_model)))
cat(sprintf("AIC GWR: %.2f\\n", gwr_model$GW.diagnostic$AICc))
cat(sprintf("R² Global OLS: %.4f\\n", summary(ols_model)$r.squared))
cat(sprintf("R² moyen GWR: %.4f\\n", mean(result_sf$GWR_localR2, na.rm=TRUE)))
cat(sprintf("R² médian GWR: %.4f\\n", median(result_sf$GWR_localR2, na.rm=TRUE)))
cat(sprintf("R² min GWR: %.4f\\n", min(result_sf$GWR_localR2, na.rm=TRUE)))
cat(sprintf("R² max GWR: %.4f\\n", max(result_sf$GWR_localR2, na.rm=TRUE)))

# Sauvegarde
cat("\\nSauvegarde des résultats...\\n")
st_write(result_sf, "{output_path}", delete_dsn=TRUE, quiet=TRUE)

cat("\\n=== ANALYSE TERMINEE AVEC SUCCES ===\\n")
"""
        
        with open(script_path, 'w', encoding='utf-8') as f:
            f.write(script_content)

    @staticmethod
    def run_analysis(r_path, layer, dependent_var, independent_vars, kernel_type,
                    bandwidth_type, bandwidth_value, adaptive, neighbors,
                    standardize, robust):
        """
        Exécute l'analyse GWR avec R
        
        Args:
            r_path: Chemin vers Rscript
            layer: Couche QGIS à analyser
            dependent_var: Nom de la variable dépendante
            independent_vars: Liste des variables indépendantes
            kernel_type: Type de kernel
            bandwidth_type: 1=CV, 2=AIC, 3=manuel
            bandwidth_value: Valeur manuelle de bande passante
            adaptive: True pour bande passante adaptative
            neighbors: Nombre de voisins (si adaptative)
            standardize: True pour standardiser
            robust: True pour GWR robuste
            
        Returns:
            tuple: (result_layer, message) ou (None, error_message)
        """
        temp_dir = tempfile.mkdtemp()
        input_shp = os.path.join(temp_dir, "input.shp")
        output_shp = os.path.join(temp_dir, "output.shp")
        script_path = os.path.join(temp_dir, "gwr_analysis.R")

        try:
            # Créer un mapping sécurisé des noms de champs
            all_selected_fields = [dependent_var] + independent_vars
            field_mapping = GWRAnalysisModule.create_safe_field_mapping(
                layer.fields(), 
                all_selected_fields
            )
            
            # Adapter les noms de variables pour le shapefile
            dependent_var_shp = field_mapping[dependent_var]
            independent_vars_shp = [field_mapping[v] for v in independent_vars]
            
            # Afficher le mapping
            print("\n" + "="*60)
            print("MAPPING DES NOMS DE COLONNES (SHAPEFILE) - GWR")
            print("="*60)
            print(f"Variable dépendante:")
            print(f"  '{dependent_var}' → '{dependent_var_shp}'")
            print(f"\nVariables indépendantes:")
            for orig, short in zip(independent_vars, independent_vars_shp):
                print(f"  '{orig}' → '{short}'")
            print("="*60 + "\n")
            
            # Exporter la couche
            success, error_msg = GWRAnalysisModule.export_layer_with_field_mapping(
                layer, input_shp, field_mapping
            )
            
            if not success:
                return None, f"Erreur lors de l'export: {error_msg}"

            # Déterminer l'approche
            if bandwidth_type == 1:
                approach = "CV"
            elif bandwidth_type == 2:
                approach = "AIC"
            else:
                approach = None

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
            GWRAnalysisModule.write_r_script_to_file(
                script_path, input_shp, output_shp, dependent_var_shp, independent_vars_shp,
                kernel_name, approach, adaptive, bandwidth_value, neighbors, standardize, robust
            )

            # Exécuter le script R
            result = subprocess.run(
                [r_path, script_path],
                capture_output=True,
                text=True,
                timeout=1800  # 30 minutes
            )

            if result.returncode != 0:
                return None, f"Erreur R :\n{result.stderr}\n\nSortie stdout:\n{result.stdout}"

            # Charger le résultat
            result_layer = QgsVectorLayer(output_shp, "GWR_Results", "ogr")

            if not result_layer.isValid():
                return None, "Erreur lors du chargement des résultats"

            return result_layer, result.stdout, temp_dir

        except subprocess.TimeoutExpired:
            return None, "L'analyse a dépassé le temps limite (30 minutes)"
        except Exception as e:
            return None, f"Erreur: {str(e)}"
