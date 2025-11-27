# -*- coding: utf-8 -*-

"""
Module d'analyse LISA (Local Indicators of Spatial Association)
Utilise le package R rgeoda
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


class LISAAnalysisModule:
    """Classe pour effectuer l'analyse LISA avec R - rgeoda"""

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
    def write_r_script_to_file(script_path, shapefile_path, output_path, analysis_type, 
                           variable, variable2, contiguity_type, order, 
                           standardize_weights, significance, standardize_var):
        """
        Écrit le script R pour LISA dans un fichier
        """
    
        # Convertir les chemins Windows en format R
        shapefile_path = shapefile_path.replace("\\", "/")
        output_path = output_path.replace("\\", "/")
    
        # Code de standardisation
        if standardize_var:
            if analysis_type == "univariate":
                standardize_code = f"""
# Standardisation de la variable sur une COPIE
sf_data_work <- sf_data
sf_data_work[["{variable}"]] <- as.numeric(scale(sf_data_work[["{variable}"]]))
"""
            else:
                standardize_code = f"""
# Standardisation des variables sur une COPIE
sf_data_work <- sf_data
sf_data_work[["{variable}"]] <- as.numeric(scale(sf_data_work[["{variable}"]]))
sf_data_work[["{variable2}"]] <- as.numeric(scale(sf_data_work[["{variable2}"]]))
"""
        else:
            standardize_code = """
# Pas de standardisation - créer une copie de travail
sf_data_work <- sf_data
"""
    
        # Code de vérification des variables secondaires
        if analysis_type == "bivariate":
            var2_check = f'''
if (!"{variable2}" %in% names(sf_data_work)) {{
  stop("Variable '{variable2}' non trouvée dans les données")
}}
'''
            var2_na_check = f'''
n_na_var2 <- sum(is.na(sf_data_work[["{variable2}"]]))
if (n_na_var2 > 0) {{
  cat(sprintf("ATTENTION: %d valeurs manquantes pour {variable2}\\n", n_na_var2))
}}
'''
        else:
            var2_check = ""
            var2_na_check = ""
        
        # Code de vérification des cas complets
        if analysis_type == "univariate":
            complete_cases_code = f"complete_cases <- complete.cases(sf_data_work[['{variable}']])"
        else:
            complete_cases_code = f"complete_cases <- complete.cases(sf_data_work[['{variable}']], sf_data_work[['{variable2}']])"
        
        # Code d'extraction des variables
        if analysis_type == "univariate":
            var_extract_code = f'''
var_data <- sf_data_work[["{variable}"]]
cat(sprintf("Variable {variable} extraite: %d valeurs\\n", length(var_data)))
'''
        else:
            var_extract_code = f'''
var_data1 <- sf_data_work[["{variable}"]]
var_data2 <- sf_data_work[["{variable2}"]]
cat(sprintf("Variables {variable} et {variable2} extraites: %d valeurs chacune\\n", length(var_data1)))
'''
        
        # Code de calcul LISA
        if analysis_type == "univariate":
            lisa_calc_code = f'''
# LISA univarié
cat("Analyse univariée de la variable {variable}\\n")
lisa_result <- local_moran(weights, var_df)
cat("Calcul LISA terminé\\n")
'''
        else:
            lisa_calc_code = f'''
# LISA bivarié
cat("Analyse bivariée des variables {variable} et {variable2}\\n")
lisa_result <- local_bimoran(weights, var_df$var1, var_df$var2)
cat("Calcul LISA terminé\\n")
'''
        
        script_content = f"""
# Suppression des warnings
options(warn = -1)

# Chargement des bibliothèques
suppressPackageStartupMessages({{
  library(rgeoda)
  library(sf)
}})

cat("=== DEBUT DE L'ANALYSE LISA ===\\n\\n")

# Lecture des données
cat("Lecture des données...\\n")
sf_data <- st_read("{shapefile_path}", quiet=TRUE)
cat(sprintf("Nombre d'entités: %d\\n", nrow(sf_data)))
cat(sprintf("Nombre de colonnes: %d\\n", ncol(sf_data)))

# Vérification de la variable
cat("\\nVérification de la variable...\\n")
if (!"{variable}" %in% names(sf_data)) {{
  stop("Variable '{variable}' non trouvée dans les données")
}}

{var2_check}

cat("Variable(s) présente(s)\\n")

{standardize_code}

# Vérification des NA
cat("\\nVérification des valeurs manquantes...\\n")
n_na_var1 <- sum(is.na(sf_data_work[["{variable}"]]))
if (n_na_var1 > 0) {{
  cat(sprintf("ATTENTION: %d valeurs manquantes pour {variable}\\n", n_na_var1))
}}

{var2_na_check}

# Suppression des lignes avec NA
{complete_cases_code}
if (sum(!complete_cases) > 0) {{
  cat(sprintf("Suppression de %d lignes avec valeurs manquantes\\n", sum(!complete_cases)))
  sf_data_work <- sf_data_work[complete_cases, ]
  sf_data <- sf_data[complete_cases, ]
}}

# Création de la matrice de poids spatiaux
cat("\\nCréation de la matrice de poids spatiaux...\\n")
cat(sprintf("Type de contiguïté: {contiguity_type}\\n"))
cat(sprintf("Ordre: {order}\\n"))
cat(sprintf("Standardisation: {str(standardize_weights).upper()}\\n"))

# Créer les poids avec rgeoda
if ("{contiguity_type}" == "queen") {{
  weights <- queen_weights(sf_data_work, order = {order}, 
                          include_lower_order = TRUE)
}} else {{
  weights <- rook_weights(sf_data_work, order = {order},
                         include_lower_order = TRUE)
}}

cat(sprintf("Matrice de poids créée avec succès\\n"))

# Préparation des données pour LISA - CRÉER UN DATA.FRAME EXPLICITE
cat("\\nPréparation des données pour LISA...\\n")

# Retirer la géométrie et créer un pur data.frame
data_no_geom <- st_drop_geometry(sf_data_work)

{var_extract_code}

# Créer le data.frame pour rgeoda
{"var_df <- data.frame(var = data_no_geom[['" + variable + "']])" if analysis_type == "univariate" 
 else "var_df <- data.frame(var1 = data_no_geom[['" + variable + "']], var2 = data_no_geom[['" + variable2 + "']])"}

# Calcul LISA
cat("\\n=== CALCUL LISA ===\\n")

{lisa_calc_code}

# Extraction des résultats
cat("\\n=== EXTRACTION DES RESULTATS ===\\n")
cat(sprintf("Seuil de significativité utilisé: {significance}\\n"))

# Extraire directement avec les fonctions rgeoda
lisa_i <- lisa_values(lisa_result)
pvals <- lisa_pvalues(lisa_result)
clusters <- lisa_clusters(lisa_result, cutoff = {significance})
labels <- lisa_labels(lisa_result)

cat(sprintf("Statistiques extraites pour %d observations\\n", length(lisa_i)))
cat(sprintf("Nombre de p-values: %d\\n", length(pvals)))
cat(sprintf("Clusters identifiés\\n"))

# Préparer les résultats - COPIE PROPRE de sf_data d'origine
result_sf <- sf_data

# Ajouter les résultats LISA
result_sf$LISA_I <- lisa_i
result_sf$LISA_pvalue <- pvals
result_sf$LISA_cluster <- clusters

# Créer une colonne catégorielle pour les clusters
cluster_labels <- c("Not significant", "High-High", "Low-Low", "Low-High", "High-Low", "Undefined", "Isolated")
result_sf$LISA_category <- factor(clusters, 
                                   levels = 0:6,
                                   labels = cluster_labels)

# Statistiques descriptives
cat("\\n=== STATISTIQUES DESCRIPTIVES ===\\n")
cat(sprintf("I de Moran moyen: %.4f\\n", mean(lisa_i, na.rm=TRUE)))
cat(sprintf("I min: %.4f\\n", min(lisa_i, na.rm=TRUE)))
cat(sprintf("I max: %.4f\\n", max(lisa_i, na.rm=TRUE)))

cat("\\n=== REPARTITION DES CLUSTERS ===\\n")
cluster_table <- table(result_sf$LISA_category)
for (i in seq_along(cluster_table)) {{
  cat(sprintf("%s: %d (%.1f%%)\\n", 
              names(cluster_table)[i], 
              cluster_table[i],
              100 * cluster_table[i] / nrow(result_sf)))
}}

# Statistiques par type de cluster
cat("\\n=== STATISTIQUES PAR CLUSTER ===\\n")
for (cat_name in names(cluster_table)) {{
  if (cluster_table[cat_name] > 0) {{
    subset_data <- result_sf[result_sf$LISA_category == cat_name, ]
    cat(sprintf("\\n%s:\\n", cat_name))
    cat(sprintf("  I moyen: %.4f\\n", mean(subset_data$LISA_I, na.rm=TRUE)))
    cat(sprintf("  p-value moyenne: %.4f\\n", mean(subset_data$LISA_pvalue, na.rm=TRUE)))
  }}
}}

# Sauvegarde
cat("\\nSauvegarde des résultats...\\n")
st_write(result_sf, "{output_path}", delete_dsn=TRUE, quiet=TRUE)

cat("\\n=== ANALYSE LISA TERMINEE AVEC SUCCES ===\\n")
"""
    
        with open(script_path, 'w', encoding='utf-8') as f:
            f.write(script_content)
    
    @staticmethod
    def run_analysis(layer, analysis_type, variable, variable2, contiguity_type, 
                    order, standardize_weights, significance, standardize_var):
        """
        Exécute l'analyse LISA avec R
        
        Args:
            layer: Couche QGIS à analyser
            analysis_type: "univariate" ou "bivariate"
            variable: Nom de la variable d'analyse
            variable2: Nom de la deuxième variable (pour bivarié)
            contiguity_type: "queen" ou "rook"
            order: Ordre de contiguïté
            standardize_weights: True pour standardiser les poids
            significance: Seuil de significativité
            standardize_var: True pour standardiser la variable
            
        Returns:
            tuple: (result_layer, message) ou (None, error_message)
        """
        temp_dir = tempfile.mkdtemp()
        input_shp = os.path.join(temp_dir, "input.shp")
        output_shp = os.path.join(temp_dir, "output.shp")
        script_path = os.path.join(temp_dir, "lisa_analysis.R")

        try:
            # Créer un mapping sécurisé des noms de champs
            selected_fields = [variable]
            if analysis_type == "bivariate" and variable2:
                selected_fields.append(variable2)
            
            field_mapping = LISAAnalysisModule.create_safe_field_mapping(
                layer.fields(), 
                selected_fields
            )
            
            # Adapter les noms de variables pour le shapefile
            variable_shp = field_mapping[variable]
            variable2_shp = field_mapping[variable2] if variable2 else None
            
            # Afficher le mapping
            print("\n" + "="*60)
            print("MAPPING DES NOMS DE COLONNES (SHAPEFILE) - LISA")
            print("="*60)
            print(f"Variable principale:")
            print(f"  '{variable}' → '{variable_shp}'")
            if variable2:
                print(f"Variable secondaire:")
                print(f"  '{variable2}' → '{variable2_shp}'")
            print("="*60 + "\n")
            
            # Exporter la couche
            success, error_msg = LISAAnalysisModule.export_layer_with_field_mapping(
                layer, input_shp, field_mapping
            )
            
            if not success:
                return None, f"Erreur lors de l'export: {error_msg}"

            # Créer le script R
            LISAAnalysisModule.write_r_script_to_file(
                    script_path, input_shp, output_shp, analysis_type,
                    variable_shp, variable2_shp, contiguity_type, order,
                    standardize_weights, significance, standardize_var
                )

            # Exécuter le script R
            result = subprocess.run(
                ['Rscript', script_path],
                capture_output=True,
                text=True,
                timeout=600  # 10 minutes pour LISA
            )

            if result.returncode != 0:
                return None, f"Erreur R :\n{result.stderr}\n\nSortie stdout:\n{result.stdout}"

            # Charger le résultat
            result_layer = QgsVectorLayer(output_shp, "LISA_Results", "ogr")

            if not result_layer.isValid():
                return None, "Erreur lors du chargement des résultats"

            return result_layer, result.stdout

        except subprocess.TimeoutExpired:
            return None, "L'analyse a dépassé le temps limite (10 minutes)"
        except Exception as e:
            return None, f"Erreur: {str(e)}"
