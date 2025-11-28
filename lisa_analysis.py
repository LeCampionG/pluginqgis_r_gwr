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
    def write_args(
        shapefile_path: str,
        output_path: str,
        variable: str,
        variable_2: str,
        analysis_type: str,
        standardize: bool,
        contiguity_type: str,
        order: int,
        standardize_weights: float,
        significance: float,
    ) -> list[str]:
        """
        Écrit le script R dans un fichier

        Args:
            shapefile_path: Chemin du shapefile d'entrée
            output_path: Chemin du shapefile de sortie
            analysis_type: type d'analyse ("univariate", "bivariate")
            variable: Première variable
            variable_2: Deuxième variable
            standardize: True pour standardiser les variables
            contiguity_type: type de contiguïté ("queen", "rock")
            order: ordre de la contiguïté
            standardize_weights: True pour standardiser les poids
            significance: valeur alpha de la p-value
        """
        # Convertir les chemins en format lisible par R
        shapefile_path = shapefile_path.replace("\\", "/")
        output_path = output_path.replace("\\", "/")

        # Créer la liste des arguments.
        # WARN: l'ordre importe.
        args = [
            shapefile_path,
            output_path,
            variable,
            variable_2,
            analysis_type,
            standardize,
            contiguity_type,
            order,
            standardize_weights,
            significance,
        ]
        # Convertir les bool en int
        args = [int(x) if isinstance(x, int) else x for x in args]
        # Convertir l'ensemble en string
        args = [str(v) for v in args]

        return args


    @staticmethod
    def run_analysis(
        r_path,
        layer,
        analysis_type,
        variable,
        variable2,
        contiguity_type,
        order,
        standardize_weights,
        significance,
        standardize_var,
    ):
        """
        Exécute l'analyse LISA avec R
        
        Args:
            r_path: Chemin vers Rscript
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

        # Le chemin vers le script R responsable des LISA
        r_script = os.path.join(os.path.dirname(__file__), "r_scripts/lisa.R")

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
                msg = f"Erreur lors de l'export: {error_msg}"
                return None, error_msg, temp_dir

            # Créer le script R
            cmd_args = LISAAnalysisModule.write_args(
                input_shp,
                output_shp,
                variable_shp,
                variable2_shp,
                analysis_type,
                standardize_var,
                contiguity_type,
                order,
                standardize_weights,
                significance,
            )

            # Exécuter le script R
            result = subprocess.run(
                [r_path, r_script] + cmd_args,
                capture_output=True,
                text=True,
                timeout=600  # 10 minutes pour LISA
            )

            if result.returncode != 0:
                msg = f"Erreur R :\n{result.stderr}\n\nSortie stdout:\n{result.stdout}"
                return None, error_msg, temp_dir

            # Charger le résultat
            result_layer = QgsVectorLayer(output_shp, "LISA_Results", "ogr")

            if not result_layer.isValid():
                msg = "Erreur lors du chargement des résultats"
                return None, msg, temp_dir

            return result_layer, result.stdout, temp_dir

        except subprocess.TimeoutExpired:
            return None, "L'analyse a dépassé le temps limite (10 minutes)"
        except Exception as e:
            return None, f"Erreur: {str(e)}"
