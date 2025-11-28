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
    def write_args(
        input_path: str,
        output_path: str,
        dependent_var: str,
        independent_vars: str,
        kernel_name: str,
        adaptive: bool,
        standardize: bool,
        bw_approach: str,
        max_iter: int,
        tolerance: float,
        criterion_mgwr: str = "dCVR"
    ) -> list[str]:
        """
        Écrit le script R dans un fichier

        Args:
            input_path: Chemin du shapefile d'entrée
            output_path: Chemin du shapefile de sortie
            dependent_var: Nom de la variable dépendante
            independent_vars: Liste des noms des variables indépendantes
            kernel_name: Type de kernel
            adaptive: True si bande passante adaptative
            standardize: True pour standardiser les variables
            bw_approach: Critère d'optimisation pour la bandwidth ("AIC", "CV")
            max_iter: Nombre maximal d'itérations
            tolerance: Tolérance de convergence
            criterion_mgwr: Critère pour déterminer la convergence de la
                            procédure de back-ftting ("CVR" or "dCVR")
        """
        # Convertir les chemins en format lisible par R
        input_path = input_path.replace("\\", "/")
        output_path = output_path.replace("\\", "/")

        # Si AIC alors AICc
        bw_approach = "AICc" if bw_approach == "AIC" else bw_approach

        # Convertir les X en string
        independent_vars = ",".join(independent_vars)

        # Créer la liste des arguments.
        # WARN: l'ordre importe.
        args = [input_path, output_path, dependent_var, independent_vars,
                kernel_name, adaptive, standardize, bw_approach, criterion_mgwr,
                max_iter, tolerance]
        # Convertir les bool en int
        args = [int(x) if isinstance(x, int) else x for x in args]
        # Convertir l'ensemble en string
        args = [str(v) for v in args]

        return args

    @staticmethod
    def run_analysis(
        r_path,
        layer,
        dependent_var,
        independent_vars,
        kernel_type,
        adaptive,
        standardize,
        criterion,
        max_iter,
        tolerance,
    ):
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

        # Le chemin vers le script R responsable de la GWR
        r_script = os.path.join(os.path.dirname(__file__), "r_scripts/mgwr.R")

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
            # MGWRAnalysisModule.write_r_script_to_file(
                # script_path, input_shp, output_shp, dependent_var_shp, independent_vars_shp,
                # kernel_name, adaptive, standardize, criterion, max_iter, tolerance
            # )
            cmd_args = MGWRAnalysisModule.write_args(
                input_shp,
                output_shp,
                dependent_var_shp,
                independent_vars_shp,
                kernel_name,
                adaptive,
                standardize,
                criterion,
                max_iter,
                tolerance,
            )

            # Exécuter le script R
            result = subprocess.run(
                [r_path, r_script] + cmd_args,
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

            return result_layer, result.stdout, temp_dir

        except subprocess.TimeoutExpired:
            return None, "L'analyse a dépassé le temps limite (60 minutes)"
        except Exception as e:
            return None, f"Erreur: {str(e)}"
