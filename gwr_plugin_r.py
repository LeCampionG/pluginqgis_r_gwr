# -*- coding: utf-8 -*-

"""
Plugin QGIS pour la Régression Géographiquement Pondérée (GWR), la MGWR et les LISA en utilisant R
ARCHITECTURE DU PLUGIN:
1. ModelSelectionDialog() pour choisir entre les analyses
2. GWRAnalysisDialog() pour l'interface GWR
3. MGWRAnalysisDialog() pour l'interface MGWR
4. LISAAnalysisDialog() pour l'interface LISA
5. run_gwr_analysis() pour l'exécution GWR (appelle gwr_analysis.py)
6. run_mgwr_analysis() pour l'exécution MGWR (appelle mgwr_analysis.py)
7. run_lisa_analysis() pour l'exécution LISA (appelle lisa_analysis.py)
8. GWRPlugin : Point d'entrée QGIS
"""

from qgis.PyQt.QtCore import QVariant
from qgis.PyQt.QtWidgets import (
    QAction, QDialog, QVBoxLayout, QComboBox, QPushButton, QLabel,
    QSpinBox, QDoubleSpinBox, QMessageBox, QRadioButton, QButtonGroup,
    QGroupBox, QHBoxLayout, QCheckBox, QProgressDialog, QListWidget,
    QAbstractItemView, QFileDialog
)

from qgis.core import (
    QgsProject, 
    QgsVectorLayer, 
    QgsVectorFileWriter, 
    QgsField
)

import subprocess
import os
import shutil
import configparser

# Importer les modules d'analyse séparés
from .gwr_analysis import GWRAnalysisModule
from .mgwr_analysis import MGWRAnalysisModule
from .lisa_analysis import LISAAnalysisModule

CONFIG_FILE = os.path.join(os.path.dirname(__file__), "config.ini")

# =====================================================
# Dialogue de sélection des méthodes
# =====================================================

class ModelSelectionDialog(QDialog):
    """Dialogue initial pour choisir entre les différentes méthodes (GWR/LISA/MGWR)"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Sélection du modèle de régression")
        self.setMinimumWidth(400)
        self.selected_model = None
        
        layout = QVBoxLayout()
        
        title = QLabel("<h2>Choisissez le type d'analyse</h2>")
        layout.addWidget(title)
        layout.addSpacing(10)
        
        # Description GWR
        gwr_group = QGroupBox("GWR - Régression Géographiquement Pondérée")
        gwr_layout = QVBoxLayout()
        gwr_desc = QLabel(
            "Utilise une seule bande passante pour toutes les variables.\n"
            "Utilise le package GWmodel"
        )
        gwr_desc.setWordWrap(True)
        gwr_layout.addWidget(gwr_desc)
        
        self.gwr_button = QPushButton("Lancer GWR")
        self.gwr_button.setMinimumHeight(40)
        self.gwr_button.clicked.connect(self.select_gwr)
        gwr_layout.addWidget(self.gwr_button)
        
        gwr_group.setLayout(gwr_layout)
        layout.addWidget(gwr_group)
        
        layout.addSpacing(10)
        
        # Description MGWR
        mgwr_group = QGroupBox("MGWR - Régression Multi-échelle")
        mgwr_layout = QVBoxLayout()
        mgwr_desc = QLabel(
            "Calcule une bande passante optimale pour chaque variable.\n"
            "Permet de capturer des processus à différentes échelles spatiales.\n"
            "Utilise le package GWmodel"
        )
        mgwr_desc.setWordWrap(True)
        mgwr_layout.addWidget(mgwr_desc)
        
        self.mgwr_button = QPushButton("Lancer MGWR")
        self.mgwr_button.setMinimumHeight(40)
        self.mgwr_button.clicked.connect(self.select_mgwr)
        mgwr_layout.addWidget(self.mgwr_button)
        
        mgwr_group.setLayout(mgwr_layout)
        layout.addWidget(mgwr_group)
        
        layout.addSpacing(10)
        
        # Description LISA
        lisa_group = QGroupBox("LISA - Indicateurs Locaux d'Association Spatiale")
        lisa_layout = QVBoxLayout()
        lisa_desc = QLabel(
            "Identifie les clusters spatiaux et les outliers spatiaux.\n"
            "Analyse univariée ou bivariée de l'autocorrélation spatiale locale.\n"
            "Nécessite le package rgeoda dans R."
        )
        lisa_desc.setWordWrap(True)
        lisa_layout.addWidget(lisa_desc)
        
        self.lisa_button = QPushButton("Lancer LISA")
        self.lisa_button.setMinimumHeight(40)
        self.lisa_button.clicked.connect(self.select_lisa)
        lisa_layout.addWidget(self.lisa_button)
        
        lisa_group.setLayout(lisa_layout)
        layout.addWidget(lisa_group)
        
        layout.addSpacing(10)
        
        cancel_button = QPushButton("Annuler")
        cancel_button.clicked.connect(self.reject)
        layout.addWidget(cancel_button)
        
        self.setLayout(layout)
    
    def select_gwr(self):
        self.selected_model = "GWR"
        self.accept()
    
    def select_mgwr(self):
        self.selected_model = "MGWR"
        self.accept()
    
    def select_lisa(self):
        self.selected_model = "LISA"
        self.accept()

# =====================================================
# Dialogue GWR (code original)
# =====================================================

class GWRAnalysisDialog(QDialog):
    """Interface pour configurer l'analyse GWR"""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Régression Géographiquement Pondérée (GWR) avec R - GWmodel")
        self.setMinimumWidth(500)
        self.setMinimumHeight(650)
        self.r_path = self.load_r_path()

        layout = QVBoxLayout()

        # Chemin Rscript
        layout.addWidget(QLabel("Chemin vers Rscript :"))
        path_layout = QHBoxLayout()
        self.r_path_label = QLabel(self.r_path if self.r_path else "(non défini)")
        self.choose_r_path_button = QPushButton("Sélectionner Rscript")
        self.choose_r_path_button.clicked.connect(self.choose_r_path)
        path_layout.addWidget(self.r_path_label)
        path_layout.addWidget(self.choose_r_path_button)
        layout.addLayout(path_layout)
        layout.addSpacing(10)

        # Sélection couche
        layout.addWidget(QLabel("Couche vectorielle:"))
        self.layer_combo = QComboBox()
        self.populate_layers()
        layout.addWidget(self.layer_combo)

        # Variable dépendante
        layout.addWidget(QLabel("Variable dépendante (Y):"))
        self.dependent_combo = QComboBox()
        self.layer_combo.currentIndexChanged.connect(self.populate_fields)
        self.populate_fields()
        layout.addWidget(self.dependent_combo)

        # Variables indépendantes
        layout.addWidget(QLabel("Variables indépendantes (X) - sélection multiple:"))
        self.independent_list = QListWidget()
        self.independent_list.setSelectionMode(QAbstractItemView.MultiSelection)
        self.dependent_combo.currentIndexChanged.connect(self.populate_independent_vars)
        self.populate_independent_vars()
        layout.addWidget(self.independent_list)

        # Kernel
        kernel_group = QGroupBox("Fonction de pondération (kernel)")
        kernel_layout = QVBoxLayout()
        self.kernel_button_group = QButtonGroup()

        self.gaussian_radio = QRadioButton("Gaussien (recommandé)")
        self.gaussian_radio.setChecked(True)
        self.kernel_button_group.addButton(self.gaussian_radio, 1)
        kernel_layout.addWidget(self.gaussian_radio)

        self.bisquare_radio = QRadioButton("Bisquare (moins sensible aux outliers)")
        self.kernel_button_group.addButton(self.bisquare_radio, 2)
        kernel_layout.addWidget(self.bisquare_radio)

        self.tricube_radio = QRadioButton("Tricube")
        self.kernel_button_group.addButton(self.tricube_radio, 3)
        kernel_layout.addWidget(self.tricube_radio)

        self.exponential_radio = QRadioButton("Exponentiel")
        self.kernel_button_group.addButton(self.exponential_radio, 4)
        kernel_layout.addWidget(self.exponential_radio)

        self.boxcar_radio = QRadioButton("Boxcar")
        self.kernel_button_group.addButton(self.boxcar_radio, 5)
        kernel_layout.addWidget(self.boxcar_radio)

        kernel_group.setLayout(kernel_layout)
        layout.addWidget(kernel_group)

        # Bande passante
        bandwidth_group = QGroupBox("Bande passante (bandwidth)")
        bandwidth_layout = QVBoxLayout()
        self.bandwidth_button_group = QButtonGroup()

        self.auto_bandwidth_radio = QRadioButton("Automatique (validation croisée CV)")
        self.auto_bandwidth_radio.setChecked(True)
        self.bandwidth_button_group.addButton(self.auto_bandwidth_radio, 1)
        bandwidth_layout.addWidget(self.auto_bandwidth_radio)

        self.aic_bandwidth_radio = QRadioButton("Automatique (AIC)")
        self.bandwidth_button_group.addButton(self.aic_bandwidth_radio, 2)
        bandwidth_layout.addWidget(self.aic_bandwidth_radio)

        self.manual_bandwidth_radio = QRadioButton("Manuelle")
        self.bandwidth_button_group.addButton(self.manual_bandwidth_radio, 3)
        bandwidth_layout.addWidget(self.manual_bandwidth_radio)

        manual_layout = QHBoxLayout()
        manual_layout.addWidget(QLabel("  Valeur:"))
        self.bandwidth_spin = QDoubleSpinBox()
        self.bandwidth_spin.setMinimum(0.1)
        self.bandwidth_spin.setMaximum(1000000)
        self.bandwidth_spin.setValue(1000)
        self.bandwidth_spin.setDecimals(2)
        manual_layout.addWidget(self.bandwidth_spin)
        manual_layout.addStretch()
        bandwidth_layout.addLayout(manual_layout)

        bandwidth_group.setLayout(bandwidth_layout)
        layout.addWidget(bandwidth_group)

        # Options
        options_group = QGroupBox("Options")
        options_layout = QVBoxLayout()

        self.adaptive_check = QCheckBox("Bande passante adaptative (nombre de voisins si coché)")
        options_layout.addWidget(self.adaptive_check)

        adaptive_layout = QHBoxLayout()
        adaptive_layout.addWidget(QLabel("  Nombre de voisins si adaptative et manuelle:"))
        self.neighbors_spin = QSpinBox()
        self.neighbors_spin.setMinimum(10)
        self.neighbors_spin.setMaximum(1000)
        self.neighbors_spin.setValue(50)
        adaptive_layout.addWidget(self.neighbors_spin)
        adaptive_layout.addStretch()
        options_layout.addLayout(adaptive_layout)

        self.standardize_check = QCheckBox("Standardiser les variables")
        self.standardize_check.setChecked(True)
        options_layout.addWidget(self.standardize_check)

        self.robust_check = QCheckBox("GWR robuste (plus résistant aux outliers)")
        self.robust_check.setChecked(False)
        options_layout.addWidget(self.robust_check)

        options_group.setLayout(options_layout)
        layout.addWidget(options_group)

        # Boutons
        button_layout = QHBoxLayout()
        self.test_r_button = QPushButton("Tester R")
        self.test_r_button.clicked.connect(self.test_r_installation)
        button_layout.addWidget(self.test_r_button)
        button_layout.addStretch()

        self.ok_button = QPushButton("Lancer l'analyse")
        self.ok_button.clicked.connect(self.accept)
        button_layout.addWidget(self.ok_button)

        layout.addLayout(button_layout)
        self.setLayout(layout)

    def load_r_path(self):
        if os.path.exists(CONFIG_FILE):
            config = configparser.ConfigParser()
            config.read(CONFIG_FILE)
            if "R" in config and "path" in config["R"]:
                return config["R"]["path"]
        return None

    def save_r_path(self, path):
        config = configparser.ConfigParser()
        config["R"] = {"path": path}
        with open(CONFIG_FILE, "w", encoding="utf-8") as f:
            config.write(f)

    def choose_r_path(self):
        path, _ = QFileDialog.getOpenFileName(self, "Sélectionnez Rscript", "", "")
        if path:
            self.r_path = path
            self.r_path_label.setText(path)
            self.save_r_path(path)

    def test_r_installation(self):
        r_exe = self.r_path or "Rscript"
        try:
            result = subprocess.run(
                [r_exe, '--version'],
                capture_output=True,
                text=True,
                timeout=20,
                shell=False
            )

            if result.returncode == 0:
                version = result.stdout.strip() or result.stderr.strip() or "(aucune sortie)"
                
                test_script = """
packages <- c('GWmodel', 'sp', 'sf')
installed <- packages %in% rownames(installed.packages())
cat(paste(packages[!installed], collapse=', '))
"""
                result = subprocess.run(
                    [r_exe, '-e', test_script],
                    capture_output=True,
                    text=True,
                    timeout=20,
                    shell=False
                )

                missing = result.stdout.strip()

                if missing:
                    QMessageBox.warning(self, "Packages manquants",
                        f"R est installé ({version}) mais les packages suivants sont manquants:\n"
                        f"{missing}\n\nInstallez-les avec:\n"
                        f"install.packages(c('{missing}'))")
                else:
                    QMessageBox.information(self, "Test réussi",
                        f"R est installé ({version}) et tous les packages nécessaires sont présents!")
            else:
                QMessageBox.critical(self, "Erreur",
                    f"Échec de l'exécution de Rscript ({r_exe}).\nCode retour: {result.returncode}")

        except subprocess.TimeoutExpired:
            QMessageBox.critical(self, "Erreur", f"Le test de R a pris trop de temps (>20s).")
        except FileNotFoundError:
            QMessageBox.critical(self, "Erreur",
                f"R n'est pas installé ou Rscript n'est pas accessible.\n\nChemin testé : {r_exe}")
        except Exception as e:
            QMessageBox.critical(self, "Erreur", f"Erreur lors du test: {str(e)}")

    def populate_layers(self):
        self.layer_combo.clear()
        for layer in QgsProject.instance().mapLayers().values():
            if isinstance(layer, QgsVectorLayer):
                self.layer_combo.addItem(layer.name(), layer)

    def populate_fields(self):
        self.dependent_combo.clear()
        layer = self.layer_combo.currentData()
        if layer:
            for field in layer.fields():
                if field.type() in [QVariant.Int, QVariant.Double, QVariant.LongLong]:
                    self.dependent_combo.addItem(field.name())

    def populate_independent_vars(self):
        self.independent_list.clear()
        layer = self.layer_combo.currentData()
        dependent_var = self.dependent_combo.currentText()
        if layer:
            for field in layer.fields():
                if field.type() in [QVariant.Int, QVariant.Double, QVariant.LongLong]:
                    if field.name() != dependent_var:
                        self.independent_list.addItem(field.name())

    def get_selected_layer(self):
        return self.layer_combo.currentData()

    def get_dependent_var(self):
        return self.dependent_combo.currentText()

    def get_independent_vars(self):
        return [i.text() for i in self.independent_list.selectedItems()]

    def get_kernel_type(self):
        kernels = {1: "gaussian", 2: "bisquare", 3: "tricube", 4: "exponential", 5: "boxcar"}
        return kernels.get(self.kernel_button_group.checkedId(), "gaussian")

    def get_bandwidth_type(self):
        return self.bandwidth_button_group.checkedId()

    def get_bandwidth_value(self):
        return self.bandwidth_spin.value()

    def get_adaptive(self):
        return self.adaptive_check.isChecked()

    def get_neighbors(self):
        return self.neighbors_spin.value()

    def get_standardize(self):
        return self.standardize_check.isChecked()

    def get_robust(self):
        return self.robust_check.isChecked()

# =====================================================
# Dialogue MGWR
# =====================================================

class MGWRAnalysisDialog(QDialog):
    """Interface pour configurer l'analyse MGWR"""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Régression Multi-échelle (MGWR) avec R - GWmodel")
        self.setMinimumWidth(500)
        self.setMinimumHeight(500)  # Réduit car moins d'options
        self.r_path = self.load_r_path()

        layout = QVBoxLayout()

        # Chemin Rscript
        layout.addWidget(QLabel("Chemin vers Rscript :"))
        path_layout = QHBoxLayout()
        self.r_path_label = QLabel(self.r_path if self.r_path else "(non défini)")
        self.choose_r_path_button = QPushButton("Sélectionner Rscript")
        self.choose_r_path_button.clicked.connect(self.choose_r_path)
        path_layout.addWidget(self.r_path_label)
        path_layout.addWidget(self.choose_r_path_button)
        layout.addLayout(path_layout)
        layout.addSpacing(10)

        # Sélection couche
        layout.addWidget(QLabel("Couche vectorielle:"))
        self.layer_combo = QComboBox()
        self.populate_layers()
        layout.addWidget(self.layer_combo)

        # Variable dépendante
        layout.addWidget(QLabel("Variable dépendante (Y):"))
        self.dependent_combo = QComboBox()
        self.layer_combo.currentIndexChanged.connect(self.populate_fields)
        self.populate_fields()
        layout.addWidget(self.dependent_combo)

        # Variables indépendantes
        layout.addWidget(QLabel("Variables indépendantes (X) - sélection multiple:"))
        self.independent_list = QListWidget()
        self.independent_list.setSelectionMode(QAbstractItemView.MultiSelection)
        self.dependent_combo.currentIndexChanged.connect(self.populate_independent_vars)
        self.populate_independent_vars()
        layout.addWidget(self.independent_list)

        # Kernel
        kernel_group = QGroupBox("Fonction de pondération (kernel)")
        kernel_layout = QVBoxLayout()
        self.kernel_button_group = QButtonGroup()

        self.gaussian_radio = QRadioButton("Gaussien (recommandé)")
        self.gaussian_radio.setChecked(True)
        self.kernel_button_group.addButton(self.gaussian_radio, 1)
        kernel_layout.addWidget(self.gaussian_radio)

        self.bisquare_radio = QRadioButton("Bisquare")
        self.kernel_button_group.addButton(self.bisquare_radio, 2)
        kernel_layout.addWidget(self.bisquare_radio)

        self.tricube_radio = QRadioButton("Tricube")
        self.kernel_button_group.addButton(self.tricube_radio, 3)
        kernel_layout.addWidget(self.tricube_radio)

        self.exponential_radio = QRadioButton("Exponentiel")
        self.kernel_button_group.addButton(self.exponential_radio, 4)
        kernel_layout.addWidget(self.exponential_radio)

        self.boxcar_radio = QRadioButton("Boxcar")
        self.kernel_button_group.addButton(self.boxcar_radio, 5)
        kernel_layout.addWidget(self.boxcar_radio)

        kernel_group.setLayout(kernel_layout)
        layout.addWidget(kernel_group)

        # Options MGWR
        options_group = QGroupBox("Options MGWR")
        options_layout = QVBoxLayout()

        self.adaptive_check = QCheckBox("Bande passante adaptative (recommandé pour MGWR)")
        self.adaptive_check.setChecked(True)
        options_layout.addWidget(self.adaptive_check)

        # Critère d'optimisation
        criterion_layout = QHBoxLayout()
        criterion_layout.addWidget(QLabel("Critère d'optimisation:"))
        self.criterion_combo = QComboBox()
        self.criterion_combo.addItem("AICc (recommandé)", "AICc")
        self.criterion_combo.addItem("AIC", "AIC")
        self.criterion_combo.addItem("CV (validation croisée)", "CV")
        criterion_layout.addWidget(self.criterion_combo)
        criterion_layout.addStretch()
        options_layout.addLayout(criterion_layout)

        self.standardize_check = QCheckBox("Standardiser les variables")
        self.standardize_check.setChecked(True)
        options_layout.addWidget(self.standardize_check)

        options_group.setLayout(options_layout)
        layout.addWidget(options_group)

        # Paramètres avancés
        advanced_group = QGroupBox("Paramètres avancés")
        advanced_layout = QVBoxLayout()

        # Nombre d'itérations
        iter_layout = QHBoxLayout()
        iter_layout.addWidget(QLabel("Nombre max d'itérations:"))
        self.max_iter_spin = QSpinBox()
        self.max_iter_spin.setMinimum(10)
        self.max_iter_spin.setMaximum(500)
        self.max_iter_spin.setValue(200)
        iter_layout.addWidget(self.max_iter_spin)
        iter_layout.addStretch()
        advanced_layout.addLayout(iter_layout)

        # Tolérance
        tol_layout = QHBoxLayout()
        tol_layout.addWidget(QLabel("Tolérance de convergence:"))
        self.tolerance_spin = QDoubleSpinBox()
        self.tolerance_spin.setMinimum(0.00001)
        self.tolerance_spin.setMaximum(0.1)
        self.tolerance_spin.setValue(0.01)
        self.tolerance_spin.setDecimals(5)
        tol_layout.addWidget(self.tolerance_spin)
        tol_layout.addStretch()
        advanced_layout.addLayout(tol_layout)

        advanced_group.setLayout(advanced_layout)
        layout.addWidget(advanced_group)

        # Note informative
        info_label = QLabel("Note : MGWR calcule automatiquement une bande passante optimale pour chaque variable indépendante.")
        info_label.setWordWrap(True)
        info_label.setStyleSheet("color: gray; font-style: italic; padding: 5px;")
        layout.addWidget(info_label)

        # Boutons
        button_layout = QHBoxLayout()
        self.test_r_button = QPushButton("Tester R et GWmodel")
        self.test_r_button.clicked.connect(self.test_r_installation)
        button_layout.addWidget(self.test_r_button)
        button_layout.addStretch()

        self.ok_button = QPushButton("Lancer l'analyse MGWR")
        self.ok_button.clicked.connect(self.accept)
        button_layout.addWidget(self.ok_button)

        layout.addLayout(button_layout)
        self.setLayout(layout)

    def load_r_path(self):
        if os.path.exists(CONFIG_FILE):
            config = configparser.ConfigParser()
            config.read(CONFIG_FILE)
            if "R" in config and "path" in config["R"]:
                return config["R"]["path"]
        return None

    def save_r_path(self, path):
        config = configparser.ConfigParser()
        config["R"] = {"path": path}
        with open(CONFIG_FILE, "w", encoding="utf-8") as f:
            config.write(f)

    def choose_r_path(self):
        path, _ = QFileDialog.getOpenFileName(self, "Sélectionnez Rscript", "", "")
        if path:
            self.r_path = path
            self.r_path_label.setText(path)
            self.save_r_path(path)

    def test_r_installation(self):
        r_exe = self.r_path or "Rscript"
        try:
            result = subprocess.run(
                [r_exe, '--version'],
                capture_output=True,
                text=True,
                timeout=20,
                shell=False
            )

            if result.returncode == 0:
                version = result.stdout.strip() or result.stderr.strip() or "(aucune sortie)"
                
                test_script = """
packages <- c('GWmodel', 'sp', 'sf')
installed <- packages %in% rownames(installed.packages())
cat(paste(packages[!installed], collapse=', '))
"""
                result = subprocess.run(
                    [r_exe, '-e', test_script],
                    capture_output=True,
                    text=True,
                    timeout=20,
                    shell=False
                )

                missing = result.stdout.strip()

                if missing:
                    QMessageBox.warning(self, "Packages manquants",
                        f"R est installé ({version}) mais les packages suivants sont manquants:\n"
                        f"{missing}\n\nInstallez-les avec:\n"
                        f"install.packages(c('{missing}'))")
                else:
                    QMessageBox.information(self, "Test réussi",
                        f"R est installé ({version}) et tous les packages nécessaires sont présents!\n\n"
                        f"Note : GWR et MGWR utilisent tous deux le package GWmodel.")
            else:
                QMessageBox.critical(self, "Erreur",
                    f"Échec de l'exécution de Rscript ({r_exe}).\nCode retour: {result.returncode}")

        except subprocess.TimeoutExpired:
            QMessageBox.critical(self, "Erreur", f"Le test de R a pris trop de temps (>20s).")
        except FileNotFoundError:
            QMessageBox.critical(self, "Erreur",
                f"R n'est pas installé ou Rscript n'est pas accessible.\n\nChemin testé : {r_exe}")
        except Exception as e:
            QMessageBox.critical(self, "Erreur", f"Erreur lors du test: {str(e)}")

    def populate_layers(self):
        self.layer_combo.clear()
        for layer in QgsProject.instance().mapLayers().values():
            if isinstance(layer, QgsVectorLayer):
                self.layer_combo.addItem(layer.name(), layer)

    def populate_fields(self):
        self.dependent_combo.clear()
        layer = self.layer_combo.currentData()
        if layer:
            for field in layer.fields():
                if field.type() in [QVariant.Int, QVariant.Double, QVariant.LongLong]:
                    self.dependent_combo.addItem(field.name())

    def populate_independent_vars(self):
        self.independent_list.clear()
        layer = self.layer_combo.currentData()
        dependent_var = self.dependent_combo.currentText()
        if layer:
            for field in layer.fields():
                if field.type() in [QVariant.Int, QVariant.Double, QVariant.LongLong]:
                    if field.name() != dependent_var:
                        self.independent_list.addItem(field.name())

    def get_selected_layer(self):
        return self.layer_combo.currentData()

    def get_dependent_var(self):
        return self.dependent_combo.currentText()

    def get_independent_vars(self):
        return [i.text() for i in self.independent_list.selectedItems()]

    def get_kernel_type(self):
        kernels = {1: "gaussian", 2: "bisquare", 3: "tricube", 4: "exponential", 5: "boxcar"}
        return kernels.get(self.kernel_button_group.checkedId(), "gaussian")

    def get_adaptive(self):
        return self.adaptive_check.isChecked()

    def get_standardize(self):
        return self.standardize_check.isChecked()

    def get_criterion(self):
        return self.criterion_combo.currentData()

    def get_max_iter(self):
        return self.max_iter_spin.value()

    def get_tolerance(self):
        return self.tolerance_spin.value()

# =====================================================
# Dialogue LISA
# =====================================================

class LISAAnalysisDialog(QDialog):
    """Interface pour configurer l'analyse LISA"""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Analyse LISA (Local Indicators of Spatial Association)")
        self.setMinimumWidth(500)
        self.setMinimumHeight(500)  # Réduit de 600 à 500
        self.r_path = self.load_r_path()

        layout = QVBoxLayout()

        # Chemin Rscript
        layout.addWidget(QLabel("Chemin vers Rscript :"))
        path_layout = QHBoxLayout()
        self.r_path_label = QLabel(self.r_path if self.r_path else "(non défini)")
        self.choose_r_path_button = QPushButton("Sélectionner Rscript")
        self.choose_r_path_button.clicked.connect(self.choose_r_path)
        path_layout.addWidget(self.r_path_label)
        path_layout.addWidget(self.choose_r_path_button)
        layout.addLayout(path_layout)
        layout.addSpacing(10)

        # Sélection couche
        layout.addWidget(QLabel("Couche vectorielle:"))
        self.layer_combo = QComboBox()
        self.populate_layers()
        layout.addWidget(self.layer_combo)

        # Type d'analyse LISA
        type_group = QGroupBox("Type d'analyse LISA")
        type_layout = QVBoxLayout()
        self.type_button_group = QButtonGroup()

        self.univariate_radio = QRadioButton("Univariée (une seule variable)")
        self.univariate_radio.setChecked(True)
        self.type_button_group.addButton(self.univariate_radio, 1)
        type_layout.addWidget(self.univariate_radio)

        self.bivariate_radio = QRadioButton("Bivariée (deux variables)")
        self.type_button_group.addButton(self.bivariate_radio, 2)
        type_layout.addWidget(self.bivariate_radio)

        type_group.setLayout(type_layout)
        layout.addWidget(type_group)

        # Variable principale
        layout.addWidget(QLabel("Variable d'analyse:"))
        self.variable_combo = QComboBox()
        self.layer_combo.currentIndexChanged.connect(self.populate_fields)
        self.populate_fields()
        layout.addWidget(self.variable_combo)

        # Variable secondaire (pour analyse bivariée)
        layout.addWidget(QLabel("Variable secondaire (uniquement pour analyse bivariée):"))
        self.variable2_combo = QComboBox()
        self.variable_combo.currentIndexChanged.connect(self.populate_variable2)
        self.populate_variable2()
        layout.addWidget(self.variable2_combo)
        self.variable2_combo.setEnabled(False)
        self.bivariate_radio.toggled.connect(lambda checked: self.variable2_combo.setEnabled(checked))

        # Matrice de poids spatiaux
        weights_group = QGroupBox("Matrice de poids spatiaux")
        weights_layout = QVBoxLayout()

        # Type de contiguïté
        cont_layout = QHBoxLayout()
        cont_layout.addWidget(QLabel("Type de contiguïté:"))
        self.contiguity_combo = QComboBox()
        self.contiguity_combo.addItem("Reine (Queen)", "queen")
        self.contiguity_combo.addItem("Tour (Rook)", "rook")
        cont_layout.addWidget(self.contiguity_combo)
        cont_layout.addStretch()
        weights_layout.addLayout(cont_layout)

        # Ordre de contiguïté
        order_layout = QHBoxLayout()
        order_layout.addWidget(QLabel("Ordre de contiguïté:"))
        self.order_spin = QSpinBox()
        self.order_spin.setMinimum(1)
        self.order_spin.setMaximum(5)
        self.order_spin.setValue(1)
        order_layout.addWidget(self.order_spin)
        order_layout.addStretch()
        weights_layout.addLayout(order_layout)

        # Standardisation des poids
        self.standardize_weights_check = QCheckBox("Standardiser les poids (row-standardized)")
        self.standardize_weights_check.setChecked(True)
        weights_layout.addWidget(self.standardize_weights_check)

        weights_group.setLayout(weights_layout)
        layout.addWidget(weights_group)

        # Options (seuil de significativité uniquement)
        options_group = QGroupBox("Options")
        options_layout = QVBoxLayout()

        # Seuil de significativité
        sig_layout = QHBoxLayout()
        sig_layout.addWidget(QLabel("Seuil de significativité:"))
        self.significance_combo = QComboBox()
        self.significance_combo.addItem("0.01 (1%)", 0.01)
        self.significance_combo.addItem("0.05 (5%)", 0.05)
        self.significance_combo.addItem("0.10 (10%)", 0.10)
        self.significance_combo.setCurrentIndex(1)  # 0.05 par défaut
        sig_layout.addWidget(self.significance_combo)
        sig_layout.addStretch()
        options_layout.addLayout(sig_layout)

        self.standardize_var_check = QCheckBox("Standardiser la variable")
        self.standardize_var_check.setChecked(False)
        options_layout.addWidget(self.standardize_var_check)

        options_group.setLayout(options_layout)
        layout.addWidget(options_group)

        # Note informative
        info_label = QLabel("Note : Les p-values sont calculées automatiquement par rgeoda avec la fonction lisa_pvalues()")
        info_label.setWordWrap(True)
        info_label.setStyleSheet("color: gray; font-style: italic; padding: 5px;")
        layout.addWidget(info_label)

        # Boutons
        button_layout = QHBoxLayout()
        self.test_r_button = QPushButton("Tester R et rgeoda")
        self.test_r_button.clicked.connect(self.test_r_installation)
        button_layout.addWidget(self.test_r_button)
        button_layout.addStretch()

        self.ok_button = QPushButton("Lancer l'analyse LISA")
        self.ok_button.clicked.connect(self.accept)
        button_layout.addWidget(self.ok_button)

        layout.addLayout(button_layout)
        self.setLayout(layout)

    def load_r_path(self):
        if os.path.exists(CONFIG_FILE):
            config = configparser.ConfigParser()
            config.read(CONFIG_FILE)
            if "R" in config and "path" in config["R"]:
                return config["R"]["path"]
        return None

    def save_r_path(self, path):
        config = configparser.ConfigParser()
        config["R"] = {"path": path}
        with open(CONFIG_FILE, "w", encoding="utf-8") as f:
            config.write(f)

    def choose_r_path(self):
        path, _ = QFileDialog.getOpenFileName(self, "Sélectionnez Rscript", "", "")
        if path:
            self.r_path = path
            self.r_path_label.setText(path)
            self.save_r_path(path)

    def test_r_installation(self):
        r_exe = self.r_path or "Rscript"
        try:
            result = subprocess.run(
                [r_exe, '--version'],
                capture_output=True,
                text=True,
                timeout=20,
                shell=False
            )

            if result.returncode == 0:
                version = result.stdout.strip() or result.stderr.strip() or "(aucune sortie)"
                
                test_script = """
packages <- c('rgeoda', 'sf')
installed <- packages %in% rownames(installed.packages())
cat(paste(packages[!installed], collapse=', '))
"""
                result = subprocess.run(
                    [r_exe, '-e', test_script],
                    capture_output=True,
                    text=True,
                    timeout=20,
                    shell=False
                )

                missing = result.stdout.strip()

                if missing:
                    QMessageBox.warning(self, "Packages manquants",
                        f"R est installé ({version}) mais les packages suivants sont manquants:\n"
                        f"{missing}\n\nInstallez-les avec:\n"
                        f"install.packages(c('{missing}'))")
                else:
                    QMessageBox.information(self, "Test réussi",
                        f"R est installé ({version}) et tous les packages nécessaires sont présents!")
            else:
                QMessageBox.critical(self, "Erreur",
                    f"Échec de l'exécution de Rscript ({r_exe}).\nCode retour: {result.returncode}")

        except subprocess.TimeoutExpired:
            QMessageBox.critical(self, "Erreur", f"Le test de R a pris trop de temps (>20s).")
        except FileNotFoundError:
            QMessageBox.critical(self, "Erreur",
                f"R n'est pas installé ou Rscript n'est pas accessible.\n\nChemin testé : {r_exe}")
        except Exception as e:
            QMessageBox.critical(self, "Erreur", f"Erreur lors du test: {str(e)}")

    def populate_layers(self):
        self.layer_combo.clear()
        for layer in QgsProject.instance().mapLayers().values():
            if isinstance(layer, QgsVectorLayer):
                self.layer_combo.addItem(layer.name(), layer)

    def populate_fields(self):
        self.variable_combo.clear()
        layer = self.layer_combo.currentData()
        if layer:
            for field in layer.fields():
                if field.type() in [QVariant.Int, QVariant.Double, QVariant.LongLong]:
                    self.variable_combo.addItem(field.name())

    def populate_variable2(self):
        self.variable2_combo.clear()
        layer = self.layer_combo.currentData()
        variable1 = self.variable_combo.currentText()
        if layer:
            for field in layer.fields():
                if field.type() in [QVariant.Int, QVariant.Double, QVariant.LongLong]:
                    if field.name() != variable1:
                        self.variable2_combo.addItem(field.name())

    def get_selected_layer(self):
        return self.layer_combo.currentData()

    def get_analysis_type(self):
        return "univariate" if self.type_button_group.checkedId() == 1 else "bivariate"

    def get_variable(self):
        return self.variable_combo.currentText()

    def get_variable2(self):
        return self.variable2_combo.currentText() if self.bivariate_radio.isChecked() else None

    def get_contiguity_type(self):
        return self.contiguity_combo.currentData()

    def get_order(self):
        return self.order_spin.value()

    def get_standardize_weights(self):
        return self.standardize_weights_check.isChecked()

    def get_significance(self):
        return self.significance_combo.currentData()

    def get_standardize_variable(self):
        return self.standardize_var_check.isChecked()

# =====================================================
# Plugin principal QGIS
# =====================================================

class GWRPlugin:
    """Point d'entrée du plugin pour QGIS"""

    def __init__(self, iface):
        self.iface = iface
        self.action = None

    def initGui(self):
        self.action = QAction("Analyse GWR / MGWR / LISA", self.iface.mainWindow())
        self.action.triggered.connect(self.run)
        self.iface.addPluginToMenu("&Analyse Spatiale", self.action)
        self.iface.addToolBarIcon(self.action)

    def unload(self):
        self.iface.removePluginMenu("&Analyse Spatiale", self.action)
        self.iface.removeToolBarIcon(self.action)

    def run(self):
        """Fonction principale - affiche d'abord le choix du modèle"""
        
        # Étape 1 : Choix du modèle
        model_dialog = ModelSelectionDialog()
        if not model_dialog.exec_():
            return
        
        selected_model = model_dialog.selected_model
        
        if selected_model == "GWR":
            self.run_gwr_analysis()
        elif selected_model == "MGWR":
            self.run_mgwr_analysis()
        elif selected_model == "LISA":
            self.run_lisa_analysis()
    
    def run_gwr_analysis(self):
        """Lance l'analyse GWR - appelle gwr_analysis.py"""
        dialog = GWRAnalysisDialog()
        
        if dialog.exec_():
            layer = dialog.get_selected_layer()
            dependent_var = dialog.get_dependent_var()
            independent_vars = dialog.get_independent_vars()
            kernel_type = dialog.get_kernel_type()
            bandwidth_type = dialog.get_bandwidth_type()
            bandwidth_value = dialog.get_bandwidth_value()
            adaptive = dialog.get_adaptive()
            neighbors = dialog.get_neighbors()
            standardize = dialog.get_standardize()
            robust = dialog.get_robust()
            r_path = dialog.r_path
            
            if not layer or not dependent_var or not independent_vars:
                QMessageBox.warning(None, "Erreur",
                    "Veuillez sélectionner une couche, une variable dépendante et au moins une variable indépendante")
                return
            
            if dialog.r_path:
                os.environ["PATH"] += ";" + os.path.dirname(dialog.r_path)
            
            progress = QProgressDialog("Analyse GWR en cours...", "Annuler", 0, 0)
            progress.setWindowTitle("Calcul GWR")
            progress.setModal(True)
            progress.show()

            # Appel du module GWR séparé
            result_layer, message, temp_folder = GWRAnalysisModule.run_analysis(
                r_path, layer, dependent_var, independent_vars, kernel_type,
                bandwidth_type, bandwidth_value, adaptive, neighbors, standardize, robust
            )

            progress.close()

            if result_layer:
                self.process_results(result_layer, layer, message, "GWR", temp_folder)
            else:
                QMessageBox.critical(None, "Erreur GWR", message)
    
    def run_mgwr_analysis(self):
        """Lance l'analyse MGWR - appelle mgwr_analysis.py"""
        dialog = MGWRAnalysisDialog()
        
        if dialog.exec_():
            layer = dialog.get_selected_layer()
            dependent_var = dialog.get_dependent_var()
            independent_vars = dialog.get_independent_vars()
            kernel_type = dialog.get_kernel_type()
            adaptive = dialog.get_adaptive()
            standardize = dialog.get_standardize()
            criterion = dialog.get_criterion()
            max_iter = dialog.get_max_iter()
            tolerance = dialog.get_tolerance()
            r_path = dialog.r_path
            
            if not layer or not dependent_var or not independent_vars:
                QMessageBox.warning(None, "Erreur",
                    "Veuillez sélectionner une couche, une variable dépendante et au moins une variable indépendante")
                return
            
            if dialog.r_path:
                os.environ["PATH"] += ";" + os.path.dirname(dialog.r_path)
            
            progress = QProgressDialog("Analyse MGWR en cours...", "Annuler", 0, 0)
            progress.setWindowTitle("Calcul MGWR")
            progress.setModal(True)
            progress.show()

            # Appel du module MGWR séparé
            result_layer, message, temp_folder = MGWRAnalysisModule.run_analysis(
                r_path, layer, dependent_var, independent_vars, kernel_type,
                adaptive, standardize, criterion, max_iter, tolerance
            )

            progress.close()

            if result_layer:
                self.process_results(result_layer, layer, message, "MGWR")
            else:
                QMessageBox.critical(None, "Erreur MGWR", message)
    
    def run_lisa_analysis(self):
        """Lance l'analyse LISA - appelle lisa_analysis.py"""
        dialog = LISAAnalysisDialog()
        
        if dialog.exec_():
            layer = dialog.get_selected_layer()
            analysis_type = dialog.get_analysis_type()
            variable = dialog.get_variable()
            variable2 = dialog.get_variable2()
            contiguity_type = dialog.get_contiguity_type()
            order = dialog.get_order()
            standardize_weights = dialog.get_standardize_weights()
            significance = dialog.get_significance()
            standardize_var = dialog.get_standardize_variable()
            r_path = dialog.r_path
            
            if not layer or not variable:
                QMessageBox.warning(None, "Erreur",
                    "Veuillez sélectionner une couche et une variable d'analyse")
                return
            
            if analysis_type == "bivariate" and not variable2:
                QMessageBox.warning(None, "Erreur",
                    "Veuillez sélectionner une deuxième variable pour l'analyse bivariée")
                return
            
            if dialog.r_path:
                os.environ["PATH"] += ";" + os.path.dirname(dialog.r_path)
            
            progress = QProgressDialog("Analyse LISA en cours...", "Annuler", 0, 0)
            progress.setWindowTitle("Calcul LISA")
            progress.setModal(True)
            progress.show()

            try:
                # Appel du module LISA séparé - SANS permutations
                result_layer, message, temp_folder = LISAAnalysisModule.run_analysis(
                    r_path, layer, analysis_type, variable, variable2,
                    contiguity_type, order, standardize_weights,
                    significance, standardize_var
                )

                progress.close()

                if result_layer:
                    self.process_results(result_layer, layer, message, "LISA")
                else:
                    QMessageBox.critical(None, "Erreur LISA", message)
            
            except Exception as e:
                progress.close()
                import traceback
                error_msg = f"Erreur lors de l'analyse LISA:\n{str(e)}\n\n{traceback.format_exc()}"
                QMessageBox.critical(None, "Erreur LISA", error_msg)
                print(error_msg)
    
    def process_results(self, result_layer, original_layer, message, model_type, temp_folder):
        """Traite les résultats (commun à GWR, MGWR et LISA) - VERSION SANS MODIFICATION DE LA COUCHE ORIGINALE"""
        
        # IMPORTANT: Cloner la couche en mémoire pour libérer les fichiers temporaires
        # Nom de la nouvelle couche
        layer_name = f"{original_layer.name()}_{model_type}"
        
        # Créer une URI pour une couche en mémoire
        geom_type = result_layer.geometryType()
        geom_type_str = {0: "Point", 1: "LineString", 2: "Polygon", 3: "Unknown", 4: "NoGeometry"}
        geom_str = geom_type_str.get(geom_type, "Point")
        
        crs = result_layer.crs().authid()
        uri = f"{geom_str}?crs={crs}"
        
        # Créer la couche en mémoire
        memory_layer = QgsVectorLayer(uri, layer_name, "memory")
        memory_provider = memory_layer.dataProvider()
        
        # Copier les champs
        memory_provider.addAttributes(result_layer.fields())
        memory_layer.updateFields()
        
        # Copier les entités
        features = list(result_layer.getFeatures())
        memory_provider.addFeatures(features)
        
        # Ajouter la couche en mémoire au projet
        QgsProject.instance().addMapLayer(memory_layer)
        
        # Maintenant on peut supprimer les fichiers temporaires en toute sécurité
        try:
            # Forcer la fermeture de la couche temporaire
            result_layer = None
            
            # Supprimer le dossier temporaire
            if temp_folder and os.path.exists(temp_folder):
                print(f"Suppression du dossier temporaire: {temp_folder}")
                shutil.rmtree(temp_folder, ignore_errors=True)
        except Exception as e:
            print(f"Avertissement: Impossible de supprimer le dossier temporaire: {e}")
            print("Le dossier sera supprimé automatiquement par le système")
        
        # Message de succès simple
        success_msg = (
            f"Analyse {model_type} terminée avec succès!\n\n"
            f"Une nouvelle couche '{memory_layer.name()}' a été créée avec les résultats.\n"
            f"La couche originale '{original_layer.name()}' n'a pas été modifiée."
        )
        
        # Vérifier si le nombre d'entités diffère
        if memory_layer.featureCount() != original_layer.featureCount():
            success_msg += (
                f"\n\n⚠ Attention: Le nombre d'entités diffère:\n"
                f"• Couche originale: {original_layer.featureCount()}\n"
                f"• Résultats {model_type}: {memory_layer.featureCount()}"
            )
        
        QMessageBox.information(None, "Succès", success_msg)

  
        # Afficher des infos dans la console
        print(f"\n=== RÉSULTATS {model_type} ===")
        print(f"✓ Nouvelle couche créée: {memory_layer.name()}")
        print(f"✓ Nombre d'entités: {memory_layer.featureCount()}")
        print(f"✓ Nombre de champs: {len(memory_layer.fields())}")
        print(f"✓ Couche originale '{original_layer.name()}' préservée")
        
        # Lister les champs de résultats
        prefix = f"{model_type}_"
        result_fields = [f.name() for f in memory_layer.fields() if f.name().startswith(prefix)]
        if result_fields:
            print(f"✓ Champs de résultats ajoutés ({len(result_fields)}):")
            for field_name in result_fields[:10]:  # Afficher max 10
                print(f"  - {field_name}")
            if len(result_fields) > 10:
                print(f"  ... et {len(result_fields) - 10} autres champs")
        
        print(f"{'='*50}\n")
        
        # Afficher les résultats détaillés
        from qgis.PyQt.QtWidgets import QTextEdit, QVBoxLayout, QDialog, QPushButton
        from qgis.PyQt.QtGui import QFont
        
        result_dialog = QDialog()
        result_dialog.setWindowTitle(f"Résultats de l'analyse {model_type}")
        result_dialog.setMinimumWidth(700)
        result_dialog.setMinimumHeight(500)
        
        layout = QVBoxLayout()
        
        text_edit = QTextEdit()
        text_edit.setReadOnly(True)
        text_edit.setPlainText(message)
        
        font = QFont("Courier New", 9)
        text_edit.setFont(font)
        
        text_edit.setStyleSheet("""
            QTextEdit {
                background-color: #2b2b2b;
                color: #f0f0f0;
                border: 1px solid #555;
                padding: 10px;
            }
        """)
        
        layout.addWidget(text_edit)
        
        ok_button = QPushButton("OK")
        ok_button.setMinimumHeight(30)
        ok_button.clicked.connect(result_dialog.accept)
        layout.addWidget(ok_button)
        
        result_dialog.setLayout(layout)
        result_dialog.exec_()


def classFactory(iface):
    """Fonction requise par QGIS pour charger le plugin"""
    return GWRPlugin(iface)
