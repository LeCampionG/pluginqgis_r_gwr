# **Plugin Qgis pour réaliser des statistiques spatiales (GWR, MGWR et LISA) à l'aide de R**

## Données tests

**Des données test sont dans le dossier "data_test"**

1. Au format shape et geopackage
2. Des polygones et des points.


## Structure du Plugin


```
gwr_plugin/
│
├── __init__.py                 # Point d'entrée du plugin
├── gwr_plugin_r.py            # Fichier principal - Interface utilisateur
├── gwr_analysis.py            # Module d'analyse GWR
├── mgwr_analysis.py           # Module d'analyse MGWR
├── lisa_analysis.py           # Module d'analyse LISA
├── config.ini                 # Fichier de configuration (crée lorsque le chemin R est spécifié suite à la première utilisation)
└── metadata.txt               # Métadonnées du plugin
```

## Objectif du plugin
Réaliser une **Régression Géographiquement Pondérée (GWR)**, une **MGWR** ou une analyse de **LISA** directement depuis QGIS en utilisant R et les package `GWmodel` et `rgeoda`.

---

## Le processus


```

  1. UTILISATEUR                                             
     └─> Sélection du modèle (GWR/MGWR/LISA)                
     └─> Configuration des paramètres                        
     └─> Validation                                          

                          ↓

  2. PYTHON (gwr_plugin_r.py)                                
     └─> Récupération de la couche QGIS                     
     └─> Appel du module d'analyse approprié                

                          ↓

  3. MODULE D'ANALYSE (gwr/mgwr/lisa_analysis.py)           
     └─> Création d'un répertoire temporaire                
     └─> Mapping des noms de champs (limite 10 caractères)  
     └─> Export de la couche en geopackage                   
     └─> Génération du script R                             

                          ↓

  4. EXÉCUTION R                                 
     └─> Lecture du shapefile                              
     └─> Création d'une copie des données (sp_data_work)    
     └─> Transformations sur la copie uniquement            
     └─> Calculs statistiques                               
     └─> Création de result_sf à partir des données         
         originales non modifiées                            
     └─> Ajout des résultats à result_sf                    
     └─> Export du gpkg de résultats                   

                          ↓

  5. RETOUR À PYTHON                                         
     └─> Chargement du shapefile de résultats               
     └─> Création d'une QgsVectorLayer                      
     └─> Nettoyage des fichiers temporaires                 

                          ↓

  6. AFFICHAGE DANS QGIS                   
     └─> Ajout de la nouvelle couche au projet              
     └─> La couche originale reste intacte                  
     └─> Affichage du message de succès                     
     └─> Affichage des résultats détaillés                  

```

### Prérequis

1. **QGIS** (version 3.x)
2. **R** (version 4.0+)
3. **Packages R requis :**
   - Pour GWR/MGWR : `GWmodel`, `sp`, `sf`
   - Pour LISA : `rgeoda`, `sf`

### Installation des packages R

```r
# Dans R ou RStudio
install.packages("GWmodel")
install.packages("sp")
install.packages("sf")
install.packages("rgeoda")
```

### Installation du plugin

1. Copiez le dossier du plugin dans le répertoire des plugins QGIS :
   - Windows : `C:\Users\[utilisateur]\AppData\Roaming\QGIS\QGIS3\profiles\default\python\plugins\`
   - Linux : `~/.local/share/QGIS/QGIS3/profiles/default/python/plugins/`
   - Mac : `~/Library/Application Support/QGIS/QGIS3/profiles/default/python/plugins/`

**ou**

1bis. Ouvrir Qgis :
      - Menu `Extensions` → `Installer depuis un ZIP`

2. Redémarrez QGIS

3. Activez le plugin :
   - Menu `Extensions` → `Installer/Gérer les extensions`
   - Recherchez "GWR Analysis R"
   - Cochez la case pour l'activer

4. Configurez le chemin R :
   - Lancez le plugin
   - Cliquez sur "Sélectionner Rscript.exe"
   - Naviguez jusqu'à votre installation R (sur windows: `C:\Program Files\R\R-4.3.1\bin\Rscript.exe`)

## L'utilisation

### Workflow général

1. **Charger votre couche vectorielle dans QGIS**
   - Formats supportés : Shapefile, GeoPackage, etc.
   - Doit contenir des attributs numériques

2. **Lancer le plugin**
   - Icône `Analyse GWR / MGWR / LISA` dans la barre d'outils

3. **Choisir le type d'analyse**
   - GWR
   - MGWR
   - LISA

4. **Configurer les paramètres**
   - Sélectionner la couche
   - Choisir les variables
   - Ajuster les paramètres

5. **Lancer l'analyse**
   - Cliquer sur "Lancer l'analyse"
   - Patienter (peut prendre plusieurs minutes)

6. **Consulter les résultats**
   - Une nouvelle couche est créée automatiquement
   - Contient toutes les colonnes originales + les résultats
   - La couche d'origine reste intacte