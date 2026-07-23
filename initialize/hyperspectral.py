import logging
import os
import sys
# import whitebox
import zipfile
import ray

import pandas as pd
import numpy as np
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, ConfusionMatrixDisplay
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Patch
import rasterio
from rasterio.plot import show
from rasterio.features import geometry_mask
import fiona
from sklearn.preprocessing import OneHotEncoder
# from keras.models import Sequential
# from keras.layers import Conv2D, MaxPooling2D, Flatten, Dense
import geopandas as gpd
from shapely.geometry import Point
from shapely.geometry import mapping
from shapely.geometry import box
from collections import Counter
import rpy2.robjects as ro
from collections import Counter
import seaborn as sns

# data preparation
from sklearn.decomposition import PCA
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, f1_score

# model selection
from sklearn.model_selection import (
    GridSearchCV,
    GroupKFold,
    GroupShuffleSplit,
    KFold,
    ShuffleSplit,
    StratifiedGroupKFold,
    StratifiedKFold,
    StratifiedShuffleSplit,
    TimeSeriesSplit,
    train_test_split,
)

# model evaluation
from sklearn.metrics import f1_score

# machine learning
from sklearn.ensemble import RandomForestClassifier as RF

# Save files
import joblib
from joblib import dump, load

from utils.download_functions import download_aop_files
from utils.apply_brdf_corrections import move_red_yellow_subset, \
                                        convert_hdf_to_envi, \
                                        create_config, \
                                        implement_brdf_correction , \
                                        convert_envi_to_tif
from pathlib import Path

# add environ
conda_env_path = Path(sys.executable).parent.parent
os.environ['PROJ_LIB'] = str(conda_env_path/'share'/'proj')

log = logging.getLogger(__name__)

def download_hyperspectral(site, year, data_raw_aop_path, hs_type):
    """Download raw, uncorrected hyperspectral data (raw and raster) for site-year
    Analog to 04-download_aop_imagery.R from https://github.com/earthlab/neon-veg

    Parameters
    ----------
    site : str
        Site name
    year : str
        Year of aop collection
    data_raw_aop_path : str
        Path to store the downloaded data
    hs_type : str
        "tile" (L3 1km x 1km tile) OR "flightline" (L1 flightline)

    Returns
    -------
    (str, str)
        Path to the result folder for hyperspectral uncorrected imagery (L3 tiles) ('*/hs_tile' or '*/hs_flightline')
        and the result raster folder ('*/tif')
    """
    path = Path(data_raw_aop_path)/site/year

    product_code = 'DP3.30010.001' #rgb
    p = path/'tif'
    download_aop_files(product_code,
                        site,
                        year,
                        str(p),
                        match_string='323000_4097000',
                        check_size=False)

    product_code = 'DP3.30026.001' #veg indices
    download_aop_files(product_code,
                       site,
                       year,
                       str(p),
                       match_string='323000_4097000',
                       check_size=False)
    zip_files = [file for file in os.listdir(p) if file.endswith('.zip')] # get the list of files
    for zip_file in zip_files:  #for each zipfile
        with zipfile.ZipFile(Path(p/zip_file)) as item: # treat the file as a zip
            item.extractall(p)  # extract it in the working directory
            item.close()
    vi_error_files = [file for file in os.listdir(p) if file.endswith('_error.tif')]
    for vi_error_file in vi_error_files: 
        os.remove(Path(p/vi_error_file)) # remove error files

    if hs_type=="tile":
        product_code = 'DP3.30006.001'  #DP3.30006.002
        p = path/'hs_tile'
    elif hs_type=="flightline": 
        product_code = 'DP1.30006.001' #DP1.30006.002
        p = path/'hs_flightline'
    else:
        print("must specify hs_type argument")
    download_aop_files(product_code,
                       site,
                       year,
                       str(p),
                       match_string='323000_4097000',
                       check_size=False)
    
    return str(p)



def correct_flightlines(site, year_inv, year_aop, data_raw_aop_path, data_int_path):
    """Correct raw L1 NEON AOP flightlines.
        1) Apply BRDF/topographic corrections
            output format: *.envi 
        2)  Convert envi to tif
            output format: *.tif
        3) Crop flightlines to AOP tile bounds
            output format: *.tif
        4) Merge tiffs that overlap spatially
            output format: *.tif
    """
    r_source = ro.r['source']
    r_source(str(Path(__file__).resolve().parent/'hyperspectral_helper.R'))

    flightline_h5_path = os.path.join(data_raw_aop_path,site,year_aop,"hs_flightline")
    
    # 1) Apply BRDF/topographic corrections
    log.info(f'Applying BRDF/topo corrections for: {site} {year_inv}')
    # move_red_yellow_subset(site, flightline_h5_path) 
    # flightline_envi_path = convert_hdf_to_envi(flightline_h5_path, 
    #                                            os.path.join(data_int_path,site,year_inv,
    #                                                         'hs_flightline/'))
    flightline_envi_path = os.path.join(data_int_path,site,year_inv, 'hs_flightline/')
    # config_path = create_config(flightline_envi_path)
    # implement_brdf_correction(config_path) 

    # 2)  Convert corrected envi to tif
    log.info(f'Converting from envi to tif format for: {site} {year_inv}')
    convert_envi_to_tif(site, flightline_h5_path, flightline_envi_path) 

    # 3) Crop flightlines to AOP tile bounds
    log.info(f'Cropping flightline tifs to AOP tile bounds for: {site} {year_inv}')
    crop_flightlines_to_tiles = ro.r('crop_flightlines_to_tiles')
    # ???? = crop_flightlines_to_tiles(site, year_inv, data_raw_aop_path, flightline_envi_path) #ais

    # 4) Merge tiffs that overlap spatially
    log.info(f'Merging spatially overlapping tiffs for: {site} {year_inv}')
    

def prep_manual_training_data(site, year, data_raw_inv_path, data_int_path, biomass_path):
    """
    Clean and organize manually delineated tree crowns with PFT labels
    """

    # Create features (points or polygons) for each tree 
    log.info(f'Creating tree crown training data features for: {site} {year}')
    r_source = ro.r['source']
    # r_source(str(Path(__file__).resolve().parent/'inventory_helper.R'))
    # match_species_to_pft = ro.r('match_species_to_pft') 
    r_source(str(Path(__file__).resolve().parent/'hyperspectral_helper.R'))
    # list_tiles_w_veg = ro.r('list_tiles_w_veg') 
    
    # Create tree polygons
    prep_manual_crown_delineations = ro.r('prep_manual_crown_delineations') 
    training_shp_path = prep_manual_crown_delineations(site, year, data_raw_inv_path, data_int_path, 
                                                   biomass_path)
    
    #ais abandoned converting this function from R to python because loading shapefiles with gpd lost the column data types. also tough to use list_tiles_w_veg function in python

    log.info('Clipped tree crown polygons data saved at: '
             f'{training_shp_path}')
    
    return training_shp_path


def prep_aop_imagery(site, year, hs_type, hs_path, tif_path, data_int_path, use_tiles_w_veg):
    """ Prepare aop imagery for extracting descriptive features for classification by 
            (1) clean/mask out shadow and non-veg
            (2) stacking imagery 
     Analog to 05-prep_aop_imagery.R from https://github.com/earthlab/neon-veg 

     ais to do: Then plots one tile of the prepped imagery. Analog to 06-plot_aop_imagery.R

    Calls a function written in R

    Parameters
    ----------
    site : str
        Site name
    year : str
        Year of aop collection
    hs_path : str
        Path where either raw hyperspectral h5 files OR BRDF-corrected tifs are stored
    tif_path : str
        Path where the raw aop tif files are stored

    Returns
    -------
    (str)
        Path to the result folder for stacked HS data ready for classification
    """
    log.info(f'Prepping aop data for: {site} {year}')
    r_source = ro.r['source']
    r_source(str(Path(__file__).resolve().parent/'hyperspectral_helper.R'))

    # Clean/ mask imagery

    # Stack imagery
    prep_aop_imagery = ro.r('prep_aop_imagery')
    stacked_aop_path = prep_aop_imagery(site, year, hs_type, hs_path, tif_path, data_int_path, 
                                        use_tiles_w_veg)
    log.info('Stacked AOP data and saved at: '
             f'{stacked_aop_path}')
    return stacked_aop_path


def extract_spectra_from_polygon(site, year, shp_path, data_int_path, data_final_path, stacked_aop_path, 
                         use_case, aggregate_from_1m_to_2m_res, ic_type):
    """Create geospatial features (points, polygons with half the maximum crown diameter) 
        for every tree in the NEON woody vegetation data set. 
        Analog to 02-create_tree_features.R from https://github.com/earthlab/neon-veg 

        Generate polygons that intersect with independent pixels in the AOP data 
        Analog to 03-process_tree_features.R from https://github.com/earthlab/neon-veg 

        Extract features (remote sensing data) for each sample (pixel) within the 
        specified shapefile (containing polygons that correspond to trees at the NEON site)
        Analog to 07-extract_training_features.R from https://github.com/earthlab/neon-veg 
    """

    # # Create features (points or polygons) for each tree 
    # log.info(f'Extracting spectra from tree crowns for: {site} {year}')

    # training_data_dir = os.path.join(data_int_path, site, year, "training")
    # ic_type_path = os.path.join(data_final_path,site,year,ic_type)

    # # get a description of the shapefile to use for naming outputs
    # shapefile_description = os.path.splitext(os.path.basename(shp_path))[0]

    # # Specify destination for extracted features
    # if use_case == "train":
    #     extracted_features_path     = os.path.join(training_data_dir)
    #     extracted_features_filename = os.path.join(extracted_features_path,
    #                                              shapefile_description+"-extracted_features_inv.tif")                
    # elif use_case=="predict":
    #     extracted_features_path     = os.path.join(ic_type_path) 
    #     extracted_features_filename = os.path.join(extracted_features_path,
    #                                              shapefile_description+"-extracted_features.tif")
    # else:
    #     print("need to specify use_case")

    # # Only run if the extracted features do not exist
    # if not os.path.exists(extracted_features_filename):
    #     # Load shapefile
    #     shp_gdf = gpd.read_file(shp_path)

    #     # Compute centroid for each geometry
    #     shp_gdf["center_X"] = shp_gdf.geometry.centroid.x
    #     shp_gdf["center_Y"] = shp_gdf.geometry.centroid.y

    #     # List all .tif raster files in the directory
    #     stacked_aop_list = [os.path.join(stacked_aop_path, f) for f in os.listdir(stacked_aop_path) if f.endswith(".tif")]

    #     # create column to track shape ID 
    #     # if training, this is the tree crown boundary
    #     # if predicting, this is the plot boundary)
    #     if use_case == "train":
    #         shp_gdf["shapeID"] = [f"tree_crown_{i}" for i in range(len(shp_gdf))]
    #     elif use_case == "predict":
    #         if "PLOTID" in shp_gdf.columns:
    #             shp_gdf = shp_gdf.rename(columns={"PLOTID": "shapeID"})
    #         elif "plotID" in shp_gdf.columns:
    #             shp_gdf = shp_gdf.rename(columns={"plotID": "shapeID"})
    #         else:
    #             print("Trying to predict without training data. Need to first switch use_case to train from predict")

    #     # Filter to tiles containing veg to speed up the next for-loop
    #     # ais does this code include tiles with only a crown from the ucla field data though?
    #     if use_case == 'train' or ic_type=='rs_inv_plots':
    #         # Read tiles_w_veg.txt and extract the first column as a list
    #         tiles_w_veg = pd.read_csv(os.path.join(training_data_dir, "tiles_w_veg.txt"), header=None)[0].tolist()
                    
    #         # Filter stacked_aop_list based on the presence of any tile name
    #         stacked_aop_list = [path for path in stacked_aop_list if any(tile in path for tile in tiles_w_veg)]

    #     # Loop through AOP tiles
    #     for stacked_aop_filename in stacked_aop_list:

    #         # Read current tile of stacked AOP data
    #         with rasterio.open(stacked_aop_filename) as src:
    #             stacked_aop_data = src.read()
    #             raster_transform = src.transform
    #             raster_crs = src.crs
    #             raster_bounds = src.bounds  # Get raster bounds

    #             # Construct the easting northing string for naming outputs
    #             east_north_string = f"{round(raster_bounds[0])}_{round(raster_bounds[1])}"
                
    #             east_north_tif_path = os.path.join(extracted_features_path, f"extracted_features_mask_{east_north_string}_{shapefile_description}.tif")
                
    #             # If CSV file already exists, skip
    #             if os.path.exists(east_north_tif_path):
    #                 continue

    #             # Check which polygons overlap with the raster
    #             shp_gdf = shp_gdf[shp_gdf.intersects(box(*raster_bounds))]

    #             # If no polygons overlap, exit
    #             if shp_gdf.empty:
    #                 print("No overlapping polygons found.")
    #                 continue
    #             else:
    #                 # Convert overlapping polygons to GeoJSON format
    #                 shapes_geojson = [mapping(geom) for geom in shp_gdf.geometry]
                    
    #                 # Mask the raster using the overlapping polygons
    #                 with rasterio.open(stacked_aop_filename) as src:
    #                     masked_raster, masked_transform = rasterio.mask.mask(dataset=src, shapes=shapes_geojson, crop=False, nodata=np.nan)
                    
    #                 # Save the masked raster as a new GeoTIFF file
    #                 with rasterio.open(
    #                     east_north_tif_path, "w", driver="GTiff",
    #                     height=masked_raster.shape[1], width=masked_raster.shape[2],
    #                     count=masked_raster.shape[0], dtype=str(masked_raster.dtype),
    #                     crs=raster_crs, transform=masked_transform
    #                 ) as dst:
    #                     dst.write(masked_raster)
                    
    #                 print(f"Masked raster saved to {east_north_tif_path}")

    #     # combine all extracted features into a single .tif
    #     paths_ls = glob.glob(os.path.join(extracted_features_path, "*.tif"))
        
    #     # refine the output csv selection 
    #     tifs = [path for path in paths_ls if f"000_{shapefile_description}.tif" in path]
        
    #     # Open and merge TIF files
    #     src_files_to_mosaic = [rasterio.open(tif) for tif in tifs]
    #     mosaic, out_trans = merge(src_files_to_mosaic)
    #     # Open and merge TIF files
    #     src_files_to_mosaic = [rasterio.open(tif) for tif in tifs]
    #     mosaic, out_trans = merge(src_files_to_mosaic)

    #     # Get metadata from the first file
    #     out_meta = src_files_to_mosaic[0].meta.copy()   
    #     # Get metadata from the first file
    #     out_meta = src_files_to_mosaic[0].meta.copy()   

    #     # Update metadata for the merged file
    #     out_meta.update({
    #         "driver": "GTiff",
    #         "height": mosaic.shape[1],
    #         "width": mosaic.shape[2],
    #         "transform": out_trans
    #     })

    #     # Write the merged TIF file
    #     with rasterio.open(extracted_features_filename, "w", **out_meta) as dest:
    #         dest.write(mosaic)

    #     # Close input files
    #     for src in src_files_to_mosaic:
    #         src.close()

    #     # Delete the individual TIF files for each tile
    #     for tif in tifs:
    #         os.remove(tif)

    r_source = ro.r['source']
    r_source(str(Path(__file__).resolve().parent/'hyperspectral_helper.R'))
    extract_spectra_from_polygon_r = ro.r('extract_spectra_from_polygon_r')

    # Extract training data from AOP data with tree polygons
    training_spectra_path = extract_spectra_from_polygon_r(site=site, 
                                                         year=year, 
                                                         data_int_path=data_int_path, 
                                                         data_final_path=data_final_path, 
                                                         stacked_aop_path=stacked_aop_path, 
                                                         shp_path=shp_path,
                                                         use_case=use_case, 
                                                         aggregate_from_1m_to_2m_res=aggregate_from_1m_to_2m_res,
                                                         ic_type=ic_type)

    log.info('Spectral features for training data saved at: '
             f'{training_spectra_path}')
    
    return training_spectra_path



def extract_spectra_3darray(raster_path, training_shp):
    # generated by cborg, adapted by AIS March 2025

    # Load the multiband raster
    with rasterio.open(raster_path) as src:
        # Prepare a 3D array to store training data (bands, height, width)
        training_array = []

        # Create a combined mask for all polygons
        combined_mask = np.zeros(src.shape, dtype=bool)

        for _, row in training_shp.iterrows():
            geometry = row['geometry']
            # Create a mask for the raster using the polygon
            mask = geometry_mask([geometry], invert=True, transform=src.transform, out_shape=src.shape, all_touched=True)

            # Combine the mask
            combined_mask |= mask  # Logical OR to combine masks

        # Read the raster values
        raster_values = src.read()

        # Extract pixel values for the masked area
        for band in range(raster_values.shape[0]):
            band_values = raster_values[band][combined_mask]
            training_array.append(band_values)

        # Convert the list of arrays into a 3D numpy array
        training_array = np.array(training_array)  # Shape will be (bands, num_samples)

        return training_array, combined_mask



def _cm_analysis(y_true, y_pred, labels, savefile):
    cm = confusion_matrix(y_true, y_pred, labels=labels, normalize='true')
    disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=labels)
    fig, ax = plt.subplots(figsize=(8, 8))
    disp.plot(ax=ax, cmap=plt.cm.Blues)
    plt.savefig(str(savefile) + '_CMnorm_test.png', dpi=300, bbox_inches='tight')
    plt.close()


def _class_report(y_true, y_pred, labels, savefile):
    report = classification_report(y_true, y_pred, labels=labels, output_dict=True)
    report_df = pd.DataFrame(report).transpose()

    # Per-class specificity: TN / (TN + FP)
    cm = confusion_matrix(y_true, y_pred, labels=labels)
    total = cm.sum()
    spec_per_class = {}
    for i, label in enumerate(labels):
        tp = cm[i, i]
        fp = cm[:, i].sum() - tp
        fn = cm[i, :].sum() - tp
        tn = total - tp - fp - fn
        spec_per_class[str(label)] = tn / (tn + fp) if (tn + fp) > 0 else np.nan

    supports = report_df.loc[[str(l) for l in labels], 'support'].values
    macro_spec    = np.nanmean(list(spec_per_class.values()))
    weighted_spec = np.average(list(spec_per_class.values()), weights=supports)

    report_df['specificity'] = np.nan
    for label_str, spec in spec_per_class.items():
        report_df.loc[label_str, 'specificity'] = spec
    report_df.loc['macro avg',    'specificity'] = macro_spec
    report_df.loc['weighted avg', 'specificity'] = weighted_spec

    report_df.to_csv(str(savefile) + '_ClassReport_test.csv')


def _plot_featimportance(importance, feat_names, savefile):
    fi_df = pd.DataFrame({'feature': feat_names, 'importance': importance})
    fi_df = fi_df.sort_values('importance', ascending=False)
    plt.figure(figsize=(10, 6))
    sns.barplot(x='importance', y='feature', data=fi_df, color='gray')
    sns.despine(offset=10)
    plt.xlabel('Feature Importance')
    plt.ylabel('Feature')
    plt.savefig(str(savefile) + '_FeatImp.png', dpi=300, bbox_inches='tight')
    plt.close()


def train_pft_classifier(sites, data_int_path, pcaInsteadOfWavelengths, ntree,
                         randomMinSamples, independentValidationSet):    
    """Train a Random Forest (RF) model to classify tree PFT using in-situ tree
        measurements for PFT labels and remote sensing data as descriptive features
        Analog to 08-classify_species.R from https://github.com/earthlab/neon-veg 

        Train random forest model and assess PFT classification accuracy.
        Model outputs will be written to a folder within the training directory 
        starting with "rf_" followed by a description of each shapefile 
        containing points or polygons per tree. 

        Evaluate classifier performance: Out-of-Bag accuracy, independent
        validation set accuracy, and Cohen's kappa. Generate confusion matrix
        to show the classification accuracy for each PFT. 
        Analog to 09-assess_accuracy.R from https://github.com/earthlab/neon-veg 

        ntree 
          RF tuning parameter, number of trees to grow. default value 500
    
        randomMinSamples 
          To reduce bias, this boolean variable indicates whether the same number 
          of samples are randomly selected per PFT class. Otherwise,
          all samples per class are used for training the classifier. 
    
        independentValidationSet=T
          if TRUE, keep separate set for validation
          ais this doesnt work if for F - troubleshoot this later

    """
    
    log.info(f'Training model on sites: {sites}')

    # Directory to save the classification outputs 
    rf_output_dir = os.path.join(data_int_path, "rf_dir")
    if not os.path.exists(rf_output_dir):
        os.makedirs(rf_output_dir) 

    rf_model_path = os.path.join(rf_output_dir, f"rf_model.joblib")
    if not os.path.exists(rf_model_path):
    
        features_df = pd.DataFrame()
        for site in sites:
            site_year_paths = [os.path.join(data_int_path, site, year) for year in os.listdir(os.path.join(data_int_path, site))]
            for site_year_path in site_year_paths:

                if not os.path.exists(os.path.join(site_year_path, "stacked_aop/wavelengths.txt")):
                    continue

                # Load data
                wavelengths = pd.read_csv(os.path.join(site_year_path, "stacked_aop/wavelengths.txt"))
                # Stacked AOP layer names
                stacked_aop_layer_names = pd.read_csv(os.path.join(site_year_path, "stacked_aop/stacked_aop_layer_names.txt"))
                # Labelled, half-diam crown polygons to be used for training/validation
                shapefile_description = os.path.basename(os.path.join(site_year_path, "training/ref_labelled_crowns.shp")).split('.')[0]
                # Csv file containing extracted features
                extracted_features_filename = os.path.join(site_year_path, "training/ref_labelled_crowns-extracted_features_inv.csv")
                if not os.path.exists(extracted_features_filename):
                    continue

                # Filter out unwanted wavelengths
                wavelength_lut = filter_out_wavelengths(wavelengths=wavelengths['wavelengths'].tolist(), 
                                                        layer_names=stacked_aop_layer_names['stacked_aop_layer_names'].tolist())
                
                # features and label to use in the RF models
                featureNames = ["shapeID", "pft"] + wavelength_lut['xwavelength'].tolist() + [name for name in stacked_aop_layer_names['stacked_aop_layer_names'].tolist() if not name.isdigit()]
            
                # Prep extracted features csv for RF model
                extracted_features_X_df = pd.read_csv(extracted_features_filename)
                extracted_features_X_df = extracted_features_X_df[featureNames] 
                # Remove X from wavelength column names
                features_temp = extracted_features_X_df.rename(columns=lambda x: x[1:] if x.startswith('X') else x)                  
                # filter the data to contain only the features of interest 
                # features_temp = extracted_features_df[featureNames] 
                # features_temp.columns = featureNames[len(wavelengths):] + [f"X{i+1}" for i in range(len(wavelengths))] 
                
                # make sure that column names are the same so we can rbind them across site-year
                numeric_columns = [col for col in features_temp.columns if col.isdigit()]
                features_temp.columns = features_temp.columns.map(lambda x: f'X{numeric_columns.index(x)+1}' if x in numeric_columns else x)

                features_df = pd.concat([features_df, features_temp])

        # Remove any rows with NA   
        features_df.dropna(inplace=True)
            
        # Prepare arrays — shapeID as groups so all pixels from one crown
        # stay on the same side of every train/test split
        feature_cols  = [c for c in features_df.columns if c not in ('shapeID', 'pft')]
        x_col_idx     = [i for i, c in enumerate(feature_cols) if c.startswith('X')]
        other_col_idx = [i for i, c in enumerate(feature_cols) if not c.startswith('X')]
        X             = features_df[feature_cols].values
        # X      = features_df.loc[:, features_df.columns.str.startswith('X')].values
        y      = features_df['pft'].values
        groups = features_df['shapeID'].values

        print("Class distribution:")
        print(pd.Series(y).value_counts(normalize=True))

        if randomMinSamples:
            min_count = pd.Series(y).value_counts().min()
            keep_idx = (pd.DataFrame({'y': y, 'g': groups})
                        .groupby('y')
                        .apply(lambda df: df.sample(min_count, random_state=42))
                        .index.get_level_values(1))
            X, y, groups = X[keep_idx], y[keep_idx], groups[keep_idx]

        # Pipeline keeps preprocessing inside folds — no leakage
        preprocessor = ColumnTransformer(transformers=[
            ('spectral', Pipeline([
                ('scaler', StandardScaler()),
                ('pca',    PCA(n_components=0.99, whiten=True)),
            ]), x_col_idx),
            ('other', StandardScaler(), other_col_idx),
        ])
        pipeline = Pipeline(steps=[
            ('preprocessor', preprocessor),
            ('classifier',   RF(class_weight='balanced', random_state=42)),
        ])
        # pipeline = Pipeline(steps=[
        #     ('scaler',     StandardScaler()),
        #     ('pca',        PCA(n_components=0.99, whiten=True)),
        #     ('classifier', RF(class_weight='balanced', random_state=42)),
        # ])
        param_grid = {
            'classifier__n_estimators':  [50, 100, 500, 1000],
            'classifier__max_features':  ['sqrt', 'log2', None],
            'classifier__max_depth':     [None, 10, 20, 30],
            'classifier__min_samples_leaf': [1, 2, 4],
        }

        # 30 rounds of random 80/20 group splits — each shapeID can appear in
        # multiple test sets, giving multiple predictions per pixel for majority voting
        # Inner: StratifiedGroupKFold — group-aware hyperparameter tuning per round
        cv_shuffle = GroupShuffleSplit(n_splits=30, test_size=0.2, random_state=42)
        cv_inner   = StratifiedGroupKFold(n_splits=3)

        predicted_targets = np.array([])
        actual_targets    = np.array([])
        outer_results     = []
        best_params_list  = []
        pixel_preds       = {i: [] for i in range(len(X))}

        for rnd, (train_idx, test_idx) in enumerate(cv_shuffle.split(X, y, groups=groups)):
            X_tr, X_te = X[train_idx], X[test_idx]
            y_tr, y_te = y[train_idx], y[test_idx]
            g_tr       = groups[train_idx]

            search = GridSearchCV(estimator=pipeline, param_grid=param_grid,
                                  scoring='f1_weighted', cv=cv_inner, refit=True)
            result = search.fit(X_tr, y_tr, groups=g_tr)
            p_te   = result.best_estimator_.predict(X_te)

            for idx, pred in zip(test_idx, p_te):
                pixel_preds[idx].append(pred)

            predicted_targets = np.append(predicted_targets, p_te)
            actual_targets    = np.append(actual_targets, y_te)
            round_f1 = f1_score(y_te, p_te, average='weighted')
            outer_results.append(round_f1)
            best_params_list.append(result.best_params_)
            log.info(f'Round {rnd+1}/30: f1={round_f1:.3f}, inner_best={result.best_score_:.3f}, '
                     f'params={result.best_params_}')

        log.info(f'Repeated CV F1: {np.mean(outer_results):.3f} +/- {np.std(outer_results):.3f}')

        # Majority vote and agreement fraction per pixel
        records = []
        for idx, preds in pixel_preds.items():
            if not preds:
                continue
            c         = Counter(preds)
            majority  = c.most_common(1)[0][0]
            agreement = c[majority] / len(preds)
            records.append({'array_idx': idx, 'shapeID': groups[idx],
                            'majority_pft': majority,
                            'agreement': agreement, 'n_predictions': len(preds)})
        uncertainty_df = pd.DataFrame(records)

        cv_stats = pd.DataFrame(best_params_list)
        cv_stats.insert(0, 'round', range(1, len(outer_results) + 1))
        cv_stats.insert(1, 'round_f1', outer_results)
        # cv_stats = pd.DataFrame({
        #     'round':             list(range(1, len(outer_results) + 1)),
        #     'round_f1':          outer_results,
        #     'best_n_estimators': [p['classifier__n_estimators'] for p in best_params_list],
        # })
        cv_stats.to_csv(os.path.join(rf_output_dir, 'cv_results_test.csv'), index=False)

        # Fit final model on all data with the same group-aware inner CV
        final_search = GridSearchCV(estimator=pipeline, param_grid=param_grid,
                                    scoring='f1_weighted', cv=cv_inner, refit=True)
        final_search.fit(X, y, groups=groups)
        rf_model = final_search.best_estimator_
        log.info(f'Final model best params: {final_search.best_params_}')
        dump(rf_model, rf_model_path)
        log.info(f'Trained PFT classifier saved: {rf_model_path}')

        # PCA visualisation using the fitted pipeline
        X_transformed = rf_model.named_steps['preprocessor'].transform(X)
        n_pcs         = rf_model.named_steps['preprocessor'].transformers_[0][1].named_steps['pca'].n_components_
        X_pca         = X_transformed[:, :n_pcs]
        pca_df   = pd.DataFrame(X_pca[:, :min(4, n_pcs)],
                                columns=[f'PC{i+1}' for i in range(min(4, n_pcs))])
        pca_df['pft'] = y
        for pc_a, pc_b in [('PC1', 'PC2'), ('PC2', 'PC3'), ('PC3', 'PC4')]:
            if pc_a not in pca_df.columns or pc_b not in pca_df.columns:
                break
            plt.figure(figsize=(8, 6))
            for pft_val in pca_df['pft'].unique():
                sub = pca_df[pca_df['pft'] == pft_val]
                plt.scatter(sub[pc_a], sub[pc_b], label=pft_val, alpha=0.6)
            plt.xlabel(pc_a); plt.ylabel(pc_b)
            plt.title(f'{pc_a} vs {pc_b}')
            plt.legend()
            plt.savefig(os.path.join(rf_output_dir, f'{pc_a.lower()}_vs_{pc_b.lower()}.png'))
            plt.close()

        # Uncertainty visualization — agreement histogram (always produced)
        plt.figure(figsize=(10, 5))
        for pft_val in sorted(uncertainty_df['majority_pft'].unique()):
            sub = uncertainty_df[uncertainty_df['majority_pft'] == pft_val]
            plt.hist(sub['agreement'], bins=20, alpha=0.5, label=pft_val, range=(0, 1))
        plt.xlabel('Agreement fraction across CV rounds')
        plt.ylabel('Pixel count')
        plt.title('Per-pixel prediction agreement (30 CV rounds, majority vote)')
        plt.legend()
        plt.savefig(os.path.join(rf_output_dir, 'uncertainty_agreement_hist_test.png'), dpi=300, bbox_inches='tight')
        plt.close()

        # Spatial uncertainty map — produced only if easting/northing exist in features_df
        east_col  = next((c for c in features_df.columns if c.lower() in ('easting', 'x')),  None)
        north_col = next((c for c in features_df.columns if c.lower() in ('northing', 'y')), None)
        if east_col and north_col:
            coord_src = features_df.iloc[keep_idx] if randomMinSamples else features_df
            coords    = coord_src[[east_col, north_col]].reset_index(drop=True)
            idxs      = uncertainty_df['array_idx'].values
            uncertainty_df['easting']  = coords.iloc[idxs][east_col].values
            uncertainty_df['northing'] = coords.iloc[idxs][north_col].values

            pft_classes = sorted(uncertainty_df['majority_pft'].unique())
            pft_to_int  = {p: i for i, p in enumerate(pft_classes)}
            cmap_pft    = plt.cm.get_cmap('tab10', len(pft_classes))

            fig, axes = plt.subplots(1, 2, figsize=(16, 6))
            sc0 = axes[0].scatter(uncertainty_df['easting'], uncertainty_df['northing'],
                                  c=uncertainty_df['agreement'], cmap='RdYlGn', s=1, vmin=0, vmax=1)
            plt.colorbar(sc0, ax=axes[0], label='Agreement fraction')
            axes[0].set_title('Prediction agreement')
            axes[0].set_xlabel('Easting'); axes[0].set_ylabel('Northing')

            sc1 = axes[1].scatter(uncertainty_df['easting'], uncertainty_df['northing'],
                                  c=uncertainty_df['majority_pft'].map(pft_to_int),
                                  cmap=cmap_pft, s=1, vmin=0, vmax=len(pft_classes) - 1)
            cb  = plt.colorbar(sc1, ax=axes[1], ticks=range(len(pft_classes)), label='PFT')
            cb.set_ticklabels(pft_classes)
            axes[1].set_title('Majority-vote PFT')
            axes[1].set_xlabel('Easting'); axes[1].set_ylabel('Northing')

            plt.tight_layout()
            plt.savefig(os.path.join(rf_output_dir, 'uncertainty_spatial_map_test.png'), dpi=300, bbox_inches='tight')
            plt.close()
        else:
            log.info('Easting/northing columns not found in features; spatial uncertainty map skipped.')

        uncertainty_df.to_csv(os.path.join(rf_output_dir, 'uncertainty_per_pixel_test.csv'), index=False)

        # Evaluation on accumulated held-out CV predictions
        savefile   = os.path.join(rf_output_dir, 'rf')
        class_list = np.unique(actual_targets)
        _cm_analysis(actual_targets, predicted_targets, class_list, savefile)
        _class_report(actual_targets, predicted_targets, class_list, savefile)

        n_pcs       = rf_model.named_steps['preprocessor'].transformers_[0][1].named_steps['pca'].n_components_
        pc_names    = [f'PC{i+1}' for i in range(n_pcs)]
        other_names = [feature_cols[i] for i in other_col_idx]
        feat_names  = pc_names + other_names
        importances = rf_model.named_steps['classifier'].feature_importances_
        _plot_featimportance(importances, feat_names, savefile)

        overall_acc = accuracy_score(actual_targets, predicted_targets)
        nir = pd.Series(y).value_counts(normalize=True).max()

        summary_lines = [
            f"Repeated CV F1 (30 rounds, mean +/- std): {np.mean(outer_results):.4f} +/- {np.std(outer_results):.4f}",
            f"Accuracy (held-out CV predictions): {overall_acc:.4f}",
            f"No-information rate (majority class frequency): {nir:.4f}",
            f"Final model best params: {final_search.best_params_}",
            "\nPer-round results:",
            cv_stats.to_string(index=False),
            "\nClassification report (held-out CV predictions):",
            classification_report(actual_targets, predicted_targets),
        ]
        with open(os.path.join(rf_output_dir, 'rf_summary_statistics_test.txt'), 'w') as f:
            f.write('\n'.join(summary_lines))
    
    # else:
    #     rf_model = load(rf_model_path)

    return rf_model_path   



def plot_cv_indices(cv, X, y, groups, ax, n_splits, cmap_data, cmap_cv, lw=10):
    """Create a sample plot for indices of a cross-validation object.
    from: https://scikit-learn.org/stable/auto_examples/model_selection/plot_cv_indices.html"""
    #use_groups = "Group" in type(cv).__name__
    # groups = group if use_groups else None
    # Generate the training/testing visualizations for each CV split
    for ii, (tr, tt) in enumerate(cv.split(X=X, y=y, groups=groups)):
        # Fill in indices with the training/test groups
        indices = np.array([np.nan] * len(X))
        indices[tt] = 1
        indices[tr] = 0

        # Visualize the results
        ax.scatter(
            range(len(indices)),
            [ii + 0.5] * len(indices), 
            c=indices,
            marker="_",
            lw=lw,
            cmap=cmap_cv,
            vmin=-0.2,
            vmax=1.2,
        )

    # Plot the data classes and groups at the end
    ax.scatter(
        range(len(X)), [ii + 1.5] * len(X), c=y, marker="_", lw=lw, cmap=cmap_data
    )
    ax.scatter(
        range(len(X)), [ii + 2.5] * len(X), c=groups, marker="_", lw=lw, cmap=cmap_data
    )

    # Formatting
    yticklabels = list(range(n_splits)) + ["class", "group"]
    ax.set(
        yticks=np.arange(n_splits + 2) + 0.5,
        yticklabels=yticklabels,
        xlabel="Sample index",
        ylabel="CV iteration",
        ylim=[n_splits + 2.2, -0.2],
        # xlim=[0, 100],
    )
    ax.set_title("{}".format(type(cv).__name__), fontsize=15)
    return ax



def fit_RF_CV_class(X, y, groups, k_fold, pca, savefile):
    # Written by Nicola Falco
    # Adapted by Anna Spiers Nov 2024
  
    ##################### Define evaluation procedure for CV
    cv = StratifiedKFold(n_splits=k_fold)#, shuffle=True, random_state=10210)

    ##################### test/training 
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=111)

    #############################################################################
    #####################  cross-validation
    if pca == True:
        pipeline = Pipeline(steps=[
            ['pca', PCA(n_components=0.99, whiten= True)],
            # ['scaler', MinMaxScaler()],
            ['classifier', RF(class_weight='balanced')] 
            ])
        savefile = (str(savefile) + '_PCA')

    else:
        pipeline = Pipeline(steps=[
            # ['scaler', MinMaxScaler()],
            ['classifier', RF(class_weight='balanced')]
            ])
        savefile = (str(savefile))
        
    ##################### Define grid search RF Parameter optimization
    # Number of trees in random forest
    n_estimators = [int(x) for x in np.linspace(start = 1000, stop = 10000, num = 4)] #5

    # Number of features to consider at every split
    max_features = ['auto', 'sqrt']

    # Maximum number of levels in tree
    max_depth = [int(x) for x in np.linspace(10, 110, num = 4)] #5
    max_depth.append(None)

    # Minimum number of samples required to split a node
    # min_samples_split = [2, 5, 10]

    # Minimum number of samples required at each leaf node
    min_samples_leaf = [1, 2, 4]

    # Method of selecting samples for training each tree
    bootstrap = [True]#, False]

    param_rf ={'classifier__n_estimators': n_estimators,
               'classifier__max_features': max_features,
               'classifier__min_samples_leaf': min_samples_leaf,
               'classifier__max_depth': max_depth,
               # 'classifier__min_samples_split': min_samples_split,
               # 'classifier__min_samples_leaf': min_samples_leaf,
               'classifier__bootstrap': bootstrap
               }

    search = GridSearchCV(estimator=pipeline, param_grid=param_rf,
                           scoring='f1_weighted', cv=cv, refit=True, verbose=3, n_jobs = 128)
  
    ##################### Fit the search to the data
    # Fit the grid search to the data
    search.fit(X_train, y_train)
    search.best_params_

    ##################### get the best performing model fit on the whole training set
    best_model = search.best_estimator_
  
    ##################### Prediction
    p_test = best_model.predict(X_test)

    ##################### Report performance
    #### average CM and classification report
    class_list= np.unique(y_test)

    cm_analysis(y_test, p_test, class_list, savefile)
    class_report(y_test, p_test, class_list,savefile)
    
    return best_model, X_train, X_test



def cm_analysis(y_true, y_pred, labels,savefile): 

    SMALL_SIZE = 12
    MEDIUM_SIZE = 15
    BIGGER_SIZE = 18
   
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
 
    sns.set(style="whitegrid")

    ##############################
 
    # compute the CM
    cm = confusion_matrix(y_true, y_pred, labels=labels, normalize='true')
    disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=labels)
    fig, ax = plt.subplots(figsize=(8,8))
    disp.plot(ax=ax, cmap=plt.cm.Blues)

    save2png= (str(savefile) + '_CMnorm_test.png')
    plt.savefig(save2png, dpi=300, bbox_inches='tight')

    

def class_report(y_test, y_pred, labels,savefile):
    report = classification_report(y_test,y_pred,labels=labels,output_dict=True)
    df_report = pd.DataFrame(report).transpose()
    save2csv = (str(savefile) + '_ClassReport_test.csv')
    df_report.to_csv(save2csv)




def serial_index_to_coordinates(serial_index, raster_width=1000):
    """
    Translates a serial index of a raster into coordinates of the raster.

    Parameters:
    serial_index (int): The serial index of the pixel in the raster.
    raster_width (int): The width of the raster.

    Returns:
    tuple: A tuple containing the x and y coordinates of the pixel in the raster.
    """
    x = (serial_index // raster_width)-1
    y = 1000 - (serial_index % raster_width) 
    return x, y





def filter_out_wavelengths(wavelengths, layer_names):
    # define the "bad bands" wavelength ranges in nanometers, where atmospheric 
    # absorption creates unreliable reflectance values. 
    # bad_band_window_1 = (1340, 1445)
    # bad_band_window_2 = (1790, 1955)
    wavelengths_np = np.array(wavelengths)

    # remove the bad bands from the list of wavelengths 
    remove_bands = wavelengths_np[(wavelengths_np > 1340) & 
                               (wavelengths_np < 1445) | 
                               (wavelengths_np > 1790) & 
                               (wavelengths_np < 1955)]
    
    # Make sure printed wavelengths and stacked AOP wavelengths match
    if not np.allclose(np.round(wavelengths_np), np.array([float(name) for name in layer_names[:len(wavelengths_np)]])):
        print("wavelengths do not match between wavelength.txt and the stacked imagery")
    
    # create a LUT that matches actual wavelength values with the column names,
    # X followed by the rounded wavelength values. 
    # Remove the rows that are within the bad band ranges. 
    wavelength_lut = pd.DataFrame({'wavelength': wavelengths_np,
                                  'xwavelength': ['X' + str(round(w)) for w in wavelengths_np]})[~np.in1d(wavelengths_np, remove_bands)]
    
    return wavelength_lut


