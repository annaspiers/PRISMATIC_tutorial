import logging
import os
import sys
import glob
import re

import pandas as pd
import geopandas as gpd
import rasterio
from rasterio import mask
from rasterio.plot import show
from rasterio.features import geometry_mask
from rasterio.transform import rowcol
import numpy as np
from shapely.geometry import Point
from shapely.geometry import box
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec #plot predicted classes next to rgb
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
import joblib
from joblib import dump, load

import rpy2.robjects as ro

from pathlib import Path

# add environ
conda_env_path = Path(sys.executable).parent.parent
os.environ['PROJ_LIB'] = str(conda_env_path/'share'/'proj')

# custom functions in another script
from initialize.lad import prep_lad
from initialize.hyperspectral import filter_out_wavelengths, \
                                        extract_spectra_from_polygon

log = logging.getLogger(__name__)

def generate_initial_conditions(site, year_inv, year_aop, data_raw_aop_path, data_int_path,

                                data_final_path, rf_model_path, stacked_aop_path, biomass_path,
                                use_case, ic_type, ic_type_path, n_plots, min_distance, plot_length, 
                                aggregate_from_1m_to_2m_res, pcaInsteadOfWavelengths, multisite):                                            
    """Generate FATES initial conditions from remote sensing data in three scenarios:
            1) inventory-based IC's over inventory plots (ic_type = field_inv_plots)
            2) remote-sensing-based IC's from plots used for forest inventory 
                    (ic_type = rs_inv_plots)
            3) remote-sensing-based IC's from plots gridded wall to wall across 
                    the eddy covariance tower spatial extent (ic_type = rs_tower_ftpt)
            4) remote-sensing-based IC's from plots randomly generated across 
                the AOP spatial extent (ic_type = rs_random_plots)
            5) remote-sensing-based IC's from plots gridded wall to wall across 
                the AOP spatial extent (ic_type = rs_wall2wall)

    We generate files formatted to initialize FATES: cohort (.css) and patch (.pss) files
    """

    if not os.path.exists(Path(data_final_path,site,year_inv)):
        os.makedirs(Path(data_final_path,site,year_inv))   

    log.info(f'Generating initial conditions for: {site} {year_inv} {ic_type}')
    r_source = ro.r['source']
    r_source(str(Path(__file__).resolve().parent/'initial_conditions_helper.R'))
    generate_pcss_field = ro.r('generate_pcss_field') 
    generate_pcss_rs = ro.r('generate_pcss_rs') 
    
    if (ic_type=="field_inv_plots"):
        cohort_path, patch_path = generate_pcss_field(site=site, 
                                    year=year_inv, 
                                    data_int_path=data_int_path, 
                                    biomass_path=biomass_path,
                                    ic_type=ic_type, 
                                    ic_type_path=ic_type_path)

    else :  #any of the rs_* ic_types      
        ### 1) Specify spatial extent of initial conditions and prepare spatial data by
        # clipping point clouds to and extracting raster data from plots
        plots_shp_path, plots_laz_path, extracted_features_csv = prepare_spatial_data_for_ic(
                                            site, year_inv, year_aop, ic_type, data_raw_aop_path, 
                                            data_int_path, stacked_aop_path, data_final_path, 
                                            ic_type_path, use_case, n_plots, min_distance, 
                                            plot_length, aggregate_from_1m_to_2m_res) 
        log.info('Plots, clipped laz, and clipped stacked rasters prepared for generating initial conditions')

        ### 2) Predict pfts on new plots
        # proxy to neon-veg-SOAPpfts/11-predict-pft.R
        # if ic_type==rs_inv_plots, use known PFTs from NEON inventory data 
        # if ic_type==rs_random_plots or rs_tower_ftpt, predict PFTs with RF classifier        
        classified_plot_PFTs_path = predict_pixel_PFTs(site=site, 
                                                       year_inv=year_inv, 
                                                       year_aop=year_aop,
                                                       data_raw_aop_path=data_raw_aop_path,
                                                       data_int_path=data_int_path,
                                                       stacked_aop_path=stacked_aop_path, 
                                                        ic_type=ic_type, 
                                                        plots_shp_path=plots_shp_path,
                                                        ic_type_path=ic_type_path, 
                                                        rf_model_path=rf_model_path,
                                                        pcaInsteadOfWavelengths=pcaInsteadOfWavelengths)
        # classified_plot_PFTs_path = predict_plot_level_PFTs(stacked_aop_path=stacked_aop_path, 
        #                                             ic_type=ic_type, ic_type_path=ic_type_path, 
        #                                             rf_model_path=rf_model_path,
        #                                             extracted_features_filename=extracted_features_csv[0],
        #                                             pcaInsteadOfWavelengths=pcaInsteadOfWavelengths)
        log.info('Classified PFTs for FATES initialization stored in: '
                f'{classified_plot_PFTs_path}')

        ### 3) Generate lad profiles on new plots 
        # if ic_type=="rs_inv_plots", use existing LAD profiles
        # or if ic_type=="rs_tower_ftpt","rs_random_plots","rs_wall2wall" generate lad profiles...
        if ic_type!="rs_inv_plots":
            if not glob.glob(str(ic_type_path+'/clipped_to_plots/'+'*lad.json') ):
                # Generate new LAD profiles
                prep_lad(laz_path=plots_laz_path, 
                            inventory_path=None, 
                            site=site, 
                            year=year_inv, 
                            output_path=ic_type_path, 
                            use_case="predict")
                #ais why can LAD not be processed for some sites?

        ### 4) Generate cohort and patch files
               
        # # Check for compatability between PFTs and cohorts in each patch
        # ?? = ??(patch_struct_summ_path)
        # Generate pcss files by assigining pfts across cohorts
        cohort_path, patch_path = generate_pcss_rs(site=site, 
                                        year=year_inv, 
                                        data_int_path=data_int_path, 
                                        biomass_path=biomass_path,
                                        ic_type=ic_type, 
                                        ic_type_path=ic_type_path, 
                                        plots_shp_path=plots_shp_path,
                                        classified_plot_PFTs_path=classified_plot_PFTs_path,
                                        plots_laz_path=plots_laz_path,
                                        multisite=multisite)
    
    log.info('Cohort and patch files for FATES initialization generated')

    return cohort_path, patch_path


def prepare_spatial_data_for_ic(site, year_inv, year_aop, ic_type, data_raw_aop_path, 
                                data_int_path, stacked_aop_path, data_final_path, ic_type_path,
                                use_case, n_plots, min_distance, plot_length, 
                                aggregate_from_1m_to_2m_res):
    # This function generates random plots across the remote sensing 
    # spatial extent and clips laz and rasters 
    r_source = ro.r['source']
    r_source(str(Path(__file__).resolve().parent/'initial_conditions_helper.R'))
    generate_random_plots = ro.r('generate_random_plots') 
    generate_wall2wall_plots = ro.r('generate_wall2wall_plots') 
    clip_norm_laz_to_shp = ro.r('clip_norm_laz_to_shp') 
    # clip_lidar_to_polygon_lidR = ro.r('clip_lidar_to_polygon_lidR') 
    
    ### Specify spatial extent of initial conditions
    if ic_type == "rs_inv_plots": 
        # use existing inventory plot extents
        aoi_shp_dir = os.path.join(data_raw_aop_path, site, year_aop, "shape")
        aoi_shp_matches = glob.glob(os.path.join(aoi_shp_dir, "*.shp"))
        if not aoi_shp_matches:
            raise FileNotFoundError(f"No .shp file found in {aoi_shp_dir}")
        aoi_shp_path = aoi_shp_matches[0]
        plots_shp_path = os.path.join(data_int_path, site, str(year_inv), "inventory_plots", "plots.shp")
        plots_laz_path = os.path.join(data_int_path, site, str(year_inv), "clipped_to_plots")

        plots_shp = gpd.read_file(plots_shp_path)
        aoi_shp   = gpd.read_file(aoi_shp_path)
        # Save figure with plots and AOI
        fig, ax = plt.subplots()
        aoi_shp.plot(ax=ax, color='none', edgecolor='black', linewidth=1)
        plots_shp.plot(ax=ax, color='red', edgecolor='red', linewidth=.75)
        plt.title('AOI and Plots for initial condition type ' + ic_type)
        plt.savefig(os.path.join(ic_type_path,'plot_coverage_across_AOI.png'), dpi=300, bbox_inches='tight')
                
    else:         
        plots_shp_path = os.path.join(ic_type_path, "ic_type_plots.shp")
        plots_laz_path = os.path.join(ic_type_path, "clipped_to_plots")

        if not os.path.exists(plots_shp_path):
            if not os.path.exists(plots_laz_path): os.makedirs(plots_laz_path)
            # Generate plots across AOI

            if ic_type == "rs_random_plots":
                aoi_shp_dir = os.path.join(data_raw_aop_path, site, year_aop, "shape")
                pattern = os.path.join(aoi_shp_dir, "*merged_tiles.shp")
                aoi_shp_path = glob.glob(pattern)[0]
                generate_random_plots(aoi_shp_path, ic_type,ic_type_path, n_plots, min_distance, plot_length)
                    
            if ic_type == "rs_tower_ftpt":
                aoi_shp_dir = os.path.join(data_raw_aop_path, site)                
                pattern = os.path.join(aoi_shp_dir, "*0p_flux_ftpt.shp")
                aoi_shp_path = glob.glob(pattern)[0]
                generate_wall2wall_plots(aoi_shp_path,plot_length,ic_type,ic_type_path)

            elif ic_type == "rs_wall2wall":
                tile_bbox = box(321000, 4097000, 322000, 4098000)
                aoi_gdf = gpd.GeoDataFrame(geometry=[tile_bbox], crs="EPSG:32611")
                aoi_shp_dir = os.path.join(data_raw_aop_path, site, year_aop, "shape")
                aoi_shp_path = os.path.join(aoi_shp_dir, "tile_321000_4097000.shp")
                aoi_gdf.to_file(aoi_shp_path)
                generate_wall2wall_plots(aoi_shp_path,plot_length,ic_type,ic_type_path)
            
        # Clip normalized point clouds to randomized plots    
        #if the number of clipped plots does not exqual the number of plot shapes, continue clipping laz
       # if len(os.listdir(plots_laz_path)) != len(gpd.read_file(plots_shp_path)): 
            clip_norm_laz_to_shp(site, year_inv, data_int_path, ic_type_path, plots_shp_path)
        #ais sometimes thereare going to be some patches dropped - what to do?
       
    # ais tried changing this function to python, but it ran incredibly slowly and didn't work 
    # put too much time into trying to make it work (see scrap below)
    # extracted_features_csv = extract_spectra_from_polygon(site=site,
    #                                                        year=year_inv,
    #                                                        data_int_path=data_int_path,
    #                                                        data_final_path=data_final_path,
    #                                                        stacked_aop_path=stacked_aop_path,
    #                                                        shp_path=plots_shp_path,
    #                                                        use_case=use_case,
    #                                                        aggregate_from_1m_to_2m_res=aggregate_from_1m_to_2m_res,                                                           
    #                                                        ic_type=ic_type)
    extracted_features_csv = None

    return (plots_shp_path, plots_laz_path, extracted_features_csv)



def predict_plot_level_PFTs(data_int_path, stacked_aop_path, ic_type_path, rf_model_path, 
                            extracted_features_filename,
                            ic_type, pcaInsteadOfWavelengths):

    classified_plot_PFTs_path = os.path.join(ic_type_path,"pfts_per_plot.csv")

    if not os.path.exists(classified_plot_PFTs_path):

        wavelengths = pd.read_csv(os.path.join(stacked_aop_path, 
                                               "wavelengths.txt"))['wavelengths'].tolist()
        stacked_aop_layer_names = pd.read_csv(os.path.join(stacked_aop_path, 
                                                           "stacked_aop_layer_names.txt"))['stacked_aop_layer_names'].tolist()
        #ais remove the following rows earlier in the workflow so I don't need to do it here
        stacked_aop_layer_names = [x for x in stacked_aop_layer_names]
            
        # Filter out unwanted wavelengths
        wavelength_lut = filter_out_wavelengths(wavelengths=wavelengths, layer_names=stacked_aop_layer_names)
        
        # features to use in the RF models
        featureNames = ["shapeID"] + wavelength_lut['xwavelength'].tolist() + [name for name in stacked_aop_layer_names if not name.isdigit()]
        
        # Prep extracted features csv for RF model
        extracted_features_df = pd.read_csv(extracted_features_filename)
        # filter the data to contain only the features of interest 
        features_df = extracted_features_df[featureNames] 

        # Remove any rows with NA   
        features_df.dropna(inplace=True)
        #ais why are there rows with NAs?
        
        ### Predict PFT ID for FATES patches 
        print("Predicting plot-level PFTs")

         # Load RF model
        rf_model = joblib.load(rf_model_path) 
                
        if pcaInsteadOfWavelengths:
            # remove the individual spectral reflectance bands from the training data
            features_noWavelengths = features_df.drop([col for col in features_df.columns if col.startswith("X")], axis=1)
            
            # Define evaluation procedure for CV
            # cv = StratifiedKFold(n_splits = 3) #ais make k_fold a global param #, shuffle=True, random_state=10210)
            
            # Load scaler and PCA used during training
            loaded_scaler = joblib.load(os.path.join(data_int_path, "rf_dir",'scaler.joblib'))
            loaded_pca = joblib.load(os.path.join(data_int_path, "rf_dir",'pca.joblib'))

            # PCA
            features_scaled = loaded_scaler.transform(features_df[[col for col in features_df.columns if col.startswith("X")]])
            # nPCs = sum(1 for item in rf_model.feature_names_in_ if item.startswith('PC'))
            # pca = PCA(n_components=nPCs) 
            # #ais ^ I may have 9 PCs used in training the RF model, but to get to 99% variance for only rs_inv plots, 
            # # one PC may be sufficient, so I'm forcing the same number of PCs as model
            features_pca = loaded_pca.transform(features_scaled) 
            features = pd.concat([features_noWavelengths.reset_index(drop=True), pd.DataFrame(features_pca, columns=[f"PC{i+1}" for i in range(nPCs)])],axis=1)
        else:
            features = features_df
            
        # Predict PFT ID
        features['pred_PFT'] = rf_model.predict(features.drop('shapeID',axis=1)) 
        
        # Summarize as percentages per plot
        features_pft_by_pct = features.groupby(['shapeID', 'pred_PFT']).size().reset_index(name='count')
        features_pft_by_pct['pct'] = features_pft_by_pct['count'] / features_pft_by_pct.groupby('shapeID')['count'].transform('sum')
        
        # Write to CSV
        features_pft_by_pct.to_csv(classified_plot_PFTs_path, index=False)

        # Notes w Nicola
        # save geometry of each pixel in addition to the coordinates of the plot (for plotting RGB behind later)
        # then each point predicted has a coordinate and I can map this point-based
    
    return classified_plot_PFTs_path


def read_raster(input_image):
    #author: Nicola Falco
    
    from osgeo import gdal
    
    # Tell GDAL to throw Python exceptions, and register all drivers
    gdal.UseExceptions()
    gdal.AllRegister()
    
    # Read in image
    img_gdal_obj = gdal.Open(input_image, gdal.GA_ReadOnly)
    
    # initialize img to zeros
    img = img_gdal_obj.ReadAsArray()
    
    return img, img_gdal_obj;



def image_vectorization_3d(img):
    #author: Nicola Falco
    
    new_shape = (img.shape[0], img.shape[1] * img.shape[2])
    img_vec = img[:, :, :].reshape(new_shape)
    img_vec = img_vec.T
    
    return img_vec



def array2gtiff(image_reference, map_img, format_img, savefile):
    #author: Nicola Falco

    from osgeo import gdal
    
    datatype = gdal.GDT_Float32
    
    img = gdal.Open(image_reference)
    cols = img.RasterXSize
    rows = img.RasterYSize

    driver = gdal.GetDriverByName(format_img)

    outDataRaster = driver.Create(savefile, cols, rows, 1, datatype)
    # sets same geotransform as input
    outDataRaster.SetGeoTransform(img.GetGeoTransform())
    # sets same projection as input
    outDataRaster.SetProjection(img.GetProjection())

    outDataRaster.GetRasterBand(1).WriteArray(map_img)
    outDataRaster.GetRasterBand(1).SetNoDataValue(-9999)

    outDataRaster.FlushCache()  # remove from memory
    del outDataRaster  # delete the data (not the actual geotiff)


def extract_coords(filename):
    match = re.search(r'(\d{6}_\d{7})', filename)
    return match.group(1) if match else None



def predict_pixel_PFTs(site, year_inv, year_aop, data_raw_aop_path, data_int_path, stacked_aop_path, plots_shp_path, ic_type_path, 
                            rf_model_path, ic_type, pcaInsteadOfWavelengths):

    classified_tiles_dir = os.path.join(data_int_path, site, year_inv, "rf_classified_tiles")
    if not os.path.exists(classified_tiles_dir):
        os.makedirs(classified_tiles_dir)
    classified_files = glob.glob(os.path.join(classified_tiles_dir, '*.tif'))
    classified_coords = {extract_coords(os.path.basename(path)) for path in classified_files if extract_coords(os.path.basename(path))}

    # Load relevant data
    stacked_aop_list = glob.glob(os.path.join(stacked_aop_path,"*.tif"))
    stacked_aop_layer_names = pd.read_csv(os.path.join(stacked_aop_path, 
                                                            "stacked_aop_layer_names.txt"))['stacked_aop_layer_names'].tolist()
    stacked_aop_coords = {extract_coords(os.path.basename(path)) for path in stacked_aop_list if extract_coords(os.path.basename(path))}

    ic_type_plots_gpd = gpd.read_file(plots_shp_path)  
    wavelengths = pd.read_csv(os.path.join(stacked_aop_path,  "wavelengths.txt"))['wavelengths'].tolist()
    
    # Filter to only stacked AOP tiles that overlap with ic_type plots
    overlapping_tiles = []
    overlapping_coordinates = []
    for raster_path in stacked_aop_list:
        with rasterio.open(raster_path) as src:
            bounds = src.bounds
            raster_box = box(bounds.left, bounds.bottom, bounds.right, bounds.top)
            if ic_type_plots_gpd.geometry.intersects(raster_box).any():
                overlapping_tiles.append(raster_path)
                overlapping_coordinates.append(f'{int(bounds.left)}_{int(bounds.bottom)}')

    # Identify tiles that still need to be classified
    remaining_tiles = [path for path in overlapping_tiles if extract_coords(os.path.basename(path)) not in classified_coords]

    # Define class labels and corresponding colors
    predicted_classes = ["rock_bare", "grass", "pine", "cedar", "fir", "shrub", "oak"]
    class_colors = ["gray", "green", "blue", "brown", "lightgreen", "orange", "purple"]

    # Create a mapping from unique strings to unique numeric values
    string_to_number = {string: i for i, string in enumerate(predicted_classes)} 
    number_to_color = {i: class_colors[i] for i in range(len(predicted_classes))}

    # Create a colormap using the specific colors
    cmap = plt.cm.colors.ListedColormap(class_colors)
                
    # Only run if we haven't finished classifying all the tiles
    if remaining_tiles:
        rf_model = joblib.load(rf_model_path) 

        # Filter out unwanted wavelengths    
        wavelength_lut = filter_out_wavelengths(wavelengths=wavelengths, layer_names=stacked_aop_layer_names)

        # Classify pixels in each overlapping tile
        print(f"Predicting pixel-level PFTs for {site} {year_inv}, saving to {classified_tiles_dir}")
        for stacked_aop_filename in remaining_tiles:

            # Load NEON tile
            tile_img, tile_gdal_obj = read_raster(stacked_aop_filename)
            [bands, rows, cols] = tile_img.shape

            # Extract relevant info from tile for naming
            xmin = tile_gdal_obj.GetGeoTransform()[0] 
            ymin = tile_gdal_obj.GetGeoTransform()[3] - 1000 
            east_north_string = f"{int(xmin)}_{int(ymin)}"

            # Save classification of tif file in intermediate data in new classified_pixels/ folder
            east_north_classified_tif_path = os.path.join(classified_tiles_dir,
                                                f"classified_{east_north_string}.tif")
                            
            # If classified tile already exists, skip
            if not os.path.exists(east_north_classified_tif_path):
                print(f"Predicting pixel-level PFTs across tile: {east_north_string}")

                # Image vectorization
                tile_vec = image_vectorization_3d(tile_img)

                # Convert the NumPy array to a Pandas DataFrame and assign column names
                tile_df = pd.DataFrame(tile_vec, columns=stacked_aop_layer_names)

                # Features to use in the RF models
                rf_features_noWavelengths = [x for x in stacked_aop_layer_names if not x.isdigit()]
                rf_features_wavelengths = round(wavelength_lut['wavelength']).astype(int).astype(str).tolist()
                rf_features_all = rf_features_noWavelengths + rf_features_wavelengths

                # Filter image to selected bands
                tile_df_features = tile_df[rf_features_all]

                # Remove rows with NaN values
                tile_df_features_noNA = tile_df_features.dropna()
                    
                if pcaInsteadOfWavelengths:
                    # The pipeline preprocessor applies scaler+PCA to spectral columns and a
                    # separate scaler to other columns. It selects by integer index, so columns
                    # must be in the same order as training: wavelengths first (as X1, X2, ...),
                    # then other features. Don't apply scaler/PCA manually — the pipeline does it.
                    rename_map   = {wl: f'X{i+1}' for i, wl in enumerate(rf_features_wavelengths)}
                    spec_renamed = tile_df_features_noNA[rf_features_wavelengths].rename(columns=rename_map)
                    other_cols   = tile_df_features_noNA[rf_features_noWavelengths]
                    features = pd.concat([spec_renamed.reset_index(drop=True),
                                          other_cols.reset_index(drop=True)], axis=1)
                else:
                    features = tile_df_features_noNA
                
                # Predict PFT ID
                features = features.loc[:, ~features.columns.isin(['eastingIDs', 'northingIDs', 'pixelNumber'])]
                features['pred_PFT'] = rf_model.predict(features) 

                # Convert the array of strings to an array of numbers using the mapping
                pred_img_num = np.vectorize(string_to_number.get)(features['pred_PFT'])

                # Replace NaN values in the original DataFrame with predictions
                pred_PFT_withNAs = pd.Series(index=tile_df_features.index, dtype='object')
                pred_PFT_withNAs[tile_df_features_noNA.index] = pred_img_num  # Fill in predictions
                
                # Reshape the vector to image shape
                pred_img = pred_PFT_withNAs.to_numpy().reshape(rows, cols)
                
                # Save to geotiff
                format_img = "GTiff"    
                array2gtiff(stacked_aop_filename, pred_img, format_img, east_north_classified_tif_path)

    # Clip classified TIFFs to the ic_type plots and save them
    # Also plot each one side-by-side with RGB
    for classified_tile_path in glob.glob(os.path.join(classified_tiles_dir, '*.tif')):

        # Check if the tif file overlaps with the polygons
        if any(coord in classified_tile_path for coord in overlapping_coordinates):
            
            with rasterio.open(classified_tile_path) as src:
                print(f"Clipping {classified_tile_path} to overlapping polygons...")

                # Filter to only polygons that overlap this tif
                classified_tile_bounds = src.bounds
                raster_box = box(classified_tile_bounds.left, classified_tile_bounds.bottom, classified_tile_bounds.right, classified_tile_bounds.top)
                overlapping_polygons = ic_type_plots_gpd[ic_type_plots_gpd.geometry.within(raster_box)]
                
                # Clip the raster to the extent of the polygons
                for _, polygon in overlapping_polygons.iterrows():

                    # Create a filename for the clipped raster
                    if "PLOTID" in polygon.index:
                        polygon_name = f"{polygon['PLOTID']}_{int(polygon.geometry.bounds[0])}_{int(polygon.geometry.bounds[1])}"  
                    elif "plotID" in polygon.index:
                        polygon_name = f"{polygon['plotID']}_{int(polygon.geometry.bounds[0])}_{int(polygon.geometry.bounds[1])}"  
                    final_clipped_path = os.path.join(ic_type_path, "clipped_to_plots")
                    os.makedirs(final_clipped_path, exist_ok=True)
                    clipped_plot_tif_path = os.path.join(final_clipped_path, f"{polygon_name}_classified_clipped.tif")
                    
                    if not os.path.exists(clipped_plot_tif_path):

                        # Clip to the current polygon (exact clip for saving)
                        out_image, out_transform = mask.mask(src, [polygon['geometry']], crop=True)

                        with rasterio.open(clipped_plot_tif_path, 'w', driver='GTiff', height=out_image.shape[1],
                                        width=out_image.shape[2], count=1, dtype='int32',
                                        crs=src.crs, transform=out_transform) as dst:
                            dst.write(out_image.astype(rasterio.uint8))

                        # Buffered read for plotting: 3 pixels of context beyond polygon boundary
                        buf_px_w = out_transform.a
                        buf_px_h = abs(out_transform.e)
                        bounds = polygon['geometry'].bounds
                        buf_box = box(bounds[0] - 3 * buf_px_w, bounds[1] - 3 * buf_px_h,
                                      bounds[2] + 3 * buf_px_w, bounds[3] + 3 * buf_px_h)
                        out_image_plot, out_transform_plot = mask.mask(src, [buf_box], crop=True)

                        # Clip RGB + CHM to polygon and plot classified / RGB / CHM
                        for rgb_image_path in glob.glob(os.path.join(data_raw_aop_path, site, year_aop, 'tif', '*0_image.tif')):
                            rgb_filename = os.path.basename(rgb_image_path)
                            rgb_coords = rgb_filename.split('_')[-3:-1]

                            # Check for overlap based on bounding box
                            if (int(rgb_coords[0]) == classified_tile_bounds[0] and int(rgb_coords[1]) == classified_tile_bounds[1]):

                                # Load corresponding RGB image (buffered)
                                with rasterio.open(rgb_image_path) as src_rgb:
                                    rgb_buf_box = box(bounds[0] - 3 * src_rgb.res[1], bounds[1] - 3 * src_rgb.res[0],
                                                      bounds[2] + 3 * src_rgb.res[1], bounds[3] + 3 * src_rgb.res[0])
                                    out_image_rgb, out_transform_rgb = rasterio.mask.mask(src_rgb, [rgb_buf_box], crop=True)

                                # Load corresponding CHM tile (same easting/northing, _CHM.tif suffix)
                                chm_candidates = glob.glob(os.path.join(
                                    data_raw_aop_path, site, year_aop, 'tif',
                                    f'*{rgb_coords[0]}_{rgb_coords[1]}_CHM.tif'))
                                out_image_chm, out_transform_chm = None, None
                                if chm_candidates:
                                    with rasterio.open(chm_candidates[0]) as src_chm:
                                        chm_buf_box = box(bounds[0] - 3 * src_chm.res[1], bounds[1] - 3 * src_chm.res[0],
                                                          bounds[2] + 3 * src_chm.res[1], bounds[3] + 3 * src_chm.res[0])
                                        out_image_chm, out_transform_chm = rasterio.mask.mask(
                                            src_chm, [chm_buf_box], crop=True)

                                # Load training polygons and filter to those overlapping this clip extent
                                training_shp = os.path.join(data_int_path, site, year_inv,
                                                             "training", "ref_labelled_crowns.shp")
                                training_polys = None
                                if os.path.exists(training_shp):
                                    all_training = gpd.read_file(training_shp)
                                    clip_geom = box(*polygon['geometry'].bounds)
                                    training_polys = all_training[all_training.geometry.intersects(clip_geom)]

                                def _overlay_polys(ax, transform, polys, label_col=None):
                                    if polys is None or polys.empty:
                                        return
                                    for _, row in polys.iterrows():
                                        geom = row.geometry
                                        parts = geom.geoms if geom.geom_type == 'MultiPolygon' else [geom]
                                        for part in parts:
                                            xs, ys = zip(*part.exterior.coords)
                                            pixel_rows, pixel_cols = rowcol(transform, xs, ys)
                                            ax.plot(pixel_cols, pixel_rows, 'w-', linewidth=0.8)
                                        if label_col and label_col in polys.columns:
                                            cx, cy = geom.centroid.x, geom.centroid.y
                                            (crow,), (ccol,) = rowcol(transform, [cx], [cy])
                                            ax.text(ccol, crow, str(row[label_col]),
                                                    fontsize=5, color='black', ha='center', va='center',
                                                    bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                                                              edgecolor='none', alpha=0.8))

                                def _draw_plot_outline(ax, transform, geom):
                                    parts = geom.geoms if geom.geom_type == 'MultiPolygon' else [geom]
                                    for part in parts:
                                        xs, ys = zip(*part.exterior.coords)
                                        pixel_rows, pixel_cols = rowcol(transform, xs, ys)
                                        ax.plot(pixel_cols, pixel_rows, 'r-', linewidth=1.5)

                                # Build figure: 3 panels when CHM available, else 2
                                has_chm = out_image_chm is not None
                                n_img_cols = 3 if has_chm else 2
                                fig = plt.figure(figsize=(6 * n_img_cols, 6))
                                gs = gridspec.GridSpec(2, n_img_cols, height_ratios=[10, 1], hspace=0.3)

                                # Order: RGB (col 0), CHM (col 1, if present), Classified (last col)
                                classified_col = 2 if has_chm else 1
                                ax_rgb        = fig.add_subplot(gs[0, 0])
                                ax_classified = fig.add_subplot(gs[0, classified_col])

                                # RGB panel (buffered)
                                ax_rgb.imshow(out_image_rgb.transpose(1, 2, 0))
                                ax_rgb.set_title('RGB')
                                ax_rgb.axis('off')
                                _overlay_polys(ax_rgb, out_transform_rgb, training_polys)
                                _draw_plot_outline(ax_rgb, out_transform_rgb, polygon['geometry'])

                                # CHM panel + colorbar (buffered)
                                if has_chm:
                                    ax_chm = fig.add_subplot(gs[0, 1])
                                    chm_data = out_image_chm.squeeze().astype(float)
                                    chm_data[chm_data <= -9999] = np.nan
                                    im_chm = ax_chm.imshow(chm_data, cmap='viridis')
                                    ax_chm.set_title('CHM')
                                    ax_chm.axis('off')
                                    _overlay_polys(ax_chm, out_transform_chm, training_polys)
                                    _draw_plot_outline(ax_chm, out_transform_chm, polygon['geometry'])
                                    chm_cbar_ax = fig.add_subplot(gs[1, 1])
                                    chm_cbar = plt.colorbar(im_chm, cax=chm_cbar_ax, orientation='horizontal')
                                    chm_cbar.set_label('Canopy height (m)')

                                # Classified panel (buffered)
                                classified_data = out_image_plot.squeeze().astype(float)
                                classified_data[classified_data < 0] = np.nan
                                im_classified = ax_classified.imshow(
                                    np.ma.masked_invalid(classified_data), cmap=cmap, vmin=0, vmax=len(string_to_number) - 1)
                                ax_classified.set_title('Classified')
                                ax_classified.axis('off')
                                _overlay_polys(ax_classified, out_transform_plot, training_polys, label_col='pft')
                                _draw_plot_outline(ax_classified, out_transform_plot, polygon['geometry'])

                                # PFT colorbar under Classified panel
                                cbar_ax = fig.add_subplot(gs[1, classified_col])
                                cbar = plt.colorbar(im_classified, cax=cbar_ax, orientation='horizontal')
                                cbar.set_ticks(range(len(string_to_number)))
                                cbar.set_ticklabels(predicted_classes)
                                cbar.set_label('Predicted PFTs')

                                # Save the comparison plot
                                comparison_plot_path = os.path.join(ic_type_path, "clipped_to_plots", f"{polygon_name}_comparison.png")
                                fig.suptitle(f"{site} {year_aop} — {polygon_name}", fontsize=14, y=1.02)
                                plt.tight_layout()
                                plt.savefig(comparison_plot_path, bbox_inches='tight')
                                plt.close(fig)

    return os.path.join(ic_type_path, "clipped_to_plots")
