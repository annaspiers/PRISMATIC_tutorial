# install.packages("bvls")  # r-bvls incompatible with R 4.4 on conda; use limSolve instead
.this_dir <- dirname(normalizePath(sys.frame(1)$ofile))
library(lidR)
library(sf)
library(randomForest)
library(ggplot2)
library(limSolve)
library(tidyr) #pivot_longer()
library(assertthat)
library(dplyr)
library(rjson)
library(future)
library(purrr)
library(terra)
plan(multisession, workers = future::availableCores()/4 ) # Use all but a few

source(file.path(.this_dir, "inventory_helper.R"))

# Drop-in replacement for bvls() using limSolve::lsei()
bvls <- function(A, b, bl, bu) {
    n <- ncol(A)
    if (is.null(n) || n == 0 || nrow(A) == 0 || length(b) == 0) {
        return(list(x=numeric(0), deviance=0))
    }
    sol <- limSolve::lsei(A=A, B=as.numeric(b),
                          G=rbind(diag(n), -diag(n)),
                          H=c(bl, -bu))
    list(x=sol$X, deviance=sum((A %*% sol$X - b)^2))
}

#ais replace this function below with sy-toan's clip_lidar_to_polygon_lidR
clip_norm_laz_to_shp <- function(site, year_inv, data_int_path, ic_type_path,
                                 plots_shp_path) {
    # Function to clip normalized lidar point cloud to plots being  
    # used to generate initial conditions
    
    laz_clipped_path <- file.path(ic_type_path, "clipped_to_plots")
    
    # Create directory if it doesn't exist                                     
    if (!dir.exists(laz_clipped_path)) {
        dir.create(laz_clipped_path) 
    }
    
    las_ctg <- lidR::readLAScatalog(file.path(data_int_path, site, year_inv,
                                        "normalized_lidar_tiles"))
    opt_output_files(las_ctg) <- file.path(laz_clipped_path, #labelled as NEON tiles are (left, bottom)
                                            "plot_{PLOTID}_{XLEFT}_{YBOTTOM}") #XCENTER}_{YCENTER}")
    opt_laz_compression(las_ctg) <- TRUE 
    plots_sf <- sf::st_read(plots_shp_path)

    # Extract the plots. Returns a list of extracted point clouds in R
    random_plots_ctg <- lidR::clip_roi(las_ctg, plots_sf)
    #ais can I parallelize this? ^
    #ais address differences in tile count between normalized_lidar_tiles, aop_tiles, and this

}


generate_random_plots <- function(aoi_shp_path, ic_type,ic_type_path, n_plots, min_distance, plot_length) { #site, year_aop, year_inv, data_raw_aop_path
                                  
    # Function to generate random points within a flight boundary
    
    random_plots_path <- file.path(ic_type_path, "ic_type_plots.shp")
    
    if (file.exists(random_plots_path)) {
        message(paste0("file with ",n_plots," random plots already generated"))
    } else {
        aop_tiles <- sf::st_read(aoi_shp_path)
        
        # Create or load a polygon
        # Limited by hyperspectral extent with "green" cloud conditions
        flight_bdry_sf <- aop_tiles$geometry
        
        #extent_temp <- read_sf("/Users/AISpiers/Downloads/QGIS_LDRD_StudyArea/NEON_raster_CHM/NEON_SOAP_A01_2019_Boundary.shp")
        #^ais what is this? Also, will need to filter this based on what tiles
        #we filter to based on the 'green' cloud cover flag, etc.
        
        points <- st_sfc(crs=st_crs(flight_bdry_sf))
        plots <- st_sfc(crs=st_crs(flight_bdry_sf))
        n <- 0
        
        while (n < n_plots) {
            point <- st_sample(flight_bdry_sf, size=1)
            
            # Check that point sufficiently far from existing points and polygon edge
            if (n == 0 ||
                as.numeric(min(st_distance(point, points))) >= min_distance) {
                if (as.numeric(min(st_geometry(obj = flight_bdry_sf) %>%
                                   st_cast(to = 'LINESTRING') %>%
                                   st_distance(y = point))) >= min_distance) {
                    
                    # Convert to 20m x 20m
                    plot <- point %>%
                        st_buffer(dist = plot_length/2) %>%
                        st_bbox() %>%
                        sf::st_as_sfc() %>%
                        sf::st_as_sf()
                    
                    plots <- rbind(plots, plot) %>%
                        as.data.frame %>%
                        sf::st_as_sf()
                    st_crs(plots) <- st_crs(flight_bdry_sf)
                    
                    n <- n+1
                }
            }
        }
        
        plots_final <- plots %>%
            tibble::rowid_to_column() %>%
            rename(PLOTID=rowid) %>%
            mutate(subplotID = "1") 
        
        # save as shp
        st_write(plots_final, random_plots_path, delete_layer = TRUE)

        p <- ggplot() +
            geom_sf(data = aoi_shp, color = "black", fill=NA) +  
            geom_sf(data = plots_shp, color = "red", fill="red") +    
            labs(title = paste0('AOI and Plots for initial condition type ')) +
            coord_sf(default_crs = sf::st_crs(32611)) +
            theme_minimal()
        ggsave(file.path(ic_type_path,'plot_coverage_across_AOI.pdf'), plot = p, width = 10, height = 8, dpi = 300)
    }
    
    return(random_plots_path)
}

generate_wall2wall_plots <- function(aoi_shp_path, plot_length,ic_type,ic_type_path) {
    # Function to generate plots wall to wall across an AOI 
        # aoi_shp_path (.shp):  path to AOI shapefile that will have plots gridded across it
        # plot_length (number): dimension of square plots in meters
        # ic_type_path (folder path):  path to save gridded plots 
        
    # Filename to save gridded plots
    gridded_plots_path <- file.path(ic_type_path, "ic_type_plots.shp")
    
    if (file.exists(gridded_plots_path)) {
        message("file with gridded plots already generated")
    } else {
        # Read in AOI shapefile         
        aop_shp <- sf::st_read(aoi_shp_path)

        grid <- sf::st_make_grid(aop_shp, cellsize = plot_length, 
                                 what = "polygons", square = TRUE,
                                 offset=c(floor(ext(aop_shp)[1]/1000)*1000,
                                          floor(ext(aop_shp)[3]/1000)*1000))

        clipped_grid <- sf::st_intersection(grid, aop_shp)

        plots_final <- clipped_grid %>% 
            # Filter to only polygons
            st_sf() %>%
            filter(grepl("POLYGON", st_geometry_type(geometry))) %>%
            # Filter to only plots that are 95% or more of the max area
            mutate(area = as.numeric(st_area(geometry))) %>%
            filter(area/(plot_length^2) > .99) %>% #ais chose this arbitrarily
            data.frame() %>%
            tibble::rowid_to_column() %>%
            dplyr::rename(PLOTID=rowid) %>%
            dplyr::mutate(subplotID = "1") 

        for (i in 1:nrow(plots_final)) {
            plots_final$e_min[i] = min(plots_final$geometry[[i]][[1]][,1])
            plots_final$n_min[i] = min(plots_final$geometry[[i]][[1]][,2])
        }

        # save as shp
        st_write(plots_final, gridded_plots_path, delete_layer = TRUE)
        
        p = ggplot2::ggplot() +
            ggplot2::geom_sf(data = aop_shp, color = "black", alpha = 0.5) +
            ggplot2::geom_sf(data = plots_final$geometry, color = "red", alpha = 0.7) +
            ggplot2::labs(title = paste0('AOI and Plots for initial condition type ', ic_type)) +
            ggplot2::theme_minimal()
        ggplot2::ggsave(file.path(ic_type_path,'patch_coverage_across_AOI.pdf'), 
                        plot = p, width = 10, height = 8, dpi = 300)
        ggplot2::ggsave(file.path(ic_type_path,'patch_coverage_across_AOI.png'), 
                        plot = p, width = 10, height = 8, dpi = 300)
    }
        
    return(gridded_plots_path)    
}




divide_patch_into_cohorts <- function(plots_laz_path, ic_type) {
    
    lad_laz_plot_paths <- list.files(file.path(plots_laz_path), 
                                     pattern="*.laz", full.names = T) %>%
        tools::file_path_sans_ext()
    
    # Initialize dataframe
    by_patch_df <- data.frame()
    #height_step <- 0.5 #m #for height step in lad_csv
    
    for (c in 1:length(lad_laz_plot_paths)) {
        print(c)
        
        # Load leaf area density files
        lad_json <- rjson::fromJSON(file=paste0(lad_laz_plot_paths[c],"_lad.json"))
        lad_csv <- read.csv(paste0(lad_laz_plot_paths[c],"_lad.csv")) 
        
        if ( length(lad_json$layer_height)==0 & nrow(lad_csv)==0 ) {
            # in a previous step, but why?
            
            # In this case, the patch is flat and has no cohorts
            patch_temp = data.frame(z=1.25,
                                    lad=0,
                                    # Then this patch has no canopy, so cohort height is 0
                                    cohort_height = 1,
                                    # take the minimum height (for calculating LAI when there is canopy)
                                    diff_z = 1.25,
                                    patch = ifelse(grepl("central", basename(lad_laz_plot_paths[c])),
                                                   sub(".*/(.*)$", "\\1", lad_laz_plot_paths[c]), #sub(".*/plot_(\\d+)_.*", "\\1", lad_laz_plot_paths[c]),
                                                   ifelse(ic_type=="rs_inv_plots",sub(".*/(.*)$", "\\1", lad_laz_plot_paths[c]),
                                                          sub(".*/plot_(\\d+)_.*","\\1",lad_laz_plot_paths[c]))),
                                    patch_bottom_left = ifelse(grepl("central", basename(lad_laz_plot_paths[c])) | ic_type=="rs_inv_plots",
                                                               NA,sub(".*?(\\d{6}_\\d{7}).*","\\1",lad_laz_plot_paths[c])),
                                    # This patch has no cohorts
                                    cohort_idx = 0) 
            print(paste0("patch ",c," has nrow(lad_csv)==0 & length(lad_json$layer_height)==0"))
            
        } else if ( length(lad_json$layer_height)==0 ) {
            #If lad_json is empty, but lad_csv has values, then we can assume 
            # there is one cohort. Pool all z heights into the tallest cohort
            
            patch_temp = data.frame(z=lad_csv$z,
                                    lad=lad_csv$lad,
                                    cohort_height = max(lad_csv$z),
                                    # take the minimum height (for calculating LAI when there is canopy)
                                    diff_z = 1.25,
                                    patch = ifelse(grepl("central", basename(lad_laz_plot_paths[c])),
                                                   sub(".*/(.*)$", "\\1", lad_laz_plot_paths[c]), #sub(".*/plot_(\\d+)_.*", "\\1", lad_laz_plot_paths[c]),
                                                   ifelse(ic_type=="rs_inv_plots",sub(".*/(.*)$", "\\1", lad_laz_plot_paths[c]),
                                                          sub(".*/plot_(\\d+)_.*","\\1",lad_laz_plot_paths[c]))),
                                    patch_bottom_left = ifelse(grepl("central", basename(lad_laz_plot_paths[c])) | ic_type=="rs_inv_plots",
                                                               NA,sub(".*?(\\d{6}_\\d{7}).*","\\1",lad_laz_plot_paths[c])),
                                    # This patch has no cohorts
                                    cohort_idx = 1) 
            print(paste0("patch ",c," has length(lad_json$layer_height)==0"))
        } else {
            # lad units are m2/m3
            if (is.null(lad_csv$X)) { #ais what is this column X?
                lad_csv <- lad_csv %>%
                    #add lad for rows beyond highest layer_height to highest layer_height
                    mutate(too_tall_lad_sum = sum(lad[z > max(lad_json$layer_height)])) %>%
                    filter(z <= max(lad_json$layer_height)) %>%
                    # Add the lad_sum to the lad of the row with the highest remaining z
                    mutate(lad = if_else(z == max(z), lad + too_tall_lad_sum, lad)) %>%
                    select(-too_tall_lad_sum)
            } else {
                lad_csv <- lad_csv %>%
                    dplyr::select(-c("X")) %>%
                    #add lad for rows beyond highest layer_height to highest layer_height
                    mutate(too_tall_lad_sum = sum(lad[z > max(lad_json$layer_height)])) %>%
                    filter(z <= max(lad_json$layer_height)) %>%
                    # Add the lad_sum to the lad of the row with the highest remaining z
                    mutate(lad = if_else(z == max(z), lad + too_tall_lad_sum, lad)) %>%
                    select(-too_tall_lad_sum)
            }
            
            #ais figure out how to vectorize this rather than for-loop
            for (l in 1:nrow(lad_csv)) {
                # identify the closest layer_height for each row
                # choosing reference height as top height in cohort
                lad_csv$cohort_height[l] <- min(lad_json$layer_height[
                    lad_json$layer_height >= lad_csv$z[l]])
                
                # calculate the difference in height from one row to the next
                # (for calculating LAI)
                lad_csv$diff_z[l] <- ifelse(l==1,lad_csv$z[l],
                                            lad_csv$z[l]-lad_csv$z[l-1])
            }
            
            cohort_idx_df <- lad_csv %>%
                filter(lad>0) %>%
                distinct(cohort_height) %>%
                arrange(desc(cohort_height)) %>%
                dplyr::mutate(cohort_idx = row_number())
            
            # ais when would 'central' be part of filename now?
            if (grepl("central", basename(lad_laz_plot_paths[c]))) { 
                patch_temp <- lad_csv %>%
                    filter(lad>0) %>%
                    # assign patch and cohort index
                    mutate(patch = basename(lad_laz_plot_paths[c]),
                           patch_bottom_left=NA) %>%
                    left_join(cohort_idx_df)
            } else {
                patch_temp <- lad_csv %>%
                    filter(lad>0) %>%
                    # assign patch # and cohort index
                    mutate(patch = ifelse(ic_type=="rs_inv_plots",
                                          sub(".*/(.*)$", "\\1", lad_laz_plot_paths[c]),
                                          sub(".*/plot_(\\d+)_.*","\\1",lad_laz_plot_paths[c])),
                           ,
                           patch_bottom_left=ifelse(ic_type=="rs_inv_plots",NA,
                                                    sub(".*?(\\d{6}_\\d{7}).*", "\\1", 
                                                        lad_laz_plot_paths[c]))) %>%
                    left_join(cohort_idx_df)
            }
        }
        # assign things back into rows of cohort_df
        by_patch_df <- rbind(by_patch_df, patch_temp)
    }
    
    return(by_patch_df)
}



assign_pft_across_cohorts <- function(by_patch_df, allom_params,
                                      perc_pfts_per_patch_df) {
    
    veg_rel_area = perc_pfts_per_patch_df %>% 
        mutate(rel_area=rowSums(select(., where(is.numeric)))) %>%
        dplyr::select(patch, rel_area) %>% distinct()
    
    breakdown_to_pft_df <- by_patch_df %>%
        # Add the rel_area column to adjust LAI
        left_join(veg_rel_area) %>% 
        
        # lai = sum(layer thickness (diff_z) * LAD of layer / av relative vegetated area in patch)
        # lai (layer 1) = stem density (layer 1) * ILA(as a function of z)
        # ILA = Bleaf * SLA , SLA is specific to height and PFT - for now just pft
        
        # assign cohort LAI
        dplyr::mutate(lai = lad*diff_z/rel_area) %>%
        dplyr::group_by(patch,cohort_idx,cohort_height,layer) %>% #ais do I mean layer
        dplyr::reframe(lai_cohort = sum(lai)) %>%
        dplyr::ungroup() %>%
        
        # assign layer LAI
        dplyr::group_by(patch, layer) %>%
        dplyr::mutate(lai_layer = sum(lai_cohort)) %>%
        dplyr::ungroup() %>%
        
        # assign total LAI for patch
        dplyr::group_by(patch) %>%
        dplyr::mutate( lai_patch_T = sum(lai_cohort) ) %>%
        dplyr::ungroup() %>%
        
        # percent that each PFT takes up LAI per patch (drop grass)
        dplyr::left_join(perc_pfts_per_patch_df %>%
                             dplyr::select(-p_g_T) %>%
                             # fill out PFT % so that they sum to 1
                             rowwise() %>%                                         
                             mutate(across(starts_with("p_"),                      
                                           ~ . / sum(c_across(starts_with("p_"))))) %>%
                             ungroup() %>% distinct()) %>%
        
        # known fractions
        #f_o_c = fraction of oaks in layer c
        #f_o_b = fraction of oaks in layer b
        #f_o_a = fraction of oaks in layer a
        mutate(f_s = ifelse(layer=="a", (p_s_T*lai_patch_T)/lai_layer, 0))
    
    # Estimate fraction of each PFT (f_p,f_c,f_o,f_f) for each layer
    message("Estimating PFT fractions across cohorts")
    
    # initialize
    params_est <- data.frame()
    for (p in unique(breakdown_to_pft_df$patch)) {
        
        temp_df <- breakdown_to_pft_df %>%
            dplyr::filter(patch == p)
        
        # Isolate the lai total for each layer
        lai_a <- { vals <- temp_df %>% dplyr::filter(layer=="a") %>% dplyr::distinct(lai_layer) %>% dplyr::pull(lai_layer); if (length(vals) > 0 && !is.na(vals[1])) vals[1] else 0 }
        lai_b <- { vals <- temp_df %>% dplyr::filter(layer=="b") %>% dplyr::distinct(lai_layer) %>% dplyr::pull(lai_layer); if (length(vals) > 0 && !is.na(vals[1])) vals[1] else 0 }
        lai_c <- { vals <- temp_df %>% dplyr::filter(layer=="c") %>% dplyr::distinct(lai_layer) %>% dplyr::pull(lai_layer); if (length(vals) > 0 && !is.na(vals[1])) vals[1] else 0 }
        lai_patch_T <- temp_df %>% dplyr::distinct(lai_patch_T) %>% dplyr::pull(lai_patch_T)
        p_p_T <- { vals <- if("p_p_T" %in% colnames(temp_df)) temp_df %>% dplyr::distinct(p_p_T) %>% dplyr::pull(p_p_T) else numeric(0); if (length(vals)>0 && !is.na(vals[1])) vals[1] else 0 }
        p_c_T <- { vals <- if("p_c_T" %in% colnames(temp_df)) temp_df %>% dplyr::distinct(p_c_T) %>% dplyr::pull(p_c_T) else numeric(0); if (length(vals)>0 && !is.na(vals[1])) vals[1] else 0 }
        p_f_T <- { vals <- if("p_f_T" %in% colnames(temp_df)) temp_df %>% dplyr::distinct(p_f_T) %>% dplyr::pull(p_f_T) else numeric(0); if (length(vals)>0 && !is.na(vals[1])) vals[1] else 0 }
        p_o_T <- { vals <- if("p_o_T" %in% colnames(temp_df)) temp_df %>% dplyr::distinct(p_o_T) %>% dplyr::pull(p_o_T) else numeric(0); if (length(vals)>0 && !is.na(vals[1])) vals[1] else 0 }
        p_s_T <- { vals <- if("p_s_T" %in% colnames(temp_df)) temp_df %>% dplyr::distinct(p_s_T) %>% dplyr::pull(p_s_T) else numeric(0); if (length(vals)>0 && !is.na(vals[1])) vals[1] else 0 }
        
        f_s_vals <- temp_df %>% dplyr::filter(layer=="a") %>% dplyr::distinct(f_s) %>% dplyr::pull(f_s)
        f_s_a <- if (isTRUE(lai_a > 0) && length(f_s_vals) > 0 && !is.na(f_s_vals[1])) f_s_vals[1] else 0
        
        # system of equations
        # 1 - f_s_a = f_c_a + f_p_a + f_f_a + f_o_a  (f_s_a is known)
        #         1 = f_c_b + f_p_b + f_f_b + f_o_b
        #         1 = f_c_c + f_p_c + f_f_c
        # p_c_T*lai_patch_T = f_c_a*lai_a + f_c_b*lai_b + f_c_c*lai_c
        # p_p_T*lai_patch_T = f_p_a*lai_a + f_p_b*lai_b + f_p_c*lai_c
        # p_f_T*lai_patch_T = f_f_a*lai_a + f_f_b*lai_b + f_f_c*lai_c
        # p_o_T*lai_patch_T = f_o_a*lai_a + f_o_b*lai_b 
        
        # also an equation, but not included above
        # p_s_T*lai_patch_T = f_s_a*lai_a + f_s_b*lai_b
        # but since there is no shrub beyond layer a, f_o_b=f_o_c=0
        # p_s_T*lai_patch_T = f_s_a*lai_a
        
        # Write in matrix form
        # ais need to scale this out to include more PFTs too
        # column order: f_c_a,f_c_b,f_c_c,f_p_a,f_p_b,f_p_c,f_f_a,f_f_b,f_f_c,f_o_a,f_o_b
        full_X <- data.frame(matrix(c(1,0,0,1,0,0,1,0,0,1,0, #layer a
                                      0,1,0,0,1,0,0,1,0,0,1, #layer b
                                      0,0,1,0,0,1,0,0,1,0,0, #layer c
                                      lai_a,lai_b,lai_c,0,0,0,0,0,0,0,0,  #cedar
                                      0,0,0,lai_a,lai_b,lai_c,0,0,0,0,0,  #pine
                                      0,0,0,0,0,0,lai_a,lai_b,lai_c,0,0,  #fir
                                      0,0,0,0,0,0,0,0,0,lai_a,lai_b),     #oak
                                    7,11, byrow=TRUE))
        colnames(full_X) <- c("f_c_a", "f_c_b", "f_c_c", 
                              "f_p_a", "f_p_b", "f_p_c", 
                              "f_f_a", "f_f_b", "f_f_c", 
                              "f_o_a", "f_o_b")
        rownames(full_X) <- c("layer a", "layer b", "layer c", 
                              "cedar", "pine", "fir", "oak")
        # Create a vector of rows to exclude based on conditions
        rows_to_exclude <- c(
            ifelse(lai_a == 0, "layer a", ""),
            ifelse(lai_b == 0, "layer b", ""),
            ifelse(lai_c == 0, "layer c", ""),
            ifelse(p_p_T == 0, "pine", ""),
            ifelse(p_c_T == 0, "cedar", ""),
            ifelse(p_f_T == 0, "fir", ""),
            ifelse(p_o_T == 0, "oak", ""))
        
        # Create a vector of column patterns to exclude based on conditions
        columns_to_exclude <- c(
            ifelse(lai_a == 0, "_a$", ""),
            ifelse(lai_b == 0, "_b$", ""),
            ifelse(lai_c == 0, "_c$", ""),
            ifelse(p_p_T == 0, "_p_", ""),
            ifelse(p_c_T == 0, "_c_", ""),
            ifelse(p_f_T == 0, "_f_", ""),
            ifelse(p_o_T == 0, "_o_", ""))        
        columns_to_exclude <- columns_to_exclude[!columns_to_exclude == ""]

        # Filter rows and columns in fewer lines
        X <- full_X %>%
            filter(!rownames(.) %in% rows_to_exclude) %>%
            select(-matches(paste(columns_to_exclude, collapse = "|"))) %>%
            as.matrix()

        if (nrow(X) == 0 || ncol(X) == 0) next

        if (f_s_a > 1) { #minority of cases
            #e.g., "SOAP_021_central"
            # f_s_a should be between 0 and 1, but sometimes it is greater than 1
            # when f_s_a > 1, just set all of layer a to shrub,
            # and split up layer b between the rest of the PFTs 
            f_s_a <- 1
            # which then forces f_p_a = f_c_a = f_f_a = f_o_a = 0, and we need 
            # to recalibrate. Since shrub can only be in layer a, but the 
            # proportion of lai in layer a is less than p_s_T, then we set 
            # p_s_T = lai_a/lai_patch_T,and recalibrate the other PFTs accordingly
            rescal <- (1 - lai_a/lai_patch_T)/(p_p_T + p_c_T + p_f_T + p_o_T)
            if (rescal!="Inf") {
                if ("p_p_T" %in% colnames(temp_df)) { p_p_T <- rescal*p_p_T } 
                if ("p_c_T" %in% colnames(temp_df)) { p_c_T <- rescal*p_c_T } 
                if ("p_f_T" %in% colnames(temp_df)) { p_f_T <- rescal*p_f_T } 
                if ("p_o_T" %in% colnames(temp_df)) { p_o_T <- rescal*p_o_T } 
            }
            # so now the system of equations above should balance out
        } else {
            rescal <- 1 
        }
        
        # see system of equations above - the left hand side
        y <- matrix(c(1 - f_s_a, 1, 1, p_c_T*lai_patch_T, p_p_T*lai_patch_T,
                      p_f_T*lai_patch_T, p_o_T*lai_patch_T), 7, 1)
        # Get rid of equations 
        incl_eqns <- c(   if (lai_a > 0) 1,
                          if (lai_b > 0) 2, if (lai_c > 0) 3,
                          if (p_c_T > 0) 4, if (p_p_T > 0) 5,
                          if (p_f_T > 0) 6, if (p_o_T > 0) 7 )
        y = y[incl_eqns]
        #assertion: ncol(X)==length(y)
        
        # If key > 0, istate is as follows: the last contains the total
        # number of components at their bounds (the bound variables). 
        
        # The absolute values of the first nbound <- tail(istate,1) entries of
        # istate are the indices of these bound components of x. The sign of
        # istate[1:nbound] indicates whether x(abs(istate[1:nbound])) is at
        # its upper or lower bound. istate[1:nbound] is positive if the
        # component is at its upper bound, negative if the component is at
        # its lower bound. istate[(nbound+1):ncol(A)] contain the indices of
        # the components of x that are active (i.e. are expected to lie
        # strictly within their bounds). When key > 0, the routine initially
        # sets the active components to the averages of their upper and
        # lower bounds.
        
        b_l <- rep(0,ncol(X))
        b_u <- rep(1,ncol(X)) 
        sol_temp <- bvls(X, y, b_l, b_u) 
        dev_df <- cbind(as.data.frame(setNames(as.list(sol_temp$x), colnames(X))),
                        data.frame(setNames(as.list(b_l), paste0("b_l_", seq_len(length(b_l)))),
                                   dev=sol_temp$deviance))
        
        i <- 1
        while (max(b_l) < 1) {
            ind_fixed <- which(dev_df[i,1:ncol(X)] == min(dev_df[i,1:ncol(X)]))
            b_l[ind_fixed] <- b_l[ind_fixed] + 0.01
            sol_temp <- bvls(X, y, b_l, b_u) #solve(X,y)
            # Add next row
            dev_df <- rbind(dev_df,c(sol_temp$x,b_l,sol_temp$deviance))
            i <- i+1
        }
        
        params_temp <- cbind(p, f_s_a, dev_df %>%
                                 # Select row with lowest deviance
                                 dplyr::filter(dev==min(dev_df$dev)) %>%
                                 dplyr::slice_head(n=1)  )
        
        params_est <- bind_rows(params_est ,params_temp)
    }

    if (nrow(params_est) == 0 || !"p" %in% colnames(params_est)) {
        warning("No valid patches found for PFT estimation — all patches skipped. Check LAD/PFT data.")
        return(list(character(0), character(0)))
    }

    params_est_f_s <- params_est %>%
        dplyr::select(patch=p, f_s_a) %>%
        mutate(f_s_b = 0) %>%
        tidyr::pivot_longer(cols=c(f_s_a,f_s_b), names_prefix = "f_s_", names_to="layer",
                            values_to="f_s")
    params_est_f_o <- if(any(colnames(params_est) %in% c("f_o_a", "f_o_b"))) params_est %>%
        dplyr::select(patch=p, any_of(c("f_o_a", "f_o_b"))) %>%
        tidyr::pivot_longer(cols=c(f_o_a,f_o_b), names_prefix = "f_o_", names_to="layer",
                            values_to="f_o")
    params_est_f_c <- if(any(colnames(params_est) %in% c("f_c_a", "f_c_b", "f_c_c"))) params_est %>%
        dplyr::select(patch=p, any_of(c("f_c_a", "f_c_b", "f_c_c"))) %>%
        tidyr::pivot_longer(cols=c(f_c_a,f_c_b,f_c_c), names_prefix = "f_c_", names_to="layer",
                            values_to="f_c")
    params_est_f_p <- if(any(colnames(params_est) %in% c("f_p_a", "f_p_b", "f_p_c"))) params_est %>%
        dplyr::select(patch=p, any_of(c("f_p_a", "f_p_b", "f_p_c"))) %>%
        tidyr::pivot_longer(cols=c(f_p_a, f_p_b, f_p_c), names_prefix = "f_p_", names_to="layer",
                            values_to="f_p")
    params_est_f_f <- if(any(colnames(params_est) %in% c("f_f_a", "f_f_b", "f_f_c"))) params_est %>%
        dplyr::select(patch=p, any_of(c("f_f_a", "f_f_b", "f_f_c"))) %>%
        tidyr::pivot_longer(cols=c(f_f_a, f_f_b, f_f_c), names_prefix = "f_f_", names_to="layer",
                            values_to="f_f")
    # Filter out NULL dfs
    valid_dfs <- Filter(Negate(is.null), list(breakdown_to_pft_df,params_est_f_s,params_est_f_c, 
                                              params_est_f_p,params_est_f_f,params_est_f_o))
    
    # # Sanity check: system of equations
    # # 1 = f_o_a + f_c_a + f_p_a and 1 = f_c_b + f_p_b
    # test1 <- reduce(valid_dfs, left_join) %>%
    #     dplyr::mutate(should_be_1 = f_o + f_c + f_f + f_p + f_s) %>%
    #     dplyr::select(patch,layer,should_be_1) %>% distinct() %>%
    #     tidyr::drop_na() %>% data.frame()
    # # right column should be 1's
    # assertthat::assert_that(round(floor(sum(test1$should_be_1)))==nrow(test1),
    #                         msg="pft fractions do not sum to 1 for at least one patch")
    
    # # p_c_T*lai_patch_T = f_c_a*lai_a + f_c_b*lai_b
    # test2 <- breakdown_to_pft_df %>%
    #     dplyr::left_join(params_est_f_c) %>%
    #     dplyr::group_by(patch) %>%
    #     dplyr::reframe(cedar_lai_patch_Total = p_c_T*lai_patch_T) %>% distinct() %>%
    #     dplyr::left_join(breakdown_to_pft_df %>%
    #                          dplyr::left_join(params_est_f_c) %>%
    #                          dplyr::group_by(patch,layer) %>%
    #                          dplyr::reframe(cedar_by_layer = f_c*lai_layer) %>%
    #                          dplyr::distinct() %>%
    #                          dplyr::group_by(patch) %>%
    #                          dplyr::summarize(cedar_tot = sum(cedar_by_layer)) ) %>%
    #     tidyr::drop_na() %>%
    #     dplyr::mutate(diff=cedar_tot - cedar_lai_patch_Total)
    # # right column should be as close to 0 as possible
    # message(paste0("largest LAI difference between allometrically and statistically derived cedar ",max(abs(test2$diff))))
    # # assertthat::assert_that(max(abs(test2$diff)) < 0.1, #ais arbitrarily picked this threshold
    # #                         msg="difference between allometrically and statistically derived cedar is too large for at least one patch")
    
    
    # # p_p_T*lai_patch_T = f_p_a*lai_a + f_p_b*lai_b
    # test3 <- breakdown_to_pft_df %>%
    #     dplyr::left_join(params_est_f_p) %>%
    #     dplyr::group_by(patch) %>%
    #     dplyr::reframe(pine_lai_patch_Total = p_p_T*lai_patch_T) %>% distinct() %>%
    #     dplyr::left_join(breakdown_to_pft_df %>%
    #                          dplyr::left_join(params_est_f_p) %>%
    #                          dplyr::group_by(patch,layer) %>%
    #                          dplyr::reframe(pine_by_layer = f_p*lai_layer) %>%
    #                          dplyr::distinct() %>%
    #                          dplyr::group_by(patch) %>%
    #                          dplyr::summarize(pine_tot = sum(pine_by_layer)) ) %>%
    #     tidyr::drop_na() %>%
    #     dplyr::mutate(diff=pine_tot - pine_lai_patch_Total)
    # # right column should be as close to 0 as possible
    # message(paste0("largest LAI difference between allometrically and statistically derived pine ",max(abs(test3$diff))))
    # assertthat::assert_that(max(abs(test3$diff)) < 0.1, #ais arbitrarily picked this threshold
    #                 msg="difference between allometrically and statistically derived cedar is too large for at least one patch")
    
    # ais are not all on par (diff is not 0), but good enough
    
    patch_all_df <- reduce(valid_dfs, left_join)  %>%
        dplyr::select(-c(lai_layer)) %>%
        tidyr::pivot_longer(cols=any_of(c("f_o","f_c","f_p","f_s","f_f")), names_to="pft", values_to="frac") %>%
        dplyr::mutate(lai = frac*lai_cohort) %>%
        dplyr::filter(lai > 0) %>%
        dplyr::group_by(patch) %>%
        dplyr::mutate(cohort_idx_new = row_number()) %>% 
        dplyr::ungroup() %>%
        dplyr::select(c(patch,pft,cohort_height,cohort_idx=cohort_idx_new,lai,lai_patch_T)) %>%
        dplyr::mutate(pft = case_when(
            pft == "f_p" ~ "pine",
            pft == "f_f" ~ "fir",
            pft == "f_c" ~ "cedar",
            pft == "f_o" ~ "oak",
            pft == "f_s" ~ "shrub") ) %>%
        dplyr::left_join(allom_params) %>%
        dplyr::group_by(patch,lai_patch_T,cohort_idx,cohort_height,pft) %>%
        dplyr::reframe(lai_cohort = lai,
                       dbh = d1 * pmin(cohort_height,Hmax)^d2,
                       agb = a1 * dbh^a2,
                       leaf_biom  = x1 * pmin(cohort_height,Hmax)^x2,
                       n_stemdens = lai_cohort / (leaf_biom  * SLA),
                       ba = n_stemdens*pi/4*dbh^2) %>%
        dplyr::ungroup() %>% distinct() %>% arrange(patch,cohort_idx) %>%
        # keep track of vegetated relative area
        dplyr::left_join(veg_rel_area)
    
    # Sanity check: LAI
    breakdown_to_pft_df %>% group_by(patch) %>%
        dplyr::summarize(lai_patch_Tot_a = sum(lai_cohort)) %>%
        dplyr::left_join(patch_all_df %>% group_by(patch) %>% 
                             dplyr::summarize(lai_patch_Tot_b = sum(lai_cohort)) ) %>%
        tidyr::drop_na() %>%
        dplyr::mutate(diff = lai_patch_Tot_a - lai_patch_Tot_b) %>%
        arrange(desc(diff))
    #ais almost spot on (diff==0) - good enough!
    
    return(patch_all_df)
}




extract_perc_pfts <- function(p, classified_tifs, predicted_classes) {
    # Find TIFFs that intersect this plot
    overlapping_tifs <- classified_tifs[sapply(classified_tifs, function(tif_path) {
        r <- terra::rast(tif_path)
        rast_bbox <- sf::st_as_sfc(st_bbox(r))
        any(st_intersects(p, rast_bbox, sparse = FALSE))
    })]
    
    if (length(overlapping_tifs) == 0) {
        warning(paste("No overlapping TIFFs found for plot", i))
        return(NULL)
    }
    
    # Merge rasters if needed
    classified_rast <- if (length(overlapping_tifs) == 1) {
        terra::rast(overlapping_tifs[1])
    } else {
        do.call(terra::merge, lapply(overlapping_tifs, terra::rast))
    }
    
    # Clip raster to plot
    p_vect <- vect(p)
    clipped_rast <- terra::crop(classified_rast, p_vect) %>% terra::mask(p_vect)
    
    # Skip if clipped raster is empty
    if (is.na(minmax(clipped_rast)[1])) {
        warning(paste("Tif is empty where plot overlaps. Skipping plot", i))
        return(NULL)
    }
    
    # Get frequency and convert to percent
    freq_table <- terra::freq(clipped_rast) %>%
        as.data.frame() %>%
        filter(!is.na(value))
    
    total_pixels <- sum(freq_table$count)
    freq_table$percent <- freq_table$count / total_pixels * 100
    freq_table$class <- predicted_classes[freq_table$value + 1]  # 0-indexed
    freq_table$plotID <- p %>% data.frame() %>% select(any_of(c("plotID", "PLOTID"))) %>% pull()
    freq_table$subplotID <- p %>% data.frame() %>% select(any_of(c("subplotID"))) %>% pull()
    
    return(freq_table[, c("plotID", "subplotID", "value", "class", "percent")])
}



generate_pcss_field <- function(site, year, data_int_path, biomass_path, ic_type, ic_type_path, 
                          multisite=FALSE, plots_shp_path=NULL, classified_plot_PFTs_path=NULL,
                          plots_laz_path=NULL) {
# Generate cohort (.css) and patch (.pss) files for FATES initialization for field inventory data

    # allom_path <- file.path(data_int_path,"AllometricParameters_SOAPeqns_byMarcos_ais.csv") # ais TODO: find this file
    allom_path <- file.path(data_int_path,"AllometricParameters_Fates_PB+AHB - SOAP equations.csv") # workaround
    warning("Using workaround allometric parameters file. Replace with AllometricParameters_SOAPeqns_byMarcos_ais.csv when found.")
    allom_params <- read.csv(allom_path)
    
    # Load cleaned individual-level data
    veg <- read.csv(biomass_path)%>% 
        rename(SLA_sytoan = SLA) %>% #to differentiate from SLA from allometry file
        
        # If an individual is multistem, select only the largest stem (usually oak)
        dplyr::group_by(individualID) %>% 
        dplyr::slice(which.max(used_diameter)) %>% 
        ungroup() %>%
        
        #Assign pft
        # ais how does this work for non-NEON inventory data?
        mutate(pft = match_species_to_pft(growthForm, taxonID)) %>% 
        # Filter out pft 'other'
        filter(pft != "other") %>% #ais do I want to do this?
        # Assign unknown plants to PFT - ais do this - all of them were shrubs anyway
        #filter(scientificName != "Unknown plant") %>% #20 plants
        left_join(allom_params) %>%
        dplyr::mutate(time = year, 
                      patch = ifelse(is.na(subplotID),paste0(plotID,"_central"), paste0(plotID,"_",subplotID)),
                      index = row_number(),
                      dbh = used_diameter,
                      n = individualStemNumberDensity, 
                      agb = a1 * dbh^a2, #ais use biomass function
                      leaf_biom = x1 * pmin(height,Hmax)^x2, #ais where is height coming from?
                      lai = n * SLA * leaf_biom)
    # lai (layer 1) = stem density (layer 1) * ILA(as a funciton of z)
    # ILA = Bleaf * SLA , SLA is specific to height and PFT - for now just pft
    
    # Visualize / sanity check
    # by patch and pft
        
    veg %>% group_by(patch) %>% dplyr::reframe(sum_lai = sum(lai,na.rm=T)) %>% 
        left_join( veg %>% group_by(patch, pft) %>%
                       dplyr::reframe(sum_lai = sum(lai, na.rm=T)) %>%
                        group_by(patch) %>%
                       slice_max(sum_lai) %>% dplyr::select(-c(sum_lai)) %>%
                       rename(pft_dom=pft) ) %>%
        ggplot(aes(sum_lai,fill=pft_dom)) + 
        geom_histogram() + xlab("LAI (m2/m2, inventory)") +
        ylab("Patch count") 
    ggsave(file.path(ic_type_path,paste0(ic_type,"_lai.png"))) 
    veg %>% group_by(patch) %>% dplyr::reframe(sum_n = sum(n,na.rm=T)) %>% 
        left_join(  veg %>% group_by(patch, pft) %>%
                        dplyr::reframe(sum_lai = sum(lai, na.rm=T)) %>%
                        group_by(patch) %>%
                        slice_max(sum_lai) %>% dplyr::select(-c(sum_lai)) %>%
                        rename(pft_dom=pft) ) %>%
        ggplot(aes(sum_n,fill=pft_dom)) + geom_histogram() + 
        xlab("Stem density (1/m2, inventory)") + ylab("Patch count") 
    ggsave(file.path(ic_type_path,paste0(ic_type,"_stemdens.png"))) 
    veg %>% group_by(patch) %>% dplyr::reframe(sum_lb = sum(leaf_biom*n,na.rm=T)) %>% 
        left_join(  veg %>% group_by(patch, pft) %>%
                        dplyr::reframe(sum_lai = sum(lai, na.rm=T)) %>%
                        group_by(patch) %>%
                        slice_max(sum_lai) %>% dplyr::select(-c(sum_lai)) %>%
                        rename(pft_dom=pft) ) %>%
        ggplot(aes(sum_lb,fill=pft_dom)) + geom_histogram() + 
        xlab("Leaf biomass (kg/m2, inventory)") +
        ylab("Patch count") 
    ggsave(file.path(ic_type_path,paste0(ic_type,"_leafbiom.png"))) 
    veg %>% group_by(patch) %>% dplyr::reframe(sum_agb = sum(agb*n,na.rm=T)) %>% 
        left_join( veg %>% group_by(patch, pft) %>%
                       dplyr::reframe(sum_lai = sum(lai, na.rm=T)) %>%
                       group_by(patch) %>%
                       slice_max(sum_lai) %>% dplyr::select(-c(sum_lai)) %>%
                       rename(pft_dom=pft) ) %>%
        ggplot(aes(sum_agb,fill=pft_dom)) + geom_histogram() + 
        xlab("Aboveground biomass (kg/m2, inventory)") +
        ylab("Patch count") 
    ggsave(file.path(ic_type_path,paste0(ic_type,"_agb.png"))) 
    veg %>% group_by(patch) %>% dplyr::reframe(sum_ba = sum(individualBasalArea*n,na.rm=T)) %>% 
        left_join(  veg %>% group_by(patch, pft) %>%
                        dplyr::reframe(sum_lai = sum(lai, na.rm=T)) %>%
                        group_by(patch) %>%
                        slice_max(sum_lai) %>% dplyr::select(-c(sum_lai)) %>%
                        rename(pft_dom=pft) ) %>%
        ggplot(aes(sum_ba,fill=pft_dom)) + geom_histogram() + 
        xlab("Basal area (cm2/m2) by patch (Inventory)") +
        ylab("Patch count") 
    ggsave(file.path(ic_type_path,paste0(ic_type,"_ba.png"))) 
    
    # Save structural attributes locally
    plot_structure_summ <- veg %>% 
        dplyr::group_by(patch) %>%
        dplyr::summarize(BasalArea=sum(individualBasalArea),
                         StemNumberDensity=sum(individualStemNumberDensity),
                         biomass=sum(individualBiomass_kg),
                         lai=sum(lai))      
    readr::write_csv(x=plot_structure_summ,file=file.path(ic_type_path,"plot_structure_summ.csv"))
    
    # Create final files
    cohort_df_initial <- veg %>%
        mutate(pft = case_when(
            pft == "pine" ~ 1,
            pft == "cedar" ~ 2,
            pft == "fir" ~ 3,
            pft == "shrub" ~ 4,
            pft == "oak" ~ 5 )) %>%
        dplyr::select(time, patch, index, dbh, height, pft, n)  
    
    # Patch file
    patch_df_initial <- veg %>%
        distinct(patch) %>% arrange(patch) %>%
        mutate(area = 1/(nrow(veg %>% distinct(patch))))
    
    # Validate css and pss
    pcss_list <- validate_pss_css(year, cohort_df_initial, patch_df_initial, ic_type, allom_params)
    cohort_df_final <- pcss_list[[1]]
    patch_df_final <- pcss_list[[2]]
    
    # Save final initial conditions --------------------------------------------
    # wiki: https://github.com/EDmodel/ED2/wiki/Initial-conditions#files-types-and-formats-for-nlied_init_mode6
    
    #no multisite for field application - creating for Danielle who needs it only for RS ic's
    
    cohort_path <- file.path(ic_type_path,"ic.css")
    patch_path  <- file.path(ic_type_path,"ic.pss")
    
    readr::write_delim(x=cohort_df_final,file=cohort_path,delim=" ",append=FALSE,quote="none")
    readr::write_delim(x=patch_df_final,file=patch_path,delim=" ",append=FALSE,quote="none")
    
    return(c(cohort_path, patch_path))
}



generate_pcss_rs <- function(site, year, data_int_path, biomass_path, ic_type, ic_type_path, 
                                        plots_shp_path=NULL, classified_plot_PFTs_path=NULL, 
                                        plots_laz_path=NULL, multisite=FALSE) {
    # Generate cohort (.css) and patch (.pss) files for FATES initialization for one of the rs_*_... ic_types 
    
    # allom_path <- file.path(data_int_path,"AllometricParameters_SOAPeqns_byMarcos_ais.csv") # ais TODO: find this file
    allom_path <- file.path(data_int_path,"AllometricParameters_Fates_PB+AHB - SOAP equations.csv") # workaround
    warning("Using workaround allometric parameters file. Replace with AllometricParameters_SOAPeqns_byMarcos_ais.csv when found.")
    allom_params <- read.csv(allom_path)
    # ais how to access this data reproducibly? I just manually put this file in that folder
    
    # load plots shp
    plots_sf <- sf::st_read(plots_shp_path)
    
    # load plots classified as PFT by percentage
    classified_tifs <- list.files(path = classified_plot_PFTs_path, pattern = "\\.tif$", full.names = TRUE)
    
    predicted_classes = c("rock_bare", "grass", "pine", "cedar", "fir", "shrub", "oak")
    
    # Extract % pft per plot
    summary_df <- purrr::map(seq_len(nrow(plots_sf)), function(i) {
        extract_perc_pfts(plots_sf[i,], classified_tifs, predicted_classes)
    }) %>% purrr::list_rbind() #takes a good bit of time
    
    # Summarize cohort structure in each patch as 
    by_patch_df = divide_patch_into_cohorts(plots_laz_path, ic_type) 
    
    perc_pfts_per_patch <- summary_df %>%
        dplyr::rename(pft = class) %>%
        dplyr::mutate(percent = percent/100,
                      patch = paste0(plotID,"_",subplotID)) %>%
                      #patch=as.character(patch)) %>% 
        
        # Remove NA classes and recalibrate percent values
        dplyr::filter(!is.na(pft)) %>%
        dplyr::group_by(patch) %>%
        dplyr::mutate(percent = percent / sum(percent)) %>%
        dplyr::ungroup() %>%
        
        dplyr::select(-c(value)) %>% 
        tidyr::pivot_wider(names_from = pft, values_from = percent) %>%
        
        # Keep patch spatial coordinates
        left_join(by_patch_df %>% select(patch,patch_bottom_left)) %>%
        distinct()
    perc_pfts_per_patch <- perc_pfts_per_patch %>%
        dplyr::mutate(dplyr::across(where(is.numeric), ~tidyr::replace_na(., 0)))
    perc_pfts_per_patch$patch = as.character(perc_pfts_per_patch$patch)
    
    # Check for plot-level compatibility between cohorts and PFT distribution. 
    combined_struct_comp_df <- by_patch_df %>%
        # percent of each pft per patch
        dplyr::left_join(perc_pfts_per_patch, by="patch") 
    
    mismatch_patches <- data.frame()
    for (p in unique(combined_struct_comp_df$patch)) {
        temp_df <- combined_struct_comp_df %>%
            dplyr::filter(patch == p)
        
        # If a pft is classified in plot, but no cohort exists that is shorter than 
        # the hmax for that pft, then make a cohort at the hmax height 
        # this will happen mostly for shrub and oak
        for (pft in allom_params$pft) {
            if (pft %in% colnames(temp_df) & length(unique(temp_df[[pft]]))>0 & 
                    min(temp_df$cohort_height) > allom_params$Hmax[allom_params$pft==pft]) {
                mismatch_patches <- rbind(mismatch_patches, c(p, pft, allom_params$Hmax[allom_params$pft==pft]))
            }
        }
    }    
    if (nrow(mismatch_patches) > 0) {
        colnames(mismatch_patches) = c("patch","pft","new_cohort_height")
        mismatch_patches$new_cohort_height = as.numeric(mismatch_patches$new_cohort_height)
    } else {
        mismatch_patches <- data.frame(patch=character(), pft=character(), new_cohort_height=numeric())
    }
    
    # add new cohorts into combined strucutre (LAD) and composition (PFT) df
    by_patch_df <- by_patch_df %>%
        # Join mismatch_patches to by_patch_df based on the patch column
        left_join(mismatch_patches, by = "patch") %>%
        # Update cohort_idx for rows where cohort_height < new_cohort
        mutate(cohort_idx = ifelse(!is.na(new_cohort_height) & z < new_cohort_height, 
                                   cohort_idx + 1, cohort_idx),
               cohort_height = ifelse(!is.na(new_cohort_height) & z < new_cohort_height, 
                                      new_cohort_height, cohort_height)) %>%
        select(-c(pft,new_cohort_height))
    # ais work for another time - fix LAD plots (plot_diagnostic) to reflect these updates
    
    # Assign layers based on number of distinct Hmax's
    hmax_values = sort(unique(allom_params$Hmax))
    by_patch_df$layer <- base::cut(by_patch_df$cohort_height,
                                   breaks = c(0,hmax_values,100), #add a cap that's arbitrarily large 
                                   # Assign any cohorts that are taller than hmax to tallest layer 
                                   labels = c(letters[1:length(hmax_values)],letters[length(hmax_values)]))
    
    # Split patches into barren (rock_bare) vs vegetated (PFTs + grass)
    patches_veg_vs_barren_df <- perc_pfts_per_patch %>%
        # Add new rows where rock_bare > 0 with `barren` appended to the patch value
        rbind(perc_pfts_per_patch %>% #dplyr::bind_rows(perc_pfts_per_patch %>%
                             dplyr::filter(rock_bare > 0) %>%
                             dplyr::mutate(patch = paste0(patch, "_barren")))
    barren_df <- patches_veg_vs_barren_df %>%
        dplyr::filter(grepl("barren$", patch)) %>%   # Filter rows where patch ends with "barren"
        dplyr::select(patch,rock_bare,patch_bottom_left)            # Keep only `patch` and `rock_bare` columns
    vegetated_df <- patches_veg_vs_barren_df %>%
        dplyr::filter(!grepl("barren$", patch)) %>%  # Filter rows where patch does not end with "barren"
        dplyr::select(-rock_bare) %>%           # Drop the `rock_bare` column
        rename_with( ~ ifelse(.x %in% predicted_classes, paste0("p_", substr(.x, 1, 1), "_T"), .x) )
    
    # Assign PFTs across cohorts using only vegetated patches
    cohorts_all_df <- assign_pft_across_cohorts(by_patch_df, allom_params, 
                                              vegetated_df) 
    
    # Sanity check plots by patch and pft
    cohorts_all_df %>% dplyr::group_by(patch) %>% 
        # dplyr::reframe(sum_lai = sum(lai_cohort)) %>% 
        dplyr::reframe(lai_per_m2_patch = lai_patch_T ) %>% 
        ggplot(aes(lai_per_m2_patch)) + 
        geom_histogram() + xlab("LAI (m2/m2, lidar)") +#LAI should be less than 10 for the most part
        ylab("Patch count") 
    ggsave(file.path(ic_type_path,paste0(ic_type,"_lai.png"))) 
    cohorts_all_df %>% dplyr::group_by(patch) %>% #,pft) %>% 
        dplyr::reframe(sum_n = sum(n_stemdens)) %>% 
        ggplot(aes(sum_n)) + #,fill=pft)) +  
        geom_histogram() + xlab("Stem density (1/m2, lidar)") +
        ylab("Patch count") 
    ggsave(file.path(ic_type_path,paste0(ic_type,"_stemdens.png"))) 
    cohorts_all_df %>% dplyr::group_by(patch) %>% #,pft) %>% 
        dplyr::reframe(sum_lb = sum(leaf_biom*n_stemdens)) %>% 
        ggplot(aes(sum_lb)) + #,fill=pft)) +  
        geom_histogram() + xlab("Leaf biomass (kg/m2, lidar)") +
        ylab("Patch count") 
    ggsave(file.path(ic_type_path,paste0(ic_type,"_leafbiom.png"))) 
    cohorts_all_df %>% dplyr::group_by(patch) %>% #,pft) %>% 
        dplyr::reframe(sum_agb = sum(agb*n_stemdens)) %>% 
        ggplot(aes(sum_agb)) + #,fill=pft)) +  
        geom_histogram() + xlab("Aboveground biomass (kg/m2, lidar)") +
        ylab("Patch count") 
    ggsave(file.path(ic_type_path,paste0(ic_type,"_agb.png")))  
    cohorts_all_df %>% dplyr::group_by(patch) %>% #,pft) %>% 
        dplyr::reframe(sum_ba = sum(ba)) %>% 
        ggplot(aes(sum_ba)) + #,fill=pft)) +  
        geom_histogram() + xlab("Basal area (cm2/m2, lidar)") +
        ylab("Patch count") 
    ggsave(file.path(ic_type_path,paste0(ic_type,"_ba.png"))) 
    
    # Save structural attributes locally
    veg_plot_structure_summ <- cohorts_all_df %>% 
        group_by(patch) %>%
        summarize(n_stemdens_sum=sum(n_stemdens),
                  agb_sum=sum(agb)) %>%
        left_join(by_patch_df %>%
                      dplyr::mutate(lai_cohort = lad*diff_z) %>%
                      dplyr::group_by(patch) %>%
                      dplyr::summarize(lai = sum(lai_cohort))) %>%
        select(patch,n_stemdens_sum,agb_sum,lai)#??basalarea,
    readr::write_csv(x=veg_plot_structure_summ,file=file.path(ic_type_path,"plot_structure_summ.csv"))
    
    # Cohort file - no cohorts in barren file
    cohort_df_initial <- cohorts_all_df %>%
        dplyr::select(-c(lai_cohort,agb,leaf_biom,rel_area)) %>%
        dplyr::rename(index = cohort_idx,
                      height= cohort_height,
                      n = n_stemdens) %>%
        dplyr::mutate(pft = case_when(
            pft == "pine" ~ 1,
            pft == "cedar" ~ 2,
            pft == "fir" ~ 3,
            pft == "shrub" ~ 4,
            pft == "oak" ~ 5 )) 
    
    # Patch file
    patch_area_df = cohorts_all_df %>%
        select(patch,rel_area) %>%
        distinct() %>%
        bind_rows(barren_df %>% rename(rel_area=rock_bare))
    patch_df_initial <- data.frame( #each row is one subplot, or patch
        # Include both barren and vegetated patch
        patch = sort(c(unique(cohort_df_initial$patch),barren_df$patch))) %>%  #ais find realistic value
        left_join(patch_area_df) %>%
        dplyr::mutate(area = rel_area/n())
    
    # Validate css and pss
    pcss_list <- validate_pss_css(year, cohort_df_initial, patch_df_initial, ic_type, allom_params)
    cohort_df_final <- pcss_list[[1]]
    patch_df_final <- pcss_list[[2]]
    
    # Save final initial conditions --------------------------------------------
    # wiki: https://github.com/EDmodel/ED2/wiki/Initial-conditions#files-types-and-
    # formats-for-nlied_init_mode6
    # FATES initial conditions are generated the same way as ED2
    
    if (multisite){  # so separate each patch into its own site
        
        if (!dir.exists(file.path(ic_type_path,"multisite"))) {
            dir.create(file.path(ic_type_path,"multisite")) }
        
        for (p in unique(vegetated_df$patch)){
            # Filter to specific site/patch
            cohort_df_temp = cohort_df_final %>% filter(patch==p) %>%
                # add patch_bottom_left. only whole number (vegetated) patches in cohort df
                left_join(patches_veg_vs_barren_df %>% select(patch, patch_bottom_left) %>% distinct())
            patch_df_temp  = patch_df_final %>% filter(patch %in% c(p,paste0(p,"_barren"))) %>%
                # add patch_bottom_left. only whole number (vegetated) patches in cohort df
                left_join(patches_veg_vs_barren_df %>% 
                              filter(patch %in% c(p,paste0(p,"_barren"))) %>%
                              # Group the barren and vegetated patches into one 'site'
                              # Think about Danielle's use case. Both of these patches
                              # take up the same spatial area
                              #mutate(patch = sub("_barren$", "", patch)) %>% 
                              select(patch, patch_bottom_left)  %>% distinct())
            
            # Site/patch file names
            if (ic_type=="rs_inv_plots") {
                cohort_path <- file.path(ic_type_path,"multisite",
                                         paste0("ic_",unique(patch_df_temp$patch),".css"))
                patch_path  <- file.path(ic_type_path,"multisite",
                                         paste0("ic_",unique(patch_df_temp$patch),".pss"))
            } else {
                cohort_path <- file.path(ic_type_path,"multisite",
                                         paste0("ic_",unique(patch_df_temp$patch_bottom_left),".css"))
                                                             #"_",unique(patch_df_temp$patch),".css"))
                patch_path  <- file.path(ic_type_path,"multisite",
                                         paste0("ic_",unique(patch_df_temp$patch_bottom_left),".pss"))
                                                             #"_",unique(patch_df_temp$patch),".pss"))
            }
            
            # Save patch/cohort file for each site
            readr::write_delim(x=cohort_df_temp %>% select(-patch_bottom_left),
                               file=cohort_path,delim=" ",append=FALSE,quote="none")
            readr::write_delim(x=patch_df_temp %>% select(-patch_bottom_left),
                               file=patch_path,delim=" ",append=FALSE,quote="none")
        }
        
    } else { # Then single site, so save as-is
        
        cohort_path <- file.path(ic_type_path,"ic.css")
        patch_path  <- file.path(ic_type_path,"ic.pss")
        
        readr::write_delim(x=cohort_df_final,file=cohort_path,delim=" ",append=FALSE,quote="none")
        readr::write_delim(x=patch_df_final,file=patch_path,delim=" ",append=FALSE,quote="none")
    }
    return(c(cohort_path, patch_path))
}



validate_pss_css <- function(year, i_css, i_pss, ic_type, allom_params){
    
    initial_list = 
        tibble::tribble( ~desc               ,              ~in_suffix         ,  ~out_suffix       ,     ~year,              ~colour  , ~f_net_def
                         , "Inventory (Field)" ,            "ic_inv"           , "inventory_neon"  ,      as.numeric(year), "#CC8829", 1.0
                         , "Inventory (RS)"    ,            "ic_rs_over_inv"   , "rsinventory_neon",      as.numeric(year), "#47B2AA", 1.0
                         , "Entire site random (RS)",       "ic_rs_over_random", "prismatic_random_neon", as.numeric(year), "#AC5CE5", 1.0 # 2.018285
                         , "Tower footprint (RS)"  ,        "ic_rs_tower_ftpt",  "prismatic_twr_neon"  ,  as.numeric(year), "#F5E940", 1.0 # 2.018285
                         , "Entire site wall to wall (RS)", "ic_rs_wall2wall",   "prismatic_wall_neon"  , as.numeric(year), "#6A9A55", 1.0 # 2.018285
        )
    pft_lookup = allom_params  %>%
        dplyr::mutate(fates_id = case_when(pft=="pine" ~ 1, pft=="cedar" ~ 2,
                                           pft=="fir" ~ 3, pft=="shrub" ~ 4,
                                           pft=="oak" ~ 5, .default = NA),
                      colour = case_when(pft=="pine" ~ "#66CCAE", pft=="cedar" ~ "#1B9E77",
                                         pft=="fir" ~ "#095941", pft=="shrub" ~ "#D95F02",
                                         pft=="oak" ~ "#7570B3", .default = "other"),
                      desc = case_when(pft=="pine" ~ "Pine", pft=="cedar" ~ "Cedar",
                                       pft=="fir" ~ "Fir", pft=="shrub" ~ "Shrub",
                                       pft=="oak" ~ "Oak", .default = "other"))
    
    #   Loop through output types, and fix files.
    i = ifelse(ic_type == "field_inv_plots", 1,
               ifelse(ic_type == "rs_inv_plots", 2,
                      ifelse(ic_type == "rs_random_plots", 3, 
                             ifelse(ic_type == "rs_tower_ftpt", 4, 
                                    ifelse(ic_type == "rs_wall2wall", 5, 0 )))))
    
    #   Load information.
    io_desc   = initial_list$desc      [i]
    io_year   = initial_list$year      [i]
    f_net_def = initial_list$f_net_def [i]
    message(" + Process data for ",tolower(io_desc)," initialisation.")    
    
    #   Validate files
    message("   - Validate files.")
    is_pcss_match = all(i_css$patch %in% i_pss$patch)
    if (! is_pcss_match){
        invalid_patch = i_css$patch[! (i_css$patch %in% i_pss$patch)]
        
        message("---~---"                                                                 )
        message(" FATAL ERROR!"                                                           )
        message("---~---"                                                                 )
        message("  The following patches exist in cohort files but not in the patch file:")
        for (p in seq_along(invalid_patch)){
            message(" + ",invalid_patch[p],".")
        }#end for (p in seq_along(invalid_patch))
        message("---~---"                                                                 )
        stop(" All patches in cohort file must appear in the patch file.")
    }#end if (! is_pcss_match)
    
    
    #   Validate and standardise PFTs.
    message("   - Validate PFTs.")
    if (is.character(i_css$pft)){
        is_pft_fine = all(i_css$pft %in% pft_lookup$as_name)
        
        #   Report invalid PFTs.
        if (! is_pft_fine){
            invalid_pft = sort(unique(i_css$pft[! (i_css$pft %in% pft_lookup$as_name)]))
            
            message("---~---"                                                           )
            message(" FATAL ERROR!"                                                     )
            message("---~---"                                                           )
            message("  The following PFTs exist in cohort files but are not recognised:")
            for (p in seq_along(invalid_pft)){
                message(" + ",invalid_pft[p],".")
            }#end for (p in seq_along(invalid_pft))
            message("---~---"                                                           )
            stop(" All PFTs in cohort file must match a value in \"pft_lookup\"."    )
        }
        
        #   Substitute names with the FATES id
        i_css$pft = pft_lookup$fates_id[match(i_css$pft,pft_lookup$as_name)]
        
    }else{
        is_pft_fine = all(i_css$pft %in% pft_lookup$fates_id)
        if (! is_pft_fine){
            invalid_pft = sort(unique(i_css$pft[! (i_css$pft %in% pft_lookup$fates_id)]))
            
            message("---~---"                                                           )
            message(" FATAL ERROR!"                                                     )
            message("---~---"                                                           )
            message("  The following PFTs exist in cohort files but are not recognised:")
            for (p in seq_along(invalid_pft)){
                message(" + ",invalid_pft[p],".")
            }#end for (p in seq_along(invalid_pft))
            message("---~---"                                                           )
            stop(" All PFTs in cohort file must match a value in \"pft_lookup\"."    )
        }#end if (! is_pft_fine)
    }#end if (is.character(icss$pft))
    
    
    #   Validate and standardise PFTs.
    message("   - Validate cohorts.")
    if (any(i_css$n == 0.)){
        #---~---
        #   Find invalid.
        #---~---
        sel_invalid = which(i_css$n == 0.)
        i_invalid   = i_css[sel_invalid,,drop=""]
        #---~---
        
        
        #---~---
        #   Report invalid PFTs.
        #---~---
        message ("---~---"                                                            )
        message (" FATAL ERROR!"                                                      )
        message ("---~---"                                                            )
        message ("   There are cohorts with zero stem density, and this is not valid.")
        message (" The first few lines are printed here. Check variable \"i_invalid\"")
        message (" for the full list."                                                )
        message ("---~---"                                                            )
        print(i_invalid)
        message ("---~---"                                                            )
        stop(" Cohort cannot have zero density."                                   )
    }
    
    #   Make sure that the patch area adds to 1
    message("   - Standardise patch area.")
    i_pss$area = i_pss$area / sum(i_pss$area)
    
    #   Standardise names...
    message("   - Standardise names.")
    if (! "nplant" %in% names(i_css)) i_css = i_css %>% rename( nplant = n     )
    if ("index"    %in% names(i_css)) i_css = i_css %>% rename( cohort = index )
    
    #   For patches and cohorts, we use either the original name or the hexadecimal code.
    # if (is.numeric(i_pss$patch )) i_pss$patch  = sprintf("0x%3.3X",i_pss$patch )
    # if (is.numeric(i_css$patch )) i_css$patch  = sprintf("0x%3.3X",i_css$patch )
    if (is.numeric(i_pss$patch )) i_pss$patch  = as.character(i_pss$patch)
    if (is.numeric(i_css$patch )) i_css$patch  = as.character(i_css$patch)
    if (is.numeric(i_css$cohort)) i_css$cohort = sprintf("0x%3.3X",i_css$cohort)
    
    #   Apply correction factor if needed.
    i_css$nplant = i_css$nplant * f_net_def
    
    
    #   Discard patches with tiny populations.
    no_css = i_css %>% filter(i_css$nplant <  1.e-8)
    i_css  = i_css %>% filter(i_css$nplant >= 1.e-8)
    
    #   Sort patches and cohorts
    i_pss = i_pss %>% arrange(patch)
    i_css = i_css %>% arrange(patch,cohort, pft, -dbh)
    
    #   Create output patch and cohort structures.
    message("   - Create output structures.")
    o_pss = tibble::tibble( time  = sprintf("%4.4i"  , io_year    )
                            , patch = sprintf("%s"     , i_pss$patch)
                            , trk   = sprintf("%5i"    , 2L         )
                            , age   = sprintf("%6.1f"  , 0.         )
                            , area  = sprintf("%16.14f", i_pss$area )
                            , water = sprintf("%5i    ", 0L         )
                            , fsc   = sprintf("%10.5f" , 0.         )
                            , stsc  = sprintf("%10.5f" , 0.         )
                            , stsl  = sprintf("%10.5f" , 0.         )
                            , ssc   = sprintf("%10.5f" , 0.         )
                            , psc   = sprintf("%10.5f" , 0.         )
                            , msn   = sprintf("%10.5f" , 0.         )
                            , fsn   = sprintf("%10.5f" , 0.         )
    )#end tibble::tibble
    
    #   Create output cohort structure.
    o_css = tibble::tibble( time   = sprintf("%4.4i"  , io_year     )
                            , patch  = sprintf("%s"     , i_css$patch )
                            , cohort = sprintf("%s"     , i_css$cohort)
                            , dbh    = sprintf("%9.3f"  , i_css$dbh   )
                            , height = sprintf("%9.3f"  , 0) #i_css$height) since EITHER dbh OR H must be 0 
                            , pft    = sprintf("%5i"    , i_css$pft   )
                            , nplant = sprintf("%16.10f", i_css$nplant)
                            , bdead  = sprintf("%9.3f"  , 0.          )
                            , balive = sprintf("%9.3f"  , 0.          )
                            , avgRg  = sprintf("%9.3f"  , 0.          )
    )#end tibble:tibble    
    
    #   Write the output files
    message("   - Write formatted files.")
    return(list(o_css, o_pss))
}