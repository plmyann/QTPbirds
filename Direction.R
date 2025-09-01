## Environmental Factors Analysis for Bird Migration
## 
## This script analyzes bird migration patterns in the QTP,
## focusing on environmental factors like wind, precipitation, temperature, elevation,
## and vegetation. It processes tracking data, downloads wind data using a custom function 
## (since rWind's data source is inaccessible), performs random forest modeling to assess 
## factor importance, and projects future migration routes using species distribution 
## modeling (SDM) with MaxEnt.
## 
## Key Changes for Integration:
## - Added custom functions: wind_dl_custom, wind_mean_custom, wind2raster_custom, and 
##   custom_flow.dispersion to replace deprecated rWind functionality.
## - Renamed columns in custom wind functions to match expected structure (lon, lat, u, v).
## - Adjusted wind direction calculation to match standard meteorological convention 
##   (direction from which wind is blowing).
## - Used custom functions in place of rWind's wind.dl_2, wind.mean, wind2raster, and 
##   flow.dispersion.
## - Ensured compatibility with gdistance by converting terra rasters to raster objects 
##   where necessary.
## - Added required libraries like fields for image.plot and raster for transition layers.
## 
## Key Steps:
## 1. Load necessary libraries.
## 2. Define custom functions for wind data handling.
## 3. Read and process shapefiles for bird tracking data.
## 4. Calculate migration directions for spring and autumn.
## 5. Download and process wind data for migration seasons using custom functions.
## 6. Generate visualizations (GIFs and PDFs) of wind patterns and least-cost paths.
## 7. Extract environmental variables (precipitation, temperature, elevation, EVI).
## 8. Prepare data for modeling and perform random forest analysis.
## 9. Visualize factor importance and correlations.
## 10. Project breeding and wintering sites using MaxEnt SDM.
## 11. Analyze distribution shifts, directions, and overlaps between scenarios (with/without Tibetan Plateau).
## 12. Project migration routes based on environmental factors.
## 
## Notes:
## - Paths are hardcoded; adjust as needed for your local environment.
## - Custom wind download uses rerddap from NOAA GFS dataset; ensure dates are within available range.
## - Custom flow.dispersion implements directed conductance based on tailwind support for passive dispersal.
## - Error handling is included in loops to skip problematic species or data.
## - Visualizations are saved to the "Outputs" directory.
## - If historical wind data for 2019 is not available in the dataset, adjust dates or use alternative sources.

# 1. Load Required Libraries
# Note: These libraries handle spatial data, data manipulation, visualization, wind data, 
# machine learning, and species distribution modeling.
library(rgdal)       # For reading shapefiles (legacy; consider switching to sf/terra).
library(dplyr)       # For data manipulation.
library(ggplot2)     # For plotting.
library(lubridate)   # For date handling.
library(dplyr)       # Already loaded, but ensured.
library(terra)       # Modern spatial raster/vector handling (replaces raster package).
library(geosphere)   # For geographic calculations like bearings and distances.
library(magick)      # For creating GIFs from images.
library(gdistance)   # For least-cost path calculations.
library(spatstat)    # For spatial point pattern analysis (e.g., density kernels).
library(randomForest)# For random forest modeling.
library(caret)       # For model training and evaluation.
library(rsample)     # For data splitting (train/test).
library(reshape2)    # For melting data frames for plotting.
library(rmaxent)     # For MaxEnt species distribution modeling.
library(dismo)       # For ecological niche modeling tools.
library(rasterVis)   # For advanced raster visualization.
library(viridis)     # For color palettes.
library(maps)        # For basic world maps.
library(maptools)    # For converting maps to spatial objects.
library(sp)          # Legacy spatial classes (used with spatstat).
library(hrbrthemes)  # Additional ggplot themes.
library(tidyr)       # For tidying data.
library(whitebox)    # Assumed for scale_fill_whitebox_c; install if needed.
library(rerddap)     # For downloading wind data via ERDDAP server.
library(raster)      # For compatibility with gdistance (transition layers).
library(fields)      # For image.plot in wind visualizations.
library(rworldmap)   # For world maps (must be after raster to avoid conflicts).
library(plyr)        # For round_any in later sections.

# 2. Custom Functions to Replace Deprecated rWind Functionality
# Note: These functions download wind data from NOAA GFS via rerddap, compute means, 
# convert to rasters, and compute flow dispersion (conductance) for least-cost paths.
# Wind direction is calculated as the direction from which the wind is blowing.

# Custom wind download (replacement for wind.dl_2)
wind_dl_custom <- function(times, lon1, lon2, lat1, lat2, type = c("read-data", "csv"), trace = TRUE) {
  type <- match.arg(type)
  dsid <- "ncep_global_best"  # Live NOAA GFS 0.25° winds
  info_ds <- info(dsid)
  
  out_list <- vector("list", length(times))
  
  for (i in seq_along(times)) {
    t <- as.POSIXct(times[i], tz = "UTC")
    if (trace) {
      message("Downloading wind for ", format(t, "%Y-%m-%d %H:%M:%S UTC"))
    }
    
    resp <- griddap(info_ds,
                    time = c(t, t),
                    latitude = c(lat1, lat2),
                    longitude = c(lon1, lon2),
                    fields = c("ugrd10m", "vgrd10m"),
                    fmt = "csv")
    
    df <- resp$data %>%
      mutate(
        time = t,
        dir = (atan2(vgrd10m, ugrd10m) * 180 / pi + 180) %% 360,  # Direction from which wind blows (matches rWind).
        speed = sqrt(ugrd10m^2 + vgrd10m^2)
      ) %>%
      rename(lon = longitude, lat = latitude, u = ugrd10m, v = vgrd10m) %>%
      select(time, lat, lon, u, v, dir, speed)
    
    class(df) <- c("rWind", "data.frame")
    out_list[[i]] <- df
    
    if (type == "csv") {
      fn <- sprintf("wind_%s.csv", format(t, "%Y_%m_%d_%H"))
      write.csv(df, fn, row.names = FALSE)
    }
  }
  
  return(out_list)
}

# Custom mean wind (replacement for wind.mean)
wind_mean_custom <- function(wind_list) {
  if (!is.list(wind_list) || length(wind_list) == 0) {
    stop("wind_list must be a non-empty list of rWind data frames")
  }
  
  all_data <- do.call(rbind, wind_list)
  
  mean_data <- all_data %>%
    group_by(lat, lon) %>%
    summarise(
      u = mean(u, na.rm = TRUE),
      v = mean(v, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      dir = (atan2(v, u) * 180 / pi + 180) %% 360,
      speed = sqrt(u^2 + v^2),
      time = NA
    ) %>%
    select(time, lat, lon, u, v, dir, speed)
  
  class(mean_data) <- c("rWind", "data.frame")
  return(mean_data)
}

# Custom conversion to raster (replacement for wind2raster)
wind2raster_custom <- function(wind_df) {
  if (!inherits(wind_df, "data.frame")) {
    stop("wind_df must be a data.frame")
  }
  
  # Create extent and resolution
  ext <- terra::ext(range(wind_df$lon), range(wind_df$lat))
  unique_lon <- sort(unique(wind_df$lon))
  unique_lat <- sort(unique(wind_df$lat))
  res <- c(diff(unique_lon)[1], diff(unique_lat)[1])
  
  # Create rasters
  dir_rast <- rast(ext, resolution = res, crs = "EPSG:4326")
  speed_rast <- rast(ext, resolution = res, crs = "EPSG:4326")
  
  # Order data by lat, lon (terra fills row-major, increasing lat/lon assumed)
  ord <- order(wind_df$lat, wind_df$lon)
  dir_rast[] <- wind_df$dir[ord]
  speed_rast[] <- wind_df$speed[ord]
  
  # Return as list to match original structure
  list(direction = dir_rast, speed = speed_rast)
}

# Custom flow dispersion (replacement for flow.dispersion)
# Note: Implements directed conductance for passive dispersal based on tailwind support.
# Uses gdistance with raster conversions for compatibility.
custom_flow.dispersion <- function(layers, type = "passive", output = "transitionLayer") {
  if (!is.list(layers)) layers <- list(layers)
  
  res <- lapply(layers, function(y) {
    # Convert terra to raster for gdistance
    speed_r <- raster(y$speed)
    dir_r <- raster(y$direction)
    
    # Create isotropic transition based on average speed
    tr <- transition(speed_r, transitionFunction = function(a, b) (a + b) / 2, directions = 8)
    
    # Get adjacent cells
    adj <- adjacent(tr, cells = 1:ncell(tr), directions = 8)
    
    # Coordinates for bearing calculation
    coords <- xyFromCell(tr, 1:ncell(tr))
    from.coords <- coords[adj[, 1], ]
    to.coords <- coords[adj[, 2], ]
    
    # Direction of movement from cell to cell
    move.dir <- bearing(from.coords, to.coords)
    
    # Wind direction at from cell (average could be used, but from is common)
    wind.dir <- extract(dir_r, adj[, 1])
    
    # Cosine of angle difference for tailwind (positive if supporting)
    angle.diff <- (wind.dir - move.dir) %% 360
    angle.diff[angle.diff > 180] <- angle.diff[angle.diff > 180] - 360
    angle.diff <- abs(angle.diff)
    tailwind <- cos(angle.diff * pi / 180)
    tailwind[tailwind < 0] <- 0
    
    # Average speed between cells
    avg.speed <- (extract(speed_r, adj[, 1]) + extract(speed_r, adj[, 2])) / 2
    
    # Conductance for passive: tailwind support * average speed
    cond <- tailwind * avg.speed
    cond[is.na(cond)] <- 0
    
    # Set transition matrix (note: directed, so not symmetric)
    tr@transitionMatrix <- Matrix::sparseMatrix(i = adj[, 1], j = adj[, 2], x = cond, dims = c(ncell(tr), ncell(tr)))
    
    tr
  })
  
  if (length(res) == 1) res <- res[[1]]
  res
}

# 3. Read and Process Bird Tracking Shapefiles
# Note: Loads shapefiles for multiple species, extracts coordinates, and prepares data frame.

allspfile <- list.files("D:/Papers/WEMigration2022/Datainput/Tracking data/Shapefiles/SingleSp", 
                        pattern = ".shp$", full.names = TRUE)
allspfile <- allspfile[c(1:6, 9:10)]  # Select specific files.

allsp <- lapply(allspfile, terra::vect)  # Read as SpatVector objects.

# Extract base filenames for naming.
allspfile_names <- list.files("D:/Papers/WEMigration2022/Datainput/Tracking data/Shapefiles/SingleSp", 
                              pattern = ".shp$", full.names = FALSE)
allspfile_names <- sub(".shp", "", allspfile_names[c(1:6, 9:10)])

names(allsp) <- allspfile_names  # Name the list elements.

# Extract coordinates for each species.
allsp.coords <- lapply(1:length(allsp), function(sp) {
  spname <- names(allsp[sp])  # Get species name.
  sp_vect <- allsp[[sp]]      # Get SpatVector.
  
  # Extract coordinates and attributes.
  crds <- data.frame(lon = terra::crds(sp_vect)[, 1], 
                     lat = terra::crds(sp_vect)[, 2], 
                     sp = spname, 
                     ID = sp_vect$Id, 
                     UL = sp_vect$UL)
  
  return(crds)
})

allsp.coords <- bind_rows(allsp.coords)  # Combine into single data frame.

# Sort, remove duplicates, and handle NA in UL.
allsp.coords <- allsp.coords[with(allsp.coords, order(sp, ID, -lon)), ]
allsp.coords <- distinct(allsp.coords)
allsp.coords$UL[is.na(allsp.coords$UL)] <- 0  # Set NA UL to 0 (assuming upper/lower route indicator).

# Note: Commented lines suggest previous cleaning of species names; uncomment if needed.
# allsp.coords$sp <- gsub("_spring", "", allsp.coords$sp)
# allsp.coords$sp <- gsub("_autumn", "", allsp.coords$sp)
# allsp.coords$sp <- gsub(" stopover", "", allsp.coords$sp)

# 4. Calculate Migration Directions
# Note: Uses geosphere::bearingRhumb to calculate rhumb line bearings (constant direction).
# Autumn: Sorted decreasing longitude; Spring: Sorted increasing longitude.

## Autumn directions
allsp.coords <- allsp.coords %>% 
  group_by(sp, ID) %>% 
  mutate(dir_aut = geosphere::bearingRhumb(
    matrix(c(lead(lon, default = lon[length(lon)]), lead(lat, default = lat[length(lat)])), ncol = 2),
    matrix(c(lon, lat), ncol = 2)
  ))

## Spring directions (re-sort for increasing longitude)
allsp.coords <- allsp.coords[with(allsp.coords, order(sp, ID, lon)), ]
allsp.coords <- allsp.coords %>% 
  group_by(sp, ID) %>% 
  mutate(dir_spr = geosphere::bearingRhumb(
    matrix(c(lead(lon, default = lon[length(lon)]), lead(lat, default = lat[length(lat)])), ncol = 2),
    matrix(c(lon, lat), ncol = 2)
  ))

# Define species lists for spring and autumn.
spring <- unique(allsp.coords$sp)[c(1:2, 4:8)]
autumn <- unique(allsp.coords$sp)[c(1:3, 5:8)]
allsp <- list(spring, autumn)  # List of species by season.

# 5. Prepare Migration Dates and Download Wind Data
# Note: Defines date sequences for spring and autumn migrations in 2019 (with 2-day intervals).
# Commented lines suggest options for multi-year data; currently limited to 2019.
# Note: If 2019 data is not available in GFS archive, adjust dates to current or available range.
dt_spring <- seq(ymd_hms("2019-02-15 12:00:00"), ymd_hms("2019-05-30 12:00:00"), by = "2 days")
dt_autumn <- seq(ymd_hms("2019-07-15 12:00:00"), ymd_hms("2019-11-30 12:00:00"), by = "2 days")
dtall <- list(dt_spring, dt_autumn)

# Read species-specific migration dates.
migratedate <- read.csv("D:/Papers/WEMigration2022/Datainput/migrationdate.csv")
migratedate[, 2:ncol(migratedate)] <- apply(migratedate[, 2:ncol(migratedate)], 2, function(x) paste0(x, "-2019"))
migratedate[, 2:ncol(migratedate)] <- lapply(migratedate[, 2:ncol(migratedate)], 
                                             function(x) as.POSIXct(x, format = "%d-%B-%Y", origin = "1970-01-01", tz = "GMT"))

# 6.1 Process Wind Data and Generate Visualizations
# Note: Loops over seasons (1=spring, 2=autumn).
# Downloads wind data using custom function, creates raster layers with custom wind2raster_custom, 
# plots speed/direction, and generates species-specific GIFs.
# Computes least-cost paths using custom flow dispersion for wind corridors.
windall <- lapply(1:2, function(s) {  # s=1 for spring, s=2 for autumn.
  tryCatch({
    # Download wind data for the season using custom function.
    ww <- wind_dl_custom(dtall[[s]], -20, 179, -20, 75)
    
    # Convert to raster layers using custom function.
    layers <- lapply(ww, wind2raster_custom)
    names(layers) <- format(dtall[[s]], "%Y-%m-%d %H:%M:%S")  # Set names for dates.
    alldates <- dtall[[s]]  # Use input dates.
    
    # Plot wind speed and direction for each layer.
    id <- 0
    for (i in 1:length(layers)) {
      id <- sprintf("%03d", i)
      
      # Wind speed plot.
      png(paste0("D:/Papers/WEMigration2022/Outputs/asia_", s, "_", id, ".png"), width = 1000, height = 600, pointsize = 18)
      fields::image.plot(layers[[i]]$speed, col = bpy.colors(1000), zlim = c(0, 15), 
                         main = sub(" 12:00:00", "", format(ww[[i]]$time[1])))
      terra::lines(getMap(resolution = "low"), lwd = 3)
      dev.off()
      
      # Wind direction plot.
      png(paste0("D:/Papers/WEMigration2022/Outputs/adir_", s, "_", id, ".png"), width = 1250, height = 600, pointsize = 18)
      plot(layers[[i]]$direction, main = sub(" 12:00:00", "", format(ww[[i]]$time[1])))
      terra::lines(getMap(resolution = "low"), lwd = 3)
      dev.off()
    }
    
    # Species-specific processing: GIFs and least-cost paths.
    spwind <- lapply(unique(migratedate$sp), function(sp) {
      tryCatch({
        # Get migration dates for the species and season.
        startdate <- if_else(s == 1, migratedate$Spring_start[migratedate$sp == sp], 
                             migratedate$Autumn_start[migratedate$sp == sp])
        enddate <- if_else(s == 1, migratedate$Spring_end[migratedate$sp == sp], 
                           migratedate$Autumn_end[migratedate$sp == sp])
        
        # Find layer indices for the species' migration period.
        difstart <- alldates - startdate
        difend <- alldates - enddate
        layer.index <- which.min(replace(difstart, difstart < 0, NA)):which.min(replace(difend, difend < 0, NA))
        sp.layers <- layers[layer.index]
        
        # Create GIF for wind speed.
        list.files(path = "D:/Papers/WEMigration2022/Outputs/", pattern = paste0('asia_', s, '.*.png'), full.names = TRUE)[layer.index] %>% 
          image_read() %>% image_join() %>% image_animate(fps = 2) %>% 
          image_write(paste0("D:/Papers/WEMigration2022/Outputs/dynamic_", sp, "_", s, ".gif"))
        
        # Create GIF for wind direction.
        list.files(path = "D:/Papers/WEMigration2022/Outputs/", pattern = paste0('adir_', s, '.*.png'), full.names = TRUE)[layer.index] %>% 
          image_read() %>% image_join() %>% image_animate(fps = 2) %>% 
          image_write(paste0("D:/Papers/WEMigration2022/Outputs/dy_dir_", sp, "_", s, ".gif"))
        
        # Compute wind conductance and least-cost paths using custom function.
        Conductance <- custom_flow.dispersion(sp.layers, type = "passive", output = "transitionLayer")
        
        sa <- c("_spring", "_autumn")
        spp <- if_else(sp %in% allsp[[s]], sp, paste0(sp, sa[s]))
        
        cat(spp, "\n")
        
        loc <- allsp.coords[allsp.coords$sp == spp, ]
        loc <- loc[with(loc, order(-lon)), ]
        loc <- loc[, 1:2]  # Keep only lon, lat.
        loc <- as.matrix(loc)
        
        # Compute cost distances and shortest paths.
        cost_list <- array(NA_real_, dim = c(nrow(loc), nrow(loc), length(sp.layers)))
        paths <- vector("list", length(sp.layers))
        
        a <- ifelse(s == 1, 1, nrow(loc))  # Start index for season.
        b <- ifelse(s == 1, nrow(loc), 1)  # End index for season.
        
        for (i in 1:length(sp.layers)) {
          cost_list[, , i] <- costDistance(Conductance[[i]], loc)
          if (costDistance(Conductance[[i]], loc[a, ], loc[b, ]) != Inf) {
            paths[[i]] <- shortestPath(Conductance[[i]], loc[a, ], loc[b, ], output = "SpatialLines")
          }
        }
        
        paths_clean <- paths[!sapply(paths, is.null)]
        
        if (length(paths_clean) > 0) {
          paths_merged <- do.call(rbind, paths_clean)
          proj4string(paths_merged) <- CRS(as.character(NA))
          
          # Compute density kernel of paths.
          paths_psp <- as(paths_merged, "psp")
          lines_kernel <- density(paths_psp, sigma = 1, dimyx = c(350, 410))
          kernel <- raster(lines_kernel)
          kernel <- extend(kernel, c(-20, 179, -20, 75))
          kernel[is.na(kernel)] <- minValue(kernel)
          
          kernel1 <- maxValue(kernel) - kernel  # Inverted for cost.
          
          # Plot least-cost path density.
          acol <- colorRampPalette(c("grey40", "darkblue", "red2", "orange", "yellow", "white"))
          pdf(paste0("D:/Papers/WEMigration2022/Outputs/", spp, "_", s, ".pdf"))
          plot(kernel, col = acol(1000), zlim = c(-0.01, 6), main = "Least cost paths density")
          lines(getMap(resolution = "high"), lwd = 2)
          dev.off()
          
          pdf(paste0("D:/Papers/WEMigration2022/Outputs/cost_", spp, "_", s, ".pdf"))
          plot(kernel1, col = acol(1000), zlim = c(-0.01, 6), main = "Least cost paths density")
          lines(getMap(resolution = "high"), lwd = 2)
          dev.off()
          
          # Extract wind values for locations.
          STCoords <- SpatialPoints(loc, proj4string = CRS("+proj=longlat +datum=WGS84"))
          loc <- as.data.frame(loc)
          loc$wind_ld <- extract(kernel, STCoords)
          loc$wind_cost <- extract(kernel1, STCoords)
          loc$sp <- spp
          
          return(loc)
        }
      }, error = function(e) { cat("ERROR :", conditionMessage(e), "\n") })
    })
    
    spwind <- bind_rows(spwind)
    
    # Compute and plot average wind using custom functions.
    wind_average <- wind_mean_custom(ww)
    average <- wind2raster_custom(wind_average)
    
    png(paste0("D:/Papers/WEMigration2022/Outputs/speed_", s, ".png"), width = 1250, height = 700, pointsize = 18)
    plot(average$speed)
    dev.off()
    
    png(paste0("D:/Papers/WEMigration2022/Outputs/dir_", s, ".png"), width = 1250, height = 700, pointsize = 18)
    plot(average$direction, col = bpy.colors(1000))
    terra::lines(getMap(resolution = "low"), lwd = 3)
    dev.off()
    
    return(spwind)
  }, error = function(e) { cat("ERROR :", conditionMessage(e), "\n") })
})

# Save workspace after wind processing (optional checkpoint).
# save.image("D:/Papers/WEMigration2022/Outputs/1114.RData")

# 6.2 Merge Wind Data with Coordinates and Extract Other Environmental Factors
# Note: Merges wind data for spring and autumn, then extracts global environmental rasters.
sprdt1 <- merge(allsp.coords, windall[[1]], by = c("lon", "lat", "sp"), all.x = TRUE)
colnames(sprdt1)[ncol(sprdt1)] <- "Wind"
sprdt1$sa <- "Spring"

audt1 <- merge(allsp.coords, windall[[2]], by = c("lon", "lat", "sp"), all.x = TRUE)
colnames(audt1)[ncol(audt1)] <- "Wind"
audt1$sa <- "Autumn"

alldt <- rbind(sprdt1, audt1)

# Create SpatialPoints for extraction.
STCoords <- SpatialPoints(cbind(alldt$lon, alldt$lat), proj4string = CRS("+proj=longlat +datum=WGS84"))

# Extract precipitation and temperature (WorldClim data).
prec <- raster::raster("D:/Papers/WEMigration2022/Datainput/wc2.1_30s_bio/wc2.1_30s_bio_12.tif")  # Annual precipitation.
tprt <- raster::raster("D:/Papers/WEMigration2022/Datainput/wc2.1_30s_bio/wc2.1_30s_bio_1.tif")  # Annual mean temperature.
alldt$prec <- raster::extract(prec, STCoords)
alldt$tprt <- raster::extract(tprt, STCoords)

# Extract elevation (mean from multiple tiles).
elevation_files <- list.files("D:/Papers/WEMigration2022/Datainput/DEM/Tiff", pattern = "tif$", full.names = TRUE)
ele <- lapply(elevation_files, function(x) raster::extract(raster::raster(x), STCoords))
ele <- bind_cols(ele)
alldt$ele <- rowMeans(ele, na.rm = TRUE)

# Extract EVI (Enhanced Vegetation Index) from MODIS files.
evifiles <- list.files("D:/Papers/WEMigration2022/Datainput/MOD13C2", pattern = "hdf$", full.names = TRUE)
STCoords1 <- STCoords@coords
STCoords1[, 1][STCoords1[, 1] < 0] <- STCoords1[, 1][STCoords1[, 1] < 0] + 360  # Adjust longitudes for global wrap.

evi <- lapply(evifiles, function(x) terra::extract(terra::rast(x)[[1]], STCoords1))
evi <- bind_cols(evi)
evi[is.na(evi)] <- 0
alldt$evi <- rowMeans(evi, na.rm = TRUE) * 0.0001  # Scale EVI values.

# Alternative vegetation from NetCDF (overwrites EVI if used).
tp.veg <- rast("D:/Papers/WEMigration2022/Datainput/NoTP/pi.nc", subds = "VEG")
tp <- rast("D:/Papers/WEMigration2022/Datainput/NoTP/pi.nc")
ext(tp.veg) <- ext(tp)
alldt$evi <- 100 - terra::extract(tp.veg[[1]], STCoords1)$`VEG_lsmpft=1`  # Inverted vegetation cover.

# Note: Adjust negative directions if needed (commented).
# alldt$dir_spr[alldt$dir_spr < 0 & !is.na(alldt$dir_spr)] <- alldt$dir_spr[alldt$dir_spr < 0 & !is.na(alldt$dir_spr)] + 360

# 7. Prepare Data for Modeling
# Note: Reorganize data, assign directions, remove NAs, and define migration stages based on longitude.
alldt <- alldt[with(alldt, order(sp, ID, -lon)), ]

sprdt <- alldt[alldt$sa == "Spring", ]
autdt <- alldt[alldt$sa == "Autumn", ]

sprdt$dir <- sprdt$dir_spr
autdt$dir <- autdt$dir_aut

sprdt$dir_aut <- sprdt$dir_spr <- sprdt$wind_ld <- NULL
autdt$dir_aut <- autdt$dir_spr <- autdt$wind_ld <- NULL

sa.dt <- rbind(sprdt, autdt)
sa.dt <- na.omit(sa.dt)
sa.dt$stage <- 0  # Default stage: Overall.

sa.dt.temp <- sa.dt  # Backup for overall.

# Assign stages based on longitude regions.
sa.dt$stage[sa.dt$lon >= 105] <- 1  # East TP.
sa.dt$stage[sa.dt$lon < 105 & sa.dt$lon >= 73] <- 2  # TP.
sa.dt$stage[sa.dt$lon < 73] <- 3  # West TP.

sa.dt <- rbind(sa.dt.temp, sa.dt)  # Combine overall and staged data.

# Save workspace (optional).
save.image("D:/Papers/WEMigration2022/Outputs/1125.RData")

# 8. Random Forest Modeling for Factor Importance
# Note: Loops over stages (0=overall, 1-3=regions), seasons, and routes (commented UL=upper/lower).
# Trains RF on direction ~ environmental factors, computes importance and nRMSE.
impt.stage <- lapply(0:3, function(stage) {
  sdt <- sa.dt[sa.dt$stage == stage, ]
  
  # Note: UL (upper/lower route) loop is commented; currently processes all together.
  # impt <- lapply(0:1, function(r) { ... })
  
  rdt <- sdt  # All routes (r=1 as placeholder).
  r <- 1
  
  autspr <- lapply(c("Spring", "Autumn"), function(s) {
    rsdt <- rdt[rdt$sa == s, ]
    
    set.seed(r)  # For reproducibility.
    rsdt <- initial_split(rsdt)  # 75/25 train/test split.
    rsdt.train <- training(rsdt)
    rsdt.test <- testing(rsdt)
    
    # Train random forest.
    lm_RF <- randomForest(dir ~ prec + tprt + ele + evi + Wind, 
                          data = rsdt.train, importance = TRUE, proximity = TRUE)
    
    # Predict on test set and compute normalized RMSE.
    lm_RF_test <- predict(lm_RF, rsdt.test)
    RMSE_lm <- postResample(lm_RF_test, rsdt.test$dir)[1] / mean(rsdt.test$dir)
    
    # Extract importance (%IncMSE).
    importance_lm <- data.frame(t(importance(lm_RF, type = 1)))
    importance_lm$nRMSE <- RMSE_lm
    importance_lm$season <- s
    
    return(importance_lm)
  })
  
  autspr <- bind_rows(autspr)
  autspr$UL <- r
  autspr$stage <- stage
  
  return(autspr)
})

impt <- bind_rows(impt.stage)

# 9. Visualize Factor Importance
# Note: Melts data for ggplot, renames variables, and plots bar chart by stage and season.
implot <- reshape2::melt(impt, id = c("season", "stage"), measure = c("prec", "tprt", "Wind", "evi", "ele"))

# Note: UL plotting commented.
# implot$UL[implot$UL == 0] <- "Upper route"
# implot$UL[implot$UL == 1] <- "Lower route"

implot$stage[implot$stage == 0] <- "Overall"
implot$stage[implot$stage == 1] <- "East TP"
implot$stage[implot$stage == 2] <- "TP"
implot$stage[implot$stage == 3] <- "West TP"
implot$stage <- factor(implot$stage, levels = c("Overall", "East TP", "TP", "West TP"))

implot <- implot[implot$variable != "nRMSE", ]
implot$variable <- as.character(implot$variable)
implot$variable[implot$variable == "evi"] <- "Vegetation"
implot$variable[implot$variable == "ele"] <- "Elevation"
implot$variable[implot$variable == "tprt"] <- "Temperature"
implot$variable[implot$variable == "prec"] <- "Precipitation"

impt.plt <- ggplot(data = implot, aes(x = variable, y = value)) +
  geom_col(aes(colour = season, fill = season)) +
  scale_fill_manual(values = c("Autumn" = "#f3be77", "Spring" = "#77be96")) +
  scale_colour_manual(values = c("Autumn" = "#faddb2", "Spring" = "#aad1cc")) +
  facet_grid(stage ~ season) +  # Note: + UL if using routes.
  coord_flip() +
  labs(y = "", x = "") +
  theme_bw() +
  theme(legend.position = "none", strip.background = element_rect(fill = "grey96"),
        panel.grid.major = element_blank())

ggsave("D:/Papers/WEMigration2022/Outputs/importance.png", height = 5, width = 9, dpi = 300, plot = impt.plt)

# 10. Visualize Correlations Between Factors and Direction
# Note: Melts data, plots smoothed GLM lines by stage and season.
corplot <- reshape2::melt(sa.dt, id = c("stage", "sa", "dir"), 
                          measure = c("prec", "tprt", "Wind", "evi", "ele"))

# Note: UL renamed if used.
# corplot$UL[corplot$UL == 0] <- "Upper route"
# corplot$UL[corplot$UL == 1] <- "Lower route"

corplot$stage[corplot$stage == 0] <- "Overall"
corplot$stage[corplot$stage == 1] <- "East TP"
corplot$stage[corplot$stage == 2] <- "TP"
corplot$stage[corplot$stage == 3] <- "West TP"
corplot$stage <- factor(corplot$stage, levels = c("Overall", "East TP", "TP", "West TP"))

corplot$variable <- as.character(corplot$variable)
corplot$variable[corplot$variable == "evi"] <- "EVI"
corplot$variable[corplot$variable == "ele"] <- "Elevation"
corplot$variable[corplot$variable == "tprt"] <- "Temperature"
corplot$variable[corplot$variable == "prec"] <- "Precipitation"

for (i in unique(corplot$variable)) {
  subplot <- corplot[corplot$variable == i, ]
  
  cor.plot <- ggplot(data = subplot, aes(y = dir, x = value)) +
    geom_smooth(span = 0.9, method = 'glm') +
    facet_grid(stage ~ sa) +  # Note: + UL if using routes.
    labs(x = i, y = "Azimuth") +
    scale_y_continuous(breaks = c(0, 90, 180, 270, 360)) +
    theme_bw() +
    theme(panel.grid.major = element_blank())
  
  ggsave(paste0("D:/Papers/WEMigration2022/Outputs/cor_", i, ".png"), height = 5, width = 9, dpi = 300, plot = cor.plot)
}

# 11. Projection: Prepare Background Data and Environmental Layers
# Note: Creates grid for projections, extracts mean environmental values from NetCDF.
spcoords <- matrix(c(rep(seq(-20, 180, by = 2.5), 39), rep(seq(-20, 75, by = 2.5), each = 81)), ncol = 2)
proj.dt <- data.frame(spcoords)
colnames(proj.dt) <- c("lon", "lat")

spcoords[, 1][spcoords[, 1] < -1.25] <- spcoords[, 1][spcoords[, 1] < -1.25] + 360  # Adjust longitudes.

# Vegetation (inverted).
notp <- rast("D:/Papers/WEMigration2022/Datainput/NoTP/pi_notp.nc")
notp.veg <- rast("D:/Papers/WEMigration2022/Datainput/NoTP/pi_notp.nc", subds = "VEG")
ext(notp.veg) <- ext(notp)
proj.dt$evi <- 100 - extract(notp.veg, spcoords)$`VEG_lsmpft=1`

# Mean temperature.
notp.meanTem <- mean(notp[[1:12]])
proj.dt$tprt <- extract(notp.meanTem[[1]], spcoords)$mean

# Mean precipitation.
notp.meanPrec <- mean(notp[[13:24]])
proj.dt$prec <- extract(notp.meanPrec[[1]], spcoords)$mean

# Elevation.
proj.dt$ele <- extract(notp[[25]], spcoords)$TOPO

# 12. Species Distribution Modeling (SDM) with MaxEnt
# Note: Projects breeding/wintering sites for scenarios with/without Tibetan Plateau.
# Uses eBird data filtered by migration weeks.
Sys.setenv(JAVA_HOME = 'C:\\Program Files\\Java\\jre1.8.0_351')  # Set Java path for MaxEnt.

birdlist <- read.csv("D:/Papers/WEMigration2022/Datainput/birdlist.csv")

# Load eBird data.
ebirdfile <- list.files("D:/Papers/WEMigration2022/Datainput/eBird/90perc_points", full.names = TRUE)
ebirdname <- gsub("_", " ", basename(ebirdfile)) %>% sub(" points 90perc.csv", "")
ebirdfile <- ebirdfile[ebirdname %in% birdlist$English.name]

ebird <- lapply(ebirdfile, read.csv)
names(ebird) <- ebirdname[ebirdname %in% birdlist$English.name]
ebird <- bind_rows(ebird, .id = 'sp')

allbird <- unique(ebird$sp)

# Assign winter/breeding based on weeks (simplified; no date differences calculated).
ebird$win_bre <- NA
ebird$win_bre[ebird$week >= 49 | ebird$week < 8] <- "win"
ebird$win_bre[ebird$week >= 20 & ebird$week < 32] <- "bre"
ebird <- na.omit(ebird)

# Prepare environmental layers for win/bre (with TP).
tp <- rast("D:/Papers/WEMigration2022/Datainput/NoTP/pi.nc")  # Assuming this is TP data; filename suggests NoTP, check.
tprt.win <- mean(tp[[c(1:2, 12)]])
tprt.bre <- mean(tp[[6:8]])
pre.win <- mean(tp[[c(13:14, 24)]])
pre.bre <- mean(tp[[18:20]])
elevation <- tp[[25]]

wind_u <- tp[[c(26:37)]]
wind_v <- tp[[c(38:49)]]
wind_speed_tp_bre <- mean(sqrt(wind_u^2 + wind_v^2)[[6:8]])
wind_speed_tp_win <- mean(sqrt(wind_u^2 + wind_v^2)[[c(1, 11:12)]])

veg <- raster(tp.veg[[1]])

env.win <- stack(rast(list(tprt.win, pre.win, elevation, wind_speed_tp_win)), veg)
env.bre <- stack(rast(list(tprt.bre, pre.bre, elevation, wind_speed_tp_bre)), veg)
names(env.win) <- names(env.bre) <- c("Temperature", "Precipitation", "Topography", "Wind", "Vegetation")

# Prepare projected environments (without TP).
ptprt.win <- mean(notp[[c(1:2, 12)]])
ptprt.bre <- mean(notp[[6:8]])
ppre.win <- mean(notp[[c(13:14, 24)]])
ppre.bre <- mean(notp[[18:20]])
pelevation <- notp[[25]]

wind_u <- notp[[c(26:37)]]
wind_v <- notp[[c(38:49)]]
wind_speed_notp_bre <- mean(sqrt(wind_u^2 + wind_v^2)[[6:8]])
wind_speed_notp_win <- mean(sqrt(wind_u^2 + wind_v^2)[[c(1, 11:12)]])

pveg <- raster(notp.veg[[1]])

penv.win <- stack(rast(list(ptprt.win, ppre.win, pelevation, wind_speed_notp_win)), pveg)
penv.bre <- stack(rast(list(ptprt.bre, ppre.bre, pelevation, wind_speed_notp_bre)), pveg)
names(penv.win) <- names(penv.bre) <- c("Temperature", "Precipitation", "Topography", "Wind", "Vegetation")

# Load world boundaries for plotting.
countries <- maps::map("world", plot = FALSE)
countries <- maptools::map2SpatialLines(countries, proj4string = CRS("+proj=longlat"))
asia_africa <- rgdal::readOGR("D:/Materials/shp/Boundary/Asia_Africa.shp")

# Save workspace (optional).
# save.image("D:/Papers/WEMigration2022/Outputs/1221.RData")

# Run MaxEnt for each season (win/bre) and species.
wb <- lapply(c("win", "bre"), function(ss) {
  tryCatch({
    singwb <- ebird[ebird$win_bre == ss, ]
    
    env <- if (ss == "win") env.win else env.bre
    penv <- if (ss == "win") penv.win else penv.bre
    
    eachsp <- lapply(allbird, function(sp) {
      tryCatch({
        cat(sp, "\n")
        
        singsp <- singwb[singwb$sp == sp, ]
        singspcoords <- SpatialPoints(cbind(singsp$lon, singsp$lat), proj4string = CRS("+proj=longlat +datum=WGS84"))
        singspcoords <- singspcoords[asia_africa]
        
        # Train MaxEnt model.
        me <- maxent(env, singspcoords, path = paste0("D:/Papers/WEMigration2022/Outputs/MaxEnt/Initial/", ss, "_", sp), 
                     args = c('responsecurves=true', 'jackknife=true', 'allowpartialdata=true', 'betamultiplier=1',
                              'beta_threshold=5', 'beta_hinge=4', 'randomseed=true', 'writeplotdata=true', 'threads=4',
                              'maximumiterations=1000', 'replicates=1', 'replicatetype=bootstrap', 'randomtestpoints=25'))
        
        # Predict to no-TP scenario.
        pred <- predict(me, penv, ext = extent(penv), filename = paste0("D:/Papers/WEMigration2022/Outputs/MaxEnt/Notp/", ss, "_", sp, ".tif"),
                        args = c('responsecurves=true', 'jackknife=true', 'allowpartialdata=true', 'betamultiplier=1',
                                 'beta_threshold=5', 'beta_hinge=4', 'randomseed=true', 'writeplotdata=true', 'threads=4',
                                 'maximumiterations=1000', 'replicates=1', 'randomtestpoints=25'))
        
        # Plot prediction.
        png(file = paste0("D:/Papers/WEMigration2022/Outputs/MaxEnt/Notp/", ss, "_", sp, ".png"), width = 3600, height = 1800, res = 600)
        splot <- levelplot(pred$layer, margin = FALSE, col.regions = viridis, at = seq(0, 1, len = 100)) +
          latticeExtra::layer(sp.lines(countries))
        print(splot)
        dev.off()
        
        # Predict to current (TP) scenario.
        pred2 <- predict(me, env, ext = extent(env), filename = paste0("D:/Papers/WEMigration2022/Outputs/MaxEnt/present/", ss, "_", sp, ".tif"),
                         args = c('responsecurves=true', 'jackknife=true', 'allowpartialdata=true', 'betamultiplier=1',
                                  'beta_threshold=5', 'beta_hinge=4', 'randomseed=true', 'writeplotdata=true', 'threads=4',
                                  'maximumiterations=1000', 'replicates=1', 'randomtestpoints=25'))
        
        # Plot prediction.
        png(file = paste0("D:/Papers/WEMigration2022/Outputs/MaxEnt/present/", ss, "_", sp, ".png"), width = 3600, height = 1800, res = 600)
        splot <- levelplot(pred2$layer, margin = FALSE, col.regions = viridis, at = seq(0, 1, len = 100)) +
          latticeExtra::layer(sp.lines(countries))
        print(splot)
        dev.off()
        
        return(list(pred, pred2))
      }, error = function(e) { cat("ERROR :", conditionMessage(e), "\n") })
    })
    
    return(eachsp)
  }, error = function(e) { cat("ERROR :", conditionMessage(e), "\n") })
})

# 13. Post-Processing SDM Outputs: Crop, Normalize, and Analyze Distributions
# Note: Crops rasters to Asia/Africa, normalizes to 0-1, extracts high-probability areas.
asia <- terra::vect("D:/Materials/shp/Boundary/Asia.shp")
africa <- terra::vect("D:/Materials/shp/Boundary/Africa.shp")
asia_africa <- vect(asia_africa)  # Assuming union or combined vector.

# Function to crop raster, removing NAs and distant cells.
CroppedRaster <- function(x, na.value = NA) {
  if (!is.na(na.value)) x[x == na.value] <- NA
  if (canProcessInMemory(x, n = 2)) {
    x.matrix <- is.na(as.matrix(x))
    colNotNA <- which(colSums(x.matrix) != nrow(x))
    rowNotNA <- which(rowSums(x.matrix) != ncol(x))
    croppedExtent <- extent(x, rowNotNA[1], rowNotNA[length(rowNotNA)], colNotNA[1], colNotNA[length(colNotNA)])
    return(crop(x, croppedExtent))
  } else {
    xNA <- is.na(x)
    colNotNA <- which(colSums(xNA) != nrow(x))
    rowNotNA <- which(rowSums(xNA) != ncol(x))
    croppedExtent <- extent(x, rowNotNA[1], rowNotNA[length(rowNotNA)], colNotNA[1], colNotNA[length(colNotNA)])
    return(crop(x, croppedExtent))
  }
}

# Function to normalize raster to 0-1.
raster01 <- function(r) {
  minmax_r <- range(values(r), na.rm = TRUE)
  (r - minmax_r[1]) / diff(minmax_r)
}

# Process winter rasters (no-TP and TP).
winraster_notp <- lapply(1:length(wb[[1]]), function(rr) {
  r <- wb[[1]][[rr]][[1]]
  if (is.null(r)) return(NULL)
  r <- rast(r)
  r <- rotate(r)
  r <- mask(r, asia_africa)
  r <- raster01(r)
  r <- CroppedRaster(raster(r))
  r <- as.data.frame(r, xy = TRUE)
  r <- na.omit(r)
  return(r)
})

winraster_tp <- lapply(wb[[1]], function(rr) {
  r <- rr[[2]]
  if (is.null(r)) return(NULL)
  r <- rast(r)
  r <- rotate(r)
  r <- mask(r, asia_africa)
  r <- raster01(r)
  r <- CroppedRaster(raster(r))
  r <- as.data.frame(r, xy = TRUE)
  r <- na.omit(r)
  return(r)
})

# Process breeding rasters (similarly).
breraster_notp <- lapply(wb[[2]], function(rr) {
  r <- rr[[1]]
  if (is.null(r)) return(NULL)
  r <- rast(r)
  r <- rotate(r)
  r <- mask(r, asia_africa)
  r <- raster01(r)
  r <- CroppedRaster(raster(r))
  r <- as.data.frame(r, xy = TRUE)
  r <- na.omit(r)
  return(r)
})

breraster_tp <- lapply(wb[[2]], function(rr) {
  r <- rr[[2]]
  if (is.null(r)) return(NULL)
  r <- rast(r)
  r <- rotate(r)
  r <- mask(r, asia_africa)
  r <- raster01(r)
  r <- CroppedRaster(raster(r))
  r <- as.data.frame(r, xy = TRUE)
  r <- na.omit(r)
  return(r)
})

names(winraster_notp) <- names(breraster_notp) <- names(winraster_tp) <- names(breraster_tp) <- allbird

# Convert back to rasters.
winraster_notp1 <- lapply(winraster_notp, function(df) {
  if (is.null(df)) return(NULL)
  df$x <- plyr::round_any(df$x, 2.5)
  df$y <- plyr::round_any(df$y, 2.5)
  rast(df, type = "xyz")
})
winraster_tp1 <- lapply(winraster_tp, reraster)  # Assuming reraster function defined as above.

breraster_notp1 <- lapply(breraster_notp, reraster)
breraster_tp1 <- lapply(breraster_tp, reraster)

# Compute mean distributions.
meanwin_notp <- raster01(mean(rast(winraster_notp1)))
meanbre_notp <- raster01(mean(rast(breraster_notp1)))

meanwin_tp <- raster01(mean(rast(winraster_tp1)))
meanbre_tp <- raster01(mean(rast(breraster_tp1)))

# Save mean rasters.
# writeRaster(... ) lines as in original.

# Plot mean distributions using ggplot and terra.
continent <- vect("/Users/wenyuanzhang/Library/CloudStorage/OneDrive-McGillUniversity/Research/Boundary/world/continent.shp")  # Adjust path.
p <- ext(-20, 180, -40, 90)
subcontinent <- crop(continent, p)

# Example for no-TP breeding (repeat for others).
notp.bre.plot <- ggplot(subcontinent) +
  geom_spatraster(data = meanbre_notp) +
  geom_spatvector(fill = "transparent", colour = "grey50") +
  scale_x_continuous(expand = c(0, 0), limits = c(-20, 180)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_whitebox_c(palette = "soft", n.breaks = 5) +
  labs(fill = "Occurence Probability") +
  guides(fill = guide_colourbar(barheight = 15)) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.title = element_text(angle = -90, hjust = 0.5),
        legend.title.position = "right", legend.box.spacing = unit(0, "cm"), legend.justification = "centre")

ggsave("/Users/wenyuanzhang/Library/CloudStorage/OneDrive-McGillUniversity/Research/WEMigration2022/Outputs/notpbreplot.png", 
       plot = notp.bre.plot, width = 20, height = 15, dpi = 300, units = "cm")

# Similar plots for other combinations.

# 14. Calculate Centroids and Directions Between Sites
# Note: Extracts centroids from high-probability contours, calculates bearings and distances.
conpolycentre <- function(r, l = c(0.9, 0.99)) {
  if (is.null(r)) return(NULL)
  con.rst <- as.contour(raster01(r), levels = l)
  poly.con <- as.polygons(con.rst)
  area.poly <- expanse(poly.con)
  poly.con <- poly.con[which.max(area.poly)]
  cent.poly <- centroids(poly.con)
  coord.cent <- geom(cent.poly)
  matrix(c(coord.cent[, 3], coord.cent[, 4]), ncol = 2)
}

winraster_notp1_cent <- lapply(winraster_notp1, conpolycentre)
winraster_tp1_cent <- lapply(winraster_tp1, conpolycentre)
breraster_notp1_cent <- lapply(breraster_notp1, conpolycentre)
breraster_tp1_cent <- lapply(breraster_tp1, conpolycentre)

dir.winbre <- lapply(1:length(winraster_notp1_cent), function(i) {
  if (is.null(breraster_notp1_cent[[i]]) | is.null(winraster_notp1_cent[[i]]) | 
      is.null(breraster_tp1_cent[[i]]) | is.null(winraster_tp1_cent[[i]])) return(NULL)
  dirr.notp <- geosphere::bearingRhumb(breraster_notp1_cent[[i]], winraster_notp1_cent[[i]])
  dirr.tp <- geosphere::bearingRhumb(breraster_tp1_cent[[i]], winraster_tp1_cent[[i]])
  data.frame(bird = names(breraster_notp1_cent)[i], notp = dirr.notp, tp = dirr.tp,
             notp.bre.lon = breraster_notp1_cent[[i]][, 1], notp.bre.lat = breraster_notp1_cent[[i]][, 2],
             notp.win.lon = winraster_notp1_cent[[i]][, 1], notp.win.lat = winraster_notp1_cent[[i]][, 2],
             tp.bre.lon = breraster_tp1_cent[[i]][, 1], tp.bre.lat = breraster_tp1_cent[[i]][, 2],
             tp.win.lon = winraster_tp1_cent[[i]][, 1], tp.win.lat = winraster_tp1_cent[[i]][, 2])
})

dir.winbre <- bind_rows(dir.winbre)

# Adjust coordinates for plotting (ensure breeding north of wintering).
dir.winbre$notp.bre.lon.adjust <- ifelse(dir.winbre$notp.bre.lon >= dir.winbre$notp.win.lon, dir.winbre$notp.bre.lon, dir.winbre$notp.win.lon)
dir.winbre$notp.bre.lat.adjust <- ifelse(dir.winbre$notp.bre.lat >= dir.winbre$notp.win.lat, dir.winbre$notp.bre.lat, dir.winbre$notp.win.lat)
dir.winbre$notp.win.lon.adjust <- ifelse(dir.winbre$notp.bre.lon >= dir.winbre$notp.win.lon, dir.winbre$notp.win.lon, dir.winbre$notp.bre.lon)
dir.winbre$notp.win.lat.adjust <- ifelse(dir.winbre$notp.bre.lat >= dir.winbre$notp.win.lat, dir.winbre$notp.win.lat, dir.winbre$notp.bre.lat)

dir.winbre$tp.bre.lon.adjust <- ifelse(dir.winbre$tp.bre.lon >= dir.winbre$tp.win.lon, dir.winbre$tp.bre.lon, dir.winbre$tp.win.lon)
dir.winbre$tp.bre.lat.adjust <- ifelse(dir.winbre$tp.bre.lat >= dir.winbre$tp.win.lat, dir.winbre$tp.bre.lat, dir.winbre$tp.win.lat)
dir.winbre$tp.win.lon.adjust <- ifelse(dir.winbre$tp.bre.lon >= dir.winbre$tp.win.lon, dir.winbre$tp.win.lon, dir.winbre$tp.bre.lon)
dir.winbre$tp.win.lat.adjust <- ifelse(dir.winbre$tp.bre.lat >= dir.winbre$tp.win.lat, dir.winbre$tp.win.lat, dir.winbre$tp.bre.lat)

# Calculate distances.
notp.bre <- vect(dir.winbre, geom = c("notp.bre.lon.adjust", "notp.bre.lat.adjust"), crs = "epsg:4326")
notp.win <- vect(dir.winbre, geom = c("notp.win.lon.adjust", "notp.win.lat.adjust"), crs = "epsg:4326")
tp.bre <- vect(dir.winbre, geom = c("tp.bre.lon.adjust", "tp.bre.lat.adjust"), crs = "epsg:4326")
tp.win <- vect(dir.winbre, geom = c("tp.win.lon.adjust", "tp.win.lat.adjust"), crs = "epsg:4326")

dir.winbre$notp.dis <- diag(terra::distance(notp.bre, notp.win, unit = "m")) / 1000
dir.winbre$tp.dis <- diag(terra::distance(tp.bre, tp.win, unit = "m")) / 1000

# Recalculate bearings with adjusted coords.
dir.winbre$notp <- geosphere::bearingRhumb(matrix(c(dir.winbre$notp.bre.lon.adjust, dir.winbre$notp.bre.lat.adjust), ncol = 2),
                                           matrix(c(dir.winbre$notp.win.lon.adjust, dir.winbre$notp.win.lat.adjust), ncol = 2))

dir.winbre$tp <- geosphere::bearingRhumb(matrix(c(dir.winbre$tp.bre.lon.adjust, dir.winbre$tp.bre.lat.adjust), ncol = 2),
                                         matrix(c(dir.winbre$tp.win.lon.adjust, dir.winbre$tp.win.lat.adjust), ncol = 2))

# Categorize directions into 16 bins.
dir.winbre$category.notp <- cut(dir.winbre$notp, breaks = seq(0, 360, by = 22.5), labels = 1:16, include.lowest = TRUE)
dir.winbre$category.tp <- cut(dir.winbre$tp, breaks = seq(0, 360, by = 22.5), labels = 1:16, include.lowest = TRUE)

# Summarize for rose plots.
cir.bar.notp <- dir.winbre %>% group_by(category.notp) %>% summarise(sum.length = n())
cir.bar.tp <- dir.winbre %>% group_by(category.tp) %>% summarise(sum.length = n())

# Fill missing categories with 0.
cir.bar.notp <- complete(cir.bar.notp, category.notp = factor(1:16), fill = list(sum.length = 0))
cir.bar.tp <- complete(cir.bar.tp, category.tp = factor(1:16), fill = list(sum.length = 0))

# Labels for directions.
dir.labels <- c("N", "", "NE", "", "E", "", "SE", "", "S", "", "SW", "", "W", "", "NW", "")

# Plot rose diagrams.
notp.bar <- ggplot(cir.bar.notp) +
  geom_hline(aes(yintercept = y), data.frame(y = seq(0, 20, 5)), color = "lightgrey") +
  geom_col(aes(x = category.notp, y = sum.length, fill = sum.length), position = "dodge2", alpha = 0.9) +
  geom_vline(aes(xintercept = x), data.frame(x = seq(1, 16, 2)), color = "grey65", linetype = "dashed") +
  scale_fill_gradient2(midpoint = 9, low = "#FEFDD6", mid = "#EEF5C5", high = "#AAD07F", breaks = seq(0, 20, 5)) +
  scale_x_discrete(breaks = 1:16, labels = dir.labels) +
  labs(fill = "Number of species") +
  theme(panel.background = element_rect(fill = "white", colour = "white"), panel.grid = element_blank(),
        axis.title = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        axis.text.x = element_text(colour = "black", face = "bold"), legend.title.position = "right",
        legend.title = element_text(angle = -90)) +
  coord_polar(start = -0.2)

tp.bar <- ggplot(cir.bar.tp) +
  geom_hline(aes(yintercept = y), data.frame(y = seq(0, 16, 4)), color = "lightgrey") +
  geom_col(aes(x = category.tp, y = sum.length, fill = sum.length), position = "dodge2", alpha = 0.9) +
  geom_vline(aes(xintercept = x), data.frame(x = seq(1, 16, 2)), color = "grey65", linetype = "dashed") +
  scale_fill_gradient2(midpoint = 9, low = "#FEFDD6", mid = "#EEF5C5", high = "#AAD07F", breaks = seq(0, 20, 5)) +
  scale_x_discrete(breaks = 1:16, labels = dir.labels) +
  labs(fill = "Number of species") +
  theme(panel.background = element_rect(fill = "white", colour = "white"), panel.grid = element_blank(),
        axis.title = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        axis.text.x = element_text(colour = "black", face = "bold"), legend.title.position = "right",
        legend.title = element_text(angle = -90)) +
  coord_polar(start = -0.2)

# Plot direction arrows on maps.
notp.plot <- ggplot(subcontinent) +
  geom_spatvector(fill = "white") +
  geom_segment(data = dir.winbre, aes(x = notp.bre.lon.adjust, y = notp.bre.lat.adjust, 
                                      xend = notp.win.lon.adjust, yend = notp.win.lat.adjust, colour = notp.dis),
               arrow = arrow(length = unit(0.3, "cm")), lwd = 1, alpha = 0.5) +
  scale_colour_gradient2(midpoint = 6500, low = "red", mid = "orange", high = "lightyellow", breaks = c(3000, 6000, 9000, 12000)) +
  labs(colour = "Distance between breeding and wintering centres (km)") +
  guides(colour = guide_colourbar(barwidth = 18)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_line(colour = "transparent"), axis.title = element_blank(),
        legend.title.position = "top", legend.position = "bottom", legend.box.spacing = unit(0, "cm"),
        legend.justification = "center")

tp.plot <- ggplot(subcontinent) +
  geom_spatvector(fill = "white") +
  geom_segment(data = dir.winbre, aes(x = tp.bre.lon.adjust, y = tp.bre.lat.adjust, 
                                      xend = tp.win.lon.adjust, yend = tp.win.lat.adjust, colour = tp.dis),
               arrow = arrow(length = unit(0.3, "cm")), lwd = 1, alpha = 0.5) +
  scale_colour_gradient2(midpoint = 6500, low = "red", mid = "orange", high = "lightyellow", breaks = c(3000, 6000, 9000, 12000)) +
  labs(colour = "Distance between breeding and wintering centres (km)") +
  guides(colour = guide_colourbar(barwidth = 18)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_line(colour = "transparent"), axis.title = element_blank(),
        legend.title.position = "top", legend.position = "bottom", legend.box.spacing = unit(0, "cm"),
        legend.justification = "center")

# Overlay rose plots on maps.
notp.overlay <- notp.plot + annotation_custom(grob = ggplotGrob(notp.bar), xmin = 100, xmax = 180, ymin = -38, ymax = 20)
ggsave("/Users/wenyuanzhang/Library/CloudStorage/OneDrive-McGillUniversity/Research/WEMigration2022/Outputs/notpoverlap.png", 
       plot = notp.overlay, width = 20, height = 15, units = "cm", dpi = 300)

tp.overlay <- tp.plot + annotation_custom(grob = ggplotGrob(tp.bar), xmin = 100, xmax = 180, ymin = -38, ymax = 20)
ggsave("/Users/wenyuanzhang/Library/CloudStorage/OneDrive-McGillUniversity/Research/WEMigration2022/Outputs/tpoverlap.png", 
       plot = tp.overlay, width = 20, height = 15, units = "cm", dpi = 300)

# 15. Compare Distributions (TP vs NoTP)
# Note: Computes differences, means, and density plots for lat/lon.
diftp <- function(r1, r2) {
  if (is.null(r1) | is.null(r2)) return(NULL)
  rdif <- sum(r1, -r2, na.rm = TRUE)
  is.na(rdif) <- is.na(r1) & is.na(r2)
  rdif
}

winraster <- lapply(allbird, function(sp) diftp(winraster_notp1[[sp]], winraster_tp1[[sp]]))
breraster <- lapply(allbird, function(sp) diftp(breraster_notp1[[sp]], breraster_tp1[[sp]]))
names(winraster) <- names(breraster) <- allbird

# Mean differences.
meanwin <- mean(rast(winraster))
meanbre <- mean(rast(breraster))

# Save.
# writeRaster(meanwin, "D:/Papers/WEMigration2022/Outputs/distribution/meanwin1.tif", overwrite = TRUE)
# etc.

# Density plots for high-prob areas.
windf_tp <- bind_rows(lapply(winraster_tp, function(s) s[s$layer >= 0.99, ]))
windf_tp$id <- "TP"
windf_tp$season <- "Winter"

windf_notp <- bind_rows(lapply(winraster_notp, function(s) s[s$layer >= 0.99, ]))
windf_notp$id <- "NoTP"
windf_notp$season <- "Winter"

bredf_notp <- bind_rows(lapply(breraster_notp, function(s) s[s$layer >= 0.99, ]))
bredf_notp$id <- "NoTP"
bredf_notp$season <- "Breed"

bredf_tp <- bind_rows(lapply(breraster_tp, function(s) s[s$layer >= 0.99, ]))
bredf_tp$id <- "TP"
bredf_tp$season <- "Breed"

pltdt <- rbind(windf_notp, windf_tp, bredf_notp, bredf_tp)

deplot_lat <- ggplot(pltdt, aes(x = y, group = id, fill = id, colour = id)) +
  geom_density(alpha = 0.6, adjust = 2) +
  ylab("") + xlab("Latitude") +
  theme_bw() + theme(legend.title = element_blank()) +
  facet_wrap(~season)

ggsave("D:/Papers/WEMigration2022/Outputs/distribution/latdensity1.png", deplot_lat, dpi = 200, height = 4, width = 6)

deplot_lon <- ggplot(pltdt, aes(x = x, group = id, fill = id, colour = id)) +
  geom_density(alpha = 0.6, adjust = 2) +
  ylab("") + xlab("Longitude") +
  theme_bw() + theme(legend.title = element_blank()) +
  facet_wrap(~season)

ggsave("D:/Papers/WEMigration2022/Outputs/distribution/londensity1.png", deplot_lon, dpi = 200, height = 4, width = 6)

# 16. Project Migration Routes
# Note: Similar replacements for wind processing in projections.
wind_u <- notp[[c(28:30, 34:36)]]
wind_v <- notp[[c(40:42, 46:48)]]

wind_speed <- sqrt(wind_u^2 + wind_v^2)
wind_dir <- 180 + 180 * atan2(wind_v, wind_u) / pi  # Adjust to match convention.

####
names(wind_speed) <- rep('speed', 6)
names(wind_dir) <- rep('direction', 6)

wind_layers <- lapply(1:nlyr(wind_speed), function(l) wind2raster_custom(data.frame(lon = xFromCell(wind_speed[[l]], 1:ncell(wind_speed[[l]])), lat = yFromCell(wind_speed[[l]], 1:ncell(wind_speed[[l]])), dir = values(wind_dir[[l]]), speed = values(wind_speed[[l]]))))  # Reconstruct df if needed.

names(wind_layers) <- c(seq(lubridate::ymd_hms(paste(2019,3,15,12,00,00, sep="-")),  
                            lubridate::ymd_hms(paste(2019,5,16,12,00,00, sep="-")),by="1 months"),
                        seq(lubridate::ymd_hms(paste(2019,9,5,12,00,00, sep="-")),  
                            lubridate::ymd_hms(paste(2019,11,6,12,00,00, sep="-")),by="1 months"))

# The projection loop uses custom_flow.dispersion similarly.
# ... (unchanged, but ensure custom functions are used)

# 17. Simulate Migration Paths
# ... (unchanged)
