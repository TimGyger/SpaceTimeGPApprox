#######################################################
### Packages
#######################################################

### Install/Load Packages
# Package names
packages <- c("rnoaa","dplyr", "arrow", "lubridate", "sf", "terra", "elevatr",
              "tidyr", "readr", "aws.s3", "purrr", "stringr")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# Random Fields
# Function to check and install a package if not already installed
install_if_missing <- function(package, version = NULL) {
  if (!requireNamespace(package, quietly = TRUE)) {
    if (!is.null(version)) {
      remotes::install_version(package, version)
    } else {
      install.packages(package)
    }
  }
}

###################
## Station Data
stations <- ghcnd_stations()
US_stations <- subset(stations,grepl("^US", stations$id) & !stations$state %in% c("AK","HI") & 
                        stations$last_year >= 2025 & stations$first_year <= 2024 & element %in% c("PRCP","TMAX"))


station_id <- unique(US_stations$id)[sample(1:30000,10)]

# download precipitation data (PRCP)
data <- meteo_pull_monitors(monitors = station_id, var = c("PRCP","TMAX"),keep_flags  = T,date_min = "2024-01-01",
                             date_max = "2025-08-31")

data <- data %>% filter(!is.na(tmax) & qflag_tmax == " " & qflag_prcp == " ") %>% select(id,date,tmax,prcp) %>% mutate(prcp = as.numeric(prcp)/10,
                                                                                                                       tmax = as.numeric(tmax)/10)
data <- data %>% right_join(US_stations, by = "id") %>% select(id,latitude, longitude, elevation, date) %>% mutate(year = year(date),doy  = yday(date))


###################
## Temperature Normals
# List all objects in the 1991–2020 daily normals folder
objs <- get_bucket(
  bucket = "noaa-normals-pds",
  prefix = "normals-daily/1991-2020/",
  region = "us-east-1",
  anonymous = TRUE,
  max = Inf
)

df <- data.frame(
  Key = sapply(objs, function(x) x$Key),
  LastModified = sapply(objs, function(x) x$LastModified),
  Size = sapply(objs, function(x) x$Size),
  stringsAsFactors = FALSE
)
df <- df %>%
  mutate(station_id = str_remove(basename(Key), "\\.csv$"))

head(df)
stations <- ghcnd_stations()
US_stations <- subset(stations,grepl("^US", stations$id) & !stations$state %in% c("AK","HI") & 
                        stations$last_year >= 2025 & stations$first_year <= 2024 & element == "TMAX")


station_id <- unique(US_stations$id)

df <- df %>%
  filter(station_id %in% station_id)


# Filter only CSV files
keys <- df$Key[sample(1:300,10)]

has_temperature <- function(key) {
  tryCatch({
    df <- s3read_using(
      FUN = read_csv,
      object = key,
      bucket = "noaa-normals-pds"
    )
    # Check if tmax exists in the "datatype" column
    df
  }, error = function(e) {
    FALSE
  })
}

results <- map(seq_along(keys), function(i) {
  key <- keys[i]
  res <- has_temperature(key)
  if ("DLY-TMAX-NORMAL" %in% names(res)){
    res <- res %>% select(STATION,DATE,LATITUDE,LONGITUDE,ELEVATION,month,day,`DLY-TMAX-NORMAL`) %>% drop_na()
    message(sprintf("Checked %d of %d keys", i, length(keys)))
    if (nrow(res) > 0) {
      print("good")
      res
    }
  } else {
    print("bad")
  }
  
})
results_clean <- keep(results, ~inherits(.x, "data.frame"))

# Combine them
df_all <- bind_rows(results_clean)

df_all$`DLY-TMAX-NORMAL` <- (df_all$`DLY-TMAX-NORMAL`-32)*5/9

###################
## Distance to coast
stations_sf <- sf::st_as_sf(df_all %>% distinct(STATION,.keep_all = T),
                            coords = c("LONGITUDE", "LATITUDE"),
                            crs = 4326)

# Load Natural Earth coastline (e.g. 1:50m resolution)
coastline <- rnaturalearth::ne_download(scale = 50,
                                        type = "coastline",
                                        category = "physical",
                                        returnclass = "sf")

# Transform to projected CRS (meters) for distances
stations_proj <- sf::st_transform(stations_sf, 3857)
coast_proj <- sf::st_transform(coastline, 3857)

# Compute distance to nearest coastline
stations_sf$distance_to_sea <- sf::st_distance(stations_proj, coast_proj) %>% 
  apply(1, min) # min distance across coastline segments

dat_final <- df_all %>% left_join(stations_sf %>% select(STATION,distance_to_sea),by = c("STATION"))

#################
## Slope and Aspect
stations <- dat_final %>%
  distinct(STATION, LATITUDE, LONGITUDE)

# Convert to sf points in WGS84
stations_sf <- sf::st_as_sf(stations, coords = c("LONGITUDE", "LATITUDE"), crs = 4326)

# --- 1. Build a bounding box around stations
bbox <- sf::st_bbox(stations_sf)
bbox_poly <- sf::st_as_sfc(bbox)
bbox_sf <- sf::st_sf(geometry = bbox_poly)

# --- 2. Download DEM for that area from AWS terrain tiles
# z = 10 or 12 for higher resolution (may take longer)
dem <- elevatr::get_elev_raster(locations = bbox_sf, z = 5, src = "aws")

# Convert to terra raster
dem_terra <- terra::rast(dem)

# --- 3. Compute slope and aspect
slope_r  <- terra::terrain(dem_terra, v = "slope", unit = "degrees", neighbors = 8)
aspect_r <- terra::terrain(dem_terra, v = "aspect", unit = "degrees", neighbors = 8)

# --- 4. Extract slope & aspect for each station
stations_vect <- terra::vect(stations_sf)
vals <- terra::extract(c(slope_r, aspect_r), stations_vect)

# Add slope/aspect back to stations data frame
stations_out <- stations %>%
  bind_cols(vals[, -1]) %>%   # remove ID column
  rename(slope_deg = slope, aspect_deg = aspect)

# --- 5. Join to your full dataset
dat_final_all <- dat_final %>%
  left_join(stations_out %>% select(STATION,slope_deg, aspect_deg), by = "STATION")

dat_final_all <- dat_final_all %>%
  mutate(
    month = as.numeric(month),
    day = as.numeric(day),
    # create day of year
    day_of_year = lubridate::yday(as.Date(paste0("2020-", month, "-", day))))

data_Temp_joined <- data %>%
  left_join(dat_final_all %>%
              select(STATION, day_of_year, `DLY-TMAX-NORMAL`, distance_to_sea, slope_deg, aspect_deg),
            by = c("id" = "STATION", "doy" = "day_of_year"))
data_Temp_joined <- data_Temp_joined %>% drop_na() %>% distinct(id,date,.keep_all = T)


##################
## Transform Coordinates
dat_final_all_sf <- sf::st_as_sf(
  data_Temp_joined,
  coords = c("longitude", "latitude"),
  crs = 4326
)

# Now transform to a projected CRS in meters (good for US)
dat_final_all_proj <- sf::st_transform(dat_final_all_sf, crs = 5070)

# Extract projected coordinates if you want numeric X/Y for GP
coords <- sf::st_coordinates(dat_final_all_proj)
dat_final_all_proj$X <- coords[,1]
dat_final_all_proj$Y <- coords[,2]

head(dat_final_all_proj)
dat_final_no_geom <- dat_final_all_proj %>%
  st_drop_geometry()

##################
## KGZones
dat_final_no_geom_sub <- dat_final_no_geom %>% distinct(id,.keep_all = T)
pts_sf <- st_as_sf(dat_final_no_geom_sub,coords = c("X","Y"), crs = 5070, remove = FALSE) 
pts_sf_ll <- sf::st_transform(pts_sf, crs = 4326)
KGZones <- kgc::LookupCZ(pts_sf_ll %>% select(id, geometry) %>% 
                           mutate(Longitude = st_coordinates(.)[,1],
                                                    Latitude  = st_coordinates(.)[,2],
                                                    rndCoord.lon = RoundCoordinates(Longitude),
                                                    rndCoord.lat  = RoundCoordinates(Latitude)) %>% select(id,Longitude,Latitude,rndCoord.lon,rndCoord.lat) )
pts_sf_ll <- data.frame(pts_sf_ll %>% select(id),KGZones)
dat_final_no_geom_zone <- dat_final_no_geom %>% left_join(pts_sf_ll, by = c("id"))
