library(data.table)
library(dplyr)
library(sf)
library(terra)
library(curl)
library(lubridate)
library(purrr)
library(arrow)

# ============================================================
# Load station/data
# ============================================================

data_Temp <- read_parquet(
  "https://raw.githubusercontent.com/TimGyger/SpaceTimeGPApprox/refs/heads/main/Data/Real_World_TEMP_PRCP.parquet"
)

data <- data_Temp

dat <- data %>%
  select(id, X, Y) %>%
  distinct()

# X/Y assumed to be NAD83 / Conus Albers, EPSG:5070
stations_sf <- st_as_sf(
  dat,
  coords = c("X", "Y"),
  crs = 5070
)

stations_ll <- st_transform(stations_sf, 4326)
coords <- st_coordinates(stations_ll)

stations_lonlat <- stations_ll |>
  mutate(
    lon = coords[, 1],
    lat = coords[, 2]
  ) |>
  st_drop_geometry() |>
  select(station_id = id, lat, lon)

stations <- stations_lonlat

# ============================================================
# USER INPUT
# ============================================================

start_date <- as.Date("2023-12-31")
end_date   <- as.Date("2025-08-31")

run_hour <- "00"

# Next-24h forecast from the last evening before/at target date:
# For target date D, use HRRR initialized at D, 00Z,
# and forecast hours f01...f48.
# For US stations, 00Z on D is usually the evening of D-1 local time.
fhrs <- 1:48

out_file <- "hrrr_daily_tmax_precip_2025.csv"

tmpdir <- file.path(tempdir(), "hrrr_grib_extract")
dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# HRRR URL
# ============================================================

hrrr_url <- function(init_date, fhr, run_hour = "00") {
  init_date <- as.Date(init_date, origin = "1970-01-01")
  ymd <- format(init_date, "%Y%m%d")
  
  sprintf(
    "https://noaa-hrrr-bdp-pds.s3.amazonaws.com/hrrr.%s/conus/hrrr.t%sz.wrfsfcf%02d.grib2",
    ymd, run_hour, fhr
  )
}

# ============================================================
# Download one GRIB message using byte ranges
# ============================================================

download_grib_message <- function(grib_url, pattern, out_grib, which_match = 1) {
  idx_url <- paste0(grib_url, ".idx")
  
  idx_txt <- tryCatch(
    readLines(idx_url, warn = FALSE),
    error = function(e) {
      stop("Could not read index file: ", idx_url)
    }
  )
  
  idx <- data.table::tstrsplit(idx_txt, ":", fixed = TRUE)
  
  idx_dt <- data.table(
    line = seq_along(idx_txt),
    msg = as.integer(idx[[1]]),
    byte_start = as.numeric(idx[[2]]),
    desc = idx_txt
  )
  
  hit <- idx_dt[grepl(pattern, desc)]
  
  if (nrow(hit) < 1) {
    stop(
      "No GRIB message found for pattern: ", pattern,
      "\nURL: ", grib_url
    )
  }
  
  if (nrow(hit) > 1) {
    message("Multiple matches for pattern: ", pattern)
    print(hit[, .(line, desc)])
    hit <- hit[which_match]
    message("Using match ", which_match, ": ", hit$desc)
  }
  
  next_start <- idx_dt[line == hit$line + 1, byte_start]
  
  range <- if (length(next_start) == 0) {
    paste0(hit$byte_start, "-")
  } else {
    paste0(hit$byte_start, "-", next_start - 1)
  }
  
  h <- curl::new_handle()
  curl::handle_setheaders(h, Range = paste0("bytes=", range))
  
  curl::curl_download(
    url = grib_url,
    destfile = out_grib,
    handle = h,
    quiet = TRUE
  )
  
  out_grib
}

# ============================================================
# Extract nearest HRRR grid cell values
# ============================================================

extract_nearest_hrrr <- function(grib_file, stations) {
  r <- terra::rast(grib_file)
  
  pts <- data.table(
    station_id = stations$station_id,
    lon = stations$lon,
    lat = stations$lat
  )
  
  v <- terra::vect(
    pts[, .(lon, lat)],
    geom = c("lon", "lat"),
    crs = "EPSG:4326"
  )
  
  vals <- terra::extract(r, v, method = "simple")[, 2]
  
  data.table(
    station_id = stations$station_id,
    value = vals
  )
}

# ============================================================
# Main extraction loop
# ============================================================

dates <- seq.Date(start_date, end_date, by = "day")

all_daily <- list()

for (d_raw in dates) {
  
  target_date <- as.Date(d_raw, origin = "1970-01-01")
  
  # HRRR 00Z initialization.
  # In the US, this corresponds to the previous local evening.
  init_date <- target_date
  
  message(
    "Processing target date ",
    format(target_date, "%Y-%m-%d"),
    " using HRRR init ",
    format(init_date, "%Y-%m-%d"),
    " ",
    run_hour,
    "Z"
  )
  
  hourly_temp <- list()
  hourly_prec <- list()
  
  for (fh in fhrs) {
    
    url <- hrrr_url(init_date, fh, run_hour)
    
    temp_file <- file.path(
      tmpdir,
      sprintf(
        "tmp_init_%s_target_%s_f%02d.grib2",
        format(init_date, "%Y%m%d"),
        format(target_date, "%Y%m%d"),
        fh
      )
    )
    
    prec_file <- file.path(
      tmpdir,
      sprintf(
        "apcp_init_%s_target_%s_f%02d.grib2",
        format(init_date, "%Y%m%d"),
        format(target_date, "%Y%m%d"),
        fh
      )
    )
    
    # --------------------------------------------------------
    # Temperature
    #
    # IMPORTANT:
    # In your current terra/HRRR setup, extracted TMP values
    # are already in Celsius, so do NOT subtract 273.15.
    # --------------------------------------------------------
    
    tryCatch({
      
      download_grib_message(
        grib_url = url,
        pattern = ":TMP:2 m above ground:.*fcst",
        out_grib = temp_file,
        which_match = 1
      )
      
      temp_vals <- extract_nearest_hrrr(temp_file, stations)
      
      temp_vals[, `:=`(
        date = target_date,
        init_date = init_date,
        run_hour = run_hour,
        fhr = fh,
        temp_C = value
      )]
      
      hourly_temp[[length(hourly_temp) + 1]] <-
        temp_vals[, .(station_id, date, init_date, run_hour, fhr, temp_C)]
      
    }, error = function(e) {
      message(
        "Temperature failed for target ",
        target_date,
        " init ",
        init_date,
        " f",
        fh,
        ": ",
        e$message
      )
    })
    
    # --------------------------------------------------------
    # Precipitation
    #
    # Use hourly accumulation:
    # f24 -> 23-24 hour acc fcst
    # f25 -> 24-25 hour acc fcst
    # ...
    # f47 -> 46-47 hour acc fcst
    # --------------------------------------------------------
    
    tryCatch({
      
      prec_pattern <- sprintf(
        ":APCP:surface:%d-%d hour acc fcst:",
        fh - 1,
        fh
      )
      
      download_grib_message(
        grib_url = url,
        pattern = prec_pattern,
        out_grib = prec_file,
        which_match = 1
      )
      
      prec_vals <- extract_nearest_hrrr(prec_file, stations)
      
      prec_vals[, `:=`(
        date = target_date,
        init_date = init_date,
        run_hour = run_hour,
        fhr = fh,
        precip_mm = value
      )]
      
      hourly_prec[[length(hourly_prec) + 1]] <-
        prec_vals[, .(station_id, date, init_date, run_hour, fhr, precip_mm)]
      
    }, error = function(e) {
      message(
        "Precip failed for target ",
        target_date,
        " init ",
        init_date,
        " f",
        fh,
        ": ",
        e$message
      )
    })
    
    if (file.exists(temp_file)) file.remove(temp_file)
    if (file.exists(prec_file)) file.remove(prec_file)
  }
  
  if (length(hourly_temp) == 0 || length(hourly_prec) == 0) {
    message("Skipping ", target_date, " because no data were extracted.")
    next
  }
  
  temp_dt <- rbindlist(hourly_temp, fill = TRUE)
  prec_dt <- rbindlist(hourly_prec, fill = TRUE)
  
  temp_dt[, lead_day := fifelse(fhr <= 24, "day1", "day2")]
  prec_dt[, lead_day := fifelse(fhr <= 24, "day1", "day2")]
  
  daily_temp <- dcast(
    temp_dt[
      ,
      .(
        hrrr_tmax_C = max(temp_C, na.rm = TRUE),
        n_temp_hours = sum(!is.na(temp_C))
      ),
      by = .(station_id, date, init_date, run_hour, lead_day)
    ],
    station_id + date + init_date + run_hour ~ lead_day,
    value.var = c("hrrr_tmax_C", "n_temp_hours")
  )
  
  daily_prec <- dcast(
    prec_dt[
      ,
      .(
        hrrr_precip_mm = sum(precip_mm, na.rm = TRUE),
        n_precip_hours = sum(!is.na(precip_mm))
      ),
      by = .(station_id, date, init_date, run_hour, lead_day)
    ],
    station_id + date + init_date + run_hour ~ lead_day,
    value.var = c("hrrr_precip_mm", "n_precip_hours")
  )
  
  daily <- merge(
    daily_temp,
    daily_prec,
    by = c("station_id", "date", "init_date", "run_hour"),
    all = TRUE
  )
  
  print(head(daily))
  
  all_daily[[length(all_daily) + 1]] <- daily
}

result <- rbindlist(all_daily, fill = TRUE)

result <- merge(
  result,
  as.data.table(stations),
  by = "station_id",
  all.x = TRUE
)

setorder(result, station_id, date)

result_day1 <- result %>%
  select(
    station_id,
    date,
    init_date_day1 = init_date,
    run_hour,
    hrrr_tmax_C_day1,
    hrrr_precip_mm_day1,
    n_temp_hours_day1,
    n_precip_hours_day1,
    lat,
    lon
  )

result_day2 <- result %>%
  mutate(date = date + 1) %>%
  select(
    station_id,
    date,
    init_date_day2 = init_date,
    hrrr_tmax_C_day2,
    hrrr_precip_mm_day2,
    n_temp_hours_day2,
    n_precip_hours_day2
  )

result_aligned <- result_day1 %>%
  left_join(result_day2, by = c("station_id", "date")) %>%
  arrange(station_id, date)

# Optional: keep only original date range
result_aligned <- result_aligned %>%
  filter(date >= start_date, date <= end_date)