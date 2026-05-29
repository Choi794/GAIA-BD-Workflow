# GAIA-BD Hexagon-Level Environmental Covariate Integration and Public Release Formatting
#
# This script extracts climatic, topographic, edaphic, and anthropogenic
# environmental covariates using the centroid coordinates of occupied
# 0.5° hexagonal grid cells. WWF terrestrial biome codes are assigned
# through spatial joins with the WWF biome layer. The resulting dataset
# is then standardized by applying concise variable names and a consistent
# column structure for public release.
#
# Input:
#   hexIDlevel_biodiversity_with_centroid_coordinates.csv
#
# Intermediate Output:
#   hexlevel_biodiversity_final.csv
#
# Final Output:
#   hexlevel_biodiversity.csv

rm(list = ls())

library(data.table)
library(raster)
library(sf)
library(dplyr)

# Load hexagon-level biodiversity summaries with centroid coordinates
setwd("F:/GAIA_datapaper/Incidence3/NumObs")

data <- fread("hexIDlevel_biodiversity_with_centroid_coordinates.csv")

# Ensure consistent column names
setnames(data, c("hexID", "LAT", "LON", "N", "K", "R"))

# Check data structure
str(data)
summary(data)

# Set working directory for environmental raster covariates
setwd("E:/FACAI/Covariate/Covraiate_data")

# List raster covariate files
cov_files <- list.files()

# Define covariate names
cov_names <- c(
  "aridity", "aspectcosine", "aspectsine",
  paste0("bio", 1:19),
  "dx", "dxx", "dy", "dyy",
  "elevation", "evapotspr", "humanfp",
  "pcurv", "roughness", "slope", "tcurv",
  "tpi", "tri"
)

# Initialize covariate columns
data[, (cov_names) := NA]

# Extract raster covariate values using hexagon centroid coordinates
for (i in seq_along(cov_files)) {
  
  raster_file <- raster(cov_files[i])
  
  # Reproject Human Footprint raster to WGS84 if needed
  if (i == 29) {
    raster_file <- projectRaster(
      raster_file,
      crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    )
  }
  
  extracted_values <- extract(
    x = raster_file,
    y = data[, .(LON, LAT)],
    method = "simple",
    df = TRUE
  )
  
  data[, (cov_names[i]) := extracted_values[, 2]]
  
  cat("Extracted covariate:", cov_names[i], "|", date(), "\n")
  
  rm(raster_file, extracted_values)
}

# Extract soil attributes from SoilGrids / WISE30sec data
soil_raster <- raster("E:/FACAI/Bigdata/Soil_WISE30sec/GISfiles/wise30sec_fin")

soil_extract <- extract(
  x = soil_raster,
  y = data[, .(LON, LAT)],
  method = "simple",
  df = TRUE
)

soil_attributes <- data.frame(soil_raster@data@attributes)[-1, ]
soil_codes <- soil_attributes[match(soil_extract[, 2], soil_attributes[, 9]), 7]

soil_table <- read.table(
  "E:/FACAI/Bigdata/Soil_WISE30sec/Interchangeable_format/HW30s_wD1.txt",
  sep = ",",
  header = TRUE
)

# Add soil covariates
data[, BD := soil_table[match(soil_codes, soil_table[, 1]), 24]]    # Bulk density
data[, Clay := soil_table[match(soil_codes, soil_table[, 1]), 21]]  # Clay content
data[, OCC := soil_table[match(soil_codes, soil_table[, 1]), 28]]   # Organic carbon content
data[, pH := soil_table[match(soil_codes, soil_table[, 1]), 34]]    # Soil pH in water
data[, EC := soil_table[match(soil_codes, soil_table[, 1]), 54]]    # Electrical conductivity
data[, CN := soil_table[match(soil_codes, soil_table[, 1]), 32]]    # Carbon-to-nitrogen ratio
data[, N2 := soil_table[match(soil_codes, soil_table[, 1]), 30]]    # Total nitrogen

# Convert hexagon centroid coordinates to spatial points
data_sf <- st_as_sf(
  data,
  coords = c("LON", "LAT"),
  crs = 4326,
  remove = FALSE
)

# Load WWF terrestrial biome shapefile
biome <- st_read("E:/FACAI/GlobalData/ecoregion/official/wwf_terr_ecos.shp")

# Reproject biome layer to WGS84
biome <- st_transform(biome, 4326)

# Retain only the WWF biome code
biome <- biome %>%
  dplyr::select(BIOME)

# Disable s2 geometry processing if needed for polygon joins
sf::sf_use_s2(FALSE)

# Assign WWF biome codes to each hexagon centroid
data_final <- st_join(data_sf, biome, left = TRUE)

# Check final data summary
summary(data_final)

# Export final hexagon-level biodiversity and covariate dataset
fwrite(
  data_final,
  "F:/GAIA_datapaper/Incidence3/NumObs/hexlevel_biodiversity_final.csv",
  row.names = FALSE
)



# ------------------------------------------------------------
# Standardize variable names and column structure for public release
# ------------------------------------------------------------                   

library(data.table)

rm(list = ls())

setwd("F:/GAIA_datapaper/Incidence3/NumObs")

# Load intermediate hexagon-level biodiversity dataset
dt <- fread("hexlevel_biodiversity_final.csv")

# Remove geometry column if present
if ("geometry" %in% names(dt)) {
  dt[, geometry := NULL]
}

# Rename variables for public release
rename_map <- c(
  aridity = "BC20",
  evapotspr = "BC21",
  aspectcosine = "AC",
  aspectsine = "AS",
  elevation = "EL",
  humanfp = "HF",
  pcurv = "PC",
  roughness = "RO",
  slope = "SL",
  tcurv = "TC",
  tpi = "TP",
  tri = "TR",
  Clay = "CL",
  OCC = "OC"
)

for (old_name in names(rename_map)) {
  if (old_name %in% names(dt)) {
    setnames(dt, old = old_name, new = rename_map[[old_name]])
  }
}

# Rename bioclimatic variables
for (i in 1:20) {
  
  old_name <- paste0("bio", i)
  new_name <- paste0("BC", i)
  
  if (old_name %in% names(dt)) {
    setnames(dt, old = old_name, new = new_name)
  }
}

# Define preferred final column order
preferred_order <- c(
  "hexID", "LAT", "LON", "N", "K", "R",
  "AC", "AS",
  paste0("BC", 1:21),
  "dx", "dxx", "dy", "dyy",
  "EL",
  "HF",
  "PC", "RO", "SL", "TC", "TP", "TR",
  "BD", "CL", "OC", "pH", "EC", "CN", "N2",
  "BIOME"
)

# Keep only columns that exist in the dataset
preferred_order <- preferred_order[preferred_order %in% names(dt)]

# Reorder columns
setcolorder(
  dt,
  c(preferred_order, setdiff(names(dt), preferred_order))
)

# Export finalized public-release dataset
fwrite(
  dt,
  "hexlevel_biodiversity.csv",
  row.names = FALSE
)