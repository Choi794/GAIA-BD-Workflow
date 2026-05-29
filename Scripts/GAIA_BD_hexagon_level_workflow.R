# Step 1. Generate species-level observation counts within each hexID
# Example workflow using a representative subset of the Plantae records.
# Because the complete Plantae dataset is extremely large, a sample dataset
# is used here to demonstrate the processing workflow. The same procedures
# were applied to the full dataset to generate the final GAIA-BD products.

library(data.table)

setwd(".")

rm(list = ls())

# Load species-level biodiversity observation data
dat <- fread("biodiversity_data_Plantae_all_sample.csv")

# Count the number of observations for each species within each hexID
result <- dat[
  ,
  .(NumObs = .N),
  by = .(hexID, kingdom, species)
]

# Preview the result
head(result)

# Export the summarized dataset
fwrite(result, "biodiversity_data_Plantae_all_NumObs.csv")



# Step 2. Merge kingdom-level observation summaries into an integrated biodiversity table

library(data.table)

rm(list = ls())

setwd(".")

# Load kingdom-level biodiversity observation summaries
dat1 <- fread("biodiversity_data_Animalia_all_NumObs.csv")
dat2 <- fread("biodiversity_data_Archaea_all_NumObs.csv")
dat3 <- fread("biodiversity_data_Bacteria_all_NumObs.csv")
dat4 <- fread("biodiversity_data_Chromista_all_NumObs.csv")
dat5 <- fread("biodiversity_data_Fungi_all_NumObs.csv")
dat6 <- fread("biodiversity_data_incertae sedis_all_NumObs.csv")
dat7 <- fread("biodiversity_data_Plantae_all_NumObs.csv")
dat8 <- fread("biodiversity_data_Protozoa_all_NumObs.csv")

# Combine all kingdom-level summaries
dat_total <- rbind(
  dat1, dat2, dat3, dat4,
  dat5, dat6, dat7, dat8
)

# Export integrated biodiversity observation summary
fwrite(dat_total, "biodiversity_data_all_NumObs.csv")



# Step 3. Generate hexagon-level biodiversity summary metrics

library(data.table)

rm(list = ls())

setwd(".")

# Load integrated species-level observation summary
dat <- fread("biodiversity_data_all_NumObs.csv")

# Check column names
names(dat)

# Calculate hexagon-level biodiversity metrics
result <- dat[
  ,
  .(
    # Total number of biodiversity observations
    N = sum(NumObs, na.rm = TRUE),
    
    # Number of taxonomic kingdoms
    K = uniqueN(kingdom),
    
    # Number of unique species
    R = uniqueN(species)
  ),
  by = hexID
]

# Preview the result
head(result)

# Export hexagon-level biodiversity summary
fwrite(result, "hexIDlevel_biodiversity.csv")



# Step 4. Extract centroid coordinates from the hexagon shapefile

library(sf)
library(data.table)

# Set working directory
setwd(".")

# Load the hexagon shapefile
hex <- st_read("hex0.5.shp")

# Check the coordinate reference system
st_crs(hex)

# Calculate centroid coordinates for each hexagon
hex_centroid <- st_centroid(hex)

# Extract centroid coordinates
centroid_coords <- st_coordinates(hex_centroid)

# Add centroid longitude and latitude to the hexagon attribute table
hex$centroid_LON <- centroid_coords[, 1]
hex$centroid_LAT <- centroid_coords[, 2]

# Remove geometry and export the centroid coordinate table
hex_centroid_table <- st_drop_geometry(hex)

fwrite(
  hex_centroid_table,
  "data/output/hex0.5_centroid_coordinates.csv"
)



# Step 5. Add hexagon centroid coordinates to hexagon-level biodiversity summaries

library(data.table)

rm(list = ls())

# Load hexagon centroid coordinates and biodiversity summary data
centroid <- fread("data/input/hex0.5_centroid_coordinates.csv")
bio <- fread("data/input/hexIDlevel_biodiversity.csv")

# Check column names
names(centroid)
names(bio)

# Select required coordinate columns and rename them
coord <- centroid[
  ,
  .(
    hexID,
    LAT = centroid_LAT,
    LON = centroid_LON
  )
]

# Merge centroid coordinates with biodiversity summaries by hexID
bio_with_coord <- merge(
  bio,
  coord,
  by = "hexID",
  all.x = TRUE
)

# Reorder columns so that LAT and LON follow hexID
setcolorder(
  bio_with_coord,
  c(
    "hexID",
    "LAT",
    "LON",
    setdiff(names(bio_with_coord), c("hexID", "LAT", "LON"))
  )
)

# Export hexagon-level biodiversity summary with centroid coordinates
fwrite(
  bio_with_coord,
  "hexIDlevel_biodiversity_with_centroid_coordinates.csv"
)



