library(sf)
library(data.table)
library(dplyr)

world <- st_read("data/input/continent_GFB2.shp")
world$CONTINENT <- NULL
world <- st_transform(world, crs = 4326)

# Hexagon data
hex <- st_make_grid(world, cellsize = c(0.5, 0.5), square = FALSE)
hex <- st_as_sf(hex)
hex$hexID <- 1:nrow(hex)
if(length(which(st_is_valid(hex) == FALSE)) > 0){
  hex <- hex[-which(st_is_valid(hex) == FALSE),]
}


# Save shapefile
st_write(hex, "data/output/hex0.5.shp", delete_layer = TRUE)