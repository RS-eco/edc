`## EBBA2-Grid

# Load sf package
library(sf)

# Read data
ebba2_shp <- st_read("rawdata/ebba2_50x50_grid.shp")

# Plot data
plot(st_geometry(ebba2_shp))

# Save to file
save(ebba2_shp, file="data/ebba2_shp.rda", compress="xz")
`