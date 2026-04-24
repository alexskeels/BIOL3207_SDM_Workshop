# =============================================================================
# BIOL3207 — Week 7 Workshop Solutions
# Frilled Lizard (Chlamydosaurus kingii) — spatial data, cleaning, geometric SDMs
#
# This script contains every code chunk from workshop_instructions_wk7.qmd
# plus the code needed to answer each "Ask yourself" question. Questions are
# included as comments, with short notes on plausible answers.
#
# Run from the project root (the folder containing data/).
# =============================================================================


# -----------------------------------------------------------------------------
# Packages
# -----------------------------------------------------------------------------
# install.packages(c(
#   "tidyverse", "knitr", "sf", "leaflet", "rnaturalearth",
#   "alphahull", "sp", "terra", "remotes"
# ))
# remotes::install_github("babichmorrowc/hull2spatial")

library(tidyverse)
library(knitr)
library(sf)
library(leaflet)
library(rnaturalearth)
library(alphahull)
library(hull2spatial)
library(sp)
library(terra)


# =============================================================================
# Task 1. Exploring the data
# =============================================================================
# Q: How much missing data is there?
# Q: What sorts of information are reliably recorded?
# Q: What sorts of information are often not recorded?
# Q: What shape is the distribution of missingness across columns?
#
# Notes for answers:
#  - There is a LOT of missing data overall: many columns are >90% NA
#    (a long tail of mostly-empty fields like measurement traits, identifier
#    remarks, sex, life stage, behaviour, habitat, depth, etc.).
#  - Reliably recorded: core taxonomy (kingdom -> species, scientificName),
#    coordinates (decimalLatitude/Longitude), eventDate/year, basisOfRecord,
#    dataResourceName, recordID. These are required by ALA's data pipeline.
#  - Often missing: ecological metadata (habitat, behaviour, sex, lifeStage),
#    measurement fields, identification remarks, georeferencing precision.
#  - The distribution is strongly RIGHT-SKEWED (a "J shape"): most variables
#    have few or no missing values, but a long tail of variables are nearly
#    100% missing. This is typical of biodiversity occurrence datasets.
# -----------------------------------------------------------------------------

# read csv
occ <- read_csv("data/frilled_lizard_ALA.csv")

# inspect
occ |>
  head() |>
  kable()

# how much missing data for each variable?
tibble(
  column      = names(occ),
  n_missing   = colSums(is.na(occ)),
  pct_missing = colMeans(is.na(occ)) * 100
) |>
  arrange(desc(pct_missing)) |>
  head(100) |>
  kable()

# Solution chunk: histogram of missingness
missingness_table <- tibble(
  column      = names(occ),
  n_missing   = colSums(is.na(occ)),
  pct_missing = colMeans(is.na(occ)) * 100
)

p_missing <- ggplot(missingness_table, aes(x = n_missing)) +
  geom_histogram() +
  theme_classic()
print(p_missing)


# =============================================================================
# Making Maps
# =============================================================================
# Q: Are records evenly distributed across the species range?
# Q: Do you see clustering near roads, cities, or coastlines?
# Q: Are some regions heavily sampled while others have very few records?
# Q: Why might these patterns occur?
#
# Notes for answers:
#  - Records are very UNEVEN. Strong clustering around Darwin, Cairns,
#    Townsville, Kakadu, and along major highways (Stuart Hwy, Bruce Hwy).
#  - Vast inland gaps in the Top End and Cape York where the species likely
#    occurs but human observers / road access are sparse.
#  - Drivers: sampling effort follows roads, settlements, national parks,
#    and tourist areas; this is "accessibility bias", not biological absence.
#  - Some offshore/coastal points hint at coordinate errors that we will
#    address in Task 2 with buffers.
# -----------------------------------------------------------------------------

# Load world coastline polygon (used later for cropping)
world <- ne_countries(scale = "medium", returnclass = "sf")

# Prepare the spatial data
occ <- occ |>
  filter(!is.na(decimalLongitude),
         !is.na(decimalLatitude))

occ_sf <- st_as_sf(
  occ,
  coords = c("decimalLongitude", "decimalLatitude"),
  crs    = 4326
)


##### !!!!! OPTIONAL - remove all points below -27 degrees latitude
occ_sf <- occ_sf[st_coordinates(occ_sf)[, "Y"] >= -27, ]


# Build an interactive map of occurrence records
map_occ <- leaflet(occ_sf) |>
  addProviderTiles(providers$OpenStreetMap.Mapnik) |>
  addCircleMarkers(
    radius      = 5,
    fillOpacity = 0.7,
    stroke      = FALSE,
    color       = "#2C7FB8"
  )
print(map_occ)


# =============================================================================
# Task 2. Clean occurrence records with buffers
# =============================================================================
# Q: Which dataset (point or polygon) likely overestimates the species range?
# Q: Which might underestimate it?
# Q: Where do you see obvious errors?
# Q: Are there plausible records just outside the polygon?
# Q: Why might some areas inside the polygon have no observations?
#
# Notes for answers:
#  - The IUCN polygon likely OVERESTIMATES: it is drawn as a smooth envelope
#    that includes unsuitable habitat (e.g. interior arid zones). Expert maps
#    are conservative outer bounds.
#  - The point data likely UNDERESTIMATES the realised range: many remote
#    parts of the Top End and Cape York are under-sampled (see Making Maps).
#  - Obvious errors: any points in the ocean, in southern Australia, or
#    overseas (museum collections geocoded to institutions, etc.).
#  - Plausible edge records: clusters just outside the southern/western
#    boundary of the polygon — these likely reflect real range extension or
#    polygon being too coarse.
#  - Areas inside the polygon with no observations usually = sampling gaps
#    (no roads, no people), not absence.
# -----------------------------------------------------------------------------

# Load the IUCN spatial polygon (shapefile)
iucn_sf <- st_read("data/frilled_lizard.shp")

# Map: occurrences + IUCN polygon
map_occ_iucn <- leaflet(occ_sf) |>
  addProviderTiles(providers$OpenStreetMap.Mapnik) |>
  addCircleMarkers(
    radius = 5, fillOpacity = 0.7, stroke = FALSE, color = "#2C7FB8"
  ) |>
  addPolygons(
    data = iucn_sf,
    fillColor   = "#009E73",
    color       = "#007F5F",
    weight      = 1,
    fillOpacity = 0.4,
    popup       = "IUCN Distribution"
  )
print(map_occ_iucn)


# Cropping to the Expert Range — buffers
# Q: Which buffer removes obvious errors?
# Q: Which keeps plausible edge populations?
# Q: At what point does the buffer become biologically unrealistic?
#
# Notes for answers:
#  - 50 km removes far-flung outliers but is tight enough to exclude some
#    plausible edge records.
#  - 100–200 km is a good compromise: keeps believable edge populations and
#    accounts for IUCN polygon coarseness.
#  - 500 km is biologically unrealistic — it sweeps in much of arid central
#    Australia and the southern coast where the species clearly does not occur.
#  - Default choice below: 200 km (matches the qmd solution chunk).
# -----------------------------------------------------------------------------

# Project to an Australian equal-area CRS (units = metres)
iucn_proj <- st_transform(iucn_sf, crs = 3577)
occ_proj  <- st_transform(occ_sf,  crs = 3577)

# Worked example: 50 km buffer
buffer_50km <- st_buffer(iucn_proj, dist = 50000)
buffer_50km <- st_transform(buffer_50km, crs = 4326)

sf_use_s2(FALSE)
buffer_50km <- st_intersection(st_make_valid(buffer_50km), world)
sf_use_s2(TRUE)

# Solution: build several buffers
buffer_100km <- st_buffer(iucn_proj, dist = 100000)
buffer_200km <- st_buffer(iucn_proj, dist = 285000)
buffer_500km <- st_buffer(iucn_proj, dist = 500000)

buffer_100km <- st_transform(buffer_100km, 4326)
buffer_200km <- st_transform(buffer_200km, 4326)
buffer_500km <- st_transform(buffer_500km, 4326)

sf_use_s2(FALSE)
buffer_100km <- st_intersection(st_make_valid(buffer_100km), world)
buffer_200km <- st_intersection(st_make_valid(buffer_200km), world)
buffer_500km <- st_intersection(st_make_valid(buffer_500km), world)
sf_use_s2(TRUE)

# Plot all buffers together
map_buffers <- leaflet(occ_sf) |>
  addProviderTiles(providers$OpenStreetMap.Mapnik) |>
  addCircleMarkers(
    radius = 5, fillOpacity = 0.7, stroke = FALSE, color = "#2C7FB8"
  ) |>
  addPolygons(
    data = iucn_sf,
    fillColor = "#009E73", color = "#007F5F",
    weight = 2, fillOpacity = 0.3, popup = "IUCN Range"
  ) |>
  addPolygons(
    data = buffer_100km,
    fillColor = "#FFD700", color = "#B8860B",
    weight = 1, fillOpacity = 0.1, popup = "100 km buffer"
  ) |>
  addPolygons(
    data = buffer_200km,
    fillColor = "#FFA500", color = "#FF8C00",
    weight = 1, fillOpacity = 0.1, popup = "200 km buffer"
  ) |>
  addPolygons(
    data = buffer_500km,
    fillColor = "#FF4500", color = "#B22222",
    weight = 1, fillOpacity = 0.1, popup = "500 km buffer"
  )
print(map_buffers)

# Pick the preferred buffer — used downstream
buffer_iucn <- buffer_200km

# Subset occurrences to those inside the chosen buffer
occ_cropped <- occ_sf[lengths(st_within(occ_sf, buffer_iucn)) > 0, ]

# Save cleaned data for next week's prac
occ_cropped_df <- as_tibble(st_coordinates(occ_cropped))
write_csv(occ_cropped_df, file = "data/frilled_lizard_ALA_cropped.csv")


# =============================================================================
# Task 3. Build geometric SDMs
# =============================================================================
# Notes:
#  - A convex hull is the smallest convex polygon containing all points. It
#    is fast to compute but tends to OVERESTIMATE range because it cannot
#    "carve out" unsuitable interior gaps.
#  - An alpha hull allows concavity (controlled by alpha): smaller alpha =
#    tighter, more concave; larger alpha -> approaches the convex hull.
# -----------------------------------------------------------------------------

# Convex hull
coords <- st_coordinates(occ_cropped)
lon <- coords[, "X"]
lat <- coords[, "Y"]

hull_indices <- chull(lon, lat)
hull_indices <- c(hull_indices, hull_indices[1])  # close polygon

hull_coords  <- cbind(lon[hull_indices], lat[hull_indices])
hull_polygon <- st_polygon(list(hull_coords))
hull_sf      <- st_sfc(hull_polygon, crs = 4326)
hull_sf      <- st_sf(geometry = hull_sf)

sf_use_s2(FALSE)
hull_sf <- st_intersection(st_make_valid(hull_sf), world)
sf_use_s2(TRUE)


# Alpha hull
alpha <- 3

xy_df <- data.frame(x = lon, y = lat)
xy_df <- unique(round(xy_df, 4))

ahull_obj <- ahull(xy_df, alpha = alpha)

# Convert ahull -> sp -> sf (boilerplate from the qmd)
ahull_spatial <- ahull2poly(ahull_obj)
sub_polys     <- slot(ahull_spatial, "polygons")

poly_list <- lapply(seq_along(sub_polys), function(i) {
  SpatialPolygons(list(sub_polys[[i]]),
                  proj4string = CRS(proj4string(ahull_spatial)))
})
sf_list  <- lapply(poly_list, st_as_sf)
ahull_sf <- do.call(rbind, sf_list)
ahull_sf$id <- seq_len(nrow(ahull_sf))
st_crs(ahull_sf) <- 4326

sf_use_s2(FALSE)
ahull_sf <- st_intersection(st_make_valid(ahull_sf), world)
sf_use_s2(TRUE)


# =============================================================================
# Task 4. Compare all layers
# =============================================================================
# Q: Which method (buffer, convex hull, alpha hull) best represents the
#    known distribution?
# Q: How do assumptions differ between layers?
# Q: What are the trade-offs between simplicity (convex hull) and
#    realism (alpha hull)?
#
# Notes for answers:
#  - Best fit: alpha hull is usually the closest match to the IUCN polygon
#    because it can carve out unsampled interior gaps, but it is sensitive
#    to alpha and to sampling bias (under-sampled regions get excluded).
#  - The IUCN polygon is built from expert knowledge (smooth, conservative,
#    can include unsuitable habitat).
#  - The buffer assumes the IUCN polygon is approximately correct and just
#    pads it; useful for cleaning, not for SDM prediction.
#  - The convex hull assumes the range is convex (rarely true for real
#    species) and ignores interior gaps -> overestimates.
#  - Trade-off: convex hull is simple and stable but biologically naive;
#    alpha hull is more realistic but depends on alpha and on having dense,
#    unbiased sampling.
# -----------------------------------------------------------------------------

map_compare <- leaflet(occ_cropped) |>
  addProviderTiles(providers$OpenStreetMap.Mapnik) |>
  addCircleMarkers(
    radius = 5, fillOpacity = 0.7, stroke = FALSE, color = "#2C7FB8"
  ) |>
  addPolygons(
    data = iucn_sf,
    fillColor = "#009E73", color = "#007F5F",
    weight = 2, fillOpacity = 0.1, popup = "IUCN Range"
  ) |>
  addPolygons(
    data = hull_sf,
    fillColor = "#FFD700", color = "#B8860B",
    weight = 1, fillOpacity = 0.1, popup = "convex hull"
  ) |>
  addPolygons(
    data = ahull_sf,
    fillColor = "#FFA500", color = "#FF8C00",
    weight = 1, fillOpacity = 0.1, popup = "alpha hull"
  ) |>
  addPolygons(
    data = buffer_iucn,
    fillColor = "#FF4500", color = "#B22222",
    weight = 1, fillOpacity = 0.1, popup = "your buffer"
  )
print(map_compare)


#### Cleaner map with just alpha and IUCN

map_compare <- leaflet(occ_cropped) |>
  addProviderTiles(providers$OpenStreetMap.Mapnik) |>
  addCircleMarkers(
    radius = 5, fillOpacity = 0.7, stroke = FALSE, color = "#2C7FB8"
  ) |>
  addPolygons(
    data = iucn_sf,
    fillColor = "#009E73", color = "#007F5F",
    weight = 2, fillOpacity = 0.1, popup = "IUCN Range"
  ) |>
  addPolygons(
    data = ahull_sf,
    fillColor = "#FFA500", color = "#FF8C00",
    weight = 1, fillOpacity = 0.1, popup = "alpha hull"
  ) 


##### Map with included and excluded points

occ_excluded <- occ_sf[lengths(st_within(occ_sf, buffer_iucn)) == 0, ]

leaflet() |>
  
  # Add a detailed basemap
  addProviderTiles(providers$OpenStreetMap.Mapnik) |>
  
  # Excluded occurrence records (outside buffer)
  addCircleMarkers(
    data = occ_excluded,
    radius = 5,
    fillOpacity = 0.8,
    stroke = FALSE,
    color = "#D7301F",   # red
    group = "Excluded points"
  ) |>
  
# Included occurrence records (inside buffer = occ_cropped)
addCircleMarkers(
  data = occ_cropped,
  radius = 5,
  fillOpacity = 0.8,
  stroke = FALSE,
  color = "#2C7FB8",   # blue
  group = "Included points"
) |>
  
  # IUCN polygon
  addPolygons(
    data = iucn_sf,
    fillColor = "#009E73",
    color = "#007F5F",
    weight = 2,
    fillOpacity = 0.1,
    popup = "IUCN Range"
  ) |>
# Convex hull
addPolygons(
  data = hull_sf,
  fillColor = "#FFD700",
  color = "#B8860B",
  weight = 1,
  fillOpacity = 0.1,
  popup = "convex hull"
) |>
  
  # Alpha hull
  addPolygons(
    data = ahull_sf,
    fillColor = "#FFA500",
    color = "#FF8C00",
    weight = 1,
    fillOpacity = 0.1,
    popup = "alpha hull"
  ) |>

# Your chosen buffer (from Task 4)
addPolygons(
  data = buffer_iucn,
  fillColor = "#FF4500",
  color = "#B22222",
  weight = 1,
  fillOpacity = 0.1,
  popup = "your buffer"
) |>
  
  # Legend for points
  addLegend(
    position = "bottomright",
    colors   = c("#2C7FB8", "#D7301F"),
    labels   = c("Included (occ_cropped)", "Excluded"),
    title    = "Occurrence points",
    opacity  = 0.8
  )

print(map_compare)

# Save the alpha hull for next week
sf::st_write(ahull_sf, "data/frilled_lizard_alpha_hull.shp", delete_layer = TRUE)
# If sf write fails, fall back to terra:
# terra::writeVector(terra::vect(ahull_sf),
#                    "data/frilled_lizard_alpha_hull.shp",
#                    overwrite = TRUE)


# =============================================================================
# Calculate Range Overlap With The IUCN (just like in Aronsson reading)
# =============================================================================

# Project to an equal-area CRS (Australia Albers, EPSG:3577 — good for AU)
iucn_eq  <- st_transform(iucn_sf,  3577)
ahull_eq <- st_transform(ahull_sf, 3577)
hull_eq <- st_transform(hull_sf, 3577)

# Make geometries valid before intersecting (avoids topology errors)
iucn_eq  <- st_make_valid(iucn_eq)
ahull_eq <- st_make_valid(ahull_eq)

# Areas
overlap       <- st_intersection(iucn_eq, ahull_eq)
overlap_area  <- sum(st_area(overlap))
iucn_area     <- sum(st_area(iucn_eq))
ahull_area    <- sum(st_area(ahull_eq))
hull_area    <- sum(st_area(hull_eq))

# Overlap / total area of the two polygons
union_area    <- iucn_area + ahull_area - overlap_area
overlap_index <- as.numeric(overlap_area / union_area)
overlap_index |> round(3)

# Overlap / total area of the two polygons
union_area    <- iucn_area + hull_area - overlap_area
overlap_index <- as.numeric(overlap_area / union_area)
overlap_index |> round(3)







# =============================================================================
# Going Further (1) — Range change through time (pre/post 2000)
# =============================================================================
# Q: Do we see a difference in the geographic range prior to and post 2000?
#
# Notes for answers:
#  - Historical (pre-2000) records are sparser and more concentrated near
#    long-established settlements / museum collection localities.
#  - Contemporary (2000+) records are far more numerous and broader,
#    largely driven by citizen science (e.g. iNaturalist into ALA) rather
#    than a true biological range expansion.
#  - Apparent "range expansion" therefore mostly reflects increased
#    sampling effort, not necessarily a shifting distribution. Caveat any
#    inference about climate-driven range shifts accordingly.
# -----------------------------------------------------------------------------

# Split the cleaned (cropped) occurrences into pre/post 2000
occ_hist <- occ_cropped |> filter(!is.na(year), year <  2000)
occ_curr <- occ_cropped |> filter(!is.na(year), year >= 2000)

cat("Pre-2000 records: ", nrow(occ_hist), "\n",
    "2000+  records: ",  nrow(occ_curr), "\n", sep = "")

# ----- Pre-2000 alpha hull --------------------------------------------------
coords_hist <- st_coordinates(occ_hist)
xy_hist <- unique(round(data.frame(x = coords_hist[, "X"],
                                   y = coords_hist[, "Y"]), 4))

ahull_hist_obj <- ahull(xy_hist, alpha = 5)
ahull_hist_sp  <- ahull2poly(ahull_hist_obj)
polys_hist     <- slot(ahull_hist_sp, "polygons")

sf_list_hist <- lapply(seq_along(polys_hist), function(i) {
  st_as_sf(SpatialPolygons(list(polys_hist[[i]]),
                           proj4string = CRS(proj4string(ahull_hist_sp))))
})
ahull_hist <- do.call(rbind, sf_list_hist)
st_crs(ahull_hist) <- 4326

sf_use_s2(FALSE)
ahull_hist <- st_intersection(st_make_valid(ahull_hist), world)
sf_use_s2(TRUE)

# ----- 2000+ alpha hull -----------------------------------------------------
coords_curr <- st_coordinates(occ_curr)
xy_curr <- unique(round(data.frame(x = coords_curr[, "X"],
                                   y = coords_curr[, "Y"]), 4))

ahull_curr_obj <- ahull(xy_curr, alpha = 5)
ahull_curr_sp  <- ahull2poly(ahull_curr_obj)
polys_curr     <- slot(ahull_curr_sp, "polygons")

sf_list_curr <- lapply(seq_along(polys_curr), function(i) {
  st_as_sf(SpatialPolygons(list(polys_curr[[i]]),
                           proj4string = CRS(proj4string(ahull_curr_sp))))
})
ahull_curr <- do.call(rbind, sf_list_curr)
st_crs(ahull_curr) <- 4326

sf_use_s2(FALSE)
ahull_curr <- st_intersection(st_make_valid(ahull_curr), world)
sf_use_s2(TRUE)

# ----- Compare on a leaflet map --------------------------------------------
map_time <- leaflet() |>
  addProviderTiles(providers$OpenStreetMap.Mapnik) |>
  addPolygons(data = ahull_hist,
              fillColor = "#1F78B4", color = "#1F78B4",
              weight = 1, fillOpacity = 0.2,
              popup = "Pre-2000 alpha hull") |>
  addPolygons(data = ahull_curr,
              fillColor = "#E31A1C", color = "#E31A1C",
              weight = 1, fillOpacity = 0.2,
              popup = "2000+ alpha hull")

print(map_time)

# =============================================================================
# End of Week 7 solutions
# =============================================================================
