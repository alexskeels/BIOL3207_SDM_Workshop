# =============================================================================
# BIOL3207 — Week 8 Workshop Solutions
# Frilled Lizard (Chlamydosaurus kingii) — statistical SDMs (GLM + GAM)
#
# This script contains every code chunk from workshop_instructions_wk8.qmd
# and workshop_instructions_wk8_extended.qmd, with answers to the
# "Ask yourself" questions written as inline comments.
#
# Unlike the workshop (where each student does ONE background method), this
# script runs all FOUR background methods so you can compare them directly.
# Object names use suffixes so nothing is overwritten:
#   _rand  = random background
#   _acc   = accessible-area buffer background
#   _env   = environmental envelope background
#   _bias  = bias-corrected (Agamid effort) background
#
# Run from the project root (the folder containing data/).
# =============================================================================


# -----------------------------------------------------------------------------
# Packages
# -----------------------------------------------------------------------------
# install.packages(c(
#   "tidyverse", "terra", "sf", "raster",
#   "skimr", "leaflet", "leafsync",
#   "caret", "mgcv", "dismo"
# ))

library(tidyverse)
library(terra)
library(sf)
library(raster)     # needed for leaflet raster display
library(skimr)
library(leaflet)
library(leafsync)
library(caret)
library(mgcv)       # for GAMs
library(dismo)      # for evaluate(), threshold(), kfold()


# -----------------------------------------------------------------------------
# Shared inputs — occurrence points, alpha hull, and the three climate stacks
# -----------------------------------------------------------------------------

# Cropped occurrence points written at the end of week 7 (columns: X, Y)
occ <- read_csv("data/frilled_lizard_ALA_cropped.csv")

# Alpha hull written at the end of week 7
ahull_v  <- vect("data/frilled_lizard_alpha_hull.shp")
ahull_sf <- st_as_sf(ahull_v)

# Present-day climate (WorldClim, Australia) — Bio01..Bio19
bioclim_data <- rast("data/worldclim_australia.tif")

# Future climate (WorldClim, MIROC6 SSP585, 2100)
future_bioclim_data <- rast("data/worldclim_australia_MIROC6_2100.tif")

# Last Glacial Maximum (CHELSA, ~21,000 years ago)
lgm_bioclim_data <- rast("data/chelsa_australia_LGM.tif")

# Each stack stores Bio01..Bio19 but in DIFFERENT band orders:
#   - WorldClim present is alphabetical (Bio01, Bio10, Bio11, ..., Bio09)
#   - WorldClim CMIP6 future + CHELSA LGM are numerical (Bio01, Bio02, ..., Bio19)
# Re-order all three to numerical Bio01..Bio19 so that band i means the same
# variable in every stack. This is essential — terra::predict matches
# predictors to model terms by NAME, so any mismatch silently feeds the wrong
# variable into the model and the LGM/2100 predictions collapse to zero.
bio_order           <- sprintf("Bio%02d", 1:19)
bioclim_data        <- bioclim_data[[bio_order]]
future_bioclim_data <- future_bioclim_data[[bio_order]]
lgm_bioclim_data    <- lgm_bioclim_data[[bio_order]]


# =============================================================================
# Task 1. Inspect the predictors
# =============================================================================
# Q: Can you guess what units the variables are in?
# Q: Do variable ranges look reasonable?
# Q: Do any of the variables appear correlated?
#
# Notes for answers:
#  - WorldClim units: temperatures (Bio01–Bio11) are °C. Precipitation
#    (Bio12–Bio19) is mm.
#  - Ranges look reasonable for Australia: Bio01 (annual mean temp) ~10–30 °C,
#    Bio12 (annual precipitation) ~50–4000 mm.
#  - Many bioclim variables ARE highly correlated by construction — e.g.
#    Bio01/Bio05/Bio06 all describe temperature levels; Bio12/Bio13/Bio16 all
#    describe wet-season rainfall. We address this in Task 3.
# -----------------------------------------------------------------------------

values(bioclim_data) |>
  as_tibble() |>
  skim()

# Define one palette per layer and reuse it for BOTH the raster image and the
# legend, so the colours on the map match the colours in the legend.
# "Spectral" reversed gives blue (cold/dry) -> red (hot/wet) — intuitive for
# temperature, and matches the Spectral palette used later for suitability.
pal_bio01 <- colorNumeric("Spectral", values(bioclim_data$Bio01),
                          na.color = "transparent", reverse = TRUE)
pal_bio12 <- colorNumeric("Spectral", values(bioclim_data$Bio12),
                          na.color = "transparent", reverse = TRUE)

bio01 <- leaflet() |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addRasterImage(raster::raster(bioclim_data$Bio01),
                 colors = pal_bio01, opacity = 0.8) |>
  addLegend(pal    = pal_bio01,
            values = values(bioclim_data$Bio01),
            title  = "Bio01 (°C)")

bio12 <- leaflet() |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addRasterImage(raster::raster(bioclim_data$Bio12),
                 colors = pal_bio12, opacity = 0.8) |>
  addLegend(pal    = pal_bio12,
            values = values(bioclim_data$Bio12),
            title  = "Bio12 (mm)")

sync(bio01, bio12)


# Same comparison for two seasonality variables:
#   Bio04 = temperature seasonality (SD of monthly temps × 100)
#   Bio15 = precipitation seasonality (CV of monthly precip)
pal_bio04 <- colorNumeric("Spectral", values(bioclim_data$Bio04),
                          na.color = "transparent", reverse = TRUE)
pal_bio15 <- colorNumeric("Spectral", values(bioclim_data$Bio15),
                          na.color = "transparent", reverse = TRUE)

bio04 <- leaflet() |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addRasterImage(raster::raster(bioclim_data$Bio04),
                 colors = pal_bio04, opacity = 0.8) |>
  addLegend(pal    = pal_bio04,
            values = values(bioclim_data$Bio04),
            title  = "Bio04 (temp SD ×100)")

bio15 <- leaflet() |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addRasterImage(raster::raster(bioclim_data$Bio15),
                 colors = pal_bio15, opacity = 0.8) |>
  addLegend(pal    = pal_bio15,
            values = values(bioclim_data$Bio15),
            title  = "Bio15 (precip CV)")

sync(bio04, bio15)


# -----------------------------------------------------------------------------
# Build the presence dataframe (one record per grid cell)
# -----------------------------------------------------------------------------

# SpatVector from occurrence csv
occ_v <- vect(occ, geom = c("X", "Y"), crs = crs(bioclim_data))

# Keep one point per raster cell (de-duplicates spatial sampling bias)
occ_v$cell_id <- cellFromXY(bioclim_data, crds(occ_v))
occ_pruned    <- occ_v[!duplicated(occ_v$cell_id), ]

# Extract climate at presences
occ_env <- terra::extract(bioclim_data, occ_pruned, ID = FALSE) |>
  as_tibble()

occ_model <- bind_cols(
  as_tibble(crds(occ_pruned)) |>
    setNames(c("longitude", "latitude")),
  occ_env
) |>
  mutate(presence = 1)


# =============================================================================
# Task 2. Generate background points — ALL four methods
# =============================================================================
# Each section produces a uniquely-named bg_model_<method> dataframe with
# columns longitude, latitude, Bio01..Bio19, presence = 0.
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Method 1: Random background (bg_model_rand)
# -----------------------------------------------------------------------------
# Assumption: no geographic or environmental restriction. Backgrounds are
# uniformly sampled across the entire Australian raster mask.
# Effect on the model: produces broad, generic backgrounds. The model is
# trained to separate presences from "anywhere in Australia" — usually leads
# to over-confident, overly-broad predicted ranges.

set.seed(123)

bg_rand <- spatSample(
  bioclim_data,
  size      = 5000,
  method    = "random",
  na.rm     = TRUE,
  as.points = TRUE
)

bg_env_rand <- terra::extract(bioclim_data, bg_rand, ID = FALSE) |>
  as_tibble()

bg_model_rand <- bind_cols(
  st_coordinates(st_as_sf(bg_rand)) |>
    as_tibble() |>
    rename(longitude = X, latitude = Y),
  bg_env_rand
) |>
  mutate(presence = 0)

# Quick visual check
leaflet() |>
  addProviderTiles(providers$OpenStreetMap.Mapnik) |>
  addRasterImage(raster::raster(bioclim_data$Bio01), opacity = 0.5) |>
  addCircleMarkers(data = bg_model_rand,
                   lng = ~longitude, lat = ~latitude,
                   radius = 3, fillOpacity = 0.4,
                   stroke = FALSE, color = "grey") |>
  addCircleMarkers(data = occ_model,
                   lng = ~longitude, lat = ~latitude,
                   radius = 5, fillOpacity = 0.7,
                   stroke = FALSE, color = "#2C7FB8")


# -----------------------------------------------------------------------------
# Method 2: Accessible-area background (bg_model_acc)
# -----------------------------------------------------------------------------
# Assumption: species can only realistically reach areas within ~300 km of
# known occurrences. We mask the raster by a buffer around the points and
# sample within that buffer.
# Effect on the model: tightens the predicted distribution towards the known
# range. Useful when dispersal is genuinely limiting.

occ_sf   <- st_as_sf(occ, coords = c("X", "Y"), crs = 4326)
occ_proj <- st_transform(occ_sf, 3577)        # equal-area CRS, units = m

buffer_acc      <- st_buffer(occ_proj, dist = 300000) |> st_union()
buffer_acc_vect <- vect(st_transform(buffer_acc, 4326))

set.seed(123)
bg_acc <- spatSample(
  mask(bioclim_data, buffer_acc_vect),
  size      = 5000,
  method    = "random",
  na.rm     = TRUE,
  as.points = TRUE
)

bg_env_acc <- terra::extract(bioclim_data, bg_acc, ID = FALSE) |>
  as_tibble()

bg_model_acc <- bind_cols(
  st_coordinates(st_as_sf(bg_acc)) |>
    as_tibble() |>
    rename(longitude = X, latitude = Y),
  bg_env_acc
) |>
  mutate(presence = 0)

# Quick visual check
leaflet() |>
  addProviderTiles(providers$OpenStreetMap.Mapnik) |>
  addRasterImage(raster::raster(bioclim_data$Bio01), opacity = 0.5) |>
  addCircleMarkers(data = bg_model_acc,
                   lng = ~longitude, lat = ~latitude,
                   radius = 3, fillOpacity = 0.4,
                   stroke = FALSE, color = "darkgreen") |>
  addCircleMarkers(data = occ_model,
                   lng = ~longitude, lat = ~latitude,
                   radius = 5, fillOpacity = 0.7,
                   stroke = FALSE, color = "#2C7FB8")


# -----------------------------------------------------------------------------
# Method 3: Environmental envelope background (bg_model_env)
# -----------------------------------------------------------------------------
# Assumption: backgrounds should occupy environments similar to those at
# presences — sampling within (min - 1 SD) to (max + 1 SD) for every bioclim
# variable measured at presence locations.
# Effect on the model: sharpens the climatic contrast (presences vs nearby
# climates), which can make the model look very confident, but is structurally
# biased toward saying "yes" inside the envelope.

env_ranges <- occ_env |>
  summarise(
    across(
      everything(),
      list(min = ~min(.x, na.rm = TRUE) - sd(.x, na.rm = TRUE),
           max = ~max(.x, na.rm = TRUE) + sd(.x, na.rm = TRUE))
    )
  )

clim_df <- terra::as.data.frame(bioclim_data, xy = TRUE, na.rm = TRUE)

ranges_long <- env_ranges |>
  pivot_longer(everything(),
               names_to  = c("variable", ".value"),
               names_sep = "_(?=[^_]+$)")

filtered_cells <- clim_df |>
  pivot_longer(cols = -c(x, y), names_to = "variable", values_to = "value") |>
  left_join(ranges_long, by = "variable") |>
  filter(value >= min, value <= max) |>
  group_by(x, y) |>
  summarise(n = n(), .groups = "drop") |>
  filter(n == nrow(ranges_long)) |>
  dplyr::select(x, y)   # dplyr:: needed because raster masks select()

set.seed(123)
bg_env_pts <- filtered_cells |> slice_sample(n = 5000)

bg_env_vect <- vect(as.data.frame(bg_env_pts),
                    geom = c("x", "y"),
                    crs  = crs(bioclim_data))
bg_env_env  <- terra::extract(bioclim_data, bg_env_vect, ID = FALSE) |>
  as_tibble()

bg_model_env <- bind_cols(
  bg_env_pts |> rename(longitude = x, latitude = y),
  bg_env_env
) |>
  mutate(presence = 0)

# Quick visual check
leaflet() |>
  addProviderTiles(providers$OpenStreetMap.Mapnik) |>
  addRasterImage(raster::raster(bioclim_data$Bio01), opacity = 0.5) |>
  addCircleMarkers(data = bg_model_env,
                   lng = ~longitude, lat = ~latitude,
                   radius = 3, fillOpacity = 0.4,
                   stroke = FALSE, color = "darkorange") |>
  addCircleMarkers(data = occ_model,
                   lng = ~longitude, lat = ~latitude,
                   radius = 5, fillOpacity = 0.7,
                   stroke = FALSE, color = "#2C7FB8")


# -----------------------------------------------------------------------------
# Method 4: Bias-corrected (Agamid effort) background (bg_model_bias)
# -----------------------------------------------------------------------------
# Assumption: where Agamid lizards (the family containing the Frilled Lizard)
# have been recorded, people HAVE looked. Sampling backgrounds proportional to
# Agamid density conditions on observer effort, so the model contrasts
# presences with "places that were searched but where the species was not
# found" rather than with arbitrary unsampled wilderness.
# Effect on the model: usually the most defensible. Tightens range to
# accessible, surveyed regions and avoids treating unsurveyed land as absence.

agamid_occ <- read.csv("data/agamidae.csv")

agamid_occ <- terra::vect(
  agamid_occ,
  geom = c("decimalLongitude", "decimalLatitude"),
  crs  = "EPSG:4326"
)

agamid_density <- terra::rasterize(
  x   = agamid_occ,
  y   = bioclim_data[[1]],
  fun = "sum"
)

agamid_density[is.na(agamid_density)] <- 0
agamid_weights <- agamid_density / max(values(agamid_density), na.rm = TRUE)
agamid_weights <- mask(agamid_weights, bioclim_data[[1]])

set.seed(123)
bck_bias <- terra::spatSample(
  x         = agamid_weights,
  size      = 5000,
  method    = "weights",
  as.points = TRUE,
  na.rm     = TRUE
)

bg_bias <- as_tibble(terra::crds(bck_bias)) |>
  rename(longitude = x, latitude = y)

bg_bias_vect <- vect(bg_bias,
                     geom = c("longitude", "latitude"),
                     crs  = crs(bioclim_data))
bg_env_bias  <- terra::extract(bioclim_data, bg_bias_vect, ID = FALSE) |>
  as_tibble()

bg_model_bias <- bind_cols(bg_bias, bg_env_bias) |>
  mutate(presence = 0)

# Quick visual check
leaflet() |>
  addProviderTiles(providers$OpenStreetMap.Mapnik) |>
  addRasterImage(raster::raster(bioclim_data$Bio01), opacity = 0.5) |>
  addCircleMarkers(data = bg_model_bias,
                   lng = ~longitude, lat = ~latitude,
                   radius = 3, fillOpacity = 0.4,
                   stroke = FALSE, color = "purple") |>
  addCircleMarkers(data = occ_model,
                   lng = ~longitude, lat = ~latitude,
                   radius = 5, fillOpacity = 0.7,
                   stroke = FALSE, color = "#2C7FB8")


# Faceted ggplot showing all four background clouds side by side, with
# the presence points (red) overlaid for context.
bg_all <- bind_rows(
  bg_model_rand |> mutate(method = "Random"),
  bg_model_acc  |> mutate(method = "Accessible area (300 km)"),
  bg_model_env  |> mutate(method = "Environmental envelope"),
  bg_model_bias |> mutate(method = "Bias-corrected (Agamid)")
) |>
  mutate(method = factor(method, levels = c(
    "Random", "Accessible area (300 km)",
    "Environmental envelope", "Bias-corrected (Agamid)"
  )))

# Australia outline derived from the climate raster mask (no extra deps)
aus_mask    <- !is.na(bioclim_data[[1]])
names(aus_mask) <- "is_land"
aus_outline <- terra::as.polygons(aus_mask, dissolve = TRUE) |>
  sf::st_as_sf() |>
  dplyr::filter(is_land == 1)

ggplot() +
  geom_sf(data = aus_outline,
          fill = "grey96", colour = "grey60", linewidth = 0.2) +
  geom_point(data = bg_all,
             aes(x = longitude, y = latitude, colour = method),
             size = 0.5, alpha = 0.4) +
  geom_point(data = occ_model,
             aes(x = longitude, y = latitude),
             colour = "red", size = 0.5, alpha = 0.7) +
  facet_wrap(~ method, nrow = 2) +
  scale_colour_manual(values = c(
    "Random"                   = "grey30",
    "Accessible area (300 km)" = "darkgreen",
    "Environmental envelope"   = "darkorange",
    "Bias-corrected (Agamid)"  = "purple"
  ), guide = "none") +
  coord_sf() +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank()) +
  labs(x = NULL, y = NULL,
       title    = "Background point clouds — four methods",
       subtitle = "Red points = presences (occ_model)")


# =============================================================================
# Helper — drop highly correlated bioclim variables (used in every block below)
# =============================================================================
# Pre-defined here once so the four method blocks below stay easy to read.
# Each block calls drop_correlated() with its own sdm_data so the retained
# predictor set is method-specific (the bg cloud changes which bioclim pairs
# correlate).

drop_correlated <- function(sdm_data, cutoff = 0.7) {
  bio_cols <- grep("^Bio", colnames(sdm_data), value = TRUE)
  cor_mat  <- cor(sdm_data[, bio_cols], use = "complete.obs")
  drop_idx <- findCorrelation(cor_mat, cutoff = cutoff)
  drop_nm  <- bio_cols[drop_idx]
  list(
    data           = sdm_data[, !colnames(sdm_data) %in% drop_nm],
    predictor_vars = setdiff(bio_cols, drop_nm),
    dropped        = drop_nm
  )
}


# =============================================================================
# Task 3–5 (block 1 of 4): Random background pipeline
# =============================================================================
# Q: Which variables were removed due to colinearity?
# Q: Do the remaining predictors make ecological sense for the Frilled Lizard?
# Q: Which model has lower AIC / higher D² ?
#
# Notes for answers (this paragraph applies to all four blocks):
#  - The set of removed variables differs slightly between background methods,
#    because the bg cloud changes the joint distribution. This is one of the
#    quietly-important effects of background choice.
#  - Variables that typically survive: Bio04 (temp seasonality), Bio15 (precip
#    seasonality), Bio02 (mean diurnal range). These make ecological sense
#    for a tropical lizard whose distribution is bounded by seasonality.
#  - GAM almost always has lower AIC and higher deviance explained than GLM
#    on the training data, but that is partly over-fitting — see CV below.
# -----------------------------------------------------------------------------

sdm_data_rand   <- bind_rows(occ_model, bg_model_rand)

drop_rand                <- drop_correlated(sdm_data_rand)
sdm_data_uncor_rand      <- drop_rand$data
predictor_vars_rand      <- drop_rand$predictor_vars
dropped_rand             <- drop_rand$dropped
print(predictor_vars_rand)
print(dropped_rand)

formula_glm_rand <- formula(paste0("presence ~ ",
                                   paste0(predictor_vars_rand, collapse = " + ")))
formula_gam_rand <- formula(paste0("presence ~ ",
                                   paste0("s(", predictor_vars_rand, ", k = 4)",
                                          collapse = " + ")))

glm_sdm_rand <- glm(formula_glm_rand,
                    data = sdm_data_uncor_rand,
                    family = binomial(link = "logit"))
gam_sdm_rand <- mgcv::gam(formula_gam_rand,
                          data = sdm_data_uncor_rand,
                          family = binomial(link = "logit"))


summary(glm_sdm_rand)
summary(gam_sdm_rand)

par(mfrow = c(2, 3), mar = c(4, 4, 1, 1))
termplot(glm_sdm_rand,
         se            = TRUE,    # confidence ribbons
         rug           = TRUE,    # data rug at the bottom
         partial.resid = FALSE)
par(mfrow = c(1, 1))

plot(gam_sdm_rand,
     pages = 1,
     scale = 0)

# Project onto present, 2100, LGM
pred_glm_rand_present <- terra::predict(bioclim_data,         glm_sdm_rand, type = "response")
pred_gam_rand_present <- terra::predict(bioclim_data,         gam_sdm_rand, type = "response")
pred_glm_rand_2100    <- terra::predict(future_bioclim_data,  glm_sdm_rand, type = "response")
pred_gam_rand_2100    <- terra::predict(future_bioclim_data,  gam_sdm_rand, type = "response")
pred_glm_rand_lgm     <- terra::predict(lgm_bioclim_data,     glm_sdm_rand, type = "response")
pred_gam_rand_lgm     <- terra::predict(lgm_bioclim_data,     gam_sdm_rand, type = "response")

# Evaluate at present-day to get a threshold
p_rand <- na.omit(sdm_data_uncor_rand[sdm_data_uncor_rand$presence == 1, ])
a_rand <- na.omit(sdm_data_uncor_rand[sdm_data_uncor_rand$presence == 0, ])

glm_eval_rand <- evaluate(
  p = as.numeric(predict(glm_sdm_rand, newdata = p_rand, type = "response")),
  a = as.numeric(predict(glm_sdm_rand, newdata = a_rand, type = "response"))
)
gam_eval_rand <- evaluate(
  p = as.numeric(predict(gam_sdm_rand, newdata = p_rand, type = "response")),
  a = as.numeric(predict(gam_sdm_rand, newdata = a_rand, type = "response"))
)

thr_glm_rand_specsens <- threshold(glm_eval_rand, stat = "spec_sens")
thr_gam_rand_specsens <- threshold(gam_eval_rand, stat = "spec_sens")

# Binary distributions at the spec_sens threshold
pred_glm_rand_present_bin <- pred_glm_rand_present > thr_glm_rand_specsens
pred_gam_rand_present_bin <- pred_gam_rand_present > thr_gam_rand_specsens
pred_glm_rand_2100_bin    <- pred_glm_rand_2100    > thr_glm_rand_specsens
pred_gam_rand_2100_bin    <- pred_gam_rand_2100    > thr_gam_rand_specsens
pred_glm_rand_lgm_bin     <- pred_glm_rand_lgm     > thr_glm_rand_specsens
pred_gam_rand_lgm_bin     <- pred_gam_rand_lgm     > thr_gam_rand_specsens


# =============================================================================
# Task 3–5 (block 2 of 4): Accessible-area background pipeline
# =============================================================================

sdm_data_acc          <- bind_rows(occ_model, bg_model_acc)

drop_acc              <- drop_correlated(sdm_data_acc)
sdm_data_uncor_acc    <- drop_acc$data
predictor_vars_acc    <- drop_acc$predictor_vars
dropped_acc           <- drop_acc$dropped
print(predictor_vars_acc)
print(dropped_acc)

formula_glm_acc <- formula(paste0("presence ~ ",
                                  paste0(predictor_vars_acc, collapse = " + ")))
formula_gam_acc <- formula(paste0("presence ~ ",
                                  paste0("s(", predictor_vars_acc, ", k = 4)",
                                         collapse = " + ")))

glm_sdm_acc <- glm(formula_glm_acc,
                   data = sdm_data_uncor_acc,
                   family = binomial(link = "logit"))
gam_sdm_acc <- mgcv::gam(formula_gam_acc,
                         data = sdm_data_uncor_acc,
                         family = binomial(link = "logit"))

pred_glm_acc_present <- terra::predict(bioclim_data,        glm_sdm_acc, type = "response")
pred_gam_acc_present <- terra::predict(bioclim_data,        gam_sdm_acc, type = "response")
pred_glm_acc_2100    <- terra::predict(future_bioclim_data, glm_sdm_acc, type = "response")
pred_gam_acc_2100    <- terra::predict(future_bioclim_data, gam_sdm_acc, type = "response")
pred_glm_acc_lgm     <- terra::predict(lgm_bioclim_data,    glm_sdm_acc, type = "response")
pred_gam_acc_lgm     <- terra::predict(lgm_bioclim_data,    gam_sdm_acc, type = "response")

p_acc <- na.omit(sdm_data_uncor_acc[sdm_data_uncor_acc$presence == 1, ])
a_acc <- na.omit(sdm_data_uncor_acc[sdm_data_uncor_acc$presence == 0, ])

glm_eval_acc <- evaluate(
  p = as.numeric(predict(glm_sdm_acc, newdata = p_acc, type = "response")),
  a = as.numeric(predict(glm_sdm_acc, newdata = a_acc, type = "response"))
)
gam_eval_acc <- evaluate(
  p = as.numeric(predict(gam_sdm_acc, newdata = p_acc, type = "response")),
  a = as.numeric(predict(gam_sdm_acc, newdata = a_acc, type = "response"))
)

thr_glm_acc_specsens <- threshold(glm_eval_acc, stat = "spec_sens")
thr_gam_acc_specsens <- threshold(gam_eval_acc, stat = "spec_sens")

pred_glm_acc_present_bin <- pred_glm_acc_present > thr_glm_acc_specsens
pred_gam_acc_present_bin <- pred_gam_acc_present > thr_gam_acc_specsens
pred_glm_acc_2100_bin    <- pred_glm_acc_2100    > thr_glm_acc_specsens
pred_gam_acc_2100_bin    <- pred_gam_acc_2100    > thr_gam_acc_specsens
pred_glm_acc_lgm_bin     <- pred_glm_acc_lgm     > thr_glm_acc_specsens
pred_gam_acc_lgm_bin     <- pred_gam_acc_lgm     > thr_gam_acc_specsens


# =============================================================================
# Task 3–5 (block 3 of 4): Environmental envelope background pipeline
# =============================================================================

sdm_data_env          <- bind_rows(occ_model, bg_model_env)

drop_env              <- drop_correlated(sdm_data_env)
sdm_data_uncor_env    <- drop_env$data
predictor_vars_env    <- drop_env$predictor_vars
dropped_env           <- drop_env$dropped
print(predictor_vars_env)
print(dropped_env)

formula_glm_env <- formula(paste0("presence ~ ",
                                  paste0(predictor_vars_env, collapse = " + ")))
formula_gam_env <- formula(paste0("presence ~ ",
                                  paste0("s(", predictor_vars_env, ", k = 4)",
                                         collapse = " + ")))

glm_sdm_env <- glm(formula_glm_env,
                   data = sdm_data_uncor_env,
                   family = binomial(link = "logit"))
gam_sdm_env <- mgcv::gam(formula_gam_env,
                         data = sdm_data_uncor_env,
                         family = binomial(link = "logit"))

pred_glm_env_present <- terra::predict(bioclim_data,        glm_sdm_env, type = "response")
pred_gam_env_present <- terra::predict(bioclim_data,        gam_sdm_env, type = "response")
pred_glm_env_2100    <- terra::predict(future_bioclim_data, glm_sdm_env, type = "response")
pred_gam_env_2100    <- terra::predict(future_bioclim_data, gam_sdm_env, type = "response")
pred_glm_env_lgm     <- terra::predict(lgm_bioclim_data,    glm_sdm_env, type = "response")
pred_gam_env_lgm     <- terra::predict(lgm_bioclim_data,    gam_sdm_env, type = "response")

p_env <- na.omit(sdm_data_uncor_env[sdm_data_uncor_env$presence == 1, ])
a_env <- na.omit(sdm_data_uncor_env[sdm_data_uncor_env$presence == 0, ])

glm_eval_env <- evaluate(
  p = as.numeric(predict(glm_sdm_env, newdata = p_env, type = "response")),
  a = as.numeric(predict(glm_sdm_env, newdata = a_env, type = "response"))
)
gam_eval_env <- evaluate(
  p = as.numeric(predict(gam_sdm_env, newdata = p_env, type = "response")),
  a = as.numeric(predict(gam_sdm_env, newdata = a_env, type = "response"))
)

thr_glm_env_specsens <- threshold(glm_eval_env, stat = "spec_sens")
thr_gam_env_specsens <- threshold(gam_eval_env, stat = "spec_sens")

pred_glm_env_present_bin <- pred_glm_env_present > thr_glm_env_specsens
pred_gam_env_present_bin <- pred_gam_env_present > thr_gam_env_specsens
pred_glm_env_2100_bin    <- pred_glm_env_2100    > thr_glm_env_specsens
pred_gam_env_2100_bin    <- pred_gam_env_2100    > thr_gam_env_specsens
pred_glm_env_lgm_bin     <- pred_glm_env_lgm     > thr_glm_env_specsens
pred_gam_env_lgm_bin     <- pred_gam_env_lgm     > thr_gam_env_specsens


# =============================================================================
# Task 3–5 (block 4 of 4): Bias-corrected background pipeline
# =============================================================================

sdm_data_bias          <- bind_rows(occ_model, bg_model_bias)

drop_bias              <- drop_correlated(sdm_data_bias)
sdm_data_uncor_bias    <- drop_bias$data
predictor_vars_bias    <- drop_bias$predictor_vars
dropped_bias           <- drop_bias$dropped
print(predictor_vars_bias)
print(dropped_bias)

formula_glm_bias <- formula(paste0("presence ~ ",
                                   paste0(predictor_vars_bias, collapse = " + ")))
formula_gam_bias <- formula(paste0("presence ~ ",
                                   paste0("s(", predictor_vars_bias, ", k = 4)",
                                          collapse = " + ")))

glm_sdm_bias <- glm(formula_glm_bias,
                    data = sdm_data_uncor_bias,
                    family = binomial(link = "logit"))
gam_sdm_bias <- mgcv::gam(formula_gam_bias,
                          data = sdm_data_uncor_bias,
                          family = binomial(link = "logit"))

pred_glm_bias_present <- terra::predict(bioclim_data,        glm_sdm_bias, type = "response")
pred_gam_bias_present <- terra::predict(bioclim_data,        gam_sdm_bias, type = "response")
pred_glm_bias_2100    <- terra::predict(future_bioclim_data, glm_sdm_bias, type = "response")
pred_gam_bias_2100    <- terra::predict(future_bioclim_data, gam_sdm_bias, type = "response")
pred_glm_bias_lgm     <- terra::predict(lgm_bioclim_data,    glm_sdm_bias, type = "response")
pred_gam_bias_lgm     <- terra::predict(lgm_bioclim_data,    gam_sdm_bias, type = "response")

p_bias <- na.omit(sdm_data_uncor_bias[sdm_data_uncor_bias$presence == 1, ])
a_bias <- na.omit(sdm_data_uncor_bias[sdm_data_uncor_bias$presence == 0, ])

glm_eval_bias <- evaluate(
  p = as.numeric(predict(glm_sdm_bias, newdata = p_bias, type = "response")),
  a = as.numeric(predict(glm_sdm_bias, newdata = a_bias, type = "response"))
)
gam_eval_bias <- evaluate(
  p = as.numeric(predict(gam_sdm_bias, newdata = p_bias, type = "response")),
  a = as.numeric(predict(gam_sdm_bias, newdata = a_bias, type = "response"))
)

thr_glm_bias_specsens <- threshold(glm_eval_bias, stat = "spec_sens")
thr_gam_bias_specsens <- threshold(gam_eval_bias, stat = "spec_sens")

pred_glm_bias_present_bin <- pred_glm_bias_present > thr_glm_bias_specsens
pred_gam_bias_present_bin <- pred_gam_bias_present > thr_gam_bias_specsens
pred_glm_bias_2100_bin    <- pred_glm_bias_2100    > thr_glm_bias_specsens
pred_gam_bias_2100_bin    <- pred_gam_bias_2100    > thr_gam_bias_specsens
pred_glm_bias_lgm_bin     <- pred_glm_bias_lgm     > thr_glm_bias_specsens
pred_gam_bias_lgm_bin     <- pred_gam_bias_lgm     > thr_gam_bias_specsens



# Palette for the binary raster: 1 = suitable -> blue, 0 = unsuitable -> transparent
pal_bin <- colorNumeric(
  palette  = c("transparent", "#2C7FB8"),
  domain   = c(0, 1),
  na.color = "transparent"
)

leaflet() |>
  addProviderTiles(providers$CartoDB.Positron) |>
  
  # Thresholded suitability — GAM × accessible-area, present day
  addRasterImage(
    raster::raster(pred_gam_rand_present_bin),
    colors  = pal_bin,
    opacity = 0.7
  ) |>
  
  # Alpha hull from week 7 (red outline, no fill)
  addPolygons(
    data        = ahull_sf,
    fillColor   = "transparent",
    color       = "red",
    weight      = 2,
    opacity     = 0.9,
    popup       = "Alpha hull (week 7)"
  ) |>
  
  # Occurrence points (one per grid cell)
  addCircleMarkers(
    data        = occ_model,
    lng         = ~longitude,
    lat         = ~latitude,
    radius      = 4,
    fillOpacity = 0.8,
    stroke      = FALSE,
    color       = "black"
  ) |>
  
  # Legend
  addLegend(
    position = "bottomright",
    colors   = c("#2C7FB8", "red", "black"),
    labels   = c("Suitable (GAM × random)",
                 "Alpha hull (week 7)",
                 "Occurrence points"),
    title    = "Present-day suitability",
    opacity  = 0.8
  )



# -----------------------------------------------------------------------------
# Summary table — one row per (background method × model)
# -----------------------------------------------------------------------------

summary_tbl <- tribble(
  ~background, ~model, ~auc,                ~threshold_specsens,        ~n_predictors,
  "rand",      "GLM",  glm_eval_rand@auc,   thr_glm_rand_specsens,      length(predictor_vars_rand),
  "rand",      "GAM",  gam_eval_rand@auc,   thr_gam_rand_specsens,      length(predictor_vars_rand),
  "acc",       "GLM",  glm_eval_acc@auc,    thr_glm_acc_specsens,       length(predictor_vars_acc),
  "acc",       "GAM",  gam_eval_acc@auc,    thr_gam_acc_specsens,       length(predictor_vars_acc),
  "env",       "GLM",  glm_eval_env@auc,    thr_glm_env_specsens,       length(predictor_vars_env),
  "env",       "GAM",  gam_eval_env@auc,    thr_gam_env_specsens,       length(predictor_vars_env),
  "bias",      "GLM",  glm_eval_bias@auc,   thr_glm_bias_specsens,      length(predictor_vars_bias),
  "bias",      "GAM",  gam_eval_bias@auc,   thr_gam_bias_specsens,      length(predictor_vars_bias)
)
print(summary_tbl)


# =============================================================================
# Task 5 (continued). K-fold cross-validation
# =============================================================================
# Run on the random-bg pipeline (the qmd's "main path"). The structure is the
# same for any of the four bg methods — just swap sdm_data_uncor_<method>,
# glm_sdm_<method> and gam_sdm_<method>.
#
# Notes for answers:
#  - Training AUC is almost always over-optimistic. CV AUC drops by 0.01–0.10
#    typically, and the GAM tends to drop more than the GLM (it overfits more).
#  - High variability across folds = the model is sensitive to which records
#    you held out — usually a sign of small/biased presence data.
# -----------------------------------------------------------------------------

K <- 5
k_fold <- kfold(sdm_data_uncor_rand, k = K)

auc_df <- data.frame(k = 1:K, glm_AUC = NA, gam_AUC = NA)

for (k in 1:K) {

  train <- na.omit(sdm_data_uncor_rand[k_fold != k, ])
  test  <- na.omit(sdm_data_uncor_rand[k_fold == k, ])

  test_p <- test[test$presence == 1, ]
  test_a <- test[test$presence == 0, ]

  mod_glm <- update(glm_sdm_rand, data = train)
  mod_gam <- update(gam_sdm_rand, data = train)

  glm_eval_tmp <- evaluate(
    p = as.numeric(predict(mod_glm, newdata = test_p, type = "response")),
    a = as.numeric(predict(mod_glm, newdata = test_a, type = "response"))
  )
  gam_eval_tmp <- evaluate(
    p = as.numeric(predict(mod_gam, newdata = test_p, type = "response")),
    a = as.numeric(predict(mod_gam, newdata = test_a, type = "response"))
  )

  auc_df[k, "glm_AUC"] <- glm_eval_tmp@auc
  auc_df[k, "gam_AUC"] <- gam_eval_tmp@auc
}

print(auc_df)
cat("Mean CV AUC — GLM:", round(mean(auc_df$glm_AUC), 3),
    "  GAM:", round(mean(auc_df$gam_AUC), 3), "\n")


# =============================================================================
# Going Further (1) — Threshold sensitivity
# =============================================================================
# Same model (GAM, random bg) — different threshold metric. This isolates the
# effect of *threshold choice* on the predicted range.
#
# Notes for answers:
#  - no_omission is forced to include every presence -> tends to predict the
#    LARGEST range.
#  - sensitivity (=TPR) is similar to no_omission in spirit.
#  - spec_sens balances false positives and false negatives -> middle ground.
#  - kappa penalises false positives more in datasets with low prevalence ->
#    tends to predict the SMALLEST range.
# Range size can vary by ~2–3× across these on the same fitted model.
# -----------------------------------------------------------------------------

thr_rand_gam_specsens     <- threshold(gam_eval_rand, stat = "spec_sens")
thr_rand_gam_kappa        <- threshold(gam_eval_rand, stat = "kappa")
thr_rand_gam_no_omission  <- threshold(gam_eval_rand, stat = "no_omission")
thr_rand_gam_sensitivity  <- threshold(gam_eval_rand, stat = "sensitivity")

# Present
pred_gam_rand_present_bin_specsens     <- pred_gam_rand_present > thr_rand_gam_specsens
pred_gam_rand_present_bin_kappa        <- pred_gam_rand_present > thr_rand_gam_kappa
pred_gam_rand_present_bin_no_omission  <- pred_gam_rand_present > thr_rand_gam_no_omission
pred_gam_rand_present_bin_sensitivity  <- pred_gam_rand_present > thr_rand_gam_sensitivity

# 2100
pred_gam_rand_2100_bin_specsens        <- pred_gam_rand_2100 > thr_rand_gam_specsens
pred_gam_rand_2100_bin_kappa           <- pred_gam_rand_2100 > thr_rand_gam_kappa
pred_gam_rand_2100_bin_no_omission     <- pred_gam_rand_2100 > thr_rand_gam_no_omission
pred_gam_rand_2100_bin_sensitivity     <- pred_gam_rand_2100 > thr_rand_gam_sensitivity

# LGM
pred_gam_rand_lgm_bin_specsens         <- pred_gam_rand_lgm > thr_rand_gam_specsens
pred_gam_rand_lgm_bin_kappa            <- pred_gam_rand_lgm > thr_rand_gam_kappa
pred_gam_rand_lgm_bin_no_omission      <- pred_gam_rand_lgm > thr_rand_gam_no_omission
pred_gam_rand_lgm_bin_sensitivity      <- pred_gam_rand_lgm > thr_rand_gam_sensitivity


# =============================================================================
# Going Further (2 + 3) — Final comparison figures
# =============================================================================
# Three figures (present, 2100, LGM), each showing 4 bg methods × {GLM, GAM}
# at the spec_sens threshold. Plus a fourth figure varying the threshold.
#
# We use base terra::plot with par(mfrow) here because leafsync::sync() of
# 8+ leaflet panels is unreadable and slow. The wk7 solution drops into
# whichever plotting tool fits the comparison best — same idea here.
# -----------------------------------------------------------------------------

dir.create("figures", showWarnings = FALSE)

plot_panels_present <- function() {
  par(mfrow = c(2, 4), mar = c(2, 2, 2, 2))
  plot(pred_glm_rand_present_bin, main = "GLM × random",        legend = FALSE)
  plot(pred_glm_acc_present_bin,  main = "GLM × accessible",    legend = FALSE)
  plot(pred_glm_env_present_bin,  main = "GLM × env. envelope", legend = FALSE)
  plot(pred_glm_bias_present_bin, main = "GLM × bias-corrected",legend = FALSE)
  plot(pred_gam_rand_present_bin, main = "GAM × random",        legend = FALSE)
  plot(pred_gam_acc_present_bin,  main = "GAM × accessible",    legend = FALSE)
  plot(pred_gam_env_present_bin,  main = "GAM × env. envelope", legend = FALSE)
  plot(pred_gam_bias_present_bin, main = "GAM × bias-corrected",legend = FALSE)
}

plot_panels_2100 <- function() {
  par(mfrow = c(2, 4), mar = c(2, 2, 2, 2))
  plot(pred_glm_rand_2100_bin, main = "GLM × random — 2100",         legend = FALSE)
  plot(pred_glm_acc_2100_bin,  main = "GLM × accessible — 2100",     legend = FALSE)
  plot(pred_glm_env_2100_bin,  main = "GLM × env. envelope — 2100",  legend = FALSE)
  plot(pred_glm_bias_2100_bin, main = "GLM × bias-corrected — 2100", legend = FALSE)
  plot(pred_gam_rand_2100_bin, main = "GAM × random — 2100",         legend = FALSE)
  plot(pred_gam_acc_2100_bin,  main = "GAM × accessible — 2100",     legend = FALSE)
  plot(pred_gam_env_2100_bin,  main = "GAM × env. envelope — 2100",  legend = FALSE)
  plot(pred_gam_bias_2100_bin, main = "GAM × bias-corrected — 2100", legend = FALSE)
}

plot_panels_lgm <- function() {
  par(mfrow = c(2, 4), mar = c(2, 2, 2, 2))
  plot(pred_glm_rand_lgm_bin, main = "GLM × random — LGM",         legend = FALSE)
  plot(pred_glm_acc_lgm_bin,  main = "GLM × accessible — LGM",     legend = FALSE)
  plot(pred_glm_env_lgm_bin,  main = "GLM × env. envelope — LGM",  legend = FALSE)
  plot(pred_glm_bias_lgm_bin, main = "GLM × bias-corrected — LGM", legend = FALSE)
  plot(pred_gam_rand_lgm_bin, main = "GAM × random — LGM",         legend = FALSE)
  plot(pred_gam_acc_lgm_bin,  main = "GAM × accessible — LGM",     legend = FALSE)
  plot(pred_gam_env_lgm_bin,  main = "GAM × env. envelope — LGM",  legend = FALSE)
  plot(pred_gam_bias_lgm_bin, main = "GAM × bias-corrected — LGM", legend = FALSE)
}

plot_panels_threshold <- function() {
  par(mfrow = c(3, 4), mar = c(2, 2, 2, 2))
  # rows = time; cols = threshold metric
  plot(pred_gam_rand_lgm_bin_specsens,        main = "LGM — spec_sens",       legend = FALSE)
  plot(pred_gam_rand_lgm_bin_kappa,           main = "LGM — kappa",           legend = FALSE)
  plot(pred_gam_rand_lgm_bin_no_omission,     main = "LGM — no_omission",     legend = FALSE)
  plot(pred_gam_rand_lgm_bin_sensitivity,     main = "LGM — sensitivity",     legend = FALSE)
  plot(pred_gam_rand_present_bin_specsens,    main = "Present — spec_sens",   legend = FALSE)
  plot(pred_gam_rand_present_bin_kappa,       main = "Present — kappa",       legend = FALSE)
  plot(pred_gam_rand_present_bin_no_omission, main = "Present — no_omission", legend = FALSE)
  plot(pred_gam_rand_present_bin_sensitivity, main = "Present — sensitivity", legend = FALSE)
  plot(pred_gam_rand_2100_bin_specsens,       main = "2100 — spec_sens",      legend = FALSE)
  plot(pred_gam_rand_2100_bin_kappa,          main = "2100 — kappa",          legend = FALSE)
  plot(pred_gam_rand_2100_bin_no_omission,    main = "2100 — no_omission",    legend = FALSE)
  plot(pred_gam_rand_2100_bin_sensitivity,    main = "2100 — sensitivity",    legend = FALSE)
}

# Show on screen
plot_panels_present()
plot_panels_2100()
plot_panels_lgm()
plot_panels_threshold()

# Save to disk
png("figures/wk8_present_bg_x_model.png", width = 1600, height = 900, res = 120)
plot_panels_present(); dev.off()

png("figures/wk8_2100_bg_x_model.png",    width = 1600, height = 900, res = 120)
plot_panels_2100();    dev.off()

png("figures/wk8_lgm_bg_x_model.png",     width = 1600, height = 900, res = 120)
plot_panels_lgm();     dev.off()

png("figures/wk8_threshold_sensitivity.png", width = 1600, height = 1200, res = 120)
plot_panels_threshold(); dev.off()


# -----------------------------------------------------------------------------
# Final discussion
# -----------------------------------------------------------------------------
# Q: Which background selection step gave you the best estimate of the
#    present-day range?
# Q: How did the methods differ from one another?
# Q: Why might some methods perform better than others?
#
# Notes for answers:
#  - "Best" depends on what you trust as ground truth. Compared to the alpha
#    hull from week 7 (a presence-only summary), the bias-corrected GAM panel
#    is usually the closest match — it concentrates suitable habitat in the
#    Top End and Cape York and excludes most of arid central Australia.
#  - Random-background panels usually predict the broadest range (the model
#    has no penalty for "yes" in unsurveyed inland areas).
#  - Environmental-envelope panels can look surprisingly tight even though
#    they are structurally biased — backgrounds drawn from inside the niche
#    envelope are easier to separate from presences than they should be.
#  - 2100 panels generally show suitable habitat shifting / contracting in
#    northern Australia under MIROC6 SSP585.
#  - LGM panels can look very different and are less reliable — climates at
#    LGM lie outside the training range, and GAM smooths in particular can
#    swing wildly under extrapolation. Treat LGM panels as indicative.
# -----------------------------------------------------------------------------


# =============================================================================
# End of Week 8 solutions
# =============================================================================
