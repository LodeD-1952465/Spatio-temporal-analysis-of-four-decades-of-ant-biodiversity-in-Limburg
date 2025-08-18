###################### Spatio-temporal analysis of ants ########################
#                                                                              #
#Authors: Lode Deckers                                                         #
#Date: 16-08-2025                                                              #
#Course: Thesis                                                                #
################################################################################

################################## Packages ####################################

library(tidyverse)
library(ggplot2)
library(viridis)
library(tmap)
library(readxl)
library(INLA)
library(eurostat)
library(giscoR)
library(dplyr)
library(reshape2)
library(gstat)
library(sp)
library(sf)

############################### General setup ##################################
# Load directionary ------------------------------------------------------------
setwd(dir = "your_working_directory")
# reproducibility measures -----------------------------------------------------
set.seed(1952465) 

# General functions ------------------------------------------------------------
## Habitat classifier ----------------------------------------------------------
classify_habitat <- function(bwk_label) {
  codes <- strsplit(trimws(bwk_label), "\\+")  # Split by '+'
  
  habitats <- sapply(codes, function(code_set) {
    first_code <- trimws(code_set)[1] 
    first_letter <- tolower(substring(first_code, 1, 1))
    second_letter <- tolower(substring(first_code, 2, 2))  
    full_code <- tolower(first_code)
    habitat_classification <- case_when(
      first_letter == "c" & second_letter %in% 
        c("d", "v", "g") ~ "Droge heide",
      first_letter == "c" & second_letter == "e" ~ "Natte heide",
      first_letter == "c" ~ "Overige heide",
      first_letter %in% c("f", "l", "s", "n", "q", "r", "v") ~ "Loofbos",
      first_letter == "p" ~ "Naaldbos",
      first_letter == "b" ~ "Landbouw",
      first_letter == "h" ~ "Grasland",
      first_letter %in% c("u", "w") ~ "Woonomgeving",
      full_code == "ka" ~ "Overige",
      full_code == "kc" ~ "Landbouw", 
      full_code == "kd" ~ "Overige",
      full_code == "kf" ~ "Woonomgeving",
      full_code == "kf-" ~ "Woonomgeving",
      full_code == "kg" ~ "Terril",
      full_code == "ki" ~ "Overige",
      full_code == "kj" ~ "Hoogstamboomgaard",
      full_code == "kk" ~ "Overige",
      full_code == "kl" ~ "Landbouw",
      full_code == "ko" ~ "Woonomgeving",
      full_code == "kp" ~ "Woonomgeving",
      full_code == "kpa" ~ "Woonomgeving",
      full_code == "kpk" ~ "Woonomgeving",
      full_code == "kr" ~ "Overige",
      full_code == "ks" ~ "Woonomgeving",
      full_code == "kt" ~ "Overige",
      full_code == "kw" ~ "Overige",
      full_code == "kz" ~ "Overige",
      full_code == "wat" ~ "Overige",
      full_code == "ng" ~ "Overige",
      TRUE ~ "Overige"
    )
    return(habitat_classification)
  })
  
  return(habitats)
}

## Habitat translator (dutch to english) ---------------------------------------
group_and_translate_habitat <- function(dutch_habitat) {
  english_habitat <- case_when(
    dutch_habitat == "Droge heide" ~ "Dry heathland",
    dutch_habitat == "Natte heide" ~ "Wet heathland", 
    dutch_habitat %in% c("Loofbos", "Naaldbos") ~ "Forest",
    dutch_habitat == "Grasland" ~ "Grassland",
    dutch_habitat %in% c("Woonomgeving", "Landbouw") ~ "Human influence", 
    dutch_habitat %in% c("Terril", "Hoogstamboomgaard") ~ "Other",
    dutch_habitat %in% c("Overige", "Overige heide") ~ "Other",
    TRUE ~ "Other"
  )
  
  return(english_habitat)
}
## Existence of a model verifieer ----------------------------------------------
check_existing_model <- function(species, model_name, family) {
  file_path <- file.path(dir_path, paste0(species, "_", model_name, "_", 
                                          family, ".rds"))
  return(file.exists(file_path))
}


############################## General set up ##################################
# Load the data ----------------------------------------------------------------
mieren <- read.csv("mieren_cluster_only_pitfalls.csv")
colnames(mieren)[colnames(mieren) == "LATIJNSE.NAAM"] <- "LATIJNSE NAAM"
colnames(mieren)[colnames(mieren) == "Rode.Lijst.Vlaanderen"] <- 
  "Rode Lijst Vlaanderen"

# General parameters -----------------------------------------------------------
year_jumps <- 6 # amount of jumps between a year is made when more than 3
                # observation years are present
indangerd_species <- "" # Default is "". The value of this parameter selects the 
                        # most present species in that category of Rode Lijst
                        # Vlaanderen
width_inches <- 706 / 96  # width saved figures
height_inches <- 415 / 96  # hieiht saved figures
pixel_size <- 8  # Pixel size (units depend on CRS e.g., meters in UTM)
resolution_grid <- 500  # ...m grid cells
threshold_species <- 0.004 # Threshold for selection of ant species
threshold_hab <- 0.05 # Threshold for selection of habitat type 
n_perm <- 500 #number of permutations for the envelopes of the variograms
hab_of_interest <- c("Droge heide", "Grasland", "Landbouw", "Loofbos", 
                     "Naaldbos", "Natte heide", "Overige heide", 
                     "Woonomgeving", "Terril", "Hoogstamboomgaard")

# Information of inference model -----------------------------------------------
best_model_name <- "formula20_binomial"
best_model_formula <- "formula20"
best_model_family <- "binomial"

# Final prior set -------------------------------------------------------------- 
# (only applied when original prior set is not used in inferecence)
spde_range <- c(15000, 0.2)  
spde_sigma <- c(3, 0.2)   
ar1_param = c(0, 0.3)
prior_name = "vague"

# Set up: Visualization --------------------------------------------------------

## Transform the mieren coordinate system to the Limburg province coordinates
mieren_sf <- st_as_sf(mieren, coords = c("X", "Y"), crs = 31370)
mieren_sf$X <- st_coordinates(mieren_sf)[,1]
mieren_sf$Y <- st_coordinates(mieren_sf)[,2]

## Obtain the Limburg province 
search_eurostat("nuts")
nuts_sf <- get_eurostat_geospatial(output_class = "sf", resolution = "01", 
                                   nuts_level = 2)
belgium_nuts <- nuts_sf[nuts_sf$CNTR_CODE == "BE", ]
limburg <- belgium_nuts[belgium_nuts$NUTS_NAME == "Prov. Limburg (BE)", ] 
limburg_utm <- st_transform(limburg, crs = 31370)

## Vizualization using 
mieren_sf_longlat <- st_transform(mieren_sf, crs = 4326)
limburg_longlat <- st_transform(limburg, crs = 4326)

########################## Statistical analaysis ###############################
# Create folder to store results of statistical inference ----------------------
if (length(hab_of_interest) == 1){
  dir_path <- paste0("statistical inference results for ", hab_of_interest, 
                     " habitat")
} else {
  dir_path <- "statistical inference results for all habitats"
}
if (!dir.exists(dir_path)) {
  dir.create(dir_path)
  cat("Folder created successfully at:", dir_path, "\n")
} else {
  cat("Folder already exists at:", dir_path, "\n")
}

# Generate the group of species ------------------------------------------------
threshold_range <- seq(0.001, 0.01, by = 0.0005)
species_counts_by_threshold <- tibble(
  threshold = numeric(), 
  n_species = integer(),
  per_species = numeric(),
  species_names = list()
)

total_ant_species <- mieren %>%
  filter(Soortgroep == 14) %>%
  pull(`LATIJNSE NAAM`) %>%
  unique() %>%
  length()

for (curr_threshold in threshold_range) {
  threshold <- nrow(mieren) * curr_threshold
  common_species_df <- mieren %>%
    count(`LATIJNSE NAAM`, Soortgroep, name = "count") %>% 
    filter(Soortgroep == 14) %>% 
    filter(count >= threshold)
  n_species <- n_distinct(common_species_df$`LATIJNSE NAAM`)
  species_names_vec <- common_species_df$`LATIJNSE NAAM`
  species_counts_by_threshold <- add_row(species_counts_by_threshold,  
                                         threshold = curr_threshold,  
                                         n_species = n_species,  
                                         per_species = n_species / 
                                           total_ant_species, 
                                         species_names = list(
                                           species_names_vec))
}

saveRDS(species_counts_by_threshold, 
        file = file.path(dir_path,
                         "species_group_per_tresholds.rds"))
species_group <- species_counts_by_threshold %>% 
  filter(threshold == threshold_species) %>% pull(species_names) %>% .[[1]]

if (length(hab_of_interest) == 1){
  if (hab_of_interest == hab1){
    species_group <- species_list_droge
  }
  
  if (hab_of_interest == hab2){
    species_group <- species_list_natte
  }
}

# Filter on habitat types that have sufficient observations --------------------
mieren_specific <- mieren %>% filter(`LATIJNSE NAAM` %in% species_group)
tot_obs <- nrow(mieren_specific)
presence_per_hab_specific <- mieren_specific %>%
  group_by(habitat_type) %>%summarise(n_present = sum(is_present_bin == 1, 
                                                      na.rm = TRUE),
                                      presence_rate = (n_present / tot_obs) * 
                                        100)
write.csv(presence_per_hab_specific, 
          file.path(dir_path, "percence_rate_habitat_species_analysed.csv"),
          row.names = FALSE)

hab_analysis <- presence_per_hab_specific %>%
  filter(presence_rate > 5 | habitat_type == 'Wet heathland' |
           habitat_type == "Forest" | habitat_type == "Human influence" |
           habitat_type == "Grassland" | habitat_type == "Dry heathland") %>% 
  pull(habitat_type)


mieren_data <- mieren %>%
  mutate(habitat_type_analysis = case_when(
    habitat_type %in% hab_analysis ~ habitat_type,
    TRUE ~ 'Other'))
mieren_data$habitat_type <- mieren_data$habitat_type_analysis

# Ensure proper habitat grouping -----------------------------------------------
pred_grid_path <- file.path(dir_path, "pred_grid_list.rds")
if (file.exists(pred_grid_path)) {
  file.remove(pred_grid_path)
}
# Fitting a single species occupancy model for every species -------------------
results_list <- list()
for (curr_species in species_group) {
  # Data preparation -----------------------------------------------------------
  mieren_specific <- mieren %>% filter(`LATIJNSE NAAM` %in% species_group)
  tot_obs <- nrow(mieren_specific)
  presence_per_hab_specific <- mieren_specific %>%
    group_by(habitat_type) %>%summarise(n_present = sum(is_present_bin == 1, 
                                                        na.rm = TRUE),
                                        presence_rate = (n_present / tot_obs) * 
                                          100)
  # Filter for habitat types with sufficient presence OR important habitat types
  hab_analysis <- presence_per_hab_specific %>%
    filter(presence_rate > 5 | habitat_type == 'Wet heathland' 
           | habitat_type == 'Human influence' | habitat_type == 'Grassland' 
           | habitat_type == 'Forest' | habitat_type == 'Dry heathland') %>%
    pull(habitat_type) 
  
  mieren_data <- mieren %>%
    mutate(habitat_type_analysis = case_when(
      habitat_type %in% hab_analysis ~ habitat_type,
      TRUE ~ 'Other'))
  mieren_data$habitat_type <- mieren_data$habitat_type_analysis
  
  is_present <- numeric(nrow(mieren_data))
  for (i in 1:nrow(mieren_data)){
    if (mieren_data$`LATIJNSE NAAM`[i] == curr_species){
      if (mieren_data$Aantal[i] > 0){
        is_present[i] <- 1
      }
    }
  }
  
  mieren_data$is_present <- is_present
  mieren_data$jaar1_num <- as.numeric(mieren_data$jaar1)
  mieren_data$year_group <- as.integer(factor(mieren_data$jaar1))
  mieren_data$presence <- factor(ifelse(mieren_data$is_present == 1, 
                                        "present", "absent"), 
                                 levels = c("absent", "present"))
  mieren_data <- mieren_data %>% filter(!is.na(habitat_type))
  mieren_data$habitat_type <- as.factor(mieren_data$habitat_type)
  mieren_data$hab_type_numeric <- as.numeric(as.factor(mieren_data$
                                                         habitat_type))
  mieren_data$is_winter <- factor(ifelse(mieren_data$is_winter_bin == 1, 
                                         "winter", "no winter"), 
                                  levels = c("no winter", "winter"))
  mieren_data$jaar1_hab <- interaction(mieren_data$year_group, 
                                       mieren_data$habitat_type, drop = TRUE)
  mieren_data$jaar1_num <- as.numeric(mieren_data$jaar1)
  mieren_data$year_group_global <- mieren_data$year_group
  mieren_data$year_group_habitat <- mieren_data$year_group
  
  # Data vizualization 
  if (!exists("limburg_utm_cached")) {
    search_eurostat("nuts")
    nuts_sf <- get_eurostat_geospatial(output_class = "sf", resolution = "01", 
                                       nuts_level = 2)
    belgium_nuts <- nuts_sf[nuts_sf$CNTR_CODE == "BE", ]
    limburg <- belgium_nuts[belgium_nuts$NUTS_NAME == "Prov. Limburg (BE)", ] 
    limburg_utm_cached <<- st_transform(limburg, crs = 31370)
    rm(nuts_sf, belgium_nuts, limburg)
    gc()
  }
  
  limburg_utm <- limburg_utm_cached
  print(paste0("data preparation of species ", curr_species, " finished."))
  
  # Statistical analysis: set up -----------------------------------------------
  ## Mesh construction
  coo <- cbind(mieren_data$X, mieren_data$Y)
  mesh_file_path <- file.path(dir_path, "mesh.rds")
  if (file.exists(mesh_file_path)) {
    mesh_list <- readRDS(mesh_file_path)
    bnd <- mesh_list$bnd
    mesh <- mesh_list$mesh
    spde <- mesh_list$spde
  } else {
    bnd <- inla.nonconvex.hull(st_coordinates(limburg_utm)[,1:2])
    mesh <- inla.mesh.2d(loc = coo, boundary = bnd, max.edge = c(20000, 40000),  
                         cutoff = 2000)
    spde <- inla.spde2.pcmatern(mesh = mesh, alpha = 2, constr = TRUE,
                                prior.range = c(5000, 0.5), 
                                prior.sigma = c(1, 0.01))
    
    saveRDS(list(bnd = bnd, mesh = mesh, spde = spde), mesh_file_path)
  }
  ## Index set construction
  timesn <- length(unique(mieren_data$jaar1))
  indexs <- inla.spde.make.index("s", n.spde = spde$n.spde, n.group = timesn)
  lengths(indexs)
  
  ## Projection matrix construction 
  A <- inla.spde.make.A(
    mesh = mesh,
    loc = cbind(mieren_data$X, mieren_data$Y),
    group = as.integer(factor(mieren_data$jaar1))
  )
  
  ## Prediction location construction
  pred_grid_path <- file.path(dir_path, "pred_grid_list.rds")
  if (file.exists(pred_grid_path)) {
    pred_grid <- readRDS(pred_grid_path)
    grid_habitats_clean <- pred_grid$grid_habitats_clean
    n_cells_clean <- pred_grid$n_cells_clean
  } else {
    bb <- st_bbox(limburg_utm)
    grid <- st_make_grid(
      limburg_utm,
      cellsize = resolution_grid,
      what = "polygons",
      square = TRUE) |> st_as_sf() |> st_intersection(limburg_utm)
    
    grid_predictions <- grid %>% mutate(cell_id = 1:n()) %>% 
      st_set_crs(st_crs(limburg_utm))
    
    # Process dominant habitats
    grid_predictions <- st_transform(grid_predictions, st_crs(limburg_utm))
    BWK_map <- st_transform(BWK_map, st_crs(limburg_utm))
    
    dominant_habitats <- grid_predictions %>% st_intersection(BWK_map) %>%
      mutate(habitat_type_dutch = sapply(BWKLABEL, function(x) classify_habitat(x)[1]),
             habitat_type = group_and_translate_habitat(habitat_type_dutch),
             fragment_area = st_area(.)) %>%  st_drop_geometry() %>%
      group_by(cell_id, habitat_type) %>% summarise(area = sum(fragment_area), 
                                                    .groups = "drop_last") %>%
      mutate(percentage = as.numeric(area / sum(area) * 100)) %>% ungroup() %>%
      group_by(cell_id) %>% slice_max(percentage, n = 1, with_ties = FALSE) %>%
      dplyr::select(cell_id, dominant_habitat = habitat_type, 
                    dominant_percentage = percentage)
    
    grid_habitats <- grid_predictions %>%
      left_join(dominant_habitats, by = "cell_id")
    na_cells <- is.na(grid_habitats$dominant_habitat)
    grid_habitats_clean <- grid_habitats[!na_cells, ]
    
    ## Filter for habitat types that are observed
    unique_obs_habitats <- levels(mieren_data$habitat_type)
    valid_habitat <- grid_habitats_clean$dominant_habitat %in% 
      unique_obs_habitats
    grid_habitats_clean <- grid_habitats_clean[valid_habitat, ]
    n_cells_clean <- nrow(grid_habitats_clean)
    grid_habitats_clean$
      dominant_habitat <- factor(grid_habitats_clean$dominant_habitat, 
                                 levels = levels(mieren_data$habitat_type))
    saveRDS(list(grid_habitats_clean = grid_habitats_clean,
                 n_cells_clean = n_cells_clean), pred_grid_path)
    rm(grid, grid_predictions, dominant_habitats, grid_habitats)
    gc()
  }
  
  ## Prediction habitat specification
  observed_years <- sort(unique(mieren_data$jaar1))
  observed_year_groups <- sort(unique(mieren_data$year_group))
  n_years <- length(observed_years)
  
  expanded_grid_clean <- do.call(rbind, lapply(1:n_years, function(i) {
    cbind(grid_habitats_clean, 
          year = observed_years[i],
          year_group = observed_year_groups[i])
  }))
  
  
  ## Spatial index construction 
  spatial.index <- inla.spde.make.index(name = "s", 
                                        n.spde = spde$n.spde,       
                                        n.group = length(unique(mieren_data$
                                                                  jaar1)))
  
  spatial.index.only <- inla.spde.make.index(name = "spatial_field", 
                                             n.spde = spde$n.spde)
  A_spatial <- inla.spde.make.A(mesh, loc = cbind(mieren_data$X, mieren_data$Y))
  
  
  ### Prepare prediction variables
  pred_hab_type <- factor(rep(grid_habitats_clean$dominant_habitat,
                              n_years), 
                          levels = levels(mieren_data$habitat_type))
  pred_hab_type_numeric <- as.numeric(as.factor(pred_hab_type))
  pred_year_group <- rep(observed_year_groups, each = n_cells_clean)
  pred_jaar1 <- rep(observed_years, each = n_cells_clean)
  pred_jaar1_hab <- interaction(pred_year_group, pred_hab_type, drop = TRUE)
  pred_is_winter <- factor(rep("no winter", n_cells_clean * n_years), 
                           levels = levels(mieren_data$is_winter))
  mieren_data$year_group_global <- mieren_data$year_group
  mieren_data$year_group_habitat <- mieren_data$year_group
  pred_year_group_global <- pred_year_group
  pred_year_group_habitat <- pred_year_group
  ### Create temporal indices, data frame and projection matrix
  indexs_pred <- inla.spde.make.index("s", n.spde = mesh$n, n.group = n_years)
  prediction_points_clean <- st_centroid(expanded_grid_clean)
  saveRDS(prediction_points_clean, 
          file = file.path(dir_path, paste0(curr_species, 
                                            "_prediction_points.rds")))
  year_indices <- match(expanded_grid_clean$year, observed_years)
  Ap <- inla.spde.make.A(mesh = mesh,
                         loc = st_coordinates(prediction_points_clean),
                         group = year_indices, fast = TRUE)
  
  indexs_pred_spatial <- inla.spde.make.index(name = "spatial_field", 
                                              n.spde = mesh$n)
  Ap_spatial <- inla.spde.make.A(mesh = mesh, 
                                 loc = st_coordinates(prediction_points_clean),
                                 fast = TRUE)
  print(paste0("Statistical set up of species ", curr_species, " finished."))
  
  # statistical analysis: model specification ----------------------------------
  rprior <- list(theta = list(prior = "pccor1", param = c(0, 0.9)))
  
  formula1 <- y ~ is_winter + hab_type + f(spatial_field, model = spde) 
  
  formula2 <- y ~ is_winter + hab_type + f(jaar1, model = "ar1", hyper = rprior)
  
  formula3 <- y ~ is_winter + hab_type + 
    f(jaar1, model = "ar1", hyper = rprior) + f(spatial_field, model = spde) 
  
  formula4 <- y ~ is_winter + hab_type + 
    f(jaar1, model = "ar1", hyper = rprior) + f(spatial_field, model = spde) + 
    f(s, model = "iid")
  
  formula5 <- y ~ is_winter + hab_type + 
    f(jaar1, model = "ar1", hyper = rprior) + f(spatial_field, model = spde) + 
    f(s, model = "iid", group = s.group, control.group = list(model = "ar1", 
                                                              hyper = rprior))
  
  formula6 <- y ~ is_winter + hab_type + 
    f(jaar1, model = "ar1", hyper = rprior) + f(spatial_field, model = spde) + 
    f(s, model = spde, group = s.group, control.group = list(model = "iid"))
  
  formula7 <- y ~ is_winter + hab_type + 
    f(jaar1, model = "ar1", hyper = rprior) + f(spatial_field, model = spde) + 
    f(s, model = spde, group = s.group, control.group = list(model = "ar1", 
                                                             hyper = rprior))
  
  # Same model building approach but now with habitat and year interaction 
  formula8 <- y ~ is_winter + hab_type + jaar1:hab_type 
  
  formula9 <- y ~ is_winter + hab_type + jaar1:hab_type + 
    f(spatial_field, model = spde)
  
  formula10 <- y ~ is_winter + hab_type + jaar1:hab_type + 
    f(jaar1, model = "ar1", hyper = rprior)
  
  formula11 <- y ~ is_winter + hab_type +  jaar1:hab_type +
    f(spatial_field, model = spde) + f(jaar1, model = "ar1", hyper = rprior)
  
  formula12 <- y ~ is_winter + hab_type + jaar1:hab_type +
    f(jaar1, model = "ar1", hyper = rprior) + f(spatial_field, model = spde) + 
    f(s, model = "iid")
  
  formula13 <- y ~ is_winter + hab_type + jaar1:hab_type +
    f(jaar1, model = "ar1", hyper = rprior) + f(spatial_field, model = spde) + 
    f(s, model = "iid", group = s.group, control.group = list(model = "ar1", 
                                                              hyper = rprior))
  
  formula14 <- y ~ is_winter + hab_type + jaar1:hab_type +
    f(jaar1, model = "ar1", hyper = rprior) + f(spatial_field, model = spde) + 
    f(s, model = spde, group = s.group, control.group = list(model = "iid"))
  
  formula15 <- y ~ is_winter + hab_type + jaar1:hab_type + 
    f(jaar1, model = "ar1", hyper = rprior) + f(spatial_field, model = spde) + 
    f(s, model = spde, group = s.group, control.group = list(model = "ar1", 
                                                             hyper = rprior))
  
  # Same model building approach but now with habitat specific temporal AR1 
  # random effect 
  formula16 <- y ~ is_winter + hab_type + 
    f(inla.group(year_group, n = length(unique(mieren_data$year_group))), 
      model = "ar1", replicate = hab_type_numeric, hyper = rprior)
  
  formula17 <- y ~ is_winter + hab_type + 
    f(inla.group(year_group, n = length(unique(mieren_data$year_group))), 
      model = "ar1", replicate = hab_type_numeric, hyper = rprior) + 
    f(spatial_field, model = spde)
  
  formula18 <- y ~ is_winter + hab_type + 
    f(inla.group(year_group, n = length(unique(mieren_data$year_group))), 
      model = "ar1", replicate = hab_type_numeric, hyper = rprior) + 
    f(spatial_field, model = spde) + 
    f(s, model = "iid")
  
  formula19 <- y ~ is_winter + hab_type + 
    f(inla.group(year_group, n = length(unique(mieren_data$year_group))), 
      model = "ar1", replicate = hab_type_numeric, hyper = rprior) + 
    f(spatial_field, model = spde) + 
    f(s, model = "iid", group = s.group, 
      control.group = list(model = "ar1", hyper = rprior))
  
  formula20 <- y ~ is_winter + hab_type + 
    f(inla.group(year_group, n = length(unique(mieren_data$year_group))), 
      model = "ar1", replicate = hab_type_numeric, hyper = rprior) + 
    f(spatial_field, model = spde) + 
    f(s, model = spde, group = s.group, 
      control.group = list(model = "iid"))
  
  formula21 <- y ~ is_winter + hab_type + 
    f(inla.group(year_group, n = length(unique(mieren_data$year_group))), 
      model = "ar1", replicate = hab_type_numeric, hyper = rprior) + 
    f(spatial_field, model = spde) + 
    f(s, model = spde, group = s.group, 
      control.group = list(model = "ar1", hyper = rprior))
  
  # Same model building approach but now with habitat specific temporal AR1 
  # random effect and interaction effect between year and habitat type 
  formula22 <- y ~ is_winter + hab_type + jaar1:hab_type + 
    f(inla.group(year_group, n = length(unique(mieren_data$year_group))), 
      model = "ar1", replicate = hab_type_numeric, hyper = rprior)
  
  formula23 <- y ~ is_winter + hab_type + jaar1:hab_type + 
    f(inla.group(year_group, n = length(unique(mieren_data$year_group))), 
      model = "ar1", replicate = hab_type_numeric, hyper = rprior) + 
    f(spatial_field, model = spde)
  
  formula24 <- y ~ is_winter + hab_type + jaar1:hab_type + 
    f(inla.group(year_group, n = length(unique(mieren_data$year_group))), 
      model = "ar1", replicate = hab_type_numeric, hyper = rprior) + 
    f(spatial_field, model = spde) + 
    f(s, model = "iid")
  
  formula25 <- y ~ is_winter + hab_type + jaar1:hab_type + 
    f(inla.group(year_group, n = length(unique(mieren_data$year_group))), 
      model = "ar1", replicate = hab_type_numeric, hyper = rprior) + 
    f(spatial_field, model = spde) + 
    f(s, model = "iid", group = s.group, control.group = list(model = "ar1", 
                                                              hyper = rprior))
  
  formula26 <- y ~ is_winter + hab_type + jaar1:hab_type + 
    f(inla.group(year_group, n = length(unique(mieren_data$year_group))), 
      model = "ar1", replicate = hab_type_numeric, hyper = rprior) + 
    f(spatial_field, model = spde) + 
    f(s, model = spde, group = s.group, control.group = list(model = "iid"))
  
  formula27 <- y ~ is_winter + hab_type + jaar1:hab_type + 
    f(inla.group(year_group, n = length(unique(mieren_data$year_group))), 
      model = "ar1", replicate = hab_type_numeric, hyper = rprior) + 
    f(spatial_field, model = spde) + 
    f(s, model = spde, group = s.group, 
      control.group = list(model = "ar1", hyper = rprior))
  
  # Same model building approach but now with habitat specific temporal AR1 
  # random effect
  formula28 <- y ~ is_winter + hab_type + 
    f(inla.group(year_group, n = length(unique(mieren_data$year_group))), 
      model = "ar1", group = hab_type_numeric, hyper = rprior)
  
  formula29 <- y ~ is_winter + hab_type + 
    f(inla.group(year_group, n = length(unique(mieren_data$year_group))), 
      model = "ar1", group = hab_type_numeric, hyper = rprior) + 
    f(spatial_field, model = spde)
  
  formula30 <- y ~ is_winter + hab_type + 
    f(inla.group(year_group, n = length(unique(mieren_data$year_group))), 
      model = "ar1", group = hab_type_numeric, hyper = rprior) + 
    f(spatial_field, model = spde) + 
    f(s, model = "iid")
  
  formula31 <- y ~ is_winter + hab_type + 
    f(inla.group(year_group, n = length(unique(mieren_data$year_group))), 
      model = "ar1", group = hab_type_numeric, hyper = rprior) + 
    f(spatial_field, model = spde) + 
    f(s, model = "iid", group = s.group, 
      control.group = list(model = "ar1", hyper = rprior))
  
  formula32 <- y ~ is_winter + hab_type + 
    f(inla.group(year_group, n = length(unique(mieren_data$year_group))), 
      model = "ar1", group = hab_type_numeric, hyper = rprior) + 
    f(spatial_field, model = spde) + 
    f(s, model = spde, group = s.group, 
      control.group = list(model = "iid"))
  
  formula33 <- y ~ is_winter + hab_type + 
    f(inla.group(year_group, n = length(unique(mieren_data$year_group))), 
      model = "ar1", group = hab_type_numeric, hyper = rprior) + 
    f(spatial_field, model = spde) + 
    f(s, model = spde, group = s.group, 
      control.group = list(model = "ar1", hyper = rprior))
  
  # Stack for all different models 
  ## Stack for formula1 (spatial random effect)
  stk.est_1 <- inla.stack(
    tag = "est",
    data = list(y = mieren_data$is_present),
    A = list(1, A_spatial),
    effects = list(
      data.frame(
        hab_type = mieren_data$habitat_type,
        is_winter = mieren_data$is_winter),
      list(spatial_field = spatial.index.only$spatial_field)))
  
  stk.pred_1 <- inla.stack(
    tag = "pred",
    data = list(y = NA),
    A = list(1, Ap_spatial),
    effects = list(
      data.frame(
        hab_type = pred_hab_type,
        is_winter = pred_is_winter),
      list(spatial_field = indexs_pred_spatial$spatial_field)))
  
  stk.full_1 <- inla.stack(stk.est_1, stk.pred_1)
  
  ## Stack for formula2 (temporal AR1 random effect)
  stk.est_2 <- inla.stack(
    tag = "est", 
    data = list(y = mieren_data$is_present),
    A = list(1),
    effects = list(
      data.frame(
        hab_type = mieren_data$habitat_type,
        is_winter = mieren_data$is_winter,
        jaar1 = mieren_data$year_group)))
  
  stk.pred_2 <- inla.stack(
    tag = "pred",
    data = list(y = NA), 
    A = list(1),
    effects = list(
      data.frame(
        hab_type = pred_hab_type,
        is_winter = pred_is_winter,
        jaar1 = pred_year_group)))
  
  stk.full_2 <- inla.stack(stk.est_2, stk.pred_2)
  
  ## Stack for formula3 (temporal AR1 + spatial random effect)
  stk.est_3 <- inla.stack(
    tag = "est",
    data = list(y = mieren_data$is_present),
    A = list(1, A_spatial),
    effects = list(
      data.frame(
        hab_type = mieren_data$habitat_type,
        is_winter = mieren_data$is_winter,
        jaar1 = mieren_data$year_group),
      list(spatial_field = spatial.index.only$spatial_field)))
  
  stk.pred_3 <- inla.stack(
    tag = "pred",
    data = list(y = NA),
    A = list(1, Ap_spatial),
    effects = list(
      data.frame(
        hab_type = pred_hab_type,
        is_winter = pred_is_winter,
        jaar1 = pred_year_group),
      list(spatial_field = indexs_pred_spatial$spatial_field)))
  
  stk.full_3 <- inla.stack(stk.est_3, stk.pred_3)
  
  # Stack for formula4 (spatio-temporal random effect type 1)
  stk.est_4 <- inla.stack(
    tag = "est",
    data = list(y = mieren_data$is_present),
    A = list(1, A_spatial, A),
    effects = list(
      data.frame(
        hab_type = mieren_data$habitat_type,
        is_winter = mieren_data$is_winter,
        jaar1 = mieren_data$year_group),
      list(spatial_field = spatial.index.only$spatial_field),
      list(s = indexs$s, s.group = indexs$s.group)))
  
  stk.pred_4 <- inla.stack(
    tag = "pred",
    data = list(y = NA),
    A = list(1, Ap_spatial, Ap),
    effects = list(
      data.frame(
        hab_type = pred_hab_type,
        is_winter = pred_is_winter,
        jaar1 = pred_year_group),
      list(spatial_field = indexs_pred_spatial$spatial_field),
      list(s = indexs_pred$s, s.group = indexs_pred$s.group)))
  
  stk.full_4 <- inla.stack(stk.est_4, stk.pred_4)
  
  # Stack for formula5 (spatio-temporal random effect type 2)
  stk.full_5 <- stk.full_4
  
  # Stack for formula6 (spatio-temporal random effect type 3)
  stk.full_6 <- stk.full_4
  
  # Stack for formula7 (spatio-temporal random effect type 4)
  stk.full_7 <- stk.full_4
  
  # Stack for formula8 (interaction habitat and year)
  stk.est_8 <- inla.stack(
    tag = "est", 
    data = list(y = mieren_data$is_present),
    A = list(1),
    effects = list(
      data.frame(
        hab_type = mieren_data$habitat_type,
        is_winter = mieren_data$is_winter,
        jaar1 = mieren_data$year_group)))
  
  stk.pred_8 <- inla.stack(
    tag = "pred",
    data = list(y = NA), 
    A = list(1),
    effects = list(
      data.frame(
        hab_type = pred_hab_type,
        is_winter = pred_is_winter,
        jaar1 = pred_year_group)))
  
  stk.full_8 <- inla.stack(stk.est_8, stk.pred_8)
  
  # Stack for formula9 (interaction + spatial random effect)
  stk.est_9 <- inla.stack(
    tag = "est", 
    data = list(y = mieren_data$is_present),
    A = list(1, A_spatial),
    effects = list(
      data.frame(
        hab_type = mieren_data$habitat_type,
        is_winter = mieren_data$is_winter,
        jaar1 = mieren_data$year_group),
      list(spatial_field = spatial.index.only$spatial_field)))
  
  stk.pred_9 <- inla.stack(
    tag = "pred",
    data = list(y = NA), 
    A = list(1, Ap_spatial),
    effects = list(
      data.frame(
        hab_type = pred_hab_type,
        is_winter = pred_is_winter,
        jaar1 = pred_year_group),
      list(spatial_field = indexs_pred_spatial$spatial_field)))
  
  stk.full_9 <- inla.stack(stk.est_9, stk.pred_9)
  
  # Stack for formula10 (interaction + temporal AR1)
  stk.est_10 <- inla.stack(
    tag = "est",
    data = list(y = mieren_data$is_present),
    A = list(1),
    effects = list(
      data.frame(
        hab_type = mieren_data$habitat_type,
        is_winter = mieren_data$is_winter,
        jaar1 = mieren_data$year_group)))
  
  stk.pred_10 <- inla.stack(
    tag = "pred",
    data = list(y = NA),
    A = list(1),
    effects = list(
      data.frame(
        hab_type = pred_hab_type,
        is_winter = pred_is_winter,
        jaar1 = pred_year_group)))
  
  stk.full_10 <- inla.stack(stk.est_10, stk.pred_10)
  
  
  # Stack for formula11 (interaction + temporal AR1 + spatial random effect)
  stk.est_11 <- inla.stack(
    tag = "est", 
    data = list(y = mieren_data$is_present),
    A = list(1, A_spatial),
    effects = list(
      data.frame(
        hab_type = mieren_data$habitat_type,
        is_winter = mieren_data$is_winter,
        jaar1 = mieren_data$year_group),
      list(spatial_field = spatial.index.only$spatial_field)))
  
  stk.pred_11 <- inla.stack(
    tag = "pred",
    data = list(y = NA), 
    A = list(1, Ap_spatial),
    effects = list(
      data.frame(
        hab_type = pred_hab_type,
        is_winter = pred_is_winter,
        jaar1 = pred_year_group),
      list(spatial_field = indexs_pred_spatial$spatial_field)))
  
  stk.full_11 <- inla.stack(stk.est_11, stk.pred_11)
  
  # Stack for formula12 (interaction + spatio-temporal type 1)
  stk.est_12 <- inla.stack(
    tag = "est", 
    data = list(y = mieren_data$is_present),
    A = list(1, A_spatial, A),
    effects = list(
      data.frame(
        hab_type = mieren_data$habitat_type,
        is_winter = mieren_data$is_winter,
        jaar1 = mieren_data$year_group),
      list(spatial_field = spatial.index.only$spatial_field),
      list(s = indexs$s, s.group = indexs$s.group)))
  
  stk.pred_12 <- inla.stack(
    tag = "pred",
    data = list(y = NA), 
    A = list(1, Ap_spatial, Ap),
    effects = list(
      data.frame(
        hab_type = pred_hab_type,
        is_winter = pred_is_winter,
        jaar1 = pred_year_group),
      list(spatial_field = indexs_pred_spatial$spatial_field),
      list(s = indexs_pred$s, s.group = indexs_pred$s.group)))
  
  stk.full_12 <- inla.stack(stk.est_12, stk.pred_12)
  
  # Stack for formula13 (interaction + spatio-temporal type 2)
  stk.full_13 <- stk.full_12
  
  # Stack for formula14 (interaction + spatio-temporal type 3)
  stk.full_14 <- stk.full_12
  
  # Stack for formula15 (interaction + spatio-temporal type 4)
  stk.full_15 <- stk.full_12
  
  stk.est_16 <- inla.stack(
    tag = "est",
    data = list(y = mieren_data$is_present),
    A = list(1),
    effects = list(
      data.frame(
        hab_type = mieren_data$habitat_type,
        is_winter = mieren_data$is_winter,
        year_group = mieren_data$year_group,
        hab_type_numeric = mieren_data$hab_type_numeric)))
  
  stk.pred_16 <- inla.stack(
    tag = "pred",
    data = list(y = NA),
    A = list(1),
    effects = list(
      data.frame(
        hab_type = pred_hab_type,
        is_winter = pred_is_winter,
        year_group = pred_year_group,
        hab_type_numeric = pred_hab_type_numeric)))
  
  stk.full_16 <- inla.stack(stk.est_16, stk.pred_16)
  
  # Stack for formula17 (habitat-specific temporal AR1 + spatial)
  stk.est_17 <- inla.stack(
    tag = "est",
    data = list(y = mieren_data$is_present),
    A = list(1, A_spatial),
    effects = list(
      data.frame(
        hab_type = mieren_data$habitat_type,
        is_winter = mieren_data$is_winter,
        year_group = mieren_data$year_group,
        hab_type_numeric = mieren_data$hab_type_numeric),
      list(spatial_field = spatial.index.only$spatial_field)))
  
  stk.pred_17 <- inla.stack(
    tag = "pred",
    data = list(y = NA),
    A = list(1, Ap_spatial),
    effects = list(
      data.frame(
        hab_type = pred_hab_type,
        is_winter = pred_is_winter,
        year_group = pred_year_group,
        hab_type_numeric = pred_hab_type_numeric),
      list(spatial_field = indexs_pred_spatial$spatial_field)))
  
  stk.full_17 <- inla.stack(stk.est_17, stk.pred_17)
  
  # Stack for formula18 (habitat-specific temporal AR1 + spatial + spatio-temporal type 1)
  stk.est_18 <- inla.stack(
    tag = "est",
    data = list(y = mieren_data$is_present),
    A = list(1, A_spatial, A),
    effects = list(
      data.frame(
        hab_type = mieren_data$habitat_type,
        is_winter = mieren_data$is_winter,
        year_group = mieren_data$year_group,
        hab_type_numeric = mieren_data$hab_type_numeric),
      list(spatial_field = spatial.index.only$spatial_field),
      list(s = indexs$s, s.group = indexs$s.group)))
  
  stk.pred_18 <- inla.stack(
    tag = "pred",
    data = list(y = NA),
    A = list(1, Ap_spatial, Ap),
    effects = list(
      data.frame(
        hab_type = pred_hab_type,
        is_winter = pred_is_winter,
        year_group = pred_year_group,
        hab_type_numeric = pred_hab_type_numeric),
      list(spatial_field = indexs_pred_spatial$spatial_field),
      list(s = indexs_pred$s, s.group = indexs_pred$s.group)))
  
  stk.full_18 <- inla.stack(stk.est_18, stk.pred_18)
  
  # Stack for formula19 (habitat-specific temporal AR1 + spatial + spatio-temporal type 2)
  stk.full_19 <- stk.full_18
  
  # Stack for formula20 (habitat-specific temporal AR1 + spatial + spatio-temporal type 3)
  stk.full_20 <- stk.full_18
  
  # Stack for formula21 (habitat-specific temporal AR1 + spatial + spatio-temporal type 4)
  stk.full_21 <- stk.full_18
  
  # Stacks for habitat-specific models with interaction (formulas 22-27)
  
  # Stack for formula22 (habitat-specific temporal AR1 + interaction)
  stk.est_22 <- inla.stack(
    tag = "est",
    data = list(y = mieren_data$is_present),
    A = list(1),
    effects = list(
      data.frame(
        hab_type = mieren_data$habitat_type,
        is_winter = mieren_data$is_winter,
        jaar1 = mieren_data$year_group,
        year_group = mieren_data$year_group,
        hab_type_numeric = mieren_data$hab_type_numeric)))
  
  stk.pred_22 <- inla.stack(
    tag = "pred",
    data = list(y = NA),
    A = list(1),
    effects = list(
      data.frame(
        hab_type = pred_hab_type,
        is_winter = pred_is_winter,
        jaar1 = pred_year_group,
        year_group = pred_year_group,
        hab_type_numeric = pred_hab_type_numeric)))
  
  stk.full_22 <- inla.stack(stk.est_22, stk.pred_22)
  
  # Stack for formula23 (habitat-specific temporal AR1 + interaction + spatial)
  stk.est_23 <- inla.stack(
    tag = "est",
    data = list(y = mieren_data$is_present),
    A = list(1, A_spatial),
    effects = list(
      data.frame(
        hab_type = mieren_data$habitat_type,
        is_winter = mieren_data$is_winter,
        jaar1 = mieren_data$year_group,
        year_group = mieren_data$year_group,
        hab_type_numeric = mieren_data$hab_type_numeric),
      list(spatial_field = spatial.index.only$spatial_field)))
  
  stk.pred_23 <- inla.stack(
    tag = "pred",
    data = list(y = NA),
    A = list(1, Ap_spatial),
    effects = list(
      data.frame(
        hab_type = pred_hab_type,
        is_winter = pred_is_winter,
        jaar1 = pred_year_group,
        year_group = pred_year_group,
        hab_type_numeric = pred_hab_type_numeric),
      list(spatial_field = indexs_pred_spatial$spatial_field)))
  
  stk.full_23 <- inla.stack(stk.est_23, stk.pred_23)
  
  # Stack for formula24 (habitat-specific temporal AR1 + interaction + spatial + spatio-temporal type 1)
  stk.est_24 <- inla.stack(
    tag = "est",
    data = list(y = mieren_data$is_present),
    A = list(1, A_spatial, A),
    effects = list(
      data.frame(
        hab_type = mieren_data$habitat_type,
        is_winter = mieren_data$is_winter,
        jaar1 = mieren_data$year_group,
        year_group = mieren_data$year_group,
        hab_type_numeric = mieren_data$hab_type_numeric),
      list(spatial_field = spatial.index.only$spatial_field),
      list(s = indexs$s, s.group = indexs$s.group)))
  
  stk.pred_24 <- inla.stack(
    tag = "pred",
    data = list(y = NA),
    A = list(1, Ap_spatial, Ap),
    effects = list(
      data.frame(
        hab_type = pred_hab_type,
        is_winter = pred_is_winter,
        jaar1 = pred_year_group,
        year_group = pred_year_group,
        hab_type_numeric = pred_hab_type_numeric),
      list(spatial_field = indexs_pred_spatial$spatial_field),
      list(s = indexs_pred$s, s.group = indexs_pred$s.group)))
  
  stk.full_24 <- inla.stack(stk.est_24, stk.pred_24)
  
  # Stack for formula25 (habitat-specific temporal AR1 + interaction + spatial + spatio-temporal type 2)
  stk.full_25 <- stk.full_24
  
  # Stack for formula26 (habitat-specific temporal AR1 + interaction + spatial + spatio-temporal type 3)
  stk.full_26 <- stk.full_24
  
  # Stack for formula27 (habitat-specific temporal AR1 + interaction + spatial + spatio-temporal type 4)
  stk.full_27 <- stk.full_24
  
  # Stack for formula28 (habitat-specific temporal AR1 exchangeable)
  stk.est_28 <- inla.stack(
    tag = "est",
    data = list(y = mieren_data$is_present),
    A = list(1),
    effects = list(
      data.frame(
        hab_type = mieren_data$habitat_type,
        is_winter = mieren_data$is_winter,
        year_group = mieren_data$year_group,
        hab_type_numeric = mieren_data$hab_type_numeric)))
  
  stk.pred_28 <- inla.stack(
    tag = "pred",
    data = list(y = NA),
    A = list(1),
    effects = list(
      data.frame(
        hab_type = pred_hab_type,
        is_winter = pred_is_winter,
        year_group = pred_year_group,
        hab_type_numeric = pred_hab_type_numeric)))
  
  stk.full_28 <- inla.stack(stk.est_28, stk.pred_28)
  
  # Stack for formula29 (habitat-specific temporal AR1 exchangeable + spatial)
  stk.est_29 <- inla.stack(
    tag = "est",
    data = list(y = mieren_data$is_present),
    A = list(1, A_spatial),
    effects = list(
      data.frame(
        hab_type = mieren_data$habitat_type,
        is_winter = mieren_data$is_winter,
        year_group = mieren_data$year_group,
        hab_type_numeric = mieren_data$hab_type_numeric),
      list(spatial_field = spatial.index.only$spatial_field)))
  
  stk.pred_29 <- inla.stack(
    tag = "pred",
    data = list(y = NA),
    A = list(1, Ap_spatial),
    effects = list(
      data.frame(
        hab_type = pred_hab_type,
        is_winter = pred_is_winter,
        year_group = pred_year_group,
        hab_type_numeric = pred_hab_type_numeric),
      list(spatial_field = indexs_pred_spatial$spatial_field)))
  
  stk.full_29 <- inla.stack(stk.est_29, stk.pred_29)
  
  # Stack for formula30 (habitat-specific temporal AR1 exchangeable + spatial + spatio-temporal type 1)
  stk.est_30 <- inla.stack(
    tag = "est",
    data = list(y = mieren_data$is_present),
    A = list(1, A_spatial, A),
    effects = list(
      data.frame(
        hab_type = mieren_data$habitat_type,
        is_winter = mieren_data$is_winter,
        year_group = mieren_data$year_group,
        hab_type_numeric = mieren_data$hab_type_numeric),
      list(spatial_field = spatial.index.only$spatial_field),
      list(s = indexs$s, s.group = indexs$s.group)))
  
  stk.pred_30 <- inla.stack(
    tag = "pred",
    data = list(y = NA),
    A = list(1, Ap_spatial, Ap),
    effects = list(
      data.frame(
        hab_type = pred_hab_type,
        is_winter = pred_is_winter,
        year_group = pred_year_group,
        hab_type_numeric = pred_hab_type_numeric),
      list(spatial_field = indexs_pred_spatial$spatial_field),
      list(s = indexs_pred$s, s.group = indexs_pred$s.group)))
  
  stk.full_30 <- inla.stack(stk.est_30, stk.pred_30)
  
  # Stack for formula31 (habitat-specific temporal AR1 + interaction + spatial + spatio-temporal type 2)
  stk.full_31 <- stk.full_30
  # Stack for formula31 (habitat-specific temporal AR1 + interaction + spatial + spatio-temporal type 3)
  stk.full_32 <- stk.full_30
  # Stack for formula33 (habitat-specific temporal AR1 + interaction + spatial + spatio-temporal type 4)
  stk.full_33 <- stk.full_30
  
  # statistical analysis: fit the models ---------------------------------------
  fit_models <- function(formula, stack, model_name, 
                         families = c("binomial")) {
    models <- list()
    
    for(family in families) {
      model_key <- paste(model_name, family, sep = "_")
      
      if(check_existing_model(curr_species, model_name, family)) {
        cat(paste("Model", model_key, "for species", curr_species, 
                  "already exists"))
        
        tryCatch({
          existing_model <- readRDS(file.path(dir_path, 
                                              paste0(curr_species, "_", 
                                                     model_name, "_", family, 
                                                     ".rds")))
          models[[model_key]] <- existing_model
          cat(paste("Successfully loaded existing model:", model_key, "\n"))
          next  
        }, error = function(e) {
          cat(paste("Error loading existing model", model_key, ":", e$message, 
                    "\n"))
          cat("Will attempt to refit the model")
        })
      }
      
      
      cat(paste("Fitting model:", model_key, "for species", curr_species, "\n"))
      
      tryCatch({
        model <- inla(formula = formula, 
                      family = family,  
                      data = inla.stack.data(stack),
                      control.predictor = list(compute = TRUE, 
                                               A = inla.stack.A(stack)),
                      control.compute = list(dic = TRUE, waic = TRUE, 
                                             cpo = TRUE,
                                             openmp.strategy = "huge"),
                      safe = FALSE, num.threads = parallel::detectCores() - 1)
        
        models[[model_key]] <- model
        saveRDS(model, file = file.path(dir_path, paste0(curr_species, "_", 
                                                         model_name, "_",
                                                         family, ".rds")))
        cat(paste("Completed", model_name, "with", family, "family", "\n"))
        
      }, error = function(e) {
        model <- NULL
        cat(paste("Error fitting model", model_key, ":", e$message, "\n"))
      })
    }
    
    return(models)
  }
  
  all_models <- list()
  all_models <- c(all_models, fit_models(formula1, stk.full_1, "formula1"))
  all_models <- c(all_models, fit_models(formula2, stk.full_2, "formula2"))
  all_models <- c(all_models, fit_models(formula3, stk.full_3, "formula3"))
  all_models <- c(all_models, fit_models(formula4, stk.full_4, "formula4"))
  all_models <- c(all_models, fit_models(formula5, stk.full_5, "formula5"))
  all_models <- c(all_models, fit_models(formula6, stk.full_6, "formula6"))
  all_models <- c(all_models, fit_models(formula7, stk.full_7, "formula7"))
  all_models <- c(all_models, fit_models(formula8, stk.full_8, "formula8"))
  all_models <- c(all_models, fit_models(formula9, stk.full_9, "formula9"))
  all_models <- c(all_models, fit_models(formula10, stk.full_10, "formula10"))
  all_models <- c(all_models, fit_models(formula11, stk.full_11, "formula11"))
  all_models <- c(all_models, fit_models(formula12, stk.full_12, "formula12"))
  all_models <- c(all_models, fit_models(formula13, stk.full_13, "formula13"))
  all_models <- c(all_models, fit_models(formula14, stk.full_14, "formula14"))
  all_models <- c(all_models, fit_models(formula15, stk.full_15, "formula15"))
  all_models <- c(all_models, fit_models(formula16, stk.full_16, "formula16"))
  all_models <- c(all_models, fit_models(formula17, stk.full_17, "formula17"))
  all_models <- c(all_models, fit_models(formula18, stk.full_18, "formula18"))
  all_models <- c(all_models, fit_models(formula19, stk.full_19, "formula19"))
  all_models <- c(all_models, fit_models(formula20, stk.full_20, "formula20"))
  all_models <- c(all_models, fit_models(formula21, stk.full_21, "formula21"))
  all_models <- c(all_models, fit_models(formula22, stk.full_22, "formula22"))
  all_models <- c(all_models, fit_models(formula23, stk.full_23, "formula23"))
  all_models <- c(all_models, fit_models(formula24, stk.full_24, "formula24"))
  all_models <- c(all_models, fit_models(formula25, stk.full_25, "formula25"))
  all_models <- c(all_models, fit_models(formula26, stk.full_26, "formula26"))
  all_models <- c(all_models, fit_models(formula27, stk.full_27, "formula27"))
  all_models <- c(all_models, fit_models(formula28, stk.full_28, "formula28"))
  all_models <- c(all_models, fit_models(formula29, stk.full_29, "formula29"))
  all_models <- c(all_models, fit_models(formula30, stk.full_30, "formula30"))
  all_models <- c(all_models, fit_models(formula31, stk.full_31, "formula31"))
  all_models <- c(all_models, fit_models(formula32, stk.full_32, "formula32"))
  all_models <- c(all_models, fit_models(formula33, stk.full_33, "formula33"))
  
  converged_models <- all_models[!sapply(all_models, is.null)]
  
  results_list[[curr_species]] <- list(
    models = converged_models,
    stacks = list(
      formula1 = stk.full_1,
      formula2 = stk.full_2,
      formula3 = stk.full_3,
      formula4 = stk.full_4,
      formula5 = stk.full_5,
      formula6 = stk.full_6,
      formula7 = stk.full_7,
      formula8 = stk.full_8,
      formula9 = stk.full_9,
      formula10 = stk.full_10,
      formula11 = stk.full_11,
      formula12 = stk.full_12,
      formula13 = stk.full_13,
      formula14 = stk.full_14,
      formula15 = stk.full_15,
      formula16 = stk.full_16,
      formula17 = stk.full_17,
      formula18 = stk.full_18,
      formula19 = stk.full_19,
      formula20 = stk.full_20,
      formula21 = stk.full_21,
      formula22 = stk.full_22,
      formula23 = stk.full_23,
      formula24 = stk.full_24,
      formula25 = stk.full_25,
      formula26 = stk.full_26,
      formula27 = stk.full_27,
      formula28 = stk.full_28,
      formula29 = stk.full_29,
      formula30 = stk.full_30,
      formula31 = stk.full_31,
      formula32 = stk.full_32,
      formula33 = stk.full_33),
    prediction_points = prediction_points_clean,
    prediction_grid = grid_habitats_clean,
    n_prediction_cells = n_cells_clean)
  
  # statistical analysis: model comparison -------------------------------------
  
  model_comparison <- data.frame(
    model = names(all_models),
    DIC = sapply(all_models, function(x) x$dic$dic),
    WAIC = sapply(all_models, function(x) x$waic$waic),
    marg_likelihood = sapply(all_models, function(x) x$mlik[1]),
    LPML = sapply(all_models, function(x) sum(log(x$cpo$cpo), na.rm = TRUE)),
    stringsAsFactors = FALSE)
  write.csv(model_comparison, 
            file = file.path(dir_path, paste0(curr_species, 
                                              "_model_comparison.csv")),
            row.names = FALSE)
  # Select best model based on WAIC (lower is better)
  if (best_model_name == ""){
    best_index <- which.min(model_comparison$WAIC)
    best_model_name <- model_comparison$model[best_index]
  } else {
    best_index <- which(model_comparison$model == best_model_name)
  }
  write.csv(model_comparison, 
            file = file.path(dir_path, paste0(curr_species, 
                                              "_model_comparison.csv")),
            row.names = FALSE)
  
  ## Select best model 
  best_model <- all_models[[best_model_name]]
  info_best_model <- strsplit(best_model_name, "_")[[1]]
  formula_name <- info_best_model[1]
  family_name <- info_best_model[2]
  best_stack <- results_list[[curr_species]]$stacks[[formula_name]]
  
  
  print(paste0("Selected best model: ", best_model_name, " for ", curr_species))
  print(paste0("Best model WAIC: ", round(model_comparison$WAIC[best_index], 
                                          2)))
  
  # statistical analysis: extracting best model results ------------------------
  best_model_info <- list(
    species = curr_species,
    best_model = best_model,
    best_model_name = best_model_name,
    formula_number = formula_name,
    family = family_name,
    waic = model_comparison$WAIC[best_index],
    dic = model_comparison$DIC[best_index], 
    marginal_likelihood = model_comparison$marg_likelihood[best_index],
    cpo = model_comparison$CPO[best_index],
    pit = model_comparison$PIT[best_index],
    best_stack = best_stack,
    mesh = mesh,
    spde_original = spde, 
    mieren_data = mieren_data,
    prediction_points = prediction_points_clean,
    original_priors = list(
      spde_range = c(30000, 0.5),
      spde_sigma = c(1.0, 0.5),
      ar1_prior = list(theta = list(prior = "pccor1", param = c(0, 0.5)))))
  saveRDS(best_model_info, file = file.path(dir_path, 
                                            paste0(curr_species, 
                                                   "_best_model_info.rds")))
  rm(all_models, mieren_data, coo, A, Ap, A_spatial, 
     Ap_spatial)
  rm(indexs, indexs_pred, indexs_pred_spatial, spatial.index.only)
  rm(pred_hab_type, pred_hab_type_numeric, pred_year_group, pred_jaar1, 
     pred_jaar1_hab, pred_is_winter)
  rm(prediction_points_clean, year_indices)
  gc()
  print(paste0("Spatio-temporal analysis for ", curr_species, " completed"))
}

saveRDS(results_list, file = file.path(dir_path, "complete_results_list.rds"))
print("Spatio-temporal analysis done for selected species")

# General setup: inference model -----------------------------------------------
## Get appropriate information of the selected model ---------------------------
result_list <- readRDS(file.path(dir_path, "complete_results_list.rds"))

for (species in species_group){
  ### Load general information
  best_info <- readRDS(file.path(dir_path, paste0(species, 
                                                  "_best_model_info.rds")))
  mesh <- best_info$mesh
  spde <- best_info$spde_original
  priors <- best_info$original_priors
  prediction_points <- best_info$prediction_points
  mieren_data <- best_info$mieren_data
  
  ### Load model specific information
  best_model <- result_list[[species]]$models[[best_model_name]]
  best_stack <- result_list[[species]]$stacks[[best_model_formula]]
  
  final_model_info <- list(
    species = species,
    best_model = best_model,
    best_model_name = best_model_name,
    formula_name = best_model_formula,
    family = best_model_family,
    best_stack = best_stack,
    mesh = mesh,
    priors = priors,
    spde_original = spde,
    prediction_points = prediction_points,
    mieren_data = mieren_data)
  saveRDS(final_model_info, file = file.path(dir_path, 
                                             paste0(species, 
                                                    "_final_model_info.rds")))
}


rm(result_list)
gc()


# Prior sensitivity analysis for the best model --------------------------------
## General setup: specification of priors --------------------------------------
all_priors <- list(
  original_info = list(
    name = "original_info",
    spde_range = c(5000, 0.5), 
    spde_sigma = c(1, 0.01),     
    ar1_param = c(0, 0.9)),
  moderate_info = list(
    name = "moderate_info",
    spde_range = c(7000, 0.5), 
    spde_sigma = c(1.5, 0.05),     
    ar1_param = c(0, 0.7)),
  weak_info = list(
    name = "weak_info", 
    spde_range = c(10000, 0.3),   
    spde_sigma = c(2, 0.1),     
    ar1_param = c(0, 0.5)),
  vague_info = list(
    name = "vague_info",
    spde_range = c(15000, 0.2),   
    spde_sigma = c(3, 0.2),     
    ar1_param = c(0, 0.3)),
  very_vague_info = list(
    name = "very_vague_info",
    spde_range = c(25000, 0.1),   
    spde_sigma = c(5, 0.3),     
    ar1_param = c(0, 0.1)),
  ultra_vague_info = list(
    name = "ultra_vague_info",
    spde_range = c(35000, 0.05),   
    spde_sigma = c(10, 0.4),     
    ar1_param = c(0, 0.05)),
  no_info = list(
    name = "no_info",
    spde_range = c(50000, 0.01),   
    spde_sigma = c(20, 0.5),     
    ar1_param = c(0, 0.01)))

## General setup: create main prior sensitivity analysis directory -------------

sensitivity_dir <- file.path(dir_path, "sensitivity_analysis")
if (!dir.exists(sensitivity_dir)) {
  dir.create(sensitivity_dir, recursive = TRUE)
}

## Prior analysis --------------------------------------------------------------

for (species in species_group) {
  fitted_models <- list() 
  waic_comparison <- data.frame()
  best_info <- readRDS(file.path(dir_path, paste0(species, 
                                                  "_final_model_info.rds")))
  curr_species <- best_info$species
  formula_name <- best_info$formula_name
  fitted_models[[curr_species]] <- list()
  for (some_prior in names(all_priors)) {
    ### Extracting the right data + updating priors
    curr_prior <- all_priors[[some_prior]]
    model_file <- file.path(sensitivity_dir, 
                            paste0(curr_species, "_", some_prior, 
                                   "_model.rds"))
    
    if (file.exists(model_file)) {
      sensitivity_model <- readRDS(model_file)
      lpml_value <- sum(log(sensitivity_model$cpo$cpo), na.rm = TRUE)
      
      fitted_models[[curr_species]][[some_prior]] <- sensitivity_model
      
      # Add to WAIC comparison table
      waic_comparison <- rbind(waic_comparison, data.frame(
        species = curr_species,
        prior_type = some_prior,
        original_formula = best_info$best_model_name,
        waic = sensitivity_model$waic$waic,
        lpml = lpml_value,
        dic = sensitivity_model$dic$dic,
        stringsAsFactors = FALSE))
      
    } else {
      new_spde <- inla.spde2.pcmatern(
        mesh = best_info$mesh, 
        alpha = 2, 
        constr = TRUE,
        prior.range = curr_prior$spde_range,
        prior.sigma = curr_prior$spde_sigma)
      new_ar1_prior <- list(theta = list(prior = "pccor1", 
                                         param = curr_prior$ar1_param))
      
      ### Data preparation -----------------------------------------------------
      mieren_data <- best_info$mieren_data
      mieren_sf <- st_as_sf(mieren_data, coords = c("X", "Y"), crs = 31370)
      mieren_sf$X <- st_coordinates(mieren_sf)[,1]
      mieren_sf$Y <- st_coordinates(mieren_sf)[,2]
      search_eurostat("nuts")
      nuts_sf <- get_eurostat_geospatial(output_class = "sf", resolution = "01", 
                                         nuts_level = 2)
      belgium_nuts <- nuts_sf[nuts_sf$CNTR_CODE == "BE", ]
      limburg <- belgium_nuts[belgium_nuts$NUTS_NAME == "Prov. Limburg (BE)", ] 
      limburg_utm <- st_transform(limburg, crs = 31370)
      
      ### Statistical analysis: set up -----------------------------------------
      #### Mesh construction 
      mesh <- best_info$mesh
      #### Index set construction
      timesn <- length(unique(mieren_data$jaar1))
      indexs <- inla.spde.make.index("s", n.spde = new_spde$n.spde, 
                                     n.group = timesn)
      lengths(indexs)
      
      #### Projection matrix construction 
      A <- inla.spde.make.A(mesh = mesh,
                            loc = cbind(mieren_data$X, mieren_data$Y),
                            group = as.integer(factor(mieren_data$jaar1)))
      
      #### Prediction location construction
      pred_grid_path <- file.path(dir_path, "pred_grid_list.rds")
      if (file.exists(pred_grid_path)) {
        pred_grid <- readRDS(pred_grid_path)
        grid_habitats_clean <- pred_grid$grid_habitats_clean
        n_cells_clean <- pred_grid$n_cells_clean
      } else {
        bb <- st_bbox(limburg_utm)
        grid <- st_make_grid(
          limburg_utm,
          cellsize = resolution_grid,
          what = "polygons",
          square = TRUE) |> st_as_sf() |> st_intersection(limburg_utm)
        
        grid_predictions <- grid %>% mutate(cell_id = 1:n()) %>% 
          st_set_crs(st_crs(limburg_utm))
        grid_predictions <- st_transform(grid_predictions, st_crs(limburg_utm))
        BWK_map <- st_transform(BWK_map, st_crs(limburg_utm))
        
        dominant_habitats <- grid_predictions %>% st_intersection(BWK_map) %>%
          mutate(habitat_type_dutch = sapply(BWKLABEL, function(x) classify_habitat(x)[1]),
                 habitat_type = group_and_translate_habitat(habitat_type_dutch),
                 fragment_area = st_area(.)) %>%  st_drop_geometry() %>%
          group_by(cell_id, habitat_type) %>% summarise(area = sum(fragment_area), 
                                                        .groups = "drop_last") %>%
          mutate(percentage = as.numeric(area / sum(area) * 100)) %>% ungroup() %>%
          group_by(cell_id) %>% slice_max(percentage, n = 1, with_ties = FALSE) %>%
          dplyr::select(cell_id, dominant_habitat = habitat_type, 
                        dominant_percentage = percentage)
        
        grid_habitats <- grid_predictions %>%
          left_join(dominant_habitats, by = "cell_id")
        na_cells <- is.na(grid_habitats$dominant_habitat)
        grid_habitats_clean <- grid_habitats[!na_cells, ]
        unique_obs_habitats <- levels(mieren_data$habitat_type)
        valid_habitat <- grid_habitats_clean$dominant_habitat %in% 
          unique_obs_habitats
        grid_habitats_clean <- grid_habitats_clean[valid_habitat, ]
        n_cells_clean <- nrow(grid_habitats_clean)
        grid_habitats_clean$
          dominant_habitat <- factor(grid_habitats_clean$dominant_habitat, 
                                     levels = levels(mieren_data$habitat_type))
        saveRDS(list(grid_habitats_clean = grid_habitats_clean,
                     n_cells_clean = n_cells_clean), pred_grid_path)
        rm(grid, grid_predictions, dominant_habitats, grid_habitats)
        gc()
      }
      observed_years <- sort(unique(mieren_data$jaar1))
      observed_year_groups <- sort(unique(mieren_data$year_group))
      n_years <- length(observed_years)
      
      expanded_grid_clean <- do.call(rbind, lapply(1:n_years, function(i) {
        cbind(grid_habitats_clean, 
              year = observed_years[i],
              year_group = observed_year_groups[i])
      }))
      
      
      #### Spatial index construction 
      spatial.index <- inla.spde.make.index(name = "s", 
                                            n.spde = new_spde$n.spde,       
                                            n.group = length(unique(mieren_data$
                                                                      jaar1)))
      
      spatial.index.only <- inla.spde.make.index(name = "spatial_field", 
                                                 n.spde = new_spde$n.spde)
      A_spatial <- inla.spde.make.A(mesh, loc = cbind(mieren_data$X,
                                                      mieren_data$Y))
      
      
      ### Prepare prediction variables
      pred_hab_type <- factor(rep(grid_habitats_clean$dominant_habitat,
                                  n_years), 
                              levels = levels(mieren_data$habitat_type))
      pred_hab_type_numeric <- as.numeric(as.factor(pred_hab_type))
      pred_year_group <- rep(observed_year_groups, each = n_cells_clean)
      pred_jaar1 <- rep(observed_years, each = n_cells_clean)
      pred_jaar1_hab <- interaction(pred_year_group, pred_hab_type, drop = TRUE)
      pred_is_winter <- factor(rep("no winter", n_cells_clean * n_years), 
                               levels = levels(mieren_data$is_winter))
      
      ### Create temporal indices, data frame and projection matrix
      indexs_pred <- inla.spde.make.index("s", n.spde = mesh$n, 
                                          n.group = n_years)
      prediction_points_clean <- st_centroid(expanded_grid_clean)
      saveRDS(prediction_points_clean, 
              file = file.path(dir_path, paste0(curr_species, 
                                                "_prediction_points.rds")))
      year_indices <- match(expanded_grid_clean$year, observed_years)
      Ap <- inla.spde.make.A(mesh = mesh,
                             loc = st_coordinates(prediction_points_clean),
                             group = year_indices, fast = TRUE)
      
      indexs_pred_spatial <- inla.spde.make.index(name = "spatial_field", 
                                                  n.spde = mesh$n)
      Ap_spatial <- inla.spde.make.A(mesh = mesh, 
                                     loc = st_coordinates(prediction_points_clean),
                                     fast = TRUE)
      
      ### statistical analysis: model specification ----------------------------
      rprior <- new_ar1_prior
      
      
      if(formula_name == "formula20"){
        curr_formula <- y ~ is_winter + hab_type + 
          f(inla.group(year_group, n = length(unique(mieren_data$year_group))), 
            model = "ar1", replicate = hab_type_numeric, hyper = rprior) + 
          f(spatial_field, model = new_spde) + 
          f(s, model = new_spde, group = s.group, 
            control.group = list(model = "iid"))
        
        stk.est_20 <- inla.stack(
          tag = "est",
          data = list(y = mieren_data$is_present),
          A = list(1, A_spatial, A),
          effects = list(
            data.frame(
              hab_type = mieren_data$habitat_type,
              is_winter = mieren_data$is_winter,
              year_group = mieren_data$year_group,
              hab_type_numeric = mieren_data$hab_type_numeric),
            list(spatial_field = spatial.index.only$spatial_field),
            list(s = indexs$s, s.group = indexs$s.group)))
        
        stk.pred_20 <- inla.stack(
          tag = "pred",
          data = list(y = NA),
          A = list(1, Ap_spatial, Ap),
          effects = list(
            data.frame(
              hab_type = pred_hab_type,
              is_winter = pred_is_winter,
              year_group = pred_year_group,
              hab_type_numeric = pred_hab_type_numeric),
            list(spatial_field = indexs_pred_spatial$spatial_field),
            list(s = indexs_pred$s, s.group = indexs_pred$s.group)))
        
        new_stk.full <- inla.stack(stk.est_20, stk.pred_20)
      }
      
      ### statistical analysis: model fitting ----------------------------------
      sensitivity_model <- try({
        inla(
          formula = curr_formula,
          family = best_info$family,
          data = inla.stack.data(new_stk.full),
          control.predictor = list(compute = TRUE, 
                                   A = inla.stack.A(new_stk.full)),
          control.compute = list(dic = TRUE, waic = TRUE, 
                                 cpo = TRUE,
                                 openmp.strategy = "huge"),
          safe = FALSE, num.threads = parallel::detectCores() - 1)
      })
      
      if (inherits(sensitivity_model, "try-error")) {
        cat(" FAILED\n")
        next
      }
      
      saveRDS(sensitivity_model, model_file)
      lpml_value <- sum(log(sensitivity_model$cpo$cpo), na.rm = TRUE)
      
      fitted_models[[curr_species]][[some_prior]] <- sensitivity_model
      
      # Add to WAIC comparison table
      waic_comparison <- rbind(waic_comparison, data.frame(
        species = curr_species,
        prior_type = some_prior,
        original_formula = best_info$best_model_name,
        waic = sensitivity_model$waic$waic,
        lpml = lpml_value,
        dic = sensitivity_model$dic$dic,
        stringsAsFactors = FALSE))
    }
  }
  write.csv(waic_comparison, file.path(sensitivity_dir,paste0(curr_species, 
                                                              "gof_table.csv")), 
            row.names = FALSE)
  
}

# General setup: inference model (only if prior set changed) -------------------
## General setup: fit inference model ------------------------------------------
results_list <- list()
for (curr_species in species_group) {
  # Data preparation -----------------------------------------------------------
  mieren_specific <- mieren %>% filter(`LATIJNSE NAAM` %in% species_group)
  tot_obs <- nrow(mieren_specific)
  presence_per_hab_specific <- mieren_specific %>%
    group_by(habitat_type) %>%summarise(n_present = sum(is_present_bin == 1, 
                                                        na.rm = TRUE),
                                        presence_rate = (n_present / tot_obs) * 
                                          100)
  # Filter for habitat types with sufficient presence OR important habitat types
  hab_analysis <- presence_per_hab_specific %>%
    filter(presence_rate > 5 | habitat_type == 'Wet heathland' 
           | habitat_type == 'Human influence' | habitat_type == 'Grassland' 
           | habitat_type == 'Forest' | habitat_type == 'Dry heathland') %>%
    pull(habitat_type) 
  
  mieren_data <- mieren %>%
    mutate(habitat_type_analysis = case_when(
      habitat_type %in% hab_analysis ~ habitat_type,
      TRUE ~ 'Other'))
  mieren_data$habitat_type <- mieren_data$habitat_type_analysis
  
  is_present <- numeric(nrow(mieren_data))
  for (i in 1:nrow(mieren_data)){
    if (mieren_data$`LATIJNSE NAAM`[i] == curr_species){
      if (mieren_data$Aantal[i] > 0){
        is_present[i] <- 1
      }
    }
  }
  
  mieren_data$is_present <- is_present
  mieren_data$jaar1_num <- as.numeric(mieren_data$jaar1)
  mieren_data$year_group <- as.integer(factor(mieren_data$jaar1))
  mieren_data$presence <- factor(ifelse(mieren_data$is_present == 1, 
                                        "present", "absent"), 
                                 levels = c("absent", "present"))
  mieren_data <- mieren_data %>% filter(!is.na(habitat_type))
  mieren_data$habitat_type <- as.factor(mieren_data$habitat_type)
  mieren_data$hab_type_numeric <- as.numeric(as.factor(mieren_data$
                                                         habitat_type))
  mieren_data$is_winter <- factor(ifelse(mieren_data$is_winter_bin == 1, 
                                         "winter", "no winter"), 
                                  levels = c("no winter", "winter"))
  mieren_data$jaar1_hab <- interaction(mieren_data$year_group, 
                                       mieren_data$habitat_type, drop = TRUE)
  mieren_data$jaar1_num <- as.numeric(mieren_data$jaar1)
  mieren_data$year_group_global <- mieren_data$year_group
  mieren_data$year_group_habitat <- mieren_data$year_group
  
  # Data vizualization 
  if (!exists("limburg_utm_cached")) {
    search_eurostat("nuts")
    nuts_sf <- get_eurostat_geospatial(output_class = "sf", resolution = "01", 
                                       nuts_level = 2)
    belgium_nuts <- nuts_sf[nuts_sf$CNTR_CODE == "BE", ]
    limburg <- belgium_nuts[belgium_nuts$NUTS_NAME == "Prov. Limburg (BE)", ] 
    limburg_utm_cached <<- st_transform(limburg, crs = 31370)
    rm(nuts_sf, belgium_nuts, limburg)
    gc()
  }
  
  limburg_utm <- limburg_utm_cached
  print(paste0("data preparation of species ", curr_species, " finished."))
  
  # Statistical analysis: set up -----------------------------------------------
  ## Mesh construction
  coo <- cbind(mieren_data$X, mieren_data$Y)
  mesh_file_path <- file.path(dir_path, "mesh.rds")
  if (file.exists(mesh_file_path)) {
    mesh_list <- readRDS(mesh_file_path)
    bnd <- mesh_list$bnd
    mesh <- mesh_list$mesh
    spde <- mesh_list$spde
  } else {
    bnd <- inla.nonconvex.hull(st_coordinates(limburg_utm)[,1:2])
    mesh <- inla.mesh.2d(loc = coo, boundary = bnd, max.edge = c(20000, 40000),  
                         cutoff = 2000)
    spde <- inla.spde2.pcmatern(mesh = mesh, alpha = 2, constr = TRUE,
                                prior.range = spde_range, 
                                prior.sigma = spde_sigma)
    
    saveRDS(list(bnd = bnd, mesh = mesh, spde = spde), mesh_file_path)
  }
  ## Index set construction
  timesn <- length(unique(mieren_data$jaar1))
  indexs <- inla.spde.make.index("s", n.spde = spde$n.spde, n.group = timesn)
  lengths(indexs)
  
  ## Projection matrix construction 
  A <- inla.spde.make.A(
    mesh = mesh,
    loc = cbind(mieren_data$X, mieren_data$Y),
    group = as.integer(factor(mieren_data$jaar1))
  )
  
  ## Prediction location construction
  pred_grid_path <- file.path(dir_path, "pred_grid_list.rds")
  if (file.exists(pred_grid_path)) {
    pred_grid <- readRDS(pred_grid_path)
    grid_habitats_clean <- pred_grid$grid_habitats_clean
    n_cells_clean <- pred_grid$n_cells_clean
  } else {
    bb <- st_bbox(limburg_utm)
    grid <- st_make_grid(
      limburg_utm,
      cellsize = resolution_grid,
      what = "polygons",
      square = TRUE) |> st_as_sf() |> st_intersection(limburg_utm)
    
    grid_predictions <- grid %>% mutate(cell_id = 1:n()) %>% 
      st_set_crs(st_crs(limburg_utm))
    
    # Process dominant habitats
    grid_predictions <- st_transform(grid_predictions, st_crs(limburg_utm))
    BWK_map <- st_transform(BWK_map, st_crs(limburg_utm))
    
    dominant_habitats <- grid_predictions %>% st_intersection(BWK_map) %>%
      mutate(habitat_type_dutch = sapply(BWKLABEL, function(x) classify_habitat(x)[1]),
             habitat_type = group_and_translate_habitat(habitat_type_dutch),
             fragment_area = st_area(.)) %>%  st_drop_geometry() %>%
      group_by(cell_id, habitat_type) %>% summarise(area = sum(fragment_area), 
                                                    .groups = "drop_last") %>%
      mutate(percentage = as.numeric(area / sum(area) * 100)) %>% ungroup() %>%
      group_by(cell_id) %>% slice_max(percentage, n = 1, with_ties = FALSE) %>%
      dplyr::select(cell_id, dominant_habitat = habitat_type, 
                    dominant_percentage = percentage)
    
    grid_habitats <- grid_predictions %>%
      left_join(dominant_habitats, by = "cell_id")
    na_cells <- is.na(grid_habitats$dominant_habitat)
    grid_habitats_clean <- grid_habitats[!na_cells, ]
    
    ## Filter for habitat types that are observed
    unique_obs_habitats <- levels(mieren_data$habitat_type)
    valid_habitat <- grid_habitats_clean$dominant_habitat %in% 
      unique_obs_habitats
    grid_habitats_clean <- grid_habitats_clean[valid_habitat, ]
    n_cells_clean <- nrow(grid_habitats_clean)
    grid_habitats_clean$
      dominant_habitat <- factor(grid_habitats_clean$dominant_habitat, 
                                 levels = levels(mieren_data$habitat_type))
    saveRDS(list(grid_habitats_clean = grid_habitats_clean,
                 n_cells_clean = n_cells_clean), pred_grid_path)
    rm(grid, grid_predictions, dominant_habitats, grid_habitats)
    gc()
  }
  
  ## Prediction habitat specification
  observed_years <- sort(unique(mieren_data$jaar1))
  observed_year_groups <- sort(unique(mieren_data$year_group))
  n_years <- length(observed_years)
  
  expanded_grid_clean <- do.call(rbind, lapply(1:n_years, function(i) {
    cbind(grid_habitats_clean, 
          year = observed_years[i],
          year_group = observed_year_groups[i])
  }))
  
  
  ## Spatial index construction 
  spatial.index <- inla.spde.make.index(name = "s", 
                                        n.spde = spde$n.spde,       
                                        n.group = length(unique(mieren_data$
                                                                  jaar1)))
  
  spatial.index.only <- inla.spde.make.index(name = "spatial_field", 
                                             n.spde = spde$n.spde)
  A_spatial <- inla.spde.make.A(mesh, loc = cbind(mieren_data$X, mieren_data$Y))
  
  
  ### Prepare prediction variables
  pred_hab_type <- factor(rep(grid_habitats_clean$dominant_habitat,
                              n_years), 
                          levels = levels(mieren_data$habitat_type))
  pred_hab_type_numeric <- as.numeric(as.factor(pred_hab_type))
  pred_year_group <- rep(observed_year_groups, each = n_cells_clean)
  pred_jaar1 <- rep(observed_years, each = n_cells_clean)
  pred_jaar1_hab <- interaction(pred_year_group, pred_hab_type, drop = TRUE)
  pred_is_winter <- factor(rep("no winter", n_cells_clean * n_years), 
                           levels = levels(mieren_data$is_winter))
  mieren_data$year_group_global <- mieren_data$year_group
  mieren_data$year_group_habitat <- mieren_data$year_group
  pred_year_group_global <- pred_year_group
  pred_year_group_habitat <- pred_year_group
  ### Create temporal indices, data frame and projection matrix
  indexs_pred <- inla.spde.make.index("s", n.spde = mesh$n, n.group = n_years)
  prediction_points_clean <- st_centroid(expanded_grid_clean)
  saveRDS(prediction_points_clean, 
          file = file.path(dir_path, paste0(curr_species, 
                                            "_prediction_points.rds")))
  year_indices <- match(expanded_grid_clean$year, observed_years)
  Ap <- inla.spde.make.A(mesh = mesh,
                         loc = st_coordinates(prediction_points_clean),
                         group = year_indices, fast = TRUE)
  
  indexs_pred_spatial <- inla.spde.make.index(name = "spatial_field", 
                                              n.spde = mesh$n)
  Ap_spatial <- inla.spde.make.A(mesh = mesh, 
                                 loc = st_coordinates(prediction_points_clean),
                                 fast = TRUE)
  print(paste0("Statistical set up of species ", curr_species, " finished."))
  
  # statistical analysis: model specification ----------------------------------
  rprior <- list(theta = list(prior = "pccor1", param = ar1_param))
  formula20 <- y ~ is_winter + hab_type + 
    f(inla.group(year_group, n = length(unique(mieren_data$year_group))), 
      model = "ar1", replicate = hab_type_numeric, hyper = rprior) + 
    f(spatial_field, model = spde) + 
    f(s, model = spde, group = s.group, 
      control.group = list(model = "iid"))
  
  # Stack for inference model
  stk.est_20 <- inla.stack(
    tag = "est",
    data = list(y = mieren_data$is_present),
    A = list(1, A_spatial, A),
    effects = list(
      data.frame(
        hab_type = mieren_data$habitat_type,
        is_winter = mieren_data$is_winter,
        year_group = mieren_data$year_group,
        hab_type_numeric = mieren_data$hab_type_numeric),
      list(spatial_field = spatial.index.only$spatial_field),
      list(s = indexs$s, s.group = indexs$s.group)))
  
  stk.pred_20 <- inla.stack(
    tag = "pred",
    data = list(y = NA),
    A = list(1, Ap_spatial, Ap),
    effects = list(
      data.frame(
        hab_type = pred_hab_type,
        is_winter = pred_is_winter,
        year_group = pred_year_group,
        hab_type_numeric = pred_hab_type_numeric),
      list(spatial_field = indexs_pred_spatial$spatial_field),
      list(s = indexs_pred$s, s.group = indexs_pred$s.group)))
  
  stk.full_20 <- inla.stack(stk.est_20, stk.pred_20)
  
  # statistical analysis: fit the models ---------------------------------------
  fit_models <- function(formula, stack, model_name, 
                         families = c("binomial")) {
    models <- list()
    
    for(family in families) {
      model_key <- paste(model_name, family, sep = "_")
      
      if(check_existing_model(curr_species, model_name, family)) {
        cat(paste("Model", model_key, "for species", curr_species, 
                  "already exists"))
        
        tryCatch({
          existing_model <- readRDS(file.path(dir_path, 
                                              paste0(curr_species, "_", 
                                                     model_name, "_", family, 
                                                     ".rds")))
          models[[model_key]] <- existing_model
          cat(paste("Successfully loaded existing model:", model_key, "\n"))
          next  
        }, error = function(e) {
          cat(paste("Error loading existing model", model_key, ":", e$message, 
                    "\n"))
          cat("Will attempt to refit the model")
        })
      }
      
      
      cat(paste("Fitting model:", model_key, "for species", curr_species, "\n"))
      
      tryCatch({
        model <- inla(formula = formula, 
                      family = family,  
                      data = inla.stack.data(stack),
                      control.predictor = list(compute = TRUE, 
                                               A = inla.stack.A(stack)),
                      control.compute = list(dic = TRUE, waic = TRUE, 
                                             cpo = TRUE,
                                             openmp.strategy = "huge"),
                      safe = FALSE, num.threads = parallel::detectCores() - 1)
        
        models[[model_key]] <- model
        saveRDS(model, file = file.path(dir_path, paste0(curr_species, "_", 
                                                         model_name, "_",
                                                         family, ".rds")))
        cat(paste("Completed", model_name, "with", family, "family", "\n"))
        
      }, error = function(e) {
        model <- NULL
        cat(paste("Error fitting model", model_key, ":", e$message, "\n"))
      })
    }
    
    return(models)
  }
  
  all_models <- list()
  all_models <- c(all_models, fit_models(formula20, stk.full_20, "formula20"))
  
  converged_models <- all_models[!sapply(all_models, is.null)]
  
  results_list[[curr_species]] <- list(
    models = converged_models,
    stacks = list(
      formula20 = stk.full_20),
    prediction_points = prediction_points_clean,
    prediction_grid = grid_habitats_clean,
    n_prediction_cells = n_cells_clean)
  
  # statistical analysis: model comparison -------------------------------------
  
  model_comparison <- data.frame(
    model = names(all_models),
    DIC = sapply(all_models, function(x) x$dic$dic),
    WAIC = sapply(all_models, function(x) x$waic$waic),
    marg_likelihood = sapply(all_models, function(x) x$mlik[1]),
    LPML = sapply(all_models, function(x) sum(log(x$cpo$cpo), na.rm = TRUE)),
    stringsAsFactors = FALSE)
  write.csv(model_comparison, 
            file = file.path(dir_path, paste0(curr_species, 
                                              "_model_comparison.csv")),
            row.names = FALSE)
  # Select best model based on WAIC (lower is better)
  if (best_model_name == ""){
    best_index <- which.min(model_comparison$WAIC)
    best_model_name <- model_comparison$model[best_index]
  } else {
    best_index <- which(model_comparison$model == best_model_name)
  }
  write.csv(model_comparison, 
            file = file.path(dir_path, paste0(curr_species, 
                                              "_model_comparison.csv")),
            row.names = FALSE)
  
  ## Select best model 
  best_model <- all_models[[best_model_name]]
  info_best_model <- strsplit(best_model_name, "_")[[1]]
  formula_name <- info_best_model[1]
  family_name <- info_best_model[2]
  best_stack <- results_list[[curr_species]]$stacks[[formula_name]]
  
  
  print(paste0("Selected best model: ", best_model_name, " for ", curr_species))
  print(paste0("Best model WAIC: ", round(model_comparison$WAIC[best_index], 
                                          2)))
  
  # statistical analysis: extracting best model results ------------------------
  best_model_info <- list(
    species = curr_species,
    best_model = best_model,
    best_model_name = best_model_name,
    formula_number = formula_name,
    family = family_name,
    waic = model_comparison$WAIC[best_index],
    dic = model_comparison$DIC[best_index], 
    marginal_likelihood = model_comparison$marg_likelihood[best_index],
    cpo = model_comparison$CPO[best_index],
    pit = model_comparison$PIT[best_index],
    best_stack = best_stack,
    mesh = mesh,
    spde_original = spde, 
    mieren_data = mieren_data,
    prediction_points = prediction_points_clean,
    original_priors = list(
      spde_range = c(30000, 0.5),
      spde_sigma = c(1.0, 0.5),
      ar1_prior = list(theta = list(prior = "pccor1", param = c(0, 0.5)))))
  saveRDS(best_model_info, file = file.path(dir_path, 
                                            paste0(curr_species, 
                                                   "_best_model_info.rds")))
  rm(all_models, mieren_data, coo, A, Ap, A_spatial, 
     Ap_spatial)
  rm(indexs, indexs_pred, indexs_pred_spatial, spatial.index.only)
  rm(pred_hab_type, pred_hab_type_numeric, pred_year_group, pred_jaar1, 
     pred_jaar1_hab, pred_is_winter)
  rm(prediction_points_clean, year_indices)
  gc()
  print(paste0("Spatio-temporal analysis for ", curr_species, " completed"))
}

saveRDS(results_list, file = file.path(dir_path, "complete_results_prior.rds"))

## General setup: information of inference model -------------------------------
results_list_prior <- readRDS(file.path(dir_path, "complete_results_prior.rds"))
for (species in species_group){
  best_info <- readRDS(file.path(dir_path, paste0(species, 
                                                  "_best_model_info.rds")))
  priors <- best_info$original_priors
  mesh <- best_info$mesh
  spde <- best_info$spde_original
  prediction_points <- best_info$prediction_points
  mieren_data <- best_info$mieren_data
  
  # Load model specific information
  best_model <- results_list_prior[[species]]$models[[best_model_name]]
  best_stack <- results_list_prior[[species]]$stacks[[best_model_formula]]
  
  final_model_info <- list(
    species = species,
    best_model = best_model,
    best_model_name = best_model_name,
    formula_name = best_model_formula,
    family = best_model_family,
    best_stack = best_stack,
    mesh = mesh,
    priors = priors,
    spde_original = spde,
    prediction_points = prediction_points,
    mieren_data = mieren_data)
  saveRDS(final_model_info, file = file.path(dir_path, 
                                             paste0(species, 
                                                    "_final_model_info.rds")))
}


rm(result_list_prior)
gc()


# statistical analysis: model validation ---------------------------------------
## Data preparation -----------------------------------------------------------
mieren_specific <- mieren %>% filter(`LATIJNSE NAAM` %in% species_group)
tot_obs <- nrow(mieren_specific)
presence_per_hab_specific <- mieren_specific %>%
  group_by(habitat_type) %>%summarise(n_present = sum(is_present_bin == 1, 
                                                      na.rm = TRUE),
                                      presence_rate = (n_present / tot_obs) * 
                                        100)
hab_analysis <- presence_per_hab_specific %>%
  filter(presence_rate > 5 | habitat_type == 'Wet heathland' 
         | habitat_type == 'Human influence' | habitat_type == 'Grassland' 
         | habitat_type == 'Forest' | habitat_type == 'Dry heathland') %>%
  pull(habitat_type) 

mieren_data <- mieren %>%
  mutate(habitat_type_analysis = case_when(
    habitat_type %in% hab_analysis ~ habitat_type,
    TRUE ~ 'Other'))
mieren_data$habitat_type <- mieren_data$habitat_type_analysis

is_present <- numeric(nrow(mieren_data))
curr_species <- species_group
for (i in 1:nrow(mieren_data)){
  if (mieren_data$`LATIJNSE NAAM`[i] %in% curr_species){
    if (mieren_data$Aantal[i] > 0){
      is_present[i] <- 1
    }
  }
}

mieren_data$is_present <- is_present
mieren_data$jaar1_num <- as.numeric(mieren_data$jaar1)
mieren_data$year_group <- as.integer(factor(mieren_data$jaar1))
mieren_data$presence <- factor(ifelse(mieren_data$is_present == 1, 
                                      "present", "absent"), 
                               levels = c("absent", "present"))
mieren_data <- mieren_data %>% filter(!is.na(habitat_type))
mieren_data$habitat_type <- as.factor(mieren_data$habitat_type)
mieren_data$hab_type_numeric <- as.numeric(as.factor(mieren_data$
                                                       habitat_type))
mieren_data$is_winter <- factor(ifelse(mieren_data$is_winter_bin == 1, 
                                       "winter", "no winter"), 
                                levels = c("no winter", "winter"))
mieren_data$jaar1_hab <- interaction(mieren_data$year_group, 
                                     mieren_data$habitat_type, drop = TRUE)
mieren_data$jaar1_num <- as.numeric(mieren_data$jaar1)



## spatial autocorrelation assessment (empirical variogram) --------------------
limburg_bbox <- st_bbox(limburg_utm)
max_x_distance <- limburg_bbox$xmax - limburg_bbox$xmin
max_y_distance <- limburg_bbox$ymax - limburg_bbox$ymin
max_diagonal <- sqrt(max_x_distance^2 + max_y_distance^2)
threshold_distance <- max_diagonal / 3

for (species in species_group) {
  all_vario_plots <- list()
  
  ### Load the necessary data 
  best_info <- readRDS(file.path(dir_path, paste0(species, 
                                                  "_final_model_info.rds")))
  curr_species <- best_info$species
  curr_formula <- best_info$formula_name
  curr_family <- best_info$family
  mieren_sp <- best_info$mieren_data
  
  ### Add jaar1_groep variable
  mieren_sp <- mieren_sp %>%
    mutate(jaar1_groep = cut(jaar1, breaks = seq(1995, 2025, by = 5),
                             include.lowest = TRUE, right = FALSE,
                             labels = c("1995-1999","2000-2004","2005-2009",
                                        "2010-2014", "2015-2019", "2020-2024")))
  
  ### Fit simple INLA model without spatial/temporal random effects
  best_model <- inla(is_present ~ is_winter + habitat_type,
                     family = "binomial",
                     data = mieren_sp,
                     control.predictor = list(compute = TRUE))
  
  
  ### Calculate residuals
  pred_fixed_only <- best_model$summary.fitted.values$mean
  
  mieren_sp$prediction <- pred_fixed_only  
  mieren_sp$residuals <- (mieren_sp$is_present - pred_fixed_only) / 
    sqrt(pred_fixed_only * (1 - pred_fixed_only))
  
  ### Calculate empirical variogram 
  vario_list <- list()
  all_jaar1_groep <- unique(mieren_sp$jaar1_groep)
  
  for (curr_groep in all_jaar1_groep) {
    curr_mieren_groep <- mieren_sp %>% 
      filter(jaar1_groep == curr_groep, !is.na(residuals))
    
    if (nrow(curr_mieren_groep) < 5) {
      next
    } else { 
      coordinates(curr_mieren_groep) <- ~ X + Y
      curr_vario <- variogram(residuals ~ 1, data = curr_mieren_groep)
      curr_vario_list <- list()  
      
      for (i in 1:n_perm) {
        curr_mieren_groep$resid_perm <- sample(curr_mieren_groep$residuals)
        vario_perm <- variogram(resid_perm ~ 1, data = curr_mieren_groep)
        curr_vario_list[[i]] <- vario_perm$gamma
      }
      
      curr_vario_matrix <- do.call(cbind, curr_vario_list)
      envelope_low <- apply(curr_vario_matrix, 1, quantile, 0.025, na.rm = TRUE)
      envelope_high <- apply(curr_vario_matrix, 1, quantile, 0.975, na.rm = TRUE)
      
      curr_vario_df <- data.frame(dist = curr_vario$dist, 
                                  gamma = curr_vario$gamma, 
                                  envelope_low = envelope_low, 
                                  envelope_high = envelope_high, 
                                  jaar1_groep = curr_groep, 
                                  species = curr_species)
      vario_list[[as.character(curr_groep)]] <- curr_vario_df
    }
  }
  
  ### Combine all results and plot
  if (length(vario_list) > 0) {
    all_vario <- do.call(rbind, vario_list)
    
    fig <- ggplot(all_vario, aes(x = dist, y = gamma)) +
      geom_ribbon(aes(ymin = envelope_low, ymax = envelope_high), 
                  fill = "gray80", alpha = 0.6, color = NA) +
      geom_line(color = "black", linewidth = 0.8) +
      facet_wrap(~ paste(jaar1_groep, "years"), 
                 scales = "free", ncol = 4) +
      xlim(0, threshold_distance) +
      labs(x = "Spatial distance (m)", 
           y = "Variogram",
           title = paste("Variograms for", curr_species)) + 
      theme_minimal()
    
    ggsave(filename = file.path(dir_path, paste0("All_variograms_of_", 
                                                 curr_species, ".png")), 
           plot = fig, width = width_inches, height = height_inches, dpi = 300)
  } else {
    warning(paste("No variogram data generated for", curr_species))
  }
}

## zero inflation assessment ---------------------------------------------------

zero_inflation_results <- data.frame()
for (species in species_group) {
  best_info <- readRDS(file.path(dir_path, paste0(species, 
                                                  "_final_model_info.rds")))
  curr_species <- best_info$species
  curr_formula <- best_info$formula_number
  curr_family <- best_info$family
  mieren_data <- best_info$mieren_data
  
  best_model <- readRDS(file.path(dir_path, paste0(species, "_", 
                                                   best_model_name, ".rds")))
  stack_data <- inla.stack.data(best_info$best_stack)
  est_indices <- inla.stack.index(best_info$best_stack, "est")
  pred_best_model <- best_model$summary.fitted.values$mean[est_indices$data]
  prob_absence <- 1 - pred_best_model
  expected_zeros <- sum(prob_absence, na.rm = TRUE) 
  all_observed_zeros <- mieren_data %>% filter(is_present_bin == 0)
  observed_zeros <- nrow(all_observed_zeros)
  zero_ratio <- observed_zeros / expected_zeros
  
  curr_zero_info <- data.frame(species = curr_species,
                               best_model = curr_formula,
                               best_family = curr_family, 
                               zero_ratio = round(zero_ratio, 3), 
                               observed_zeros = observed_zeros, 
                               expected_zeros = expected_zeros)
  zero_inflation_results <- rbind(zero_inflation_results, curr_zero_info)
}

write.csv(zero_inflation_results, 
          file = file.path(dir_path, "zero_inflation_results_best_models.csv"), 
          row.names = FALSE)

## Predictive performance: overall (posterior predictive distribution) ---------
for (species in species_group) {
  best_info <- readRDS(file.path(dir_path, paste0(species, 
                                                  "_final_model_info.rds")))
  curr_species <- best_info$species
  curr_formula <- best_info$formula_number
  curr_family <- best_info$family
  mieren_sp <- best_info$mieren_data
  
  best_model <- readRDS(file.path(dir_path, paste0(species, "_", 
                                                   best_model_name, ".rds")))
  cpo <- best_model$cpo$cpo
  pit <- best_model$cpo$pit
  cpo_failure <- best_model$cpo$failure
  
  diagnostics_df <- data.frame(species = curr_species,
                               observation_id = 1:length(cpo), cpo = cpo,
                               pit = pit, cpo_failure = cpo_failure)
  
  write.csv(diagnostics_df, 
            file = file.path(dir_path, paste0(curr_species, 
                                              "_diagnostics_data.csv")),
            row.names = FALSE)
  
  ## Vizualization
  cpo_plot <- ggplot(diagnostics_df, aes(x = cpo)) +
    geom_histogram(bins = 30, fill = "gray", alpha = 0.7) +
    labs(title = paste("CPO Distribution -", curr_species),
         x = "CPO Values", y = "Frequency") +
    theme_minimal()
  
  pit_plot <- ggplot(diagnostics_df, aes(x = pit)) +
    geom_histogram(bins = 20, fill = "gray", alpha = 0.7) +
    geom_hline(yintercept = nrow(diagnostics_df$pit)/20, 
               linetype = "dashed", color = "red") +
    labs(x = "PIT Values", y = "Frequency") + theme_minimal()
  
  # Save plots
  ggsave(filename = file.path(dir_path, paste0(curr_species, 
                                               "_cpo_histogram.png")), 
         plot = cpo_plot, width = width_inches, height = height_inches)
  
  ggsave(filename = file.path(dir_path, paste0(curr_species, 
                                               "_pit_histogram.png")), 
         plot = pit_plot, width = width_inches, height = height_inches)
}

## Predictive performance: overall (prediction accuracy) -----------------------
for (species in species_group) {
  best_info <- readRDS(file.path(dir_path, paste0(species, 
                                                  "_final_model_info.rds")))
  curr_species <- best_info$species
  curr_formula <- best_info$formula_number
  curr_family <- best_info$family
  mieren_sp <- best_info$mieren_data
  
  best_model <- readRDS(file.path(dir_path, paste0(species, "_", 
                                                   best_model_name, ".rds")))
  best_stack <- best_info$best_stack
  stack_indices <- inla.stack.index(best_stack, "est")
  fitted_values <- best_model$summary.fitted.values$mean[stack_indices$data]
  fitted_lower <- best_model$summary.fitted.values$
    `0.025quant`[stack_indices$data]
  fitted_upper <- best_model$summary.fitted.values$
    `0.975quant`[stack_indices$data]
  observed <- mieren_sp$is_present
  mape_probability <- mean(abs(observed - fitted_values)) * 100
  
  
  mieren_sp$fitted_mean <- fitted_values
  mieren_sp$fitted_lower <- fitted_lower
  mieren_sp$fitted_upper <- fitted_upper
  mieren_sp$residual <- mieren_sp$is_present - mieren_sp$fitted_mean
  mieren_sp$abs_residual <- abs(mieren_sp$residual)
  
  saveRDS(mieren_sp, file = file.path(dir_path, 
                                      paste0(curr_species, 
                                             "_data_abs_residuals.rds")))
}
# statistical analysis: predictions --------------------------------------------
## General setup: load prediction grid -----------------------------------------
pred_grid_path <- file.path(dir_path, "pred_grid_list.rds")
if (file.exists(pred_grid_path)) {
  pred_grid <- readRDS(pred_grid_path)
  grid_habitats_clean <- pred_grid$grid_habitats_clean
  n_cells_clean <- pred_grid$n_cells_clean
} else {
  bb <- st_bbox(limburg_utm)
  grid <- st_make_grid(
    limburg_utm,
    cellsize = resolution_grid,
    what = "polygons",
    square = TRUE) |> st_as_sf() |> st_intersection(limburg_utm)
  
  grid_predictions <- grid %>% mutate(cell_id = 1:n()) %>% 
    st_set_crs(st_crs(limburg_utm))
  
  # Process dominant habitats
  grid_predictions <- st_transform(grid_predictions, st_crs(limburg_utm))
  BWK_map <- st_transform(BWK_map, st_crs(limburg_utm))
  
  dominant_habitats <- grid_predictions %>% st_intersection(BWK_map) %>%
    mutate(habitat_type_dutch = sapply(BWKLABEL, function(x) classify_habitat(x)[1]),
           habitat_type = group_and_translate_habitat(habitat_type_dutch),
           fragment_area = st_area(.)) %>%  st_drop_geometry() %>%
    group_by(cell_id, habitat_type) %>% summarise(area = sum(fragment_area), 
                                                  .groups = "drop_last") %>%
    mutate(percentage = as.numeric(area / sum(area) * 100)) %>% ungroup() %>%
    group_by(cell_id) %>% slice_max(percentage, n = 1, with_ties = FALSE) %>%
    dplyr::select(cell_id, dominant_habitat = habitat_type, 
                  dominant_percentage = percentage)
  
  grid_habitats <- grid_predictions %>%
    left_join(dominant_habitats, by = "cell_id")
  na_cells <- is.na(grid_habitats$dominant_habitat)
  grid_habitats_clean <- grid_habitats[!na_cells, ]
  
  ## Filter for habitat types that are observed
  unique_obs_habitats <- levels(mieren_data$habitat_type)
  valid_habitat <- grid_habitats_clean$dominant_habitat %in% 
    unique_obs_habitats
  grid_habitats_clean <- grid_habitats_clean[valid_habitat, ]
  n_cells_clean <- nrow(grid_habitats_clean)
  grid_habitats_clean$
    dominant_habitat <- factor(grid_habitats_clean$dominant_habitat, 
                               levels = levels(mieren_data$habitat_type))
}
## Generate prediction grid map (same for every species) -----------------------
plot_data <- grid_habitats_clean %>% 
  dplyr::select(dominant_habitat, dominant_percentage, cell_id)

# Ensure dominant_habitat is a factor
plot_data$dominant_habitat <- as.factor(plot_data$dominant_habitat)

# Create color palette for habitats
habitat_colors <- scales::hue_pal()(length(levels(plot_data$dominant_habitat)))
names(habitat_colors) <- levels(plot_data$dominant_habitat)

fig4 <- ggplot() + 
  geom_sf(data = limburg_utm, fill = "white", color = "black") +
  geom_sf(data = plot_data, aes(fill = dominant_habitat), 
          color = "white", size = 0.1, alpha = 0.8) + 
  scale_fill_manual(values = habitat_colors, name = "Habitat type", 
                    drop = FALSE) + theme_minimal() + 
  labs(title = "Prediction grid by habitat type", x = "", y = "") +
  theme(legend.position = "right", 
        plot.title = element_text(hjust = 0.5),
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank())

ggsave(filename = file.path(dir_path, "prediction_grid_habitats.png"), 
       plot = fig4, width = width_inches, height = height_inches, dpi = 300)
## Generate predictions for each species ---------------------------------------
pred_grid_path <- file.path(dir_path, "pred_grid_list.rds")
pred_grid <- readRDS(pred_grid_path)
grid_habitats_clean <- pred_grid$grid_habitats_clean
all_years <- sort(unique(mieren$jaar1))  
selected_years <- seq(min(all_years), max(all_years), by = year_jumps)
if (!2024 %in% selected_years){
  selected_years <- c(selected_years, 2024)
}
predictions_list <- list()
for (curr_species in species_group) {
  best_info <- readRDS(file.path(dir_path, paste0(curr_species, 
                                                  "_final_model_info.rds")))
  best_model <- readRDS(file.path(dir_path, paste0(curr_species, "_", 
                                                   best_model_name, ".rds")))
  
  index_pred <- inla.stack.index(best_info$best_stack, "pred")$data
  fitted_means <- best_model$summary.fitted.values[index_pred, "mean"]
  
  pred_years <- sort(unique(best_info$mieren_data$jaar1))
  n_years <- length(pred_years)
  n_cells_clean <- pred_grid$n_cells_clean
  unique_cell_ids <- grid_habitats_clean$cell_id
  
  pred_df <- data.frame(
    cell_id = rep(unique_cell_ids, n_years),
    year = rep(pred_years, each = n_cells_clean),
    species = curr_species,
    model_used = best_model_name,
    pred_prob = plogis(fitted_means),
    pred_logit = fitted_means
  )
  
  predictions_list[[curr_species]] <- pred_df
  
  plot_df <- grid_habitats_clean %>%
    left_join(pred_df, by = "cell_id") %>%
    filter(year %in% selected_years)
  
  if (nrow(plot_df) > 0) {
    range_pred <- range(plot_df$pred_prob, na.rm = TRUE)
    
    fig <- ggplot() + 
      geom_sf(data = limburg_utm, fill = NA) +
      geom_sf(data = plot_df, aes(fill = pred_prob), color = NA) +
      scale_fill_viridis_c(
        name = "Presence Probability", 
        option = "H", 
        trans = "sqrt", 
        limits = c(0, max(range_pred)),
        breaks = seq(0, max(range_pred), length.out = 5), 
        guide = guide_colorbar(barheight = unit(2, "cm"),
                               barwidth = unit(0.5, "cm"))) + 
      theme_minimal() + 
      facet_wrap(~year) + 
      ggtitle(paste("Presence probability -", curr_species)) + 
      theme(axis.text = element_blank(), axis.ticks = element_blank(), 
            axis.title = element_blank())
    
    ggsave(filename = file.path(dir_path, paste0("pred_pres_plot_", 
                                                 curr_species, ".png")), 
           plot = fig, width = width_inches, height = height_inches, dpi = 300)
  }
  
  print(paste("Completed predictions and visualization for ", curr_species))
}

# Save all predictions at the end
saveRDS(predictions_list, file = file.path(dir_path, "predictions_list.rds"))
all_preds <- do.call(rbind, predictions_list)
saveRDS(all_preds, file = file.path(dir_path, "all_preds.rds"))



## summary statistics of the habitats of the grid (using BWK_map) --------------
summary_hab <- table(grid_habitats_clean$dominant_habitat)
write.csv(summary_hab, 
          file = file.path(dir_path, "summary_habitate_predictions.csv"), 
          row.names = TRUE)