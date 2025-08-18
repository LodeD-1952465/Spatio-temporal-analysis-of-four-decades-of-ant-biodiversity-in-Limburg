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
mieren_traps <- read.csv("mieren_clustered.csv") 

# BWK map limburg
BWK_map <- st_read("BWK_map2/Shapefile/BwkHab70000.shp")

# General parameters -----------------------------------------------------------
specific_hab <- "" # leave empty when you want no specific habitat type
species_name <- ""  # leave empty when you want the most present species
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
best_model_name <- "formula20_binomial" # Name of the model with lowest WAIC/LMPL. Default is "". 
# Set up: data for specific species (only to study 1 specific species) ---------
## Select most observed species if no specific species or habitat was given
species_highest <- species_name
if (specific_hab == ""){
  if (species_name == ""){
    all_species <- mieren %>% group_by(`LATIJNSE NAAM`) %>%
      summarise(total_obs = sum(is_present_bin, na.rm = TRUE),
                total_years = length(unique(jaar2)),
                years_observed = list(unique(jaar2))) %>% 
      arrange(desc(total_obs)) %>% dplyr::select(`LATIJNSE NAAM`, total_obs, 
                                                 years_observed, total_years)
    species_highest <- all_species$`LATIJNSE NAAM`[1]
    years_observed <- all_species$years_observed[1]
  }
}

## Select most observed species if specific habitat but no specific species is 
## given
if (!(specific_hab == "")){
  mieren_specific_hab <- mieren %>% filter(habitat_type == specific_hab)
  all_species_spec_hab <- mieren_specific_hab %>% group_by(`LATIJNSE NAAM`) %>%
    summarise(total_obs = sum(is_present_bin, na.rm = TRUE),
              total_years = length(unique(jaar2)),
              years_observed = list(unique(jaar2))) %>% 
    arrange(desc(total_obs)) %>% dplyr::select(`LATIJNSE NAAM`, total_obs, 
                                               years_observed, total_years)
  species_highest <- all_species_spec_hab$`LATIJNSE NAAM`[1]
  years_observed <- all_species_spec_hab$years_observed[1]
}
## Select most observed species data when species name is given
if (!(species_name == "")){
  spec_species <- mieren %>% filter(`LATIJNSE NAAM` == species_name)
  spec_species_data <- spec_species %>%
    summarise(total_obs = sum(is_present_bin, na.rm = TRUE),
              total_years = length(unique(jaar2)),
              years_observed = list(unique(jaar2)))
  years_observed <- spec_species_data$years_observed
}

## Select most observed species data for given level of indangerd species
if (!(indangerd_species == "")){
  indangerd_mieren <- mieren %>% 
    filter(`Rode Lijst Vlaanderen` == indangerd_species )
  indangerd_data <- indangerd_mieren %>%  group_by(`LATIJNSE NAAM`) %>%
    summarise(total_obs = sum(is_present_bin, na.rm = TRUE),
              total_years = length(unique(jaar2)),
              years_observed = list(unique(jaar2))) %>%
    dplyr::select(`LATIJNSE NAAM`, total_obs, years_observed, total_years)
  species_highest <- indangerd_data$`LATIJNSE NAAM`[1]
  years_observed <- indangerd_data$years_observed[1]
}
mieren_data <- mieren
mier_data <- mieren %>% filter(`LATIJNSE NAAM` == species_highest)

## filter on habitat type that occurs for specified species and year 
all_habitats <- unique(mier_data$habitat_type[!is.na(mieren$habitat_type)])
table_habitat <- table(mier_data$habitat_type)
table_present <- table(mieren_data$is_present)
table_present
## filter on habitats that occur more than 10 times
table_habitat_big <- table_habitat[table_habitat > 10]
selected_habitats <- names(table_habitat_big)

## summary statistics regarding presence and habitat type
table(mier_data$habitat_type)
table(mier_data$is_present)

## create presence/absence for highest year for species specified above
mieren_data <- mieren_data %>% filter(!is.na(Aantal))
is_present_bin <- numeric(nrow(mieren_data))
for (i in 1:nrow(mieren_data)){
  if (mieren_data$`LATIJNSE NAAM`[i] == species_highest){
    if (mieren_data$Aantal[i] > 0){
      is_present_bin[i] <- 1
    }
  }
}

mieren_data <- mieren_data %>% dplyr::select(-is_present_bin, -presence)
mieren_data$is_present_bin <- is_present_bin
mieren_data$presence <- factor(ifelse(mieren_data$is_present_bin == 1, "present",
                                      "absent"), level= c("absent", "present"))
mieren_data <- mieren_data %>% filter(!is.na(habitat_type))
mieren_data$habitat_type <- factor(mieren_data$habitat_type, 
                                   levels = all_habitats)

mieren_data$is_winter <- factor(ifelse(mieren_data$is_winter_bin == 1, "winter",
                                       "no winter"), level= c("no winter", 
                                                              "winter"))
mieren_data$jaar1_num <- as.numeric(mieren_data$jaar1)
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

########################## Data Exploration: Global ############################
# Create folder to store data exploration results ------------------------------
dir_path_exploration <- "data exploration results"

if (!dir.exists(dir_path_exploration)) {
  dir.create(dir_path_exploration)
  cat("Folder created successfully at:", dir_path_exploration, "\n")
} else {
  cat("Folder already exists at:", dir_path_exploration, "\n")
} 

## Transform the mieren coordinate system to the Limburg province coordinates
mieren_data <- mieren
mieren_sf <- st_as_sf(mieren_data, coords = c("X", "Y"), crs = 31370)
mieren_sf$X <- st_coordinates(mieren_sf)[,1]
mieren_sf$Y <- st_coordinates(mieren_sf)[,2]

### Obtain the Limburg province 
search_eurostat("nuts")
nuts_sf <- get_eurostat_geospatial(output_class = "sf", resolution = "01", 
                                   nuts_level = 2)
belgium_nuts <- nuts_sf[nuts_sf$CNTR_CODE == "BE", ]
limburg <- belgium_nuts[belgium_nuts$NUTS_NAME == "Prov. Limburg (BE)", ] 
limburg_utm <- st_transform(limburg, crs = 31370)

# Proportion of soortgroepen data exploration ----------------------------------
species_counts <- table(mieren_data$Soortgroep)
species_proportions <- prop.table(species_counts)
species_df <- data.frame(
  Soortgroep = names(species_proportions),
  Count = as.numeric(species_counts),
  Proportion = as.numeric(species_proportions))
write.csv(species_df, file.path(dir_path_exploration, 
                                "prop_soorgroep.csv"))

# Number of different unique species per location data exploration ------------- 
loc_species <- mieren_data %>%  group_by(`Plaatscode`) %>%
  summarise(number_species = n_distinct(`LATIJNSE NAAM`), X = first(X),
            Y = first(Y),.groups = "drop")
saveRDS(loc_species, 
        file = file.path(dir_path_exploration,
                         "amount_species_per_loc.rds"))
# Total amount of observation locations ----------------------------------------
tot_obs_loc <- mieren_data %>% summarise(n_count = n_distinct(Plaatscode)) 
tot_obs_loc <- tot_obs_loc$n_count[1]

# Catching method data exploration ---------------------------------------------
mieren_traps$Methode <- tolower(mieren_traps$Methode)
catch_types <- table(mieren_traps$Methode)
prop_catch_types <- prop.table(catch_types) * 100
prop_catch_types <- sort(prop_catch_types, decreasing = TRUE)
trap_percentages <- data.frame(trap_type = names(prop_catch_types),percentage = as.numeric(prop_catch_types))
write.csv(trap_percentages, file.path(dir_path_exploration, 
                                 "amount_catch_types.csv"))

# Yearly amount of traps set data exploration ----------------------------------
## overall amount of yearly traps
num_traps_yearly <- mieren_data %>% group_by(jaar1) %>% 
  summarise(total_traps = n_distinct(Plaatscode))
mieren_only <- mieren_data %>% filter(Soortgroep == 14)
ant_traps_yearly <- mieren_only %>% group_by(jaar1) %>% 
  summarise(ant_traps = n_distinct(Plaatscode))
all_traps_yearly <- num_traps_yearly %>%
  left_join(ant_traps_yearly, by = "jaar1", suffix = c("_overall", "_ants"))
colnames(all_traps_yearly)[colnames(all_traps_yearly) == "jaar1"] <- "Jaar"
all_traps_yearly <- all_traps_yearly %>% 
  mutate( ant_traps = ifelse(is.na(ant_traps), 0, ant_traps))
write.csv(all_traps_yearly, file.path(dir_path_exploration, 
                                      "amount_traps_yearly.csv"))
plot_traps <- all_traps_yearly %>%
  pivot_longer(cols = c(total_traps, ant_traps), 
               names_to = "Type", values_to = "Aantal")

fig0 <- ggplot(plot_traps, aes(x = Jaar, y = Aantal, fill = Type)) +
  geom_col(position = "dodge") +
  labs(title = "", x = "Year", 
       y = "Number of unique traps", fill = "Type") + 
  theme_minimal()

ggsave(filename = file.path(dir_path_exploration, 
                            "temp_dist_overall.png"), 
       plot = fig0, width = width_inches, height = height_inches, dpi = 300)

## Dry and wet heath specific yearly amount of traps
dry_heath_traps_yearly <- mieren_data %>%
  filter(habitat_type == "Droge heide") %>%
  group_by(jaar1) %>%
  summarise(dry_heath_traps = n_distinct(Plaatscode))

wet_heath_traps_yearly <- mieren_data %>%
  filter(habitat_type == "Natte heide") %>%
  group_by(jaar1) %>%
  summarise(wet_heath_traps = n_distinct(Plaatscode))
heath_traps_yearly <- dry_heath_traps_yearly %>%
  full_join(wet_heath_traps_yearly, by = "jaar1") %>%
  mutate(dry_heath_traps = ifelse(is.na(dry_heath_traps), 0, dry_heath_traps),
         wet_heath_traps = ifelse(is.na(wet_heath_traps), 0, wet_heath_traps))

plot_heath_traps <- heath_traps_yearly %>%
  pivot_longer(cols = c(dry_heath_traps, wet_heath_traps),
               names_to = "Type", values_to = "Aantal")

fig1 <- ggplot(plot_heath_traps, aes(x = jaar1, y = Aantal, fill = Type)) +
  geom_col(position = "dodge") + labs(title = "", x = "Year", 
                                      y = "Number of unique traps", 
                                      fill = "Habitat type") +
  theme_minimal()

ggsave(filename = file.path(dir_path_exploration, "temp_dist_heath.png"),
       plot = fig1, width = width_inches, height = height_inches, 
       dpi = 300)

# Amount of different habitat types data exploration ---------------------------
am_hab_type <- table(mieren_data$habitat_type)
prop_hab_type <- prop.table(am_hab_type)
per_hab_type <- prop_hab_type * 100
write.csv(per_hab_type, file.path(dir_path_exploration, 
                                  "percentage_habitat_type_observations.csv"))

# Percentage of different habitats prediction cells (overall) ------------------
BWK_map_class <- BWK_map %>% mutate(habitat_type = classify_habitat(BWK_map$
                                                                      BWKLABEL))
habitat_area_summary <- BWK_map_class %>% mutate(area = st_area(.)) %>%
  mutate(area = as.numeric(area)) %>% group_by(habitat_type) %>%  
  summarise(total_area = sum(area, na.rm = TRUE)) %>%
  mutate(percentage = (total_area / sum(total_area)) * 100,
         area_km2 = total_area / 1000000) %>% arrange(desc(percentage)) 
habitat_area_summary_table <- habitat_area_summary %>% 
  dplyr::select(habitat_type, total_area, percentage, area_km2) %>%
  sf::st_drop_geometry()
write.csv(habitat_area_summary_table, 
          file.path(dir_path_exploration, "percentage_habitat_type_all.csv"))

# Trap percentages per habitat observations (overall)---------------------------
trap_per_hab <- table(mieren_sf$habitat_type)
per_traps_hab <- prop.table(trap_per_hab)
write.csv(per_traps_hab, file.path(dir_path_exploration, 
                                   "percent_traps_per_habitat.csv"))

# Trap locations per habitat type ----------------------------------------------

habitat_colors <- scales::hue_pal()(length(unique(mieren_sf$habitat_type)))
names(habitat_colors) <- unique(mieren_sf$habitat_type)
fig2 <- ggplot() + geom_sf(data = limburg_utm, fill = "white", 
                           color = "black") +
  geom_point(data = mieren_sf, aes(x = X, y = Y, color = habitat_type), 
             size = 1, alpha = 0.7) + 
  scale_color_manual(values = habitat_colors, name = "Habitat Type", 
                     drop = FALSE) + theme_minimal() + 
  labs(title = "Trap locations by Habitat Type", x = "", y = "", 
       caption = paste("Total observations:", nrow(mieren_sf))) +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(filename = file.path(dir_path_exploration, 
                            "trap_locations_ants.png"), plot = fig2,
       width = width_inches, height = height_inches, dpi = 300)

# Temporal summary for all habitats of interest (ants only) --------------------
for (habitat in hab_of_interest) {
  mieren_habitat <- mieren_only %>% filter(habitat_type == habitat)
  all_species_habitat <- mieren_habitat %>% 
    group_by(`LATIJNSE NAAM`) %>%
    summarise(total_obs = sum(is_present_bin, na.rm = TRUE),
              total_years = length(unique(jaar2)),
              years_observed = toString(unique(jaar2))) %>% 
    arrange(desc(total_obs)) %>% 
    dplyr::select(`LATIJNSE NAAM`, total_obs, years_observed, total_years)
  write.csv(all_species_habitat, 
            file.path(dir_path_exploration, paste0("temporal_summary_",
                                                   habitat, ".csv")))
}

# Amount of trap locations winter/no winter (overall) --------------------------
overall_traps_winter <- table(mieren_data$is_winter)
overall_df <- data.frame(category = names(overall_traps_winter), 
                         overall_count = as.numeric(overall_traps_winter),
                         overall_prop = as.numeric(prop.table(
                           overall_traps_winter)))
ant_traps_winter <- table(mieren_only$is_winter)
ant_df <- data.frame(category = names(ant_traps_winter),
                     ant_count = as.numeric(ant_traps_winter),
                     ant_prop = as.numeric(prop.table(ant_traps_winter)))
all_winter_df <- merge(overall_df, ant_df, by = "category", all = TRUE)
write.csv(all_winter_df, file.path(dir_path_exploration, 
                                   "percent_traps_winter.csv"))
# Amount of trap locations winter/no winter (for each habitat) -----------------
for (habitat in hab_of_interest) {
  mieren_habitat <- mieren_data %>% filter(habitat_type == habitat)
  am_trap_habitat_winter <- table(mieren_habitat$is_winter)
  prop_trap_habitat_winter <- prop.table(am_trap_habitat_winter)
  write.csv(prop_trap_habitat_winter, file.path(dir_path_exploration,
                                                paste0("prop_trap_winter_",
                                                       habitat, ".csv")))
}
# Trap locations by winter/no winter (overall) ---------------------------------
selected_years <- seq(min(mieren_data$jaar1), max(mieren_data$jaar1), 
                      by = year_jumps)
selected_years <- c(selected_years, 2024)
mieren_sf_viz <- mieren_sf %>% filter(jaar1 %in% selected_years)
mieren_sf_viz <- mieren_sf_viz %>% arrange(is_winter)

fig3 <- ggplot() + geom_sf(data = limburg_utm, fill = "white", 
                           color = "black") +
  geom_point(data = mieren_sf_viz, aes(x = X, y = Y, color = is_winter), 
             size = 1, alpha = 0.7) + theme_minimal() +  
  facet_wrap(~jaar1) + 
  labs(title = "", x = "", y = "") +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(filename = file.path(dir_path_exploration, 
                            "is_winter_loc_overall.png"), 
       plot = fig3, width = width_inches, height = height_inches, dpi = 300)

# Trap locations by winter/no winter (for each habitat) ------------------------
mieren_sf_viz$jaar1_fac <- as.factor(mieren_sf_viz$jaar1)
for (habitat in hab_of_interest){
  mieren_sf_hab <- mieren_sf_viz %>% filter(habitat_type == habitat)
  mieren_sf_hab$jaar1_fac <- droplevels(mieren_sf_hab$jaar1_fac)
  if (nrow(mieren_sf_hab) == 0) {
    message(paste("No data for habitat:", habitat, "- skipping plot."))
    next
  }
  fig4 <- ggplot() + geom_sf(data = limburg_utm, fill = "white", 
                             color = "black") +
    geom_point(data = mieren_sf_hab, aes(x = X, y = Y, color = is_winter), 
               size = 1, alpha = 0.7) + theme_minimal() +  
    facet_wrap(~jaar1_fac) + 
    labs(title = "", x = "", y = "") +
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(filename = file.path(dir_path_exploration, paste0("is_winter_loc_",
                                                           habitat, ".png")), 
         plot = fig4, width = width_inches, height = height_inches, dpi = 300)
}

# Indangerd species histogram (overall) ----------------------------------------
mieren$`Rode Lijst Vlaanderen` <- as.factor(mieren_data$`Rode Lijst Vlaanderen`)
mieren_hist <- mieren_data %>% filter(!is.na(`Rode Lijst Vlaanderen`))
mieren_hist <- mieren_hist %>% filter(Soortgroep == 14)
selected_years <- seq(min(mieren_hist$jaar1), max(mieren_hist$jaar1), 
                      by = year_jumps) 
mieren_hist <- mieren_hist %>% filter(jaar1 %in% selected_years)
total_observations_per_year <- mieren_hist %>% group_by(jaar1) %>%
  summarise(total_observations = n(), .groups = "drop")
mieren_summary_hist <- mieren_hist %>% 
  group_by(jaar1, `Rode Lijst Vlaanderen`) %>%
  summarise(presence_count = sum(is_present_bin == 1),.groups = "drop") %>%
  left_join(total_observations_per_year, by = "jaar1") %>%
  mutate(presence_rate = presence_count / total_observations)
mieren_summary_hist <- mieren_summary_hist %>% 
  filter(!is.na(mieren_summary_hist$`Rode Lijst Vlaanderen`))

# Define threat categories in Dutch and their English translations
threat_order_dutch <- c("Met uitsterven bedreigd", "Sterk bedreigd", "Bedreigd", 
                        "Kwetsbaar", "Momenteel niet bedreigd", "Zeldzaam")
threat_order_english <- c("Critically Endangered", "Endangered", "Threatened", 
                          "Vulnerable", "Not Currently Threatened", "Rare")

# Create translation mapping
threat_translation <- setNames(threat_order_english, threat_order_dutch)

# Apply factor levels and translate
mieren_summary_hist$`Rode Lijst Vlaanderen` <- factor(
  mieren_summary_hist$`Rode Lijst Vlaanderen`, levels = threat_order_dutch)

# Create English version for plotting
mieren_summary_hist$RedList_English <- threat_translation[mieren_summary_hist$`Rode Lijst Vlaanderen`]
mieren_summary_hist$RedList_English <- factor(mieren_summary_hist$RedList_English, 
                                              levels = threat_order_english)

fig5 <- ggplot(mieren_summary_hist, 
               aes(x = RedList_English, 
                   y = presence_rate)) +
  geom_col(position = "dodge", width = 0.7) + 
  geom_text(aes(label = paste0(round(presence_rate, 1))), 
            position = position_dodge(width = 1),
            vjust = -0.5, size = 3.5) +  
  facet_wrap(~jaar1, ncol = 1, scales = "free_y") +  
  theme_bw() +
  labs(y = "Presence rate (%)", x = "Red List Flanders") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), 
        axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(), 
        strip.text = element_text(size = 12), 
        panel.spacing = unit(1, "lines")) +  
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, 1.5))

ggsave(filename = file.path(dir_path_exploration, 
                            "endangered_species_hist.png"), 
       plot = fig5, width = width_inches, height = height_inches + 1,
       dpi = 300)


#################### Data Exploration: species specific ########################
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


# Overview of number of observations for every selected species ----------------

num_obs_per_species <- mieren %>% filter(`LATIJNSE NAAM` %in% species_group) %>% 
  group_by(`LATIJNSE NAAM`) %>%
  summarise(n_observations = n(), .groups = "drop") %>% 
  rename(species_name = `LATIJNSE NAAM`) %>%
  arrange(desc(n_observations))

write.csv(num_obs_per_species, 
          file = file.path(dir_path, "species_group_observation_counts.csv"), 
          row.names = FALSE)

# Number of observations per covariates for every species analyzed -------------
for (curr_species in species_group) {
  # Data preparation -----------------------------------------------------------
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
  
  # Data exploration: number of observations for every habitate type -----------
  species_specific <- mieren_data %>% filter(`LATIJNSE NAAM` == curr_species)
  presence_counts <- species_specific %>%
    filter(is_present_bin == 1) %>%
    count(habitat_type, name = "presence_count")
  
  percentage_presence <- species_specific %>%
    group_by(habitat_type) %>%
    summarise(
      presence_count = n(),
      total_obs = sum(is_present == 1),
      presence_percentage = 100 * (presence_count / total_obs)
    )
  
  write.csv(percentage_presence, 
            file = file.path(dir_path, paste0(curr_species,
                                              "_presence_per_hab.csv")))
  rm(presence_counts, percentage_presence)
  gc()
  
  # Data exploration: number of observations for every winter/no winter --------
  species_specific <- mieren_data %>% filter(`LATIJNSE NAAM` == curr_species)
  obs_is_winter <- species_specific %>% group_by(is_winter_bin) %>%
    summarize(num_winter_obs = sum(is_winter_bin == 1),
              num_no_winter_obs = sum(is_winter_bin == 0),
              total_obs = n())
  table_is_winter <- table(species_specific$is_winter)
  prop_is_winter <- prop.table(table_is_winter)
    
  saveRDS(prop_is_winter, file = file.path(dir_path, 
                                           paste0(curr_species,
                                                  "prop_is_winter.rds")))
  write.csv(obs_is_winter, 
            file = file.path(dir_path, paste0(curr_species,
                                              "_is_winter_obs.csv")))
  rm(obs_is_winter)
  gc()
  print(paste0("data exploration of species ", curr_species, " finished."))
}
# Species specific data exploration for every species analysed -----------------
for (species in species_group){
  # Create species-specific map ------------------------------------------------
  short_species <- species
  if (length(strsplit(species, ",")[[1]]) > 1){
    short_species <- strsplit(species, ",")[[1]][1]
    if (length(strsplit(short_species, "\\(")[[1]]) > 1){
      short_species <- strsplit(short_species, "\\(")[[1]][1]
    }
  }
  if (length(strsplit(short_species, "\\(")[[1]]) > 1){
    short_species <- strsplit(species, "\\(")[[1]][1]
  }
  short_species <- trimws(short_species)
  clean_short_species <- ""
  for (i in 1:nchar(short_species)) {
    char <- substr(short_species, i, i)
    if (char == " ") {
      clean_short_species <- paste0(clean_short_species, "_")
    } else {
      clean_short_species <- paste0(clean_short_species, char)
    }
  }
  
  species_dir_name <- paste0("data exploration of ", short_species)
  species_dir_path <- file.path(dir_path, species_dir_name)
  if (!dir.exists(species_dir_path)) {
    dir.create(species_dir_path, recursive = TRUE)
  }
  
  # Data preparation -------------------------------------------------------------
  mieren_data <- mieren
  is_present <- numeric(nrow(mieren_data))
  for (i in 1:nrow(mieren_data)){
    if (mieren_data$`LATIJNSE NAAM`[i] == species){
      if (mieren_data$Aantal[i] > 0){
        is_present[i] <- 1
      }
    }
  }
  mieren_data$is_present <- is_present
  mieren_data$jaar1_num <- as.numeric(mieren_data$jaar1)
  mieren_data$year_group <- as.integer(factor(mieren_data$jaar1))
  mieren_data$habitat_type <- as.factor(mieren_data$habitat_type)
  selected_years <- seq(min(mieren_data$jaar1), max(mieren_data$jaar1), 
                        by = year_jumps)
  selected_years <- c(selected_years, 2024)
  
  write.csv(mieren_data, 
            file = file.path(species_dir_path, paste0("data_", short_species, 
                                                      ".csv")), 
            row.names = FALSE)
  
  # Presence/absence distribution per habitat type -----------------------------
  mieren_pa <- mieren_data %>% filter(!is.na(habitat_type))
  mieren_pa_rates <- mieren_pa %>% group_by(jaar1, habitat_type) %>%
    summarise(presence_count = sum(is_present == 1), total_observations = n(),
              presence_rate = round(presence_count / total_observations, 4),
              .groups = "drop")
  
  mieren_pa_viz <- mieren_pa_rates %>% filter(jaar1 %in% selected_years)
  fig11 <- ggplot(mieren_pa_viz, aes(x = habitat_type, y = presence_rate)) +
    geom_col(position = "dodge", width = 0.7) + 
    geom_text(aes(label = ifelse(presence_rate > 0, round(presence_rate, 2), 
                                 "")), position = position_dodge(width = 1), 
              vjust = -0.5, size = 2.5) +  
    facet_wrap(~jaar1, ncol = 1, scales = "free_y") + theme_bw() +
    labs(y = "Presence rate", x = "Habitat Type") +
    theme(axis.text.x = element_text(size = 8), 
          axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
          strip.text = element_text(size = 8), 
          panel.spacing = unit(1, "lines")) +  
    scale_y_continuous(expand = expansion(mult = c(0, 0.25)),  
                       limits = c(0, max(mieren_pa_rates$presence_rate)))
  
  ggsave(filename = file.path(species_dir_path, paste0("hab_type_dist_", 
                                                       species, ".png")), 
         plot = fig11, width = width_inches, height = height_inches, dpi = 300)
  
  # Presence/absence locations per habitat type --------------------------------
  mieren_sf_pa <- st_as_sf(mieren_pa, coords = c("X", "Y"), 
                           crs = st_crs(limburg_utm))
  mieren_sf_pa_viz <- mieren_sf_pa %>% filter(jaar1 %in% selected_years)
  mieren_sf_pa_viz <- mieren_sf_pa_viz %>%
    mutate(
      X_orig = st_coordinates(.)[,1],
      Y_orig = st_coordinates(.)[,2],
      X_jitter = X_orig + rnorm(n(), 0, 1000),
      Y_jitter = Y_orig + rnorm(n(), 0, 1000)
    ) %>%
    st_drop_geometry() %>%
    st_as_sf(coords = c("X_jitter", "Y_jitter"), crs = st_crs(limburg_utm))
  
  mieren_sf_pa_viz <- st_filter(mieren_sf_pa_viz, limburg_utm)
  mieren_sf_pa_viz <- mieren_sf_pa_viz %>%
    mutate(presence_factor = factor(is_present, levels = c("0", "1"))) %>%
    arrange(presence_factor)
  
  presence_colors <- scales::hue_pal()(2)
  names(presence_colors) <- c("0", "1")
  fig12 <- ggplot() + 
    geom_sf(data = limburg_utm, fill = "white", color = "black") +
    geom_sf(data = mieren_sf_pa_viz, aes(color = presence_factor), 
            size = 0.5, alpha = 0.7) + 
    scale_color_manual(values = presence_colors, 
                       name = "Presence", 
                       labels = c("Absent", "Present"), 
                       drop = FALSE) + facet_wrap(~jaar1) + theme_minimal() + 
    labs(title = "", x = "", y = "") +
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          legend.position = "right", 
          plot.title = element_text(hjust = 0.5), 
          panel.spacing = unit(0.5, "lines"))
  ggsave(filename = file.path(species_dir_path, 
                              paste0("presence_absence_locations_", 
                                     short_species, ".png")), 
         plot = fig12, width = width_inches, height = height_inches, dpi = 300)
  
  # Table with presence rates by year and habitat type -------------------------
  mieren_table_data <- mieren_data %>% 
    filter(!is.na(habitat_type))
  
  presence_rate_table <- mieren_table_data %>%
    group_by(jaar1, habitat_type) %>%
    summarise(presence_count = sum(is_present == 1, na.rm = TRUE),
              total_observations = n(),
              presence_rate = round(presence_count / total_observations, 5),
              .groups = "drop") %>% 
    dplyr::select(jaar1, habitat_type, presence_rate) %>%
    pivot_wider(names_from = habitat_type, 
                values_from = presence_rate, values_fill = 0) %>% arrange(jaar1)
  
  write.csv(presence_rate_table, 
            file = file.path(species_dir_path, 
                             paste0("presence_rate_by_habitat", 
                                    clean_short_species, ".csv")), 
            row.names = TRUE)
  
  # Table with all the presence absence ----------------------------------------
  amount_presence_absence <- table(mieren_sf_pa$is_present)
  write.csv(amount_presence_absence, 
            file = file.path(species_dir_path, paste0("amount_presence_absence",
                                                      short_species, ".csv")), 
            row.names = TRUE)
}

