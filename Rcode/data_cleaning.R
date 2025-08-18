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
library(data.table)

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

############################# Loading the data #################################
# Mieren data
mieren_original <- read_excel("Ongewervelden-allegegevens.xlsx")
mieren_original <- fread("Ongewervelden-allegegevens.csv")

# Methode van vangst dataset
mieren_methode <- read_excel("Hmsc/Mieren-methode.xlsx")

# BWK volgens cristina (BWK = biologische waarderingskaart voor bodemtype)
BWK_class <- read_excel("BWK-Cristina2.xlsx") 

# BWK map limburg
BWK_map <- st_read("BWK_map2/Shapefile/BwkHab70000.shp")

######################### Modifications to the data ############################
# Filter on important variables ------------------------------------------------

mieren <- mieren %>% dplyr::select(Soortgroep, Familie, `LATIJNSE NAAM`, 
                                   maand1, jaar1, dag1, maand2, jaar2, dag2, 
                                   Aantal, Mannetjes, Vrouwtjes, Gemeente, 
                                   Plaats, Biotoop, Plaatscode, Methode, 
                                   `Rode Lijst Vlaanderen`, Code)

# Filter on dates that occur after 1994 ---------------------------------------
mieren <- mieren %>% filter(mieren$jaar1 >= 1995)

# Filter on periods of trap setting  -------------------------------------------
mieren$date1 <- make_date(mieren$jaar1, mieren$maand1, mieren$dag1)
mieren$date2 <- make_date(mieren$jaar2, mieren$maand2, mieren$dag2)
mieren$diff_days <- as.numeric(difftime(mieren$date2, mieren$date1,
                                        units = "days"))
mieren <- mieren %>% filter(diff_days >= 10 & diff_days <= 30)

# Coordinates of the different trap locations ----------------------------------
map.shp <- st_read("mappen/punten.shp")
map.shp <- map.shp[order(map.shp$Plaatscode),]
pos_loc <- unique(mieren$Plaatscode)
map.shp <- map.shp[map.shp$Plaatscode %in% pos_loc,]
geometry <- rep(NA, nrow(mieren))
X <- rep(NA, nrow(mieren))
Y <- rep(NA, nrow(mieren))
for (i in 1:nrow(mieren)) {
  if (!is.na(mieren$Plaatscode[i])) {
    for (j in 1:nrow(map.shp)) {
      if (!is.na(map.shp$Plaatscode[j]) && 
          mieren$Plaatscode[i] == map.shp$Plaatscode[j]) {
        geometry[i] <- map.shp$geometry[j]
        X[i] <- map.shp$X[j]
        Y[i] <- map.shp$Y[j]
        break
      }
    }
  }
}

mieren$X <- X
mieren$Y <- Y
mieren <- mieren %>% filter(!is.na(mieren$X) | !is.na(mieren$Y)) %>% 
  filter(X != 0 & Y != 0)


# Presence absence of all ant species construction -----------------------------
mieren <- mieren %>% filter(!(is.na(Aantal)))
is_present_bin <- numeric(nrow(mieren))
for (i in 1:nrow(mieren)){
  if (mieren$Soortgroep[i] == 14){
    if (mieren$Aantal[i] > 0){
      is_present_bin[i] <- 1
    }
  }
}
mieren$is_present_bin <- is_present_bin
mieren$presence <- factor(ifelse(mieren$is_present_bin == 1, "present",
                                 "absent"), level= c("absent", "present"))
# Habitat type construction ----------------------------------------------------
## Extract most dominant habitat type for pixels 
mieren <- mieren %>% filter(!is.na(Code) & Code != "Overige")
mieren$habitat_type <- classify_habitat(mieren$Code)

# Complete loop that does both grouping and English translation
habitat_type_grouped <- mieren$habitat_type
for (i in 1:nrow(mieren)){
  curr_hab <- mieren$habitat_type[i]
  
  # Group and translate to English
  if (curr_hab == "Woonomgeving" | curr_hab == "Landbouw"){
    habitat_type_grouped[i] <- "Human influence"
  } else if (curr_hab == "Loofbos" | curr_hab == "Naaldbos"){
    habitat_type_grouped[i] <- "Forest"
  } else if (curr_hab == "Grasland"){
    habitat_type_grouped[i] <- "Grassland"
  } else if (curr_hab == "Droge heide"){
    habitat_type_grouped[i] <- "Dry heathland"
  } else if (curr_hab == "Natte heide"){
    habitat_type_grouped[i] <- "Wet heathland"
  } else if (curr_hab == "Overige heide"){
    habitat_type_grouped[i] <- "Other heathland"
  } else if (curr_hab == "Hoogstamboomgaard"){
    habitat_type_grouped[i] <- "High-stem orchard"
  } else if (curr_hab == "Overige"){
    habitat_type_grouped[i] <- "Other"
  }
  # All other habitat types remain unchanged
}

mieren$habitat_type <- habitat_type_grouped
mieren$habitat_type <- as.factor(mieren$habitat_type)
# Seasonality construction -----------------------------------------------------
is_winter_bin <- numeric(nrow(mieren))
for (i in 1:nrow(mieren)){
  if (mieren$maand1[i] >= 12 | mieren$maand1[i] < 2){
    is_winter_bin[i] <- 1
  } 
  if (mieren$maand1[i] == 2){
    if ((mieren$dag1[i] + mieren$diff_days[i]) <= 29){
      is_winter_bin[i] <- 1
    }
  }
}
mieren$is_winter_bin <- is_winter_bin
mieren$is_winter <- factor(ifelse(mieren$is_winter_bin == 1, "winter",
                                  "no winter"), level= c("no winter", 
                                                         "winter"))
# Removing habitat types that are not properly defined -------------------------
## if for a location 2 habitats are observed, they are excluded from the study
hab_types_loc <- mieren %>%
  group_by(Plaatscode, habitat_type) %>% summarise(first_year = min(jaar1), 
                                                   .groups = "drop") %>%
  group_by(Plaatscode) %>%
  summarise(
    number_of_habitats = n_distinct(habitat_type),
    unique_habitats = paste(unique(habitat_type), collapse = ", "),
    habitat_first_years = paste(paste0(habitat_type, ": ", first_year), 
                                collapse = "; "))
multi_hab_loc <- hab_types_loc %>% filter(number_of_habitats > 1)
write.csv(multi_hab_loc, "locations_with_multiple_habitat_types.csv")
single_hab_loc <- hab_types_loc %>% filter(number_of_habitats == 1)
mieren <- mieren %>% filter(Plaatscode %in% single_hab_loc$Plaatscode)

# Grouping of clusters into 1 observation --------------------------------------
## Check if there Aantallen are the same for the traps placed close to each
## other (was done to minimize mearusemetn error)
mieren_only <- mieren %>% filter(Soortgroep == 14)
mieren_counts <- mieren_only %>% group_by(Plaatscode, date1, date2) %>%
  summarise(group_size = n(), n_unique_counts = n_distinct(Aantal),
    percent_unique = (n_distinct(Aantal) / n()) * 100, .groups = "drop")
percent_traps_unique <- mean(mieren_counts$percent_unique == 100) * 100

global_counts <- mieren %>% group_by(Plaatscode, date1, date2) %>%
  summarise(group_size = n(), n_unique_counts = n_distinct(Aantal),
            percent_unique = (n_distinct(Aantal) / n()) * 100, .groups = "drop")
percent_traps_unique_global <- mean(global_counts$percent_unique == 100) * 100

## Combine mierenset (no small scale variation present)
mieren_clustered <- mieren %>% 
  group_by(Plaatscode, date1, date2, `LATIJNSE NAAM`) %>%
  mutate(is_present_bin = max(is_present_bin, na.rm = TRUE)) %>% ungroup() %>%
  distinct(Plaatscode, date1, date2, `LATIJNSE NAAM`, .keep_all = TRUE)


# Filter on trap type ----------------------------------------------------------
traps_of_interest <- c("bodemval", "Bodemval") 
mieren_all_traps <- mieren_clustered
mieren <- mieren_clustered %>% filter(Methode %in% traps_of_interest)
# Saving the modified dataset --------------------------------------------------
write.csv(mieren_clustered, "mieren_cluster_only_pitfalls.csv", 
          row.names = FALSE)
write.csv(mieren_all_traps, "mieren_clustered.csv", 
          row.names = FALSE)

