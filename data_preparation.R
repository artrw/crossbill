library(tidyverse)
library(auk) # for reading eBird data
library(lubridate) # for extracting years from Date objects
library(sf) # for using shapefile data
library(tigris) # for state and county maps
library(tidycensus) # for county population data
library(geodist) # for fast point-to-point geospatial distance calculations
library(fixest) # for efficient fixed effect model estimation
library(stargazer) # for pretty regression output

# Arthur Wardle observer id 935644

# Census API key for package tidycensus
census_api_key("ADD API HERE")

# Location of eBird sampling data
sampling_filepath <- "ebd_sampling_relAug-2021/ebd_sampling_relAug-2021.txt"

# Forest Service "Ranger District Boundaries" shapefile, downloaded from
# https://data.fs.usda.gov/geodata/edw/datasets.php?dsetCategory=boundaries
rd4_sf <- st_read("S_USA.RangerDistrict/S_Usa.RangerDistrict.shp") %>%
  filter(REGION == "04") %>% # Region 4 only: NV, UT, parts of ID, WY
  st_make_valid() %>%
  st_transform(crs = st_crs(4326)) %>%
  select(FORESTNAME, DISTRICTOR, DISTRICTNA, GIS_ACRES, geometry) %>%
  rename(forest_name = FORESTNAME,
         district_num = DISTRICTOR,
         district_name = DISTRICTNA,
         acreage = GIS_ACRES) %>%
  mutate(color = case_when(district_name == "Minidoka Ranger District" ~ "red4",
                           TRUE ~ "slategrey"))

rd4_raw <- auk_sampling(sampling_filepath) %>%
  auk_bbox(rd4_sf) %>%
  auk_year(2008:2021) %>%
  auk_filter(file = tempfile(), overwrite = TRUE) %>%
  read_sampling(unique = FALSE)    

rd4_data <- rd4_raw %>%
  select(latitude, longitude, observation_date, observer_id) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = st_crs(4326)) %>%
  st_join(rd4_sf, join = st_within) %>%
  filter(!is.na(forest_name)) %>% # Remove obs not in a National Forest 
  # Split year at July 6 (date of Cassia Crossbill split)
  mutate(year = lubridate::year(observation_date),
         altyear = ifelse(lubridate::month(observation_date) > 7 |
                         (lubridate::month(observation_date) == 7 &
                          lubridate::mday(observation_date) >= 6),
                       lubridate::year(observation_date),
                       lubridate::year(observation_date) - 1)) %>%
  st_drop_geometry() 

annual_site_visits <- function(data){
  years <- sort(unique(data$year))
  sites <- sort(unique(data$district_name))
  output <- data.frame(district_name = sites)
  for (year in years){
    output <- left_join(output, 
                        count(as_tibble(data[data$year == year,]),
                              district_name,
                              name = paste0("year_", year)),
                        by = "district_name")
  }
  output[is.na(output)] <- 0
  return(output)
}

rd4_visits_data <- annual_site_visits(rd4_data) %>%
  pivot_longer(cols = starts_with("year_"), 
               names_to = "year", 
               names_prefix = "year_",
               values_to = "visits") %>%
  left_join(rd4_sf[,c("district_name", "district_num", "forest_name", "color")], 
            by = "district_name") %>%
  mutate(treatment = ifelse(district_name == "Minidoka Ranger District" &
                              year >= 2017, 1, 0),
         rubytreatment = ifelse(district_name == "Mountain City-Ruby Mountains-Jarbidge Ranger District" &
                                  year >= 2017, 1, 0),
         district_num = as.numeric(district_num),
         year = as.numeric(year))  

write_csv(rd4_visits_data, "rd4_visits_data.csv")

annual_site_observers <- function(data){
  data <- data %>%
    select(district_name, observer_id, year) %>%
    distinct()
  years <- sort(unique(data$year))
  sites <- sort(unique(data$district_name))
  output <- data.frame(district_name = sites)
  for (year in years){
    output <- left_join(output, 
                        count(as_tibble(data[data$year == year,]),
                              district_name,
                              name = paste0("year_", year)),
                        by = "district_name")
  }
  output[is.na(output)] <- 0
  return(output)
}

rd4_observers_data <- annual_site_observers(rd4_data) %>%
  pivot_longer(cols = starts_with("year_"), 
               names_to = "year", 
               names_prefix = "year_",
               values_to = "observers") %>%
  left_join(rd4_sf[,c("district_name", "district_num", "forest_name", "color")], 
            by = "district_name") %>%
  mutate(treatment = ifelse(district_name == "Minidoka Ranger District" &
                              year >= 2017, 1, 0),
         rubytreatment = ifelse(district_name == "Mountain City-Ruby Mountains-Jarbidge Ranger District" &
                                  year >= 2017, 1, 0),
         district_num = as.numeric(district_num),
         year = as.numeric(year)) 

write_csv(rd4_observers_data, "rd4_observers_data.csv")

minidoka_observers <- rd4_data %>%
  filter(district_name == "Minidoka Ranger District") %>%
  select(observer_id) %>%
  distinct() %>%
  mutate(observer_id = as.numeric(sub("obs", "", observer_id))) 

write_csv(minidoka_observers, "minidoka_observers.csv")

minidoka_cuts <- split(minidoka_observers$observer_id,
                       ceiling(seq_along(minidoka_observers$observer_id)/500))
print(paste(length(minidoka_cuts), "cuts to completion"))
minidoka_observers_raw <- data.frame()
for (cut in seq_along(minidoka_cuts)){
  auk_cut <- auk_sampling(sampling_filepath) %>%
    auk_observer(minidoka_cuts[[cut]]) %>%
    auk_filter(file = tempfile(), overwrite = TRUE) %>%
    read_sampling(unique = FALSE)  
  minidoka_observers_raw <- rbind(minidoka_observers_raw, auk_cut)
  print(paste("Cut", as.character(cut), "completed"))
}
rm(minidoka_cuts, auk_cut, cut)

minidoka_observers_data <- minidoka_observers_raw %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = st_crs(4326)) %>%
  st_join(rd4_sf[rd4_sf$district_name == "Minidoka Ranger District",], 
          join = st_within)

unique_observer_years_minidoka <- rd4_data %>%
  filter(district_name == "Minidoka Ranger District") %>%
  select(observer_id, year) %>%
  distinct() 

origin <- function(data, observer, year, district, buffer = 50){
  # Warning: Data MUST be arranged by observation date
  data <- data[data$observer_id == observer,] 
  first <- which(lubridate::year(data$observation_date) == year & 
                   data$district_name == district)[1]
  if (first < 10){
    return("Insufficient Prior Records")
  }
  county_counts <- table(data[max(1, (first-buffer)):(first-1),]$county_state)
  if (max(county_counts) < 10){
    return("No Good Match")
  } else
    return(names(sort(county_counts, decreasing = TRUE)[1]))
}

unique_observer_years_minidoka$origin <- map2_chr(unique_observer_years_minidoka$observer_id,
                                                  unique_observer_years_minidoka$year,
                                                  origin,
                                                  data = minidoka_observers_data,
                                                  district = "Minidoka Ranger District")

observers <- rd4_data %>%
  select(observer_id) %>%
  distinct() %>%
  mutate(observer_id = as.numeric(sub("obs", "", observer_id))) 

observer_cuts <- split(observers$observer_id,
                       ceiling(seq_along(observers$observer_id)/500))
print(paste(length(observer_cuts), "cuts to completion"))
observers_raw <- data.frame()
for (cut in seq_along(observer_cuts)){
  auk_cut <- auk_sampling(sampling_filepath) %>%
    auk_observer(observer_cuts[[cut]]) %>%
    auk_filter(file = tempfile(), overwrite = TRUE) %>%
    read_sampling(unique = FALSE)  
  observers_raw <- rbind(observers_raw, auk_cut)
  print(paste("Cut", as.character(cut), "completed"))
}
rm(observer_cuts, auk_cut, cut)

observers_data <- observers_raw %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = st_crs(4326)) %>%
  st_join(rd4_sf, join = st_within) %>%
  arrange(observation_date) %>%
  mutate(county_state = paste0(county, ", ", state))

unique_observer_years <- rd4_data %>%
  filter(!is.na(district_name)) %>%
  select(observer_id, year, district_name) %>%
  distinct() 

unique_observer_years$origin <- pmap_chr(list(unique_observer_years$observer_id,
                                              unique_observer_years$year,
                                              unique_observer_years$district_name),
                                         origin,
                                         data = observers_data)

write_csv(unique_observer_years, "unique_observer_years_origins.csv")

unique_observer_years <- read_csv("unique_observer_years_origins.csv")

us_sf <- tigris::states()
  
us_county_sf <- tigris::counties() %>%
  st_transform(4326) %>%
  rename(COUNTYNAME = NAME) %>%
  left_join(st_drop_geometry(us_sf)[,c("STATEFP", "NAME")], by = "STATEFP") %>%
  rename(STATENAME = NAME) %>%
  mutate(NAME = paste0(COUNTYNAME, ", ", STATENAME),
         FIP = paste0(STATEFP, COUNTYFP)) %>%
  {.[.$FIP == "29510", "NAME"] <- "St. Louis City, Missouri"; .} %>%
  {.[.$FIP == "51770", "NAME"] <- "Roanoke City, Virginia"; .} %>%
  {.[.$FIP == "51760", "NAME"] <- "Richmond City, Virginia"; .} %>%
  {.[.$FIP == "51600", "NAME"] <- "Fairfax City, Virginia"; .} %>%
  {.[.$FIP == "51515", "NAME"] <- "Bedford City, Virginia"; .} %>%
  {.[.$FIP == "51620", "NAME"] <- "Franklin City, Virginia"; .} %>%
  {.[.$FIP == "24510", "NAME"] <- "Baltimore City, Maryland"; .} %>%
  {.[.$FIP == "35013", "NAME"] <- "Dona Ana, New Mexico"; .} %>%
  {.[.$FIP == "17099", "NAME"] <- "La Salle, Illinois"; .} %>%
  {.[.$FIP == "02105", "NAME"] <- "Skagway-Hoonah-Angoon, Alaska"; .} %>% # Judgement call, boundary change
  {.[.$FIP == "02198", "NAME"] <- "Prince of Wales-Outer Ketchikan, Alaska"; .} %>% # Judgement call, boundary change
  select(NAME, geometry) %>%
  rename(county = NAME)
  
# Canadian Census District shapefile from 
# https://www12.statcan.gc.ca/census-recensement/2011/geo/bound-limit/bound-limit-2016-eng.cfm

canada_cd_sf <- st_read("lcd_000b16a_e/lcd_000b16a_e.shp") %>%
  st_transform(4326) %>%
  filter(!(PRNAME %in% c("Northwest Territories / Territoires du Nord-Ouest",
                         "Nunavut", "Yukon"))) %>%
  {.[.$PRNAME == "British Columbia / Colombie-Britannique", "PRNAME"] <- 
    "British Columbia"; .} %>%
  {.[.$PRNAME == "Newfoundland and Labrador / Terre-Neuve-et-Labrador", "PRNAME"] <- 
    "Newfoundland and Labrador"; .} %>%
  {.[.$PRNAME == "Quebec / Québec", "PRNAME"] <- 
    "Quebec"; .} %>%
  {.[.$PRNAME == "Prince Edward Island / Île-du-Prince-Édouard", "PRNAME"] <- 
    "Prince Edward Island"; .} %>%
  {.[.$PRNAME == "Nova Scotia / Nouvelle-Écosse", "PRNAME"] <- 
    "Nova Scotia"; .} %>%
  {.[.$PRNAME == "New Brunswick / Nouveau-Brunswick", "PRNAME"] <- 
    "New Brunswick"; .} %>%
  rename(county = CDNAME, state = PRNAME) %>%
  mutate(county = paste0(county, ", ", state)) %>%
  select(county, geometry)

county_sf <- rbind(us_county_sf, canada_cd_sf)
rm(us_county_sf, canada_cd_sf)

unique_observer_years_filtered <- unique_observer_years %>%
  mutate(origin = ifelse(origin == "NA, Guam", "Guam, Guam", origin)) %>%
  mutate(origin = ifelse(origin == "NA, American Samoa", 
                         "Western, American Samoa", 
                         origin)) %>% # Judgement call, highest population Samoa county
  mutate(origin_state = sub('.*,\\s*', '', origin)) %>%
  mutate(origin_state = ifelse(origin_state %in% c(us_sf$NAME, "Alberta",
                                                   "British Columbia",
                                                   "Manitoba", "New Brunswick",
                                                   "Newfoundland and Labrador",
                                                   "Nova Scotia", "Ontario",
                                                   "Prince Edward Island",
                                                   "Quebec", "Saskatchewan"),
                                                   origin_state, 
                               ifelse(origin_state %in% c("Insufficient Prior Records",
                                                          "No Good Match"), 
                                      origin_state, "Foreign"))) %>%
  filter(!(origin_state %in% c("Insufficient Prior Records", "No Good Match", "Foreign"))) 

write_csv(unique_observer_years_filtered, "unique_observer_years_filtered.csv")
                                      
county_centroids <- county_sf %>%
  mutate(centroid = st_centroid(geometry)) %>%
  cbind(as.data.frame(st_coordinates(.$centroid)[,c("X", "Y")]))

rd4_centroids <- rd4_sf %>%
  mutate(centroid = st_centroid(geometry)) %>%
  cbind(as.data.frame(st_coordinates(.$centroid)[,c("X", "Y")]))

us_county_populations <- get_decennial(geography = "county",
                                       variables = c("P001001"),
                                       geometry = FALSE) %>%
  rename(county = NAME, 
         population = value) %>%
  mutate(county = sub(" City and County,", ",", county)) %>%
  mutate(county = sub(" County,", ",", county)) %>%
  mutate(county = sub(" City and Borough,", ",", county)) %>%
  mutate(county = sub(" Borough,", ",", county)) %>%
  mutate(county = sub(" Parish,", ",", county)) %>%
  mutate(county = sub(" Census Area,", ",", county)) %>%
  mutate(county = sub(" Municipality,", ",", county)) %>%
  {.[.$county == "Baltimore city, Maryland", "county"] <- "Baltimore City, Maryland"; .} %>%
  {.[.$county == "St. Louis city, Missouri", "county"] <- "St. Louis City, Missouri"; .} %>%
  {.[.$county == "Fairfax city, Virginia", "county"] <- "Fairfax City, Virginia"; .} %>%
  {.[.$county == "Franklin city, Virginia", "county"] <- "Franklin City, Virginia"; .} %>%
  {.[.$county == "Richmond city, Virginia", "county"] <- "Richmond City, Virginia"; .} %>%  
  {.[.$county == "Roanoke city, Virginia", "county"] <- "Roanoke City, Virginia"; .} %>%
  mutate(county = sub(" city,", ",", county)) %>%
  {.[.$county == "LaSalle, Illinois", "county"] <- "La Salle, Illinois"; .} %>%
  {.[.$county == "La Salle, Louisiana", "county"] <- "LaSalle, Louisiana"; .} %>%
  {.[.$county == "Shannon, South Dakota", "county"] <- "Oglala Lakota, South Dakota"; .} %>% # Name change 2015
  {.[.$county == "Skagway, Alaska", "county"] <- 
    "Skagway-Hoonah-Angoon, Alaska"; .} %>% # Judgement call, boundary change
  {.[.$county == "Prince of Wales-Hyder, Alaska", "county"] <- 
    "Prince of Wales-Outer Ketchikan, Alaska"; .} %>% # Judgement call, boundary change
  select(county, population) %>%
  rbind(c("Guam, Guam", 159358)) %>% # https://www.census.gov/newsroom/releases/archives/2010_census/cb11-cn179.html
  rbind(c("Western, American Samoa", 55519)) %>% # https://www.census.gov/newsroom/releases/archives/2010_census/cb11-cn177.html
  mutate(population = as.integer(population))

# 2011 Canadian Census population data from 
# https://www12.statcan.gc.ca/census-recensement/2011/dp-pd/prof/details/download-telecharger/comprehensive/comp-csv-tab-dwnld-tlchrgr.cfm?Lang=E#tabs2011

canada_cd_populations <- read_csv("98-316-XWE2011001-701_CSV/98-316-XWE2011001-701.csv",
                                  skip = 1,
                                  col_names = TRUE,
                                  col_types = "-cc--c-n-----",
                                  locale = locale(encoding = "latin1")) %>%
  filter(Characteristic == "Population in 2011") %>%
  mutate(state = Prov_Name, 
         county = CD_Name, 
         population = as.integer(Total)) %>%
  filter(!(state %in% c("Northwest Territories", "Nunavut", "Yukon"))) %>%
  {.[.$county == "Lajemmerais", "county"] <- "Marguerite-D'Youville"; .} %>%
  mutate(county = paste0(county, ", ", state)) %>%
  select(county, population)

district_county_distances <- expand.grid(district_name = unique(rd4_sf$district_name),
                                         origin = unique(county_centroids$county)) %>%
  mutate(distance = map2_dbl(district_name,
                             origin,
                             ~ as.double(geodist(rd4_centroids[rd4_centroids$district_name == as.character(.x), c("X", "Y")],
                                                 county_centroids[county_centroids$county == as.character(.y), c("X", "Y")],
                                                 measure = "haversine"))/1000)) %>%
  left_join(us_county_populations,
            by = c("origin" = "county")) %>%
  left_join(canada_cd_populations,
            by = c("origin" = "county")) %>%
  mutate(population = ifelse(is.na(population.x), 
                             population.y, 
                             population.x)) %>%
  select(-c(population.x, population.y)) %>% #CHECK
  filter(!is.na(population))

rd4_representative_distances <- district_county_distances %>%
  mutate(bin = case_when(distance <= 100 ~ 100,
                         distance <= 500 ~ 500,
                         distance <= 1000 ~ 1000,
                         distance <= 3000 ~ 3000,
                         TRUE ~ 11000)) %>%
  group_by(district_name, bin) %>%
  mutate(representative_dist = sum(population * distance)/sum(population)) %>%
  select(district_name, bin, representative_dist) %>%
  distinct()

write_csv(rd4_representative_distances, "rd4_representative_distances.csv")

rd4_zone_populations <- district_county_distances %>%
  mutate(bin = case_when(distance <= 100 ~ 100,
                         distance <= 500 ~ 500,
                         distance <= 1000 ~ 1000,
                         distance <= 3000 ~ 3000,
                         TRUE ~ 11000)) %>%
  group_by(district_name, bin) %>%
  summarise(bin_pop = sum(population))

write_csv(rd4_zone_populations, "rd4_zone_populations.csv")

rd4_zone_counts <- unique_observer_years_filtered %>%
  left_join(district_county_distances, by = c("district_name", "origin")) %>%
  mutate(bin = case_when(distance <= 100 ~ 100,
                         distance <= 500 ~ 500,
                         distance <= 1000 ~ 1000,
                         distance <= 3000 ~ 3000,
                         TRUE ~ 11000)) %>%
  group_by(district_name, year, bin) %>%
  mutate(bin_count = n()) %>%
  ungroup() %>%
  select(district_name, year, bin, bin_count) %>%
  distinct() %>%
  right_join(expand.grid(district_name = rd4_sf$district_name, 
                         year = 2008:2021,
                         bin = c(100, 500, 1000, 3000, 11000)),
             by = c("district_name", "year", "bin")) %>%
  mutate(bin_count = ifelse(is.na(bin_count), 0, bin_count),
         treat = ifelse(district_name == "Minidoka Ranger District" &
                          year >= 2017, 1, 0),
         ruby_treat = ifelse(district_name == "Mountain City-Ruby Mountains-Jarbidge Ranger District" &
                              year >= 2017, 1, 0)) %>%
  mutate(treat_bins = as.factor(ifelse(treat == 1, bin, 0)),
         ruby_treat_bins = as.factor(ifelse(ruby_treat == 1, bin, 0))) %>%
  left_join(rd4_zone_populations, by = c("district_name", "bin")) %>%
  mutate(bin_pop_100k = bin_pop/100000) %>%
  mutate(bin_per_100k = ifelse(bin_count == 0, 0, bin_count/bin_pop_100k))

write_csv(rd4_zone_counts, "rd4_zone_counts.csv")

rd4_zone_counts_wide <- rd4_zone_counts %>%
  select(district_name, year, bin, bin_per_100k) %>%
  pivot_wider(names_from = bin, values_from = bin_per_100k) %>%
  rename(bin100 = `100`,
         bin500 = `500`,
         bin1000 = `1000`,
         bin3000 = `3000`,
         bin11000 = `11000`)

