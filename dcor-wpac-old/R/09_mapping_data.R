# Load necessary libraries
library(tidyverse)    # For data manipulation (dplyr) and plotting (ggplot2)
library(sf)           # For handling spatial data
library(rnaturalearth) # For world map data
library(ggspatial)    # For scale bars and north arrows
library(ggrepel)      # For non-overlapping labels
library(patchwork)    # For combining main map and inset map
library(tidygeocoder) # For geocoding the sample locations
library(ggsn)

#Load non-geocoded sample data
all_samples <- read_csv("data-raw/final_sample_map_data.csv")
head(all_samples)

sample_counts<- all_samples %>%
  count(Locality, Region, Country, name="SampleCount")

# Create a new dataframe with only the unique locations
unique_locations <- all_samples %>%
  distinct(Locality, Region, Country, Stratum)

# Combine location columns into a single address string for better accuracy.
# Sometimes just a 'Locality' and 'Country' is enough. Test what works best.
unique_locations <- unique_locations %>%
  mutate(full_address = paste(Locality, Region, Country, sep = ", "))

# Geocode the unique addresses
# Export unique addresses to CSV to add Lat/Long to
#write.csv(unique_locations, "data-raw/geocoded_locations.csv")

# Load your geocoded and summarized data
locations <- read_csv("data-raw/geocoded_locations.csv")
colnames(locations)<-c("Ind"  ,       "Locality"  ,   "Region"    ,   "Country"   ,   "Stratum"  ,    "full_address" ,"Latitude",     "Longitude")
head(locations)
locations_with_counts<-locations%>%
  left_join(sample_counts, by=c("Locality", "Region", "Country"))

country_summary <- locations_with_counts %>%
  group_by(Country) %>%
  summarise(
    # Sum the samples for a total count per country
    TotalSampleCount = sum(SampleCount),
    # Calculate a central point for plotting
    Latitude = mean(Latitude),
    Longitude = mean(Longitude)
  )


library(ggOceanMaps)
library(ggspatial)
map <- basemap(limits = c(90, -100, -30, 50), rotate = TRUE)

# Build the final map
main_map <- map +
  geom_spatial_point(
    data = as.data.frame(country_summary), 
    aes(x = Longitude, y = Latitude, fill = Country, size = TotalSampleCount), # Changed to 'fill'
    shape = 21,         # Use a circle shape with an outline
    color = "black",    # Set the outline color to black
    alpha = 0.8
  ) +
  scale_size_continuous(range = c(5, 15), name = "Number of Samples") + # Increased the size range
  
  # Add all labels here
  labs(
    title = "Sampling Locations",
    fill = NULL,  # Changed to 'fill' to match the aesthetic
    x = "Longitude", 
    y = "Latitude"
  ) +
  
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    panel.grid = element_blank()
  )

# Display the map
main_map


birds_head_data <- locations_with_counts %>%
  filter(Stratum %in% c("Bird's Head-Summer", "Bird's Head-Winter"))%>%
  group_by(Stratum, Latitude, Longitude) %>%
  summarise(SampleCount = sum(SampleCount, na.rm = TRUE),
  .groups = 'drop' # This tells dplyr to ungroup after summarizing
  )

head(birds_head_data)

inset_map <- basemap(limits = c(130, 138, -3, 1)) +
  
  # Layer 1: Draw large circles
  geom_spatial_point(
    data = as.data.frame(birds_head_data),
    aes(x = Longitude, y = Latitude, fill = Stratum),
    size = 15,
    shape = 21,
    color = "black",
    crs = "EPSG:4269"
  ) +
  
  # Layer 2: Add the sample count numbers
  geom_text(
    data = as.data.frame(birds_head_data),
    aes(x = Longitude, y = Latitude, label = SampleCount),
    color = "white",
    size = 5,
    fontface = "bold"
  ) +
  
  # Change the legend title
  labs(fill = NULL) +
  
  # Theme adjustments
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    
    legend.position = c(0.95, 0.95), 
    # Set the legend's own anchor point to its top-right corner
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = alpha("white", 0.5)),
    panel.grid = element_blank()
  ) 

# Display the map
inset_map

main_map
inset_map
