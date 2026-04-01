# Tidy tuesdays data exploration

library(tidyverse)


# --------------------
# Olympics Milano
# --------------------
dataset_url <- "https://raw.githubusercontent.com/chendaniely/olympics-2026/refs/heads/main/data/final/olympics/olympics_events.csv"
df <- readr::read_csv(dataset_url)

summary(df)

events_per_sport <- df %>%
  group_by(discipline_name) %>%
  summarize(n = n())

ggplot(events_per_sport, aes (x = reorder(discipline_name, n), y = n))+
  geom_col()


# Explorations
# how many events of each sport give a medal
# what venue was the most used
# what day had the most events


# --------------------
# Ireland grants
# --------------------
install.packages("tidytuesdayR")
eire <- tidytuesdayR::tt_load('2026-02-24')

sfi_grants <- eire$sfi_grants


# --------------------
# The Simpsons
# --------------------
#https://toddwschneider.com/posts/the-simpsons-by-the-data/
  
tuesdata <- tidytuesdayR::tt_load('2025-02-04')

simpsons_characters <- tuesdata$simpsons_characters
simpsons_episodes <- tuesdata$simpsons_episodes
simpsons_locations <- tuesdata$simpsons_locations
simpsons_script_lines <- tuesdata$simpsons_script_lines



# --------------------
# Pokemon
# --------------------
tuesdata <- tidytuesdayR::tt_load('2025-04-01')
pokemon_df <- tuesdata$pokemon_df


# --------------------
# British library funding
# --------------------
# https://anjackson.net/2024/11/27/updating-the-data-on-british-library-funding/

tuesdata <- tidytuesdayR::tt_load('2025-07-15')
bl_funding <- tuesdata$bl_funding


# --------------------
# IKEA
# --------------------
tuesdata <- tidytuesdayR::tt_load(2020, week = 45)
ikea.csv <- tuesdata$ikea.csv
