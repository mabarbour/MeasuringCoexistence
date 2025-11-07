
## Setup ----

# load libraries
library(RCurl)
library(tidyverse)

# download data
aphid_url <- getURL("https://raw.githubusercontent.com/mabarbour/EcoEvoStability/refs/heads/main/data/aphid_data_final.csv")

# manipulate data for analysis. Modified from:
# https://github.com/mabarbour/EcoEvoStability/blob/main/code/manage_data.R
aphid_data <- read.csv(text = aphid_url) %>%
  mutate(# Estimate aphid abundance for the entire cage. First calculate aphids per plant then multiply by 16 (total number of plants per cage)
    R_cage = as.integer((R*aphids_per_click)/plants*16),
    Y_cage = as.integer((Y*aphids_per_click)/plants*16),
    G_cage = as.integer((G*aphids_per_click)/plants*16),
    Aphids_cage = R_cage + Y_cage + G_cage, 
    Treatment = ifelse(ID %in% c("3RP", "8RP", "13RP", "18RP", "23RP", "28RP"), "RP",
                       ifelse(ID %in% c("4YP", "9YP", "14YP", "19YP", "24YP", "29YP"), "YP",
                              ifelse(ID %in% c("5GP", "10GP", "15GP", "20GP", "25GP", "30GP"), "GP",
                                     ifelse(ID %in% c("2RYGP", "7RYGP", "12RYGP", "17RYGP", "22RYGP", "27RYGP"), "RYGP",
                                            ifelse(ID %in% c("1RYG", "6RYG", "11RYG", "16RYG", "21RYG", "26RYG"), "RYG", NA)))))) %>%
  group_by(ID) %>%
  mutate(# determine lagged abundances for MAR(1) models
    R_cage_lag1 = lag(R_cage, n = 1, order_by = week),
    Y_cage_lag1 = lag(Y_cage, n = 1, order_by = week),
    G_cage_lag1 = lag(G_cage, n = 1, order_by = week),
    Aphids_cage_lag1 = R_cage_lag1 + Y_cage_lag1 + G_cage_lag1) %>%
  ungroup() %>%
  #filter(week > 0, week < 4 | week > 6) %>% # exclude weeks of non-sampling and heat-wave
  # mutate and rename variables to facilitate model processing
  mutate(ptoid_in = ifelse(Treatment %in% c("RYGP","RP","YP","GP") & week > 3, 1, 0),
         Aphids_sub = ifelse(Aphids_cage_lag1 > 0 | aphids_extinct == 0, 1, 0), 
         # log1p abundance variables
         lnRt1 = log1p(R_cage),
         lnRt = log1p(R_cage_lag1),
         lnYt1 = log1p(Y_cage),
         lnYt = log1p(Y_cage_lag1),
         lnGt1 = log1p(G_cage),
         lnGt = log1p(G_cage_lag1)) %>% 
  filter(ptoid_in == 0, Aphids_sub == 1) %>%
  # focus on select data columns for subsequent analyses
  select(week, ID, Treatment, 
         R_cage, R_cage_lag1, Y_cage, Y_cage_lag1, G_cage, G_cage_lag1,
         lnRt1, lnRt, lnYt1, lnYt, lnGt1, lnGt) %>%
  drop_na()

# view data
aphid_data
