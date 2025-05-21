library(dplyr)
library(tidyr)
library(ggplot2)
library(ambiR)
library(readr)
library(httr)
library(jsonlite) 
library(sf)
library(fuzzyjoin)
library(rgbif)
library(worrms)
library(patchwork)

source("R/utils.R")


bf <- read_delim("data/HentData_WDXSYOYGJE.csv", 
                 delim = ";", escape_double=FALSE, 
                 locale=locale(decimal_mark = ",", 
                               grouping_mark = ".", 
                               encoding = "WINDOWS-1252"), 
                 trim_ws = TRUE, 
                 show_col_types = FALSE)


res <- match_species(bf, column_species = "Artsnavn")

