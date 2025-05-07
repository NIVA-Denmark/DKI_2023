library(dplyr)
library(tidyr)
library(ggplot2)
library(ambiR)
library(readr)
library(httr)
library(jsonlite) 
library(sf)

# pak::pak("niva-denmark/ambiR", upgrade=T)
# install.packages("C:/Data/repositories/ambiR_0.1.0.tar.gz")

source("R/utils.R")

# read salinity mean values (calculate_salinity.R)

psal <- readRDS("data/psal_mean.Rds")

psal <- psal %>% 
  sf::st_as_sf(coords=c("lon","lat"), crs=sf::st_crs(4326)) %>%
  sf::st_transform(crs=sf::st_crs(25833))

  
## DKI

# Read species observation data from ODA.

bf <- read_delim("data/HentData_WDXSYOYGJE.csv", 
                 delim = ";", escape_double=FALSE, 
                 locale=locale(decimal_mark = ",", 
                               grouping_mark = ".", 
                               encoding = "WINDOWS-1252"), 
                 trim_ws = TRUE, 
                 show_col_types = FALSE)

#  check sampling methods used and numbers of records approved

stats <- bf %>%
  group_by(`Prøveredskabsareal (cm2)`,Kvalitet) %>% 
  summarise(n=n(), .groups="drop") %>%
  pivot_wider(names_from = Kvalitet, values_from = n)

print(stats)

#' filter data to exclude non-approved species registrations
#' this is not run unless you change the conditional expression

if(FALSE){
  bf <- bf %>%
    filter(Kvalitet == "Godkendt")
}


# get unique locations for bottom fauna samples
bf_posn <- bf %>%
  select(ov_id=ObservationsStedNr,
         ov_navn=ObservationsStedNavn, 
         lon = `Prøvetagnings Længdegrad`,
         lat = `Pprøvtagnings Breddegrad`) %>%
  group_by(ov_id, ov_navn) %>%
  summarise(lon = median(lon),
            lat=median(lat), .groups="drop")
         
bf_posn <- bf_posn %>%
  sf::st_as_sf(coords=c("lon","lat"), crs=sf::st_crs(4326)) %>%
  sf::st_transform(crs=sf::st_crs(25833))

# find closest salinity point for each bottom fauna position
closest <- sf::st_nearest_feature(bf_posn, psal)

bf_posn$psal <- psal$psal[closest]

bf_psal <- bf_posn %>%
  sf::st_drop_geometry() %>%
  select(ov_id, psal)

# back to the sample data :-)
# ... select relevant columns

bf <- bf %>%
  select(ov_id=ObservationsStedNr,
         ov_navn=ObservationsStedNavn, Dato, Tid, 
         sample_nr=Prøvetagningsnummer, sample_method=Prøvetagningsudstyr,
         sample_area_cm2 =  `Prøveredskabsareal (cm2)`, n=`Antal (stk)`, 
         Artsrække, Artsnavn)

# replace_names

bf <- bf %>%
  mutate(Artsnavn=stringr::str_replace_all(Artsnavn, "indet.", "sp.")) %>%
  mutate(Artsnavn=.fix_name(Artsnavn)) %>%
  mutate(Artsnavn = .valid_species(Artsnavn))


# ------------ DKI with pooled samples --------

# For DKI v2, samples are pooled


bf_p <- bf %>%
  mutate(samples_per_group = round(1000/sample_area_cm2, digits = 0))


# add pooled sample ID
bf_p <- bf_p %>%
  DKI_group_samples(
    group_vars=c("ov_navn","ov_id","Dato","sample_method"),
    column_sample_id = "sample_nr",
    column_samples_per_group="samples_per_group")

# summarise within pooled samples
bf_p <- bf_p %>%
  group_by(ov_id, ov_navn, Dato, Tid, 
           sample_method, grp_id,
           Artsrække, Artsnavn) %>%
  summarise(n=sum(n), .groups="drop")


# first run with AMBI

ambi_res_p <- ambiR::AMBI(bf_p, 
                    by=c("ov_id", "ov_navn", "Dato", "Tid"),
                    var_species = "Artsnavn",
                    var_count = "n",
                    var_rep = "grp_id",
                    interactive = F)

not_found <- ambi_res_p$matched %>%
  filter(is.na(group)) %>%
  distinct(Artsrække, Artsnavn)

ArtsnavnWoRMS <- not_found$Artsnavn %>%
  purrr::map(.GetSpeciesName, 
             .progress ="Checking WoRMS species names")

not_found$Artsnavn2 <- ArtsnavnWoRMS

not_found <- not_found %>%
  mutate(Artsnavn2 = .fix_name(Artsnavn2))

not_found <- not_found %>%
  mutate(Artsnavn2 = ifelse(Artsnavn==Artsnavn2,
                            NA_character_,
                            Artsnavn2))

new_names <- not_found %>%
  filter(!is.na(Artsnavn2)) %>%
  select(Artsnavn, Artsnavn2)

# replace names where new ones were found

bf_p <- bf_p %>%
  left_join(new_names, by="Artsnavn") %>%
  mutate(Artsnavn=ifelse(is.na(Artsnavn2),Artsnavn,Artsnavn2)) %>%
  select(-Artsnavn2)

ambi_res_p <- ambiR::AMBI(bf_p, 
                     by=c("ov_id", "ov_navn", "Dato", "Tid"),
                     var_species = "Artsnavn",
                     var_count = "n",
                     var_rep = "grp_id",
                     interactive = F)


ambi_p <- ambi_res_p$AMBI


dki_p <- ambi_p %>%
  select(ov_id, ov_navn, Dato, Tid, AMBI, H, S, fNA, N) %>%
  left_join(bf_psal, by="ov_id") %>%
  mutate(DKI = DKI2(AMBI, H, N, psal))

# ------------ DKI with individual samples --------

# reuse the corrected names found earlier
bf <- bf %>%
  left_join(new_names, by="Artsnavn") %>%
  mutate(Artsnavn=ifelse(is.na(Artsnavn2),Artsnavn,Artsnavn2)) %>%
  select(-Artsnavn2)

ambi_res <- ambiR::AMBI(bf, 
                     by=c("ov_id", "ov_navn", "Dato", "Tid"),
                     var_species = "Artsnavn",
                     var_count = "n",
                     var_rep = "sample_nr",
                     interactive = F)


ambi <- ambi_res$AMBI


dki <- ambi %>%
  select(ov_id, ov_navn, Dato, Tid, AMBI, H, S, fNA, N) %>%
  left_join(bf_psal, by="ov_id") %>%
  mutate(DKI = DKI2(AMBI, H, N, psal))

# ---- compare pooled / individual results ------

# compare DKI calculated from individual samples with DKI from pooled samples

cf <- dki %>%
  select(ov_id, ov_navn, Dato, Tid, psal, DKI) %>%
  left_join(dki_p %>%
              select(ov_id, Dato, Tid, DKI_pooled=DKI), 
            by=c("ov_id", "Dato", "Tid"))

p <- ggplot() +
  geom_point(data=cf, 
             aes(x=DKI, y=DKI_pooled, colour = psal)) +
  scale_color_distiller(palette = "Spectral", 
                        limits=c(5,35), breaks=c(5,15,25,35)) +
  theme_minimal() +
  theme(legend.position = "inside", 
        legend.position.inside = c(0.8,0.1),
        legend.justification = c(0,0))

ggsave(p, filename="compare_DKI.png", height=10, width=10, 
       units="cm", dpi=300, bg="white")
