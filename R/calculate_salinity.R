# get the mean bottom salinity to be used in DKI calculations

library(dplyr)
library(tidyr)
library(readr)

# read salinity data (not saved in github repository)
psal <- read_delim("data/HentData_LFLXTQXRGW.csv", 
                   delim = ";", escape_double=FALSE, 
                   locale=locale(decimal_mark = ",", 
                                 grouping_mark = ".", 
                                 encoding = "WINDOWS-1252"), 
                   trim_ws = TRUE, 
                   show_col_types = FALSE)

psal <- psal %>% 
  select(ov_id=ObservationsStedNr,
         ov_navn=ObservationsStedNavn, 
         Længde, Bredde,
         Dato,
         z = `Dybde (m)`,
         psal = KorrigeretResultat)

# find deepest measurement for each profile
psal <- psal %>% 
  group_by(ov_id, ov_navn, Længde, Bredde, Dato) %>%
  arrange(desc(z)) %>%
  slice(1) %>%
  ungroup()

# find average for each ov
psal <- psal %>% 
  group_by(ov_id, ov_navn) %>%
  summarise(lon=median(Længde), 
            lat=median(Bredde), 
            psal=mean(psal, na.rm=T), .groups = "drop")

# save the salinity results

saveRDS(psal, file="data/psal_mean.Rds")



