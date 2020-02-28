# =-------------------------------------------------------------------
# Climate Indexes Project
# Alejandra E. - Harold A. 
# Februrary 2020
# =-------------------------------------------------------------------

# =--------------------
# Packages 
library(tidyverse)
library(raster)
library(ncdf4)
library(sf)
library(future)
library(furrr)
library(lubridate)
library(glue)
library(cowsay)
# library(fs)
# =--------------------

# Countries: 
# Pakistan - Benin, Burkina Faso, Ivory Coast, Mali, Mali, Cameroon, Togo, 
# Nigeria, Ethiopia, Kenya, Malawi, Mozambique, Zambia, Tunisia, India, Vietnam


all_countries <- getData("ISO3")%>%
  as.tibble() %>% 
  filter(NAME %in% c('Pakistan', 'Benin', 'Burkina Faso', 'Ivory Coast', 'Mali', 'Mali', 'Cameroon', 'Togo', 
                     'Nigeria', 'Ethiopia', 'Kenya', 'Malawi', 'Mozambique', 'Zambia', 'Tunisia', 'India', 'Vietnam'))

# Descargas de los shp. 
tictoc::tic()
ctry_shps <-  do.call("bind", lapply(all_countries$NAME, 
                                   function(x) getData('GADM', country=x, level=0))) %>% 
  sf::st_as_sf()
tictoc::toc()

colombia <- getData('GADM', country= 'Colombia', level=0) %>% sf::st_as_sf()
 
all_world <- shapefile('//dapadfs/workspace_cluster_9/Coffee_Cocoa2/_cam/_shp/_base/all_countries.shp') %>% 
  sf::st_as_sf()




# =--------------------
path_Chirps <- '//dapadfs/data_cluster_4/observed/gridded_products/chirps/daily/'
path_CHITS <- '//catalogue/BaseLineData_cluster04/GLOBAL/Climate/CHIRTS/'



list.files(path_CHITS)

# tibble(name = list.files(path_Chirps, pattern = 'chirps-v2.0.')) %>% 
#   mutate(new = str_remove(name, 'chirps-v2.0.') %>% str_remove(., '.tif')) %>%
#   filter()

test <- tibble(Date = seq(ymd('1985-01-01'),ymd('2015-12-31'),by='day') ) %>% 
  mutate(date_mod = str_replace_all(Date, '-','.'), 
         path_rain = glue::glue('{path_Chirps}chirps-v2.0.{date_mod}.tif')) %>% 
  dplyr::select(-date_mod)


raster(test$path_rain[1]) %>% crop(ctry_shps) %>% mask(ctry_shps)  %>% plot()


# tictoc::tic()
# chirps <- stack(test$path_rain) %>% 
#   crop(ctry_shps)
# tictoc::toc()



head(test)





