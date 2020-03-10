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
library(fst)
# =--------------------


# Ruta Principal para guardados: 
principal_path <- '//dapadfs/workspace_cluster_8/climateriskprofiles/data/'


# \\dapadfs\workspace_cluster_8\climateriskprofiles\data\Chirps_Chirts


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



test %>% filter(Date < '1986-01-01') 
  

path_test <- filter(test, Date == '1985-01-01') %>% 
  dplyr::select(path_rain) %>% 
  pull(path_rain)



# =- Aqui estamos creando la función de leida de los datos. 


date <- "1985.01.01"


do_raster_to_table <- function(date, path_Chirps, path_CHITS){
  
  a <- raster(glue::glue('{path_Chirps}chirps-v2.0.{date}.tif')) %>% crop(ctry_shps) 
  bc <- stack(glue::glue('{path_CHITS}Tmax.{date}.tif'), glue::glue('{path_CHITS}Tmin.{date}.tif')) %>% crop(ctry_shps) 
  
  abc <- stack(a, bc) %>% mask(ctry_shps)
  
  names(abc) <- c('prec', 'tmax', 'tmin')
  points <- abc %>% rasterToPoints() %>% as.tibble() %>% 
    mutate(Date = date %>% lubridate::as_date()) %>% dplyr::select(Date, everything()) %>% 
    mutate_each(funs(replace(., . == -9999, NA_real_)), -Date, -x, -y) %>% 
    drop_na() %>% 
    mutate(id = 1:nrow(.)) %>% 
    dplyr::select(id, everything())
  
  fst::write.fst(x = points, path = glue::glue('D:/OneDrive - CGIAR/Desktop/P_indices_H/fst/D_{date}.fst'))

return(points)}


tictoc::tic()
plan(multiprocess)
options(future.globals.maxSize= 891289600)
ten_days <- tibble(Date = seq(ymd('1985-01-01'),ymd('1985-01-03'),by='day') %>% str_replace_all('-', '.')) %>% 
  mutate(climate = furrr::future_map(.x = Date, .f = do_raster_to_table, path_Chirps = path_Chirps, path_CHITS = path_CHITS))
gc()
gc(reset = T)
tictoc::toc() # 6.57


# st_union(ctry_shps, points)
# ggplot() + 
#   geom_raster(data = points, aes(x = x, y = y, fill =  prec)) + 
#   geom_sf(data = ctry_shps, color = gray(.5), fill = NA) + 
#   theme_bw()
# 
# 
# ggplot() + 
#   geom_raster(data = points, aes(x = x, y = y, fill =  tmax)) + 
#   geom_sf(data = ctry_shps, color = gray(.5), fill = NA) + 
#   theme_bw()
# 
# 
# ggplot() + 
#   geom_raster(data = points, aes(x = x, y = y, fill =  tmin)) + 
#   geom_sf(data = ctry_shps, color = gray(.5), fill = NA) + 
#   theme_bw()





# Para Colombia... (comparación). 


catalog <- read_csv('D:/OneDrive - CGIAR/Desktop/P_indices_H/Profiles_AA/st_cuencas.csv') 




list.files(path_CHITS)


# tibble(Date = seq(ymd('1985-01-02'),ymd('2015-12-31'),by='day') %>% str_replace_all('-', '.')) 



date <- '1985.01.11'
coord <- catalog %>% dplyr::select( LONG, LAT) %>% raster::coordinates()
parallel::detectCores()

extract_stations <- function(date, coord, path_CHITS){
  Temp <- stack(glue::glue('{path_CHITS}Tmax.{date}.tif'), glue::glue('{path_CHITS}Tmin.{date}.tif'))
  
  data <- raster::extract(Temp, coord) %>% 
    as_tibble() %>% 
    set_names(c('tmax', 'tmin')) %>% 
    mutate(id = 1:nrow(.)) %>% 
    bind_cols(as_tibble(coord) , . ) %>% 
    dplyr::select(id, everything())
  
  print(date)
return(data)}




1500000
# raster('//catalogue/BaseLineData_cluster04/GLOBAL/Climate/CHIRTS/Tmax.1985.09.07.tif')



raster('//catalogue/BaseLineData_cluster04/GLOBAL/Climate/CHIRTS/Tmax.1999.10.16.tif')

raster('//catalogue/BaseLineData_cluster04/GLOBAL/Climate/CHIRTS/Tmax.2000.10.18.tif')



x <- list.files('//catalogue/BaseLineData_cluster04/GLOBAL/Climate/CHIRTS/', pattern ='Tmax',full.names = TRUE)


# file.size(x[2])
# y <- x[sapply(x, file.size) < 74901352]
# fileList <- x [1:5]
# sapply(x [1:5], file.size)

tictoc::tic()
x1 <- x %>% 
  as_tibble() %>% 
  mutate(size = purrr::map(.x = value, .f = file.size))
tictoc::toc()

datos <- x1 %>% 
  unnest() %>% 
  filter(size < 74901352) %>% 
  mutate(name = str_remove(value, '//catalogue/BaseLineData_cluster04/GLOBAL/Climate/CHIRTS/')) %>% 
  dplyr::select(-value)


write_csv(datos, 'sizes.csv')

datos_t <- datos %>% pull(name) %>% str_remove('Tmax.') %>% str_remove('.tif')


# =-------------------------------------------------------

tictoc::tic()
# plan(multiprocess)
# options(future.globals.maxSize= 891289600)
points <- tibble(Date = seq(ymd('1985-01-01'),ymd('2015-12-31'),by='day') %>% str_replace_all('-', '.')) %>% 
  filter(!(Date %in% datos_t)) %>% 
  mutate(climate = purrr::map(.x = Date, .f = extract_stations, coord = coord, path_CHITS = path_CHITS))
# gc()
# gc(reset = T)
tictoc::toc() # 3.21 horas (estaba corriendo otra cosa a la vez).


points_mod <- points %>% 
  unnest() 


IDEAM <- readRDS('D:/OneDrive - CGIAR/Desktop/P_indices_H/Profiles_AA/IDEAM_temp.rds') %>% 
  nest(-codS)  %>% 
  mutate(codS = as.numeric(codS)) %>% 
  rename('IDEAM' = 'data')

all_Date <- IDEAM %>% 
  dplyr::select(IDEAM) %>% filter(row_number() == 1) %>% unnest() %>% dplyr::select(Date)




satelital_data <- points_mod %>% 
  mutate(Date = as_date(Date)) %>% 
  rename('x' = 'LONG' ,   'y' = 'LAT' ,'tmax_C' = 'tmax', 'tmin_C' = 'tmin') %>% 
  full_join(dplyr::select(catalog, ID,NAME) %>% rename('id' =  'ID', 'codS' = 'NAME'), .) %>% 
  nest(-id, -codS, -x, -y) %>% 
  mutate(data = purrr::map(.x = data, .f = function(x, y){full_join(y, x)}, y = all_Date))




all_data_for_ts <-  inner_join(satelital_data %>% unnest, IDEAM %>% unnest) 
write_csv(x = all_data_for_ts, path = 'data_for_compT.csv')


data <- all_data_for_ts %>%
  nest(-id, -codS, -x , -y) %>% dplyr::select(data) %>% filter(row_number() == 1) %>% unnest()


control <- function(data, id){
  # Tmax...
  tmax <- data %>% 
    dplyr::select( Date, contains('tmax')) %>% 
    rename( 'tmax_CHIRTS' = 'tmax_C', 'tmax_Ideam'  = 'tmax_I') %>% 
    mutate(year = lubridate::year(Date), month = lubridate::month(Date), day = lubridate::day(Date))
  
  
  tmax_g <- tmax %>% dplyr::select(-year, -month, -day) %>% 
    gather(key = var, value = value, -Date) %>% 
    ggplot(aes(x = Date, value, colour =  var)) + 
    geom_line() + 
    scale_x_date(date_breaks = "2 year", date_labels = "%Y") +
    labs(title = 'Tmax', x = NULL, y = 'Temperature (°C)', colour = NULL) + 
    theme_bw() + 
    theme(legend.position = 'top')
  
  
  wna_tmax <- drop_na(tmax) 
  
  # Calculo de indices. 
  Index_tmax <- wna_tmax %>% # Calculo con la base de datos... en general. 
    group_by(month) %>% 
    mutate(less = tmax_Ideam-tmax_CHIRTS) %>% 
    summarise(Kolmogorov = ks.test(tmax_CHIRTS, tmax_Ideam)$p.value,
              Wald_Wolfowitz = DescTools::RunsTest(tmax_CHIRTS, tmax_Ideam, exact=FALSE)$p.value, 
              cor_p = cor(tmax_CHIRTS, tmax_Ideam),
              cor_k = cor(tmax_CHIRTS, tmax_Ideam, method = 'kendall') )
  
  RMSE_tmax <- wna_tmax %>% group_by(month) %>% yardstick::rmse(tmax_CHIRTS, tmax_Ideam) %>% 
    dplyr::select(month, .estimate) %>% rename(RMSE = '.estimate')
  
  tmax_I <- inner_join(Index_tmax, RMSE_tmax) %>% mutate(var = 'tmax')
  
  
  Index_tmax_y <- wna_tmax %>% # Calculo con la base de datos... en general. 
    group_by(year) %>% 
    mutate(less = tmax_Ideam-tmax_CHIRTS) %>% 
    summarise(Kolmogorov = ks.test(tmax_CHIRTS, tmax_Ideam)$p.value,
              Wald_Wolfowitz = DescTools::RunsTest(tmax_CHIRTS, tmax_Ideam, exact=FALSE)$p.value, 
              cor_p = cor(tmax_CHIRTS, tmax_Ideam),
              cor_k = cor(tmax_CHIRTS, tmax_Ideam, method = 'kendall') ) 
  
  
  RMSE_tmax_y <- wna_tmax %>% group_by(year) %>% yardstick::rmse(tmax_CHIRTS, tmax_Ideam) %>% 
    dplyr::select(year, .estimate) %>% rename(RMSE = '.estimate')
  
  tmax_I_y <- inner_join(Index_tmax_y, RMSE_tmax_y) %>% mutate(var = 'tmax')
  
  
  # Tmin ... 
  tmin <- data%>% 
    dplyr::select( Date, contains('tmin')) %>% 
    rename( 'tmin_CHIRTS' = 'tmin_C', 'tmin_Ideam'  = 'tmin_I') %>% 
    mutate(year = lubridate::year(Date), month = lubridate::month(Date), day = lubridate::day(Date))
  
  
  tmin_g <- tmin %>% dplyr::select(-year, -month, -day) %>% 
    gather(key = var, value = value, -Date) %>% 
    ggplot(aes(x = Date, value, colour =  var)) + 
    geom_line() + 
    scale_x_date(date_breaks = "2 year", date_labels = "%Y") +
    labs(title = 'Tmin', x = NULL, y = 'Temperature (°C)', colour = NULL) + 
    theme_bw() + 
    theme(legend.position = 'top')
  
  
  wna_tmin <- drop_na(tmin) 
  # =---------
  Index_tmin <- wna_tmin %>% # Calculo con la base de datos... en general. 
    group_by(month) %>% 
    mutate(less = tmin_Ideam-tmin_CHIRTS) %>% 
    summarise(Kolmogorov = ks.test(tmin_CHIRTS, tmin_Ideam)$p.value,
              Wald_Wolfowitz = DescTools::RunsTest(tmin_CHIRTS, tmin_Ideam, exact=FALSE)$p.value, 
              cor_p = cor(tmin_CHIRTS, tmin_Ideam),
              cor_k = cor(tmin_CHIRTS, tmin_Ideam, method = 'kendall') )
  
  RMSE_tmin <- wna_tmin %>% group_by(month) %>% yardstick::rmse(tmin_CHIRTS, tmin_Ideam) %>% 
    dplyr::select(month, .estimate) %>% rename(RMSE = '.estimate')
  
  tmin_I <- inner_join(Index_tmin, RMSE_tmin) %>% mutate(var = 'tmin')
  
  # =------------
  Index_tmin_y <- wna_tmin %>% # Calculo con la base de datos... en general. 
    group_by(year) %>% 
    mutate(less = tmin_Ideam-tmin_CHIRTS) %>% 
    summarise(Kolmogorov = ks.test(tmin_CHIRTS, tmin_Ideam)$p.value,
              Wald_Wolfowitz = DescTools::RunsTest(tmin_CHIRTS, tmin_Ideam, exact=FALSE)$p.value, 
              cor_p = cor(tmin_CHIRTS, tmin_Ideam),
              cor_k = cor(tmin_CHIRTS, tmin_Ideam, method = 'kendall') ) 
  
  RMSE_tmin_y <- wna_tmin %>% group_by(year) %>% yardstick::rmse(tmin_CHIRTS, tmin_Ideam) %>% 
    dplyr::select(year, .estimate) %>% rename(RMSE = '.estimate')
  
  tmin_I_y <- inner_join(Index_tmin_y, RMSE_tmin_y) %>% mutate(var = 'tmin')
  
  a <- full_join(tmin_I, tmax_I) %>% dplyr::select(var, everything()) %>% nest(-var) %>% rename('Month' = 'data')
  b <- full_join(tmin_I_y, tmax_I_y) %>% dplyr::select(var, everything()) %>% nest(-var) %>% rename('Year' = 'data')
  
  final_data <- inner_join(a, b) %>% mutate(id) %>% nest(-id)
  
  
  png(paste0('D:/OneDrive - CGIAR/Desktop/P_indices_H/Profiles_AA/Comp_st/' ,id,'.png'), width = 900, height = 720, res = 120)
  gridExtra::grid.arrange(tmax_g, tmin_g)
  dev.off()
  
return(final_data)}


data_each_st <- all_data_for_ts %>%
  nest(-id, -codS, -x , -y) %>% # dplyr::select(data) %>% filter(row_number() == 1) %>% unnest()
  mutate(Comp = purrr::map2(.x = data, .y = codS, .f = control))


# =- 
data_each_st %>% dplyr::select(Comp) %>% unnest() %>% 
  unnest() %>% filter(var == 'tmin') %>%  dplyr::select(id, Month) %>% unnest() %>% 
  write_csv(. , path = 'D:/OneDrive - CGIAR/Desktop/P_indices_H/Profiles_AA/Comp_st/tmin_year.csv')
  








# =--------------------------------------------------
# NASA POWER
# =--------------------------------------------------

# library(nasapower)
# 
# # for(i in 1:9){
# #   Sys.sleep(5)
# #   print(i)
# # }
# 
# tictoc::tic()
# daily_region_ag <- get_power(
#   community = "AG",
#   lonlat = c(150.5, -28.5 , 153.5, -25.5),
#   pars = c("RH2M", "T2M"),
#   dates = c("1985-01-01", "2015-12-31"),
#   temporal_average = "DAILY"
# )
# tictoc::toc() # 3.298
# 
# daily_region_ag %>% dplyr::select(LON, LAT) %>% group_by(LON, LAT) %>% count() %>% View()



#    id       Date         x      y prec      tmax       tmin
test <- fst::read.fst(path = '//dapadfs/workspace_cluster_8/climateriskprofiles/data/Chirps_Chirts/D_1985.01.01.fst')
test %>% as.tibble() %>% filter(row_number() == 1) %>% dplyr::select(x, y)

# Para la descarga de NASA necesito tener este formato. 
#    id       Date         x      y prec      tmax       tmin

# principal_path
 
# //dapadfs/workspace_cluster_8/climateriskprofiles/data/NASA

Nasa_download <- function(id, coords){
  
  lat =  coords$y
  lon =  coords$x
  
  json_file <- paste0("https://power.larc.nasa.gov/cgi-bin/v1/DataAccess.py?&request=execute&identifier=SinglePoint&parameters=ALLSKY_SFC_SW_DWN,T2M_MAX,T2M_MIN&startDate=19850101&endDate=",format(as.Date("2015/12/31"),"%Y%m%d"),"&userCommunity=AG&tempAverage=DAILY&outputList=ASCII&lat=",lat,"&lon=",lon)
  # Esta mostrando un error que no conozco.
  json_data <- jsonlite::fromJSON(json_file)
  
  
  var_NASA <- tibble(id = id, Date = seq(ymd('1985-01-01'),ymd('2015-12-31'),by='day'), x = lon, y = lat)  %>%  
    mutate(tmin_N = json_data$features$properties$parameter$T2M_MIN %>% unlist, 
           tmax_N = json_data$features$properties$parameter$T2M_MAX %>% unlist, 
           srad = json_data$features$properties$parameter$ALLSKY_SFC_SW_DWN %>% unlist) %>% 
    na_if(-99)
  
  fst::write.fst(x = var_NASA, path = glue::glue('//dapadfs/workspace_cluster_8/climateriskprofiles/data/NASA/{id}.fst'))
  
return(var_NASA)}


by_list_nasa <- function(list_n){

    x <- purrr::walk2(.x =  list_n$id, .y = list_n$coords, .f = Nasa_download)
    Sys.sleep(60)

return(x)}


# j <- por_ahora %>% nest(-id) %>% rename(coords = data) %>%
#      mutate(group =  rep(1:2, each = 10) ) %>%
#      group_split(group) %>% purrr::map(.f = by_list_nasa)

id_complete <-test %>% 
  as.tibble() %>% dplyr::select(id, x, y) %>% 
  unique()  %>% nest(-id) %>% rename(coords = data) 


# nrow(id_complete )
# 375856/90
# length(c(rep(1:4176, each = 90),  rep(4177, 16) ) )- 375856

# tictoc::tic()
# NASA <- id_complete  %>% nest(-id) %>% rename(coords = data)  %>% 
#   mutate(group =  c(rep(1:4176, each = 90),  rep(4177, 16) )  ) %>% 
#      group_split(group) %>% purrr::map(.f = by_list_nasa)
# tictoc::toc()



countries_id <- read_csv('id_country.csv')
id_run <- list.files('//dapadfs/workspace_cluster_8/climateriskprofiles/data/NASA/') %>% 
  str_remove('.fst') %>% 
  as.numeric()

#  'Mali','Ethiopia'


id_complete %>% right_join(countries_id) %>% filter(id == 95)


tictoc::tic()
Pakistan <- id_complete %>% right_join(countries_id)%>%
  filter(Country %in% c('Pakistan')) %>% 
  filter(!(id %in% id_run))   %>%  filter(row_number() > 15923) %>% 
  mutate(group =  c(rep(178:353, each = 90), rep(354, 83) )  ) %>% 
  group_split(group) %>% purrr::map(.f = by_list_nasa)
tictoc::toc()


# 
# tictoc::tic()
# Mali <- id_complete %>% inner_join(countries_id) %>%
#   filter(Country %in% c('Mali')) %>% 
#   filter(!(id %in% id_run))  %>%  # nrow()
#   mutate(group =  c(rep(1:473, each = 90), rep(474, 69) )  ) %>% 
#   group_split(group) %>% purrr::map(.f = by_list_nasa)
# tictoc::toc()
# 
# 
# tictoc::tic()
# Ethiopia <- id_complete %>% inner_join(countries_id) %>%
#   filter(Country %in% c('Ethiopia')) %>% 
#   filter(!(id %in% id_run))  %>%  # nrow()
#   mutate(group =  c(rep(1:410, each = 90), rep(411, 31) )  ) %>% 
#   group_split(group) %>% purrr::map(.f = by_list_nasa)
# tictoc::toc()
# 
# 








  




# =--------------------------------------------------------------
#    id       Date         x      y prec      tmax       tmin
test <- fst::read.fst(path = '//dapadfs/workspace_cluster_8/climateriskprofiles/data/Chirps_Chirts/D_1985.01.01.fst')
test <- test %>% as.tibble() #%>% filter(row_number() == 1) # %>% dplyr::select(x, y)

spg <- test %>% 
  dplyr::select(x, y, id)

coordinates(spg) <- ~ x + y
gridded(spg) <- TRUE
# coerce to raster
rasterDF <- raster(spg)
rasterDF
# plot(rasterDF)
#  Por ahora =  'Burkina Faso'
countries <-  c('Pakistan', 'Mali','Ethiopia' , 'Burkina Faso', 'Benin', 'Ivory Coast', 'Cameroon', 'Togo', 
                'Nigeria',  'Kenya', 'Malawi', 'Mozambique', 'Zambia', 'Tunisia', 'India', 'Vietnam')

special_countries_cod <- getData("ISO3")%>%
  as.tibble() %>% 
  filter(NAME %in% countries)

probando = list()
for(i in 1:16){
  # =- 
  special_countries <-  special_countries_cod[i, ]
  
  # Descargas de los shp. 
  tictoc::tic()
  special_shps <-  do.call("rbind", lapply(special_countries$NAME, 
                                           function(x) getData('GADM', country=x, level=0))) %>% 
    sf::st_as_sf()
  tictoc::toc()
  
  # ggplot() + 
  #   geom_raster(data = test, aes(x = x, y = y, fill = id)) + 
  #   geom_sf(data = special_shps, color = 'red', fill = NA) +
  #   theme_bw()
  
  
  probando[[i]] <- mask(rasterDF, special_shps) %>% rasterToPoints() %>% 
    as.tibble() %>% mutate(ISO3 = special_countries$ISO3, Country = special_countries$NAME)
}




probando %>% bind_rows() %>% write_csv('id_country.csv')