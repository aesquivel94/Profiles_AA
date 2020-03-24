# =-------------------------------------------------------------------
# Climate Indexes Project
# Alejandra E. - Harold A. 
# Februrary 2020
# =-------------------------------------------------------------------

rm(list = ls())
gc(reset = TRUE)

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


options(warn = -1, scipen = 999)

suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, raster, vroom, sp, fst))
# =--------------------

# Ruta Principal para guardados: 
root <- '//dapadfs/workspace_cluster_8/climateriskprofiles/data/'


# =------ Datos procesados...

# Datos observacionales:
datos_Obs <- '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/data/observational_data'

# Datos de modelos de clima sin corrección de sesgo:
d_model_CSC_sesg <- '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/data/gcm_0_05deg_lat_county'

# Datos de modelos de clima con corrección de sesgo:
datos_cli_cs <- '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/data/bc_quantile_0_05deg_lat_county'

# =------------------------------------------------------------------------------

# =----------------------------------
# Identificacion de pixel para ETH
# =----------------------------------
country <- 'Ethiopia'
county  <- 'Arsi'
time <-  'past' # Historic - Future. 
country1 <- 'ethiopia'


# =----------------------------------
# Do folders. 
# =----------------------------------
save_folders <- '//dapadfs/workspace_cluster_8/climateriskprofiles/results/'

# if(dir.exists(output) == FALSE){
#   dir.create(output)
# }else{cat('Ok')}

# =----------------------------------



# Load county shapefile
shp <- raster::shapefile(paste0(root,'/shps/Ethiopia/ETH_adm2.shp'))
shp <- shp[shp@data$NAME_2 == county,]
plot(shp)

# Load id coords
crd <- vroom(paste0(root,'/id_country.csv'), delim = ',')
crd <- crd %>%
  dplyr::filter(Country == country)
pnt <- crd %>% dplyr::select(x,y) %>% sp::SpatialPoints(coords = .)
crs(pnt) <- crs(shp)
# Filter coordinates that are present in the county
pnt <- sp::over(pnt, shp) %>% data.frame %>% dplyr::select(ID_0:TYPE_2) %>% complete.cases() %>% which()
crd <- crd[pnt,]
crd <<- crd  %>% 
  mutate(county = county)


# =-------------------------------------------------------
# Nasa Download...
# =-------------------------------------------------------

Nasa_download <- function(id, coords){
  
  file_to_save <- glue::glue('//dapadfs/workspace_cluster_8/climateriskprofiles/data/NASA/{id}.fst')
  
  if(file.exists(file_to_save) == FALSE){
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
    
    fst::write.fst(x = var_NASA, path = file_to_save)
  }else{ var_NASA <- NULL}
  
  return(var_NASA)}

by_list_nasa <- function(list_n){
  
  x <- purrr::walk2(.x =  list_n$id, .y = list_n$coords, .f = Nasa_download)
  Sys.sleep(60)
  
  return(x)}

# =- 
  # crd %>% nest(-ISO3, -Country, -county, -id) %>% rename(coords = data) %>%
  #   filter(!(id %in% id_run)) %>% #filter(row_number() <= 195) %>%
  #   mutate(group =  c(rep(1:2, each = 90), rep(3, 15))) %>%
  #   group_split(group) %>% purrr::map(.f = by_list_nasa)

  

# =----------------------------------------------------------------------------
# =----------------------------------------------------------------------------
# =----------------------------------------------------------------------------

# ------------------------------------------------------------------------------------- #
# Index Functions. 
# ------------------------------------------------------------------------------------- #
rsum.lapply <- function(x, n=3L) # Calculate rollin sum
{
  lapply(1:(length(x)-n+1), function(i)
  {
    # Sum for n consecutive days
    z <- sum(x[i:(i+n-1)])
    # Indices used to calculate the sum
    seq.sum <- as.numeric(i:(i+n-1))
    # List with SUM and INDICES
    results <- list(z, seq.sum)
    return(results)
  })
}
cumulative.r.sum <- function(results){ unlist(lapply(results, function(x){z <- x[[1]]; return(z)})) } # Extract the SUM
is.leapyear <- function(year){ return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0)) } # Function to identify leap years


## CDD. Maximum number of consecutive dry days
calc_cdd <- function(PREC, p_thresh=1){
  runs <- rle(PREC < p_thresh)
  cons_days <- max(runs$lengths[runs$values==1], na.rm=TRUE)
  return(cons_days)
}
calc_cddCMP <- cmpfun(calc_cdd)

## P5D. Maximum 5-day running average precipitation
calc_p5d <- function(PREC){
  runAvg <- caTools::runmean(PREC, k=5, endrule='NA')
  runAvg <- max(runAvg, na.rm=TRUE)
  return(runAvg)
}
calc_p5dCMP <- cmpfun(calc_p5d)

## NT35. Number of days with max. temperature above 35?C
calc_hts <- function(TMAX, t_thresh=35) {
  hts <- length(which(TMAX >= t_thresh))
  return(hts)
}
calc_htsCMP <- cmpfun(calc_hts)

## P95. 95th percentile of daily precipitation
calc_p95 <- function(PREC){
  quantile(PREC, probs = .95, na.rm = T)
}
calc_p95CMP <- cmpfun(calc_p95)



# =----------------------------------------------------------
# function to run each semester. 
run_each_semester <- function(one, semester){
  
  one2 <- one %>%
    dplyr::select(id, year) %>%
    unique()
  
  
  # one2 <- one2 %>% dplyr::mutate(season = list(rsum.lapply(x = one$prec %>% as.numeric, n = 150)))
  
  one2 <- one2 %>%
    dplyr::mutate(season = purrr::map(.x =year, .f = function(y){
      oi = list(rsum.lapply(x = filter(one, year == y)$prec %>% as.numeric, n = 150))[[1]] } ) )
  
  # one2 <- mutate(one2, sum.ind = list(season[[1]][[which.max(cumulative.r.sum(season[[1]]))]]))
  
  one2 <- one2 %>%
    mutate(sum.ind = purrr::map(.x = season, function(y){list(y[[which.max(cumulative.r.sum(y[[1]]))]])[[1]]}))
  
  # data <-  one2 %>%
  #   mutate(semester = semester) %>%
  #   dplyr::mutate(CDD = calc_cddCMP(one$prec[one2$sum.ind[[1]][[2]]])) %>%
  #   dplyr::mutate(P5D = calc_p5dCMP(one$prec[one2$sum.ind[[1]][[2]]])) %>%
  #   dplyr::mutate(P95 = calc_p95CMP(one$prec[one2$sum.ind[[1]][[2]]])) %>%
  #   dplyr::mutate(NT35 = calc_htsCMP(one$tmax[one2$sum.ind[[1]][[2]]], t_thresh = 35)) %>%
  #   dplyr::select(-season) %>%
  #   dplyr::select(id, year, semester, everything())
  
  data <- one2 %>% 
    mutate(semester = semester) %>% 
    mutate(climate = purrr::map(.x = year, .f = function(x){filter(one , year == x)})) %>%
    mutate(index = purrr::map2(.x = sum.ind, .y = climate, .f = function(x, y){
      tibble(CDD = calc_cddCMP(y$prec[x[[2]]]), 
             P5D = calc_p5dCMP(y$prec[x[[2]]]), 
             P95 = calc_p95CMP(y$prec[x[[2]]]), 
             NT35 = calc_htsCMP(y$prec[x[[2]]], t_thresh = 35))
    })) %>% unnest(index)
  
  
  return(data)}
# =----------------------------------------------------------




# =----------------------------------------------------------
# Funciones para suelos... hay que tener en cuenta que... primero se debe correr 
# el balance hidrico... (creo...)







# =----------------------------------------------------------------------------
# =----------------------------------------------------------------------------
# =----------------------------------------------------------------------------

tictoc::tic()
observacional_data <- read_rds('//dapadfs/workspace_cluster_8/climateriskprofiles/data/observational_data/ethiopia/arsi.RDS')
tictoc::toc()


# Soil <- fst::read_fst(paste0('//dapadfs/workspace_cluster_8/climateriskprofiles/data/soil_data.fst'))
Soil1 <- fst::read.fst('//dapadfs/workspace_cluster_8/climateriskprofiles/data/soilcp_data.fst') %>% 
  as_tibble() %>% 
  dplyr::select(id, soilcp) %>% 
  filter(id %in%  pull(crd, id))

# =-
observacional_data <- observacional_data %>% 
  mutate(county = county, semester = 1) %>%
  mutate(id = as.integer(id))


observacional_data <- inner_join(observacional_data, Soil1) %>%
  # filter(is.na(soilcp)) %>% 
  mutate(soilcp = ifelse(is.na(soilcp), 100, soilcp))


# px <- observacional_data %>%
#   dplyr::group_split(id) %>%
#   .[[1]]

# reading_run_pys(px)
# Run... Lo correre como una lista... creo que asi sera mas facil...
  
# reading_run_pys(pixel_inf = px, climate_path = root)
reading_run_pys <- function( px){
  
  id <- px$id ; semester <- px$semester ; ISO3 <- px$ISO3 ; county <- px$county; soilcp <- px$soilcp
  
  if(semester == 2){
    # Si se tienen 2 se tiene un semestres...
    # semester = 2
  data_s <-  dplyr::select(px, Climate)  %>% unnest() %>% 
    mutate(year = lubridate::year(Date), month = lubridate::month(Date), 
                                                      semester = ifelse(month < 7, 1, 2)) %>% 
    dplyr::select(-month) 
  
    one <- data_s %>% 
      dplyr::select(id, semester,  year, prec, tmax, tmin) %>%
      nest(-semester)
    
    data_base <- purrr::map2(.x = one$data, .y = one$semester, .f = run_each_semester) %>% bind_rows()
    
  }else if(semester == 1){

    one  <- dplyr::select(px, Climate)  %>% unnest() %>% 
      mutate(year = lubridate::year(Date)) %>%
      dplyr::select(id, year, prec, tmax, tmin)


    data_base <- run_each_semester(one = one, semester =  1)

  }else{ data_base <- NULL}
  
  return(data_base)}


# Aqui vamos a correr los indices que tengan solo clima normal...

tictoc::tic()
plan(multiprocess)
options(future.globals.maxSize= 891289600)
tictoc::tic()

index_by_pixel <- observacional_data %>% # filter(row_number() < 20) %>%
  dplyr::group_split(id) %>%
  furrr::future_map(.f = reading_run_pys) %>% 
  bind_rows()

gc()
gc(reset = T)
tictoc::toc() # 8 min. 



# data_with_C_index <- observacional_data %>%inner_join(. , index_by_pixel)

# data_with_C_index <- dplyr::select(index_by_pixel,-climate)

output_P <- glue::glue('{save_folders}/{country}/{time}/')
readr::write_rds(x = index_by_pixel, path = glue::glue('{output_P}{county}_1985_2015.RDS'))  




# Mapas presente -----






# \\dapadfs.cgiarad.org\workspace_cluster_8\climateriskprofiles\data\bc_quantile_0_05deg_lat_county\ethiopia\ipsl_cm5a_mr\2021_2045
