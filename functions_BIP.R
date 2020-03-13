###  Present Basic functions. 
# By: H. Achicanoy & A. Esquivel.
# CIAT - March - 2020.

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
library(glue)

options(warn = -1, scipen = 999)
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, vroom, raster, fst, compiler, caTools))
# =--------------------

# =--------------------
# Packages 
root <- '//dapadfs/workspace_cluster_8/climateriskprofiles/data/'
climate_path <- glue::glue('{root}Chirps_Chirts/')


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


# ------------------------------------------------------------------------------------- #
# Index Functions. 
# ------------------------------------------------------------------------------------- #

# fls  <- list.files(paste0(root,'/data/Chirps_Chirts'), pattern = '*.fst$', full.names = T)

# Id - country
id_country <- read_csv('//dapadfs/workspace_cluster_8/climateriskprofiles/data/id_country.csv') 

id_country_f <- id_country %>% 
  filter(Country %in% c('Pakistan', 'Mali','Ethiopia'))

# Years
# all_years <- 1985:2015

# Toy trabajando aqui...
all_country <- id_country_f %>%
  nest(-ISO3, -Country ) %>% 
  mutate(semester = c(1, 1, 2)) # %>% 
# filter(row_number() == 3) %>% 
# dplyr::select(data) %>% 
# unnest()

# Temporalmente olvidar la temporada. 
# Para un solo ID

# =----------------------------------------------------------
# function to run each semester. 
run_each_semester <- function(one, semester){
  
  one2 <- one %>%
    dplyr::select(id, year) %>%
    unique()
  
  
  one2 <- one2 %>%
    dplyr::mutate(season = list(rsum.lapply(x = one$prec %>% as.numeric, n = 150)))
  
  
  one2 <- mutate(one2, sum.ind = list(season[[1]][[which.max(cumulative.r.sum(season[[1]]))]]))
  
  data <- one2 %>%
    mutate(semester = semester) %>%
    dplyr::mutate(CDD = calc_cddCMP(one$prec[one2$sum.ind[[1]][[2]]])) %>%
    dplyr::mutate(P5D = calc_p5dCMP(one$prec[one2$sum.ind[[1]][[2]]])) %>%
    dplyr::mutate(P95 = calc_p95CMP(one$prec[one2$sum.ind[[1]][[2]]])) %>%
    dplyr::mutate(NT35 = calc_htsCMP(one$tmax[one2$sum.ind[[1]][[2]]], t_thresh = 35)) %>%
    dplyr::select(-season) %>%
    dplyr::select(id, year, semester, everything())
  
  return(data)}
# =----------------------------------------------------------
# function to run for each id, year, semester. 
# # pixel_inf <- tibble(id = 95, year = 1986,  semester = 2)

reading_run_pys <- function( pixel_inf, climate_path){
  
  year <- pixel_inf$year; id <- pixel_inf$id ; semester <- pixel_inf$semester
  if(semester == 2){
    # Si se tienen 2 se tiene un semestres...
    # semester = 2
    data_s <- tibble(date = list.files(path = climate_path, pattern = as.character(year)) %>% str_remove('D_') %>% str_remove('.fst') %>%  lubridate::as_date(Date), 
                     files = list.files(path = climate_path, pattern = as.character(year), full.names = TRUE)) %>% 
      mutate(year = lubridate::year(date), month = lubridate::month(date), 
             semester = ifelse(month < 7, 1, 2)) %>% 
      dplyr::select(-month) %>% 
      mutate(data = purrr::map(.x = files, .f = function(x, id){ fst::read_fst(x, from = id, to = id)}, id = id)) %>% 
      dplyr::select(-files)
    
    
    one <- data_s %>% 
      unnest() %>% 
      dplyr::select(id, semester,  year, prec, tmax, tmin) %>%
      nest(-semester)
    
    
    data_base <- purrr::map2(.x = one$data, .y = one$semester, .f = run_each_semester) %>% bind_rows()
    
  }else if(semester == 1){
    
    # tictoc::tic()
    one  <- tibble(date = list.files(path = climate_path, pattern = as.character(year)) %>% str_remove('D_') %>% str_remove('.fst') %>%  lubridate::as_date(Date), 
                   files = list.files(path = climate_path, pattern = as.character(year), full.names = TRUE)) %>% 
      mutate(data = purrr::map(.x = files, .f = function(x, id){ fst::read_fst(x, from = id, to = id)}, id = id)) %>% 
      dplyr::select(-files) %>% 
      mutate(year = lubridate::year(date))
    # tictoc::toc() # 3.28
    
    
    one <- one %>% 
      unnest() %>% 
      dplyr::select(id, year, prec, tmax, tmin) 
    
    
    data_base <- run_each_semester(one = one, semester =  1)
    
  }else{ data_base <- NULL}
  
  return(data_base)}

# =----------------------------------------------------------
run_pixelY <- function(ISO3,  data, climate_path){
  
  id <- data$id; x <- data$x ; y <- data$y ; semester <- data$semester
  
  # Aqui se agrupa por year...
  # tictoc::tic() 
  run_years_ps <- tibble(semester = semester, id =  id, year = 1985:2015)  %>% 
    mutate(count = 1:nrow(.)) %>% 
    nest(-count) %>% 
    mutate(index = purrr::map(.x = data,  .f = reading_run_pys, climate_path = climate_path))
  # tictoc::toc()  # 25.51
  
  data_index <- run_years_ps %>% 
    dplyr::select(-data) %>% mutate(ISO3 = ISO3) %>% 
    unnest(index) %>% 
    inner_join(data) %>% 
    dplyr::select(ISO3, id, x, y, semester, year, CDD, P5D , P95, NT35,  sum.ind)
  
  
  sum.ind <- data_index %>% 
    dplyr::select(ISO3, id,x,y,semester, year, sum.ind) %>%
    unnest() %>% # filter(row_number() == 1) %>% dplyr::select(sum.ind) %>% unnest()
    mutate(rows = purrr::map(.x = sum.ind, .f = length)) %>% 
    unnest(rows) 
  
  
  Acum <- sum.ind %>% filter(rows == 1) %>% unnest( sum.ind) %>% dplyr::select(-rows)
  fst::write_fst(x = Acum , glue::glue('//dapadfs/workspace_cluster_8/climateriskprofiles/data/metadata/periodos/acum/Acum_{ISO3}_{id}.fst'))
  
  sum.ind %>% 
    filter(rows > 1) %>% 
    dplyr::select(-rows) %>%
    mutate(file = glue::glue('//dapadfs/workspace_cluster_8/climateriskprofiles/data/metadata/periodos/{ISO3}/{ISO3}_{id}_{year}_{semester}') ) %>% 
    unnest(sum.ind) %>% 
    nest(-file) %>% 
    mutate(save = purrr::map2(.x = data, .y = file, .f = fst::write_fst))
  
  
  index <- data_index %>% 
    dplyr::select(ISO3, id, x, y, semester, year, CDD, P5D , P95, NT35)
  fst::write_fst(x = index , glue::glue('//dapadfs/workspace_cluster_8/climateriskprofiles/data/metadata/periodos/Index_climate/I_{ISO3}_{id}.fst'))
  
  return(data_index)}

# id <- all_country %>% filter(id == 95) %>% pull(id)
# data_country <- all_country %>% filter(row_number() == 3) %>% dplyr::select(data) %>% unnest() %>% 
#   dplyr::filter(id == 95)

# =------------------ Funcion para correr todos los pixel de un pais...
# Necesito tener un objeto que se llame ISO3


# mod <- all_country  %>% 
#   mutate(data_mod = purrr::map(.x = data, .f = function(x){filter(x, row_number() < 3)})) %>% 
#   dplyr::select(ISO3, Country, semester, data_mod)

# ISO3 <- 'MLI'
# country_folder <- glue::glue('//dapadfs/workspace_cluster_8/climateriskprofiles/data/metadata/periodos/{ISO3}/')
# ISO3
# if(dir.exists(country_folder) == FALSE){dir.create(country_folder)}else{cowsay::say('ok', 'smallcat')}




# tictoc::tic()
# plan(multiprocess)
# options(future.globals.maxSize= 891289600)
# # tictoc::tic()
# all_pixles <-  all_country %>%
#   # filter(row_number() == 1) %>%
#   # dplyr::select(ISO3, data) %>%
#   unnest() %>%
#   mutate(rows = 1:nrow(.)) %>%
#   nest(- ISO3, -Country, -rows) %>%
#   # mutate(run_by_id = purrr::map2(.x = ISO3, .y = data, .f = run_pixelY, climate_path = climate_path))# mutate(run_by_id = purrr::map2(.x = ISO3, .y = data, .f = run_pixelY, climate_path = climate_path))
#   mutate(run_by_id = furrr::future_map2(.x = ISO3, .y = data, .f = run_pixelY, climate_path = climate_path))
# gc()
# gc(reset = T)
# tictoc::toc() # 6.57



# probando <- all_country %>%
#   # filter(row_number() == 1) %>%
#   # dplyr::select(ISO3, data) %>%
#   unnest() %>%
#   mutate(rows = 1:nrow(.)) %>%
#   nest(- ISO3, -Country, -rows) %>% 
#   filter(row_number()==1)


# probando$data[[1]]
# probando$ISO3
# 
# run_pixelY(ISO3 = probando$ISO3, data = probando$data[[1]], climate_path = climate_path)




# no_cores <- availableCores() - 1
plan(cluster, workers = 30)

all_pixles <-  all_country %>%
  # filter(row_number() == 1) %>%
  # dplyr::select(ISO3, data) %>%
  unnest() %>%
  mutate(rows = 1:nrow(.)) %>%
  nest(- ISO3, -Country, -rows) %>%
  # mutate(run_by_id = purrr::map2(.x = ISO3, .y = data, .f = run_pixelY, climate_path = climate_path))# mutate(run_by_id = purrr::map2(.x = ISO3, .y = data, .f = run_pixelY, climate_path = climate_path))
  mutate(run_by_id = furrr::future_map2(.x = ISO3, .y = data, .f = run_pixelY, climate_path = climate_path))

