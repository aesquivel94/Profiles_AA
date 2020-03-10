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

fls  <- list.files(paste0(root,'/data/Chirps_Chirts'), pattern = '*.fst$', full.names = T)

# Id - country
id_country <- read_csv('//dapadfs/workspace_cluster_8/climateriskprofiles/data/id_country.csv') 

id_country_f <- id_country %>% 
  filter(Country %in% c('Pakistan', 'Mali','Ethiopia'))

# Years
all_years <- 1985:2015


all_country <- id_country_f %>% 
  nest(-ISO3, -Country ) %>% 
  filter(row_number() == 3) %>% 
  dplyr::select(data) %>% 
  unnest()

# Temporalmente olvidar la temporada. 
# Para un solo ID

id <- all_country %>% filter(id == 95) %>% pull(id)


# Aqui se agrupa por year...
id_year <- tibble(id =  id, year = all_years) %>% 
  filter(year == 1985)



id_year 


# Aqui la funcion year - id 
year = 1985
id = 95
semester = 1


# Si solo se tiene un semestre...




#  Funcion interna para correr los semestres



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


# pixel_inf <- tibble(id = 95, year = 1986,  semester = 2)
reading_run_pys(pixel_inf = tibble(id = 95, year = 1986,  semester = 1), climate_path =  climate_path)

#


