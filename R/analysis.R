# goa dusky 2024 assessment converted to RTMB
# ben.williams
# 2025-01

# load ----
library(RTMB)
library(Matrix)
library(tidyverse)
library(here)
library(scico)
library(stringr)
# theme_set(afscassess::theme_report())

# functions 
source(here::here("R", "models.r"))
source(here::here("R", "utils.r"))

get_vec <- function(x, data, ind = 0, type = 1) {
  if(type!=1) {
    as.numeric(str_split(data[grep(x, data)], "\t")[[1]][1])
  } else {
  as.numeric(na.omit(as.numeric(str_split(data[grep(x, data)+ind], " ")[[1]])))
  }
}

# data ----
# to match up exactly with ADMB use values from the input/output
CTL = readLines(here::here('admb', 'goa_dusk.ctl'))
DAT = readLines(here::here('admb', 'goa_dusk_2024.dat'))
REP = readLines(here::here('admb', 'base.rep'))
PAR = readLines(here::here('admb', 'base.par'))

rec_age = get_vec("rec_age", DAT, 1)
years = get_vec("styr", DAT, 1):get_vec("endyr", DAT, 1)
ages = rec_age:(get_vec("nages_D", DAT, 1)+rec_age-1)
length_bins = get_vec("len_bin_labels", DAT, 1)
waa = get_vec("Weight", REP, 0)
maa = get_vec("Maturity", REP, 0)
wt_mature = maa * waa * 0.5
spawn_mo = 3


catch_obs = get_vec('obs_catch', DAT,ind=3)
catch_ind = rep(1, length(years))
catch_wt = ifelse(years<=1991, 2, 50)

srv_yrs = get_vec('Trawl survey years:', DAT, ind=1)
srv_ind = ifelse(years %in% srv_yrs, 1, 0)
srv_obs = get_vec('obs_srv1_biom', DAT, ind=1)
srv_sd = get_vec('obs_srv1_se', DAT, ind=1)


REP[grep("Obs_P_fish_age", REP):(grep("Pred_P_fish_age", REP)-2)] %>% 
  str_split(" ") %>%
  purrr::map_if(., is.character, as.numeric) %>% 
  purrr::map(~as.data.frame(t(.))) %>% 
  dplyr::bind_rows() %>% 
  dplyr::slice(-1) %>% 
  dplyr::select_if(~sum(!is.na(.))>0) -> obs

fish_age_yrs = obs$V1
fish_age_ind = ifelse(years%in% fish_age_yrs, 1, 0)
fish_age_iss = obs$V33
fish_age_obs = unname(t(as.matrix(dplyr::select(obs, V3:V29))))

REP[grep("Obs_P_srv1_age", REP):(grep("Pred_P_srv1_age", REP)-2)] %>% 
  str_split(" ") %>%
  purrr::map_if(., is.character, as.numeric) %>% 
  purrr::map(~as.data.frame(t(.))) %>% 
  dplyr::bind_rows() %>% 
  dplyr::slice(-1) %>% 
  dplyr::select_if(~sum(!is.na(.))>0) -> obs

srv_age_yrs = obs$V1
srv_age_ind = ifelse(years%in% srv_age_yrs, 1, 0)
srv_age_iss = obs$V33
srv_age_obs = unname(t(as.matrix(dplyr::select(obs, V3:V29))))

REP[grep("Obs_P_fish_size", REP):(grep("Pred_P_fish_size", REP)-2)] %>% 
  str_split(" ") %>%
  purrr::map_if(., is.character, as.numeric) %>% 
  purrr::map(~as.data.frame(t(.))) %>% 
  dplyr::bind_rows() %>% 
  dplyr::slice(-1) %>% 
  dplyr::select_if(~sum(!is.na(.))>0) -> obs

fish_size_yrs = obs$V1
fish_size_ind = ifelse(years%in% fish_size_yrs, 1, 0)
fish_size_iss = obs$V38
fish_size_obs = unname(t(as.matrix(dplyr::select(obs, V3:V34))))

DAT[(grep("Size-age transition matrix:", DAT)+2):(grep("age error transition matrix", DAT)-5)] %>% 
  str_split(" ") %>%
  purrr::map_if(., is.character, as.numeric) %>% 
  purrr::map(~as.data.frame(t(.))) %>% 
  dplyr::bind_rows() %>% 
  as.matrix() %>% 
  unname() -> saa

DAT[(grep("age error transition matrix", DAT)+2):(grep("end of file marker", DAT)-5)] %>% 
  str_split(" ") %>%
  purrr::map_if(., is.character, as.numeric) %>% 
  purrr::map(~as.data.frame(t(.))) %>% 
  dplyr::bind_rows() %>% 
  as.matrix() %>% 
  unname() -> ae

yield = get_vec("\\yield", CTL, type=2)
srv_wt = get_vec("wt_srv1", CTL, type=2)
fish_age_wt = get_vec("wt_fish_age", CTL, type=2)
fish_size_wt = get_vec("wt_fish_size", CTL, type=2)
srv_age_wt = get_vec("wt_srv1_age", CTL, type=2)
wt_fmort_reg = get_vec("wt_fmort_reg", CTL, type=2)
wt_rec_var = get_vec("wt_rec_var", CTL, type=2)

# pars 
log_M = get_vec("logm", PAR, ind=1)
log_a50C = log(get_vec("\\ba50:", PAR, ind=1))
deltaC = get_vec("\\bdelta:", PAR, ind=1)
log_a50S = log(get_vec("\\ba50_srv1", PAR, ind=1))
deltaS = get_vec("\\bdelta_srv1", PAR, ind=1)
log_q = get_vec("\\blog_q_srv1", PAR, ind=1)
log_mean_R = get_vec("\\blog_mean_rec", PAR, ind=1)
init_log_Rt = rev(get_vec('\\blog_rec_dev', PAR, ind=1)[1:(nrow(ae)-2)]) # ADMB brings in inital rec devs different, so flip them around
log_Rt = get_vec('\\blog_rec_dev', PAR, ind=1)[(nrow(ae)-1):(length(years)+(nrow(ae)-2))]
log_mean_F = get_vec("\\blog_avg_F", PAR, ind=1)
log_Ft = get_vec("\\blog_F_devs", PAR, ind=1)
log_F35 = log(get_vec("\\bmF35", PAR, ind=1))
log_F40= log(get_vec("\\bmF40", PAR, ind=1))
log_F50 = log(get_vec("\\bmF50", PAR, ind=1))
sigmaR = get_vec("\\bsigr", PAR, ind=1)

# data ----

data <- list(
  ages = ages,
  years = years,
  length_bins = length_bins,
  waa = waa,
  maa = maa,
  wt_mature = wt_mature,
  spawn_mo = spawn_mo,
  catch_ind = catch_ind,
  catch_obs = catch_obs,
  catch_wt = catch_wt,
  srv_obs = srv_obs,
  srv_ind = srv_ind,
  srv_sd = srv_sd,
  srv_wt = srv_wt,
  fish_age_obs = fish_age_obs,
  fish_age_ind = fish_age_ind,
  fish_age_iss = fish_age_iss,
  fish_age_wt = fish_age_wt,
  srv_age_obs = srv_age_obs,
  srv_age_ind = srv_age_ind,
  srv_age_iss = srv_age_iss,
  srv_age_wt = srv_age_wt,
  fish_size_obs = fish_size_obs,
  fish_size_ind = fish_size_ind,
  fish_size_iss = fish_size_iss,
  fish_size_wt = fish_size_wt,
  age_error = ae,
  size_age = saa,
  wt_fmort_reg = 2,
  wt_rec_var = 1,
  mean_M = 0.07,
  cv_M = 0.05,
  mean_q = 1,
  cv_q = 0.447213595,
  mean_sigmaR = 1.5,
  cv_sigmaR = 0.447213595,
  yield_ratio = yield
)
str(data)

# pars ----
pars = list(log_M = log_M,
            log_a50C = log_a50C,
            deltaC = deltaC,
            log_a50S = log_a50S,
            deltaS = deltaS,
            log_q = log_q,
            log_mean_R = log_mean_R,
            init_log_Rt = init_log_Rt,
            log_Rt = log_Rt,
            log_mean_F = log_mean_F,
            log_Ft =  log_Ft,
            log_F35 = log_F35,
            log_F40 = log_F40,
            log_F50 = log_F50,
            sigmaR = sigmaR)

map = list(log_M = factor(NA),
            log_a50C = factor(NA),
            deltaC = factor(NA),
            log_a50S = factor(NA),
            deltaS = factor(NA),
            log_q = factor(NA),
            log_mean_R = factor(NA),
            init_log_Rt = factor(rep(NA, length(init_log_Rt))),
            log_Rt = factor(rep(NA, length(log_Rt))),
            log_mean_F = factor(NA),
            log_Ft =  factor(rep(NA, length(log_Ft))),
            log_F35 = factor(NA),
            log_F40 = factor(NA),
            log_F50 = factor(NA),
            sigmaR = factor(NA))



# without running model ----
f(pars)
obj <- RTMB::MakeADFun(f,
                       pars,
                       map=map)
report <- obj$report(obj$env$last.par.best)
proj_bio(report)

# run model ----
pars = list(log_M = log(0.07),
            log_a50C = log(7),
            deltaC = 3,
            log_a50S = log(7),
            deltaS = 3,
            log_q = 0,
            log_mean_R = 3,
            init_log_Rt = rep(0, nrow(ae)-2),
            log_Rt = rep(0, length(years)),
            log_mean_F = -3,
            log_Ft =  rep(0, length(years)),
            log_F35 = 0,
            log_F40 = 0,
            log_F50 = 0,
            sigmaR = 1.5)
obj1 = RTMB::MakeADFun(f,
                       pars,
                       map=list(log_M = factor(NA)))
fit = nlminb(start=obj1$par,
             objective = obj1$fn,
             greadfient = obj1$gr,
             control = list(iter.max=100000,
                            eval.max=20000))
rep1 <- obj1$report(obj1$env$last.par.best)
proj_bio(rep1)

