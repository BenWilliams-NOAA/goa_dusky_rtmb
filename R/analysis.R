# goa dusky 2024 assessment converted to RTMB
# ben.williams
# 2025-01

# load ---
source(here::here("R", "models.r"))
source(here::here("R", "utils.r"))

# data ----
# to match up exactly with ADMB use values from the input/output
CTL = readLines(here::here('admb', 'goa_dusk.ctl'))
DAT = readLines(here::here('admb', 'goa_dusk_2024.dat'))
REP = readLines(here::here('admb', 'base.rep'))
PAR = readLines(here::here('admb', 'base.par'))

rec_age = get_vec("rec_age", DAT, 1)
years = get_vec("styr", DAT, 1):get_vec("endyr", DAT, 1)
ages = rec_age:(get_vec("nages_D", DAT, 1)+rec_age-1)
length_bins = as.integer(get_vec("len_bin_labels", DAT, 1))
waa = get_vec("Weight", REP, 0)
maa = get_vec("Maturity", REP, 0)
wt_mature = maa * waa * 0.5
spawn_mo = 3

catch_obs = get_vec('obs_catch', DAT,ind=3)
catch_ind = as.integer(rep(1, length(years)))
catch_wt = ifelse(years<=1991, 2, 50)

srv_yrs = as.integer(get_vec('Trawl survey years:', DAT, ind=1))
srv_ind = as.integer(ifelse(years %in% srv_yrs, 1, 0))
srv_obs = get_vec('obs_srv1_biom', DAT, ind=1)
srv_sd = get_vec('obs_srv1_se', DAT, ind=1)

REP[grep("Obs_P_fish_age", REP):(grep("Pred_P_fish_age", REP)-2)] %>% 
  str_split(" ") %>%
  purrr::map_if(., is.character, as.numeric) %>% 
  purrr::map(~as.data.frame(t(.))) %>% 
  dplyr::bind_rows() %>% 
  dplyr::slice(-1) %>% 
  dplyr::select_if(~sum(!is.na(.))>0) -> obs

fish_age_yrs = as.integer(obs$V1)
fish_age_ind = as.integer(ifelse(years%in% fish_age_yrs, 1, 0))
fish_age_iss = obs$V33
fish_age_obs = unname(t(as.matrix(dplyr::select(obs, V3:V29))))

REP[grep("Obs_P_srv1_age", REP):(grep("Pred_P_srv1_age", REP)-2)] %>% 
  str_split(" ") %>%
  purrr::map_if(., is.character, as.numeric) %>% 
  purrr::map(~as.data.frame(t(.))) %>% 
  dplyr::bind_rows() %>% 
  dplyr::slice(-1) %>% 
  dplyr::select_if(~sum(!is.na(.))>0) -> obs

srv_age_yrs = as.integer(obs$V1)
srv_age_ind = as.integer(ifelse(years%in% srv_age_yrs, 1, 0))
srv_age_iss = obs$V33
srv_age_obs = unname(t(as.matrix(dplyr::select(obs, V3:V29))))

REP[grep("Obs_P_fish_size", REP):(grep("Pred_P_fish_size", REP)-2)] %>% 
  str_split(" ") %>%
  purrr::map_if(., is.character, as.numeric) %>% 
  purrr::map(~as.data.frame(t(.))) %>% 
  dplyr::bind_rows() %>% 
  dplyr::slice(-1) %>% 
  dplyr::select_if(~sum(!is.na(.))>0) -> obs

fish_size_yrs = as.integer(obs$V1)
fish_size_ind = as.integer(ifelse(years%in% fish_size_yrs, 1, 0))
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

# data list ----
# data and parameters are input as lists
# i've renamed a number of inputs and added the weights to the data as well (instead of a .ctl file)
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
  srv_yrs = srv_yrs,
  srv_obs = srv_obs,
  srv_ind = srv_ind,
  srv_sd = srv_sd,
  srv_wt = srv_wt,
  fish_age_yrs = fish_age_yrs,
  fish_age_obs = fish_age_obs,
  fish_age_ind = fish_age_ind,
  fish_age_iss = fish_age_iss,
  fish_age_wt = fish_age_wt,
  srv_age_yrs = srv_age_yrs,
  srv_age_obs = srv_age_obs,
  srv_age_ind = srv_age_ind,
  srv_age_iss = srv_age_iss,
  srv_age_wt = srv_age_wt,
  fish_size_yrs = fish_size_yrs,
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
# make sure indexes are integers and matrices are unnamed
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
str(pars)
# to hold a parameter fixed it gets "mapped"
# to compare the assessments we can first fix all parameters and check results without optimizing the model
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
f(pars) # check to see that the model runs (should be a reasonable nll output)
compiler::enableJIT(0)
# build the model
obj <- RTMB::MakeADFun(f,
                       pars,
                       map=map)
# get the output from the model
report <- obj$report(obj$env$last.par.best)
# since projections cannot be auto differentiated they are moved to an external function (found in utils.r)
proj_bio(report)

# run model ----
# optime the model using the same starting values as found in ADMB
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
# this is the optimization
fit = nlminb(start=obj1$par,
             objective = obj1$fn,
             gradient = obj1$gr,
             control = list(iter.max=100000,
                            eval.max=20000))
rep1 <- obj1$report(obj1$env$last.par.best)

# admb
data.frame(year = 2025:2026,
           spawn_bio = c(get_vec("Female_Spawning Biomass for 2025", REP, 1),
                         get_vec("Female_Spawning_Biomass for 2026", REP, 1)),
           total_bio = c(get_vec("TotalBiomass for 2025", REP, 1),
                         get_vec("TotalBiomass for 2026", REP, 1)),
           catch_abc = c(get_vec("\\bABC for 2025", REP, 1),
                         get_vec("\\bABC for 2026", REP, 1)),
           catch_ofl = c(get_vec("\\bOFL for 2025", REP, 1),
                         get_vec("\\bOFL for 2026", REP, 1)),
           F40 = c(get_vec("F_ABC for 2025", REP, 1),
                   get_vec("F_ABC for 2026", REP, 1)),
           F35 = c(get_vec("F_OFL for 2025", REP, 1),
                   get_vec("F_OFL for 2026", REP, 1)))
# bridge
proj_bio(report)
# optimized
proj_bio(rep1)

# some slight round error differences, but essentially the same - can doublecheck all the loglikelihood values etc.
# ssc noted that they appreciated tables 10-23:10-26 in the northern rockfish SAFE appendix https://www.npfmc.org/wp-content/PDFdocuments/SAFE/2024/GOAnork.pdf

# review some of the outputs
# sd of key outputs ----
sdrep = sdreport(obj1, getJointPrecision=TRUE)
summary(sdrep, "report") %>% 
  as.data.frame() %>% 
  mutate(lci = Estimate - 1.96 * `Std. Error`,
         uci = Estimate + 1.96 * `Std. Error`) %>% 
  tibble::rownames_to_column('item') %>% 
  filter(!(item %in% c('q', 'M'))) %>% 
  mutate(item = gsub('\\..*', '', item),
         year = c(data$years, data$srv_yrs, rep(data$years,3))) -> m24_se

## mcmc ----
Q <- sdrep$jointPrecision
# random effects or no?
if(!is.null(Q)) {
  M <- solve(Q)
} else {
  M <- sdrep$cov.fixed
  Q <- solve(M)
}
globals <- list(data = data)
chains <- 5

# Cole recommends 5 chains, 1000 iters, 250 warmup
mcmc <- sample_sparse_tmb(obj1, iter = 1000,
                          warmup = 400, 
                          chains = 1, cores = 1,
                          metric='dense', Qinv=M, Q=Q,
                          globals = globals, skip_optimization=TRUE)

post = adnuts::extract_samples(mcmc)
mpost = as.matrix(post)
r = list()
# process posterior
for(i in 1:nrow(mpost)) {
  r[[i]] = obj1$report(mpost[i,])
}

rep_out <- function(rep, data, item = 'spawn_bio') {
  purrr::map(rep, item) %>% 
    purrr::map(., ~as.data.frame(.)) %>% 
    bind_rows(.id = 'sim') %>% 
    mutate(year = rep(data$years, length(rep)))
}

purrr::map(r, proj_bio) %>% 
  purrr::map(., ~as.data.frame(.)) %>% 
  bind_rows(.id = 'sim') -> projs

tots = rep_out(r, data, 'tot_bio')  
ssb = rep_out(r, data, 'spawn_bio')
Ft = rep_out(r, data, 'Ft')
recs = rep_out(r, data, 'recruits')  

plot_par(item ='log_a50S', post=post, rep=rep1, rep_item='a50S')
plot_par(item ='log_q', post=post, rep=rep1, rep_item='q')


# francis rewt
weights = francis_rewt(data, model=f, pars, map=list(log_M=factor(NA)))

# update data with new weights
data$fish_age_wt = last(weights$fac)
data$srv_age_wt = last(weights$sac)
data$fish_size_wt = last(weights$fsc)

# rerun model with new weights
obj2 = RTMB::MakeADFun(f,
                       pars,
                       map=list(log_M = factor(NA)))
# this is the optimization
fit2 = nlminb(start=obj2$par,
             objective = obj2$fn,
             gradient = obj2$gr,
             control = list(iter.max=100000,
                            eval.max=20000))
rep2 <- obj2$report(obj2$env$last.par.best)
proj_bio(rep2)

# plot catch
data.frame(year = data$years,
           obs = data$catch_obs,
           rep1 = rep1$catch_pred,
           rep2 = rep2$catch_pred) %>% 
  tidyr::pivot_longer(-year) %>% 
  ggplot(aes(year, value, color=name, shape=name)) + 
  geom_point() + 
  geom_line() +
  theme_minimal() +
  scico::scale_color_scico_d(palette = 'roma')



# plot survey
data.frame(year = data$srv_yrs,
           obs = data$srv_obs,
           sd = data$srv_sd,
           rep1 = rep1$srv_pred,
           rep2 = rep2$srv_pred) %>% 
  mutate(lci = obs - 1.96 * sd,
         uci = obs + 1.96 * sd) %>% 
  tidyr::pivot_longer(-c(year, sd, obs, lci, uci)) %>% 
  ggplot(aes(year, value, color=name, shape=name)) + 
  geom_point(aes(y=obs), color='gray') + 
  geom_errorbar(aes(ymin=lci, ymax=uci), color='gray') + 
  geom_line() +
  theme_minimal() +
  scico::scale_color_scico_d(palette = 'roma') +
  expand_limits(y=0)

# comp ----

plot_comp <- function(data, rep, type, plot=TRUE) {
  obs = as.data.frame(data[paste0(type, '_obs')])
  pred = as.data.frame(rep[paste0(type, '_pred')])
  names(pred) = names(obs) = data[paste0(type, '_yrs')][[1]]
  id = deparse(substitute(rep))
  if(length(grep('age', type))==0) {
    ind = data$length_bins
    x = "Length bins"
  } else {
    ind = data$ages
    x = "Ages"
  }
  obs %>% 
    tidytable::mutate(index = ind,
                      id = 'obs') %>% 
    tidytable::bind_rows(
      pred %>% 
        tidytable::mutate(index = ind, 
                          id = id)) %>% 
        tidytable::pivot_longer(-c(index, id)) -> dat
    
  if(isTRUE(plot)){
    dat %>% 
  ggplot(aes(index, value, color=id)) + 
    geom_col(data = . %>% filter(id=='obs'), color='gray', fill='lightgray') + 
    geom_line(data = . %>% filter(id!='obs')) + 
    facet_wrap(~name, dir='v') +
    theme_minimal() +
      xlab(x)
  } else {
  dat
  }
}

plot_comp(data=data, rep=rep1, type="fish_age")
plot_comp(data=data, rep=rep1, type="fish_size")
plot_comp(data=data, rep=rep1, type="srv_age")

# compare rep1 & rep2
plot_comp(data=data, rep=rep1, type="srv_age", plot=F) %>% 
  bind_rows(
plot_comp(data=data, rep=rep2, type="srv_age", plot=F) 
) %>% 
  ggplot(aes(index, value, color=id)) + 
  geom_col(data = . %>% filter(id=='obs'), color='gray', fill='lightgray', alpha = 0.4) + 
  geom_line(data = . %>% filter(id!='obs')) + 
  facet_wrap(~name, dir='v') +
  scico::scale_color_scico_d(palette = 'roma')


rep1$nll
rep2$nll

rep1$ssqcatch
rep2$ssqcatch

rep1$like_srv
rep2$like_srv

rep1$like_rec
rep2$like_srv

rep1$like_fish_age
rep2$like_fish_age

rep1$like_srv_age
rep2$like_srv_age

rep1$like_fish_size
rep2$like_fish_size
