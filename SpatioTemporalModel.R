library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

##############################################################
#### Thermal sensitivity in the phenological distirbution ####
##############################################################

## Model code
write("data{
      int<lower=0> N;   // length of data
      int<lower=0> y[N];     // cater abundance for each sample
      vector[N] date;   // ordinal date of sample
      vector[N] temp_m;   // site-year mean temp for window predicting mu
      vector[N] temp_ls;   // site-year mean temp for window predicting logsigma
      vector[N] temp_lm;   // site-year mean temp for window predicting logmax
      int<lower=0> N_site;   // no. of sites
      int<lower=0> N_year;   // no. of years
      int<lower=0> N_siteyear;   // no. of site by year combinations 
      int<lower=0> N_siteday;   // no. of day by site by year combinations
      int<lower=0> N_treeID;   // no. of unique tree identities
      int<lower=0> N_rec;   // no. of recorders
      int<lower=0> N_resid;   // no. of residuals (N)
      int<lower=0,upper=N_site> site_id[N];    // site id for each sample
      int<lower=0,upper=N_year> year_id[N];   // year id for each sample
      int<lower=0,upper=N_siteyear> siteyear_id[N];    // site-year id for each sample
      int<lower=0,upper=N_siteday> siteday_id[N];    //day-site-year id for each sample
      int<lower=0,upper=N_treeID> treeID_id[N];   // tree identity id for each sample
      int<lower=0,upper=N_rec> rec_id[N];   // recorder id for each sample
      int<lower=0,upper=N_resid> resid_id[N];   // residual id for each sample
      }
      
      parameters{
      real mu;   // mean timing intercept
      real logsigma;   // log width intercept
      real logmax;   // log height intercept
      real t_mu;   // slope for change in mean timing with temp
      real t_logsigma;   // slope for change in log width with temp
      real t_logmax;   // slope for change in log height with temp
      cholesky_factor_corr[3] L;   // cholesky decomposition L matrix for gaussian parameters across random effects (used for s, y and sy)
      matrix[3, N_site] site_scaled;   // scaled site random effects, row for each gaussian parameter
      matrix[3, N_year] year_scaled;    // scaled year random effects
      matrix[3, N_siteyear] siteyear_scaled;   // scaled site-year random effects
      vector[N_siteday] siteday_scaled;   // scaled day-site-year random effects
      vector[N_treeID] treeID_scaled;   // scaled tree ID random effects
      vector[N_rec] rec_scaled;     // scaled recorder random effects  
      vector[N_resid] resid_scaled;     // scaled residual random effects 
      vector<lower=0>[3] sd_site;     // standard deviations of cater site random effects for each gaussian parameter
      vector<lower=0>[3] sd_year;     // standard deviations of cater year random effects for each gaussian parameter
      vector<lower=0>[3] sd_siteyear;     // standard deviations of cater site-year random effects for each gaussian parameter
      real<lower=0> sd_siteday;     // standard deviation of day-site-year random effects 
      real<lower=0> sd_treeID;     // standard deviation of tree ID random effects 
      real<lower=0> sd_rec;     // standard deviation of recorder random effects 
      real<lower=0> sd_resid;     // standard deviation of residual random effects 
      }
      
      transformed parameters{
      matrix[N_site,3] site_effs;     // cater unscaled site random effects, column for each parameter
      matrix[N_year,3] year_effs;     // cater unscaled year random effects
      matrix[N_siteyear,3] siteyear_effs;     // cater unscaled site-year random effects
      
      site_effs = (diag_pre_multiply(sd_site, L) * site_scaled)';  // sds as a diag matrix multipled by L matrix then by scaled effects and transposed
      year_effs = (diag_pre_multiply(sd_year, L) * year_scaled)';
      siteyear_effs = (diag_pre_multiply(sd_siteyear, L) * siteyear_scaled)';
      }
      
      model{
      vector[N] y_1; // logscale cater estimates
      vector[N_siteday] siteday_effects = sd_siteday * siteday_scaled;   // unscaled day-site-year random effects
      vector[N_treeID] treeID_effects = sd_treeID * treeID_scaled;   // unscaled tree ID random effects 
      vector[N_rec] rec_effects = sd_rec * rec_scaled;   // unscaled recorder random effects 
      vector[N_resid] resid_effects = sd_resid * resid_scaled;   // unscaled residual random effects
      
      y_1 = logmax + t_logmax*temp_lm - 
      ((date-(mu+t_mu*temp_m+site_effs[site_id,1]+year_effs[year_id,1]+siteyear_effs[siteyear_id,1])) .* 
      (date-(mu+t_mu*temp_m+site_effs[site_id,1]+year_effs[year_id,1]+siteyear_effs[siteyear_id,1]))) ./ 
      (2.*(exp(logsigma+t_logsigma*temp_ls+site_effs[site_id,2]+year_effs[year_id,2]+siteyear_effs[siteyear_id,2])).*
      (exp(logsigma+t_logsigma*temp_ls+site_effs[site_id,2]+year_effs[year_id,2]+siteyear_effs[siteyear_id,2]))) 
      + site_effs[site_id,3] + year_effs[year_id,3] 
      + siteyear_effs[siteyear_id,3] + siteday_effects[siteday_id] 
      + treeID_effects[treeID_id] + rec_effects[rec_id] + resid_effects[resid_id];   // gaussian function with temp slopes and random terms
      
      y ~ poisson_log(y_1);
      
      // Priors
      to_vector(site_scaled) ~ normal(0,1); 
      to_vector(year_scaled) ~ normal(0,1); 
      to_vector(siteyear_scaled) ~ normal(0,1); 
      siteday_scaled ~ normal(0,1); 
      treeID_scaled ~ normal(0,1); 
      rec_scaled ~ normal(0,1); 
      resid_scaled ~ normal(0,1);
      mu ~ normal(0,10);
      logsigma ~ normal(0,10);
      logmax ~ normal(0,10);
      t_mu ~ normal(0,10);           
      t_logsigma ~ normal(0,10);     
      t_logmax ~ normal(0,10);       
      L ~ lkj_corr_cholesky(2.0);  
      to_vector(sd_site) ~ cauchy(0,10);
      to_vector(sd_year) ~ cauchy(0,10);
      to_vector(sd_siteyear) ~ cauchy(0,10);
      sd_siteday ~ cauchy(0,10);
      sd_treeID ~ cauchy(0,10);
      sd_rec ~ cauchy(0,10);
      sd_resid ~ cauchy(0,10);
      }    
      
      generated quantities {
      matrix[3, 3] omega;  // correlation matrix = L times L transposed
      omega = L * L';
      }
      ",
      
      "TempCater.stan"
)
stanc("TempCater.stan")
TempCaterMod <- stan_model("TempCater.stan")


##Stan data list
data_cater <- read.csv("~/data_cater.csv")

stan_data_Temp <-list(
  N=nrow(data_cater),
  date=data_cater$datescaled,
  temp_m=data_cater$temp_m.cent,
  temp_ls=data_cater$temp_ls.cent,
  temp_lm=data_cater$temp_lm.cent,
  N_site=length(unique(data_cater$site)),
  N_siteyear=length(unique(data_cater$siteyear)),
  N_year=length(unique(data_cater$year)),
  N_siteday=length(unique(data_cater$siteday)),
  N_treeID=length(unique(data_cater$treeID)),
  N_rec=length(unique(data_cater$recorder)),
  N_resid=length(unique(data_cater$resid)),
  y=data_cater$caterpillars,
  site_id=as.numeric(as.factor(data_cater$site)),
  siteyear_id=as.numeric(as.factor(data_cater$siteyear)),
  year_id=as.numeric(as.factor(data_cater$year)),
  siteday_id=as.numeric(as.factor(data_cater$siteday)),
  treeID_id=as.numeric(as.factor(data_cater$treeID)),
  rec_id=as.numeric(as.factor(data_cater$recorder)),
  resid_id=as.numeric(as.factor(data_cater$resid))
)


## Sampling
TempModel <- sampling(object=TempCaterMod, data=stan_data_Temp, 
                            iter=4500, warmup=2000, chains=4, thin=5, cores=4,
                            control = list(adapt_delta = 0.94)) 