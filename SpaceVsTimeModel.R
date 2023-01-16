library(lme4)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(rstantools)

####################################################
#### Spatial vs temporal thermal sensitivity in ####
####        the phenological distribution       ####
####################################################

#### Analysis ####

# Model code
write("data{
      int<lower=0> N;
      vector[N] date;
      vector[N] temp_m_s;           //
      vector[N] temp_ls_s;           //
      vector[N] temp_lm_s;           //
      vector[N] temp_m_y;           //
      vector[N] temp_ls_y;           //
      vector[N] temp_lm_y;           //
      int<lower=0> N_site;
      int<lower=0> N_year;
      int<lower=0> N_siteyear;
      int<lower=0> N_siteday;
      int<lower=0> N_treeID;
      int<lower=0> N_rec;
      int<lower=0> N_resid;
      int<lower=0> y[N];
      int<lower=0,upper=N_site> site_id[N]; 
      int<lower=0,upper=N_year> year_id[N]; 
      int<lower=0,upper=N_siteyear> siteyear_id[N]; 
      int<lower=0,upper=N_siteday> siteday_id[N]; 
      int<lower=0,upper=N_treeID> treeID_id[N]; 
      int<lower=0,upper=N_rec> rec_id[N]; 
      int<lower=0,upper=N_resid> resid_id[N];
      }
      
      parameters{
      real mu;
      real logsigma;
      real logmax;
      real t_mu_s;            //
      real t_logsigma_s;      //
      real t_logmax_s;        //
      real t_mu_y;            //
      real t_logsigma_y;      //
      real t_logmax_y;        //
      cholesky_factor_corr[3] L;
      matrix[3, N_site] site_scaled; 
      matrix[3, N_year] year_scaled; 
      matrix[3, N_siteyear] siteyear_scaled;
      vector[N_siteday] siteday_scaled; 
      vector[N_treeID] treeID_scaled; 
      vector[N_rec] rec_scaled; 
      vector[N_resid] resid_scaled;
      vector<lower=0>[3] sd_site;
      vector<lower=0>[3] sd_year; 
      vector<lower=0>[3] sd_siteyear;
      real<lower=0> sd_siteday;
      real<lower=0> sd_treeID;
      real<lower=0> sd_rec;
      real<lower=0> sd_resid;
      }
      
      transformed parameters{
      matrix[N_site,3] site_effs;
      matrix[N_year,3] year_effs;
      matrix[N_siteyear,3] siteyear_effs;  
      
      site_effs = (diag_pre_multiply(sd_site, L) * site_scaled)'; 
      year_effs = (diag_pre_multiply(sd_year, L) * year_scaled)';
      siteyear_effs = (diag_pre_multiply(sd_siteyear, L) * siteyear_scaled)';
      }
      
      model{
      vector[N] y_1;
      vector[N_siteday] siteday_effects = sd_siteday * siteday_scaled; 
      vector[N_treeID] treeID_effects = sd_treeID * treeID_scaled; 
      vector[N_rec] rec_effects = sd_rec * rec_scaled; 
      vector[N_resid] resid_effects = sd_resid * resid_scaled;
      
      y_1 = logmax + t_logmax_s*temp_lm_s + t_logmax_y*temp_lm_y - 
      ((date-(mu+t_mu_s*temp_m_s+t_mu_y*temp_m_y+site_effs[site_id,1]+year_effs[year_id,1]+siteyear_effs[siteyear_id,1])) .* 
      (date-(mu+t_mu_s*temp_m_s+t_mu_y*temp_m_y+site_effs[site_id,1]+year_effs[year_id,1]+siteyear_effs[siteyear_id,1]))) ./ 
      (2.*(exp(logsigma+t_logsigma_s*temp_ls_s+t_logsigma_y*temp_ls_y+site_effs[site_id,2]+year_effs[year_id,2]+siteyear_effs[siteyear_id,2])).*
      (exp(logsigma+t_logsigma_s*temp_ls_s+t_logsigma_y*temp_ls_y+site_effs[site_id,2]+year_effs[year_id,2]+siteyear_effs[siteyear_id,2]))) 
      + site_effs[site_id,3] + year_effs[year_id,3] 
      + siteyear_effs[siteyear_id,3] + siteday_effects[siteday_id] 
      + treeID_effects[treeID_id] + rec_effects[rec_id] + resid_effects[resid_id];
      
      
      y ~ poisson_log(y_1);
      
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
      t_mu_s ~ normal(0,10);           //
      t_logsigma_s ~ normal(0,10);     //
      t_logmax_s ~ normal(0,10);       //
      t_mu_y ~ normal(0,10);           //
      t_logsigma_y ~ normal(0,10);     //
      t_logmax_y ~ normal(0,10);       //
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
      matrix[3, 3] omega;
      omega = L * L'; 
      }
      ",
      
      "SpaceVsTimeCater.stan"
)
stanc("SpaceVsTimeCater.stan")
SpaceVsTimeCaterMod <- stan_model("SpaceVsTimeCater.stan")

# Load data
data_cater <- read.csv("~/data_cater.csv")
temps <- read.csv("~/sy_temp_by_para.csv")

# Site and year for temperatures
temps$site <- substr(temps$siteyear, 1,3)
temps$year <- substr(temps$siteyear,5,8)
temps <- temps[-which(temps$site=="IVN"),] # site with no caterpillar data

# Model site mean temperature for each phenological parameter window
mod.mu <- lmer(temp.mu~1+(1|site)+(1|year), data=temps)
mod.ls <- lmer(temp.logsigma~1+(1|site)+(1|year), data=temps)
mod.lm <- lmer(temp.logmax~1+(1|site)+(1|year), data=temps)

plot(mod.mu)
qqnorm(resid(mod.mu))
qqline(resid(mod.mu))
plot(mod.ls)
qqnorm(resid(mod.ls))
qqline(resid(mod.ls))
plot(mod.lm)
qqnorm(resid(mod.lm))
qqline(resid(mod.lm))

temps.s <- data.frame(site=rownames(coef(mod.mu)$site),
                      mod_mu=coef(mod.mu)$site[,1],
                      mod_ls=coef(mod.ls)$site[,1],
                      mod_lm=coef(mod.lm)$site[,1])

#Bring site means into temp data and calculate annual deviations
temps <- merge(temps, temps.s, by="site")
temps$temp_m_s.cent <- temps$mod_mu-mean(temps$mod_mu)
temps$temp_ls_s.cent <- temps$mod_ls-mean(temps$mod_ls)
temps$temp_lm_s.cent <- temps$mod_lm-mean(temps$mod_lm)
temps$temp_m_y <- temps$temp.mu-temps$mod_mu
temps$temp_ls_y <- temps$temp.logsigma-temps$mod_ls
temps$temp_lm_y <- temps$temp.logmax-temps$mod_lm
temps$site <- NULL
temps$year <- NULL
data_cater <- merge(data_cater, temps, by="siteyear")

# Stan data list
data_cater <- merge(data_cater, temps, by="siteyear")

stan_data_SvT <-list(
  N=nrow(data_cater),
  date=data_cater$datescaled,
  temp_m_s=data_cater$temp_m_s.cent,
  temp_ls_s=data_cater$temp_ls_s.cent,
  temp_lm_s=data_cater$temp_lm_s.cent,
  temp_m_y=data_cater$temp_m_y,
  temp_ls_y=data_cater$temp_ls_y,
  temp_lm_y=data_cater$temp_lm_y,
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

# Sampling
SpaceTimeModel <- sampling(object=SpaceVsTimeCaterMod, data=stan_data_SvT, 
                           iter=4500, warmup=2000, chains=4, thin=5, cores=2,
                           control = list(adapt_delta = 0.94))

## Results ##

#Extract the posterior
stanpost <- function(model){
  
  df <- as.data.frame(extract(model, permuted=FALSE))
  chns <- model@sim$chains
  npara <- (ncol(df))/chns
  
  # Organise + combine chains of posterior distributions
  remove.start <- function(x, n){  #function to remove first n number of characters in a string
    substr(x, nchar(x)-(nchar(x)-n-1), nchar(x))
  }
  
  fullpost <- data.frame(matrix(NA, nrow = chns*nrow(df), ncol = npara))
  
  for(i in 1:npara){
    fullpost[,i] <- stack(df[,(i*chns-(chns-1)):(i*chns)])[,1]
    colnames(fullpost)[i] <- remove.start(colnames(df)[(i*chns-(chns-1))],8)
  }
  
  return(fullpost)
}  #convert multi-chained posteriors into one column per parameter   FUNCTION
STpost <- stanpost(STTempMod)

#columns of interest
STcoefs <- STpost[,1:9]

STint <- data.frame(coef=colnames(STcoefs)[1:3],
                    mean=colMeans(STcoefs)[1:3],
                    lowci=posterior_interval(as.matrix(STcoefs), prob=0.95)[1:3,1],
                    upci=posterior_interval(as.matrix(STcoefs), prob=0.95)[1:3,2])

STslopes <- data.frame(coef=colnames(STcoefs)[4:9],
                       mean=colMeans(STcoefs)[4:9],
                       lowci=posterior_interval(as.matrix(STcoefs), prob=0.95)[4:9,1],
                       upci=posterior_interval(as.matrix(STcoefs), prob=0.95)[4:9,2])

STareacoef <- data.frame(coef=c("t_area_s","t_area_y"),
                         mean=c(mean(STcoefs$t_logsigma_s+STcoefs$t_logmax_s),mean(STcoefs$t_logsigma_y+STcoefs$t_logmax_y)),
                         lowci=c(posterior_interval(as.matrix(STcoefs$t_logsigma_s+STcoefs$t_logmax_s), prob=0.95)[,1],posterior_interval(as.matrix(STcoefs$t_logsigma_y+STcoefs$t_logmax_y), prob=0.95)[,1]),
                         upci=c(posterior_interval(as.matrix(STcoefs$t_logsigma_s+STcoefs$t_logmax_s), prob=0.95)[,2],posterior_interval(as.matrix(STcoefs$t_logsigma_y+STcoefs$t_logmax_y), prob=0.95)[,2]))
STslopes <- rbind(STslopes,STareacoef) 

# Unscaled timing slope means and CIs
Smu_unscaled <- round(c(mean(STcoefs$t_mu_s*14.08795), posterior_interval(as.matrix(STcoefs$t_mu_s*14.08795), prob=0.95)),2)
Tmu_unscaled <- round(c(mean(STcoefs$t_mu_y*14.08795), posterior_interval(as.matrix(STcoefs$t_mu_y*14.08795), prob=0.95)),2)

# Exponentiated height slope means and CIs
Smax <- round(c(mean(exp(STcoefs$t_logmax_s)), posterior_interval(as.matrix(exp(STcoefs$t_logmax_s)), prob=0.95)),2)
Tmax <- round(c(mean(exp(STcoefs$t_logmax_y)), posterior_interval(as.matrix(exp(STcoefs$t_logmax_y)), prob=0.95)),2)

# Exponentiated width slope means and CIs
Ssig <- round(c(mean(exp(STcoefs$t_logsigma_s)), posterior_interval(as.matrix(exp(STcoefs$t_logsigma_s)), prob=0.95)),2)
Tsig <- round(c(mean(exp(STcoefs$t_logsigma_y)), posterior_interval(as.matrix(exp(STcoefs$t_logsigma_y)), prob=0.95)),2)

# Exponentiated area slope means and CIs
Sarea <- round(c(mean(exp(STcoefs$t_logsigma_s+STcoefs$t_logmax_s)), posterior_interval(as.matrix(exp(STcoefs$t_logsigma_s+STcoefs$t_logmax_s)), prob=0.95)),2)
Tarea <- round(c(mean(exp(STcoefs$t_logsigma_y+STcoefs$t_logmax_y)), posterior_interval(as.matrix(exp(STcoefs$t_logsigma_y+STcoefs$t_logmax_y)), prob=0.95)),2)

# Differences between slopes in space and time
t_mu_dif <- (STcoefs$t_mu_s-STcoefs$t_mu_y)*14.08795
t_sig_dif <- exp(STcoefs$t_logsigma_s)-exp(STcoefs$t_logsigma_y)
t_max_dif <- exp(STcoefs$t_logmax_s)-exp(STcoefs$t_logmax_y)
t_area_dif <- exp(STcoefs$t_logsigma_s+STcoefs$t_logmax_s)-exp(STcoefs$t_logsigma_y+STcoefs$t_logmax_y)

STslopedif <- data.frame(coef=c("t_mu_dif","t_sig_dif","t_max_dif","t_area_dif"),
                         mean=c(mean(t_mu_dif),mean(t_sig_dif),mean(t_max_dif),mean(t_area_dif)),
                         lowci=c(posterior_interval(as.matrix(t_mu_dif), prob=0.95)[,1],posterior_interval(as.matrix(t_sig_dif), prob=0.95)[,1],posterior_interval(as.matrix(t_max_dif), prob=0.95)[,1],posterior_interval(as.matrix(t_area_dif), prob=0.95)[,1]),
                         upci=c(posterior_interval(as.matrix(t_mu_dif), prob=0.95)[,2],posterior_interval(as.matrix(t_sig_dif), prob=0.95)[,2],posterior_interval(as.matrix(t_max_dif), prob=0.95)[,2],posterior_interval(as.matrix(t_area_dif), prob=0.95)[,2]))
