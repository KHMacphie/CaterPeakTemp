library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(metafor)
library(MCMCglmm)
library(dplyr)

###############################################################
####   Sliding window to identify period when temperature  ####
#### best predicts each caterpillar phenological parameter ####
###############################################################

# Model for site by year (SY) combination predictions of caterpillar phenological parameters
write("data{
      int<lower=0> N;
      vector[N] date;
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
      vector[N_siteday] siteday_effects;
      vector[N_treeID] treeID_effects;
      vector[N_rec] rec_effects;
      vector[N_resid] resid_effects;
      siteday_effects = sd_siteday * siteday_scaled; 
      treeID_effects = sd_treeID * treeID_scaled; 
      rec_effects = sd_rec * rec_scaled; 
      resid_effects = sd_resid * resid_scaled;
      
      y_1 = logmax - ((date-(mu+site_effs[site_id,1]+year_effs[year_id,1]+siteyear_effs[siteyear_id,1])) .* 
      (date-(mu+site_effs[site_id,1]+year_effs[year_id,1]+siteyear_effs[siteyear_id,1]))) ./ 
      (2.*(exp(logsigma+site_effs[site_id,2]+year_effs[year_id,2]+siteyear_effs[siteyear_id,2])).*
      (exp(logsigma+site_effs[site_id,2]+year_effs[year_id,2]+siteyear_effs[siteyear_id,2]))) 
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
      
      "SiteYearPeaks.stan"
)
stanc("SiteYearPeaks.stan")
SiteYearPeaksMod <- stan_model("SiteYearPeaks.stan")

# Load caterpillar data
data_cater <- read.csv("~/data_cater.csv")

# Data in stan model format
stan_data_Peaks <-list(
  N=nrow(data_cater),
  date=data_cater$datescaled,
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
SY_peaks <- sampling(object=SiteYearPeaksMod, data=stan_data_Peaks, 
                                 iter=4000, warmup=1500, chains=4, thin=5, cores=4,
                                 control = list(adapt_delta = 0.9),
                                 pars=c("mu", "logsigma","logmax","sd_site","sd_year","sd_siteyear",
                                        "sd_resid","site_effs","year_effs","siteyear_effs",
                                        "site_scaled","year_scaled","siteyear_scaled",
                                        "L", "omega")) 

# Extract posterior distributions and combine chains for each parameter
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
}

post <- stanpost(SY_peaks)

# Create df for SY parameter estimates 
sy.effs <- data.frame(site_id=as.numeric(as.factor(data_cater$site)), 
                      site=data_cater$site,
                      year_id=as.numeric(as.factor(data_cater$year)), 
                      year=data_cater$year,
                      siteyear_id=as.numeric(as.factor(data_cater$siteyear)), 
                      siteyear=data_cater$siteyear)
sy.effs <- distinct(sy.effs)
sy.effs.long <- sy.effs[rep(seq_len(nrow(sy.effs)), each = 3), ]

# Extract/calculate mode and vcv matrix for each SY parameter estimate
for(i in 1:length(sy.effs$site_id)){
  
  mu.i <- mcmc(post[,1] + post[,sy.effs$site_id[i]+13] + post[,sy.effs$year_id[i]+145] + post[,sy.effs$siteyear_id[i]+169])
  logsigma.i <- mcmc(post[,2] + post[,sy.effs$site_id[i]+57] + post[,sy.effs$year_id[i]+153] + post[,sy.effs$siteyear_id[i]+462])
  logmax.i <- mcmc(post[,3] + post[,sy.effs$site_id[i]+101] + post[,sy.effs$year_id[i]+161] + post[,sy.effs$siteyear_id[i]+755])
  
  sy.effs.long$para[i*3-2] <- "mu"
  sy.effs.long$para[i*3-1] <- "logsigma"
  sy.effs.long$para[i*3] <- "logmax"
  sy.effs.long$yi[i*3-2] <- posterior.mode(mu.i)
  sy.effs.long$yi[i*3-1] <- posterior.mode(logsigma.i)
  sy.effs.long$yi[i*3] <- posterior.mode(logmax.i)
  
  df <- data.frame(mu.i=mu.i, logsigma.i=logsigma.i, logmax.i=logmax.i)
  df.m <- as.matrix(df)
  df.cov <- cov(df.m)
  
  sy.effs.long[i*3-2,9:11] <- df.cov[1,]
  sy.effs.long[i*3-1,9:11] <- df.cov[2,]
  sy.effs.long[i*3,9:11] <- df.cov[3,]
  
  mu.i <- NULL
  logsigma.i <- NULL
  logmax.i <- NULL
  df <- NULL
  df.m <- NULL
  df.cov <- NULL
}

colnames(sy.effs.long)[9:11] <- c("vi.mu","vi.logsigma","vi.logmax")
sy.effs.long$id <- rep(seq(1,length(sy.effs$siteyear),1),1, each=3)

# VCV matrix for metafor
V <- bldiag(lapply(split(sy.effs.long[,c("vi.mu","vi.logsigma","vi.logmax")], sy.effs.long$id), as.matrix))

#Column for linking yi for each para to the temperature slopes
sy.effs.long$para.mu <- ifelse(sy.effs.long$para=="mu", 1, 0)
sy.effs.long$para.logsigma <- ifelse(sy.effs.long$para=="logsigma", 1, 0)
sy.effs.long$para.logmax <- ifelse(sy.effs.long$para=="logmax", 1, 0)

# Load average daily temperatures
temp <- read.csv("~/sy_daily_temp.csv")

# Organise windows of temperature to be tested for logsigma and logmax
firstcol <- 2 #first column with a day's temp
minwin <- 28 #minimum duration
shiftdur <- 14 #change duration by
shiftst <- 14 #change start day by
maxwin <- ncol(temp)-(firstcol-1)# no maximum duration

duration <- c(seq((minwin-1),(maxwin-1),shiftdur)) #duration is -1 because will be start+duration
start <- c(seq(firstcol,length(colnames(temp)),shiftst)) #column to start in
bounds <- data.frame(start=rep(start,length(duration)))
bounds$duration <- rep(duration, 1, each=length(start))
bounds$end <- bounds$start+bounds$duration
bounds <- bounds[-c(which(bounds$end>length(colnames(temp)))),] #remove ones which over run days in temp data file
bounds$ID <- paste(colnames(temp)[bounds$start],"_",(bounds$duration+1)) #ID for each window

mv.wndws.sigmax <- data.frame(logsigma.ID=bounds$ID,logmax.ID=bounds$ID)
mv.wndws.sigmax <- expand.grid(mv.wndws.sigmax) # all window options for logsig and logmax paired up

# Organise windows of temperature to be tested for mu
firstcol.mu <- 2 #first column with a day's temp
minwin.mu <- 28 #minimum duration
shiftdur.mu <- 14 #change duration by
shiftst.mu <- 7 #change start day by
maxwin.mu <- ncol(temp)-(firstcol-1)# no maximum duration
maxst.mu <- 44 #column of latest start day 100

duration.mu <- c(seq((minwin.mu-1),(maxwin.mu-1),shiftdur.mu)) #duration is -1 becasue  will be start+duration
start.mu <- c(seq(firstcol.mu,length(colnames(temp)),shiftst.mu)) #column to start in
bounds.mu <- data.frame(start=rep(start.mu,length(duration.mu)))
bounds.mu$duration <- rep(duration.mu, 1, each=length(start.mu))
bounds.mu$end <- bounds.mu$start+bounds.mu$duration
bounds.mu <- bounds.mu[-c(which(bounds.mu$end>length(colnames(temp)))),] #remove ones which over run temp data file
bounds.mu <- bounds.mu[-c(which(bounds.mu$start>maxst.mu)),] #remove ones which over run temp data file
bounds.mu$ID <- paste(colnames(temp)[bounds.mu$start],"_",(bounds.mu$duration+1))

# All window combinations for the three parameters
mv.wndws <- data.frame(mu.ID=character(0),logsigma.ID=character(0),logmax.ID=character(0))

for(i in 1:nrow(bounds.mu)){
  df <- data.frame(mu.ID=rep(bounds.mu$ID[i],nrow(mv.wndws.sigmax)),logsigma.ID=mv.wndws.sigmax$logsigma.ID,logmax.ID=mv.wndws.sigmax$logmax.ID)
  mv.wndws <- rbind(mv.wndws, df)
}

# Putting all window IDs into one file 
bounds <- rbind(bounds,bounds.mu)
bounds <- distinct(bounds)

# Columns for info from sliding window
mv.wndws$mu.st <- rep(NA, nrow(mv.wndws))
mv.wndws$mu.end <- rep(NA, nrow(mv.wndws))
mv.wndws$mu.dur <- rep(NA, nrow(mv.wndws))
mv.wndws$logsigma.st <- rep(NA, nrow(mv.wndws))
mv.wndws$logsigma.end <- rep(NA, nrow(mv.wndws))
mv.wndws$logsigma.dur <- rep(NA, nrow(mv.wndws))
mv.wndws$logmax.st<- rep(NA, nrow(mv.wndws))
mv.wndws$logmax.end <- rep(NA, nrow(mv.wndws))
mv.wndws$logmax.dur <- rep(NA, nrow(mv.wndws))
mv.wndws$AIC <- rep(NA, nrow(mv.wndws))
mv.wndws$logLik <- rep(NA, nrow(mv.wndws))
mv.wndws$k <- rep(NA, nrow(mv.wndws))
mv.wndws$temp.mu.beta <- rep(NA, nrow(mv.wndws))
mv.wndws$temp.mu.se <- rep(NA, nrow(mv.wndws))
mv.wndws$temp.mu.p <- rep(NA, nrow(mv.wndws))
mv.wndws$temp.logsigma.beta <- rep(NA, nrow(mv.wndws))
mv.wndws$temp.logsigma.se <- rep(NA, nrow(mv.wndws))
mv.wndws$temp.logsigma.p <- rep(NA, nrow(mv.wndws))
mv.wndws$temp.logmax.beta <- rep(NA, nrow(mv.wndws))
mv.wndws$temp.logmax.se <- rep(NA, nrow(mv.wndws))
mv.wndws$temp.logmax.p <- rep(NA, nrow(mv.wndws))

# Df for temp data to store in each iteration of the loop and after
slidwin.temp <- data.frame(siteyear=temp$siteyear)

# Run the models
for(i in 1:nrow(mv.wndws)){
  #mean time window
  slidwin.temp$temp.mu <- apply(temp[,bounds$start[which(bounds$ID==mv.wndws$mu.ID[i])]
                                      :bounds$end[which(bounds$ID==mv.wndws$mu.ID[i])]
                                      ],1,mean)
  slidwin.temp$temp.logsigma <- apply(temp[,bounds$start[which(bounds$ID==mv.wndws$logsigma.ID[i])]
                                            :bounds$end[which(bounds$ID==mv.wndws$logsigma.ID[i])]
                                            ],1,mean)
  slidwin.temp$temp.logmax <- apply(temp[,bounds$start[which(bounds$ID==mv.wndws$logmax.ID[i])]
                                          :bounds$end[which(bounds$ID==mv.wndws$logmax.ID[i])]
                                          ],1,mean)
  #pair with SY cater peak
  store <- pmatch(sy.effs.long$siteyear, slidwin.temp$siteyear, duplicates.ok = T)
  sy.effs.long$temp.mu <- slidwin.temp$temp.mu[store]
  sy.effs.long$temp.logsigma <- slidwin.temp$temp.logsigma[store]
  sy.effs.long$temp.logmax <- slidwin.temp$temp.logmax[store]
  
  multivar.mod<--999
  #run model
  try(multivar.mod <- rma.mv(yi, V, mods=~ para + temp.mu:para.mu + temp.logsigma:para.logsigma + temp.logmax:para.logmax - 1, random=list(~para|siteyear, ~para|year), data=sy.effs.long, struct="UN", method="ML")) 
  
  if(multivar.mod[[1]][1]!=-999){
    #store window details
    mv.wndws$mu.st[i] <- bounds$start[which(bounds$ID==mv.wndws$mu.ID[i])]+56
    mv.wndws$mu.end[i] <- bounds$end[which(bounds$ID==mv.wndws$mu.ID[i])]+56
    mv.wndws$mu.dur[i] <- length(bounds$start[which(bounds$ID==mv.wndws$mu.ID[i])]
                                 :bounds$end[which(bounds$ID==mv.wndws$mu.ID[i])])
    mv.wndws$logsigma.st[i] <- bounds$start[which(bounds$ID==mv.wndws$logsigma.ID[i])]+56
    mv.wndws$logsigma.end[i] <- bounds$end[which(bounds$ID==mv.wndws$logsigma.ID[i])]+56
    mv.wndws$logsigma.dur[i] <- length(bounds$start[which(bounds$ID==mv.wndws$logsigma.ID[i])]
                                       :bounds$end[which(bounds$ID==mv.wndws$logsigma.ID[i])])
    mv.wndws$logmax.st[i] <- bounds$start[which(bounds$ID==mv.wndws$logmax.ID[i])]+56
    mv.wndws$logmax.end[i] <- bounds$end[which(bounds$ID==mv.wndws$logmax.ID[i])]+56
    mv.wndws$logmax.dur[i] <- length(bounds$start[which(bounds$ID==mv.wndws$logmax.ID[i])]
                                     :bounds$end[which(bounds$ID==mv.wndws$logmax.ID[i])])
    
    #store model details
    mv.wndws$AIC[i] <- AIC(multivar.mod)
    mv.wndws$logLik[i] <- logLik(multivar.mod)[1]
    mv.wndws$k[i] <- multivar.mod$k
    mv.wndws$temp.mu.beta[i] <- multivar.mod$beta[4]
    mv.wndws$temp.mu.se[i] <- multivar.mod$se[4]
    mv.wndws$temp.mu.p[i] <- multivar.mod$pval[4]
    mv.wndws$temp.logsigma.beta[i] <- multivar.mod$beta[5]
    mv.wndws$temp.logsigma.se[i] <- multivar.mod$se[5]
    mv.wndws$temp.logsigma.p[i] <- multivar.mod$pval[5]
    mv.wndws$temp.logmax.beta[i] <- multivar.mod$beta[6]
    mv.wndws$temp.logmax.se[i] <- multivar.mod$se[6]
    mv.wndws$temp.logmax.p[i] <- multivar.mod$pval[6]
    
    slidwin.temp$temp.mu <- NULL
    slidwin.temp$temp.logsigma <- NULL
    slidwin.temp$temp.logmax <- NULL
    store <- NULL
    sy.effs.long$temp.mu <- NULL
    sy.effs.long$temp.logsigma <- NULL
    sy.effs.long$temp.logmax <- NULL
    multivar.mod <- NULL
    print(i)
    }
  }



# Identify row with best fitting windows
bestwindow <- which(mv.wndws$AIC==min(mv.wndws$AIC,na.rm=T))

# Extract the mean temperatures for each SY in those windows and pair with caterpillar data
slidwin.temp$temp.mu <- apply(temp[,bounds$start[which(bounds$ID==mv.wndws$mu.ID[bestwindow])]
                                    :bounds$end[which(bounds$ID==mv.wndws$mu.ID[bestwindow])]
                                    ],1,mean)
slidwin.temp$temp.logsigma <- apply(temp[,bounds$start[which(bounds$ID==mv.wndws$logsigma.ID[bestwindow])]
                                          :bounds$end[which(bounds$ID==mv.wndws$logsigma.ID[bestwindow])]
                                          ],1,mean)
slidwin.temp$temp.logmax <- apply(temp[,bounds$start[which(bounds$ID==mv.wndws$logmax.ID[bestwindow])]
                                        :bounds$end[which(bounds$ID==mv.wndws$logmax.ID[bestwindow])]
                                        ],1,mean)
data_cater <- merge(data_cater, slidwin.temp, by="siteyear")

#mean centre each temperature variable (already in online data)
data_cater$temp_m.cent <- data_cater$temp.mu-mean(data_cater$temp.mu)
data_cater$temp_ls.cent <- data_cater$temp.logsigma-mean(data_cater$temp.logsigma)
data_cater$temp_lm.cent <- data_cater$temp.logmax-mean(data_cater$temp.logmax)
