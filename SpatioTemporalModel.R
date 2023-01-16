library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(rstantools)
library(dplyr)
library(ggplot2)
library(forcats)
library(gridExtra)
library(grid)
library(mvtnorm)

##############################################################
#### Thermal sensitivity in the phenological distribution ####
##############################################################


#### Analysis ####

# Model code
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


# Stan data list
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


# Sampling
TempModel <- sampling(object=TempCaterMod, data=stan_data_Temp, 
                            iter=4500, warmup=2000, chains=4, thin=5, cores=4,
                            control = list(adapt_delta = 0.94)) 


#### Results ####

# Extract posterior distributions from the model
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

post <- stanpost(TempMod)

# Find parameters wanted
coefs <- post[,1:6]
SYcoefs <- post[,46248:47282]
coefvarcov <- post[,c(1:6,46235:46247,47283:47291)]

# Temp data as df by SYs
temps <- read.csv("~/sy_temp_by_para.csv")               

# Make SY effects into one column per SY:para (not site, year and siteyear) ####
SY <- data.frame(siteyear_id=as.numeric(as.factor(data_cater$siteyear)),
                 siteyear=as.factor(data_cater$siteyear),
                 site_id=as.numeric(as.factor(data_cater$site)),
                 site=as.factor(data_cater$site),
                 year_id=as.numeric(as.factor(data_cater$year)),
                 year=as.factor(data_cater$year))
SY <- distinct(SY)
store <- pmatch(SY$siteyear,temps$siteyear)
SY$temp.mu <- temps$temp.mu[store]  
SY$temp.ls <- temps$temp.logsigma[store]  
SY$temp.lm <- temps$temp.logmax[store] 
SY$t.mu.cent <- SY$temp.mu-mean(data_cater$temp.mu)
SY$t.ls.cent <- SY$temp.ls-mean(data_cater$temp.logsigma)
SY$t.lm.cent <- SY$temp.lm-mean(data_cater$temp.logmax)

mu.site <- c()
ls.site <- c()
lm.site <- c()
mu.year <- c()
ls.year <- c()
lm.year <- c()
mu.siteyear <- c()
ls.siteyear <- c()
lm.siteyear <- c()
SY.mu <- c()
SY.ls <- c()
SY.lm <- c()

for(i in 1:44){
  mu.site[[i]] <- SYcoefs[,i]
  ls.site[[i]] <- SYcoefs[,i+44]
  lm.site[[i]] <- SYcoefs[,i+88]
}

for(i in 1:8){
  mu.year[[i]] <- SYcoefs[,i+132]
  ls.year[[i]] <- SYcoefs[,i+140]
  lm.year[[i]] <- SYcoefs[,i+148]
}

for(i in 1:293){
  mu.siteyear[[i]] <- SYcoefs[,i+156]
  ls.siteyear[[i]] <- SYcoefs[,i+449]
  lm.siteyear[[i]] <- SYcoefs[,i+742]
}

for(i in 1:293){
  SY.mu[[i]] <- mu.siteyear[[i]] + mu.site[[SY$site_id[which(SY$siteyear_id==i)]]] + mu.year[[SY$year_id[which(SY$siteyear_id==i)]]] 
  SY.ls[[i]] <- ls.siteyear[[i]] + ls.site[[SY$site_id[which(SY$siteyear_id==i)]]] + ls.year[[SY$year_id[which(SY$siteyear_id==i)]]]
  SY.lm[[i]] <- lm.siteyear[[i]] + lm.site[[SY$site_id[which(SY$siteyear_id==i)]]] + lm.year[[SY$year_id[which(SY$siteyear_id==i)]]]
}

# Timing by temperature- slope and siteyear points
slope.mu <- data.frame(tempcent=seq(min(SY$t.mu.cent),max(SY$t.mu.cent),0.01))
slope.mu$temp <- slope.mu$tempcent + mean(data_cater$temp.mu)

for(i in 1:nrow(slope.mu)){ 
  slope.mu$mean.mu[i] <- mean((coefs$mu+coefs$t_mu*slope.mu$tempcent[i])*14.08795 + 147.9033)
  slope.mu$lci.mu[i] <- posterior_interval(as.matrix((coefs$mu+coefs$t_mu*slope.mu$tempcent[i])*14.08795 + 147.9033),prob=0.95)[1]
  slope.mu$uci.mu[i] <- posterior_interval(as.matrix((coefs$mu+coefs$t_mu*slope.mu$tempcent[i])*14.08795 + 147.9033),prob=0.95)[2]
}

for (i in 1:nrow(SY)){
  X <- (coefs$mu+coefs$t_mu*SY$t.mu.cent[i]+SY.mu[[i]])*14.08795 + 147.9033
  SY$mu.points[i] <- mean(X)
}

temppoints.mu <- data.frame(temp=c(mean(data_cater$temp.mu)+2,mean(data_cater$temp.mu)+1,mean(data_cater$temp.mu),mean(data_cater$temp.mu)-1,mean(data_cater$temp.mu)-2),
                            y=rep(132,5))
RedBlue <- c("#2E70FF", "#599EF2", "#FFE6A3", "#F58762", "#FA4444")
(muslope <- ggplot(data=slope.mu, aes(temp, mean.mu))+
    geom_line(lwd=1, col=1)+
    geom_ribbon(aes(x=temp, ymin=lci.mu, ymax=uci.mu), alpha=0.25, fill="gray")+
    geom_point(data=SY, aes(temp.mu, mu.points), size=0.7)+
    geom_point(data=temppoints.mu, aes(temp, y), colour=rev(RedBlue), pch=15, size=2.5)+
    xlab("Temperature (°C)")+
    ylab("Mean Timing (days)")+
    theme_bw()+
    coord_cartesian(ylim=c(133.6,179), xlim=c(2.8,8.8))+
    annotate(x=8.9, y=178.5, label="b", geom="text", col=1, cex=5)+
    theme(text=element_text(size= 15))) 

# Width by temperature- slope and siteyear points
#half variance from terms marginalised over added to mean because estimated as a lognormal distribution
var.ls <- coefvarcov[,8]^2 + coefvarcov[,11]^2 + coefvarcov[,14]^2
slope.ls <- data.frame(tempcent=seq(min(SY$t.ls.cent),max(SY$t.ls.cent),0.01))
slope.ls$temp <- slope.ls$tempcent + mean(data_cater$temp.logsigma)


for(i in 1:nrow(slope.ls)){    
  slope.ls$mean.sigma[i] <- mean((exp(coefs$logsigma+coefs$t_logsigma*slope.ls$tempcent[i] + var.ls/2))*14.08795)
  slope.ls$lci.sigma[i] <- posterior_interval(as.matrix((exp(coefs$logsigma+coefs$t_logsigma*slope.ls$tempcent[i] + var.ls/2))*14.08795),prob=0.95)[1]
  slope.ls$uci.sigma[i] <- posterior_interval(as.matrix((exp(coefs$logsigma+coefs$t_logsigma*slope.ls$tempcent[i] + var.ls/2))*14.08795),prob=0.95)[2]
}

for (i in 1:nrow(SY)){
  X <- exp(coefs$logsigma+coefs$t_logsigma*SY$t.ls.cent[i]+SY.ls[[i]])*14.08795
  SY$sig.points[i] <- mean(X)
}

temppoints.sig <- data.frame(temp=c(mean(data_cater$temp.logsigma)+2,mean(data_cater$temp.logsigma)+1,mean(data_cater$temp.logsigma),mean(data_cater$temp.logsigma)-1,mean(data_cater$temp.logsigma)-2),
                             y=rep(4.5,5))

(sigslope <- ggplot(data=slope.ls, aes(temp, mean.sigma))+
    geom_line(lwd=1, col=1)+
    geom_ribbon(aes(x=temp, ymin=lci.sigma, ymax=uci.sigma), alpha=0.25, fill="gray")+
    geom_point(data=SY, aes(temp.ls, sig.points), size=0.7)+
    geom_point(data=temppoints.sig, aes(temp, y), colour=rev(RedBlue), pch=15, size=2.5)+
    xlab("Temperature (°C)")+
    ylab("Width (days)")+
    theme_bw()+
    annotate(x=10, y=25.5, label="d", geom="text", col=1, cex=5)+
    coord_cartesian(ylim=c(5.25,26), xlim=c(5.4,9.95))+
    theme(text=element_text(size= 15))) 

# Height by temperature- slope and siteyear points 
#half variance from terms marginalised over added to mean because estimated as a lognormal distribution
var.lm <- coefvarcov[,9]^2 + coefvarcov[,12]^2 + coefvarcov[,15]^2
var.basic <- coefvarcov[,16]^2 + coefvarcov[,17]^2 + coefvarcov[,18]^2 + coefvarcov[,19]^2
slope.lm <- data.frame(tempcent=seq(min(SY$t.lm.cent),max(SY$t.lm.cent),0.01))
slope.lm$temp <- slope.lm$tempcent + mean(data_cater$temp.logmax)

for(i in 1:nrow(slope.lm)){    
  slope.lm$mean.max[i] <- mean(exp(coefs$logmax+coefs$t_logmax*slope.lm$tempcent[i] + (var.lm+var.basic)/2))
  slope.lm$lci.max[i] <- posterior_interval(as.matrix(exp(coefs$logmax+coefs$t_logmax*slope.lm$tempcent[i] + (var.lm+var.basic)/2)),prob=0.95)[1]
  slope.lm$uci.max[i] <- posterior_interval(as.matrix(exp(coefs$logmax+coefs$t_logmax*slope.lm$tempcent[i] + (var.lm+var.basic)/2)),prob=0.95)[2]
  
  slope.lm$mean.lmax[i] <- mean(coefs$logmax+coefs$t_logmax*slope.lm$tempcent[i])
  slope.lm$lci.lmax[i] <- posterior_interval(as.matrix(coefs$logmax+coefs$t_logmax*slope.lm$tempcent[i]),prob=0.95)[1]
  slope.lm$uci.lmax[i] <- posterior_interval(as.matrix(coefs$logmax+coefs$t_logmax*slope.lm$tempcent[i]),prob=0.95)[2]
}

for (i in 1:nrow(SY)){
  X <- exp(coefs$logmax+coefs$t_logmax*SY$t.lm.cent[i]+SY.lm[[i]]+var.basic/2)
  SY$max.points[i] <- mean(X)
  Y <- (coefs$logmax+coefs$t_logmax*SY$t.lm.cent[i]+SY.lm[[i]])
  SY$lmax.points[i] <- mean(Y)
}

cropped <- subset(SY, max.points>=2)
(lmaxslope <- ggplot(data=slope.lm, aes(temp, mean.lmax))+
    geom_line(lwd=1, col=1)+
    geom_ribbon(aes(x=temp, ymin=lci.lmax, ymax=uci.lmax), alpha=0.25, fill="gray")+
    geom_point(data=SY, aes(temp.lm, lmax.points), size=0.4)+
    geom_point(data=cropped, aes(temp.lm, lmax.points), size=0.6, colour="red")+
    xlab("")+
    ylab("logHeight")+
    theme_bw()+
    theme(text=element_text(size= 10),axis.title.x = element_blank())) 

temppoints.max <- data.frame(temp=c(mean(data_cater$temp.logmax)+2,mean(data_cater$temp.logmax)+1,mean(data_cater$temp.logmax),mean(data_cater$temp.logmax)-1,mean(data_cater$temp.logmax)-2),
                             y=rep(-0.05,5))

(maxslope <- ggplot(data=slope.lm, aes(temp, mean.max))+
    geom_line(lwd=1, col=1)+
    geom_ribbon(aes(x=temp, ymin=lci.max, ymax=uci.max), alpha=0.25, fill="gray")+
    geom_point(data=SY, aes(temp.lm, max.points), size=0.7)+
    geom_point(data=temppoints.max, aes(temp, y), colour=rev(RedBlue), pch=15, size=2.5)+
    xlab("Temperature (°C)")+
    ylab("Height (abundance/branch)")+
    theme_bw()+
    coord_cartesian(ylim = c(0.02,2), xlim=c(5.1,11.3))+
    theme(text=element_text(size= 15))+ 
    annotate(x=11.4, y=2, label="c", geom="text", col=1, cex=5)+
    annotation_custom(ggplotGrob(lmaxslope), 
                      xmin = 4.8, xmax = 7.4,
                      ymin = 0.9, ymax = 2.05)) 

# Temp windows plot
windows <- data.frame(Parameter=c("Mean timing","Mean timing", "Height", "Height", "Width", "Width"), 
                      Date=c(65, 106, 100, 141, 58, 155)) 
windows$Parameter <- relevel(as.factor(windows$Parameter), ref="Mean timing")

windows2 <- data.frame(Parameter=c("Mean timing","Mean timing", "Height", "Height", "Width", "Width"), 
                       Date=c(65, 106, 100, 141, 58, 85)) 
windows2$Parameter <- relevel(as.factor(windows2$Parameter), ref="Mean timing")

windows3 <- data.frame(Parameter=c("Mean timing","Mean timing", "Height", "Height", "Width", "Width"), 
                       Date=c(65, 106, 100, 141, 58, 141)) 
windows3$Parameter <- relevel(as.factor(windows3$Parameter), ref="Mean timing")

windows4 <- data.frame(Parameter=c("Mean timing","Mean timing", "Height", "Height", "Width", "Width"), 
                       Date=c(65, 106, 100, 141, 72, 155)) 
windows4$Parameter <- relevel(as.factor(windows4$Parameter), ref="Mean timing")

windows5 <- data.frame(Parameter=c("Mean timing","Mean timing", "Height", "Height", "Width", "Width"), 
                       Date=c(65, 106, 100, 141, 58, 127)) 
windows5$Parameter <- relevel(as.factor(windows5$Parameter), ref="Mean timing")

windows6 <- data.frame(Parameter=c("Mean timing","Mean timing", "Height", "Height", "Width", "Width"), 
                       Date=c(65, 106, 100, 141, 58, 99)) 
windows6$Parameter <- relevel(as.factor(windows6$Parameter), ref="Mean timing")

mycolblack <- rgb(100, 100, 100, max = 250, alpha = 50, names = "blacktrans")

(windowplot <- ggplot(windows, aes(fct_rev(Parameter), Date))+
    geom_line(lwd=1.5, col=1)+
    coord_flip(ylim = c(58,161))+
    theme_bw()+
    xlab("")+
    ylab("Ordinal Date")+
    annotate(x=3.3, y=161, label="a", geom="text", col=1, cex=5)+
    theme(text=element_text(size= 15), axis.title.y = element_blank()))

Fig3 <- grid.arrange(windowplot,muslope,maxslope,sigslope, nrow = 4, heights = c(0.4,1,1,1)) #11.5x4.5



## Expectation over time for different temps- requires simulation under model ##
#random draw from random terms inc. covariance for s, y and sy
#into equations for date by temp estimates
#simulate multipel times from same variance estimates- take mean 
#repeat for every iteration

# f1 = logmax + temp*t + random(s+y+sy)
# f2 = date - mu - temp*t - random(s+y+sy)
# f3 = sqrt(2) * exp(logsig + temp*t + random(s+y+sy))
# f4 = sum other ran effs

# Group parameters for ease
ranef.cov <- coefvarcov[,7:15]
ranef.basic <- coefvarcov[,16:19]
omega <- coefvarcov[,20:28]

# Values to calculate from within sim
temp <- c(-2, -1, 0, 1, 2)
dat.cnt <- seq(-3.5,4.5,0.05)
areatemp <- seq(-2,2,0.1)
N.h <- 10000 #no. of sims

# Set up dfs for simulation
Esim <- data.frame(matrix(NA, nrow = length(dat.cnt)*length(temp), ncol = N.h+2)) #for mean expectations
Emean <- data.frame(matrix(NA, nrow = length(dat.cnt)*length(temp), ncol = nrow(coefvarcov)+2)) #for mean expectations
Asim <- data.frame(matrix(NA, nrow = length(areatemp), ncol = N.h+1))
Amean <- data.frame(matrix(NA, nrow = length(areatemp), ncol = nrow(coefvarcov)+1))
colnames(Esim)[c(1,2)] <- c("dat.cnt", "temp")
colnames(Emean)[c(1,2)] <- c("dat.cnt", "temp")
colnames(Asim)[1] <- c("temp")
colnames(Amean)[1] <- c("temp")

Esim$dat.cnt <- rep(dat.cnt,5)
Esim$temp <- rep(temp, each=length(dat.cnt),1)
Emean$dat.cnt <- rep(dat.cnt,5)
Emean$temp <- rep(temp, each=length(dat.cnt),1)

Asim$temp <- areatemp
Amean$temp <- areatemp

# Simulate predictions from the model
system.time(for (i in 1:nrow(coefvarcov)){ # for each iteration
  #one correlation matrix between the 3 parameters, difference variance for each random term
  ##covar matrix: 1=mu, 2=logsigma, 3=logmax
  cor.matrix <- matrix(c(omega[i,1], omega[i,2], omega[i,3],
                         omega[i,4], omega[i,5], omega[i,6],
                         omega[i,7], omega[i,8], omega[i,9]), nrow=3, ncol=3) #correlation matrix- called omega in model code
  
  site.sd <- diag(c(ranef.cov[i,1], ranef.cov[i,2], ranef.cov[i,3])) # for mu, logsigma and logmax respectively
  year.sd <- diag(c(ranef.cov[i,4], ranef.cov[i,5], ranef.cov[i,6]))
  styr.sd <- diag(c(ranef.cov[i,7], ranef.cov[i,8], ranef.cov[i,9]))
  
  site.covar <- site.sd%*%cor.matrix%*%site.sd  
  year.covar <- year.sd%*%cor.matrix%*%year.sd
  styr.covar <- styr.sd%*%cor.matrix%*%styr.sd
  
  for(h in 1:N.h){ #for one iteration draw 50000 times
    site.ef <- rmvnorm(n=1, mean=c(0,0,0), sigma=site.covar) #rmvnorm random number generator for multivariate with VCV
    year.ef <- rmvnorm(n=1, mean=c(0,0,0), sigma=year.covar)
    styr.ef <- rmvnorm(n=1, mean=c(0,0,0), sigma=styr.covar)
    
    stdy.ef <- rnorm(1, mean=0, sd=ranef.basic[i,1])
    tree.ef <- rnorm(1, mean=0, sd=ranef.basic[i,2])
    rcrd.ef <- rnorm(1, mean=0, sd=ranef.basic[i,3])
    resd.ef <- rnorm(1, mean=0, sd=ranef.basic[i,4])
    
    f1 <- coefs[i,3] + coefs[i,6]*Esim$temp + site.ef[3] + year.ef[3] + styr.ef[3]
    f2 <- Esim$dat.cnt - coefs[i,1] - coefs[i,4]*Esim$temp - site.ef[1] - year.ef[1] - styr.ef[1]
    f3 <- sqrt(2) * exp(coefs[i,2] + coefs[i,5]*Esim$temp + site.ef[2] + year.ef[2] + styr.ef[2])
    
    f4 <- stdy.ef + tree.ef + rcrd.ef + resd.ef
    
    Esim[,h+2] <- exp(f1+f4-(f2/f3)^2)
    
    f1 <- NULL
    f2 <- NULL
    f3 <- NULL
    f4 <- NULL
    
    lm.int <- coefs[i,3] + site.ef[3] + year.ef[3] + styr.ef[3] + stdy.ef + tree.ef + rcrd.ef + resd.ef
    ls.int <- coefs[i,2] + site.ef[2] + year.ef[2] + styr.ef[2]
    Asim[,h+1] <- exp(ls.int + lm.int + (coefs[i,5]+coefs[i,6])*Asim$temp + log(sqrt(2*pi)))
    
  } 
  
  Emean[,i+2] <- apply(Esim[,3:ncol(Esim)], 1, mean)
  Amean[,i+1] <- apply(Asim[,2:ncol(Asim)], 1, mean)
  
  Esim[,c(3:ncol(Esim))] <- NULL  
  Asim[,c(2:ncol(Asim))] <- NULL 
  print(paste("iteration",i))
  
} )

# Mean and CIs for estimate at each date or temperature
Emean$mean <- apply(Emean[,3:ncol(Emean)], 1, mean)
Emean$lwci <- posterior_interval(t(as.matrix(Emean[1:nrow(Emean),3:(ncol(Emean)-1)])), prob=0.95)[,1]
Emean$upci <- posterior_interval(t(as.matrix(Emean[1:nrow(Emean),3:(ncol(Emean)-2)])), prob=0.95)[,2]
Emean$dat <- Emean$dat.cnt*14.08795 + 147.9033

Amean$mean.unscal <- apply((Amean[,2:ncol(Amean)]*14.08795), 1, mean)
Amean$lwci.unscal <- posterior_interval(t(as.matrix(Amean[1:nrow(Amean),2:(ncol(Amean)-1)]*14.08795)), prob=0.95)[,1]
Amean$upci.unscal <- posterior_interval(t(as.matrix(Amean[1:nrow(Amean),2:(ncol(Amean)-2)]*14.08795)), prob=0.95)[,2]

# Plot peaks for abundance estimate by date at each temp
RedBlue <- c("#2E70FF", "#599EF2", "#FFE6A3", "#F58762", "#FA4444")
durline0.05 <- data.frame(x=c(100,200),y=c(0.05,0.05))
durline0.1 <- data.frame(x=c(100,200),y=c(0.1,0.1))
(peakplot <- ggplot(Emean, aes(dat, mean, col=as.factor(temp)))+
    geom_line(lwd=0.9)+
    geom_line(data=durline0.05, aes(x,y), linetype="dashed", col=1)+
    geom_line(data=durline0.1, aes(x,y), linetype="dashed", col=1)+
    annotate(x=120, y=0.38, label="a", geom="text", col=1, cex=5)+
    theme_bw()+
    xlab("Ordinal Date")+
    ylab("Abundance")+
    theme(text = element_text(size=15),
          plot.margin = unit(c(5.5,0,5.5,5.5), "pt"))+
    coord_cartesian(ylim = c(0.01,0.375), xlim=c(119,174))+
    guides(color = "none")+
    scale_colour_manual(values=RedBlue)) #6x4.5" 


# Duration CIs plot
#find dates for each iteration that are closest to 0.05 and 0.1 
#in Emean to estimate the duration

Dur.1 <- data.frame()
Dur.05 <- data.frame()
dur <- data.frame(matrix(NA,ncol=5,nrow=10))
colnames(dur) <- c("temp", "abund", "mean", "lowci", "upci")

for(j in 1:5){
  
  df <- subset(Emean, temp==unique(Emean$temp)[j])
  
  for(i in 3:2002){
    low.1 <- df$dat[which(abs(df[,i][1:which(df[,i]==max(df[,i]))]-0.1)==min(abs(df[,i][1:which(df[,i]==max(df[,i]))]-0.1)))]
    up.1 <- df$dat[which(df[,i]==max(df[,i]))-1+which(abs(df[,i][which(df[,i]==max(df[,i])):nrow(df)]-0.1)==min(abs(df[,i][which(df[,i]==max(df[,i])):nrow(df)]-0.1)))]
    Dur.1[i-2,j] <- up.1-low.1
    
    low.05 <- df$dat[which(abs(df[,i][1:which(df[,i]==max(df[,i]))]-0.05)==min(abs(df[,i][1:which(df[,i]==max(df[,i]))]-0.05)))]
    up.05 <- df$dat[which(df[,i]==max(df[,i]))-1+which(abs(df[,i][which(df[,i]==max(df[,i])):nrow(df)]-0.05)==min(abs(df[,i][which(df[,i]==max(df[,i])):nrow(df)]-0.05)))]
    Dur.05[i-2,j] <- up.05-low.05
  }
  colnames(Dur.1)[j] <- paste("temp",unique(Emean$temp)[j])
  colnames(Dur.05)[j] <- paste("temp",unique(Emean$temp)[j])
  
  dur$temp[(j*2-1):(j*2)] <- c(unique(Emean$temp)[j],unique(Emean$temp)[j])
  dur$abund[(j*2-1):(j*2)] <- c("0.1","0.05")
  dur$mean[(j*2-1):(j*2)] <- c(mean(Dur.1[,j]),mean(Dur.05[,j]))
  dur[(j*2-1),4:5] <- c(posterior_interval(as.matrix(Dur.1[,j]), prob=0.95))
  dur[(j*2),4:5] <- c(posterior_interval(as.matrix(Dur.05[,j]), prob=0.95))
  
}

# Posterior mean and CIs for difference in duration at 2 and -2 oC
Dif_0.05 <- Dur.05$`temp 2` - Dur.05$`temp -2`
Dif_0.1 <- Dur.1$`temp 2` - Dur.1$`temp -2`

mean(Dif_0.05)
posterior_interval(as.matrix(Dif_0.05), prob=0.95)
mean(Dif_0.1)
posterior_interval(as.matrix(Dif_0.1), prob=0.95)

# Duration plot
label <- data.frame(abund="0.1", lab="b")

(durplot <- ggplot(dur, aes(temp, mean, colour=as.factor(temp)))+
    geom_point(size=2)+
    geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
    geom_hline(yintercept=0, linetype="dotted", color = "black")+
    geom_text(x=2.45, y=2, aes(label=lab), data=label, col=1, cex=5)+
    theme_bw()+
    xlab("")+
    ylab("Duration (days)")+
    coord_flip()+
    xlim(-2.5,2.5)+
    guides(color = "none")+
    scale_colour_manual(values=RedBlue)+ 
    facet_grid(rows = vars(fct_rev(as.factor(abund))),switch="y")+#+
    theme(text = element_text(size=15),
          strip.background = element_blank(),
          axis.text.y = element_blank(), axis.ticks = element_blank(),
          plot.margin = unit(c(5.5,9,5.5,0), "pt")))

# Area by temperature- slope and siteyear points 
temppoints.area <- data.frame(temp=c(2,1,0,-1,-2),
                              y=rep(2.5,5))

(areaslope <- ggplot(Amean, aes(temp, mean.unscal))+
    geom_line(lwd=1, col=1)+
    geom_ribbon(aes(x=temp, ymin=lwci.unscal, ymax=upci.unscal), alpha=0.25, fill="gray")+
    geom_point(data=temppoints.area, aes(temp, y), colour=rev(RedBlue), pch=15, size=2.5)+
    annotate(x=-1.95, y=20.8, label="c", geom="text", col=1, cex=5)+
    xlab("Temperature (°C)")+
    ylab("Area (caterpillar days)")+
    theme_bw()+
    coord_cartesian(ylim=c(3.15,20.5))+
    theme(text=element_text(size= 15),
          plot.margin = unit(c(5.5,5.5,5.5,15), "pt"))) #4x3

Fig4 <- grid.arrange(peakplot,durplot,areaslope, ncol = 3, widths = c(1,0.9,1)) #4.2x10

grid.lines(x=c(0.265,0.38), y=c(0.359,0.60), gp=gpar(fill="black"), 
           arrow=arrow(type="closed", length=unit(2,"mm")))

grid.lines(x=c(0.265,0.37), y=c(0.251,0.33), gp=gpar(fill="black"), 
           arrow=arrow(type="closed", length=unit(2,"mm")))


# Mean and CIs for difference between +2 height and -2 height
E_n2 <- subset(Emean, temp==-2)
E_2 <- subset(Emean, temp==2)

max_n2 <- c()
max_2 <- c()

for(i in 1:2000){
  max_n2[i] <- max(E_n2[,i+2])
  max_2[i] <- max(E_2[,i+2])
}
hist(max_2/max_n2,200)
abline(v=mean(max_2/max_n2))

mean(max_2/max_n2)
posterior_interval(as.matrix(max_2/max_n2), prob=0.95)
