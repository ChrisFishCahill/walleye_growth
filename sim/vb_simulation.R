library(INLA)
library(TMB)
library(RandomFields)
library(raster)
library(RANN)
library(ggplot2)
library(dplyr)

rf_sim <- function(model, x, y) {
  set.seed(sample.int(1e5L, 1L))
  suppressMessages(
    RandomFields::RFsimulate(model = model, x = x, y = y)$variable1
  )
}

Sim_Fn <- function( n_years=nyears, n_stations=nstations, SpatialScale=SpatialScale, SD_O=SD_O,
 	                  rho=rho, beta1=beta1, beta0=beta0, sigma=sigma, Likelihood="LogNormal"){

  # Spatial-temporal model component:
  Loc = cbind( "x"=runif(n_stations, min=0,max=10), "y"=runif(n_stations, min=0,max=10) )
  rf_eps <- RandomFields::RMgauss(var=SD_O^2, scale=SpatialScale)
  #rf_eps <- RandomFields::RMmatern(nu=1, var=SD_O^2, scale = SpatialScale ) #SpatialScale = 1 / kappa
  eps_st <- list()
  for(t in 1:n_years){
    if(t==1){
     eps_st[[t]] <- rf_sim(rf_eps, Loc[,"x"], Loc[,"y"])
    } else {
     eps_st[[t]] <- rho*eps_st[[t-1]] + rf_sim(rf_eps, Loc[,"x"], Loc[,"y"])
    }
  }

  #Simulate a covariate x:
  x <- matrix(runif(n=n_years*n_stations, min=-1, max=1), nrow=n_stations, ncol=n_years)

  # Calculate Eta_it --> this is spatial temporal field for omega
  Eta_it   = array(NA, dim=c(n_stations,n_years))
  for(t in 1:n_years){
    Eta_it[,t] =  beta0  + beta1*x[ ,t] + eps_st[[t]] #fixed effects + Temporally evolving spatial field
  }

  #Correct for exceedingly low omegas
  Eta_it[which(Eta_it < 8.0)] <- NA

  #Simulate the assymptotic lengths:
  Epsilon_Linf <- rnorm(n_stations, mean=0, sd=SD_linf)
  Linfs <- Global_Linf + Epsilon_Linf

  #Simulate the hypothetical length at age 0's:
  Epsilon_t0  <- rnorm(n_stations, mean=0, sd=SD_t0)
  t0s <- Global_t0 + Epsilon_t0

  # Simulate data
  DF = NULL
  for(s in 1:n_stations){
    for(t in 1:n_years){
      Ages  <- sample(Age_Range, Nfish, replace=T)
      Linf  <-  Linfs[s] #Each lake gets a random deviate for Linf.
      T0    <- t0s[s]  #Each lake gets a random deviate for t0
      Omega <- Eta_it[s,t] #Predictors + Space-time for Omega

      lpreds <-   Linf*(1-exp(-(Omega/Linf) * (Ages - T0 ))) #Deterministic VB preds

      sigmas = cv*lpreds

      if(!any(is.na(lpreds))){
       if(Likelihood=="Normal"){Simulated_Length = rnorm(Nfish, mean=lpreds, sd=sigmas)}
       if(Likelihood=="Lognormal"){Simulated_Length = rlnorm0(Nfish, lpreds, cv)}
       if(Likelihood=="Gamma"){Simulated_Length = rgamma(Nfish, shape=1/cv^2, scale=lpreds*cv^2)}
      } else {
       Simulated_Length <- lpreds
      }
      if(any(is.na(lpreds))){Simulated_Length = lpreds}

      #if(min(Simulated_Length) < 0){Simulated_Length=0}

      Tmp = data.frame("Lake"=rep(s, Nfish) , "Year"=rep(t, Nfish), "x1"=rep(x[s,t], Nfish), "Age" = Ages,
                       "Simulated_Length"= Simulated_Length)
      DF = rbind(DF, Tmp)
  }}
  #plot(DF$Age,DF$Simulated_Length)

  DF = cbind( DF, 'Longitude'=Loc[DF[,'Lake'],1], 'Latitude'=Loc[DF[,'Lake'],2] )
  DF = data.frame(DF, row.names=NULL)

  # Return stuff
  Sim_List = list("DF"=DF, "Loc"=Loc, "Eta_it"=Eta_it,
                  "eps_st"=eps_st, "n_years"=n_years,
                  "n_stations"=n_stations)
  Sim_List[["Parameters"]] = c('SpatialScale'=SpatialScale, 'SigmaO'=SD_O,
                               'rho'=rho, "beta0"=beta0, "beta1"=beta1,
                               "sigma"=sigmas, 'Linf' = Global_Linf,
                               "t0" = Global_t0, 'SD_t0' = SD_t0,
                               'SD_linf'=SD_linf)
  return(Sim_List)
}

rlnorm0 <- function(n = 1, mean, coeffOfVar)
{
   sigma <- sqrt(log(coeffOfVar^2 + 1))
   mu <- log(mean) - sigma^2 / 2
   rlnorm(n, mu, sigma)
}

#---------------------------------------------------------------------------------------------------------------------------------------
#Estimate the model:
#---------------------------------------------------------------------------------------------------------------------------------------

#Simulation parameters:
#Note--the x and y ranges are 0-10

Age_Range    <-  0:25
Global_Linf  <- 55.70
SD_linf      <- 7.41
Global_t0    <- -1
SD_t0        <- 0.3
beta1        <- runif(1, -0.25, 0.25)
beta0        <- 14.79 #runif(1, 14, 16)
kappa        <- 0.2
SpatialScale <- 1/kappa
Range        <- sqrt(8) / kappa
SD_O         <- 4.66
rho          <- 0.91
cv           <-  0.069

Nfish        <- 50
nstations    <- 50
nyears       <- 15

Likelihood="Gamma" #Normal, Lognormal, or Gamma

setwd("sim/")
Version = "vb_sim_estimation"

# Compile
compile( paste0(Version,".cpp") )
dyn.load( dynlib(Version) )

Nsim = 1000
Estimates <- matrix(NA, nrow=Nsim, ncol=10)

set_par_value <- function(opt, par) {
  as.numeric(opt$par[par == names(opt$par)])
}

seed <-  sample.int(1e6, 1)
set.seed( seed )
#Run simulation:
for(i in 1:Nsim){

 Sim_List = Sim_Fn( n_years=nyears, n_stations=nstations, SpatialScale=SpatialScale, SD_O=SD_O,
 	                  rho=rho, beta1=beta1, beta0=beta0, sigma=sigma, Likelihood=Likelihood )

 DF = Sim_List[["DF"]]
 loc_xy_orig = loc_xy = Sim_List[["Loc"]]

 mesh = inla.mesh.create( loc_xy, refine=TRUE, extend=-0.5, cutoff=0.01 ) #
 spde = inla.spde2.matern( mesh, alpha=2 ) #nu=1 in matern simulation means alpha should be 2 here.
 #plot(mesh)
 #points(loc_xy_orig, pch=16, col="Steelblue")

 ggplot(DF, aes(Age, Simulated_Length)) + geom_point()
  #       facet_wrap(~Lake)

 d <- reshape2::melt(Sim_List[["eps_st"]]) %>%
        dplyr::mutate(x = rep(Sim_List[["Loc"]][,"x"], nyears),
        y = rep(Sim_List[["Loc"]][,"y"], nyears))

 ggplot(d, aes(x, y, col = value)) + geom_point() +
  facet_wrap(~L1) +
  scale_color_gradient2()

 #Extract the sparse matricies
 spdeMatrices = spde$param.inla[c("M0","M1","M2")]

 if(Likelihood=="Normal"){CTL <- 1}
 if(Likelihood=="Lognormal"){CTL <- 2}
 if(Likelihood=="Gamma"){CTL <- 3}

 # Build tagged list inputs
 Data = list("Nobs"=nrow(DF), "length_i"=DF$Simulated_Length, "age_i" = DF$Age,
             "lake_i" = DF$Lake - 1,
             "X_ij_omega"= model.matrix(~ -1 + DF$x1),
             "Nlakes" =  length(unique(DF$Lake)),
             "spdeMatrices" = spdeMatrices, "CTL" = CTL,
             "s_i" = DF$Lake-1,
             "t_i" = DF$Year-1,
             "x_s" = mesh$idx$loc-1,
             "n_t" = max(DF$Year) )

 Parameters = list("ln_global_omega" = log(beta0),
   	   					   "ln_global_linf" = log(Global_Linf),
                   "ln_sd_linf" = log(SD_linf),
  								 "global_tzero" = Global_t0,
                   "ln_sd_tzero" = log(SD_t0),
  								 "b_j_omega" = beta1,
  								 "eps_omega_st" = matrix(0,  nrow=mesh$n,ncol=Data$n_t ),
                   "eps_linf" = rep(0, Data$Nlakes ),
                   "eps_t0" = rep(0, Data$Nlakes ),
                   "ln_cv" = log(cv),
  								 "ln_kappa" = log(kappa),
	                 "ln_tau_O" = 1,
                   "rho_unscaled" = 2 * plogis(rho) - 1)

 Random = c("eps_linf", "eps_t0", "eps_omega_st")

 Data$ar1 <- 1L
 Map <- list()
 if (Data$ar1 == 0L) Map[["rho_unscaled"]] = factor(NA)
 Obj <- MakeADFun(data=Data, parameters=Parameters, random=Random,
                  hessian=FALSE, DLL=Version,map = Map)

 #Obj$fn( Obj$par )
 #Obj$gr( Obj$par )

 # Phased for a bit of a speed boost:
 # ---------------------------------
 # Phase 1: (fixed effects)

 Map <- list()
 Map[["ln_sd_linf"]] = factor(NA)
 Map[["eps_linf"]] = rep(factor(NA), length(Parameters$eps_linf))
 Map[["ln_sd_tzero"]] = factor(NA)
 Map[["eps_t0"]] = rep(factor(NA), length(Parameters$eps_t0))
 Map[["ln_kappa"]] = factor(NA)
 Map[["ln_tau_O"]] = factor(NA)
 Map[["rho_unscaled"]] = factor(NA)
 Map[["eps_omega_st"]] = factor(matrix(NA,  nrow=mesh$n,ncol=Data$n_t ))

 Obj <- MakeADFun(data=Data, parameters=Parameters, random=NULL,
                  hessian=FALSE, DLL=Version, map = Map)

 Opt = TMBhelper::Optimize(obj=Obj,
                           control=list(eval.max=1000, iter.max=1000),
                           getsd=T, newtonsteps=1, bias.correct=F)

 Parameters[["ln_global_omega"]] = set_par_value(Opt, "ln_global_omega")
 Parameters[["ln_global_linf"]] = set_par_value(Opt, "ln_global_linf")
 Parameters[["global_tzero"]] = set_par_value(Opt, "global_tzero")
 Parameters[["b_j_omega"]] = set_par_value(Opt, "b_j_omega")
 Parameters[["ln_cv"]] = set_par_value(Opt, "ln_cv")

 # ----------------------------------
 # Phase 2: (+ random effects)

 Map <- list()
 if (Data$ar1 == 0L) Map[["rho_unscaled"]] = factor(NA)

 Obj <- MakeADFun(data=Data, parameters=Parameters, random=Random,
                  hessian=FALSE, DLL=Version, map = Map)

 Opt = TMBhelper::Optimize(obj=Obj,
                           control=list(eval.max=1000, iter.max=1000),
                           getsd=T, newtonsteps=1, bias.correct=F)

 Opt

 SD = sdreport( Obj )
 final_gradient = Obj$gr( Opt$par )
 if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

 #Git yer report:
 Report = Obj$report()

 Estimates[i,1] = SD$value[["rho"]]
 Estimates[i,2] = SD$value[["SigmaO"]]
 Estimates[i,3] = exp(as.numeric(Opt$par["ln_kappa"]))
 Estimates[i,4] = exp(as.numeric(Opt$par["ln_cv"]))
 Estimates[i,5] = exp(as.numeric(Opt$par["ln_global_omega"]))
 Estimates[i,6] = exp(as.numeric(Opt$par["ln_global_linf"]))
 Estimates[i,7] = exp(as.numeric(Opt$par["ln_sd_linf"]))
 Estimates[i,8] = as.numeric(Opt$par["global_tzero"])
 Estimates[i,9] = exp(as.numeric(Opt$par["ln_sd_tzero"]))
 Estimates[i,10] = as.numeric(Opt$par["b_j_omega"])

 par(mfrow=c(4,3))
 hist(Estimates[,1], main="", xlab="Rho")
 abline(v=rho, lty=3, lwd=3, col="Steelblue")

 hist(Estimates[,2], main="", xlab="SD_O")
 abline(v=SD_O, lty=3, lwd=3, col="Steelblue")

 hist(Estimates[,3], main="", xlab="Kappa")
 abline(v=kappa, lty=3, lwd=3, col="Steelblue")

 hist(Estimates[,4], main="", xlab="CV")
 abline(v=cv, lty=3, lwd=3, col="Steelblue")

 hist(Estimates[,5], main="", xlab="Beta0" )
 abline(v=beta0, lty=3, lwd=3, col="Steelblue")

 hist(Estimates[,6], main="", xlab="Linf")
 abline(v=Global_Linf, lty=3, lwd=3, col="Steelblue")

 hist(Estimates[,7], main="", xlab="SD_Linf")
 abline(v=SD_linf, lty=3, lwd=3, col="Steelblue")

 hist(Estimates[,8], main="", xlab="Tzero")
 abline(v=Global_t0, lty=3, lwd=3, col="Steelblue")

 hist(Estimates[,9], main="", xlab="SD_t0")
 abline(v=SD_t0, lty=3, lwd=3, col="Steelblue")

 hist(Estimates[,10], main="", xlab="beta1")
 abline(v=beta1, lty=3, lwd=3, col="Steelblue")

 print(i)
}

DF$eps <- Report$eps_i
d <-   reshape2::melt(DF$eps) %>%
       mutate(WBID=DF$Lake,
              x = DF$Longitude,
              y = DF$Latitude,
              year = DF$Year)%>%
              group_by(WBID, year)

#plot eps_i through time x space
p <- group_by(d, year) %>%
     # mutate(eps = value - mean(value)) %>%
     mutate(eps = value) %>%
     ggplot(aes(x, y, col =eps )) + geom_point() +
     facet_wrap(~year) +
     xlab("Easting (km)") + ylab("Northing (km)") +
     scale_color_gradient2()
p

# ggsave("eps_i.png", p, dpi=600, width=9, height=8, units=c("in"))

png( file="100sim_NormalLike_50fish_50stations.png", width=8, height=8, res=500, units="in")
par(mfrow=c(4,3))
hist(Estimates[,1], main="", xlab="Rho", breaks=15)
abline(v=rho, lty=3, lwd=3, col="Steelblue")

hist(Estimates[,2], main="", xlab="SD_O", breaks=35)
abline(v=SD_O, lty=3, lwd=3, col="Steelblue")

hist(Estimates[,3], main="", xlab="Kappa", xlim=c(0.15, 0.45), breaks=50)
abline(v=kappa, lty=3, lwd=3, col="Steelblue")

hist(Estimates[,4], main="", xlab="CV")
abline(v=cv, lty=3, lwd=3, col="Steelblue")

hist(Estimates[,5], main="", xlab="Beta0" )
abline(v=beta0, lty=3, lwd=3, col="Steelblue")

hist(Estimates[,6], main="", xlab="Linf")
abline(v=Global_Linf, lty=3, lwd=3, col="Steelblue")

hist(Estimates[,7], main="", xlab="SD_Linf")
abline(v=SD_linf, lty=3, lwd=3, col="Steelblue")

hist(Estimates[,8], main="", xlab="Tzero")
abline(v=Global_t0, lty=3, lwd=3, col="Steelblue")

hist(Estimates[,9], main="", xlab="SD_t0")
abline(v=SD_t0, lty=3, lwd=3, col="Steelblue")

hist(Estimates[,10], main="", xlab="beta1")
abline(v=beta1, lty=3, lwd=3, col="Steelblue")

dev.off()
