library(INLA)
library(TMB)
library(RandomFields)
library(raster)
library(RANN)
library(ggplot2)
library(dplyr)
library(truncnorm)

rf_sim <- function(model, x, y) {
  set.seed(sample.int(1e5L, 1L))
  suppressMessages(
    RandomFields::RFsimulate(model = model, x = x, y = y)$variable1
  )
}

rlnorm0 <- function(n = 1, mean, coeffOfVar)
{
   sigma <- sqrt(log(coeffOfVar^2 + 1))
   mu <- log(mean) - sigma^2 / 2
   rlnorm(n, mu, sigma)
}

Sim_Fn <- function( n_years=nyears, n_stations=nstations, SpatialScale=SpatialScale, SD_O=SD_O,
 	                  rho=rho, betas=betas, beta0=beta0, cv=cv, Likelihood="LogNormal", Design="Balanced"){

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

  #Simulate covariates:
  X_ij <- array(runif(n=n_years*n_stations*length(betas), min=-1, max=1), dim=c(nrow=n_stations, ncol=n_years, length(betas)))

  Eta_it   = array(NA, dim=c(n_stations,n_years))
  for(t in 1:n_years){
    Eta_it[,t] =  beta0  + betas[1]*X_ij[ ,t, 1] + betas[2]*X_ij[ ,t, 2] + betas[3]*X_ij[ ,t, 3] + eps_st[[t]] #fixed effects + Temporally evolving spatial field
  }

  #Exclude exceedingly low simulated omegas (< 40 mm)-could instead use log-linear regression
  Eta_it[which(Eta_it < 4.0)] <- NA

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

      lpreds <-   Linf*(1-exp(-(Omega/Linf) * (Ages - T0 ))) #VB predictions

      if(!any(is.na(lpreds))){
       if(Likelihood=="Normal"){Simulated_Length = rtruncnorm(Nfish, mean=lpreds, sd=cv*lpreds, a=1, b=Inf)}
       if(Likelihood=="Lognormal"){Simulated_Length = rlnorm0(Nfish, lpreds, cv)}
       if(Likelihood=="Gamma"){Simulated_Length = rgamma(Nfish, shape=1/cv^2, scale=lpreds*cv^2)}
      } else {
       Simulated_Length <- lpreds
      }

      Tmp = data.frame("Lake"=rep(s, Nfish) , "Year"=rep(t, Nfish), "x1"=rep(X_ij[s,t,1], Nfish),
                       "x2"=rep(X_ij[s,t,2], Nfish), "x3"=rep(X_ij[s,t,3], Nfish), "Age" = Ages,
                       "Simulated_Length"= Simulated_Length)
      DF = rbind(DF, Tmp)
  }}

  DF = cbind( DF, 'Longitude'=Loc[DF[,'Lake'],1], 'Latitude'=Loc[DF[,'Lake'],2] )
  DF = data.frame(DF, row.names=NULL)

  if(Design=="Unbalanced"){
   lakes_to_keep <- sample(unique(DF$Lake), lakes_good_data*length(unique(DF$Lake)), replace=FALSE)

   Which2Drop = setdiff(1:nrow(DF),lakes_to_keep)

   DF[Which2Drop,'Simulated_Length'] = NA
   Which2Keep = sample(1:nrow(DF), size=2000*Sim_List$n_years, replace=FALSE)
   Which2Drop = setdiff(1:nrow(DF),Which2Keep)
   message("Not yet implemented"); break
  }

  # Return stuff
  Sim_List = list("DF"=DF, "Loc"=Loc, "Eta_it"=Eta_it,
                  "eps_st"=eps_st, "n_years"=n_years,
                  "n_stations"=n_stations)
  Sim_List[["Parameters"]] = c('SpatialScale'=SpatialScale, 'SigmaO'=SD_O,
                               'rho'=rho, "beta0"=beta0, "beta_ij"=betas,
                               'Linf' = Global_Linf, "t0" = Global_t0,
                               'SD_t0' = SD_t0,'SD_linf'=SD_linf)
  Sim_List
}

#--------------------
#Run the simulation:
#---------------------

#Simulation parameters:
Age_Range    <- 0:25
Global_Linf  <- 55.70
SD_linf      <- 7.41
Global_t0    <- -1
SD_t0        <- 0.3
betas        <- c(-1,0,1) #runif(1, -0.25, 0.25)
beta0        <- 14.79 #runif(1, 14, 16)
kappa        <- 0.5   #c(0.2, 0.5, 0.75)
SpatialScale <- 1/kappa
Range        <- sqrt(8) / kappa
SD_O         <- 4.66
rho          <- 0.5
cv           <-  0.069
Nfish        <- 50

Likelihood=c("Normal", "Lognormal", "Gamma")
Likelihood = Likelihood[1]

Design = c("Balanced", "Unbalanced")
Design = Design[1]
lakes_good_data <- 0.1 #number of lakes with good data
percent_of_data  <- 0.35 #percent of total dataset represented by lakes_good_data

setwd("sim/")
VersionSpatial = "vb_spdeXar1"
VersionNonSpatial = "vb_nonspatial"

# Compile
compile( paste0(VersionSpatial,".cpp") )
compile( paste0(VersionNonSpatial,".cpp") )

seed <-  sample.int(1e6, 1)
set.seed( seed )
#Run simulation:

N_years <- c(4,7,10,15,20)
N_lakes <- c(25,50,80,100)
Nsim <- 30
ptm <- proc.time()
for(nyears in unique(N_years)){
  for(Nlakes in unique(N_lakes)){
    replicate=1
    while(replicate  < Nsim ){
    Sim_List = Sim_Fn( n_years=nyears, n_stations=Nlakes, SpatialScale=SpatialScale, SD_O=SD_O,
 	                     rho=rho, betas=betas, beta0=beta0, cv=cv, Likelihood=Likelihood, Design=Design )

    DF = Sim_List[["DF"]]
    loc_xy_orig = loc_xy = unique(DF[,c("Longitude", "Latitude")])

    mesh = inla.mesh.create( loc_xy, refine=TRUE, extend=-0.5, cutoff=0.01 )
    spde = inla.spde2.matern( mesh, alpha=2 )

    #------Plots--------
    #plot(mesh)
    #points(loc_xy_orig, pch=16, col="Steelblue")
    #
    #ggplot(DF, aes(Age, Simulated_Length)) + geom_point() + facet_wrap(~Lake)
    #
    # d <- reshape2::melt(Sim_List[["eps_st"]]) %>%
    #        dplyr::mutate(x = rep(Sim_List[["Loc"]][,"x"], nyears),
    #        y = rep(Sim_List[["Loc"]][,"y"], nyears))
    #
    # ggplot(d, aes(x, y, col = value)) + geom_point() +
    #  facet_wrap(~L1) +
    #  scale_color_gradient2()
    #-------------------

    #Extract the sparse matricies
    spdeMatrices = spde$param.inla[c("M0","M1","M2")]

    if(Likelihood=="Normal"){CTL <- 1}
    if(Likelihood=="Lognormal"){CTL <- 2}
    if(Likelihood=="Gamma"){CTL <- 3}

    # Build tagged list inputs
    data_spatial = list("Nobs"=nrow(DF), "length_i"=DF$Simulated_Length, "age_i" = DF$Age,
                        "lake_i" = DF$Lake - 1,
                        "X_ij_omega"= model.matrix(~ -1 + DF$x1 + DF$x2 + DF$x3),
                        "Nlakes" =  length(unique(DF$Lake)),
                        "spdeMatrices" = spdeMatrices, "CTL" = CTL,
                        "s_i" = DF$Lake-1,
                        "t_i" = DF$Year-1,
                        "n_t" = max(DF$Year) )

     data_nonspatial = list("Nobs"=nrow(DF), "length_i"=DF$Simulated_Length, "age_i" = DF$Age,
                            "lake_i" = DF$Lake - 1,
                            "X_ij_omega"= model.matrix(~ -1 + DF$x1 + DF$x2 + DF$x3),
                            "Nlakes" =  length(unique(DF$Lake)),"CTL" = CTL )

     parameters_spatial = list("ln_global_omega" = log(beta0),
   	       					           "ln_global_linf" = log(Global_Linf),
                               "ln_sd_linf" = log(SD_linf),
  					    			         "global_tzero" = Global_t0,
                               "ln_sd_tzero" = log(SD_t0),
  								             "b_j_omega" = betas,
  								             "eps_omega_st" = matrix(0,  nrow=mesh$n,ncol=data_spatial$n_t ),
                               "eps_linf" = rep(0, data_spatial$Nlakes ),
                               "eps_t0" = rep(0, data_spatial$Nlakes ),
                               "ln_cv" = log(cv),
  								             "ln_kappa" = log(kappa),
	                             "ln_tau_O" =  -2,
                               "rho" = rho )

      parameters_nonspatial = list("ln_global_omega" = log(beta0),
                                   "ln_sd_omega" = log(SD_O),
   	   		          	             "ln_global_linf" = log(Global_Linf),
                                   "ln_sd_linf" = log(SD_linf),
  								                 "global_tzero" = Global_t0,
                                   "ln_sd_tzero" = log(SD_t0),
  								                 "b_j_omega" = betas,
  								                 "eps_omega" = rep(0, data_nonspatial$Nlakes ),
                                   "eps_linf" = rep(0, data_nonspatial$Nlakes ),
                                   "eps_t0" = rep(0, data_nonspatial$Nlakes ),
                                   "ln_cv" = log(cv) )

      random_spatial = c("eps_linf", "eps_t0", "eps_omega_st")
      random_nonspatial = c("eps_linf", "eps_t0", "eps_omega")

      dyn.load( dynlib(VersionSpatial) )
      obj_spatial <- MakeADFun(data=data_spatial, parameters=parameters_spatial,
                               random=random_spatial, hessian=FALSE, DLL=VersionSpatial)

      opt_spatial = tryCatch(TMBhelper::Optimize(obj=obj_spatial,
                             control=list(eval.max=1000, iter.max=1000),
                             getsd=T, newtonsteps=1, bias.correct=T,
                             lower=c(rep(-Inf,11),-0.999), upper=c(rep(Inf,11),0.999)),
                             error = function(e) e)

      dyn.load( dynlib(VersionNonSpatial) )
      obj_nonspatial <- MakeADFun(data=data_nonspatial, parameters=parameters_nonspatial,
                                  random=random_nonspatial, hessian=FALSE, DLL=VersionNonSpatial)

      opt_nonspatial = tryCatch(TMBhelper::Optimize(obj=obj_nonspatial,
                                control=list(eval.max=1000, iter.max=1000),
                                getsd=T, newtonsteps=1, bias.correct=F),
                                error = function(e) e)

      #if the estimation behaves, save results & advance loop
      if(!inherits(opt_spatial, "error") && !inherits(opt_nonspatial, "error")){

        SD = sdreport( obj_spatial )
        SD_nonspatial = sdreport( obj_nonspatial )

        gradients <- c(c(obj_spatial$gr( opt_spatial$par ), obj_nonspatial$gr( opt_nonspatial$par )))
        if( any(abs(gradients)>0.0001) | SD$pdHess==FALSE | SD_nonspatial$pdHess==FALSE ) stop("Not converged")

        sim_rep <- list()
        sim_rep[["sim_data"]] = Sim_List

        spatial_summary = summary(SD)
        spatial_report = obj_spatial$report()
        opt_spatial[["summary"]] = spatial_summary
        opt_spatial[["report"]] = spatial_report
        sim_rep[["opt_spatial"]] = opt_spatial

        nonspatial_summary = summary(SD_nonspatial)
        nonspatial_report = obj_nonspatial$report()
        opt_nonspatial[["summary"]] = nonspatial_summary
        opt_nonspatial[["report"]] = nonspatial_report
        sim_rep[["opt_nonspatial"]] = opt_nonspatial

        file_name <- paste(paste0(paste0(nyears, "Years"), Nlakes, "Lakes"), paste(replicate,"Replicate.RData",sep=""), sep="_")
        save(sim_rep, file=file_name)
        print(paste(paste("Simulation", replicate, sep=" "), "Complete", sep=" "))
        replicate <- replicate + 1
      } #try
    } #while
    print(paste(paste("Nlakes", Nlakes, sep=" "), "Complete", sep=" "))
  } #nlakes
} #nyears
proc.time() - ptm





sim_results <- array(NA, dim=c(2, Nsim, 10, 3), dimnames=list(c("Spatial","Nonspatial"), NULL,
                     names(opt_spatial$par),
                     c("10%","50%","90%")))

for(replicate in 1:Nsim){
  load( paste(replicate,"Simulation.RData",sep="") )

  sim_results["Spatial", replicate, names(opt_spatial$par),
              c("10%","50%","90%")] = sim_rep$opt_spatial$summary[c(names(opt_spatial$par)),'Estimate']%o%rep(1,3) +
              sim_rep$opt_spatial$summary[c(names(opt_spatial$par)), 'Std. Error']%o%qnorm(c(0.1,0.5,0.9))

  sim_results["Nonspatial", replicate, c("ln_global_omega", "ln_global_linf", "ln_sd_linf", "global_tzero", "ln_sd_tzero", "b_j_omega", "ln_cv"),
               c("10%","50%","90%")] = sim_rep$opt_nonspatial$summary[c("ln_global_omega", "ln_global_linf", "ln_sd_linf", "global_tzero", "ln_sd_tzero", "b_j_omega", "ln_cv"),
              'Estimate']%o%rep(1,3) + sim_rep$opt_nonspatial$summary[c("ln_global_omega", "ln_global_linf", "ln_sd_linf", "global_tzero", "ln_sd_tzero", "b_j_omega", "ln_cv"), 'Std. Error']%o%qnorm(c(0.1,0.5,0.9))
}

par(mfrow=c(2,1))
boxplot(sim_results["Spatial", , "b_j_omega", 2])
boxplot(sim_results["Nonspatial", , "b_j_omega", 2], add=T)

abline(h=beta1, lty=3)

d <-   reshape2::melt(sim_results[c("Spatial", "Nonspatial"), , "b_j_omega", 2])
p <- ggplot(d, aes(Var1, value))
p <- p + geom_boxplot()
p <- p +  geom_hline(aes(yintercept=beta1), colour="#990000")
p



save( tmp, file=paste(RepFile,"Save_Spatial.RData",sep=""))
load("2Simulation.RData")
tmp
names(tmp[["spatial_summary"]])


opt_spatial$summary["rho", "Estimate"]
sim_rep$opt_spatial$summary["rho", "Estimate"]
sim_rep$opt_spatial$summary["rho", "Std. Error"]

sim_rep$opt_nonspatial$summary["b_j_omega", "Estimate"]
sim_rep$opt_spatial$summary["b_j_omega", "Estimate"]

#Unbalanced vs. balanced design

#Non-spatial model for simulation
#4,8,12,16,20 years
#25,50,75,100 lakes

If possible, estimate non-spatial model in last loop,
then do data processing after below.




vector <- as.data.frame(unlist(Estimates))
vector$eps_t0

list.condition <- sapply(Estimates, function(x) class(x)=="eps_omega_st")
output.list  <- Estimates[list.condition]

DF$eps <- spatial_report$eps_i
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

png( file="100sim_Gamma_50fish_80stations_4Years.png", width=8, height=8, res=500, units="in")
par(mfrow=c(4,3))
hist(Estimates[[1:Nsim]], main="", xlab="Rho")
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

dev.off()

d2 <- reshape2::melt(Sim_List[["eps_st"]]) %>%
         dplyr::mutate(x = rep(Sim_List[["Loc"]][,"x"], nyears),
         y = rep(Sim_List[["Loc"]][,"y"], nyears))

  ggplot(d2, aes(x, y, col = value)) + geom_point() +
   facet_wrap(~L1) +
   scale_color_gradient2()







  # Build tagged list inputs
 Data = list("Nobs"=nrow(DF), "length_i"=DF$Simulated_Length, "age_i" = DF$Age,
             "lake_i" = DF$Lake - 1,
             "X_ij_omega"= model.matrix(~ -1 + DF$x1),
             "Nlakes" =  length(unique(DF$Lake)),"CTL" = CTL
             )

 Parameters = list("ln_global_omega" = log(beta0),
                   "ln_sd_omega" = log(SD_O),
   	   					   "ln_global_linf" = log(Global_Linf),
                   "ln_sd_linf" = log(SD_linf),
  								 "global_tzero" = Global_t0,
                   "ln_sd_tzero" = log(SD_t0),
  								 "b_j_omega" = beta1,
  								 "eps_omega" = rep(0, Data$Nlakes ),
                   "eps_linf" = rep(0, Data$Nlakes ),
                   "eps_t0" = rep(0, Data$Nlakes ),
                   "ln_cv" = log(cv)
                   )

 Random = c("eps_linf", "eps_t0", "eps_omega")

 Version = "vb_nonspatial"


# Compile
compile( paste0(Version,".cpp") )
dyn.load( dynlib(Version) )

 obj_nonspatial <- MakeADFun(data=Data, parameters=Parameters, random=Random,
                   hessian=FALSE, DLL=Version)

 opt_nonspatial = tryCatch(TMBhelper::Optimize(obj=obj_nonspatial,
                           control=list(eval.max=1000, iter.max=1000),
                           getsd=T, newtonsteps=1, bias.correct=T),
                           error = function(e) e)
