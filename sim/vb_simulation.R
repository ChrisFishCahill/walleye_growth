library(INLA)
library(TMB)
library(RandomFields)
library(raster)
library(RANN)
library(ggplot2)
library(dplyr)

rf_sim <- function(model, x, y) {
  # set.seed(sample.int(1e5L, 1L))
  suppressMessages(
    RandomFields::RFsimulate(model = model, x = x, y = y)$variable1
  )
}

Sim_Fn <-
 function( n_years, n_stations=100, SpatialScale=SpatialScale, SD_O=0.5, rho=0.8, beta1=-0.25, beta0=100, sigma=0.05,title=i ){

  # Spatial model
  Loc = cbind( "x"=runif(n_stations, min=0,max=1), "y"=runif(n_stations, min=0,max=1) )
  model_O <- RandomFields::RMmatern(nu=1, var=SD_O^2, scale = SpatialScale ) #SpatialScale = 1 / kappa

  # Simulate Omega
  u.t <- array(NA, dim=c(n_stations, n_years))
  #Simulate the temporally evolving spatial field as per Zuur et al. 2017 p. 321:
  for(i in 1:n_years){
   RFoptions(seed=i * 42, printlevel=0) #This was necessary post RF package update Feb 2019
   u.t[,i] <- rf_sim(model = model_O, x=Loc[,'x'], y=Loc[,'y'])
  }
  # print(u.t)

  v <- u.t
  for(i in 2:n_years){v[,i] <- rho*v[,i-1] + u.t[,i]}

  #Simulate a covariate x:
  x <- matrix(runif(n=n_years*n_stations, min=0, max=1.5), nrow=n_stations, ncol=n_years)

  # Calculate Eta_it --> this is spatial temporal field for omega
  Eta_it   = array(NA, dim=c(n_stations,n_years))
  for(t in 1:n_years){
    Eta_it[,t] =  beta0  + beta1*x[ ,t] + v[,t] #mu(i,t) = fixed effects + Temporally evolving spatial field
  }
  #
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

      Tmp = data.frame("Lake"=rep(s, Nfish) , "Year"=rep(t, Nfish), "x1"=rep(x[s,t], Nfish), "Age" = Ages,
                       "Simulated_Length"=rlnorm(Nfish, log(lpreds), cv))
     		                                 #rnorm(Nfish, mean=lpreds, sd=sigmas))
      DF = rbind(DF, Tmp)
    }}

  DF = cbind( DF, 'Longitude'=Loc[DF[,'Lake'],1], 'Latitude'=Loc[DF[,'Lake'],2] )
  DF = data.frame(DF, row.names=NULL)
  par(mfrow=c(1,1))

  # Return stuff
  Sim_List = list("DF"=DF, "Loc"=Loc, "Omega"=Omega, "v"=v, "Eta_it"=Eta_it, "n_years"=n_years, "n_stations"=n_stations)
  Sim_List[["Parameters"]] = c('SpatialScale'=SpatialScale, 'SigmaO'=SD_O,  'rho'=rho, "beta0"=beta0, "beta1"=beta1, "sigma"=sigmas,
                               'Linf' = Global_Linf, "t0" = Global_t0, 'SD_t0' = SD_t0, 'SD_linf'=SD_linf)
  return(Sim_List)
}

#---------------------------------------------------------------------------------------------------------------------------------------
#Estimate the model:
#---------------------------------------------------------------------------------------------------------------------------------------

#Simulation parameters:
SpatialScale <- 1.1
SD_O         <- 0.5
rho          <- 0.4
beta1        <- -2
beta0        <- 150/10

SD_linf     <- 7.22
Global_Linf <- 550/10

SD_t0       <- 0.54
Global_t0   <- -1

cv          <- 0.08

Age_Range = 2:25

Nfish       <- 100
nstations   <- 50
nyears      <- 20

beta0       <- 12
setwd("sim/")
Version = "vb_sim_estimation"

# Compile
compile( paste0(Version,".cpp") )
dyn.load( dynlib(Version) )

set.seed( 34 )

#Make fake data:
Sim_List = Sim_Fn( n_years=nyears, n_stations=nstations, SpatialScale=SpatialScale, SD_O=SD_O,
	                 rho=rho, beta1=beta1, beta0=beta0, sigma=sigma, title=i )

DF = Sim_List[["DF"]]
loc_xy_orig = loc_xy = Sim_List[["Loc"]]

# Build SPDE object using INLA
mesh = inla.mesh.create( loc_xy, refine=TRUE, extend=-0.5, cutoff=0.01 ) #
spde = inla.spde2.matern( mesh, alpha=2 ) #nu=1 in matern simulation means alpha should be 2 here.
plot(mesh)
points(loc_xy_orig, pch=16, col="Steelblue")

ggplot(DF, aes(Age, Simulated_Length)) + geom_point() +
       facet_wrap(~Lake)

d <- reshape2::melt(Sim_List[["v"]]) %>%
       dplyr::mutate(x = rep(Sim_List[["Loc"]][,"x"], ncol(Sim_List[["v"]])),
       y = rep(Sim_List[["Loc"]][,"y"], ncol(Sim_List[["v"]])))

ggplot(d, aes(x, y, col = value)) + geom_point() +
 facet_wrap(~Var2) +
 scale_color_gradient2()

#Extract the sparse matricies
spdeMatrices = spde$param.inla[c("M0","M1","M2")]

# Build tagged list inputs
Data = list("Nobs"=nrow(DF), "length_i"=DF$Simulated_Length, "age_i" = DF$Age,
            "lake_i" = DF$Lake - 1,
            "X_ij_omega"= model.matrix(~ -1 + DF$x1),
            "Nlakes" =  length(unique(DF$Lake)),
            "spdeMatrices" = spdeMatrices, "CTL" = 2,
            "s_i" = DF$Lake-1,
            "t_i" = DF$Year-1,
            "t_prev_i" = DF$Year-2,
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
  								"ln_kappa" = 1.5,
	                "ln_tau_O" = 1,
                  "rho_unscaled" = 2 * plogis(rho) - 1)

Random = c("eps_linf", "eps_t0", "eps_omega_st")

# Map <- list()
# Map[["ln_kappa"]] = factor(NA)
# Map[["ln_tau_O"]] = factor(NA)
# Map[["rho_unscaled"]] = factor(NA)
# Map[["eps_omega_st"]] = factor(matrix(NA,  nrow=mesh$n,ncol=Data$n_t ))

Obj <- MakeADFun(data=Data, parameters=Parameters, random=Random,
	               hessian=FALSE, DLL=Version, map=NULL)

Obj$fn( Obj$par )
Obj$gr( Obj$par )

#4 minutes on my machine
Opt = TMBhelper::Optimize( obj=Obj,
                           control=list(eval.max=1000, iter.max=1000),
                           getsd=T, newtonsteps=1, bias.correct=F)
Opt

SD = sdreport( Obj )
final_gradient = Obj$gr( Opt$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

#Git yer report:
Report = Obj$report()

Report$rho
rho

Report$Range
1/SpatialScale

DF$eps <- Report$eps_i
d <-   reshape2::melt(DF$eps) %>%
       mutate(WBID=DF$Lake,
              x = DF$Longitude,
              y = DF$Latitude,
              year = DF$Year)%>%
              group_by(WBID, year)

#plot mean-centered eps_i through time x space
p <- group_by(d, year) %>%
     mutate(eps = value - mean(value)) %>%
     ggplot(aes(x, y, col =eps )) + geom_point() +
     facet_wrap(~year) +
     xlab("Easting (km)") + ylab("Northing (km)") +
     scale_color_gradient2(limits = c(-0.01, 0.01))
p

#Well, that ain't proper..
# ggsave("eps_i.png", p, dpi=600, width=9, height=8, units=c("in"))

