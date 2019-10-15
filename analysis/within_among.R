#-------------------------------------------------------------------------------------------------------------------
#Nonspatial and Spatial-temporal von Bertalanffy with density dependent growth regression
#within-among parameterization
#Cahill Oct 2019
#-------------------------------------------------------------------------------------------------------------------

library(dplyr)
library(TMB)
library(INLA)
library(ggplot2)

data <- readRDS("C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/analysis/vB_analysis_august_2019_cahill.rds")
data <- as.data.frame(data)

#Set up the mean density for each lake for among lake effect:
newdata <-
  group_by(data, WBID) %>% #WBID = Waterbody ID
  summarize(lake_mean_dens = mean(wallEffDen)) %>%
  mutate(lake_mean_dens.std = arm::rescale(lake_mean_dens)) #Rescale it

ggplot(newdata, aes(y=lake_mean_dens.std, x=WBID)) + geom_point()

data <- left_join(data, newdata, by = "WBID")

#Set up density centered within each lake for random slopes:
newdata <- group_by(data, WBID) %>%
  mutate(lake_centered_dens = wallEffDen - mean(wallEffDen))

data$lake_centered_dens.std = arm::rescale(newdata$lake_centered_dens) #Rescale it

ggplot(data, aes(y=lake_centered_dens.std, x=WBID)) + geom_point()

#--------------
#Nonspatial model
#--------------

setwd("C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/executables")
VersionNonSpatial = "vb_nonspatial"
compile(paste0(VersionNonSpatial, ".cpp"))
dyn.load( dynlib(VersionNonSpatial) )

random_nonspatial = c("eps_linf", "eps_t0", "eps_omega", "eps_slope")

#Normal Likelihood                                                                                                                                          "ln_b_sex", "b_j_omega") )
CTL = 1

data_nonspatial = list("Nobs"=nrow(data), "length_i"=data$TL, "age_i" = data$Age,
                       "lake_i" = data$Lake - 1, "X_ij_omega"= model.matrix(~ -1 + data$lake_mean_dens.std), #among lake effect
                       "within_lake_i" = data$lake_centered_dens.std, #within lake effect
                       "sex_i" = data$SexCode,"Nlakes" = length(unique(data$Lake)),"CTL" = CTL,
                       "predTF_i"=rep(0, nrow(data)))


parameters_nonspatial = list("ln_global_omega" = log(14),
                             "ln_global_linf" = log(45), "ln_sd_linf" = log(7), "global_tzero" = -1,
                             "ln_sd_tzero" = log(3), "ln_b_sex" = log(4.760871), "b_j_omega" = rep(0, ncol(data_nonspatial$X_ij_omega)),
                             "mu_slope" = 0, #average slope term
                             "eps_omega" = rep(0, data_nonspatial$Nlakes ),
                             "eps_linf" = rep(0, data_nonspatial$Nlakes ),
                             "eps_t0" = rep(0, data_nonspatial$Nlakes ),
                             "eps_slope" = rep(0, data_nonspatial$Nlakes ),
                             "ln_cv" = log(0.2), "ln_sd_omega" = log(4.5),
                             "ln_sd_slope"=0)


obj_nonspatial <- MakeADFun(data=data_nonspatial, parameters=parameters_nonspatial,
                            random=random_nonspatial, hessian=FALSE, DLL=VersionNonSpatial)

opt_nonspatial = TMBhelper::Optimize(obj=obj_nonspatial,
                                     control=list(eval.max=1000, iter.max=1000),
                                     getsd=T, newtonsteps=1, bias.correct=F)

SD = sdreport( obj_nonspatial )
final_gradient = obj_nonspatial$gr( opt_nonspatial$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

ParHat_nonspatial = as.list( opt_nonspatial$SD, "Estimate" )
SEHat_nonspatial  = as.list( opt_nonspatial$SD, "Std. Error" )

cbind("slope"=(ParHat_nonspatial$eps_slope + ParHat_nonspatial$mu_slope),"se"=SEHat_nonspatial$eps_slope)

ParHat_nonspatial$mu_slope
SEHat_nonspatial$mu_slope

ParHat_nonspatial$b_j_omega #first value is among-lake coefficient
SEHat_nonspatial$b_j_omega


#----------
# spatial model:
#----------

loc_xy <- unique(data[ ,c("X_TTM_c","Y_TTM_c") ] )
loc_xy <- loc_xy/1000 #Put distance in kms

#mesh = inla.mesh.2d(loc=loc_xy, max.edge=c(62,1000))
mesh = inla.mesh.create( loc_xy, refine=TRUE, extend=-0.5, cutoff=0.01 )

mesh$n

#Create the inputs for spatial model
spde = inla.spde2.matern( mesh )
spdeMatrices = spde$param.inla[c("M0","M1","M2")]
random_spatial = c("eps_omega_st", "eps_linf", "eps_t0")

VersionSpatial = "spdeXAR1_v3"
compile(paste0(VersionSpatial, ".cpp"))
dyn.load( dynlib(VersionSpatial) )

#Normal Likelihood--spatial                                                                                                                                          "ln_b_sex", "b_j_omega") )
CTL = 1

data_spatial = list("Nobs" = nrow(data), "length_i"=data$TL, "age_i" = data$Age,
                    "lake_i" = data$Lake - 1, "sex_i" = data$SexCode,
                    "X_ij_omega"= model.matrix(~ -1 + data$lake_mean_dens.std), #among lake effect
                    "within_lake_i" = data$lake_centered_dens.std, #within lake effect
                    "Nlakes" = length(unique(data$Lake)), "spdeMatrices" = spdeMatrices,
                    "s_i" = data$Lake-1, "t_i" = data$Year-1, "CTL" = CTL, "predTF_i"=rep(0, nrow(data)))

parameters_spatial = list("ln_global_omega" = log(13.28954),
                          "ln_global_linf" = log(55),
                          "ln_sd_linf" = log(7.275452),
                          "ln_sd_slope" = 2.0,
                          "global_tzero" = -1.260783461,
                          "ln_sd_tzero" = log(0.538755),
                          "ln_b_sex" = log(4.760871),
                          "b_j_omega" = rep(0, ncol(data_spatial$X_ij_omega)),
                          "mu_slope" = 0,
                          "eps_omega_st" = matrix(0,  nrow=mesh$n,ncol=max(data$Year) ),
                          "eps_linf" = rep(0,  length(unique(data$Lake))),
                          "eps_t0" = rep(0, length(unique(data$Lake))),
                          "eps_slope" = rep(0, length(unique(data$Lake))),
                          "ln_cv" = -2.5043055,
                          "ln_kappa" = -2.7,
                          "ln_tau_O" = 0,
                          "rho" = 0.9 )

random_spatial = c("eps_omega_st", "eps_linf", "eps_t0", "eps_slope")

obj_spatial = MakeADFun(data=data_spatial, parameters=parameters_spatial,
                        random=random_spatial, hessian=FALSE, DLL=VersionSpatial)

#about ten minutes to run on my machine:
opt_spatial  = TMBhelper::Optimize(obj=obj_spatial,
                                   control=list(eval.max=1000, iter.max=1000),
                                   getsd=T, newtonsteps=1, bias.correct=T,
                                   lower=c(rep(-Inf,12),-0.999), upper=c(rep(Inf,12),0.999))

SD = sdreport( obj_spatial )
final_gradient = obj_spatial$gr( opt_spatial$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

ParHat = as.list( opt_spatial$SD, "Estimate" )
SEHat  = as.list( opt_spatial$SD, "Std. Error" )
ParHat$eps_slope
SEHat$eps_slope

cbind("slope"=(ParHat$eps_slope + ParHat$mu_slope),"se"=SEHat$eps_slope)
opt_spatial

ParHat$mu_slope
SEHat$mu_slope

ParHat_nonspatial$b_j_omega #first value is among-lake coefficient
SEHat_nonspatial$b_j_omega
