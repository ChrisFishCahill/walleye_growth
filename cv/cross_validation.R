#-------------------------------------------------------------------------------------------------------------------
#Spatial-temporal von Bertalanffy with density dependent growth regression
#Galluci and Quinn (1979) Parameterization
#Cross-validation
#
#Model structure:
#length_pred_i = linfinity*(1-exp(-( omega/linfinity) * (age_i (i) - t0 )))
#
# Where:
# linfinity = linfinity_global + sex_effect_i + eps_linf
# t0  = t0_global + eps_t0
# omega   = omega_global + eta_fixed_i + eps_omega_st
#
#eta_fixed_i = Effective Intraspecific (Walleye) Density +
#              Effective Interspecific Density +
#              Interaction b/n Intra & Inter Density +
#              Growing Degree Days
#
# Probability of Random Coefficients:
# epsilon_linf ~ N(0, ln_sd_linf)
# epsilon_t0 ~ N(0, ln_sd_tzero)
# eps_omega_st = rho*eps_omega_s,t-1 + u_st
# where u_st ~ N(0, SIMGA) --> SIGMA is estimated as per INLA approach
#
#Likelihoods:
# if(CTL == 1 == Normal)
# if(CTL == 2 == Lognormal)
# if(CTL == 3 == Gamma)
#
#Cahill May 2019
#-------------------------------------------------------------------------------------------------------------------

library(dplyr)
library(plot3D)
library(plotly)
library(viridis)
library(magick)
library(lattice)
library(ggmap)
library(TMB)
library(INLA)

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

#Read in the data and build a mesh
data <- readRDS("vB_analysis_march_2019_cahill.rds")
data <- as.data.frame(data)

loc_xy <- unique(data[ ,c("X_TTM_c","Y_TTM_c") ] )
loc_xy <- loc_xy/1000 #Put distance in kms

mesh = inla.mesh.create( loc_xy, refine=TRUE, extend=-0.5, cutoff=0.01 )

# png(file="plots/Mesh.png",width=9.50,height=7.00,units="in",res=600)
plot(mesh)
points(loc_xy, col="Steelblue", pch=1)
# dev.off()

#This mesh is ugly as sin, but gives the same answers as the more complex meshes
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

# Extra meshes:

# RangeGuess <- 30    #~ 1/3 study area
# MaxEdge    <- RangeGuess/ 5 #as per https://haakonbakka.bitbucket.io/btopic104.html
# #
# bound.outer <- diff(range(loc_xy[,1]))/3
# RangeGuess <- 50

# mesh2 = inla.mesh.2d(loc=loc_xy, max.edge = c(2,5)*RangeGuess )
#
# mesh3 = inla.mesh.2d(loc=loc_xy, max.edge = c(3,5)*RangeGuess )
#
# mesh4 = inla.mesh.2d(loc=loc_xy, max.edge = c(1,5)*RangeGuess )
#
# mesh4 = inla.mesh.2d(loc=loc_xy, max.edge = c(1,5)*RangeGuess,
#                      cutoff = MaxEdge,
#                      offset = c(MaxEdge, bound.outer))

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#Load the models:
VersionNonSpatial = "vb_nonspatial"
compile(paste0(VersionNonSpatial, ".cpp"))
dyn.load( dynlib(VersionNonSpatial) )

VersionSpatial = "spdeXAR1_v3"
compile(paste0(VersionSpatial, ".cpp"))
dyn.load( dynlib(VersionSpatial) )

random_nonspatial = c("eps_linf", "eps_t0", "eps_omega")
random_spatial = c("eps_omega_st", "eps_linf", "eps_t0")

CTL = 1 #1=normal, 2=lognormal, 3=gamma
K = 10

#Vectors for pred negloglike:
nonspatial <- spatial <- rep(NA, K)
#-----------------------------------------------------------------------------------------------
#leave K/100 % of lakes out cross validation:
set.seed(99)
lakes_to_predict = sample(unique(data$Lake), replace=F)
Nlakes <- length(unique(data$Lake))
by <- round((K/100)*Nlakes)

to <- seq(from=1, to=Nlakes-by, by=by)
nonspatial <- spatial <- rep(NA, K)

#Create the SPDE inputs for TMB
spde = inla.spde2.matern( mesh )
spdeMatrices = spde$param.inla[c("M0","M1","M2")]

for(k in to){
  #eight lakes for first K-1 partitions
  if(k !=max(to)) Partition_i = ifelse(data$Lake %in% lakes_to_predict[k:(k + by - 1)],1,0)
  #nine total lakes in partition == K
  if(k == max(to)) Partition_i = ifelse(data$Lake %in% lakes_to_predict[k:(k + by)],1,0)

  data_nonspatial = list("Nobs"=nrow(data), "length_i"=data$TL, "age_i" = data$Age,
                         "lake_i" = data$Lake - 1, "X_ij_omega"= model.matrix(~ -1 + data$wallEffDen.Std +
                                                                                data$compEffDen.Std +
                                                                                data$GDD.Std  +
                                                                                data$wallEffDen.Std:data$compEffDen.Std),
                         "sex_i" = data$SexCode,"Nlakes" = length(unique(data$Lake)),"CTL" = CTL,
                         "predTF_i"=Partition_i)

  data_spatial = list("Nobs"      = nrow(data), "length_i"=data$TL, "age_i" = data$Age,
                      "lake_i"    = data$Lake - 1, "sex_i" = data$SexCode,
                      "X_ij_omega"= model.matrix(~ -1 + data$wallEffDen.Std + #.Std --> already standardized
                                                   data$compEffDen.Std +
                                                   data$GDD.Std  +
                                                   data$wallEffDen.Std:data$compEffDen.Std),
                      "Nlakes" = length(unique(data$Lake)),
                      "spdeMatrices" = spdeMatrices,
                      "s_i" = data$Lake-1,
                      "t_i" = data$Year-1,
                      "CTL" = CTL,
                      "predTF_i"=Partition_i)

  parameters_nonspatial = list("ln_global_omega" = log(14),
                               "ln_global_linf" = log(45), "ln_sd_linf" = log(7), "global_tzero" = -1,
                               "ln_sd_tzero" = log(3), "ln_b_sex" = log(4.760871), "b_j_omega" = rep(0, ncol(data_nonspatial$X_ij_omega)),
                               "eps_omega" = rep(0, data_nonspatial$Nlakes ), "eps_linf" = rep(0, data_nonspatial$Nlakes ),
                               "eps_t0" = rep(0, data_nonspatial$Nlakes ),"ln_cv" = log(0.2), "ln_sd_omega" = log(4.5) )

  parameters_spatial = list("ln_global_omega" = log(13.28954),
                            "ln_global_linf" = log(55),
                            "ln_sd_linf" = log(7.275452),
                            "global_tzero" = -1.260783461,
                            "ln_sd_tzero" = log(0.538755),
                            "ln_b_sex" = log(4.760871),
                            "b_j_omega" = rep(0, ncol(data_spatial$X_ij_omega)),
                            "eps_omega_st" = matrix(0,  nrow=mesh$n,ncol=max(data$Year) ),
                            "eps_linf" = rep(0,  length(unique(data$Lake))),
                            "eps_t0" = rep(0, length(unique(data$Lake))),
                            "ln_cv" = -2.5043055,
                            "ln_kappa" = -2.7,
                            "ln_tau_O" = 0,
                            "rho" = 0.9 )

  #Fit the nonspatial model:
  obj_nonspatial <- MakeADFun(data=data_nonspatial, parameters=parameters_nonspatial,
                              random=random_nonspatial, hessian=FALSE, DLL=VersionNonSpatial)

  opt_nonspatial = TMBhelper::Optimize(obj=obj_nonspatial,
                                       control=list(eval.max=1000, iter.max=1000),
                                       getsd=T, newtonsteps=1, bias.correct=F)

  rep_nonspatial <- obj_nonspatial$report()

  #fit the spatial model:
  obj_spatial <- MakeADFun(data=data_spatial, parameters=parameters_spatial,
                           random=random_spatial, hessian=FALSE, DLL=VersionSpatial)

  opt_spatial  = TMBhelper::Optimize(obj=obj_spatial,
                                     control=list(eval.max=1000, iter.max=1000),
                                     getsd=T, newtonsteps=1, bias.correct=F,
                                     lower=c(rep(-Inf,13),-0.999), upper=c(rep(Inf,13),0.999))

  rep_spatial <- obj_spatial$report()

  #Record the results
  nonspatial[which(to==k)] <- rep_nonspatial$pred_jnll / sum(Partition_i)
  spatial[which(to==k)] <- rep_spatial$pred_jnll / sum(Partition_i)

  print(which(to==k))
}
#--------------------------------------------------------------------------------------------------

#Leave one lake out cross validation:
loo_nonspatial <- loo_spatial <- rep(NA, length(unique(data$Lake)))

for(k in unique(data$Lake)){
  #which lake to predict:
  Partition_i = ifelse(data$Lake==k,1,0)

  data_nonspatial = list("Nobs"=nrow(data), "length_i"=data$TL, "age_i" = data$Age,
                         "lake_i" = data$Lake - 1, "X_ij_omega"= model.matrix(~ -1 + data$wallEffDen.Std +
                                                                                data$compEffDen.Std +
                                                                                data$GDD.Std  +
                                                                                data$wallEffDen.Std:data$compEffDen.Std),
                         "sex_i" = data$SexCode,"Nlakes" = length(unique(data$Lake)),"CTL" = CTL,
                         "predTF_i"=Partition_i)

  data_spatial = list("Nobs"      = nrow(data), "length_i"=data$TL, "age_i" = data$Age,
                      "lake_i"    = data$Lake - 1, "sex_i" = data$SexCode,
                      "X_ij_omega"= model.matrix(~ -1 + data$wallEffDen.Std + #.Std --> already standardized
                                                   data$compEffDen.Std +
                                                   data$GDD.Std  +
                                                   data$wallEffDen.Std:data$compEffDen.Std),
                      "Nlakes" = length(unique(data$Lake)),
                      "spdeMatrices" = spdeMatrices,
                      "s_i" = data$Lake-1,
                      "t_i" = data$Year-1,
                      "CTL" = CTL,
                      "predTF_i"=Partition_i)

  parameters_nonspatial = list("ln_global_omega" = log(14),
                               "ln_global_linf" = log(45), "ln_sd_linf" = log(7), "global_tzero" = -1,
                               "ln_sd_tzero" = log(3), "ln_b_sex" = log(4.760871), "b_j_omega" = rep(0, ncol(data_nonspatial$X_ij_omega)),
                               "eps_omega" = rep(0, data_nonspatial$Nlakes ), "eps_linf" = rep(0, data_nonspatial$Nlakes ),
                               "eps_t0" = rep(0, data_nonspatial$Nlakes ),"ln_cv" = log(0.2), "ln_sd_omega" = log(4.5) )

  parameters_spatial = list("ln_global_omega" = log(13.28954),
                            "ln_global_linf" = log(55),
                            "ln_sd_linf" = log(7.275452),
                            "global_tzero" = -1.260783461,
                            "ln_sd_tzero" = log(0.538755),
                            "ln_b_sex" = log(4.760871),
                            "b_j_omega" = rep(0, ncol(data_spatial$X_ij_omega)),
                            "eps_omega_st" = matrix(0,  nrow=mesh$n,ncol=max(data$Year) ),
                            "eps_linf" = rep(0,  length(unique(data$Lake))),
                            "eps_t0" = rep(0, length(unique(data$Lake))),
                            "ln_cv" = -2.5043055,
                            "ln_kappa" = -2.7,
                            "ln_tau_O" = 0,
                            "rho" = 0.9 )

  #Fit the nonspatial model:
  obj_nonspatial <- MakeADFun(data=data_nonspatial, parameters=parameters_nonspatial,
                              random=random_nonspatial, hessian=FALSE, DLL=VersionNonSpatial)

  opt_nonspatial = TMBhelper::Optimize(obj=obj_nonspatial,
                                       control=list(eval.max=1000, iter.max=1000),
                                       getsd=T, newtonsteps=1, bias.correct=F)

  rep_nonspatial <- obj_nonspatial$report()

  #fit the spatial model:
  obj_spatial <- MakeADFun(data=data_spatial, parameters=parameters_spatial,
                           random=random_spatial, hessian=FALSE, DLL=VersionSpatial)

  opt_spatial  = TMBhelper::Optimize(obj=obj_spatial,
                                     control=list(eval.max=1000, iter.max=1000),
                                     getsd=T, newtonsteps=1, bias.correct=F,
                                     lower=c(rep(-Inf,13),-0.999), upper=c(rep(Inf,13),0.999))

  rep_spatial <- obj_spatial$report()

  #Record the results
  loo_nonspatial[k] <- rep_nonspatial$pred_jnll / sum(Partition_i)
  loo_spatial[k] <- rep_spatial$pred_jnll / sum(Partition_i)
  print(k)
}

cv_nonspatial = mean(nonspatial)
cv_spatial = mean(spatial)

cv_loo_nonspatial = mean(loo_nonspatial)
cv_loo_spatial = mean(loo_spatial)
