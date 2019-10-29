#-------------------------------------------------------------------------------------------------------------------
#Nonspatial and Spatial-temporal von Bertalanffy with density dependent growth regression
#Galluci and Quinn (1979) Parameterization
#
#Model structure:
#length_pred_i = linfinity*(1-exp(-( omega/linfinity) * (age_i (i) - t0 )))
#
# Where:
# linfinity = linfinity_global + sex_effect_i + eps_linf
# t0  = t0_global + eps_t0
# omega   = omega_global + eta_fixed_i + eps_omega_st (spatial) OR eps_omega (nonspatial)
#
#eta_fixed_i = Effective Intraspecific (Walleye) Density +
#              Effective Interspecific Density +
#              Interaction b/n Intra & Inter Density +
#              Growing Degree Days
#
# Probability of Random Coefficients:
# epsilon_linf ~ N(0, sd_linf)
# epsilon_t0 ~ N(0, sd_tzero)
# epsilon_omega ~ N(0, sd_omega) --> nonspatial model
#
# eps_omega_st = rho*eps_omega_s,t-1 + u_st --> spatial model
# where u_st ~ N(0, SIMGA) --> SIGMA is estimated as per INLA approach
#
# Likelihoods:
# if(CTL == 1 == Normal)
# if(CTL == 2 == Lognormal)
# if(CTL == 3 == Gamma)
#
#Cahill 26 Oct 2019
#-------------------------------------------------------------------------------------------------------------------
start_time <- Sys.time()

library(dplyr)
library(plot3D)
library(plotly)
library(viridis)
library(magick)
library(lattice)
library(ggmap)
library(TMB)
library(INLA)
library(ggpubr)

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#Read in the data
#-------------------------------------------------------------------------------------------------------------------

data <- readRDS("C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/analysis/vB_analysis_august_2019_cahill.rds")
data <- as.data.frame(data)

sum(data$FL < 4.0)
data <- data[which(data$FL > 4.0),]
data <- data[which(data$FL  < 80), ]
nrow(data)
#-------------------------------------------------------------------------------------------------------------------

AICs <- data.frame(Model=c("Normal_Nonspatial", "Normal_Nonspatial_Reduced",
                           "Lognormal_Nonspatial", "Lognormal_Nonspatial_Reduced",
                           "Gamma_Nonspatial", "Gamma_Nonspatial_Reduced",
                           "Normal_Spatial", "Normal_Spatial_Reduced",
                           "Lognormal_Spatial", "Lognormal_Spatial_Reduced",
                           "Gamma_Spatial", "Gamma_Spatial_Reduced"),
                   AIC_REML=NA, AIC_ML=NA, h_block_cv=NA, lolo_cv=NA)

#-------------------------------------------------------------------------------------------------------------------
#Fit the nonspatial models:

setwd("C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/executables")
VersionNonSpatial = "vb_nonspatial"
compile(paste0(VersionNonSpatial, ".cpp"))
dyn.load( dynlib(VersionNonSpatial) )

#Normal Likelihood                                                                                                                                          "ln_b_sex", "b_j_omega") )
CTL = 1
random_nonspatial = c("eps_linf", "eps_t0", "eps_omega")

data_nonspatial = list("Nobs"=nrow(data), "length_i"=data$TL, "age_i" = data$Age,
                       "lake_i" = data$Lake - 1, "X_ij_omega"= model.matrix(~ -1 + data$wallEffDen.Std +
                                                                              data$compEffDen.Std +
                                                                              data$GDD.Std  +
                                                                              data$wallEffDen.Std:data$compEffDen.Std),
                       "sex_i" = data$SexCode,"Nlakes" = length(unique(data$Lake)),"CTL" = CTL,
                       "predTF_i"=rep(0, nrow(data)))


parameters_nonspatial = list("ln_global_omega" = log(14),  "ln_sd_omega" = log(4.5),
                             "ln_global_linf" = log(45), "ln_sd_linf" = log(7), "global_tzero" = -1,
                             "ln_sd_tzero" = log(3), "ln_b_sex" = log(4.760871), "b_j_omega" = rep(0, ncol(data_nonspatial$X_ij_omega)),
                             "eps_omega" = rep(0, data_nonspatial$Nlakes ),
                             "eps_linf" = rep(0, data_nonspatial$Nlakes ),
                             "eps_t0" = rep(0, data_nonspatial$Nlakes ),"ln_cv" = log(0.2)
)

obj_nonspatial <- MakeADFun(data=data_nonspatial, parameters=parameters_nonspatial,
                            random=random_nonspatial, hessian=FALSE, DLL=VersionNonSpatial)

opt_nonspatial = TMBhelper::Optimize(obj=obj_nonspatial,
                                     control=list(eval.max=1000, iter.max=1000),
                                     getsd=T, newtonsteps=1, bias.correct=F)

SD = sdreport( obj_nonspatial )
final_gradient = obj_nonspatial$gr( opt_nonspatial$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

AICs$AIC_ML[which(AICs$Model=="Normal_Nonspatial")] <- opt_nonspatial$AIC

ParHat_nonspatial = as.list( opt_nonspatial$SD, "Estimate" )
SEHat_nonspatial  = as.list( opt_nonspatial$SD, "Std. Error" )

#Re-run with REML
Use_REML=T
if( Use_REML==TRUE ) random_nonspatial = union( random_nonspatial, c("ln_global_omega",
                                                                     "ln_global_linf", "global_tzero",
                                                                     "ln_b_sex", "b_j_omega") )

obj_nonspatial <- MakeADFun(data=data_nonspatial, parameters=parameters_nonspatial,
                            random=random_nonspatial, hessian=FALSE, DLL=VersionNonSpatial)

opt_nonspatial = TMBhelper::Optimize(obj=obj_nonspatial,
                                     control=list(eval.max=1000, iter.max=1000),
                                     getsd=T, newtonsteps=1, bias.correct=F)

SD = sdreport( obj_nonspatial )
final_gradient = obj_nonspatial$gr( opt_nonspatial$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

AICs$AIC_REML[which(AICs$Model=="Normal_Nonspatial")] <- opt_nonspatial$AIC

#Nonspatial normal reduced
CTL = 1
random_nonspatial = c("eps_linf", "eps_t0", "eps_omega")

data_nonspatial = list("Nobs"=nrow(data), "length_i"=data$TL, "age_i" = data$Age,
                       "lake_i" = data$Lake - 1, "X_ij_omega"= model.matrix(~ -1 + data$wallEffDen.Std +
                                                                              data$compEffDen.Std +
                                                                              data$GDD.Std),
                       "sex_i" = data$SexCode,"Nlakes" = length(unique(data$Lake)),"CTL" = CTL,
                       "predTF_i"=rep(0, nrow(data)))

parameters_nonspatial = list("ln_global_omega" = log(14),  "ln_sd_omega" = log(4.5),
                             "ln_global_linf" = log(45), "ln_sd_linf" = log(7), "global_tzero" = -1,
                             "ln_sd_tzero" = log(3), "ln_b_sex" = log(4.760871), "b_j_omega" = rep(0, ncol(data_nonspatial$X_ij_omega)),
                             "eps_omega" = rep(0, data_nonspatial$Nlakes ),
                             "eps_linf" = rep(0, data_nonspatial$Nlakes ),
                             "eps_t0" = rep(0, data_nonspatial$Nlakes ),"ln_cv" = log(0.2)
)

obj_nonspatial <- MakeADFun(data=data_nonspatial, parameters=parameters_nonspatial,
                            random=random_nonspatial, hessian=FALSE, DLL=VersionNonSpatial)

opt_nonspatial = TMBhelper::Optimize(obj=obj_nonspatial,
                                     control=list(eval.max=1000, iter.max=1000),
                                     getsd=T, newtonsteps=1, bias.correct=F)

SD = sdreport( obj_nonspatial )
final_gradient = obj_nonspatial$gr( opt_nonspatial$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

AICs$AIC_ML[which(AICs$Model=="Normal_Nonspatial_Reduced")] <- opt_nonspatial$AIC

ParHat_nonspatial = as.list( opt_nonspatial$SD, "Estimate" )
SEHat_nonspatial  = as.list( opt_nonspatial$SD, "Std. Error" )

#Re-run with REML
Use_REML=T
if( Use_REML==TRUE ) random_nonspatial = union( random_nonspatial, c("ln_global_omega",
                                                                     "ln_global_linf", "global_tzero",
                                                                     "ln_b_sex", "b_j_omega") )

obj_nonspatial <- MakeADFun(data=data_nonspatial, parameters=parameters_nonspatial,
                            random=random_nonspatial, hessian=FALSE, DLL=VersionNonSpatial)

opt_nonspatial = TMBhelper::Optimize(obj=obj_nonspatial,
                                     control=list(eval.max=1000, iter.max=1000),
                                     getsd=T, newtonsteps=1, bias.correct=F)

SD = sdreport( obj_nonspatial )
final_gradient = obj_nonspatial$gr( opt_nonspatial$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

AICs$AIC_REML[which(AICs$Model=="Normal_Nonspatial_Reduced")] <- opt_nonspatial$AIC

#####################################
#Nonspatial Lognormal Models
#####################################
CTL = 2

data_nonspatial = list("Nobs"=nrow(data), "length_i"=data$TL, "age_i" = data$Age,
                       "lake_i" = data$Lake - 1, "X_ij_omega"= model.matrix(~ -1 + data$wallEffDen.Std +
                                                                              data$compEffDen.Std +
                                                                              data$GDD.Std  +
                                                                              data$wallEffDen.Std:data$compEffDen.Std),
                       "sex_i" = data$SexCode,"Nlakes" = length(unique(data$Lake)),"CTL" = CTL,
                       "predTF_i"=rep(0, nrow(data)))


parameters_nonspatial = list("ln_global_omega" = log(14),  "ln_sd_omega" = log(4.5),
                             "ln_global_linf" = log(45), "ln_sd_linf" = log(7), "global_tzero" = -1,
                             "ln_sd_tzero" = log(3), "ln_b_sex" = log(4.760871), "b_j_omega" = rep(0, ncol(data_nonspatial$X_ij_omega)),
                             "eps_omega" = rep(0, data_nonspatial$Nlakes ),
                             "eps_linf" = rep(0, data_nonspatial$Nlakes ),
                             "eps_t0" = rep(0, data_nonspatial$Nlakes ),"ln_cv" = log(0.2)
)

obj_nonspatial <- MakeADFun(data=data_nonspatial, parameters=parameters_nonspatial,
                            random=random_nonspatial, hessian=FALSE, DLL=VersionNonSpatial)

opt_nonspatial = TMBhelper::Optimize(obj=obj_nonspatial,
                                     control=list(eval.max=1000, iter.max=1000),
                                     getsd=T, newtonsteps=1, bias.correct=F)

SD = sdreport( obj_nonspatial )
final_gradient = obj_nonspatial$gr( opt_nonspatial$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

AICs$AIC_ML[which(AICs$Model=="Lognormal_Nonspatial")] <- opt_nonspatial$AIC

#Re-run with REML
Use_REML=T
if( Use_REML==TRUE ) random_nonspatial = union( random_nonspatial, c("ln_global_omega",
                                                                     "ln_global_linf", "global_tzero",
                                                                     "ln_b_sex", "b_j_omega") )

obj_nonspatial <- MakeADFun(data=data_nonspatial, parameters=parameters_nonspatial,
                            random=random_nonspatial, hessian=FALSE, DLL=VersionNonSpatial)

opt_nonspatial = TMBhelper::Optimize(obj=obj_nonspatial,
                                     control=list(eval.max=1000, iter.max=1000),
                                     getsd=T, newtonsteps=1, bias.correct=F)

SD = sdreport( obj_nonspatial )
final_gradient = obj_nonspatial$gr( opt_nonspatial$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

AICs$AIC_REML[which(AICs$Model=="Lognormal_Nonspatial")] <- opt_nonspatial$AIC

#####################################
#Nonspatial Lognormal Models Reduced
#####################################

CTL = 2

data_nonspatial = list("Nobs"=nrow(data), "length_i"=data$TL, "age_i" = data$Age,
                       "lake_i" = data$Lake - 1, "X_ij_omega"= model.matrix(~ -1 + data$wallEffDen.Std +
                                                                              data$compEffDen.Std +
                                                                              data$GDD.Std ),
                       "sex_i" = data$SexCode,"Nlakes" = length(unique(data$Lake)),"CTL" = CTL,
                       "predTF_i"=rep(0, nrow(data)))


parameters_nonspatial = list("ln_global_omega" = log(14),  "ln_sd_omega" = log(4.5),
                             "ln_global_linf" = log(45), "ln_sd_linf" = log(7), "global_tzero" = -1,
                             "ln_sd_tzero" = log(3), "ln_b_sex" = log(4.760871), "b_j_omega" = rep(0, ncol(data_nonspatial$X_ij_omega)),
                             "eps_omega" = rep(0, data_nonspatial$Nlakes ),
                             "eps_linf" = rep(0, data_nonspatial$Nlakes ),
                             "eps_t0" = rep(0, data_nonspatial$Nlakes ),"ln_cv" = log(0.2)
)

obj_nonspatial <- MakeADFun(data=data_nonspatial, parameters=parameters_nonspatial,
                            random=random_nonspatial, hessian=FALSE, DLL=VersionNonSpatial)

opt_nonspatial = TMBhelper::Optimize(obj=obj_nonspatial,
                                     control=list(eval.max=1000, iter.max=1000),
                                     getsd=T, newtonsteps=1, bias.correct=F)

SD = sdreport( obj_nonspatial )
final_gradient = obj_nonspatial$gr( opt_nonspatial$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

AICs$AIC_ML[which(AICs$Model=="Lognormal_Nonspatial_Reduced")] <- opt_nonspatial$AIC

#Re-run with REML
Use_REML=T
if( Use_REML==TRUE ) random_nonspatial = union( random_nonspatial, c("ln_global_omega",
                                                                     "ln_global_linf", "global_tzero",
                                                                     "ln_b_sex", "b_j_omega") )

obj_nonspatial <- MakeADFun(data=data_nonspatial, parameters=parameters_nonspatial,
                            random=random_nonspatial, hessian=FALSE, DLL=VersionNonSpatial)

opt_nonspatial = TMBhelper::Optimize(obj=obj_nonspatial,
                                     control=list(eval.max=1000, iter.max=1000),
                                     getsd=T, newtonsteps=1, bias.correct=F)

SD = sdreport( obj_nonspatial )
final_gradient = obj_nonspatial$gr( opt_nonspatial$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

AICs$AIC_REML[which(AICs$Model=="Lognormal_Nonspatial_Reduced")] <- opt_nonspatial$AIC

###########################
#Nonspatial Gamma Models
###########################

CTL = 3

data_nonspatial = list("Nobs"=nrow(data), "length_i"=data$TL, "age_i" = data$Age,
                       "lake_i" = data$Lake - 1, "X_ij_omega"= model.matrix(~ -1 + data$wallEffDen.Std +
                                                                              data$compEffDen.Std +
                                                                              data$GDD.Std  +
                                                                              data$wallEffDen.Std:data$compEffDen.Std),
                       "sex_i" = data$SexCode,"Nlakes" = length(unique(data$Lake)),"CTL" = CTL,
                       "predTF_i"=rep(0, nrow(data)))


parameters_nonspatial = list("ln_global_omega" = log(14),  "ln_sd_omega" = log(4.5),
                             "ln_global_linf" = log(45), "ln_sd_linf" = log(7), "global_tzero" = -1,
                             "ln_sd_tzero" = log(3), "ln_b_sex" = log(4.760871), "b_j_omega" = rep(0, ncol(data_nonspatial$X_ij_omega)),
                             "eps_omega" = rep(0, data_nonspatial$Nlakes ),
                             "eps_linf" = rep(0, data_nonspatial$Nlakes ),
                             "eps_t0" = rep(0, data_nonspatial$Nlakes ),"ln_cv" = log(0.2)
)

obj_nonspatial <- MakeADFun(data=data_nonspatial, parameters=parameters_nonspatial,
                            random=random_nonspatial, hessian=FALSE, DLL=VersionNonSpatial)

opt_nonspatial = TMBhelper::Optimize(obj=obj_nonspatial,
                                     control=list(eval.max=1000, iter.max=1000),
                                     getsd=T, newtonsteps=1, bias.correct=F)

SD = sdreport( obj_nonspatial )
final_gradient = obj_nonspatial$gr( opt_nonspatial$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

AICs$AIC_ML[which(AICs$Model=="Gamma_Nonspatial")] <- opt_nonspatial$AIC

#Re-run with REML
Use_REML=T
if( Use_REML==TRUE ) random_nonspatial = union( random_nonspatial, c("ln_global_omega",
                                                                     "ln_global_linf", "global_tzero",
                                                                     "ln_b_sex", "b_j_omega") )

obj_nonspatial <- MakeADFun(data=data_nonspatial, parameters=parameters_nonspatial,
                            random=random_nonspatial, hessian=FALSE, DLL=VersionNonSpatial)

opt_nonspatial = TMBhelper::Optimize(obj=obj_nonspatial,
                                     control=list(eval.max=1000, iter.max=1000),
                                     getsd=T, newtonsteps=1, bias.correct=F)

SD = sdreport( obj_nonspatial )
final_gradient = obj_nonspatial$gr( opt_nonspatial$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

AICs$AIC_REML[which(AICs$Model=="Gamma_Nonspatial")] <- opt_nonspatial$AIC

###########################
#Nonspatial Gamma Models
###########################

CTL = 3

data_nonspatial = list("Nobs"=nrow(data), "length_i"=data$TL, "age_i" = data$Age,
                       "lake_i" = data$Lake - 1, "X_ij_omega"= model.matrix(~ -1 + data$wallEffDen.Std +
                                                                              data$compEffDen.Std +
                                                                              data$GDD.Std),
                       "sex_i" = data$SexCode,"Nlakes" = length(unique(data$Lake)),"CTL" = CTL,
                       "predTF_i"=rep(0, nrow(data)))


parameters_nonspatial = list("ln_global_omega" = log(14),  "ln_sd_omega" = log(4.5),
                             "ln_global_linf" = log(45), "ln_sd_linf" = log(7), "global_tzero" = -1,
                             "ln_sd_tzero" = log(3), "ln_b_sex" = log(4.760871), "b_j_omega" = rep(0, ncol(data_nonspatial$X_ij_omega)),
                             "eps_omega" = rep(0, data_nonspatial$Nlakes ),
                             "eps_linf" = rep(0, data_nonspatial$Nlakes ),
                             "eps_t0" = rep(0, data_nonspatial$Nlakes ),"ln_cv" = log(0.2)
)

obj_nonspatial <- MakeADFun(data=data_nonspatial, parameters=parameters_nonspatial,
                            random=random_nonspatial, hessian=FALSE, DLL=VersionNonSpatial)

opt_nonspatial = TMBhelper::Optimize(obj=obj_nonspatial,
                                     control=list(eval.max=1000, iter.max=1000),
                                     getsd=T, newtonsteps=1, bias.correct=F)

SD = sdreport( obj_nonspatial )
final_gradient = obj_nonspatial$gr( opt_nonspatial$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

AICs$AIC_ML[which(AICs$Model=="Gamma_Nonspatial_Reduced")] <- opt_nonspatial$AIC

#Re-run with REML
Use_REML=T
if( Use_REML==TRUE ) random_nonspatial = union( random_nonspatial, c("ln_global_omega",
                                                                     "ln_global_linf", "global_tzero",
                                                                     "ln_b_sex", "b_j_omega") )

obj_nonspatial <- MakeADFun(data=data_nonspatial, parameters=parameters_nonspatial,
                            random=random_nonspatial, hessian=FALSE, DLL=VersionNonSpatial)

opt_nonspatial = TMBhelper::Optimize(obj=obj_nonspatial,
                                     control=list(eval.max=1000, iter.max=1000),
                                     getsd=T, newtonsteps=1, bias.correct=F)

SD = sdreport( obj_nonspatial )
final_gradient = obj_nonspatial$gr( opt_nonspatial$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

AICs$AIC_REML[which(AICs$Model=="Gamma_Nonspatial_Reduced")] <- opt_nonspatial$AIC

#-------------------------------------------------------------------------------------------------------------------
#Fit the Spatial models:
#-------------------------------------------------------------------------------------------------------------------

#Build a mesh for the spatial models:
loc_xy <- unique(data[ ,c("X_TTM_c","Y_TTM_c") ] )
loc_xy <- loc_xy/1000 #Put distance in kms

#mesh = inla.mesh.2d(loc=loc_xy, max.edge=c(62,1000)) #Better mesh, but slower
mesh = inla.mesh.create( loc_xy, refine=TRUE, extend=-0.5, cutoff=0.01 ) #faster mesh

mesh$n

#png(file="Mesh.png",width=9.50,height=7.00,units="in",res=600)
plot(mesh)
points(loc_xy, col="Steelblue", pch=1)
#dev.off()

# Other meshes to try--some run faster, some slower, but all give same answers (Cahill pers. obs)
# RangeGuess <- 30    #~ 1/3 study area
# MaxEdge    <- RangeGuess/ 5 #as per https://haakonbakka.bitbucket.io/btopic104.html
#
# bound.outer <- diff(range(loc_xy[,1]))/3
# RangeGuess <- 50
#
# mesh1 = inla.mesh.create( loc_xy, refine=TRUE, extend=-0.5, cutoff=0.01 ) #fastest mesh
#
# mesh2 = inla.mesh.2d(loc=loc_xy, max.edge = c(2,5)*RangeGuess )
#
# mesh3 = inla.mesh.2d(loc=loc_xy, max.edge = c(3,5)*RangeGuess )
#
# mesh4 = inla.mesh.2d(loc=loc_xy, max.edge = c(1,5)*RangeGuess )
#
# mesh5 = inla.mesh.2d(loc=loc_xy, max.edge = c(1,5)*RangeGuess,
#                      cutoff = MaxEdge,
#                      offset = c(MaxEdge, bound.outer))

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

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
                    "X_ij_omega"= model.matrix(~ -1 + data$wallEffDen.Std +
                                                 data$compEffDen.Std + data$GDD.Std  +
                                                 data$wallEffDen.Std:data$compEffDen.Std),
                    "Nlakes" = length(unique(data$Lake)), "spdeMatrices" = spdeMatrices,
                    "s_i" = data$Lake-1, "t_i" = data$Year-1, "CTL" = CTL, "predTF_i"=rep(0, nrow(data)))


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
                          "logit_rho" = 3.2 ) #2 * plogis(3) - 1 = 0.9

random_spatial = c("eps_omega_st", "eps_linf", "eps_t0")

Use_REML = F
if( Use_REML==TRUE ) random_spatial = union( random_spatial, c("ln_global_omega",
                                                               "ln_global_linf", "global_tzero",
                                                               "ln_b_sex", "b_j_omega") )

obj_spatial = MakeADFun(data=data_spatial, parameters=parameters_spatial,
                        random=random_spatial, hessian=FALSE, DLL=VersionSpatial)

opt_spatial  = TMBhelper::Optimize(obj=obj_spatial,
                                   control=list(eval.max=1000, iter.max=1000),
                                   getsd=T, newtonsteps=1, bias.correct=T)

SD = sdreport( obj_spatial )
final_gradient = obj_spatial$gr( opt_spatial$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")


AICs$AIC_ML[which(AICs$Model=="Normal_Spatial")] <- opt_spatial$AIC

#Re-run with REML
Use_REML=T
if( Use_REML==TRUE ) random_spatial = union( random_spatial, c("ln_global_omega",
                                                               "ln_global_linf", "global_tzero",
                                                               "ln_b_sex", "b_j_omega") )

obj_spatial = MakeADFun(data=data_spatial, parameters=parameters_spatial,
                        random=random_spatial, hessian=FALSE, DLL=VersionSpatial)

opt_spatial  = TMBhelper::Optimize(obj=obj_spatial,
                                   control=list(eval.max=1000, iter.max=1000),
                                   getsd=T, newtonsteps=1, bias.correct=T)

SD = sdreport( obj_spatial )
final_gradient = obj_spatial$gr( opt_spatial$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

AICs$AIC_REML[which(AICs$Model=="Normal_Spatial")] <- opt_spatial$AIC

############################################
#Normal Likelihood--spatial reduced
############################################
CTL = 1

data_spatial = list("Nobs" = nrow(data), "length_i"=data$TL, "age_i" = data$Age,
                    "lake_i" = data$Lake - 1, "sex_i" = data$SexCode,
                    "X_ij_omega"= model.matrix(~ -1 + data$wallEffDen.Std +
                                                 data$compEffDen.Std + data$GDD.Std ),
                    "Nlakes" = length(unique(data$Lake)), "spdeMatrices" = spdeMatrices,
                    "s_i" = data$Lake-1, "t_i" = data$Year-1, "CTL" = CTL, "predTF_i"=rep(0, nrow(data)))


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
                          "logit_rho" = 3.2 ) #2 * plogis(3) - 1 = 0.9

random_spatial = c("eps_omega_st", "eps_linf", "eps_t0")

Use_REML = F
if( Use_REML==TRUE ) random_spatial = union( random_spatial, c("ln_global_omega",
                                                               "ln_global_linf", "global_tzero",
                                                               "ln_b_sex", "b_j_omega") )

obj_spatial = MakeADFun(data=data_spatial, parameters=parameters_spatial,
                        random=random_spatial, hessian=FALSE, DLL=VersionSpatial)

opt_spatial  = TMBhelper::Optimize(obj=obj_spatial,
                                   control=list(eval.max=1000, iter.max=1000),
                                   getsd=T, newtonsteps=1, bias.correct=T)

SD = sdreport( obj_spatial )
final_gradient = obj_spatial$gr( opt_spatial$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")


AICs$AIC_ML[which(AICs$Model=="Normal_Spatial_Reduced")] <- opt_spatial$AIC

#Re-run with REML
Use_REML=T
if( Use_REML==TRUE ) random_spatial = union( random_spatial, c("ln_global_omega",
                                                               "ln_global_linf", "global_tzero",
                                                               "ln_b_sex", "b_j_omega") )

obj_spatial = MakeADFun(data=data_spatial, parameters=parameters_spatial,
                        random=random_spatial, hessian=FALSE, DLL=VersionSpatial)

opt_spatial  = TMBhelper::Optimize(obj=obj_spatial,
                                   control=list(eval.max=1000, iter.max=1000),
                                   getsd=T, newtonsteps=1, bias.correct=T)

SD = sdreport( obj_spatial )
final_gradient = obj_spatial$gr( opt_spatial$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

AICs$AIC_REML[which(AICs$Model=="Normal_Spatial_Reduced")] <- opt_spatial$AIC

#################################################3

#Spatial Lognormal Models

CTL = 2

data_spatial = list("Nobs" = nrow(data), "length_i"=data$TL, "age_i" = data$Age,
                    "lake_i" = data$Lake - 1, "sex_i" = data$SexCode,
                    "X_ij_omega"= model.matrix(~ -1 + data$wallEffDen.Std +
                                                 data$compEffDen.Std + data$GDD.Std  +
                                                 data$wallEffDen.Std:data$compEffDen.Std),
                    "Nlakes" = length(unique(data$Lake)), "spdeMatrices" = spdeMatrices,
                    "s_i" = data$Lake-1, "t_i" = data$Year-1, "CTL" = CTL, "predTF_i"=rep(0, nrow(data)))


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
                          "logit_rho" = 3.2 ) #2 * plogis(3) - 1 = 0.9

random_spatial = c("eps_omega_st", "eps_linf", "eps_t0")

Use_REML = F
if( Use_REML==TRUE ) random_spatial = union( random_spatial, c("ln_global_omega",
                                                               "ln_global_linf", "global_tzero",
                                                               "ln_b_sex", "b_j_omega") )

obj_spatial = MakeADFun(data=data_spatial, parameters=parameters_spatial,
                        random=random_spatial, hessian=FALSE, DLL=VersionSpatial)

opt_spatial  = TMBhelper::Optimize(obj=obj_spatial,
                                   control=list(eval.max=1000, iter.max=1000),
                                   getsd=T, newtonsteps=1, bias.correct=T)

SD = sdreport( obj_spatial )
final_gradient = obj_spatial$gr( opt_spatial$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

AICs$AIC_ML[which(AICs$Model=="Lognormal_Spatial")] <- opt_spatial$AIC

#Re-run with REML
Use_REML=T
if( Use_REML==TRUE ) random_spatial = union( random_spatial, c("ln_global_omega",
                                                               "ln_global_linf", "global_tzero",
                                                               "ln_b_sex", "b_j_omega") )

obj_spatial = MakeADFun(data=data_spatial, parameters=parameters_spatial,
                        random=random_spatial, hessian=FALSE, DLL=VersionSpatial)

opt_spatial  = TMBhelper::Optimize(obj=obj_spatial,
                                   control=list(eval.max=1000, iter.max=1000),
                                   getsd=T, newtonsteps=1, bias.correct=T)

SD = sdreport( obj_spatial )
final_gradient = obj_spatial$gr( opt_spatial$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

AICs$AIC_REML[which(AICs$Model=="Lognormal_Spatial")] <- opt_spatial$AIC

#################################################
#Spatial Lognormal Models Reduced
#################################################

CTL = 2

data_spatial = list("Nobs" = nrow(data), "length_i"=data$TL, "age_i" = data$Age,
                    "lake_i" = data$Lake - 1, "sex_i" = data$SexCode,
                    "X_ij_omega"= model.matrix(~ -1 + data$wallEffDen.Std +
                                                 data$compEffDen.Std + data$GDD.Std),
                    "Nlakes" = length(unique(data$Lake)), "spdeMatrices" = spdeMatrices,
                    "s_i" = data$Lake-1, "t_i" = data$Year-1, "CTL" = CTL, "predTF_i"=rep(0, nrow(data)))


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
                          "logit_rho" = 3.2 ) #2 * plogis(3) - 1 = 0.9

random_spatial = c("eps_omega_st", "eps_linf", "eps_t0")

Use_REML = F
if( Use_REML==TRUE ) random_spatial = union( random_spatial, c("ln_global_omega",
                                                               "ln_global_linf", "global_tzero",
                                                               "ln_b_sex", "b_j_omega") )

obj_spatial = MakeADFun(data=data_spatial, parameters=parameters_spatial,
                        random=random_spatial, hessian=FALSE, DLL=VersionSpatial)

opt_spatial  = TMBhelper::Optimize(obj=obj_spatial,
                                   control=list(eval.max=1000, iter.max=1000),
                                   getsd=T, newtonsteps=1, bias.correct=T)

SD = sdreport( obj_spatial )
final_gradient = obj_spatial$gr( opt_spatial$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

AICs$AIC_ML[which(AICs$Model=="Lognormal_Spatial_Reduced")] <- opt_spatial$AIC

#Re-run with REML
Use_REML=T
if( Use_REML==TRUE ) random_spatial = union( random_spatial, c("ln_global_omega",
                                                               "ln_global_linf", "global_tzero",
                                                               "ln_b_sex", "b_j_omega") )

obj_spatial = MakeADFun(data=data_spatial, parameters=parameters_spatial,
                        random=random_spatial, hessian=FALSE, DLL=VersionSpatial)

opt_spatial  = TMBhelper::Optimize(obj=obj_spatial,
                                   control=list(eval.max=1000, iter.max=1000),
                                   getsd=T, newtonsteps=1, bias.correct=T)

SD = sdreport( obj_spatial )
final_gradient = obj_spatial$gr( opt_spatial$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

AICs$AIC_REML[which(AICs$Model=="Lognormal_Spatial_Reduced")] <- opt_spatial$AIC

####################################
#Spatial Gamma Models
###################################

CTL = 3

data_spatial = list("Nobs" = nrow(data), "length_i"=data$TL, "age_i" = data$Age,
                    "lake_i" = data$Lake - 1, "sex_i" = data$SexCode,
                    "X_ij_omega"= model.matrix(~ -1 + data$wallEffDen.Std +
                                                 data$compEffDen.Std + data$GDD.Std  +
                                                 data$wallEffDen.Std:data$compEffDen.Std),
                    "Nlakes" = length(unique(data$Lake)), "spdeMatrices" = spdeMatrices,
                    "s_i" = data$Lake-1, "t_i" = data$Year-1, "CTL" = CTL, "predTF_i"=rep(0, nrow(data)))

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
                          "logit_rho" = 3.2 ) #2 * plogis(3) - 1 = 0.9

random_spatial = c("eps_omega_st", "eps_linf", "eps_t0")

Use_REML = F
if( Use_REML==TRUE ) random_spatial = union( random_spatial, c("ln_global_omega",
                                                               "ln_global_linf", "global_tzero",
                                                               "ln_b_sex", "b_j_omega") )

obj_spatial = MakeADFun(data=data_spatial, parameters=parameters_spatial,
                        random=random_spatial, hessian=FALSE, DLL=VersionSpatial)

opt_spatial  = TMBhelper::Optimize(obj=obj_spatial,
                                   control=list(eval.max=1000, iter.max=1000),
                                   getsd=T, newtonsteps=1, bias.correct=T)

SD = sdreport( obj_spatial )
final_gradient = obj_spatial$gr( opt_spatial$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

AICs$AIC_ML[which(AICs$Model=="Gamma_Spatial")] <- opt_spatial$AIC

#Re-run with REML
Use_REML=T
if( Use_REML==TRUE ) random_spatial = union( random_spatial, c("ln_global_omega",
                                                               "ln_global_linf", "global_tzero",
                                                               "ln_b_sex", "b_j_omega") )

obj_spatial = MakeADFun(data=data_spatial, parameters=parameters_spatial,
                        random=random_spatial, hessian=FALSE, DLL=VersionSpatial)

opt_spatial  = TMBhelper::Optimize(obj=obj_spatial,
                                   control=list(eval.max=1000, iter.max=1000),
                                   getsd=T, newtonsteps=1, bias.correct=T)

SD = sdreport( obj_spatial )
final_gradient = obj_spatial$gr( opt_spatial$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

AICs$AIC_REML[which(AICs$Model=="Gamma_Spatial")] <- opt_spatial$AIC
ParHat = as.list( opt_spatial$SD, "Estimate" )

####################################
#Spatial Gamma Models Reduced
###################################

CTL = 3

data_spatial = list("Nobs" = nrow(data), "length_i"=data$TL, "age_i" = data$Age,
                    "lake_i" = data$Lake - 1, "sex_i" = data$SexCode,
                    "X_ij_omega"= model.matrix(~ -1 + data$wallEffDen.Std +
                                                 data$compEffDen.Std + data$GDD.Std ),
                    "Nlakes" = length(unique(data$Lake)), "spdeMatrices" = spdeMatrices,
                    "s_i" = data$Lake-1, "t_i" = data$Year-1, "CTL" = CTL, "predTF_i"=rep(0, nrow(data)))

parameters_spatial = list("ln_global_omega" = 2.69,
                          "ln_global_linf" = 4.021588,
                          "ln_sd_linf" = 1.990446,
                          "global_tzero" = -1.17,
                          "ln_sd_tzero" = log(0.538755),
                          "ln_b_sex" = log(1.525),
                          "b_j_omega" = rep(0, ncol(data_spatial$X_ij_omega)),
                          "eps_omega_st" = matrix(0,  nrow=mesh$n,ncol=max(data$Year) ),
                          "eps_linf" = rep(0,  length(unique(data$Lake))),
                          "eps_t0" = rep(0, length(unique(data$Lake))),
                          "ln_cv" = -2.5043055,
                          "ln_kappa" = -2.79,
                          "ln_tau_O" = -0.05,
                          "logit_rho" = 3.12 ) #2 * plogis(3) - 1 = 0.9

random_spatial = c("eps_omega_st", "eps_linf", "eps_t0")

Use_REML = F
if( Use_REML==TRUE ) random_spatial = union( random_spatial, c("ln_global_omega",
                                                               "ln_global_linf", "global_tzero",
                                                               "ln_b_sex", "b_j_omega") )

obj_spatial = MakeADFun(data=data_spatial, parameters=parameters_spatial,
                        random=random_spatial, hessian=FALSE, DLL=VersionSpatial)

opt_spatial  = TMBhelper::Optimize(obj=obj_spatial,
                                   control=list(eval.max=1000, iter.max=1000),
                                   getsd=T, newtonsteps=1, bias.correct=T)

SD = sdreport( obj_spatial )
final_gradient = obj_spatial$gr( opt_spatial$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

AICs$AIC_ML[which(AICs$Model=="Gamma_Spatial_Reduced")] <- opt_spatial$AIC

#Re-run with REML
Use_REML=T
if( Use_REML==TRUE ) random_spatial = union( random_spatial, c("ln_global_omega",
                                                               "ln_global_linf", "global_tzero",
                                                               "ln_b_sex", "b_j_omega") )

obj_spatial = MakeADFun(data=data_spatial, parameters=parameters_spatial,
                        random=random_spatial, hessian=FALSE, DLL=VersionSpatial)

opt_spatial  = TMBhelper::Optimize(obj=obj_spatial,
                                   control=list(eval.max=1000, iter.max=1000),
                                   getsd=T, newtonsteps=1, bias.correct=T)

SD = sdreport( obj_spatial )
final_gradient = obj_spatial$gr( opt_spatial$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

AICs$AIC_REML[which(AICs$Model=="Gamma_Spatial_Reduced")] <- opt_spatial$AIC

#------------------------------------------------------------------------------------------------------------
AICs$Delta_AIC_REML <- AICs$AIC_REML - min(AICs$AIC_REML)
AICs$Delta_AIC_ML <- AICs$AIC_ML - min(AICs$AIC_ML)
print(AICs)
saveRDS(AICs, file="C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/analysis/AICs")
#All spatial-temporal models ~ 9000-1000 points better than analogous nonspatial-temporal models (AIC_REML column)
#Normal likelihood appears to be most parsimonious (comparing ML AICs among spatial column)
#aics <- readRDS(file="C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/analysis/AICs")
#------------------------------------------------------------------------------------------------------------
