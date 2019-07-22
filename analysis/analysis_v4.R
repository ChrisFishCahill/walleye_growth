#-------------------------------------------------------------------------------------------------------------------
#Spatial-temporal von Bertalanffy with density dependent growth regression
#Galluci and Quinn (1979) Parameterization
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
#Cahill ^(;,;)^ July 2019
#TODO
#Clean code
#Git code
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

data <- readRDS("C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/analysis/vB_analysis_march_2019_cahill.rds")
data <- as.data.frame(data)

# convert mean effective density per net to mean effDen / area
# data$wallEffDen1 <- data$wallEffDen/data$Area
# data$wallEffDen.Std <- arm::rescale(data$wallEffDen1)
#-------------------------------------------------------------------------------------------------------------------

AICs <- data.frame(Model=c("Normal_Nonspatial", "Lognormal_Nonspatial", "Gamma_Nonspatial",
                           "Normal_Spatial", "Lognormal_Spatial", "Gamma_Spatial"), AIC_REML=NA, AIC_ML=NA, h_block_cv=NA,
                           lolo_cv=NA)

#-------------------------------------------------------------------------------------------------------------------
#Fit the nonspatial models REML:

setwd("C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/executables")
VersionNonSpatial = "vb_nonspatial"
compile(paste0(VersionNonSpatial, ".cpp"))
dyn.load( dynlib(VersionNonSpatial) )

random_nonspatial = c("eps_linf", "eps_t0", "eps_omega")

#Normal Likelihood                                                                                                                                          "ln_b_sex", "b_j_omega") )
CTL = 1

data_nonspatial = list("Nobs"=nrow(data), "length_i"=data$TL, "age_i" = data$Age,
                       "lake_i" = data$Lake - 1, "X_ij_omega"= model.matrix(~ -1 + data$wallEffDen.Std +
                                                                              data$compEffDen.Std +
                                                                              data$GDD.Std  +
                                                                              data$wallEffDen.Std:data$compEffDen.Std),
                       "sex_i" = data$SexCode,"Nlakes" = length(unique(data$Lake)),"CTL" = CTL,
                       "predTF_i"=rep(0, nrow(data)))


parameters_nonspatial = list("ln_global_omega" = log(14),
                             "ln_global_linf" = log(45), "ln_sd_linf" = log(7), "global_tzero" = -1,
                             "ln_sd_tzero" = log(3), "ln_b_sex" = log(4.760871), "b_j_omega" = rep(0, ncol(data_nonspatial$X_ij_omega)),
                             "eps_omega" = rep(0, data_nonspatial$Nlakes ), "eps_linf" = rep(0, data_nonspatial$Nlakes ),
                             "eps_t0" = rep(0, data_nonspatial$Nlakes ),"ln_cv" = log(0.2), "ln_sd_omega" = log(4.5) )


obj_nonspatial <- MakeADFun(data=data_nonspatial, parameters=parameters_nonspatial,
                            random=random_nonspatial, hessian=FALSE, DLL=VersionNonSpatial)

opt_nonspatial = TMBhelper::Optimize(obj=obj_nonspatial,
                                     control=list(eval.max=1000, iter.max=1000),
                                     getsd=T, newtonsteps=1, bias.correct=F)

AICs$AIC_ML[which(AICs$Model=="Normal_Nonspatial")] <- opt_nonspatial$AIC

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

AICs$AIC_REML[which(AICs$Model=="Normal_Nonspatial")] <- opt_nonspatial$AIC

#Nonspatial Lognormal Models
CTL = 2

data_nonspatial = list("Nobs"=nrow(data), "length_i"=data$TL, "age_i" = data$Age,
                       "lake_i" = data$Lake - 1, "X_ij_omega"= model.matrix(~ -1 + data$wallEffDen.Std +
                                                                              data$compEffDen.Std +
                                                                              data$GDD.Std  +
                                                                              data$wallEffDen.Std:data$compEffDen.Std),
                       "sex_i" = data$SexCode,"Nlakes" = length(unique(data$Lake)),"CTL" = CTL,
                       "predTF_i"=rep(0, nrow(data)))


parameters_nonspatial = list("ln_global_omega" = log(14),
                             "ln_global_linf" = log(45), "ln_sd_linf" = log(7), "global_tzero" = -1,
                             "ln_sd_tzero" = log(3), "ln_b_sex" = log(4.760871), "b_j_omega" = rep(0, ncol(data_nonspatial$X_ij_omega)),
                             "eps_omega" = rep(0, data_nonspatial$Nlakes ), "eps_linf" = rep(0, data_nonspatial$Nlakes ),
                             "eps_t0" = rep(0, data_nonspatial$Nlakes ),"ln_cv" = log(0.2), "ln_sd_omega" = log(4.5) )


obj_nonspatial <- MakeADFun(data=data_nonspatial, parameters=parameters_nonspatial,
                            random=random_nonspatial, hessian=FALSE, DLL=VersionNonSpatial)

opt_nonspatial = TMBhelper::Optimize(obj=obj_nonspatial,
                                     control=list(eval.max=1000, iter.max=1000),
                                     getsd=T, newtonsteps=1, bias.correct=F)

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

AICs$AIC_REML[which(AICs$Model=="Lognormal_Nonspatial")] <- opt_nonspatial$AIC

#Nonspatial Gamma Models
CTL = 3

data_nonspatial = list("Nobs"=nrow(data), "length_i"=data$TL, "age_i" = data$Age,
                       "lake_i" = data$Lake - 1, "X_ij_omega"= model.matrix(~ -1 + data$wallEffDen.Std +
                                                                              data$compEffDen.Std +
                                                                              data$GDD.Std  +
                                                                              data$wallEffDen.Std:data$compEffDen.Std),
                       "sex_i" = data$SexCode,"Nlakes" = length(unique(data$Lake)),"CTL" = CTL,
                       "predTF_i"=rep(0, nrow(data)))


parameters_nonspatial = list("ln_global_omega" = log(14),
                             "ln_global_linf" = log(45), "ln_sd_linf" = log(7), "global_tzero" = -1,
                             "ln_sd_tzero" = log(3), "ln_b_sex" = log(4.760871), "b_j_omega" = rep(0, ncol(data_nonspatial$X_ij_omega)),
                             "eps_omega" = rep(0, data_nonspatial$Nlakes ), "eps_linf" = rep(0, data_nonspatial$Nlakes ),
                             "eps_t0" = rep(0, data_nonspatial$Nlakes ),"ln_cv" = log(0.2), "ln_sd_omega" = log(4.5) )


obj_nonspatial <- MakeADFun(data=data_nonspatial, parameters=parameters_nonspatial,
                            random=random_nonspatial, hessian=FALSE, DLL=VersionNonSpatial)

opt_nonspatial = TMBhelper::Optimize(obj=obj_nonspatial,
                                     control=list(eval.max=1000, iter.max=1000),
                                     getsd=T, newtonsteps=1, bias.correct=F)

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

AICs$AIC_REML[which(AICs$Model=="Gamma_Nonspatial")] <- opt_nonspatial$AIC

#-------------------------------------------------------------------------------------------------------------------
#Fit the Spatial models:
#-------------------------------------------------------------------------------------------------------------------

#Build a mesh for the spatial models:
loc_xy <- unique(data[ ,c("X_TTM_c","Y_TTM_c") ] )
loc_xy <- loc_xy/1000 #Put distance in kms

mesh = inla.mesh.create( loc_xy, refine=TRUE, extend=-0.5, cutoff=0.01 )

#png(file="Mesh3.png",width=9.50,height=7.00,units="in",res=600)
plot(mesh)
points(loc_xy, col="Steelblue", pch=1)
#dev.off()

# Other meshes to try--much longer run time but same answers (Cahill pers. obs)

# RangeGuess <- 30    #~ 1/3 study area
# MaxEdge    <- RangeGuess/ 5 #as per https://haakonbakka.bitbucket.io/btopic104.html
# # #
# bound.outer <- diff(range(loc_xy[,1]))/3
# RangeGuess <- 50
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

#mesh=mesh4
mesh3 <- inla.mesh.2d(loc=loc_xy, max.edge=c(75,1000))
mesh3 <- inla.mesh.2d(loc=loc_xy, max.edge=c(62,1000))

plot(mesh3)
points(loc_xy, col="steelblue")
mesh3$n
#-------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
mesh = mesh3

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
                    "X_ij_omega"= model.matrix(~ -1 + data$wallEffDen.Std + #.Std --> already standardized
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
                          "rho" = 0.9 )

random_spatial = c("eps_omega_st", "eps_linf", "eps_t0")

Use_REML = F
if( Use_REML==TRUE ) random_spatial = union( random_spatial, c("ln_global_omega",
                                                               "ln_global_linf", "global_tzero",
                                                               "ln_b_sex", "b_j_omega") )

obj_spatial = MakeADFun(data=data_spatial, parameters=parameters_spatial,
                        random=random_spatial, hessian=FALSE, DLL=VersionSpatial)

opt_spatial  = TMBhelper::Optimize(obj=obj_spatial,
                                   control=list(eval.max=1000, iter.max=1000),
                                   getsd=T, newtonsteps=1, bias.correct=T,
                                   lower=c(rep(-Inf,13),-0.999), upper=c(rep(Inf,13),0.999))

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
                                   getsd=T, newtonsteps=1, bias.correct=T,
                                   lower=c(rep(-Inf,5),-0.999), upper=c(rep(Inf,5),0.999))

AICs$AIC_REML[which(AICs$Model=="Normal_Spatial")] <- opt_spatial$AIC

#Spatial Lognormal Models

CTL = 2

data_spatial = list("Nobs" = nrow(data), "length_i"=data$TL, "age_i" = data$Age,
                    "lake_i" = data$Lake - 1, "sex_i" = data$SexCode,
                    "X_ij_omega"= model.matrix(~ -1 + data$wallEffDen.Std + #.Std --> already standardized
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
                          "rho" = 0.9 )

random_spatial = c("eps_omega_st", "eps_linf", "eps_t0")

Use_REML = F
if( Use_REML==TRUE ) random_spatial = union( random_spatial, c("ln_global_omega",
                                                               "ln_global_linf", "global_tzero",
                                                               "ln_b_sex", "b_j_omega") )

obj_spatial = MakeADFun(data=data_spatial, parameters=parameters_spatial,
                        random=random_spatial, hessian=FALSE, DLL=VersionSpatial)

opt_spatial  = TMBhelper::Optimize(obj=obj_spatial,
                                   control=list(eval.max=1000, iter.max=1000),
                                   getsd=T, newtonsteps=1, bias.correct=T,
                                   lower=c(rep(-Inf,13),-0.999), upper=c(rep(Inf,13),0.999))

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
                                   getsd=T, newtonsteps=1, bias.correct=T,
                                   lower=c(rep(-Inf,5),-0.999), upper=c(rep(Inf,5),0.999))

AICs$AIC_REML[which(AICs$Model=="Lognormal_Spatial")] <- opt_spatial$AIC

#Spatial Gamma Models

CTL = 3

data_spatial = list("Nobs" = nrow(data), "length_i"=data$TL, "age_i" = data$Age,
                    "lake_i" = data$Lake - 1, "sex_i" = data$SexCode,
                    "X_ij_omega"= model.matrix(~ -1 + data$wallEffDen.Std + #.Std --> already standardized
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
                          "rho" = 0.9 )

random_spatial = c("eps_omega_st", "eps_linf", "eps_t0")

Use_REML = F
if( Use_REML==TRUE ) random_spatial = union( random_spatial, c("ln_global_omega",
                                                               "ln_global_linf", "global_tzero",
                                                               "ln_b_sex", "b_j_omega") )

obj_spatial = MakeADFun(data=data_spatial, parameters=parameters_spatial,
                        random=random_spatial, hessian=FALSE, DLL=VersionSpatial)

opt_spatial  = TMBhelper::Optimize(obj=obj_spatial,
                                   control=list(eval.max=1000, iter.max=1000),
                                   getsd=T, newtonsteps=1, bias.correct=T,
                                   lower=c(rep(-Inf,13),-0.999), upper=c(rep(Inf,13),0.999))

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
                                   getsd=T, newtonsteps=1, bias.correct=T,
                                   lower=c(rep(-Inf,5),-0.999), upper=c(rep(Inf,5),0.999))

AICs$AIC_REML[which(AICs$Model=="Gamma_Spatial")] <- opt_spatial$AIC

#------------------------------------------------------------------------------------------------------------
AICs$Delta_AIC_REML <- AICs$AIC_REML - min(AICs$AIC_REML)
AICs$Delta_AIC_ML <- AICs$AIC_ML - min(AICs$AIC_ML)
print(AICs)

#All spatial-temporal models ~ 8700-9000 points better than analogous nonspatial-temporal models (AIC_REML column)
#Normal likelihood appears to be most parsimonious (comparing ML AICs among spatial column)

#------------------------------------------------------------------------------------------------------------
#Cross validations
#h (block) cross validation:
if(FALSE){ #Takes ~ 50 hours to run
data$X_TTM_c <- data$X_TTM_c/1000
data$Y_TTM_c <- data$Y_TTM_c/1000
loc_xy <- unique(data[ ,c("X_TTM_c","Y_TTM_c") ] )

loc_xy <- as.data.frame(loc_xy)
names(loc_xy) <- c("X_TTM_c", "Y_TTM_c")

# Declare blocks for cross-validation
# Cahill picked these visually, could also use k-means.  These clusters seem reasonable.
loc_xy$Block = NA
loc_xy[which(loc_xy$X_TTM_c < 400 & loc_xy$Y_TTM_c > 6400),"Block" ] <- 1
loc_xy[which(loc_xy$X_TTM_c < 410 & loc_xy$Y_TTM_c < 6200 & is.na(loc_xy$Block)),"Block" ] <- 2
loc_xy[which(loc_xy$X_TTM_c < 450 & loc_xy$Y_TTM_c < 6000 & is.na(loc_xy$Block)),"Block" ] <- 2
loc_xy[which(loc_xy$X_TTM_c < 500 & loc_xy$Y_TTM_c > 6100 & is.na(loc_xy$Block)),"Block" ] <- 3
loc_xy[which(loc_xy$X_TTM_c < 660 & loc_xy$Y_TTM_c > 6000 & is.na(loc_xy$Block)),"Block" ] <- 4
loc_xy[which(loc_xy$X_TTM_c < 625 & loc_xy$Y_TTM_c < 6000 & is.na(loc_xy$Block)),"Block" ] <- 5
loc_xy[which(loc_xy$Y_TTM_c < 5800 & is.na(loc_xy$Block)),"Block" ] <- 6
loc_xy[which(loc_xy$X_TTM_c > 660 & loc_xy$Y_TTM_c > 5800 & is.na(loc_xy$Block)),"Block" ] <- 7

ggplot(loc_xy,aes(x=X_TTM_c,y=Y_TTM_c, color=as.factor(Block)))+geom_point()+
  theme_bw() + scale_colour_manual(values = RColorBrewer::brewer.pal(7, "Dark2") )

data <- dplyr::left_join(data, loc_xy, by=c("X_TTM_c", "Y_TTM_c"))
h_cvs <- matrix(NA, ncol=length(unique(data$Block)), nrow=nrow(AICs))
rownames(h_cvs) <- AICs$Model

for(model in unique(AICs$Model)){
for(k in unique(data$Block)){
  if(grepl("Norm", model)){CTL <- 1}
  if(grepl("Log", model)){CTL <- 2}
  if(grepl("Gamma", model)){CTL <- 3}

  Partition_i = ifelse(data$Block != k, 0, 1)
  if(grepl("Nonspatial", model)){
   data_nonspatial = list("Nobs" = nrow(data), "length_i" = data$TL, "age_i" = data$Age,
                          "lake_i" = data$Lake - 1,
                          "X_ij_omega" = model.matrix(~ -1 + data$wallEffDen.Std + data$compEffDen.Std +
                                                        data$GDD.Std  + data$wallEffDen.Std:data$compEffDen.Std),
                          "sex_i" = data$SexCode,"Nlakes" = length(unique(data$Lake)),"CTL" = CTL, "predTF_i"=Partition_i)

   parameters_nonspatial = list("ln_global_omega" = log(14),
                                "ln_global_linf" = log(45), "ln_sd_linf" = log(7), "global_tzero" = -1,
                                "ln_sd_tzero" = log(3), "ln_b_sex" = log(4.760871), "b_j_omega" = rep(0, ncol(data_nonspatial$X_ij_omega)),
                                "eps_omega" = rep(0, data_nonspatial$Nlakes ), "eps_linf" = rep(0, data_nonspatial$Nlakes ),
                                "eps_t0" = rep(0, data_nonspatial$Nlakes ),"ln_cv" = log(0.2), "ln_sd_omega" = log(4.5) )

   #Fit the nonspatial model:
   obj_nonspatial <- MakeADFun(data=data_nonspatial, parameters=parameters_nonspatial,
                               random=random_nonspatial, hessian=FALSE, DLL=VersionNonSpatial)

   opt_nonspatial = TMBhelper::Optimize(obj=obj_nonspatial,
                                        control=list(eval.max=1000, iter.max=1000),
                                        getsd=T, newtonsteps=1, bias.correct=F)

   rep <- obj_nonspatial$report()
  }
  if(grepl("Spatial", model)){
  #fit the spatial model:
   data_spatial = list("Nobs" = nrow(data), "length_i" = data$TL, "age_i" = data$Age,
                       "lake_i" = data$Lake - 1, "sex_i" = data$SexCode,
                       "X_ij_omega"= model.matrix(~ -1 + data$wallEffDen.Std + data$compEffDen.Std +
                                                    data$GDD.Std + data$wallEffDen.Std:data$compEffDen.Std),
                       "Nlakes" = length(unique(data$Lake)), "spdeMatrices" = spdeMatrices, "s_i" = data$Lake-1,
                       "t_i" = data$Year-1, "CTL" = CTL, "predTF_i"=Partition_i)

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

   obj_spatial <- MakeADFun(data=data_spatial, parameters=parameters_spatial,
                            random=random_spatial, hessian=FALSE, DLL=VersionSpatial)

   opt_spatial  = TMBhelper::Optimize(obj=obj_spatial,
                                      control=list(eval.max=1000, iter.max=1000),
                                      getsd=T, newtonsteps=1, bias.correct=F,
                                      lower=c(rep(-Inf,13),-0.999), upper=c(rep(Inf,13),0.999))

   rep <- obj_spatial$report()
  }
  #Record the results
  h_cvs[model, k] <- rep$pred_jnll / sum(Partition_i)
 }
  print(paste0(model , k))
}

AICs$h_block_cv <- rowMeans(h_cvs)
saveRDS(h_cvs, file="C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/analysis/h_cvs")

#--------------------------------------------------------------------------------------------------------------
#Leave one lake out cross validation:
lolo_cvs <- matrix(NA, ncol=length(unique(data$Lake)), nrow=nrow(AICs))
rownames(lolo_cvs) <- AICs$Model

for(model in unique(AICs$Model)){
  for(k in unique(data$Lake)){
    if(grepl("Norm", model)){CTL <- 1}
    if(grepl("Log", model)){CTL <- 2}
    if(grepl("Gamma", model)){CTL <- 3}

    Partition_i = ifelse(data$Lake==k,1,0)

    if(grepl("Nonspatial", model)){
      data_nonspatial = list("Nobs" = nrow(data), "length_i" = data$TL, "age_i" = data$Age,
                             "lake_i" = data$Lake - 1,
                             "X_ij_omega" = model.matrix(~ -1 + data$wallEffDen.Std + data$compEffDen.Std +
                                                           data$GDD.Std  + data$wallEffDen.Std:data$compEffDen.Std),
                             "sex_i" = data$SexCode,"Nlakes" = length(unique(data$Lake)),"CTL" = CTL, "predTF_i"=Partition_i)

      parameters_nonspatial = list("ln_global_omega" = log(14),
                                   "ln_global_linf" = log(45), "ln_sd_linf" = log(7), "global_tzero" = -1,
                                   "ln_sd_tzero" = log(3), "ln_b_sex" = log(4.760871), "b_j_omega" = rep(0, ncol(data_nonspatial$X_ij_omega)),
                                   "eps_omega" = rep(0, data_nonspatial$Nlakes ), "eps_linf" = rep(0, data_nonspatial$Nlakes ),
                                   "eps_t0" = rep(0, data_nonspatial$Nlakes ),"ln_cv" = log(0.2), "ln_sd_omega" = log(4.5) )

      #Fit the nonspatial model:
      obj_nonspatial <- MakeADFun(data=data_nonspatial, parameters=parameters_nonspatial,
                                  random=random_nonspatial, hessian=FALSE, DLL=VersionNonSpatial)

      opt_nonspatial = TMBhelper::Optimize(obj=obj_nonspatial,
                                           control=list(eval.max=1000, iter.max=1000),
                                           getsd=T, newtonsteps=1, bias.correct=F)

      rep <- obj_nonspatial$report()
    }
    if(grepl("Spatial", model)){
      #fit the spatial model:
      data_spatial = list("Nobs" = nrow(data), "length_i" = data$TL, "age_i" = data$Age,
                          "lake_i" = data$Lake - 1, "sex_i" = data$SexCode,
                          "X_ij_omega"= model.matrix(~ -1 + data$wallEffDen.Std + data$compEffDen.Std +
                                                       data$GDD.Std + data$wallEffDen.Std:data$compEffDen.Std),
                          "Nlakes" = length(unique(data$Lake)), "spdeMatrices" = spdeMatrices, "s_i" = data$Lake-1,
                          "t_i" = data$Year-1, "CTL" = CTL, "predTF_i"=Partition_i)

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

      obj_spatial <- MakeADFun(data=data_spatial, parameters=parameters_spatial,
                               random=random_spatial, hessian=FALSE, DLL=VersionSpatial)

      opt_spatial  = TMBhelper::Optimize(obj=obj_spatial,
                                         control=list(eval.max=1000, iter.max=1000),
                                         getsd=T, newtonsteps=1, bias.correct=F,
                                         lower=c(rep(-Inf,13),-0.999), upper=c(rep(Inf,13),0.999))

      rep <- obj_spatial$report()
    }
    #Record the results
    lolo_cvs[model, k] <- rep$pred_jnll / sum(Partition_i)
  }
  print(paste0(model , k))
}

AICs$lolo_cv <- rowMeans(lolo_cvs)
saveRDS(lolo_cvs, file="C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/analysis/lolo_cvs")
saveRDS(AICs, file="C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/analysis/AICs")

end_time <- Sys.time()

time <- end_time - start_time
time
}
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------

#Read in and reformat the information theoretic information for the real data:
lolo_cvs <- readRDS("C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/analysis/lolo_cvs")
h_cvs <- readRDS("C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/analysis/h_cvs")
AICs <- readRDS("C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/analysis/AICs")
AICs$h_block_cv <- rowMeans(h_cvs)
AICs$lolo_cv <- rowMeans(lolo_cvs)

#saveRDS(AICs, file="C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/analysis/InformationTheoreticsTable")
it <- readRDS("C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/analysis/InformationTheoreticsTable")
write.csv(it, "C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/analysis/InformationTheoreticsTable.csv")

#------------------------------------------------------------------------------------------------------------

#Run best nonspatial/spatial models with REML--LogNormal model (best as judged via PredNLL)
CTL <- 2 #lognormal
Partition_i <- rep(0, nrow(data)) #use all the data

Use_REML=T
if( Use_REML==TRUE ) random_spatial = union( random_spatial, c("ln_global_omega",
                                                               "ln_global_linf", "global_tzero",
                                                               "ln_b_sex", "b_j_omega") )

data_spatial = list("Nobs" = nrow(data), "length_i" = data$TL, "age_i" = data$Age,
                    "lake_i" = data$Lake - 1, "sex_i" = data$SexCode,
                    "X_ij_omega"= model.matrix(~ -1 + data$wallEffDen.Std + data$compEffDen.Std +
                                                 data$GDD.Std + data$wallEffDen.Std:data$compEffDen.Std),
                    "Nlakes" = length(unique(data$Lake)), "spdeMatrices" = spdeMatrices, "s_i" = data$Lake-1,
                    "t_i" = data$Year-1, "CTL" = CTL, "predTF_i"=Partition_i)

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


obj_spatial = MakeADFun(data=data_spatial, parameters=parameters_spatial,
                        random=random_spatial, hessian=FALSE, DLL=VersionSpatial)

opt_spatial  = TMBhelper::Optimize(obj=obj_spatial,
                                   control=list(eval.max=1000, iter.max=1000),
                                   getsd=T, newtonsteps=1, bias.correct=T,
                                   lower=c(rep(-Inf,5),-0.999), upper=c(rep(Inf,5),0.999))

SD = sdreport( obj_spatial )
final_gradient = obj_spatial$gr( opt_spatial$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

#Extract the intercept and SE
ParHat = as.list( opt_spatial$SD, "Estimate" )
SEHat  = as.list( opt_spatial$SD, "Std. Error" )

rep <- obj_spatial$report()

Range <- opt_spatial$SD["value"]$value[which(names(opt_spatial$SD["value"]$value)=="Range")]
rho <- opt_spatial$SD["value"]$value[which(names(opt_spatial$SD["value"]$value)=="rho")]
Range
rho


#nonspatial model:
dyn.load( dynlib(VersionNonSpatial) )

random_nonspatial = c("eps_linf", "eps_t0", "eps_omega")

data_nonspatial = list("Nobs"=nrow(data), "length_i"=data$TL, "age_i" = data$Age,
                       "lake_i" = data$Lake - 1, "X_ij_omega"= model.matrix(~ -1 + data$wallEffDen.Std +
                                                                              data$compEffDen.Std +
                                                                              data$GDD.Std  +
                                                                              data$wallEffDen.Std:data$compEffDen.Std),
                       "sex_i" = data$SexCode,"Nlakes" = length(unique(data$Lake)),"CTL" = CTL,
                       "predTF_i"=Partition_i)


parameters_nonspatial = list("ln_global_omega" = log(14),
                             "ln_global_linf" = log(45), "ln_sd_linf" = log(7), "global_tzero" = -1,
                             "ln_sd_tzero" = log(3), "ln_b_sex" = log(4.760871), "b_j_omega" = rep(0, ncol(data_nonspatial$X_ij_omega)),
                             "eps_omega" = rep(0, data_nonspatial$Nlakes ), "eps_linf" = rep(0, data_nonspatial$Nlakes ),
                             "eps_t0" = rep(0, data_nonspatial$Nlakes ),"ln_cv" = log(0.2), "ln_sd_omega" = log(4.5) )

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

SDnonspatial = sdreport( obj_nonspatial )
final_gradient = obj_nonspatial$gr( opt_nonspatial$par )
if( any(abs(final_gradient)>0.0001) | SDnonspatial$pdHess==FALSE ) stop("Not converged")

rep_nonspatial <- obj_nonspatial$report()


#Extract the intercept and SE
ParHat_nonspatial = as.list( opt_nonspatial$SD, "Estimate" )
SEHat_nonspatial  = as.list( opt_nonspatial$SD, "Std. Error" )

ParHat_nonspatial$b_j_omega
SEHat_nonspatial$b_j_omega


# Let's plot the results of the model, without and with
# the spatial correlation side by side.

NumberOfBetas <-  length(ParHat_nonspatial$b_j_omega)+2

Combined <- rbind(ParHat_nonspatial$b_j_omega + SEHat_nonspatial$b_j_omega%o%qnorm(c(0.025,0.50, 0.975)),
                  exp(ParHat_nonspatial$`ln_global_omega`) + exp(SEHat_nonspatial$`ln_global_omega`%o%qnorm(c(0.025,0.50, 0.975))),
                  exp(ParHat_nonspatial$ln_global_linf) + exp(SEHat_nonspatial$ln_global_linf%o%qnorm(c(0.025,0.50, 0.975))),

                  ParHat$b_j_omega + SEHat$b_j_omega%o%qnorm(c(0.025,0.50, 0.975)),
                  exp(ParHat$`ln_global_omega`) + exp(SEHat$`ln_global_omega`%o%qnorm(c(0.025,0.50, 0.975))),
                  exp(ParHat$ln_global_linf) + exp(SEHat$ln_global_linf%o%qnorm(c(0.025,0.50, 0.975))))
colnames(Combined) <- c("lower95%", "MLE", "upper95%")
Combined <- as.data.frame(Combined)

Combined$WhichModel <- rep(c("Nonspatial", "Spatial"),
                           each = NumberOfBetas)

Combined$WhichModel <- factor(Combined$WhichModel, levels=c(rep(c("Spatial", "Nonspatial"))))

Combined$WhichVariable <- rep(c("Intraspecific Density", "Interspecific Density", "GDD", "Density Interaction",
                                "Growth Rate Intercept", "Linf"), 2)

Combined$WhichVariable <- factor(Combined$WhichVariable, levels=c("Growth Rate Intercept", "Intraspecific Density",
                                                                  "Interspecific Density",  "Density Interaction", "GDD",  "Linf"))
colnames(Combined) <- c("Lo", "Mean", "Up", "WhichModel", "WhichVariable")
Combined


p <- ggplot()
p <- p + geom_point(data = Combined,
                    aes(x = WhichModel,
                        y = Mean)
)

p <- p + geom_errorbar(data = Combined,
                       aes(x = WhichModel,
                           ymax = Up,
                           ymin = Lo),
                       width=0.2)

p <- p + geom_abline(slope=0, intercept=0, linetype=3, size=1.35, colour="darkblue")

p <- p + xlab("Parameter") + ylab("Value")
p <- p + theme(text = element_text(size = 15))

p <- p + facet_wrap( ~ WhichVariable, scales = "free_y")
p <- p + theme(legend.position="none")
p

#ggsave("C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/plots/95_CI_Compare.png", p, dpi=600, width=8, height=8, units=c("in"))

#% of sites within the Range:
sum((dist(loc_xy)) <= Range)/ length(dist(loc_xy))

hist(dist(loc_xy), xlab="Distance between Lakes (km)", main="")
abline(v=Range, lty=3, lwd=3, col="Steelblue")

#Why are the walleye density coefficients different?

#Lake specific omega vs. effective density:
#declare parameters from fit
t0       = ParHat$global_tzero
global_linf = exp(ParHat$ln_global_linf)
b_sex = exp(ParHat$ln_b_sex)
b_j_omega = ParHat$b_j_omega
eps_t0 = ParHat$eps_t0
eps_omega_st_i = rep$eps_i
eps_linf = ParHat$eps_linf
global_omega = exp(ParHat$ln_global_omega)
CV = exp(ParHat$ln_cv)
data$eps <- rep$eps_i

omegas <- global_omega + b_j_omega[1]*data$wallEffDen.Std + eps_omega_st_i

data$omegasMLE <- omegas

#Make a sequence of effective densities
MyData <- data.frame(
  Density = seq(min(data$wallEffDen.Std), max(data$wallEffDen.Std),
                    length.out=100))

#Create X matrix
X <- model.matrix(~ Density, data = MyData)

Betas       <-c(global_omega, b_j_omega[1])
MyData$Pred <- X %*% Betas
vcov <- matrix(0, nrow=2, ncol=2)
diag(vcov) <- c(exp(SEHat$ln_global_omega) , SEHat$b_j_omega[1])
MyData$SE <- diag(X %*% vcov %*% t(X))

MyData$SeUp <- MyData$Pred + 1.96 * MyData$SE
MyData$SeLo <- MyData$Pred - 1.96 * MyData$SE

#Plot the predicted values
data$Block <- as.factor(data$Block)
p <- ggplot() + theme_classic()
p <- p + geom_point(data = data,
                    aes(y = omegasMLE, x = wallEffDen.Std, fill=Block),
                    size = 2, pch=21)
p <- p + ggtitle("Spatial Model") + ylim(5, 26)

p <- p + scale_fill_manual(values = RColorBrewer::brewer.pal(7, "Dark2") )

p <- p + ylab(expression(paste(paste("Growth Rate ", cm%.%year^{-1} ), (omega)))) +
  xlab("Effective Walleye Density (Standardized)")
p <- p + theme(text = element_text(size=15), plot.title = element_text(hjust = 0.5)) #, legend.position="none")
p <- p + geom_line(data = MyData,
                   aes(x = Density,
                       y = Pred),
                   linetype=2, size=2, colour="darkblue")

p <- p + geom_ribbon(data = MyData,
                     aes(x = Density,
                         ymax = SeUp,
                         ymin = SeLo ),
                     alpha = 0.6)

p <- p + theme(legend.title = element_blank(), legend.position = "none")
p

#Nonspatial model:
omegas <- exp(ParHat_nonspatial$ln_global_omega) + ParHat_nonspatial$b_j_omega[1]*data$wallEffDen.Std +
  ParHat_nonspatial$eps_omega[data$Lake]

data$OmegasMLE_nonspatial <- omegas

#Create X matrix
X <- model.matrix(~ Density, data = MyData)

Betas_nonspatial <- c(exp(ParHat_nonspatial$ln_global_omega), ParHat_nonspatial$b_j_omega[1])
MyData$Pred_nonspatial <- X %*% Betas_nonspatial
vcov <- matrix(0, nrow=2, ncol=2)
diag(vcov) <- c(exp(SEHat_nonspatial$ln_global_omega) , SEHat_nonspatial$b_j_omega[1])
MyData$SE_nonspatial <-  diag(X %*% vcov %*% t(X))

MyData$SeUp_nonspatial <- MyData$Pred_nonspatial + 1.96 * MyData$SE_nonspatial
MyData$SeLo_nonspatial <- MyData$Pred_nonspatial - 1.96 * MyData$SE_nonspatial

p1 <- ggplot() + theme_classic( )
p1 <- p1 + geom_point(data = data,
                    aes(y = OmegasMLE_nonspatial, x = wallEffDen.Std, fill=Block),
                    size = 2, pch=21)
p1 <- p1 + scale_fill_manual(values = RColorBrewer::brewer.pal(7, "Dark2") )

p1 <- p1 + ylab(expression(paste(paste("Growth Rate ", cm%.%year^{-1} ), (omega)))) +
  xlab("Effective Walleye Density (Standardized)") + ylim(5, 26)
p1 <- p1 + theme(text = element_text(size=15), plot.title = element_text(hjust = 0.5)) #, legend.position="none")
p1 <- p1 + geom_line(data = MyData,
                   aes(x = Density,
                       y = Pred_nonspatial),
                   linetype=2, size=2, colour="darkblue")
p1 <- p1 + ggtitle("Nonspatial Model")

p1 <- p1 + geom_ribbon(data = MyData,
                     aes(x = Density,
                         ymax = SeUp_nonspatial,
                         ymin = SeLo_nonspatial ),
                     alpha = 0.6)

p1 <- p1 + theme(legend.title = element_blank(), legend.position = "none")

# p1 <- p1 + facet_wrap( ~ Year)

p1

mapbox <- c(-120, 49, -109.99, 60)
glgmap   <- get_stamenmap(bbox=mapbox, maptype="toner-hybrid", zoom=6, crop=T, force=T)
map <- ggmap(glgmap) + geom_point(aes(Long_c, Lat_c, fill=Block), pch=21, size=2, data=data ) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(7, "Dark2"), name="Spatial \n Block" )
map <- map + ylab("Latitude") + xlab("Longitude")
map <- map + theme_classic( ) + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.75),
                                      axis.title=element_text(size=15))

big <- ggarrange(ggarrange(p, p1, ncol=1, nrow=2), map)
#ggsave("C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/plots/Multipanel.png", big,
#dpi=1200, width=11, height=8, units=c("in"))

p3 <- gridExtra::grid.arrange(p, p1, ncol=1)

#extract the SEs of the latent field
for(i in 1:nrow(data)){
  data$se_omega[i] <- SEHat$eps_omega_st[data_spatial$s_i[i]+1, data_spatial$t_i[i]+1]
}

d <- mutate(data, WBID=data$WBID, x = data$Long_c,
            y = data$Lat_c, year = data$Year + 2000,
            eps = data$eps,
            se=data$se_omega)%>% group_by(WBID, year)


#Plot eps_i (st raneffs of the lakes) through time x Space:
mapbox <- c(-120.5, 48.95, -109.5, 60.5)

glgmap   <- get_stamenmap(bbox=mapbox, maptype="terrain", zoom=4, crop=T, force=T)

map <- ggmap(glgmap) + geom_point(aes(Long_c, Lat_c, fill=se), pch=21, size=2, data=d ) +
  #scale_fill_gradient2(name="Random \n Effect" )
  scale_fill_gradient2(name="Standard \n Error", midpoint=median(d$se) )

map <- map + ylab("Latitude") + xlab("Longitude")
map <- map + theme_classic( ) + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.75),
                                      axis.title=element_text(size=15),
                                      panel.spacing.x = unit(7.5, "mm"))

map <- map + facet_wrap(~year, nrow=3)
map

#ggsave("C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/plots/SpaceTimeFieldTerrainSE.png",
 #      map, dpi=1200, width=13, height=10, units=c("in"))

map <- qmap('alberta', zoom = 5, maptype = 'satellite') + geom_point(aes(Long_c, Lat_c, fill=eps), pch=21, size=2, data=d ) +
  scale_fill_gradient2(name="Random \n Effect" )
map <- map + ylab("Latitude") + xlab("Longitude")
map <- map + theme_classic( ) + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.75),
                                      axis.title=element_text(size=15),
                                      panel.spacing.x = unit(7.5, "mm"))

map <- map + facet_wrap(~year, nrow=3)
map

# ggsave("C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/plots/SpaceTimeFieldGoogle.png",
#        map, dpi=1200, width=11, height=8, units=c("in"))

#png(file="C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/plots/Mesh.png",
 #   width=8,height=11,units="in",res=1200)
plot(mesh)
points(loc_xy, col="Steelblue", pch=1)
#dev.off()

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

#Generate/plot the predictive distribution conditional on the best parameter estimates a la Gelman:
rlnorm0 <- function(n = 1, mean, coeffOfVar)
{
  sigma <- sqrt(log(coeffOfVar^2 + 1))
  mu <- log(mean) - sigma^2 / 2
  rlnorm(n, mu, sigma)
}
#Takes about a minute to simulate
Nsim      <- 10000
Age_Seq <- 0:26
Stats     <- matrix(0, nrow=length(Age_Seq), ncol=5)
pred_dist <- matrix(0, nrow=nrow(data), ncol=Nsim)
set.seed(1)

#declare parameters from fit
t0       = ParHat$global_tzero
global_linf = exp(ParHat$ln_global_linf)
b_sex = exp(ParHat$ln_b_sex)
b_j_omega = ParHat$b_j_omega
eps_t0 = ParHat$eps_t0
eps_omega_st_i = rep$eps_i
eps_linf = ParHat$eps_linf
global_omega = exp(ParHat$ln_global_omega)
CV = exp(ParHat$ln_cv)

#create the deviates
for(i in 1:nrow(data)) {
  Lake <- data$Lake[i]

  tzero = t0 + eps_t0[Lake]

  omega = global_omega +
    b_j_omega[1]*data$wallEffDen.Std[i] +
    b_j_omega[2]*data$compEffDen.Std[i] +
    b_j_omega[3]*data$GDD.Std[i] +
    b_j_omega[4]*data$wallEffDen.Std[i]*data$compEffDen.Std[i] +
    eps_omega_st_i[i]

  linf = global_linf + b_sex*data$SexCode[i] + eps_linf[Lake]

  lpred = linf*(1-exp(-(omega/linf) * (data$Age[i] - t0 )))
  sigma = CV*lpred
  if(data_spatial$CTL==1) {pred_dist[i, 1:Nsim] = rnorm(Nsim, lpred, sigma)}
  if(data_spatial$CTL==2) {pred_dist[i, 1:Nsim] = rlnorm0(Nsim, lpred, CV)}
  if(data_spatial$CTL==3) {pred_dist[i, 1:Nsim] = rgamma(Nsim, shape=1/CV^2, scale=lpred*CV^2)}

  print(i)
}

#calculate the relevant quantiles for plotting
for(I in Age_Seq){
  sub.dat <- pred_dist[which(data$Age==I), ]
  Stats[I+1,] <- quantile(sub.dat,prob=c(0.025,0.25,0.5,0.75,0.975))
  print(I)
}

#Put it all in the main data frame:
Stats <- cbind(Age_Seq, Stats)
Stats <- as.data.frame(Stats)
colnames(Stats) <- c("Age", "fit.025", "fit.25", "fit.5", "fit.75", "fit.975")

#join on Age:
data <- left_join(data, Stats)

#plot the predictive distribution w/ quantiles:
p <- ggplot(data = data, aes(y = TL, x = Age, fill=Sex))
p <- p + xlab("Age") + ylab("Total Length (cm)")
p <- p + geom_point(shape = 21, size = 1.0, position="jitter") +
   scale_fill_manual(values=c("white", "black"))

p <- p + geom_ribbon(aes(x = Age, ymax = fit.975, ymin = fit.025), fill = grey(0.5), alpha = 0.5)
p <- p + geom_ribbon(aes(x = Age, ymax = fit.25, ymin = fit.75), fill = grey(1), alpha = 0.8)
p <- p + geom_line(aes(x = Age, y = fit.5), linetype = "dashed",  col="black", size=0.75)
p <- p + scale_x_continuous(limits = c(-0.5, 27), breaks=c(0,6,11,16,21,26)) +
  scale_y_continuous(limits = c(0, 86), breaks=c(0,20,40,60,80))
p <- p + theme_classic() +
  theme(axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.text=element_text(size=12, colour='black'),
        axis.title=element_text(size=15) )
p

# ggsave("C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/plots/LognormalPredictiveDistribution.png",
#        p, dpi=1200, width=11, height=8, units=c("in"))

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

#quantify coverage:

quantvec <- NULL
for(i in 0:26){
  qs <- quantile(pred_dist[which(data$Age == i)], c(0.025, 0.975))
  sub.data <- data[which(data$Age == i), "TL"]
  prop <- sum(sub.data >= qs[1] & sub.data <= qs[2])/length(sub.data)
  quantvec[i+1] <- prop
}

print(paste( round(100*sum(quantvec) / length(quantvec), 2), "% coverage within 95% C.I.s", sep=""))

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

#lake-year-specific mle predictions
#circles = girls, squares = boys
Age_Seq <- 0:26
pdf("C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/plots/LognormalLakeYearPredictions.pdf",
    width=8, height=11)
par(mfrow=c(3,3))
for(i in unique(data$WBID)){
  sub.dat <- data[which(data$WBID==i),]
  Lake <- unique(sub.dat$Lake)
  Name <- unique(sub.dat$Name)
  for(j in unique(sub.dat$Year)) {
    sub.sub.dat <- sub.dat[which(sub.dat$Year==j),]
    tzero = t0 +
      eps_t0[Lake]

    omega = global_omega +
      b_j_omega[1]*sub.sub.dat$wallEffDen.Std[1] +
      b_j_omega[2]*sub.sub.dat$compEffDen.Std[1] +
      b_j_omega[3]*sub.sub.dat$GDD.Std[1] +
      b_j_omega[4]*sub.sub.dat$wallEffDen.Std[1]*sub.sub.dat$compEffDen.Std[1] +
      sub.sub.dat$eps[1]

    linf = global_linf +
      b_sex*c(0,1) +
      eps_linf[Lake]

    lpred_m <- lpred_f <- NA
    for(a in 1:length(Age_Seq)){
      lpred_m[a] = linf[1]*(1-exp(-(omega/linf[1]) * (Age_Seq[a] - t0 )))
      lpred_f[a] = linf[2]*(1-exp(-(omega/linf[2]) * (Age_Seq[a] - t0 )))
    }

    plot(lpred_f~Age_Seq, type="l", lty=2, col="black", ylim=c(0,85), main=paste(Name, j, sep= " " ),
         xlab="Age (Years)", ylab="Total Length (cm)", cex.main=1, lwd=1.5)
    points(lpred_m~Age_Seq, type="l", col="black", lwd=1.5)
    points(sub.sub.dat$TL~sub.sub.dat$Age, pch=sub.sub.dat$SexCode)
    text(x=20, y=25, labels = paste0("Linf = ", format(round(linf[2],2), nsmall=2)))
    text(x=20, y=15, labels = paste0("Omega = ", format(round(omega,2), nsmall=2)))
    text(x=20, y=5, labels = paste0("T0 = ", format(round(t0,2), nsmall=2)))

}}
dev.off()

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

#Let's do some more residual checks:

#Some ugly residual plots for what they are worth:
rep <- obj_spatial$report()
length_pred <- rep$length_pred
data$E1 <- (data$TL - length_pred)

#png("C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/plots/residuals.png",
 #   width=8,height=11,units="in",res=1200)

#Residuals by lake:
xyplot(data$E1 ~ rep$length_pred | data$Name,
       panel=function(x, y){
         panel.xyplot(x, y)
         #panel.lmline(x, y, lty = 2)
       } )
#dev.off()

#Residuals by year:
xyplot(data$E1 ~ rep$length_pred | data$Year,
       panel=function(x, y){
         panel.xyplot(x, y)
       } )

# Homogeneity
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x=length_pred, y=data$E1, ylab="Residual", xlab="Length")
abline(h = 0, v = 1)

# Normality
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
hist(data$E1, main="", xlab="Residual")

#png("C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/plots/residuals2.png",
#    width=8,height=11,units="in",res=1200)

# Independence due to model misfit
#par(mfrow = c(3,2), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = data$wallEffDen, y = data$E1, xlab="Intraspecific Density", ylab="Residual")
abline(h = 0)

plot(x = data$compEffDen, y = data$E1, xlab="Interspecific Density", ylab="Residual")
abline(h = 0)

plot(x = data$GDD, y = data$E1, xlab="Growing Degree Days", ylab="Residual")
abline(h = 0)

plot(x = data$Year, y = data$E1, xlab="Year", ylab="Residual")
abline(h = 0)

#Residuals are getting a little broader through time, but
#also more samples collected through time...

#residuals by sex
boxplot(data$E1~data$SexCode, xlab="Sex")
abline(h = 0)

#dev.off()

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#Look at the covariance matrix:

scaledCovar <- cov2cor(opt_spatial$SD$cov.fixed) #scales the covariance matrix
colnames(scaledCovar) <- rownames(scaledCovar) <- rownames(opt_spatial$SD$cov.fixed)
print(round(scaledCovar, 2))

#-------------------------------------------------------------------------------------------------------------------
#End End End
#-------------------------------------------------------------------------------------------------------------------



