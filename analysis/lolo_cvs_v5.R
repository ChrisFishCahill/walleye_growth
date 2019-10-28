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
CVs <- data.frame(Model=c("Normal_Nonspatial", "Normal_Nonspatial_Reduced",
                          "Lognormal_Nonspatial", "Lognormal_Nonspatial_Reduced",
                          "Gamma_Nonspatial", "Gamma_Nonspatial_Reduced",
                          "Normal_Spatial", "Normal_Spatial_Reduced",
                          "Lognormal_Spatial", "Lognormal_Spatial_Reduced",
                          "Gamma_Spatial", "Gamma_Spatial_Reduced"),
                  lolo_cv=NA)

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#Create inputs for the Spatial models:
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

#Inputs for the nonspatial model
setwd("C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/executables")
VersionNonSpatial = "vb_nonspatial"
compile(paste0(VersionNonSpatial, ".cpp"))
dyn.load( dynlib(VersionNonSpatial) )

random_nonspatial = c("eps_linf", "eps_t0", "eps_omega")

#h (block) cross validation:
data$X_TTM_c <- data$X_TTM_c/1000
data$Y_TTM_c <- data$Y_TTM_c/1000
loc_xy <- unique(data[ ,c("X_TTM_c","Y_TTM_c") ] )

loc_xy <- as.data.frame(loc_xy)
names(loc_xy) <- c("X_TTM_c", "Y_TTM_c")


#Leave one lake out cross validation:
lolo_cvs <- matrix(NA, ncol=length(unique(data$Lake)), nrow=nrow(CVs))
rownames(lolo_cvs) <- CVs$Model

for(model in unique(CVs$Model)){
  for(k in unique(data$Lake)){
    if(grepl("Norm", model)){CTL <- 1}
    if(grepl("Log", model)){CTL <- 2}
    if(grepl("Gamma", model)){CTL <- 3}

    Partition_i = ifelse(data$Lake==k,1,0)
    if(grepl("Nonspatial", model)){
    if(!grepl("Reduced", model)){
      data_nonspatial = list("Nobs" = nrow(data), "length_i" = data$TL, "age_i" = data$Age,
                             "lake_i" = data$Lake - 1,
                             "X_ij_omega" = model.matrix(~ -1 + data$wallEffDen.Std + data$compEffDen.Std +
                                                           data$GDD.Std  + data$wallEffDen.Std:data$compEffDen.Std),
                             "sex_i" = data$SexCode,"Nlakes" = length(unique(data$Lake)),"CTL" = CTL, "predTF_i"=Partition_i)
    } else {
      data_nonspatial = list("Nobs" = nrow(data), "length_i" = data$TL, "age_i" = data$Age,
                             "lake_i" = data$Lake - 1,
                             "X_ij_omega" = model.matrix(~ -1 + data$wallEffDen.Std + data$compEffDen.Std +
                                                           data$GDD.Std),
                             "sex_i" = data$SexCode,"Nlakes" = length(unique(data$Lake)),"CTL" = CTL, "predTF_i"=Partition_i)
    }

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
    if(!grepl("Reduced", model)){
      #fit the spatial model:
      data_spatial = list("Nobs" = nrow(data), "length_i" = data$TL, "age_i" = data$Age,
                          "lake_i" = data$Lake - 1, "sex_i" = data$SexCode,
                          "X_ij_omega"= model.matrix(~ -1 + data$wallEffDen.Std + data$compEffDen.Std +
                                                       data$GDD.Std + data$wallEffDen.Std:data$compEffDen.Std),
                          "Nlakes" = length(unique(data$Lake)), "spdeMatrices" = spdeMatrices, "s_i" = data$Lake-1,
                          "t_i" = data$Year-1, "CTL" = CTL, "predTF_i"=Partition_i)
    } else {

      data_spatial = list("Nobs" = nrow(data), "length_i" = data$TL, "age_i" = data$Age,
                          "lake_i" = data$Lake - 1, "sex_i" = data$SexCode,
                          "X_ij_omega"= model.matrix(~ -1 + data$wallEffDen.Std + data$compEffDen.Std +
                                                       data$GDD.Std),
                          "Nlakes" = length(unique(data$Lake)), "spdeMatrices" = spdeMatrices, "s_i" = data$Lake-1,
                          "t_i" = data$Year-1, "CTL" = CTL, "predTF_i"=Partition_i)
    }


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
                                "logit_rho" = 3.12 ) #2 * plogis(3) - 1 = 0.9

      obj_spatial <- MakeADFun(data=data_spatial, parameters=parameters_spatial,
                               random=random_spatial, hessian=FALSE, DLL=VersionSpatial)

      opt_spatial  = TMBhelper::Optimize(obj=obj_spatial,
                                         control=list(eval.max=1000, iter.max=1000),
                                         getsd=T, newtonsteps=1, bias.correct=F)

      rep <- obj_spatial$report()
    }
    #Record the results
    lolo_cvs[model, k] <- rep$pred_jnll / sum(Partition_i)
  }
  print(paste0(model , k))
}

CVs$lolo_cv <- rowMeans(lolo_cvs)
saveRDS(lolo_cvs, file="C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/analysis/lolo_cvs")
#saveRDS(AICs, file="C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/analysis/AICs")

end_time <- Sys.time()

time <- end_time - start_time
time
