#-----------------------------
#Simulation testing von B with and without AR-1 spatial-temporal dependency
#Coded by Cahill April 2019
#-----------------------------

#Install and load Libraries
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
    Eta_it[,t] =  beta0  + betas[1]*X_ij[ ,t, 1] + betas[2]*X_ij[ ,t, 2] +
      betas[3]*X_ij[ ,t, 3] + eps_st[[t]] #fixed effects + Temporally evolving spatial field
  }

  #Exclude exceedingly low simulated omegas (< 5 cm/year)
  Eta_it[which(Eta_it < 5.0)] <- 5.0
  #could instead use log-linear regression, or change SD_O parameter

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
       if(Likelihood=="Gamma"){Simulated_Length = rgamma(Nfish, shape=1/cv^2, scale=lpreds*cv^2)} }  else {
       Simulated_Length <- lpreds
      }

      Tmp = data.frame("Lake"=rep(s, Nfish) , "Year"=rep(t, Nfish), "x1"=rep(X_ij[s,t,1], Nfish),
                       "x2"=rep(X_ij[s,t,2], Nfish), "x3"=rep(X_ij[s,t,3], Nfish), "Age" = Ages,
                       "Simulated_Length"= Simulated_Length, "Lpred"=lpreds)
      DF = rbind(DF, Tmp)
    }}

  DF = cbind( DF, 'Longitude'=Loc[DF[,'Lake'],1], 'Latitude'=Loc[DF[,'Lake'],2] )
  DF = data.frame(DF, row.names=NULL)

  if(Design=="Unbalanced"){
   #Select lakes to keep
   lakes_to_keep <- sample(unique(DF$Lake), lakes_good_data*length(unique(DF$Lake)), replace=FALSE)
   DF_sentinel <- DF[which(DF$Lake %in% lakes_to_keep),]

   #lakes to reduce
   lakes_to_reduce <- unique(DF$Lake)[!unique(DF$Lake) %in% lakes_to_keep]
   DF_to_reduce <- DF[which(DF$Lake %in% lakes_to_reduce),]

   #Amount to reduce
   N_data_reduced <- round(nrow(DF_sentinel)/(percent_of_data*100)*100, 0)
   reduced_to_keep <- sample(1:nrow(DF_to_reduce), N_data_reduced, replace=FALSE)
   reduced_to_drop <- setdiff(1:nrow(DF_to_reduce),reduced_to_keep)

   DF_to_reduce[reduced_to_drop,'Simulated_Length'] = NA
   DF = rbind(DF_sentinel, DF_to_reduce)
  }

  DF = filter(DF, Lpred > 0)
  DF = filter(DF, Simulated_Length > 0)

  # Return stuff
  Sim_List = list("DF"=DF, "Loc"=Loc, "Eta_it"=Eta_it,
                  "eps_st"=eps_st, "n_years"=n_years,
                  "n_stations"=n_stations, "eps_t0"=Epsilon_t0, "eps_linf"=Epsilon_Linf)
  Sim_List[["Parameters"]] = c('SpatialScale'=SpatialScale, 'SigmaO'=SD_O,
                               'rho'=rho, "beta0"=beta0, "beta_ij"=betas,
                               'Linf' = Global_Linf, "t0" = Global_t0,
                               'SD_t0' = SD_t0,'SD_linf'=SD_linf)
  Sim_List
}

#-----------------------------------------
#Run the simulation and estimate models
#-----------------------------------------

#Sampling and Likelihood parameters:
Nsim <- 100
Age_Range <- 0:25
Nfish <- 50
N_years <- c(4,7,10,15,20)
N_lakes <- 50

Likelihood=c("Normal", "Lognormal", "Gamma")
Likelihood = Likelihood[2]

Designs = c("Balanced", "Unbalanced")
Designs = Designs[2]
lakes_good_data <- 0.1 #percent of lakes with good data
percent_of_data  <- 0.35 #percent of total dataset represented by lakes_good_data

#von B parameters:
Global_Linf  <- 55.70
SD_linf      <- 7.41
Global_t0    <- -1
SD_t0        <- 0.3
betas        <- c(-1,0,1)
beta0        <- 14.79
cv           <- 0.069

#Spatial-temporal Parameters
kappa        <- 0.5
SpatialScale <- 1/kappa
Range        <- sqrt(8) / kappa
SD_O         <- 4.66
rho          <- 0.5

# Compile & load
setwd("C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/sim/")
VersionSpatial = "vb_spdeXar1"
VersionNonSpatial = "vb_nonspatial"
compile( paste0(VersionSpatial,".cpp") )
compile( paste0(VersionNonSpatial,".cpp") )
dyn.load( dynlib(VersionSpatial) )
dyn.load( dynlib(VersionNonSpatial) )

seed <-  sample.int(1e6, 1)
set.seed( seed ) #416009, 345614

#Run the simulation:
ptm <- proc.time()
for(design in Designs){
  for(nyears in unique(N_years)){
    for(nlakes in unique(N_lakes)){
    replicate=1
    while(replicate < (Nsim+1) ){
    Sim_List = Sim_Fn( n_years=nyears, n_stations=N_lakes, SpatialScale=SpatialScale, SD_O=SD_O,
 	                     rho=rho, betas=betas, beta0=beta0, cv=cv, Likelihood=Likelihood, Design=design )

    DF = Sim_List[["DF"]]
    loc_xy_orig = loc_xy = unique(DF[,c("Longitude", "Latitude")])

    mesh = inla.mesh.create( loc_xy, refine=TRUE, extend=-0.5, cutoff=0.01 )
    spde = inla.spde2.matern( mesh, alpha=2 )

    #------Plots--------
    # plot(mesh)
    # points(loc_xy_orig, pch=16, col="Steelblue")
    #
    # ggplot(DF, aes(Age, Simulated_Length)) + geom_point() + facet_wrap(~Lake)
    #
    # d <- reshape2::melt(Sim_List[["eps_st"]]) %>%
    #        dplyr::mutate(x = rep(Sim_List[["Loc"]][,"x"], nyears),
    #        y = rep(Sim_List[["Loc"]][,"y"], nyears))
    #
    # ggplot(d, aes(x, y, col = value)) + geom_point() +
    #  facet_wrap(~L1) +
    #  scale_color_gradient2()
    # -------------------

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
                        "spdeMatrices" = spdeMatrices, "CTL" = CTL, "predTF"=rep(0, nrow(DF)),
                        "s_i" = DF$Lake-1,
                        "t_i" = DF$Year-1 )

    saveRDS(data_spatial, file="checker.RDA") #record in case simulation dies

    data_nonspatial = list("Nobs"=nrow(DF), "length_i"=DF$Simulated_Length, "age_i" = DF$Age,
                           "lake_i" = DF$Lake - 1, "X_ij_omega"= model.matrix(~ -1 + DF$x1 + DF$x2 + DF$x3),
                           "Nlakes" =  length(unique(DF$Lake)),"CTL" = CTL, "predTF"=rep(0, nrow(DF)) )

    parameters_spatial = list("ln_global_omega" = log(beta0), "ln_global_linf" = log(Global_Linf),
                               "ln_sd_linf" = log(SD_linf), "global_tzero" = Global_t0, "ln_sd_tzero" = log(SD_t0),
                               "b_j_omega" = betas, "eps_omega_st" = matrix(0,  nrow=mesh$n,ncol=nyears ),
                               "eps_linf" = rep(0, data_spatial$Nlakes),"eps_t0" = rep(0, data_spatial$Nlakes),
                               "ln_cv" = log(cv), "ln_kappa" = log(kappa), "ln_tau_O" =  -2, "rho" = rho )

    parameters_nonspatial = list("ln_global_omega" = log(beta0),
                                 "ln_global_linf" = log(Global_Linf), "ln_sd_linf" = log(SD_linf), "global_tzero" = Global_t0,
                                 "ln_sd_tzero" = log(SD_t0),"b_j_omega" = betas, "eps_omega" = rep(0, data_nonspatial$Nlakes ),
                                 "eps_linf" = rep(0, data_nonspatial$Nlakes), "eps_t0" = rep(0, data_nonspatial$Nlakes),
                                 "ln_cv" = log(cv), "ln_sd_omega" = log(SD_O) )

    random_spatial = c("eps_linf", "eps_t0", "eps_omega_st")
    random_nonspatial = c("eps_linf", "eps_t0", "eps_omega")

    obj_spatial <- MakeADFun(data=data_spatial, parameters=parameters_spatial,
                             random=random_spatial, hessian=FALSE, DLL=VersionSpatial)

    opt_spatial <- tryCatch(TMBhelper::Optimize(obj=obj_spatial,control=list(eval.max=1000, iter.max=1000),
                                               getsd=T, newtonsteps=1, bias.correct=T,
                                               lower=c(rep(-Inf,11),-0.999), upper=c(rep(Inf,11),0.999)),
                            error = function(e) print(e) )

    obj_nonspatial <- MakeADFun(data=data_nonspatial, parameters=parameters_nonspatial,
                                random=random_nonspatial, hessian=FALSE, DLL=VersionNonSpatial)

    opt_nonspatial <- tryCatch(TMBhelper::Optimize(obj=obj_nonspatial,
                                                  control=list(eval.max=1000, iter.max=1000),
                                                  getsd=T, newtonsteps=1, bias.correct=F),
                              error = function(e) print(e))

      #if the estimation behaves,  save results & advance loop:
      if(!inherits(opt_spatial, "error") && !inherits(opt_nonspatial, "error") &&
         is.null(opt_spatial$h) && is.null(opt_nonspatial$h)){

        SD = sdreport( obj_spatial )
        spatial_report = obj_spatial$report()

        SD_nonspatial = sdreport( obj_nonspatial )
        nonspatial_report = obj_nonspatial$report()

        sim_rep <- list()
        sim_rep[["sim_data"]] = Sim_List

        spatial_summary = summary(SD)
        opt_spatial[["summary"]] = spatial_summary
        opt_spatial[["report"]] = spatial_report
        sim_rep[["opt_spatial"]] = opt_spatial

        nonspatial_summary = summary(SD_nonspatial)
        opt_nonspatial[["summary"]] = nonspatial_summary
        opt_nonspatial[["report"]] = nonspatial_report
        sim_rep[["opt_nonspatial"]] = opt_nonspatial

        file_name <- paste(design, paste(paste0(paste0(nyears, "Years"), nlakes, "Lakes"), paste(replicate,"Replicate.RData",sep=""), sep="_"), sep="_")
        save(sim_rep, file=file_name)
        rm("sim_rep", "data_spatial", "data_nonspatial", "obj_spatial", "obj_nonspatial", "Sim_List")
        print(paste(paste("Simulation", file_name, sep=" "), "Complete", sep=" "))
        replicate <- replicate + 1

        #remove object if replicate fishies:
        if (file.exists("checker.RDA")) file.remove("checker.RDA")
      } #try
    } #while
    } #nlakes
  } #nyears
} #designs

proc.time() - ptm

#Loop through working directory and re-organize results into an array:
Models <- c("Spatial", "Nonspatial")

sim_results <- array(NA, dim=c(length(Models), length(Designs), length(N_years), length(N_lakes), Nsim, 14, 3),
                     dimnames=list(c("Spatial","Nonspatial"), Designs, c(N_years),
                                   c(N_lakes), paste("Sim=",1:Nsim,sep=""),
                                   c("ln_global_omega", "ln_global_linf", "ln_sd_linf", "global_tzero", "ln_sd_tzero", "b_j_omega1",
                                     "b_j_omega2", "b_j_omega3", "ln_cv", "ln_kappa", "ln_tau_O", "rho", "ln_sd_omega", "sigmaO"),
                                   c("10%","50%","90%")))

for(design in Designs){
for(nyears in N_years){
for(nlakes in N_lakes){
for(rep in 1:Nsim){
  file_name <- paste(design, paste(paste0(paste0(nyears, "Years"), nlakes, "Lakes"), paste(rep,"Replicate.RData",sep=""), sep="_"), sep="_")
  load( file_name )

  #spatial model:
  sim_results[Models[1], design, as.character(nyears), as.character(nlakes),
              rep,c("ln_global_omega", "ln_global_linf", "ln_sd_linf", "global_tzero", "ln_sd_tzero", "b_j_omega1",
                          "b_j_omega2", "b_j_omega3", "ln_cv", "ln_kappa", "ln_tau_O", "rho", "sigmaO"),
              c("10%","50%","90%")] = c(sim_rep$opt_spatial$summary[1:length(sim_rep$opt_spatial$par),'Estimate'],
                                        sim_rep$opt_spatial$summary[which(row.names(sim_rep$opt_spatial$summary)=="SigmaO"),"Estimate"])%o%rep(1,3) +
              c(sim_rep$opt_spatial$summary[1:length(sim_rep$opt_spatial$par),'Std. Error'],
              sim_rep$opt_spatial$summary[which(row.names(sim_rep$opt_spatial$summary)=="SigmaO"),"Std. Error"])%o%qnorm(c(0.1,0.5,0.9))

  #nonspatial model:
  sim_results[Models[2], design, as.character(nyears), as.character(nlakes),
              rep, c("ln_global_omega", "ln_global_linf", "ln_sd_linf", "global_tzero", "ln_sd_tzero", "b_j_omega1",
                           "b_j_omega2", "b_j_omega3", "ln_cv", "ln_sd_omega"),
              c("10%","50%","90%")] = sim_rep$opt_nonspatial$summary[1:length(sim_rep$opt_nonspatial$par),'Estimate']%o%rep(1,3) +
              sim_rep$opt_nonspatial$summary[1:length(sim_rep$opt_nonspatial$par), 'Std. Error']%o%qnorm(c(0.1,0.5,0.9))
}}}}

#Plot spatial vs. nonspatial:

truth = c("ln_global_omega"=log(beta0), "ln_global_linf"=log(Global_Linf), "ln_sd_linf"=log(SD_linf),
          "global_tzero"=Global_t0, "ln_sd_tzero"=log(SD_t0), "ln_cv"=log(cv), "b_j_omega1"=betas[1],
          "b_j_omega2"=betas[2], "b_j_omega3"=betas[3], "ln_kappa"=log(kappa),
          "rho"=rho, "sigmaO"=SD_O)

d <- reshape2::melt(sim_results[c("Spatial", "Nonspatial"),Designs,as.character(N_years),
                                as.character(nlakes), 1:Nsim,
                                c("ln_global_omega", "ln_global_linf", "ln_sd_linf",
                                  "global_tzero", "ln_sd_tzero", "ln_cv", "b_j_omega1",
                                  "b_j_omega2", "b_j_omega3"), "50%"],
                    varnames=c("Model", "Design", "Nyears", "Replicate",
                               "Parameter"), value.name="MLE") %>%
              mutate(Truth=truth[Parameter]) %>%
              mutate(MLE=ifelse(grepl("ln", Parameter), exp(MLE), MLE)) %>%
              mutate(Truth=ifelse(grepl("ln", Parameter), exp(Truth), Truth)) %>%
              mutate(Parameter=stringr::str_remove(Parameter, "ln_"))
d$Nyears <- as.factor(d$Nyears)

labels <- c(global_omega = "Omega Intercept", cv = "Coefficient of Variation",
            global_linf= "Linf Intercept", b_j_omega1="Beta 1", b_j_omega2="Beta 2",
            b_j_omega3="Beta 3", sd_linf="Linf StDev", global_tzero="T0",
            sd_tzero="T0 StDev")

#Balanced data
p <- ggplot(subset(d, Design %in% "Balanced"), aes(x=Nyears, y=MLE, fill=Model)) +
     geom_boxplot(outlier.alpha=0.2) + scale_fill_manual(values = c(rgb(0,0,0,0.3), rgb(1,1,1,0.3))) +
     xlab("Number of Years") + ylab("Parameter Value") + guides(aes(shape=NA)) +
     labs(fill = "")  + facet_wrap(~Parameter, scales=c("free_y"), labeller=labeller(Parameter = labels)) +
     scale_x_discrete(limits=levels(d$Nyears)) + theme(legend.key=element_blank())

p <- p + ggtitle(paste0(N_lakes," Lakes with Balanced Sampling Program"))+ theme(plot.title = element_text(hjust = 0.5))
p <- p + geom_hline(aes(yintercept=Truth), linetype=3, size=1.35, colour="darkblue")
p

ggsave("50_Lakes_Balanced_boxplot.png", p, scale = 1, width=11, height=8, units=c("in"), dpi = 500 )

#Unbalanced data
p <- ggplot(subset(d, Design %in% "Unbalanced"), aes(x=Nyears, y=MLE, fill=Model)) +
  geom_boxplot(outlier.alpha=0.2) + scale_fill_manual(values = c(rgb(0,0,0,0.3), rgb(1,1,1,0.3))) +
  xlab("Number of Years") + ylab("Parameter Value") + guides(aes(shape=NA)) +
  labs(fill = "")  + facet_wrap(~Parameter, scales=c("free_y"), labeller=labeller(Parameter = labels)) +
  scale_x_discrete(limits=levels(d$Nyears)) + theme(legend.key=element_blank())

p <- p + ggtitle(paste0(N_lakes," Lakes with Unbalanced Sampling Program"))+ theme(plot.title = element_text(hjust = 0.5))
p <- p + geom_hline(aes(yintercept=Truth), linetype=3, size=1.35, colour="darkblue")
p

ggsave("50_Lakes_Unbalanced_boxplot.png", p, scale = 1, width=11, height=8, units=c("in"), dpi = 500 )

#-----------------------------
#Plot the spatial model parameters only for balanced and unbalanced:
d <- reshape2::melt(sim_results[c("Spatial", "Nonspatial"),Designs,as.character(N_years),
                                as.character(nlakes), 1:Nsim,
                                c("ln_global_omega", "ln_global_linf", "ln_sd_linf",
                                  "global_tzero", "ln_sd_tzero", "ln_cv", "b_j_omega1",
                                  "b_j_omega2", "b_j_omega3", "ln_kappa", "rho", "sigmaO"), "50%"],
                    varnames=c("Model", "Design", "Nyears", "Replicate",
                               "Parameter"), value.name="MLE") %>%
  mutate(Truth=truth[Parameter]) %>%
  mutate(MLE=ifelse(grepl("ln", Parameter), exp(MLE), MLE)) %>%
  mutate(Truth=ifelse(grepl("ln", Parameter), exp(Truth), Truth)) %>%
  mutate(Parameter=stringr::str_remove(Parameter, "ln_"))

d$Nyears <- as.factor(d$Nyears)

labels <- c(global_omega = "Omega Intercept", cv = "Coefficient of Variation",
            global_linf= "Linf Intercept", b_j_omega1="Beta 1", b_j_omega2="Beta 2",
            b_j_omega3="Beta 3", sd_linf="Linf StDev", global_tzero="T0",
            sd_tzero="T0 StDev", kappa="Spatial Kappa", rho="Rho", sigmaO="StDev_O")

#Balanced data
p <- ggplot(subset(d, Design %in% "Balanced" & Model %in% "Spatial"), aes(x=Nyears, y=MLE)) +
  geom_violin(fill=rgb(0,0,0,0.3)) +
  xlab("Number of Years") + ylab("Parameter Value") + guides(aes(shape=NA)) +
  labs(fill = "") + facet_wrap(~Parameter, scales=c("free_y"), labeller=labeller(Parameter = labels)) +
  scale_x_discrete(limits=levels(d$Nyears)) + theme(legend.key=element_blank())

p <- p + ggtitle(paste0(N_lakes," Lakes with Balanced Sampling Program"))+ theme(plot.title = element_text(hjust = 0.5))
p <- p + geom_hline(aes(yintercept=Truth), linetype=3, size=1.35, colour="darkblue")
p <- p + geom_jitter(width=0.06)
p

ggsave("50_Lakes_Balanced_violin_st.png", p, scale = 1, width=11, height=8, units=c("in"), dpi = 600 )

#Unbalanced data
p <- ggplot(subset(d, Design %in% "Unbalanced" & Model %in% "Spatial"), aes(x=Nyears, y=MLE)) +
  geom_violin(fill=rgb(0,0,0,0.3)) +
  xlab("Number of Years") + ylab("Parameter Value") + guides(aes(shape=NA)) +
  labs(fill = "") + facet_wrap(~Parameter, scales=c("free_y"), labeller=labeller(Parameter = labels)) +
  scale_x_discrete(limits=levels(d$Nyears)) + theme(legend.key=element_blank())

p <- p + ggtitle(paste0(N_lakes," Lakes with Unbalanced Sampling Program"))+ theme(plot.title = element_text(hjust = 0.5))
p <- p + geom_hline(aes(yintercept=Truth), linetype=3, size=1.35, colour="darkblue")
p <- p + geom_jitter(width=0.06)
p

ggsave("50_Lakes_Unbalanced_violin_st.png", p, scale = 1, width=11, height=8, units=c("in"), dpi = 600 )

#----------------------------------------------
#Does the model recover the spatial-temporal field?
#Plot it for a single replicate fit w/ 20 years of data:

#Truth:
Sim_List = sim_rep[["sim_data"]]
DF = Sim_List[["DF"]]

d <- reshape2::melt(Sim_List[["eps_st"]]) %>%
       dplyr::mutate(x = rep(Sim_List[["Loc"]][,"x"], nyears),
       y = rep(Sim_List[["Loc"]][,"y"], nyears))

p <- ggplot(d, aes(x, y, col = value)) + geom_point(size=1) +
  facet_wrap(~L1) + theme(aspect.ratio = 1) + scale_color_gradient2()

#Estimated:
DF$eps <- sim_rep$opt_spatial$report$eps_i

d <-   reshape2::melt(DF$eps) %>%
       mutate(WBID=DF$Lake,
              x = DF$Longitude,
              y = DF$Latitude,
              year = DF$Year)%>%
              group_by(WBID, year)

#plot eps_i through time x space
p1 <- group_by(d, year) %>%
     # mutate(eps = value - mean(value)) %>%
     mutate(eps = value) %>%
     ggplot(aes(x, y, col =eps )) + geom_point(size=1) +
     facet_wrap(~year) + theme(aspect.ratio = 1) +
     xlab("x") + ylab("y") +
     scale_color_gradient2()
p1

p2 <- ggpubr::ggarrange(p, p1, labels = c("Truth", "Estimated"))

ggsave("RandomField.png", p2)
