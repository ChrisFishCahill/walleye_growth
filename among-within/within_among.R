#-------------------------------------------------------------------------------------------------------------------
#Nonspatial and Spatial-temporal von Bertalanffy with density dependent growth regression
#within-among parameterization on reduced dataset
#Cahill Oct 2019
#-------------------------------------------------------------------------------------------------------------------

library(dplyr)
library(TMB)
library(INLA)
library(ggplot2)

data <- readRDS("C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/analysis/vB_analysis_august_2019_cahill.rds")
data <- as.data.frame(data)

#Select the lakes with >= 3 years of data
data <- data %>% group_by(WBID) %>%
                    filter(n_distinct(Year) >= 3) %>% ungroup()

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

#Reorder the lakes or else tmb explodes
data <- within(data, Lake <- as.numeric(interaction(data$WBID, drop=TRUE, lex.order=F)))
data <- data[order(data$Lake),] #Order data on groups for ease

mapbox <- c(-120.5, 48.95, -109.5, 60.5)

glgmap   <- ggmap::get_stamenmap(bbox=mapbox, maptype="terrain", zoom=4, crop=T, force=T)

map <- ggmap::ggmap(glgmap) + geom_point(aes(Long_c, Lat_c), pch=21, size=2, data=data )

map

#----------
# spatial model:
#----------

setwd("C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/among-within")
VersionSpatial = "among_within"
compile(paste0(VersionSpatial, ".cpp"))
dyn.load( dynlib(VersionSpatial) )

loc_xy <- unique(data[ ,c("X_TTM_c","Y_TTM_c") ] )
loc_xy <- loc_xy/1000 #Put distance in kms

mesh = inla.mesh.2d(loc=loc_xy, max.edge=c(62,1000)) #lots of minutes
#mesh = inla.mesh.create( loc_xy, refine=TRUE, extend=-0.5, cutoff=0.01 ) #2.3 minutes
mesh$n

#Create the inputs for spatial model
spde = inla.spde2.matern( mesh )
spdeMatrices = spde$param.inla[c("M0","M1","M2")]
random_spatial = c("eps_omega_st", "eps_linf", "eps_t0")

#Lognormal Likelihood--spatial                                                                                                                                          "ln_b_sex", "b_j_omega") )
CTL = 2

data_spatial = list("Nobs" = nrow(data), "length_i"=data$TL, "age_i" = data$Age,
                    "lake_i" = data$Lake - 1, "sex_i" = data$SexCode,
                    "X_ij_omega"= model.matrix(~ -1 + data$lake_mean_dens.std + #among lake effect
                                                      data$compEffDen.Std +
                                                      data$GDD.Std  +
                                                      data$wallEffDen.Std:data$compEffDen.Std),
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
                          "rho" = 2.9) # 2 * plogis(3) - 1 = 0.9

random_spatial = c("eps_omega_st", "eps_linf", "eps_t0", "eps_slope")

Use_REML = T
if( Use_REML==TRUE ) random_spatial = union( random_spatial, c("ln_global_omega",
                                                               "ln_global_linf", "global_tzero",
                                                               "ln_b_sex", "b_j_omega", "mu_slope") )

obj_spatial = MakeADFun(data=data_spatial, parameters=parameters_spatial,
                        random=random_spatial,
                        hessian=FALSE, DLL=VersionSpatial)

opt_spatial  = TMBhelper::Optimize(obj=obj_spatial,
                                    control=list(eval.max=1000, iter.max=1000),
                                    getsd=T, newtonsteps=1, bias.correct=T)

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

slopes = data.frame(ParHat$eps_slope + SEHat$eps_slope%o%qnorm(c(0.025, 0.5, 0.975)))
colnames(slopes) <- c("Lower95", "MLE", "Upper95")
slopes$Name = unique(data$Name)
print(slopes)

slopes = slopes[with(slopes, order(MLE)),]
slopes$Name <- factor(slopes$Name, as.character(slopes$Name))


p <- ggplot(data = slopes,
                    aes(x = MLE,
                        xmax = Upper95,
                        xmin = Lower95,
                        y = Name)) +  geom_point() +
  geom_segment( aes(x = Lower95, xend = Upper95,
                    y = Name, yend=Name)) +
  xlab("Slope Estimate for Walleye Density Effect") + ylab("Lake") +
  theme(text = element_text(size = 10)) +
  theme(axis.title = element_text(size=15)) +
  theme(axis.line = element_line(colour = 'black', size = 1, linetype = 'solid'))

p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
p
# ggsave(file="C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/among-within/slopes.png",
#        height=10, width=9)
