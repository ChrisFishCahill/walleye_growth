#TODO clean up code, make more sexy figures
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

#------------------------
#Read in the data
#------------------------

data <- readRDS("C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/analysis/vB_analysis_august_2019_cahill.rds")
data <- as.data.frame(data)

sum(data$FL < 4.0)
data <- data[which(data$FL > 4.0),]
data <- data[which(data$FL  < 80), ]
nrow(data)

setwd("C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/executables")
VersionSpatial = "spdeXAR1_v3"
compile(paste0(VersionSpatial, ".cpp"))
dyn.load( dynlib(VersionSpatial) )


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

#Create the inputs for spatial model
spde = inla.spde2.matern( mesh )
spdeMatrices = spde$param.inla[c("M0","M1","M2")]
random_spatial = c("eps_omega_st", "eps_linf", "eps_t0")

#Run best nonspatial/spatial models with REML--LogNormal model (best as judged via PredNLL)
CTL <- 2 #lognormal
Partition_i <- rep(0, nrow(data)) #fit to all data

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
                          "logit_rho" = 2.9 )

obj_spatial = MakeADFun(data=data_spatial, parameters=parameters_spatial,
                        random=random_spatial, hessian=FALSE, DLL=VersionSpatial)

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
Range

#----------
#Make the flying spaghetti monster (seizure) plot:
#----------

t0       = ParHat$global_tzero
global_linf = exp(ParHat$ln_global_linf)
b_sex = exp(ParHat$ln_b_sex)
b_j_omega = ParHat$b_j_omega
eps_t0 = ParHat$eps_t0
eps_omega_st_i = rep$eps_i
eps_linf = ParHat$eps_linf
global_omega = exp(ParHat$ln_global_omega)
CV = exp(ParHat$ln_cv)

data$lpred <- rep$length_pred

data$pchCode = ifelse(data$SexCode==0, 16, 1)
Age_Seq <- 0:26
data$eps = rep$eps_i
data$colCode = ifelse(data$SexCode==0, "black", "white")

t_col <- function(color, percent = 50, name = NULL) {
  #	  color = color name
  #	percent = % transparency
  #	   name = an optional name for the color
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100,
               names = name)
  ## Save the color
  invisible(t.col)
}

mycol1 <- t_col("steelblue", perc = 15, name = "lt.blue")
mycol2 <- t_col("Orange", perc = 50, name = "lt.orange")

# png("C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/plots/seizure_plot.png",
#     width=8,height=6,units="in",res=1200)

par(mar = c(4.1, 4.1, 2.1, 2.1))

plot(data$FL~data$Age, lty=2, col="white", ylim=c(0,80), xlim=c(-.5, 26.5),
     xlab="Age (Years)", ylab="Total Length (cm)", lwd=2.5,
     las=1, cex.lab=1.25, cex.axis=1.15, bty="l", yaxs="i", xaxs="i")
points(data$FL~jitter(data$Age, factor=1.5), pch=21, bg=data$colCode, cex=0.5)

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
    Age_Seq = min(sub.sub.dat$Age):max(sub.sub.dat$Age)
    lpred_m <- lpred_f <- NA
    for(a in 1:length(Age_Seq)){
      lpred_m[a] = linf[1]*(1-exp(-(omega/linf[1]) * (Age_Seq[a] - t0 )))
      lpred_f[a] = linf[2]*(1-exp(-(omega/linf[2]) * (Age_Seq[a] - t0 )))
    }
    lines(lpred_m~Age_Seq, type="l", lwd=0.5, col=mycol1)
    lines(lpred_f~Age_Seq, type="l", lwd=0.5, col=mycol2)
  }}

Ages <- 0:max(data$Age)
Lpreds = global_linf*(1-exp(-(global_omega/global_linf) * (Ages - t0 )))
lines(Ages, Lpreds, lty=1, lwd=4, col="steelblue")

Lpreds = (global_linf + b_sex)*(1-exp(-(global_omega/(global_linf+b_sex)) * (Ages - t0 )))
lines(Ages, Lpreds, lty=1, lwd=4, col="orange")

# dev.off()

#--------------
#Make a spatial-temporal plot
#--------------

data$omega <- rep$omega_i

d <- mutate(data, WBID=data$WBID, x = data$Long_c,
            y = data$Lat_c, year = data$Year + 1999,
            eps = data$eps,
            mu = data$omega,
            se=data$se_omega) %>% group_by(WBID, year) %>%
     filter(year < 2018)
nrow(d)

#----------------
colfunc <- colorRampPalette(c("darkblue","white", "darkorange2"))
colfunc(25)
plot(rep(1,25),col=colfunc(25),pch=19,cex=3)


can1<-raster::getData('GADM', country="CAN", level=1)
alta = can1[can1$NAME_1 %in% "Alberta", ]

plot(alta)
points(data$Long_c, data$Lat_c)

alta = can1[can1$NAME_1 %in% "Alberta", ]
str(alta)
class(alta)

alta <- spTransform(alta,
                       CRS("+proj=longlat +datum=WGS84"))
alta.fort = fortify(alta)

names(alta.fort)[1] <- "Long_c"
names(alta.fort)[2] <- "Lat_c"

p <- ggplot(NULL) + theme_classic( ) +
     geom_polygon(colour="black", fill="white", data=alta.fort, aes(x=Long_c, y=Lat_c, group=id)) +
     geom_point( pch=21, size=2.0, data=d, aes(Long_c, Lat_c, fill=mu) ) +
     scale_x_continuous(breaks=c(-120, -115, -110)) +
     scale_y_continuous(breaks=c(49,52,56,60)) +
     scale_fill_gradient2(low="steelblue", high="darkorange1",
                          midpoint= 15.5,
     name=bquote(atop(Growth~Rate~bold((omega)),~cm%.%year^{-1})))

p <- p + facet_wrap(~year, nrow=3) +
         ylab("Latitude") + xlab("Longitude") +
         theme(axis.title=element_text(size=15),
               strip.background = element_blank(),
               strip.text.x = element_blank(),
               panel.spacing.x=unit(1, "lines"),
               panel.spacing.y=unit(0.5, "lines"))
p <- p + geom_text(data = d,
                   mapping = aes(x = -112, y = 59.5, label = year,
                                 group=year))

p <- p + ggalt::coord_proj(
  paste0(CRS("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs ")))

ggsave("C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/plots/growth_predictions.tiff",
       p, dpi=1200, width=11, height=8, units=c("in"))

#------------
#Now, do it for four years of data
#-------------

d <- mutate(data, WBID=data$WBID, x = data$Long_c,
            y = data$Lat_c, year = data$Year + 1999,
            eps = data$eps,
            mu = data$omega,
            se=data$se_omega) %>% group_by(WBID, year) %>%
  filter(year < 2018)

target = c(2008, 2011, 2014, 2017)
d <- d %>% filter(year %in% target)

p <- ggplot(NULL) + theme_classic( ) +
  geom_polygon(colour="black", fill="white", data=alta.fort, aes(x=Long_c, y=Lat_c, group=id)) +
  geom_point( pch=21, size=2.0, data=d, aes(Long_c, Lat_c, fill=mu) ) +
  scale_x_continuous(breaks=c(-120, -115, -110)) +
  scale_y_continuous(breaks=c(49,52,56,60)) +
  scale_fill_gradient2(low="steelblue", high="darkorange1",
                       midpoint= 15.5,
                       name=bquote(atop(Growth~Rate~bold((omega)),~cm%.%year^{-1})))

p <- p + facet_wrap(~year, nrow=1) +
  ylab("Latitude") + xlab("Longitude") +
  theme(axis.title=element_text(size=15),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing.x=unit(1, "lines"),
        panel.spacing.y=unit(0.5, "lines"))
p <- p + geom_text(data = d,
                   mapping = aes(x = -112, y = 59.5, label = year,
                                 group=year))

p <- p + ggalt::coord_proj(
  paste0(CRS("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs ")))

# ggsave("C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/plots/growth_predictions_reduced.png",
#         p, dpi=1200, width=11, height=4, units=c("in"))

#------------------------------
#Plot the spatial correlation
#------------------------------

Kappa = exp(ParHat$ln_kappa)

D <- as.matrix(dist(loc_xy))

dis.vec <- seq(0, max(D), length = 1000)
Cor.M <- (Kappa * dis.vec) * besselK(Kappa * dis.vec, 1) #matern correlation
Cor.M[1] <- 1

png( file="C:/Users/Chris Cahill/Documents/GitHub/walleye_growth/plots/SpatialCorrelation.tiff",
     width=11, height=5, res=1200, units="in")
par(mfrow=c(1,2),  mar=c(3,3,2,1), mgp=c(2,0.5,0), tck=-0.02)

plot(x = dis.vec, y = Cor.M,
     type = "l", cex.lab = 1.25, cex.axis=1.15,
     xlab = "Distance (km)", yaxs="i", xaxs="i",
     ylab = "Correlation", xlim = c(0, 110), bty="l",
     las=1)

hist(dist(loc_xy), xlab="Distance between Lakes (km)", main="", cex.lab = 1.25,
     breaks=40, bty="o", cex.axis=1.15, yaxs="i", las=1)
abline(v=Range, lty=3, lwd=6, col="Steelblue")
dev.off()
