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
#Cahill ^(;,;)^ March 2019
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

png(file="plots/Mesh.png",width=9.50,height=7.00,units="in",res=600)
plot(mesh)
points(loc_xy, col="Steelblue", pch=1)
dev.off()

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

#Create the inputs for TMB and run the model

Version = "spdeXAR1_v3"
compile(paste0(Version, ".cpp"))
dyn.load( dynlib(Version) )

spde = inla.spde2.matern( mesh )

spdeMatrices = spde$param.inla[c("M0","M1","M2")]

Data = list("Nobs"      = nrow(data), "length_i"=data$TL, "age_i" = data$Age,
            "lake_i"    = data$Lake - 1, "sex_i" = data$SexCode,
            "X_ij_omega"= model.matrix(~ -1 + data$wallEffDen.Std + #.Std --> already standardized
                                              data$compEffDen.Std +
                                              data$GDD.Std  +
                                              data$wallEffDen.Std:data$compEffDen.Std),
            "Nlakes" = length(unique(data$Lake)),
            "spdeMatrices" = spdeMatrices,
            "s_i" = data$Lake-1,
            "t_i" = data$Year-1,
            "n_t" = max(data$Year),
            "t_prev_i" = data$Year-2,
            "CTL" = 2  #CTL==1 --> Normal, 2 --> Lognormal, 3 --> Gamma
)


Parameters = list("ln_global_omega"  = log(13.28954),
                  "ln_global_linf"   = log(55),
                  "ln_sd_linf"       = log(7.275452),
		  "global_tzero"     = -1.260783461,
                  "ln_sd_tzero"      = log(0.538755),
                  "ln_b_sex"         = log(4.760871),
                  "b_j_omega"        = rep(0, ncol(Data$X_ij_omega)),
                  "eps_omega_st"     = matrix(0,  nrow=mesh$n,ncol=Data$n_t ),
                  "eps_linf"         = rep(0,  length(unique(data$Lake))),
		  "eps_t0"           = rep(0, length(unique(data$Lake))),
		  "ln_cv"            = -2.5043055,
		  "ln_kappa"         = -2,
		  "ln_tau_O"         = -5,
		  "rho_unscaled"     = 2 * plogis(0.6) - 1 # --> 2 * plogis(your_rho_guess) - 1
)

Random = c("eps_omega_st", "eps_linf", "eps_t0")
Use_REML = FALSE
if( Use_REML==TRUE ) Random = union( Random, c("ln_global_omega", "ln_global_linf",
                                    "global_tzero","ln_b_sex", "b_j_omega") )

Obj <- MakeADFun(data=Data, parameters=Parameters, random=Random, hessian=FALSE, map=NULL)

Obj$fn( Obj$par )
Obj$gr( Obj$par )

#Runs in 4.3 minutes on my laptop
Opt  = TMBhelper::Optimize(obj=Obj,
			   control=list(eval.max=1000, iter.max=1000),
			   getsd=T, newtonsteps=1, bias.correct=F)

Opt
#estimates/standard errors appear reasonable

SD = sdreport( Obj )
final_gradient = Obj$gr( Opt$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

#Extract the intercept and SE
ParHat = as.list( Opt$SD, "Estimate" )
SEHat  = as.list( Opt$SD, "Std. Error" )

rep <- Obj$report()

Range <- Opt$SD["value"]$value[which(names(Opt$SD["value"]$value)=="Range")]
rho <- Opt$SD["value"]$value[which(names(Opt$SD["value"]$value)=="rho")]
#Opt$SD$sd
Range
rho

hist(dist(loc_xy), xlab="Distance between Lakes (km)", main="")
abline(v=Range, lty=3, lwd=3, col="Steelblue")
# Range (14 km) is pretty small given most sites are not within that distance:
# perhaps not that surprising...

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

#Generate/plot the predictive distribution conditional on the best parameter estimates a la Gelman:

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
data$eps <- rep$eps_i

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
	if(Data$CTL==1) {pred_dist[i, 1:Nsim] = rnorm(Nsim, lpred, sigma)}
	if(Data$CTL==2) {pred_dist[i, 1:Nsim] = rlnorm(Nsim, 	log(lpred) - sigma^2/2, sigma)}
	if(Data$CTL==3) {pred_dist[i, 1:Nsim] = rgamma(Nsim, shape=1/CV^2, scale=lpred*CV^2)}

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
p <- ggplot(data = data, aes(y = TL, x = Age))
p <- p + xlab("Age") + ylab("Total Length (cm)")
p <- p + geom_point(shape = 16, size = 1.0, position="jitter")
p <- p + geom_ribbon(aes(x = Age, ymax = fit.975, ymin = fit.025), fill = grey(0.5), alpha = 0.5)
p <- p + geom_ribbon(aes(x = Age, ymax = fit.25, ymin = fit.75), fill = grey(1), alpha = 0.8)
p <- p + geom_line(aes(x = Age, y = fit.5), linetype = "dashed",  col="black", size=0.75)
p <- p + scale_x_continuous(limits = c(-0.5, 27), breaks=c(0,6,11,16,21,26)) +
	 scale_y_continuous(limits = c(0, 86), breaks=c(0,20,40,60,80))
p <- p + ggthemes::theme_tufte() +
	       theme(axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
		     axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid'),
	             axis.text=element_text(size=14, colour="black"),
	       	     axis.title=element_text(size=20,face="bold") )
p

ggsave("plots/PredictiveDistributionFit.png", p, dpi=600, width=9, height=8, units=c("in"))

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

#quantify coverage:

quantvec <- NULL
for(i in 0:26){
	qs <- quantile(pred_dist[which(data$Age == i)], c(0.025, 0.975))
	sub.data <- data[which(data$Age == i), "TL"]
	prop <- sum(sub.data >= qs[1] & sub.data <= qs[2])/length(sub.data)
	#hist(sub.data, main=paste(paste0("Age", i), round(prop,2), sep="  "),
	#		 xlim=qs+c(-150,150), xlab="", ylab="")
	#abline(v=qs, col="Steelblue", lwd=3, lty=3)
	quantvec[i+1] <- prop
}

print(paste( round(100*sum(quantvec) / length(quantvec), 2), "% coverage within 95% C.I.s", sep=""))
#Okay-ish..
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

#lake-year-specific mle predictions
#circles = girls, squares = boys
Age_Seq <- 0:26
pdf("plots/lake_year_mles.pdf")
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
			b_j_omega[1]*sub.sub.dat$wallEffDen.Std[j] +
			b_j_omega[2]*sub.sub.dat$compEffDen.Std[j] +
			b_j_omega[3]*sub.sub.dat$GDD.Std[j] +
			b_j_omega[4]*sub.sub.dat$wallEffDen.Std[j]*sub.sub.dat$compEffDen.Std[j] +
			eps_omega_st_i[which(data$WBID==i & data$Year==j)[1]]

		linf = global_linf +
			b_sex*sub.sub.dat$SexCode[j] +
			eps_linf[Lake]

		lpred_m <- lpred_f <- NA
		for(a in 1:length(Age_Seq)){
			lpred_m[a] = linf*(1-exp(-(omega/linf) * (Age_Seq[a] - t0 )))
			lpred_f[a] = (b_sex + linf)*(1-exp(-(omega/(linf+b_sex)) * (Age_Seq[a] - t0 )))
		}

		plot(lpred_f~Age_Seq, type="l", col="black", ylim=c(0,85), main=paste(Name, j, sep= " " ),
		     xlab="Age (Years)", ylab="Total Length (cm)", cex.main=1, lwd=1.5)
		points(lpred_m~Age_Seq, type="l", col="steelblue", lwd=1.5)
		points(sub.sub.dat$TL~sub.sub.dat$Age, pch=sub.sub.dat$SexCode)
}}
dev.off()

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

#General thoughts on the patterns in the lake-year-specific mle.pdf:
#The model is fitting, but there are situations where TL-Age data in a given year/lake are
#more or less linear.  This sometimes results in ugly patterns of under/over prediction, which makes sense
#because the von B will be a garbage fit to a linear pattern. Also, given we are fitting 252 lake-year
#combinations of von bertalanffys this is probably to be ocassionally expected.

#This is also occuring because I am assuming each lake has a single Linfinity and t0--this was necessary to
#maintain estimates in the realm of reality given data limitations.  So, the curvature/growth rate (omega)
#and its spatial-temporal random effects can only compensate to a certain degree. I could put a spatial-temporal
#field on linfinity as well, but would likely need an informative prior to keep it in the realm of reality.

#I would be *very* careful trying to use this model to say something like "the growth curve for walleye in
#lake x in year y is this--> sometimes this is fine, but other times this would catastrophically fail.

#Nonetheless, it is reassuring that the general estimates of the betas are nearly identical to the spatial-only
#version of this model.  The results between these spatial/spatial-temporal models are biolgically the same.I'm
#not entirely convinced that this model is inappropriate for its original purpose, which was to
#determine the effect of a few covariates on omega, but I may be over-valuing a model that I built...

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

#Plot stuff related to the spatial temporal field

data$E2 <- (data$TL - rep$length_pred)^2

d <-   reshape2::melt(data$E2)       %>%
       mutate(WBID=data$WBID,
              x = data$X_TTM_c/1000,
              y = data$Y_TTM_c/1000,
              year = data$Year,
              eps = data$eps )	     %>%
              group_by(WBID, year)


#Plot eps_i (st raneffs of the lakes) through time x Space:
p <- ggplot(d, aes(x, y, col =eps )) + geom_point() +
     facet_wrap(~year) +
     xlab("Easting (km)") + ylab("Northing (km)") +
     scale_colour_viridis_c() #scale_color_gradient2()
p

ggsave("plots/omega_st_facet_year.png", p, dpi=600, width=11, height=8, units=c("in"))

#plot mean-centered eps_i through time x space
p <- group_by(d, year) %>%
     mutate(eps = eps - mean(eps)) %>%
     ggplot(aes(x, y, col =eps )) + geom_point() +
     facet_wrap(~year) +
     xlab("Easting (km)") + ylab("Northing (km)") +
     scale_color_gradient2(limits = c(-0.01, 0.01))

ggsave("plots/omega_st_mean_centered.png", p, dpi=600, width=11, height=8, units=c("in"))

# So, the space random effects change through time, but are rather consistent in
# a given year.  Interesting--perhaps because data don't support the more complex
# space time model? Nonetheless, much different in terms of coefficients between spatial &
# spatial-temporal versions of this model lead to identical ecological inference.
# Nice to know that the ecological inference holds up with the AR-1 ST field

#Now let's look at the SE of the field...
d <- reshape2::melt(SEHat[["eps_omega_st"]]) %>%
     mutate(x = rep(mesh$loc[,1], ncol(SEHat[["eps_omega_st"]])),
            y = rep(mesh$loc[,2], ncol(SEHat[["eps_omega_st"]])),
            year = col(SEHat[["eps_omega_st"]]))

#Plot the SEs through time-space
p <- ggplot(d, aes(x, y, col = value)) +
     geom_point() +
     facet_wrap(~Var2) +
     scale_colour_viridis_c()
p

ggsave("plots/omega_st_SE_facet_year.png", p, dpi=600, width=11, height=8, units=c("in"))

#Reconstruct the spatial temporal field for omega using inla.mesh.projector
#Same as above but visualized differently via inla mesh projector

st_field_omega <- matrix(ParHat$eps_omega_st, ncol=length(unique(data$Year)), nrow=mesh$n)
for(i in 2:ncol(st_field_omega)){
  st_field_omega[ ,i] <- rho*st_field_omega[ ,i-1] + st_field_omega[,i]
}

png( file="plots/omega_st_facet_year_base.png", width=11, height=8, res=300, units="in")
par(mfrow=c(3,6), mar=c(4.1, 4, 3.5, 3.5))
for(dtPlot in 1:(max(data$Year)))
{
  proj = inla.mesh.projector(mesh, dims=c(100, 100))

  XYUTM = mesh$loc[,1:2]
  latentFieldML = st_field_omega[,1]

  fields::image.plot(proj$x,proj$y, inla.mesh.project(proj, latentFieldML),col =  viridis(35),
                     xlab = 'Easting Km', ylab = 'Northing Km',
                     zlim = c(-3.5, 2.75),
                     main = paste("Year ", dtPlot ,sep = ""),
                     cex.lab = 1.1,cex.axis = 1.1, cex.main=1, cex.sub= 1.1)

  sub.data <- data[which(data$Year == dtPlot),]
  loc_xy_y <- unique(sub.data[ ,c("X_TTM_c","Y_TTM_c") ] )
  loc_xy_y <- as.matrix(loc_xy_y)/1000

  points(loc_xy_y, col="black", pch=16, cex=0.5)
  #contour(proj$x, proj$y,inla.mesh.project(proj, latentFieldML),nlevels = 2 ,add = T,labcex  = 1,cex = 0.5)
  #lines(borders, lwd=3)
}
dev.off()

#And now for the Standard Errors of the spatial-temporal effects
SEHat$eps_omega_st

st_field_omega_SE <- matrix(SEHat$eps_omega_st, ncol=length(unique(data$Year)), nrow=mesh$n)

png( file="plots/omega_st_SE_facet_year_base.png", width=11, height=8, res=300, units="in")
par(mfrow=c(3,6), mar=c(4.1, 4, 3.5, 3.5))
for(dtPlot in 1:(max(data$Year)))
{
  proj = inla.mesh.projector(mesh, dims=c(100, 100))

  XYUTM = mesh$loc[,1:2]
  latentFieldML = st_field_omega_SE[,i]

  fields::image.plot(proj$x,proj$y, inla.mesh.project(proj, latentFieldML),col =  viridis(35),
                     xlab = 'Easting Km', ylab = 'Northing Km',
                     zlim = c(0.15, 2.1),
                     main = paste("Year ", dtPlot ,sep = ""),
                     cex.lab = 1.1,cex.axis = 1.1, cex.main=1, cex.sub= 1.1)

  sub.data <- data[which(data$Year == dtPlot),]
  loc_xy_y <- unique(sub.data[ ,c("X_TTM_c","Y_TTM_c") ] )
  loc_xy_y <- as.matrix(loc_xy_y)/1000

  points(loc_xy_y, col="black", pch=16, cex=0.5)
  #contour(proj$x, proj$y,inla.mesh.project(proj, latentFieldML),nlevels = 2 ,add = T,labcex  = 1,cex = 0.5)
  #lines(borders, lwd=3)
}
dev.off()

# SE's higher in locations where there are lakes, unsurprisingly
# Not changing much through time, similar to the MAP estimates of the random effect (last plot)
# Maybe index this to extract only lakes similar to the plot above



#Plot RMSE through space-time --> for what they are worth

d <-   reshape2::melt(data$E2)       %>%
       mutate(WBID=data$WBID,
              x = data$X_TTM_c/1000,
              y = data$Y_TTM_c/1000,
              year = data$Year,
              eps = data$eps )	     %>%
         group_by(WBID, year)        %>%
         mutate(RMSE = sqrt(mean(value)), MSE=mean(value))

ggplot(d, aes(x, y, col = RMSE)) + geom_point() +
  scale_colour_viridis_c() + #	scale_colour_gradient2()
  xlab("Easting (km)") + ylab("Northing (km)") +
  facet_wrap(~year)

#Plot MSE through space-time
ggplot(d, aes(x, y, col = MSE)) + geom_point() +
  scale_colour_viridis_c() + #	scale_colour_gradient2()
  xlab("Easting (km)") + ylab("Northing (km)") +
  facet_wrap(~year)

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

#Let's do some more residual checks:

#Some ugly residual plots for what they are worth:
length_pred <- rep$length_pred
data$E1 <- (data$TL - length_pred)

#Residuals by year:
xyplot(data$E1 ~ rep$length_pred | data$Year,
			 panel=function(x, y){
			 panel.xyplot(x, y)
			 } )
#Okay-ish.

#Residuals by lake:
xyplot(data$E1 ~ rep$length_pred | data$Name,
			 panel=function(x, y){
			 panel.xyplot(x, y)
			 	#panel.lmline(x, y, lty = 2)
			 } )
#not the prettiest here -- definitely some patterns.

# Homogeneity
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x=length_pred, y=data$E1, ylab="Residual", xlab="Length")
abline(h = 0, v = 1)
#I've seen worse, but not the best..

# Normality
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
hist(data$E1, main="", xlab="Residual")
#meh

# Independence due to model misfit
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = data$wallEffDen, y = data$E1, xlab="Intraspecific Density", ylab="Residual")
abline(h = 0)
#Doesn't appear to be any nonlinear tomfoolery afoot for walleye density

# Independence due to model misfit
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = data$compEffDen, y = data$E1, xlab="Interspecific Density", ylab="Residual")
abline(h = 0)
#Okay

# Independence due to model misfit
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = data$GDD, y = data$E1, xlab="Growing Degree Days", ylab="Residual")
abline(h = 0)

# Independence due to model misfit
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = data$Year, y = data$E1, xlab="Year", ylab="Residual")
abline(h = 0)

#Residuals are getting a little broader through time, but
#also more samples collected through time...

#residuals by sex
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
boxplot(data$E1~data$SexCode, xlab="Sex")
abline(h = 0)

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

#This next section ignores the interaction between walleye and interspecific density & plots
#density vs. omega -- just to get a feel for what is going on

#Lake specific omega vs. effective density:
omegas <- global_omega + b_j_omega[1]*data$wallEffDen.Std +
	                 b_j_omega[2]*data$compEffDen.Std +
	                 b_j_omega[3]*data$GDD.Std +
	                 b_j_omega[4]*data$wallEffDen.Std*data$compEffDen.Std +
                         eps_omega_st_i
data$omegas <- omegas

#Omegas vs. density
p <- ggplot(data, aes(data$wallEffDen.Std, omegas)) + geom_point() +
	    ylab("Growth Rate cm/year (Omega)") + xlab("Effective Walleye Density (IGNORES INTERACTION)")
p

#Omegas vs. density by year
p1 <- ggplot(data, aes(data$wallEffDen.Std, omegas)) + geom_point() +
	    ylab("Growth Rate cm/year (Omega)") + xlab("Effective Walleye Density (IGNORES INTERACTION)") +
	    facet_wrap(~Year) +
      geom_smooth(data = data, method="lm", aes(x = wallEffDen.Std, y =  omegas))
p1
#Omegas generally decline with increased conspecific density

#Omegas vs. GDD
p2 <- ggplot(data, aes(data$GDD.Std, omegas)) + geom_point() +
	ylab("Growth Rate cm/year (Omega)") + xlab("GDD")
p2

#Omegas vs. GDD by year
p3 <- ggplot(data, aes(data$GDD.Std, omegas)) + geom_point() +
	    ylab("Growth Rate cm/year (Omega)") + xlab("GDD") +
	    geom_smooth(data = data, method="lm", aes(x = GDD.Std, y =  omegas)) +
	    facet_wrap(~Year)
p3

#Clear relationship between growth rate early in life and GDD

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

#Visualize the interaction between intra- and inter-specific density
#with *supremely* ugly code:

grid.lines <- 35
#Create an x and y sequence based on the range of your variables
x.pred <- seq(min(data$wallEffDen.Std), max(data$wallEffDen.Std), length.out = grid.lines)
y.pred <- seq(min(data$compEffDen.Std), max(data$compEffDen.Std), length.out = grid.lines)
xy <- expand.grid( x = x.pred, y = y.pred)

#Create the Z-dimension predictions holding GDD at median level:
z.pred <- NA
for(i in 1:nrow(xy)){
 z.pred[i] <- global_omega +                      #global growth rate for critters
 	            b_j_omega[1]*xy$x[i] +              #walleye density
	            b_j_omega[2]*xy$y[i] +              #interspecific density
 	            b_j_omega[3]*median(data$GDD.Std) + #hold GDD at median level
	            b_j_omega[4]*prod(xy[i,])           #interaction term
}

#Determine whether a given omega is above/below the plane:
data$pch <- NA
for(i in 1:nrow(data)){
 #Find the x y coordinates for obs i
 x <-  data$wallEffDen.Std[i]
 y <-  data$compEffDen.Std[i]

 #find the index for zpreds (omega)
 x_pred <- unique(xy[which(abs(xy[ ,1]-x)==min(abs(xy[,1]-x))), 1])
 y_pred <- unique(xy[which(abs(xy[ ,2]-y)==min(abs(xy[,2]-y))), 2])
 idx <- which(xy[ ,1] == x_pred & xy[,2] == y_pred)

 #declare point character for below
 ifelse(data$omegas[i] > z.pred[idx], data$pch[i] <- 1, data$pch[i] <- 16)
}

z.pred <- matrix(z.pred, nrow=grid.lines, ncol=grid.lines)

#make a pretty plot:
png(file="plots/InteractionPlot.png",width=9.50,height=7.00,units="in",res=600)
par(mfrow=c(1,1))
scatter3D(x=data$wallEffDen.Std, y=data$compEffDen.Std, z = data$omegas, type="n",
	  theta = 25, phi=17, ticktype = "detailed", type="h",pch = data$pch, cex = 1, col=viridis(50),
	  xlab = "Intraspecific Density", ylab = "Interspecific Density", zlab = "Omega (Growth Rate)",
          surf = list(x = x.pred, y = y.pred, z = z.pred, facets = NA), main = "")

dev.off()

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

#Visualize some of the key parameters of interest:
#beta1 -- > intraspecific effective density effect
#beta2 -- > interspecific  ''        ''      ''
#beta3 -- > Growing Degree Days
#beta4 -- > Density interaction term
#

MLE <- Opt$SD$par.fixed
SE  <- SEHat[which(names(SEHat)==names(MLE))]

d <- as.data.frame(cbind(MLE, SE=unlist(SE)))
#Declare the lower and upper 95% CIs (in case you believe in that sort of thing)
d$low <- d$MLE - 1.96*d$SE
d$up  <- d$MLE + 1.96*d$SE
head(d)
d$Parameter <- rownames(d)
rownames(d) <- NULL

d1 <- as.data.frame(cbind(MLE=Opt$SD$value, SE=Opt$SD$sd))
d1$Parameter <- rownames(d1)
rownames(d1) <- NULL
d1$low <- d1$MLE - 1.96*d1$SE
d1$up  <- d1$MLE + 1.96*d1$SE
d1

d <- as.data.frame(rbind(d, d1))

p <- ggplot()
p <- p + geom_point(data = d,
                    aes(x = Parameter, y = MLE) )

p <- p + geom_errorbar(data = d,
                       aes(x = Parameter,
                           ymax = up,
                           ymin = low),
                       width=0.2)

p <- p + xlab("Parameters") + ylab("MLE and 95% CI")
p <- p + theme(text = element_text(size=10))
p <- p + geom_hline(yintercept = 0, lty = 2)
p

d

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

#Look at the covariance matrix:

scaledCovar <- cov2cor(Opt$SD$cov.fixed) #scales the covariance matrix
colnames(scaledCovar) <- rownames(scaledCovar) <- rownames(Opt$SD$cov.fixed)
print(round(scaledCovar, 2))

#ln_kappa and ln_tau_O are crazy correlated--probably expected

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

#Profile parameters?
Profile = FALSE
if(Profile == TRUE){
	#profile all the parameters:
	profile1  <- tmbprofile(Obj, "ln_tau_O", parm.range = c(-10, 10) )
	profile2  <- tmbprofile(Obj, "ln_kappa", parm.range = c(-10, 10) )
	profile3  <- tmbprofile(Obj, "ln_global_omega", parm.range = c(-10, 10) )
	profile4  <- tmbprofile(Obj, "ln_global_linf", parm.range = c(-10, 10) )
	profile5  <- tmbprofile(Obj, "ln_sd_linf", parm.range = c(-10, 10) )
	profile6  <- tmbprofile(Obj, "global_tzero", parm.range = c(-10, 10) )
	profile7  <- tmbprofile(Obj, "ln_sd_tzero", parm.range = c(-10, 10) )
	profile8  <- tmbprofile(Obj, "ln_b_sex", parm.range = c(-10, 10) )
	profile9  <- tmbprofile(Obj, "ln_cv", parm.range = c(-10, 10) )

	par(mfrow=3,3)

	plot(profile1)
	plot(profile2)
	plot(profile3)
	plot(profile4)
	plot(profile5)
	plot(profile6)
	plot(profile7)
	plot(profile8)
	plot(profile9)
}

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

#End End End




