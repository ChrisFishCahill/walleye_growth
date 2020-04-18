#---------------------
# Load packages
#NOTE: ggplot3.3.0 seems to be incompatible with ggalt0.4.0::coord_proj()
#NOTE: here I force install of ggplot3.2.1 and things behave
#NOTE: install.packages("https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.2.1.tar.gz", repo=NULL, type="source")
#Or, install this version of ggalt that Elio Campitelli figured out cleverly:
#devtools::install_github("eliocamp/ggalt@new-coord-proj")
library(ggplot2)
theme_set(theme_light())
library(TMB)
library(future)
library(tidyr)
library(INLA)
library(purrr)
library(furrr)
library(dplyr)
library(TMBhelper)
library(arm)
library(gtable)
library(grid)
library(gridExtra)
library(ggalt)
library(ggpubr)
library(ggforce)

# Read in ml and reml results
ml <- readRDS(file = "analysis2/ML_fits.rds")
reml <- readRDS(file = "analysis2/REML_fits.rds")

# Read in the data and make it idential to analysis:
data <- readRDS("analysis2/vB_analysis_august_2019_cahill.rds")
data[, c("X_TTM_c", "Y_TTM_c")] <- data[, c("X_TTM_c", "Y_TTM_c")] / 1000 # put distance in km
Loc <- unique(data[, c("X_TTM_c", "Y_TTM_c")])

Loc$Block <- NA
Loc[which(Loc$X_TTM_c < 400 & Loc$Y_TTM_c > 6400), "Block"] <- 1
Loc[which(Loc$X_TTM_c < 410 & Loc$Y_TTM_c < 6200 & is.na(Loc$Block)), "Block"] <- 2
Loc[which(Loc$X_TTM_c < 450 & Loc$Y_TTM_c < 6000 & is.na(Loc$Block)), "Block"] <- 2
Loc[which(Loc$X_TTM_c < 500 & Loc$Y_TTM_c > 6100 & is.na(Loc$Block)), "Block"] <- 3
Loc[which(Loc$X_TTM_c < 660 & Loc$Y_TTM_c > 6000 & is.na(Loc$Block)), "Block"] <- 4
Loc[which(Loc$X_TTM_c < 625 & Loc$Y_TTM_c < 6000 & is.na(Loc$Block)), "Block"] <- 5
Loc[which(Loc$Y_TTM_c < 5800 & is.na(Loc$Block)), "Block"] <- 6
Loc[which(Loc$X_TTM_c > 660 & Loc$Y_TTM_c > 5800 & is.na(Loc$Block)), "Block"] <- 7
data <- dplyr::left_join(data, Loc, by = c("X_TTM_c", "Y_TTM_c"))

data <- data %>%
  group_by(WBID) %>%
  summarize(lake_mean_dens = mean(wallEffDen)) %>%
  mutate(lake_centered_dens.std = arm::rescale(lake_mean_dens)) %>%
  dplyr::left_join(data, .x, by = "WBID")

data <- within(data, Lake <- as.numeric(interaction(data$WBID, drop = TRUE, lex.order = F)))
data <- data[order(data$Lake), ]

mesh <- inla.mesh.create(Loc[, c("X_TTM_c", "Y_TTM_c")], refine = TRUE, extend = -0.5, cutoff = 0.01) # faster mesh
#---------------------
# make Table 1 (misery):

cv_h_block <- readRDS("analysis2/cv_h_block.rds")

unique(cv_h_block$convergence)
MyTable <- cv_h_block %>%
  group_by(sig_varies_fitted, fit_interaction) %>%
  summarize(h_block_score = sum(cv_score) / n_distinct(cv_fold))

cv_lolo <- readRDS("analysis2/cv_lolo.rds")

unique(cv_lolo$convergence)
MyTable <- cv_lolo %>%
  group_by(sig_varies_fitted, fit_interaction) %>%
  summarize(cv_lolo_score = sum(cv_score) / n_distinct(cv_fold)) %>%
  dplyr::left_join(MyTable, .x, by = c("sig_varies_fitted", "fit_interaction"))

MyTable$`REML AIC` <- NA
MyTable$`ML AIC` <- NA

remls <- map_dbl(.x = reml, .f = function(.x) {
  return(.x$opt$AIC)
})[1:nrow(MyTable)] # get rid of random slopes model
mls <- map_dbl(.x = ml, .f = function(.x) {
  return(.x$opt$AIC)
})[1:nrow(MyTable)]

which_order <- c(
  "ar1 st reduced", "ar1 st full", "both reduced",
  "both full", "by lake reduced", "by lake full",
  "by time reduced", "by time full"
)

MyTable$`REML AIC` <- remls[order(factor(names(remls), levels = which_order))]
MyTable$`ML AIC` <- mls[order(factor(names(mls), levels = which_order))]

MyTable$fit_interaction <- ifelse(MyTable$fit_interaction == "TRUE", "Yes", "No")

MyTable$`REML AIC` <- round(MyTable$`REML AIC`, 1)
MyTable$`REML AIC` <- MyTable$`REML AIC` - min(MyTable$`REML AIC`)
MyTable$`REML AIC` <- sprintf("%0.1f", MyTable$`REML AIC`)

MyTable$`ML AIC` <- round(MyTable$`ML AIC`, 1)
MyTable$`ML AIC` <- MyTable$`ML AIC` - min(MyTable$`ML AIC`)
MyTable$`ML AIC` <- sprintf("%0.1f", MyTable$`ML AIC`)

MyTable$cv_lolo_score <- round(MyTable$cv_lolo_score, 2)
MyTable$cv_lolo_score <- sprintf("%0.2f", MyTable$cv_lolo_score)
MyTable$h_block_score <- round(MyTable$h_block_score, 2)
MyTable$h_block_score <- sprintf("%0.2f", MyTable$h_block_score)

colnames(MyTable) <- c(
  "bold(Model)", "bold(Interaction~term)",
  "bold(LOLO~CV)", "bold(H~block~CV)", paste0("bold(Delta", "~REML~AIC)"),
  paste0("bold(Delta", "~ML~AIC)")
)

g <- tableGrob(MyTable,
  rows = NULL,
  theme = ttheme_minimal(
    parse = FALSE,
    colhead = list(fg_params = list(fontsize = 12, fontface = 3, parse = TRUE))
  )
)

g <- gtable_add_rows(g, heights = unit(0.15, "npc"), pos = 0)

g <- gtable_add_grob(g,
  grobs = segmentsGrob(
    x0 = unit(0, "npc"),
    y0 = unit(0, "npc"),
    x1 = unit(1, "npc"),
    y1 = unit(0, "npc"),
    gp = gpar(lwd = 3.0)
  ),
  t = 1, b = 1, l = 1, r = 6
)

g <- gtable_add_grob(g,
  grobs = segmentsGrob(
    x0 = unit(0, "npc"),
    y0 = unit(0, "npc"),
    x1 = unit(1, "npc"),
    y1 = unit(0, "npc"),
    gp = gpar(lwd = 3.0)
  ),
  t = 2, b = 2, l = 1, r = 6
)

g <- gtable_add_grob(g,
  grobs = segmentsGrob(
    x0 = unit(0, "npc"),
    y0 = unit(0, "npc"),
    x1 = unit(1, "npc"),
    y1 = unit(0, "npc"),
    gp = gpar(lwd = 3.0)
  ),
  t = 10, b = 10, l = 1, r = 6
)

# ggsave("analysis2/table_1.png", g, width = 7, height = 4.5, dpi=600, units="in")
#---------------------
# Figure 3 parameter comparison plot:
b_j_omega <- reml$`ar1 st reduced`$rep$par.random[which(names(reml$`ar1 st reduced`$rep$par.random) == "b_j_omega")]
se_s <- as.list( reml$`ar1 st reduced`$rep, "Std. Error" )$b_j_omega

int_term <-  reml$`ar1 st full`$rep$par.random[which(names(reml$`ar1 st full`$rep$par.random) == "b_j_omega")][5]
se_int <- as.list( reml$`ar1 st full`$rep, "Std. Error" )$b_j_omega[5]

all_betas = c(b_j_omega, int_term)
se_s[5] = se_int
all_ses = se_s

#by lake models:
b_j_omega <- reml$`by lake reduced`$rep$par.random[which(names(reml$`by lake reduced`$rep$par.random) == "b_j_omega")]
se_s <- as.list( reml$`by lake reduced`$rep, "Std. Error" )$b_j_omega

int_term <-  reml$`by lake full`$rep$par.random[which(names(reml$`by lake full`$rep$par.random) == "b_j_omega")][5]
se_int <- as.list( reml$`by lake full`$rep, "Std. Error" )$b_j_omega[5]

all_betas_bylake = c(b_j_omega, int_term)
se_s[5] = se_int
all_ses_bylake = se_s

Combined = rbind(all_betas + all_ses%o%qnorm(c(0.025, 0.5, 0.975)),
                 all_betas_bylake + all_ses_bylake%o%qnorm(c(0.025, 0.5, 0.975)))
Combined[1, 1:3] = exp(Combined[1, 1:3])
Combined[6, 1:3] = exp(Combined[6, 1:3])
colnames(Combined) <- c("lower95", "MLE", "upper95")
Combined <- as.data.frame(Combined)

Combined$WhichModel <- rep(c("ar1 st", "by lake"),
                           each = length(all_betas))

Combined$WhichModel <- factor(Combined$WhichModel, levels=c(rep(c("ar1 st", "by lake"))))

Combined$WhichVariable <- rep(c("omega[0]", "beta[1] ~ Intraspecific ~ Density",
                                         "beta[2] ~ Interspecific ~ Density",
                                         "beta[3] ~ GDD",
                                         "beta[4] ~Density ~ Interaction"), 2)

Combined$WhichVariable <- factor(Combined$WhichVariable, levels=c("omega[0]", "beta[1] ~ Intraspecific ~ Density",
                                         "beta[2] ~ Interspecific ~ Density",
                                         "beta[3] ~ GDD",
                                         "beta[4] ~Density ~ Interaction"))

p <- ggplot() + theme_bw()
p <- p + geom_point(data = Combined,
                    aes(x = WhichModel,
                        y = Mean)
)

p <- ggplot() + theme_bw()
p <- p + geom_point(data = Combined,
                    aes(x = WhichModel,
                        y = MLE)
)

p <- p + geom_linerange(data = Combined,
                       aes(x = WhichModel,
                           ymax = upper95,
                           ymin = lower95))

p <- p + geom_abline(slope=0, intercept=0, lty = 1, alpha = 0.6)
p <- p + xlab("Model") + ylab("Parameter Value")

p <- p + facet_wrap( ~ WhichVariable, scales = "free_y", nrow=1,
                     labeller = label_parsed)
p <- p + theme(plot.title = element_text(hjust = 0.5, size=14),
      axis.title = element_text(size=14),
      axis.text.x = element_text(size=12),
      axis.text.y = element_text(size=12),
      legend.key=element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text.x = element_text(size = 15, colour = "black") ) +
      ggsidekick::theme_sleek()

p <- p + labs(fill = "")
p

ggsave("analysis2/Fig_3.png",
        p, width=11, height=5,
        units=c("in"), dpi = 1200 )

#---------------------
# Figure 5 Flying spaghetti monster plot:

t0 <- reml$`ar1 st reduced`$rep$par.random["global_tzero"]
eps_t0 <- reml$`ar1 st reduced`$rep$par.random[which(names(reml$`ar1 st reduced`$rep$par.random) == "eps_t0")]
ln_global_linf <- reml$`ar1 st reduced`$rep$par.random["ln_global_linf"]
eps_linf <- reml$`ar1 st reduced`$rep$par.random[which(names(reml$`ar1 st reduced`$rep$par.random) == "eps_linf")]
ln_b_sex <- reml$`ar1 st reduced`$rep$par.random["ln_b_sex"]
b_j_omega <- reml$`ar1 st reduced`$rep$par.random[which(names(reml$`ar1 st reduced`$rep$par.random) == "b_j_omega")]
log_sd <- reml$`ar1 st reduced`$opt$par["log_sd"]
eps_st <- reml$`ar1 st reduced`$rep$par.random[which(names(reml$`ar1 st reduced`$rep$par.random) == "eps_omega_st")]
eps_st <- matrix(eps_st, nrow = mesh$n, ncol = length(unique(data$Year)))
data$eps_i <- NA

data$eps <- purrr::map_dbl(1:nrow(data), function(i) {
  eps_st[data$Lake[i], data$Year[i]]
})

data$pchCode <- ifelse(data$SexCode == 0, 16, 1)
Age_Seq <- 0:26
data$colCode <- ifelse(data$SexCode == 0, "black", "white")

t_col <- function(color, percent = 50, name = NULL) {
  # 	  color = color name
  # 	percent = % transparency
  # 	   name = an optional name for the color
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
    max = 255,
    alpha = (100 - percent) * 255 / 100,
    names = name
  )
  ## Save the color
  invisible(t.col)
}

mycol1 <- t_col("steelblue", perc = 15, name = "lt.blue")
mycol2 <- t_col("darkorange", perc = 50, name = "lt.orange")

# png("analysis2/Fig_5.png",
#   width = 8, height = 6, units = "in", res = 1200
# )

par(mar = c(4.1, 4.1, 2.1, 2.1))
plot(data$FL ~ data$Age,
  lty = 2, col = "white", ylim = c(0, 80), xlim = c(-.5, 26.5),
  xlab = "Age (Years)", ylab = "Total Length (cm)", lwd = 2.5,
  las = 1, cex.lab = 1.25, cex.axis = 1.15, bty = "l", yaxs = "i", xaxs = "i"
)
points(data$FL ~ jitter(data$Age, factor = 1.5), pch = 21, bg = data$colCode, cex = 0.5)
for (i in unique(data$WBID)) {
  sub.dat <- data[which(data$WBID == i), ]
  Lake <- unique(sub.dat$Lake)
  Name <- unique(sub.dat$Name)
  for (j in unique(sub.dat$Year)) {
    sub.sub.dat <- sub.dat[which(sub.dat$Year == j), ]
    tzero <- t0 +
      eps_t0[Lake]

    omega <- exp(b_j_omega[1] +
      b_j_omega[2] * sub.sub.dat$wallEffDen.Std[1] +
      b_j_omega[3] * sub.sub.dat$compEffDen.Std[1] +
      b_j_omega[4] * sub.sub.dat$GDD.Std[1] +
      # b_j_omega[4]*sub.sub.dat$wallEffDen.Std[1]*sub.sub.dat$compEffDen.Std[1] +
      sub.sub.dat$eps[1])
    linf <- exp(ln_global_linf +
      ln_b_sex * c(0, 1) +
      eps_linf[Lake])
    Age_Seq <- min(sub.sub.dat$Age):max(sub.sub.dat$Age)
    lpred_m <- lpred_f <- NA
    for (a in 1:length(Age_Seq)) {
      lpred_m[a] <- linf[1] * (1 - exp(-(omega / linf[1]) * (Age_Seq[a] - t0)))
      lpred_f[a] <- linf[2] * (1 - exp(-(omega / linf[2]) * (Age_Seq[a] - t0)))
    }
    lines(lpred_m ~ Age_Seq, type = "l", lwd = 0.5, col = mycol1)
    lines(lpred_f ~ Age_Seq, type = "l", lwd = 0.5, col = mycol2)
  }
}

Ages <- 0:max(data$Age)
Lpreds <- exp(ln_global_linf) * (1 - exp(-(exp(b_j_omega[1]) / exp(ln_global_linf)) * (Ages - t0)))
lines(Ages, Lpreds, lty = 1, lwd = 4, col = "steelblue")

Lpreds <- exp(ln_global_linf + ln_b_sex) * (1 - exp(-(exp(b_j_omega[1]) / exp(ln_global_linf + ln_b_sex)) * (Ages - t0)))
lines(Ages, Lpreds, lty = 1, lwd = 4, col = "orange")

# dev.off()

#---------------------
# Figure 6 growth predictions in space-time
data$omega <- purrr::map_dbl(1:nrow(data), function(i) {
  omega <- exp(b_j_omega[1] +
    b_j_omega[2] * data$wallEffDen.Std[i] +
    b_j_omega[3] * data$compEffDen.Std[i] +
    b_j_omega[4] * data$GDD.Std[i] +
    data$eps[i])
})

d <- mutate(data, WBID=data$WBID, x = data$Long_c,
            y = data$Lat_c, year = data$Year + 1999,
            eps = data$eps,
            mu = data$omega) %>% group_by(WBID, year) %>%
     filter(year < 2018)

colfunc <- colorRampPalette(c("darkblue","white", "darkorange2"))
colfunc(25)
plot(rep(1,25),col=colfunc(25),pch=19,cex=3)

can1<-raster::getData('GADM', country="CAN", level=1)
alta = can1[can1$NAME_1 %in% "Alberta", ]

alta <- spTransform(alta,
                       CRS("+proj=longlat +datum=WGS84"))
alta.fort = fortify(alta)

names(alta.fort)[1] <- "Long_c"
names(alta.fort)[2] <- "Lat_c"

p <- ggplot(NULL) + theme_classic() +
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
p

# ggsave("analysis2/Fig_6.png",
#        p, dpi=1200, width=11, height=8, units=c("in"))

#-------------------------------
#Figure 7--spatial range
Kappa = exp(reml$`ar1 st reduced`$opt$par["ln_kappa"])
Range = reml$`ar1 st reduced`$rep$value["Range"]

D <- as.matrix(dist(Loc))

dis.vec <- seq(0, max(D), length = 1000)
Cor.M <- (Kappa * dis.vec) * besselK(Kappa * dis.vec, 1) #matern correlation
Cor.M[1] <- 1

# png( "analysis2/Fig_7.png",
#   width=11, height=5, res=1200, units="in")
par(mfrow=c(1,2),  mar=c(3,3,2,1), mgp=c(2,0.5,0), tck=-0.02)

plot(x = dis.vec, y = Cor.M,
     type = "l", cex.lab = 1.25, cex.axis=1.15,
     xlab = "Distance (km)", yaxs="i", xaxs="i",
     ylab = "Correlation", xlim = c(0, 110), bty="l",
     las=1)

hist(dist(Loc), xlab="Distance between Lakes (km)", main="", cex.lab = 1.25,
     breaks=40, bty="o", cex.axis=1.15, yaxs="i", las=1)
abline(v=Range, lty=3, lwd=6, col="Steelblue")
# dev.off()

#-------------------------------
#Behemoth misery multipanel plot (i.e., figure 4)
#Make the ar1 st plot:
t0 <- reml$`ar1 st reduced`$rep$par.random["global_tzero"]
eps_t0 <- reml$`ar1 st reduced`$rep$par.random[which(names(reml$`ar1 st reduced`$rep$par.random) == "eps_t0")]
ln_global_linf <- reml$`ar1 st reduced`$rep$par.random["ln_global_linf"]
eps_linf <- reml$`ar1 st reduced`$rep$par.random[which(names(reml$`ar1 st reduced`$rep$par.random) == "eps_linf")]
ln_b_sex <- reml$`ar1 st reduced`$rep$par.random["ln_b_sex"]
b_j_omega <- reml$`ar1 st reduced`$rep$par.random[which(names(reml$`ar1 st reduced`$rep$par.random) == "b_j_omega")]
log_sd <- reml$`ar1 st reduced`$opt$par["log_sd"]
eps_st <- reml$`ar1 st reduced`$rep$par.random[which(names(reml$`ar1 st reduced`$rep$par.random) == "eps_omega_st")]
eps_st <- matrix(eps_st, nrow = mesh$n, ncol = length(unique(data$Year)))
data$eps <- NA

data$eps <- purrr::map_dbl(1:nrow(data), function(i) {
  eps_st[data$Lake[i], data$Year[i]]
})

y_reduced_spatial <- purrr::map_dbl(1:nrow(data), function(i) {
  b_j_omega[2]*data$wallEffDen.Std[i] + data$eps[i]
})

y_full_spatial <- purrr::map_dbl(1:nrow(data), function(i) {
    b_j_omega[1] +
    b_j_omega[2]*data$wallEffDen.Std[i] +
    b_j_omega[3]*data$compEffDen.Std[i] +
    b_j_omega[4]*data$GDD.Std[i] +
    data$eps[i]
})

data$pres_spatial = y_full_spatial - y_reduced_spatial
plot(data$pres_spatial~data$wallEffDen.Std)
data$Block <- as.factor(data$Block)
p <- ggplot() + theme_classic()
p <- p +
  ggforce::geom_mark_ellipse(data = data[data$Block %in% c("6"),],
                       aes(y = pres_spatial, x = wallEffDen.Std),
    expand = unit(3, "mm"), lty='dotted') +
    geom_point(data = data,
                      aes(y = pres_spatial, x = wallEffDen.Std, fill=Block),
                      size = 2, pch=21)

p <- p + ggtitle("ar1 st") + ylim(2.5, 2.8)
p <- p + scale_fill_manual(values = RColorBrewer::brewer.pal(7, "Dark2"), name="Spatial \n Block")
p <- p + ylab(expression(paste(paste("Log Growth Rate ", cm%.%year^{-1} ), (omega)))) +
  xlab("Intraspecific Density (Standardized)")
p <- p + theme(text = element_text(size=15), plot.title = element_text(hjust = 0.5))
p <- p +  geom_abline(slope=b_j_omega[2], intercept=b_j_omega[1],
                      linetype=1, size=0.5, colour="black")
p <- p +  geom_abline(slope=0, intercept=b_j_omega[1],
                      linetype=2, size=0.5, colour="black")

p <- p + ggsidekick::theme_sleek() + theme(legend.position = "none")

p <- p + scale_x_continuous(breaks=c(-1.0, -0.5, 0.0, 0.5, 1.0, 1.5))
p

#################################################################
#Behemoth misery multipanel plot
#Make the `by lake` plot
t0 <- reml$`by lake reduced`$rep$par.random["global_tzero"]
eps_t0 <- reml$`by lake reduced`$rep$par.random[which(names(reml$`by lake reduced`$rep$par.random) == "eps_t0")]
ln_global_linf <- reml$`by lake reduced`$rep$par.random["ln_global_linf"]
eps_linf <- reml$`by lake reduced`$rep$par.random[which(names(reml$`by lake reduced`$rep$par.random) == "eps_linf")]
ln_b_sex <- reml$`by lake reduced`$rep$par.random["ln_b_sex"]
b_j_omega <- reml$`by lake reduced`$rep$par.random[which(names(reml$`by lake reduced`$rep$par.random) == "b_j_omega")]
log_sd <- reml$`by lake reduced`$opt$par["log_sd"]
eps_omega_lake <- reml$`by lake reduced`$rep$par.random[which(names(reml$`by lake reduced`$rep$par.random) == "eps_omega_lake")]
data$eps_omega_i <- NA

data$eps_omega_i <- purrr::map_dbl(1:nrow(data), function(i) {
  eps_omega_lake[data$Lake[i]]
})

y_reduced_by_lake <- purrr::map_dbl(1:nrow(data), function(i) {
  b_j_omega[2]*data$wallEffDen.Std[i] + data$eps_omega_i[i]
})

y_full_by_lake <- purrr::map_dbl(1:nrow(data), function(i) {
    b_j_omega[1] +
    b_j_omega[2]*data$wallEffDen.Std[i] +
    b_j_omega[3]*data$compEffDen.Std[i] +
    b_j_omega[4]*data$GDD.Std[i] +
    data$eps_omega_i[i]
})

data$pres_by_lake = y_full_by_lake - y_reduced_by_lake
plot(data$pres_by_lake~data$wallEffDen.Std)
data$Block <- as.factor(data$Block)
p1 <- ggplot() + theme_classic()
p1 <- p1 + geom_point(data = data,
                      aes(y = pres_by_lake, x = wallEffDen.Std, fill=Block),
                      size = 2, pch=21)

p1 <- p1 + ggtitle("by lake") + ylim(2.5, 2.8)
p1 <- p1 + scale_fill_manual(values = RColorBrewer::brewer.pal(7, "Dark2"), name="Spatial \n Block")
p1 <- p1 + ylab(expression(paste(paste("Log Growth Rate ", cm%.%year^{-1} ), (omega)))) +
  xlab("Intraspecific Density (Standardized)")
p1 <- p1 + theme(text = element_text(size=15), plot.title = element_text(hjust = 0.5)) #, legend.position="none")
p1 <- p1 +  geom_abline(slope=b_j_omega[2], intercept=b_j_omega[1],
                      linetype=1, size=0.5, colour="black")
p1 <- p1 +  geom_abline(slope=0, intercept=b_j_omega[1],
                      linetype=2, size=0.5, colour="black")

p1 <- p1 + ggsidekick::theme_sleek() + theme(legend.position = "none")

p1 <- p1 + scale_x_continuous(breaks=c(-1.0, -0.5, 0.0, 0.5, 1.0, 1.5))
p1
p

#Plot the spatial blocks:

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

map <- ggplot(NULL) + theme_classic( ) +
  geom_polygon(colour="black", fill="white", data=alta.fort, aes(x=Long_c, y=Lat_c, group=id)) +
  ggforce::geom_mark_ellipse(data = data[data$Block %in% c("6"),],
                       aes(Long_c, Lat_c),
    expand = unit(2, "mm"), lty='dotted') +
    geom_point( pch=21, size=2.0, data=data, aes(Long_c, Lat_c, fill=as.factor(Block)) ) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(7, "Dark2"), name="Spatial \n Block" ) +
  scale_x_continuous(breaks=c(-120, -115, -110)) +
  scale_y_continuous(breaks=c(49,52,56,60))

map <- map +
  ylab("Latitude") + xlab("Longitude") +
  theme(axis.title=element_text(size=15),
        axis.text=element_text(size=12),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing.x=unit(1, "lines"),
        panel.spacing.y=unit(0.5, "lines"))

map <- map + ggalt::coord_proj(
  paste0(CRS("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs ")))

map <- map + ggsidekick::theme_sleek()

big <- ggarrange(ggarrange(p, p1, ncol=1, nrow=2), map, widths=c(1,1))

# ggsave("analysis2/Fig_4.png", big,
#        dpi=1200, width=10, height=6.5, units=c("in"))
