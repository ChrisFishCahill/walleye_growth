#---------------------------------------------------------------------------------------------
library(dplyr)
library(plot3D)
library(plotly)
library(viridis)
library(magick)
library(lattice)
library(ggmap)
#library(arm)--> conflicts with some other stuff so I call using arm::resecale below 

#                          Christopher L. Cahill ^(;,;)^ 28 Feb 2019
#----------------------------------------------------------------------------------------------

data <- readRDS("C:/Users/Chris Cahill/Documents/GitHub/walleye/data/vB_analysis_feb_2019_cahill.rds")

mean_length <- group_by(data, Age, Year) %>%  filter(Age > 0) %>% 
	             summarize(m = mean(TL), d = median(wallEffDen))

par(mfrow=c(1,1))
plot(mean_length$m ~ mean_length$d)

setwd("C:/Users/Chris Cahill/Documents/Github/walleye/growth/vb/plots/tl_age_wdens")
#png( file="3dPlot.png", width=8, height=11, res=200, units="in")
scatter3D(arm::rescale(data$wallEffDen),   data$Age, data$TL, 
					xlab="Den", ylab="Age", zlab="TL",
					theta=50, phi = 5, 
					colvar=arm::rescale(data$wallEffDen), 
					col= viridis(30), pch=16, cex=0.5)

#----------------------------------------------------------------------------------------------
#Make some gifs-- figuring this out was miserable:
#Set the path for ImageMagick:
Sys.setenv(
	PATH = paste(
		Sys.getenv("PATH"), 
		"C:\\Program Files\\ImageMagick-7.0.8-Q16\\convert.exe",
		sep = ";"
	)
)

#Age vs TL vs Effective Walleye Density:
#png(file="example%03d.png", width=480, heigh=480)
for (i in seq(0, 350 ,5)){
	print(scatter3D(arm::rescale(data$wallEffDen),   data$Age, data$TL, 
									xlab="Den", ylab="Age", zlab="TL",
						theta=i, phi = 5, colvar=arm::rescale(data$wallEffDen), 
						col= viridis(30), pch=16, cex=1, grid=T) )
}

#dev.off()
# convert pngs to one gif using ImageMagick
#system("convert -delay 40 *.png example_2_reduced.gif")
#---------------------------------------------------------------------------------------------
#Age vs. Length vs. Effective competitor density
setwd("C:/Users/Chris Cahill/Documents/Github/walleye/growth/vb/plots/tl_age_cdens")
#png(file="example%03d.png", width=480, heigh=480)
for (i in seq(0, 350 ,5)){
	print(scatter3D(arm::rescale(data$compEffDen),   data$Age, data$TL,
									xlab="Comp", ylab="Age", zlab="TL",
									theta=i, phi = 5, colvar=arm::rescale(data$compEffDen), 
									col= viridis(30), pch=16, cex=1, grid=T) )
}

#dev.off()
#system("magick -delay 60 *.png competitor.gif")
#file.remove(list.files(pattern=".png")) #This command will get rid of the static images
#--------------------------------------------------------------------------------------------
#Age vs. TL vs. GDD

setwd("C:/Users/Chris Cahill/Documents/Github/walleye/growth/vb/plots/tl_age_gdd")

#png(file="example%03d.png", width=480, heigh=480)
for (i in seq(0, 350 ,5)){
	print(scatter3D(arm::rescale(data$GDD),   data$Age, data$TL, 
									xlab="GDD", ylab="Age", zlab="TL",
									theta=i, phi = 5, colvar=arm::rescale(data$GDD), 
									col= viridis(30), pch=16, cex=1, grid=T) )
}

#dev.off()
#system("magick -delay 60 *.png GDD.gif")
#--------------------------------------------------------------------------------------------
#Split the data up 3 equal chunks based on competitor density:
# data$fCompetitorDensity <- NA
# qs <- quantile(data$compEffDen, probs=c(0.25, 0.75))
# for(i in 1:nrow(data)){
# 	if(data$compEffDen[i] < qs[1]){
# 		data$fCompetitorDensity[i] <- "Low"
# 	}
# 	if(data$compEffDen[i] >= qs[1] && data$compEffDen[i] <= qs[2]) {
# 		data$fCompetitorDensity[i] <- "Med"
# 	}
# 	if(data$compEffDen[i] > qs[2]){
# 		data$fCompetitorDensity[i] <- "High"
#   }
# }
# 
# p <- ggplot()
# p <- p + geom_point(data=data, 
# 										aes(y=TL, x=Age), 
# 										shape=1, 
# 										size=1)
# p <- p + xlab("Age") + ylab("TL")
# p <- p + theme(text= element_text(size=15))
# p <- p + geom_smooth(data = data, 
# 										 aes(x=Age, 
# 										 		 y=TL)
#                      )
# p <- p + facet_grid(.~ fCompetitorDensity)
# p
# 
# data0 <- data %>% filter(fCompetitorDensity=="Low")
# data1 <- data %>% filter(fCompetitorDensity=="Med")
# data2 <- data %>% filter(fCompetitorDensity=="High")
# 
# #-----------------------------------------------------------------------------------
#Low competitor density:
#-----------------------------------------------------------------------------------
par(mar=c(2, 1, 1.0, 2))
setwd("C:/Users/Chris Cahill/Documents/Github/walleye/growth/vb/plots/tl_age_interaction")
#png(file="example%03d.png", width=480, heigh=480)

for (i in seq(0, 350 ,5)){
	print(scatter3D(arm::rescale(data0$wallEffDen), data0$Age, data0$TL, 
				          xlab="Den", ylab="Age", zlab="TL",
				          theta=i, phi = 5, colvar=arm::rescale(data0$wallEffDen), 
				          col= viridis(30), pch=16, cex=1, main="Low"))
}

#dev.off()
# convert pngs to one gif using ImageMagick
#system("magick -delay 60 *.png low_competitor_density.gif")

#-----------------------------------------------------------------------------------
#Medium Competitor density:
#-----------------------------------------------------------------------------------

setwd("C:/Users/Chris Cahill/Documents/Github/walleye/growth/vb/plots/tl_age_interactionMed")
#png(file="example%03d.png", width=480, heigh=480)
for (i in seq(0, 350 ,5)){
	print(scatter3D(arm::rescale(data1$wallEffDen),   data1$Age, data1$TL, 
									xlab="Den", ylab="Age", zlab="TL",
									theta=i, phi = 5, colvar=arm::rescale(data1$wallEffDen), 
									col= viridis(30), pch=16, cex=1, main="Med"))
}

#dev.off()
# convert pngs to one gif using ImageMagick
#system("magick -delay 60 *.png med_competitor_density.gif")

#-----------------------------------------------------------------------------------
#High competitor density:
#-----------------------------------------------------------------------------------
par(mar=c(1.5, 1, 1.0, 1.5))
setwd("C:/Users/Chris Cahill/Documents/Github/walleye/growth/vb/plots/tl_age_interactionHigh")
#png(file="example%03d.png", width=480, heigh=480)

for (i in seq(0, 350 ,5)){
	print(scatter3D(arm::rescale(data2$wallEffDen),   data2$Age, data2$TL,
									xlab="Den", ylab="Age", zlab="TL",
									theta=i, phi = 5, colvar=arm::rescale(data2$wallEffDen), 
									col= viridis(30), pch=16, cex=1, main="High", xlim=c(-1.2, 1.2) ) )
}

#dev.off()
# convert pngs to one gif using ImageMagick
#system("magick -delay 60 *.png High_competitor_density.gif")
#-----------------------------------------------------------------------------------

#Let's explore the spatial/spatial-temporal fin data collection program

xyplot(Y_TTM_c/1000~X_TTM_c/1000|factor(Year, labels=unique(Year)), 
			 aspect="iso", 
			 col=1, 
			 pch=16, 
			 data=data, 
			 xlab="X_km", ylab="Y_km"
			 )

#looks like the program ramps up until around 2004, and then is pretty consistent
#until 2017

#Percentage of samples in each year
table(data$Year) / nrow(data)
#Again--fairly reasonable spread

range(data$Long_c)
range(data$Lat_c)
range(data$Elevation)
glgmap <- get_map(location = c(-120, 49, -109, 60), 
									zoom = 5,
									maptype = "terrain")

p <- ggmap(glgmap) 
p <- p + geom_point(aes(Long_c, Lat_c), 
										data=data)
p <- p + xlab("Longitude") + ylab("Latitude")
p <- p + theme(text = element_text(size=8))
p <- p + facet_wrap(~Year, ncol=5)
p

setwd("C:/Users/Chris Cahill/Documents/Github/walleye/growth/vb/plots/maps")
#ggsave("Facet_FIN_2000_2019.png", dpi=500, width=8, height=11, units=c("in")) 

#Are there bigger lakes in a certain area?
p <- ggmap(glgmap) 
p <- p + geom_point(aes(Long_c, Lat_c, size=Area), 
										data=data)
p <- p + xlab("Longitude") + ylab("Latitude")
p <- p + theme(text = element_text(size=8))
p 

#Are there deeper lakes in a certain area?
p <- ggmap(glgmap) 
p <- p + geom_point(aes(Long_c, Lat_c, size=MaxDepth), 
										data=data)
p <- p + xlab("Longitude") + ylab("Latitude")
p <- p + theme(text = element_text(size=8))
p

# --> some of these depth measurements are not right.  Glenifer lake = 33 m,
# --> but data say 100 maybe some people are forgetting to convert feet to meters? 
# --> Just a guess...

#Are there old fish in a certain lakes?

#Calculate max and median age per WBID, and join with data:
adat <- data %>% group_by(WBID, Year) %>% 
	      summarize("MaxAge" = max(Age), "MedianAge" = median(Age)) 
data <- left_join(data, adat)

adat <- data %>% group_by(Year, WBID) %>% 
        tally(name="Samples")
data <- left_join(data, adat)


#Maximum age:

location <- c(-120, 49, -109, 60)
glgmap1 <- get_stamenmap(location, 
									       zoom = 5,
									       source="stamen",
									       maptype = "toner" )

p <- ggmap(glgmap1) + theme_bw()
p <- p + geom_point(aes(Long_c, Lat_c, col=MaxAge), 
										data=data)
p <- p + xlab("Longitude") + ylab("Latitude")
p <- p + theme(text = element_text(size=8))
p <- p + facet_wrap(~Year)
p
# ggsave("MaxAge.png", dpi=500, width=8, height=11, units=c("in")) 

#Median Age:

glgmap1 <- get_stamenmap(location, 
												 zoom = 5,
												 source="stamen",
												 maptype = "toner" )

p <- ggmap(glgmap1) + theme_bw()
p <- p + geom_point(aes(Long_c, Lat_c, col=MedianAge), 
										data=data)
p <- p + xlab("Longitude") + ylab("Latitude")
p <- p + theme(text = element_text(size=8))
p <- p + facet_wrap(~Year)
p
#ggsave("MedianAge.png", dpi=500, width=8, height=11, units=c("in"))

#Sample sizes in time and space:
glgmap1 <- get_stamenmap(location, 
												 zoom = 5,
												 source="stamen",
												 maptype = "toner" )

p <- ggmap(glgmap1) + theme_bw()
p <- p + geom_point(aes(Long_c, Lat_c, col=Samples), 
										data=data)
p <- p + xlab("Longitude") + ylab("Latitude")
p <- p + theme(text = element_text(size=8))
p <- p + facet_wrap(~Year)
p
#ggsave("SampleSize.png", dpi=300, width=8, height=11, units=c("in"))

#I don't see much on age plots...

#-----------------------------------------------------------------------------------
#Outliers:
source("C:/Users/Chris Cahill/Documents/Github/HighstatLibV10.R")

MyVar <- c("wallEffDen", "compEffDen", "Elevation", "GDD")

Mydotplot(data[,MyVar])
#Nothing crazy in terms of outliers

Mypairs(data[,MyVar])
#Nothing crazy here in terms of correlations or outliers

corvif(data[,MyVar])
#Variance inflation factors are beautiful :)


#Do we have enough data to include an interaction:

#Split the data up 3 equal chunks based on competitor density:
data <- as.data.frame(data)
data$fCompetitorDensity <- NA
qs <- quantile(data$compEffDen, probs=c(0.25, 0.75))
for(i in 1:nrow(data)){
	if(data$compEffDen[i] < qs[1]){
		data$fCompetitorDensity[i] <- "Low"
	}
	if(data$compEffDen[i] >= qs[1] && data$compEffDen[i] <= qs[2]) {
		data$fCompetitorDensity[i] <- "Med"
	}
	if(data$compEffDen[i] > qs[2]){
		data$fCompetitorDensity[i] <- "High"
  }
}

p <- ggplot()
p <- p + geom_point(data=data,
										aes(y=TL, x=Age),
										shape=1,
										size=1)
p <- p + xlab("Age") + ylab("TL")
p <- p + theme(text= element_text(size=15))
p <- p + geom_smooth(data = data,
										 aes(x=Age,
										 		 y=TL)
                     )
p <- p + facet_grid(.~ fCompetitorDensity)
p

#Though, may be better to look at those 3d-plots.  I think the answer is "yes", it is
#okay to include an interaction for my purposes
#----------------------------------------------------------------------------------


