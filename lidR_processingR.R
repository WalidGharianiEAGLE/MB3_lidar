# Processing    : Generate the the DTM, Normalize the lidar data, Filter the point clouds from the first return, remove the outlines, segmentation of trees and Buildings
# version 1.0.0 : testing
# Author : Walid Ghariani

library(lidR)
library(rLiDAR)
library(raster)
library(rasterVis)
############## ~~Read the laz file: We choosed 4 aois with most coverage on buildings to merge and work with 
# This data are not provided in github dues to the size of teh files
lid1<- lidR::readLAS("input_lidar/laser_1_2021-02-25-10-13-37_392_RemoveOutliers.las")
lid2<- lidR::readLAS("input_lidar/laser_1_2021-02-25-10-13-37_608_RemoveOutliers.las")
lid3<- lidR::readLAS("input_lidar/laser_1_2021-02-25-10-13-37_344_RemoveOutliers.las")
lid4<- lidR::readLAS("input_lidar/laser_1_2021-02-25-10-13-37_704_RemoveOutliers.las")

############## ~~Merge the data into a single las/laz file  
lid_wurz <- rbind(lid1, lid2, lid3, lid4)
# make a copy
lidR_wurz <- lid_wurz
lid_wurz
#plot(lidR_wurz)
writeLAS(lidR_wurz, "./output_lidar/merged_lidar.laz")
############## ~~ This data will be then used with pdal 
# Use of classify_ground function to classify the point clouds into ground and non-ground
# here we are using the lidar data processed by the pdal pipeline 
#{
#"type":"filters.assign","assignment":"NumberOfReturns[:]=1"
#},
#{
#  "type":"filters.assign","assignment":"ReturnNumber[:]=1"
#}"""

#############  Testing with 10% of the data 
# https://github.com/Jean-Romain/lidR/issues/209
lidR_wurz2 <- lidR::readLAS("./input_lidar/merged_lidar_nr.laz",
                            filter = "-keep_random_fraction 0.1",
                            select = "xyzrn")
lidR_wurz2
plot(lidR_wurz2)

writeLAS(lidR_wurz2, "./output_lidar/merged_lidar_nr_10percent.laz")
# Using the Cloth Simulation Filter
# --------------------------------------

# (Parameters chosen mainly for speed)
mycsf <-  csf(rigidness = 3, cloth_resolution = 1)
lidR_wurz2 <- classify_ground(lidR_wurz2, mycsf)
lidR_wurz2
plot(lidR_wurz2, color = "Classification")
# 2 option => pmf
ws  <- seq(3,12, 3)
th  <- seq(0.1, 1.5, length.out = length(ws))

#lidR_wurz2 <- classify_ground(lidR_wurz2, pmf(ws, th))

############## ~~Generate the the DTM
DTM2 = grid_terrain(lidR_wurz2,  res = .1,algorithm = knnidw(k = 6L, p=2))
DTM2

plot(DTM2, main = "DTM at Wurzburg AOI")

gplot(DTM2)+
  geom_raster(aes(x=x, y=y, fill=value))+
  scale_fill_viridis_c()+
  coord_equal()+
  theme(text = element_text(size = 12))+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(title = "DTM at Wurzburg AOI",
       x="Longitude", y="Latitude")

############## ~~Normalize the lidar data = Remove the topography
lidar_aoi2 <- lidR_wurz2
lidar_aoi2 <- normalize_height(lidar_aoi2, DTM2)
lidar_aoi2
plot(lidar_aoi2)
############## ~~ Filter pts clouds from the first return  classified as ground
lidar_aoi2 = lasfilter(lidar_aoi2, Classification != 2L & ReturnNumber == 1L)
lidar_aoi2
############## ~~remove pts clouds lower then or eq to 0 
lidar_aoi2 <- filter_poi(lidar_aoi2, Z >= 0)
lidar_aoi2
plot(lidar_aoi2)
############## ~~generate the nDSM aka HAG
nDSM2 <- grid_canopy(lidar_aoi2, res = .1, p2r())
nDSM2
plot(nDSM2, main ="nDSM at Wurzburg AOI")

gplot(nDSM2)+
  geom_raster(aes(x=x, y=y, fill=value))+
  scale_fill_gradient(low="#FFFF80", high="#FF0000")+
  coord_equal()+
  theme(text = element_text(size = 12))+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(title = "nDSM at WurzburgAOI",
       x="Longitude", y="Latitude")
############# Segmentation of normalized data : potential trees and buildings segmentation
lidar_seg2 <- lidar_aoi2 
lidar_seg2
lidar_seg2 <- segment_shapes(lidar_seg2, shp_plane(th1 = 25, th2 = 6, k = 64), "Coplanar")#k=174 when "-keep_random_fraction 0.3" 
lidar_seg2
plot(lidar_seg2, color = "Coplanar",colorPalette = c("darkgreen", "red"))

############# Automate the whole prpocess
automated_Seg <- function (input){
  lidR_wurz2 <- lidR::readLAS(input,
                            filter = "-keep_random_fraction 0.1",
                            select = "xyzrn")
  mycsf <-  csf(rigidness = 3, cloth_resolution = 1)
  lidR_wurz2 <- classify_ground(lidR_wurz2, mycsf)
  DTM2 = grid_terrain(lidR_wurz2,  res = .1,algorithm = knnidw(k = 6L, p=2))
  writeRaster(DTM2 , "./output_lidar/DTM_lidR.tif")
  lidar_aoi2 <- lidR_wurz2
  lidar_aoi2 <- normalize_height(lidar_aoi2, DTM2)
  lidar_aoi2 = lasfilter(lidar_aoi2, Classification != 2L & ReturnNumber == 1L)
  lidar_aoi2 <- filter_poi(lidar_aoi2, Z >= 0)
  nDSM2 <- grid_canopy(lidar_aoi2, res = .1, p2r())
  writeRaster(nDSM2 , "./output_lidar/nDSM_lidR.tif")
  lidar_seg2 <- lidar_aoi2 
  lidar_seg2 <- segment_shapes(lidar_seg2, shp_plane(th1 = 25, th2 = 6, k = 64), "Coplanar")
  plot(lidar_seg2, color = "Coplanar",colorPalette = c("darkgreen", "red"))
}
automated_Seg("./input_lidar/merged_lidar_nr.laz")

########### 3Dviz with rLiDAR
# save the filtered data as .las
writeLAS(lidar_aoi2 , "./Output_lidar/lidar_filtered.las")

rLAS <- rLiDAR::readLAS("./Output_lidar/lidar_filtered.las",short=TRUE)

summary(rLAS)

# Lest Viz the data on a 3D

# Define the color ramp
# color ramp
myColorRamp <- function(colors, values) {
  v <- (values - min(values))/diff(range(values))
  x <- colorRamp(colors)(v)
  rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
}

# Color by height
col <- myColorRamp(c("blue","green","yellow","red"),rLAS[,3])

# plot 2D
plot(rLAS[,1], rLAS[,2], col=col, xlab="UTM.Easting", ylab="UTM.Northing", main="Color by height")

# plot 3D
library("rgl")
points3d(rLAS[,1:3], col=col, axes=FALSE, xlab="", ylab="", zlab="")
axes3d(c("x+", "y-", "z-"))                     # axes
grid3d(side=c('x+', 'y-', 'z'), col="gray")     # grid
title3d(xlab = "UTM.Easting", ylab = "UTM.Northing",zlab = "Height(m)", col="red") # title
planes3d(0, 0, -1, 0.001, col="gray",alpha=0.7) # terrain
