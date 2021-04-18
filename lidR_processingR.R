# Processing    : Generate the the DTM, Normalize the lidar data, Filter the point clouds from the first return, remove the outlines, segmentation of trees and Buildings
# version 1.0.0 : testing

library(lidR)
library(rLiDAR)
library(raster)
############## ~~Apply a progressive morphologiocal filter "pmf"
# Use of classify_ground function to classify the point clouds into ground and non-ground
# here we are using the lidar data processed by the pdal pipeline 
#{
#"type":"filters.assign","assignment":"NumberOfReturns[:]=1"
#},
#{
#  "type":"filters.assign","assignment":"ReturnNumber[:]=1"
#}"""
  
lidR_wurz <- lidR::readLAS("./input_lidar/merged_lidar_nr.laz", select = "xyzrn")
lidR_wurz
# Using the Progressive Morphological Filter
ws  <- seq(3,12, 3)
th  <- seq(0.1, 1.5, length.out = length(ws))

lidR_wurz <- classify_ground(lidR_wurz, pmf(ws, th))
plot(lidR_wurz, color = "Classification")
lidR_wurz
############## ~~Generate the the DTM
DTM = grid_terrain(lidR_wurz,  res = 0.1,algorithm = knnidw(k = 6L, p=2))

DTM

plot(DTM)

############## ~~Normalize the lidar data = Remove the topography
lidar_aoi <- lidR_wurz
lidar_aoi <- normalize_height(lidar_aoi, DTM)
lidar_aoi
plot(lidar_aoi)
############## ~~ Filter pts clouds from the first return  classified as ground
lidar_aoi = lasfilter(lidar_aoi, Classification != 2L & ReturnNumber == 1L)
lidar_aoi
############## ~~remove pts clouds lower then or eq to 0 
lidar_aoi <- filter_poi(lidar_aoi, Z >= 0)
lidar_aoi
plot(lidar_aoi)

############## ~~generate the nDSM aka HAG
nDSM <- grid_canopy(lidar_aoi, res = .1, p2r())
nDSM
plot(nDSM)
############# Segmentation of normalized data : potential trees and buildings segmentation
lidar_seg <- lidar_aoi 
lidar_seg <- segment_shapes(lidar_seg, shp_plane(th1 = 25, th2 = 6, k = 8), "Coplanar")
plot(lidar_seg, color = "Coplanar",colorPalette = c("darkgreen", "red"))



#############  Testing with 10% of the data 
# https://github.com/Jean-Romain/lidR/issues/209
lidR_wurz2 <- lidR::readLAS("./input_lidar/merged_lidar_nr.laz",
                            filter = "-keep_random_fraction 0.1",
                            select = "xyzrn")
lidR_wurz2
plot(lidR_wurz2)
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
plot(DTM2)

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
plot(nDSM2)
############# Segmentation of normalized data : potential trees and buildings segmentation
lidar_seg2 <- lidar_aoi2 
lidar_seg2 <- segment_shapes(lidar_seg2, shp_plane(th1 = 25, th2 = 6, k = 8), "Coplanar")
lidar_seg2
plot(lidar_seg2, color = "Coplanar",colorPalette = c("darkgreen", "red"))