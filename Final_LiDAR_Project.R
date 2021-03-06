# refrence paper: https://www.researchgate.net/publication/330203798_Optimizing_the_Remote_Detection_of_Tropical_Rainforest_Structure_with_Airborne_Lidar_Leaf_Area_Profile_Sensitivity_to_Pulse_Density_and_Spatial_Sampling
#=======================================================================#
#                           #I. canopyLazR
#=======================================================================#

# Install canopyLazR from GitHub
install_github("akamoske/canopyLazR")
library(plyr)
?canopyLazR
# Load the library
library(canopyLazR)
library(uuid)
library(rlas)
library(devtools)

library(lidR)
library(raster)
library(sp)
library(rgdal)
library(ggplot2)
library(viridisLite)
library(rasterVis)
library(grid)
library(scales)
library(RStoolbox)

setwd(choose.dir("D:/Lidar2/openTopography/pr5"))
#1. Convert .laz or .las file into a voxelized lidar array
?laz.to.array
laz.data <- laz.to.array("points.laz", 
                         voxel.resolution = 10, 
                         z.resolution = 1,
                         use.classified.returns = F)

head(laz.data)
#2. Level the voxelized array to mimic a canopy height model
?canopy.height.levelr
level.canopy <- canopy.height.levelr(lidar.array = laz.data)
#3. Estimate LAD for each voxel in leveled array
?machorn.lad()
lad.estimates <- machorn.lad(leveld.lidar.array = level.canopy, 
                             voxel.height = 1, 
                             beer.lambert.constant = NULL)

lad.estimates
#4. Convert the LAD array into a single raster stack
lad.raster <- lad.array.to.raster.stack(lad.array = lad.estimates, 
                                        laz.array = laz.data, 
                                        epsg.code = 32611)
lad.raster
plot(lad.raster)
#5. Create a single LAI raster from the LAD raster stack
lai.raster <- raster::calc(lad.raster, fun = sum, na.rm = TRUE)
lai.raster
plot(lai.raster , main="LAI"~(m^2~m^-2))
gplot(lai.raster)+
  geom_raster(aes(x=x, y=y, fill=value))+
  scale_fill_viridis_c()+
  coord_equal()+
  theme_classic()+
  theme(text = element_text(size = 14))+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(title = "LAI"~(m^2~m^-2))
#6. Convert the list of LAZ arrays into a ground and canopy height raster
?array.to.ground.and.canopy.rasters()
grd.can.rasters <- array.to.ground.and.canopy.rasters(laz.data, 32611)

#6.a/ Plot the ground raster DTM
plot(grd.can.rasters$ground.raster, main="DTM")

#6.b/ Plot the canopy height raster DSM
plot(grd.can.rasters$canopy.raster, main="DSM")
gplot(grd.can.rasters$canopy.raster)+
  geom_raster(aes(x=x, y=y, fill=value))+
  scale_fill_viridis_c()+
  coord_equal()+
  theme_classic()+
  theme(text = element_text(size = 14))+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(title = "DSM")
#6.c/ Plot the canopy height model raster CHM
chm<- grd.can.rasters$chm.raster * 0.3048
plot(chm, main= "CHM")
gplot(chm)+
  geom_raster(aes(x=x, y=y, fill=value))+
  scale_fill_viridis_c()+
  coord_equal()+
  theme_classic()+
  theme(text = element_text(size = 14))+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(title = "CHM (m)")

#7. Calculate max LAD and height of max LAD
?lad.ht.max()
max.lad <- lad.ht.max(lad.array = lad.estimates, 
                      laz.array = laz.data, 
                      ht.cut = 5, #Height that calculations will exclude
                      epsg.code = 32611) #CRS

max.lad
#7.a/ Plot the max LAD raster
LAD<- max.lad$max.lad.raster*0.3048
plot(LAD, main="Max LAD"~(m^2~m^-3))
#7.b/ Plot the height of max LAD raster
CHM_MAX_LAD<- max.lad$max.lad.ht.raster *0.3048
plot(CHM_MAX_LAD, main="Height of max LAD (m)")

#Examine the Vertical distribution of LAD
df<- data.frame(LAD=values(LAD),
                Height_MAX_LAD=values(CHM_MAX_LAD))
#VIZ 
g<-ggplot(df, aes(LAD, Height_MAX_LAD))+
  geom_point(pch=21, cex=2, fill="blue" )+
  labs(title = "Vartical distribution of LAD", 
       x="LAD"~(m^2~m^-3),y="Canopy Height of maximum LAD (m)")

g # Note that there are NA voxels which were are not considered for LAD calculation
library(ggExtra)
ggMarginal(g,type = "hist")
ggMarginal(g, type = "density", fill="lightgray")
ggMarginal(g, type = "boxplot", fill="lightgray")


#=======================================================================#
#                           #II. lidR
#=======================================================================#
library("lidR")
lidar<- readLAS("points.laz")
lidar
plot(lidar)

############################### 1. Digital Terrain Model  ##############################################
?grid_terrain()
col <- height.colors(50)
dtm1 = grid_terrain(lidar, algorithm = knnidw(k = 6L, p=2))
dtm1
dtm2 = grid_terrain(lidar, algorithm = tin())
dtm2
dtm3 = grid_terrain(lidar, algorithm = kriging(k = 10L))
dtm3

## Not run:
plot(dtm1,main="DTM1")
plot(dtm2, main="DTM2")
plot(dtm3, main="DTM3")
#3D ViZ
plot_dtm3d(dtm1)
plot_dtm3d(dtm2)
plot_dtm3d(dtm3)

################################# 2. Digital Surface Model ####################################
# the function to Create a canopy surface model (i.e. canopz height model, CHM) using a point cloud. For each pixel the function returns
#the highest point found (point-to-raster)
?grid_canopy
col <- height.colors(50)

## Points-to-raster algorithm with a resolution of 1 meter
DSM1 = grid_canopy(lidar, res= 1, algorithm = p2r())
DSM1
plot(DSM1,col=col ,main= "DSM 1")

# Local maximum algorithm with a resolution of 1 meter replacing each
# point by a 20 cm radius circle of 8 points
DSM2 = grid_canopy(lidar, res=1, algorithm = p2r(subcircle = 0.2))
DSM2
plot(DSM2,col=col, main="DSM 2")

# Points-to-raster algorithm with a resolution of 0.5 meters replacing each
# point by a 20-cm radius circle of 8 points
DSM3 <- grid_canopy(lidar, res = 0.5, p2r(0.2))
DSM3
plot(DSM3, col = col, main="DSM 3")

#3D VIZ
plot_dtm3d(DSM1)
plot_dtm3d(DSM2)
plot_dtm3d(DSM3)

# Local maximum algorithm with a resolution of 1 meter replacing each
# point by a 10 cm radius circle of 8 points and interpolating the empty
# pixels using the 3-nearest neighbours and an inverse-distance weighting.
DSM4 = grid_canopy (lidar, 1, algorithm = dsmtin(max_edge = 10))
plot(DSM4, col=col, main="DSM_4")
plot_dtm3d(DSM4)
# Basic triangulation and rasterization of first returns
DSM5 <- grid_canopy(lidar, res = 0.5, dsmtin())
plot(DSM5, col = col, main="DSM_5")
plot_dtm3d(DSM5)

# Khosravipour et al. pitfree algorithm
DSM6 <- grid_canopy(lidar, res = 0.5, pitfree(c(0,2,5,10,15), c(0, 1.5)))
plot(DSM6, col = col, main="DSM_6")
plot_dtm3d(DSM6)

################################# 3. Canopy Heigth Model (CHM) ####################################
CHM_LAZ<- (DSM1 - dtm1)*0.3048
CHM_LAZ
plot(CHM_LAZ, main=" CHM (m)")


################################# 4. Map the pulse or point density  ######################################################################
# function to spatially map the pulse or point density in a lidar point cloud
?grid_density
#Point density Map with Grid size of 5
d = grid_density(lidar, res=5)
plot(d, main="Pulse at 5m grid size")

# with another grid size
d = grid_density(lidar, 10)
plot(d, main="Pulse at 10m grid size")


################################  5. Clip LiDAR points   ##########################################################
# function to clip a given point cloud
?lasclipRectangle
extent(lidar)
subset1 = lasclipRectangle(lidar, 571900, 1013800, 572200, 1014100)
subset2 = lasclipRectangle(lidar, 571900, 1013800, 572100, 1013900)
plot(subset1)
plot(subset2)

###############################  6. Decimate a LAS object  ###########################################################
# function to thin a point cloud based on defined point density of pattern
#Decimate a LAS object: Reduce the number of points
?lasfilterdecimate
lidr = readLAS("points.laz", select = "xyz")
lidr
# Select the highest point within each pixel of an overlayed grid
thinned1 = lasfilterdecimate(lidr, algorithm=highest(5))
plot(grid_density(thinned1))
plot(thinned1)

## Select points randomly to reach an overall density of 1
# By default the method is homogenize = TRUE
thinned2 = lasfilterdecimate(lidr, algorithm = random(density = 1))
plot(grid_density(lidr))
plot(grid_density(thinned3))
## Select points randomly to reach an homogeneous density of 1
# Method homogenize = FALSE enables a global pulse density to be reached
thinned3 = lasfilterdecimate(lidr, algorithm = homogenize(density = 1,res = 5))
plot(grid_density(thinned3))
plot(thinned3)

################## 7. lasfilter:  Return points with matching conditions ############################################
?lasfilter
# Select the first returns classified as ground
firstground1 = lasfilter(lidar, Classification == 2L & ReturnNumber == 1L)
firstground1
plot(firstground1)
# Multiple arguments are equivalent to &
firstground = lasfilter(lidar, Classification == 2L, ReturnNumber == 1L)
firstground
plot(firstground)
# Multiple criteria
first_or_ground = lasfilter(lidar, Classification == 2L | ReturnNumber == 1L)
first_or_ground
plot(first_or_ground)

################################# 8. Smooth a point cloud ##########################################
?lassmooth
lidr = readLAS("points.laz", select = "xyz")
lidr
lidr <- lasfiltersurfacepoints(lidr, 1)
plot(lidr)

lidr <- lassmooth(lidr, size=5,method =  "gaussian",
                  shape = "circle", sigma = 2.5)
plot(lidr)          # =size/6  

lidr <- lassmooth(lidr, size=5, method = "gaussian", "square", sigma = 2.5)
plot(lidr)

lidr <- lasunsmooth(lidr)
plot(lidr)

#END... The rest of the script is just for trial purposes


################################ 9. Classify points as 'ground' or 'not ground' ############################################
# this function classifies the point cloud into ground and non-ground points
?lasground
rad1 = readLAS("points.laz", select = "xyzrn")
rad1
plot(rad1)


# 1/Using the Progressive Morphological Filter
# --------------------------------------

ws  <- seq(3,12, 3)
th  <- seq(0.1, 1.5, length.out = length(ws))

rad1 <- lasground(rad1, pmf(ws, th)) #pmf(): a ground-segmentation function
plot(rad1, color = "Classification")

#' 2/# Using the Cloth Simulation Filter
# --------------------------------------

# (Parameters chosen mainly for speed)
mycsf <- csf(TRUE, 1, 1, time_step = 1)
rad2 = readLAS("points.laz", select = "xyzrn")
rad2
rad2 <- lasground(rad2, mycsf)
plot(rad2, color = "Classification")








############10. Individual tree segmentation + Compute the hull of each tree##################################
?lastrees()
LASfile <- system.file("extdata", "MixedConifer.laz", package="lidR")
las <- readLAS(LASfile, select = "xyz", filter = "-drop_z_below 0")

# Using Li et al. (2012)
las <- lastrees(las, li2012(R = 3, speed_up = 5))
plot(las, color = "treeID")

laz = readLAS("points.laz", select = "xyz", filter = "-drop_z_below 0")
plot(laz)
#A/ Using Li et al. (2012)
laz <- lastrees(laz, li2012(R = 3, speed_up = 5))
plot(laz, color = "treeID")
plot(laz, color = "treeID", colorPalette = pastel.colors(200))

# Compute the hull of each segmented tree.
?tree_hulls
# Only the hulls
convex_hulls = tree_hulls(laz)
plot(convex_hulls)
# The hulls + some user-defined metrics
convex_hulls = tree_hulls(laz, func = ~list(Zmax = max(Z)))
spplot(convex_hulls, "Zmax")

# The bounding box
bbox_hulls = tree_hulls(laz, "bbox")
plot(bbox_hulls)

## Not run: 
install.packages("concaveman")
library(concaveman)
concave_hulls = tree_hulls(laz, "concave")   ##It takes forever time
sp::plot(concave_hulls)
## End(Not run)

#Silva 2016
laz_2 = readLAS("points.laz", select = "xyz", filter = "-drop_z_below 0")
laz_2
col <- pastel.colors(200)

chm <- grid_canopy(laz_2, res = 0.5, p2r(0.3))
ker <- matrix(1,3,3)
chm <- raster::focal(chm, w = ker, fun = mean, na.rm = TRUE)

ttops_2 <- tree_detection(chm, lmf(4, 2))
laz_2  <- lastrees(laz_2, silva2016(chm, ttops_2))
plot(laz_2, color = "treeID", colorPalette = col)

#dalponte2016 
laz_3 = readLAS("points.laz", select = "xyz", filter = "-drop_z_below 0")

laz_3   <- lastrees(laz_3, dalponte2016(chm, ttops_2))
plot(laz_3, color = "treeID", colorPalette = col)






# Importing LAS file:
LASfile <- system.file("extdata", "LASexample1.las", package="rLiDAR")

# Reading LAS file
LA<-readLAS(LASfile,short=TRUE)
LA
# Height subsetting the data
xyz<-subset(LAS[,1:3],LAS[,3] >= 1.37)

# Getting LiDAR clusters
set.seed(1)
clLAS<-kmeans(xyz, 32)

# Set the points id 
id<-as.factor(clLAS$cluster)

# Set the xyid input
xyid<-cbind(xyz[,1:2],id)

# Compute the LiDAR convex hull of the clusters 
chullTrees<-chullLiDAR2D(xyid)

# Plotting the LiDAR convex hull
library(sp)
plot(SpatialPoints(xyid[,1:2]),cex=0.5,col=xyid[,3])
plot(chullTrees$chullPolygon,add=TRUE, border='green')

# Get the ground-projected area of LiDAR convex hull
chullList<-chullTrees$chullArea 
summary(chullList) # summary 