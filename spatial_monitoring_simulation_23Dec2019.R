library(raster)
library(rgdal)
library(spatstat)
library(sp)
library(maptools)
library(spatialEco)
library(SDraw) #srs functions from SDraw package sample from point objects and retains attributes
library(velox) #faster raster handling

options(scipen=999)

# source('E:/Lucretia/Ferruginous_Hawk_Project/Monitoring/Raptor_Monitoring/mode_func.R')
# source('E:/Lucretia/Ferruginous_Hawk_Project/Monitoring/Raptor_Monitoring/hex_grid_func.R')
# 
# #setwd
# setwd("E:/Lucretia/Ferruginous_Hawk_Project/Monitoring")
# 
# #import species RSFs
# hawk.rsf <- raster("largehawksGLM3_7bins.tif")
# 
# #import townships
# townships <- readOGR(dsn="E:/Lucretia/Ferruginous_Hawk_Project/WY_background_layers",layer="townships_all_possible")
# #import boundaries
# state <- readOGR(dsn="E:/Lucretia/Ferruginous_Hawk_Project/WY_background_layers",layer="State_Boundary")
# hawk.bound <- readOGR(dsn="E:/Lucretia/Ferruginous_Hawk_Project/WY_background_layers/Wyo_ecoregions_EPA",layer="study_areas")


####Zach's computer
source('C:/Users/zwallac2/Dropbox (UW WYNDD)/Proj_GovESA_FEHA_GOEA_Monitoring/Data/Analyses/Raptor_Monitoring/mode_func.R')
# source('C:/Users/zwallac2/Dropbox (UW WYNDD)/Proj_GovESA_FEHA_GOEA_Monitoring/Data/Analyses/Raptor_Monitoring/hex_grid_func.R') #ZW function is faster

#ditectories and paths
setwd("C:/Users/zwallac2/Dropbox (UW WYNDD)/Proj_GovESA_FEHA_GOEA_Monitoring/Data/Analyses/Raptor_Monitoring")
shpOut <- ("C:/Users/zwallac2/Dropbox (UW WYNDD)/Proj_GovESA_FEHA_GOEA_Monitoring/Data/Analyses/Raptor_Monitoring_shpOut")

#import species RSFs
hawk.rsf <- raster("C:/Users/zwallac2/Dropbox (UW WYNDD)/Proj_GovESA_FEHA_GOEA_Monitoring/GIS/RSF/largehawksGLM3_7bins.tif")

#import townships
townships <- readOGR(dsn="C:/Users/zwallac2/Dropbox (UW WYNDD)/Proj_BI_USFS_FEHA_Movement/GIS/Base_Data/townships",layer="townships_all_possible")

#import boundaries
state <- readOGR(dsn="C:/Users/zwallac2/Dropbox (UW WYNDD)/Proj_GovESA_FEHA_GOEA_Monitoring/GIS/Base_Data/Administrative",layer="Wyoming_Boundary")
hawk.bound <- spTransform(readOGR(dsn="C:/Users/zwallac2/Dropbox (UW WYNDD)/Proj_GovESA_FEHA_GOEA_Monitoring/GIS/Base_Data/Study_area_shapefile",layer="study_areas"),proj4string(state))
#### End Zach's computer

#check crs
identicalCRS(hawk.rsf, townships)
identicalCRS(state, townships)
identicalCRS(hawk.bound, townships)

#plot
plot(hawk.rsf)
plot(state, add=T)
plot(hawk.bound, add=T)
plot(townships, add=T)

#convert raster to points
hawk.rsf.pts <- rasterToPoints(hawk.rsf, spatial=TRUE)
colnames(hawk.rsf.pts@data) <- "rsf" #rename rsf bins column
hawk.rsf.pts <- hawk.rsf.pts[complete.cases(hawk.rsf.pts@data),] #drop NA rows, if necessary
#not necessary to clip points by study area polygon, but could be for goea
saveRDS(hawk.rsf.pts, file="hawkRsfPts") #just in case because rasterToPoints is slow
# hawk.rsf.pts <- readRDS("hawkRsfPts")


# Set paramaters for simulation
## parms not yet included in code are commented out
## moved all parameter values here before simulation loop

## Simulation rep
n.reps <- 1

## Nests
reoccu.mean <- 0.60 #Average reoccupancy rate from 2011-2013 study
n.nests.mean <- round(sum(hawk.bound@data$area_km2 / c(83,127.1,85.5)) / reoccu.mean, 0) #area of each ecoregion / estimated mean density from Olson et al. 2015
# N nests at yr0 is calculated as (1) the sum of the area of each ecoregion (2) divided by the density of occupied nests (km2/nest) per region estimated from distance sampling in 2010-2011 
# (3) Divided by the average re-occupancy rate during the 2011-2013 study
# Step 3 is necessary because the density estimates from  distance sampling are for occupied nests only.
# If we started with a sample size of N occupied nests at t0 and multipled that by our chosen re-occupancy rate, we would underestimate the number of nests at t1.
# Instead, we estimate population of nests that is 0.60 occupied at t0.
# Do you agree that this makes sense?

n.nests.sd <- 100 #stdev of n nests
min.dist.mean <- 1500 # minimum distance between nests (m)
min.dist.sd <- 100 #stdev of inter-nest distance

## RSF bins
rsf.bins <- c(1,2,3,4,5,6,7) #rsf raster bins
nest.probs <- c(0,.02,.03,.09,.14,.2,.52) #proportion of nests within each bin

## Demographic rates
occupancy <- 0.6
# occupancy.sd <- 0.1
success <- 0.6
# success.sd <- 0.1
young <- 0.7
# young.sd <- 0.2

# decline <- 0.20 #over n.years

## Detection rates
occ.det <- 0.5
# occ.prod.det <- 0.80 #detection probability for occupancy with reporductive success
# occ.misclass <- 0.10 #probability of misclassifying a site with reproduction as not reproducing

## Effort
n.years <- 10
# samp.freq <- 2 #Every two years
n.samples <- 200 #number of hexagons sampled
n.vis.occ <- 4
# n.vis.prod <- 2

hex.area <- 1200 #area of hexagons (ha)
hex.cluster <- "Y"
# strat.type <- c("none","rsf")


# Make hexagon grid
##moved up to keep outside of loop
##faster than other function. 1200-ha grid took 6.7 mins with this code
grid <- st_make_grid(st_as_sf(hawk.bound), #convert to sf type
                     cellsize = 2*sqrt(hex.area*1e+4/(2*sqrt(3))), #calculate short axis (m) from area (ha)
                     square = FALSE)
##clean up attributes
hex.grid <- spTransform(as(grid, 'Spatial'), proj4string(townships)) #convert back to st and reproject
hex.ids <- sapply(slot(hex.grid, "polygons"), function(x) slot(x, "ID")) #get IDs from hex.grid
hex.df <- data.frame(ID=1:length(hex.grid), row.names=hex.ids) #make data frame of IDs
hex.grid <- SpatialPolygonsDataFrame(hex.grid, hex.df) #join dataframe to spatial polygons

##summarize rsf within hex grid -- this will be used for stratified sampling
#get average (mode) RSF value in each hex
rsf.v <- velox(hawk.rsf)
hex.modes <- rsf.v$extract(sp=hex.grid,fun=getmode) #use getmode function ##added na.omit to mode_func.R to work for hexagons that partially overlap rsf
hex.grid@data$rsfMode <- hex.modes[,1]
saveRDS(hex.grid, file="hexGrid") #save
# writeOGR(hex.grid, shpOut, "hexGridTest", driver="ESRI Shapefile")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation loop

##store all results in a big list
results <- list()

set.seed(11)

for (j in 1:n.reps){

#initial nest sample
#draw parameter values
n.nests <- round(rnorm(1,n.nests.mean,n.nests.sd),0)
min.dist <- rnorm(1,min.dist.mean,min.dist.sd)
size.samp <- round(nest.probs*n.nests,0) #sample size per bin

#loop over rsf bins
d.list <- list() #empty list
for (i in 1:length(rsf.bins)){
  d.list[[i]] <- srs.point(hawk.rsf.pts[hawk.rsf.pts@data$rsf==rsf.bins[i],],n=size.samp[i]) #srs.point from SDraw package
  }
d.list[sapply(d.list,is.null)] <- NULL #drop any bins with 0 nests
samp <- do.call(rbind,d.list) #bind into single dataframe
samp@data$sampleID <- 1:nrow(samp@data) #renumber

# #check
# length(samp);n.nests
# table(samp$rsf); size.list
# plot(samp, col=samp@data$rsf)

#for any pair of nests that are closer than the min distance, drop one
#messy code, but works
dmat <- spDists(samp,diagonal=F) #matrix of distances
pairs <- data.frame(which(dmat<=min.dist & dmat>0, arr.in=TRUE)) #pairs of points with <min distance
pairs.unique <- data.frame(unique(t(apply(pairs, 1, sort)))) #unique pairs only
colnames(pairs.unique) <- colnames(pairs) #replace column ids lost with apply
pairs.unique$rsf1 <- merge(pairs.unique,samp@data,by.x="row",by.y="sampleID",all.x=TRUE,all.y=FALSE)[,4] #add rsf value for rows
pairs.unique$rsf2 <- merge(pairs.unique,samp@data,by.x="col",by.y="sampleID",all.x=TRUE,all.y=FALSE)[,5] #add rsf value for cols
pairs.unique$diff <- pairs.unique$rsf1-pairs.unique$rsf2 #add difference in rsf values

#apply function to choose one of each pair with highest rsf
drops <- apply(pairs.unique, 1, function(x){
  if (x['diff']<0){x['row']
  } else if (x['diff']>0){x['col']
      } else if (x['diff']==0){sample(x[1:2],1)} #choose randomly if rsf values of pair are equal
})

samp.thin <- samp[!(samp$sampleID %in% drops),] #thinned nest sample spatial data file
drops.sp <- samp[(samp$sampleID %in% drops),] #dropped nests spatial data file

#check number of nests per township
nests.town <- over(samp.thin,townships)
nests.per <- table(as.numeric(table(nests.town$joinfield)))

# This approach works, but does not draw more nests to replace those dropped. This decreases the sample size  by ~3% for 1500 nests
# We could consider adding code to replace nests and match the sample size exactly


# Demographic data for simulated nests
sim.nests <- samp.thin@data #make table
sim.nests$X <- samp.thin$x
sim.nests$Y <- samp.thin$y

demo.list <- list()
for (i in 1:n.years){
  year <- i
  sampleID <- sim.nests$sampleID
  n.nests <- nrow(sim.nests)
  occ <- rbinom(n.nests,1,occupancy) #Add hyper-parameter for year-specific occupancy SE 
  succ <- rbinom(n.nests,1,success)*occ #Removed some of the if statments by multiplying by 0-1 vectors
  yng <- rpois(n.nests,1)*succ

   for (j in seq_along(1:n.vis.occ)){
     vis <- rbinom(n.nests,1,occ.det)*occ
     if(j==1){
       det.hist <- vis
     } else{
       det.hist <- cbind(det.hist, vis)
     }
   }
   colnames(det.hist) <- paste0(rep("vis",n.vis.occ),1:n.vis.occ)
  
   head(det.hist)
  
  demo.list[[i]] <- cbind(sampleID,year,occ,succ,yng,det.hist)
}
demo.dat <- data.frame(do.call(rbind,demo.list))

# Summarize results
# library(plyr)
# ddply(demo.dat, .(year), summarise, meanOccuYear=mean(occ), nOccuYear=sum(occ))



# Sample nests with hexagons
## Works with or without clusters

# # hex.cluster =="N"
# hex.hip <- hip.polygon(hex.grid, n.samples) #2.92 min. Faster than hip.point on centroids for some reason
# colnames(hex.hip@data)[1] <- "sampleOrder" #rename master sample order column
# hex.samp <- hex.grid[hex.grid@data$ID %in% hex.hip@data$ID,] #subset hex with sample points

# hex.cluster=="Y"
hex.hip <- hip.polygon(hex.grid, n.samples/7) #divide n.samples by 7 to match final sample size with 'flower' clusters
colnames(hex.hip@data)[1] <- "sampleOrder" #rename master sample order column
hex.samp <- hex.grid[hex.grid@data$ID %in% hex.hip@data$ID,] #subset hex with sample points
mat <- gIntersects(hex.samp, hex.grid, byid=T)
pairs <- data.frame(which(mat==TRUE, arr.in=TRUE))
pairs$sampHex <- rep(colnames(mat),each=7)
pairs <- pairs[order(pairs$row),]
hex.clust <- hex.grid[hex.grid@data$ID %in% pairs$row,]
hex.clust@data$cluster <- pairs$col
hex.clust$type <- rep("cluster",nrow(hex.clust@data))
hex.clust$type[hex.clust$ID %in% hex.samp$ID] <- "sample"
hex.samp <- hex.clust
# plot(hawk.bound, col="grey90")
# plot(hex.samp, border=hex.samp$cluster, add=T)

#nests in hexagons
nest.hex <- over(samp.thin, hex.samp)
nest.samp.tab <- nest.hex[!is.na(nest.hex$ID),]
nest.samp.tab$geometryID <- rownames(nest.samp.tab)
nest.samp <- samp.thin[samp.thin$geometryID %in% nest.samp.tab$geometryID,]
nest.samp <- merge(nest.samp, nest.samp.tab, by="geometryID", all.x=T, all.y=F)
colnames(nest.samp@data)[4] <- "hexID"

# plot(hex.samp)
# plot(nest.samp,col="blue",add=T)

## summarize
n.nests.samp <- nrow(nest.samp) #N nests in sample
n.nests.hex <- table(as.numeric(table(nest.hex$ID))) #table of N nests per hex
n.nests.cluster <- table(as.numeric(table(nest.hex$cluster))) #table of N nests per hex
## save
# writeOGR(samp.thin, shpOut, "nestSample", driver="ESRI Shapefile")
# writeOGR(hex.samp, shpOut, "hexSample", driver="ESRI Shapefile")


# Select demographic data for nests in hex sample
samp.dat <- demo.dat[demo.dat$sampleID %in% nest.samp$sampleID,]
samp.dat.mrg <- merge(samp.dat, nest.samp@data, by="sampleID", all.x=T, all.y=F)
samp.dat.mrg <- samp.dat.mrg[order(c(samp.dat.mrg$year, samp.dat.mrg$sampleID)),]


# Need to summarize histories for hex with >1 nest


# Use data to run occupancy models


# Save all results in big list
results[[j]] <- list(nest.sim=list(params=list(n.nests.mean=n.nests.mean,
                                               n.nests.sd=n.nests.sd,
                                               n.nests.draw=n.nests,
                                               min.dist.mean=min.dist.mean,
                                               min.dist.sd=min.dist.sd,
                                               min.dist.draw=min.dist,
                                               rsf.bins=rsf.bins,
                                               nest.probs=nest.probs,
                                               n.nests.bins=size.samp),
                                   n.drops=table(drops.sp@data$rsf),
                                   nests.per=nests.per,
                                   nest.sample=samp.thin),
                     sample.sim=list(params=list(n.years=n.years,
                                                 # samp.freq=samp.freq,
                                                 n.samples=n.samples,
                                                 hex.area.ha=hex.area,
                                                 hex.cluster=hex.cluster),
                                                 # strat.type=strat.type
                                     n.nests.samp=n.nests.samp,
                                     n.nests.hex=n.nests.hex,
                                     n.nests.cluster=n.nests.cluster,
                                     nest.samp=nest.samp,
                                     hex.samp=hex.samp
                                     ))

} #end simulation loop 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

results

#summary
##parameters
##number of nests per township
##number of nests dropped

# min(nndist(samp.thin$x, samp.thin$y))>min.dist #check min dist


#check distribution of nest datasets in bins
hawk.nests <- readOGR(dsn="E:/Lucretia/Ferruginous_Hawk_Project/RSF/Eagles RSF",layer="active_feha_nests")
hist.hawks <- readOGR(dsn="E:/Lucretia/Ferruginous_Hawk_Project/RSF/Eagles RSF",layer="historic_hawk_nests_w_yng_gte2000")
rand.nests <- readOGR(dsn="E:/Lucretia/Ferruginous_Hawk_Project/RSF/Eagles RSF",layer="random_pts_n1000p")

hawk.nests$bin <- extract(hawk.rsf,hawk.nests)
table(hawk.nests$bin)/nrow(hawk.nests)
hist.hawks$bin <- extract(hawk.rsf,hist.hawks)
table(hist.hawks$bin)/nrow(hist.hawks)
rand.nests$bin <- extract(hawk.rsf,rand.nests)
table(rand.nests$bin)/nrow(rand.nests)












