#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  Ferruginous Hawk and Golden Eagle Monitoring Program Simulations
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(raster)
library(rgdal)
library(spatstat)
library(sp)
library(sf) #st_make_grid
library(maptools)
library(spatialEco)
library(SDraw) #srs functions from SDraw package sample from point objects and retains attributes
library(velox) #faster raster handling
library(rgeos) #gIntersects
library(plyr) #summarize results
library(TSP) #optimal routes between multiple points

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

#directories and paths
setwd("C:/Users/zwallac2/Dropbox (UW WYNDD)/Proj_GovESA_FEHA_GOEA_Monitoring/Data/Analyses/Raptor_Monitoring")
shpOut <- ("C:/Users/zwallac2/Dropbox (UW WYNDD)/Proj_GovESA_FEHA_GOEA_Monitoring/Data/Analyses/Raptor_Monitoring_shpOut")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import and process files

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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set paramaters for simulation
## parms not yet included in code are commented out
## moved all parameter values here before simulation loop

## Simulation rep
n.reps <- 1

## Nests
reoccu.mean <- 0.60 #Average reoccupancy rate from 2011-2013 study
n.nests.mean <- round(sum(hawk.bound@data$area_km2 / c(83,127.1,85.5)) / reoccu.mean, 0) #area of each ecoregion / estimated mean density from Olson et al. 2015

# N nests at yr0 is calculated as:
# (1) the sum of the area of each ecoregion
# (2) divided by the density of occupied nests (km2/nest) per region estimated from distance sampling in 2010-2011 
# (3) Divided by the average re-occupancy rate during the 2011-2013 study
# Step 3 is necessary because the density estimates from  distance sampling are for occupied nests only.
# If we started with the sample size of occupied nests at t0 and multipled that by our chosen re-occupancy rate, we would underestimate the number of nests at t1.
# Instead, we estimate a population size of nests that is 60% occupied at t0.
# Do you agree that this makes sense?

n.nests.sd <- 100 #stdev of n nests
min.dist.mean <- 1500 # minimum distance between nests (m)
min.dist.sd <- 100 #stdev of inter-nest distance

## RSF bins
rsf.bins <- c(1,2,3,4,5,6,7) #rsf raster bins
nest.probs <- c(0,.02,.03,.09,.14,.2,.52) #proportion of nests within each bin from rsf results

## Demographic rates
occupancy.mean <- 0.6
occupancy.sd <- 0 #add inter-annual variation to occupancy
success <- 0.6
# success.sd <- 0.1
young <- 0.7
# young.sd <- 0.2

decline <- 0 #decline over n.years as decimal (e.g. 0.20)

## Detection rates
occ.det <- 0.5 ### This is a tough one to estimate becuase our detection rates are for RE-OCCUPANCY
# occ.prod.det <- 0.80 #detection probability for occupancy with reporductive success
# occ.misclass <- 0.10 #probability of misclassifying a site with reproduction as not reproducing

## Effort
n.years <- 10
# samp.freq <- 2 #Every two years
n.samples <- 200 #number of hexagons sampled
n.vis.occ <- 3
# n.vis.prod <- 2

hex.area <- 1200 #area of hexagons (ha)
hex.cluster <- "N"
rsf.strata <- "N" #stratify by RSF bins?
removal <- "Y" #removal sampling?


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
colnames(hex.grid@data)[1] <- "hexID"

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
set.seed(11)

occu.sim <- function(){

# results <- list()
results <- vector(mode = "list", length = n.reps)

for (j in 1:n.reps){

start_time <- Sys.time()
  
#initial nest sample
#draw parameter values
n.nests <- round(rnorm(1,n.nests.mean,n.nests.sd),0)
min.dist <- rnorm(1,min.dist.mean,min.dist.sd)
size.samp <- round(nest.probs*n.nests,0) #sample size per bin
occupancy <- rnorm(n.years,
                   seq(occupancy.mean, occupancy.mean-decline, length.out=n.years), #occupancy with declining trend
                   occupancy.sd) #sd of occupancy

# # plot
# plot(1:n.years, seq(occupancy.mean, occupancy.mean-decline, length.out=n.years), ylim=c(0,1)) #without variation
# for(i in 1:1000){
#   occupancy <- rnorm(n.years,
#                      seq(occupancy.mean, occupancy.mean-decline, length.out=n.years), #occupancy with declining trend
#                      occupancy.sd) #sd of occupancy
# points(1:n.years, occupancy, col="blue", pch=20, cex=0.5) #with sd
# }

#loop over rsf bins
d.list <- list() #empty list
for (i in 1:length(rsf.bins)){
  d.list[[i]] <- srs.point(hawk.rsf.pts[hawk.rsf.pts@data$rsf==rsf.bins[i],],n=size.samp[i]) #srs.point from SDraw package
  }
d.list[sapply(d.list,is.null)] <- NULL #drop any bins with 0 nests
samp <- do.call(rbind,d.list) #bind into single dataframe
colnames(samp@data)[1] <- "nestID" #rename for clarity
samp@data$nestID <- 1:nrow(samp@data) #renumber


# #check
# length(samp);n.nests
# table(samp$rsf); size.samp
# plot(samp, col=samp@data$rsf)

#for any pair of nests that are closer than the min distance, drop one
#messy code, but works
dmat <- spDists(samp,diagonal=F) #matrix of distances
pairs <- data.frame(which(dmat<=min.dist & dmat>0, arr.in=TRUE)) #pairs of points with <min distance
pairs.unique <- data.frame(unique(t(apply(pairs, 1, sort)))) #unique pairs only
colnames(pairs.unique) <- colnames(pairs) #replace column ids lost with apply
pairs.unique$rsf1 <- merge(pairs.unique,samp@data,by.x="row",by.y="nestID",all.x=TRUE,all.y=FALSE)[,4] #add rsf value for rows
pairs.unique$rsf2 <- merge(pairs.unique,samp@data,by.x="col",by.y="nestID",all.x=TRUE,all.y=FALSE)[,5] #add rsf value for cols
pairs.unique$diff <- pairs.unique$rsf1-pairs.unique$rsf2 #add difference in rsf values

#apply function to choose one of each pair with highest rsf
drops <- apply(pairs.unique, 1, function(x){
  if (x['diff']<0){x['row']
  } else if (x['diff']>0){x['col']
      } else if (x['diff']==0){sample(x[1:2],1)} #choose randomly if rsf values of pair are equal
})

samp.thin <- samp[!(samp$nestID %in% drops),] #thinned nest sample spatial data file
drops.sp <- samp[(samp$nestID %in% drops),] #dropped nests spatial data file

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
  nestID <- sim.nests$nestID
  n.nests <- nrow(sim.nests)
  occ <- rbinom(n.nests,1,occupancy[i])
  succ <- rbinom(n.nests,1,success)*occ
  yng <- rpois(n.nests,1)*succ

  demo.list[[i]] <- cbind(nestID,year,occ,succ,yng)
}
demo.dat <- data.frame(do.call(rbind,demo.list))
# Summarize results
# ddply(demo.dat, .(year), summarise, meanOccuYear=mean(occ), nOccuYear=sum(occ))


# Sample hexagons

if(hex.cluster=="Y"){
  hex.hip <- hip.polygon(hex.grid, n.samples/7) #divide n.samples by 7 to match final sample size with 'flower' clusters
  colnames(hex.hip@data)[1] <- "hexSampleOrder" #rename master sample order column
  hex.samp <- hex.grid[hex.grid@data$hexID %in% hex.hip@data$hexID,] #subset hex with sample points
  mat <- gIntersects(hex.samp, hex.grid, byid=T) #matrix of neighboring hexagons
  pairs <- data.frame(which(mat==TRUE, arr.in=TRUE)) #collapse to table of pairs
  sampHex <- data.frame(col=1:ncol(mat),
                        sampHex=colnames(mat))
  pairs <- pairs[order(pairs$row),] #reorder
  hex.clust <- hex.grid[hex.grid@data$hexID %in% pairs$row,] #select all hexagons in sample
  hex.clust@data$clusterID <- pairs$col #add cluster id
  hex.clust$type <- rep("cluster",nrow(hex.clust@data)) #add column indicating if hex is first of second stage cluster sample
  hex.clust$type[hex.clust$hexID %in% hex.samp$hexID] <- "sample" #indicate first stage sample hexagons
  hex.samp <- hex.clust #assign final sample to same name as unclustered method above
} else{
  hex.hip <- hip.polygon(hex.grid, n.samples) #2.92 min. Faster than hip.point on centroids for some reason
  colnames(hex.hip@data)[1] <- "hexSampleOrder" #rename master sample order column
  hex.samp <- hex.grid[hex.grid@data$hexID %in% hex.hip@data$hexID,] #subset hex with sample points
}


# #plot
# plot(hawk.bound, col="grey90")
# plot(hex.samp, border=hex.samp$clusterID, add=T)

# #cluster sample is slightly smaller because some cluster cells fall outside the study area
# #could buffer study area or draw cluster centroids only from cells with 6 neighbors
# hex.mat <- gIntersects(hex.grid, hex.grid, byid=T) #matrix of neighboring hexagons
# hex.pairs <- data.frame(which(hex.mat==TRUE, arr.in=TRUE)) #collapse to table of pairs
# n.neighbs <- ddply(hex.pairs, .(row), summarize, n=length(col))
# keeps <- n.neighbs$row[n.neighbs$n==7]
# tmp <- hex.grid[hex.grid$hexID %in% keeps,]
# tmp2 <- hex.grid[!(hex.grid$hexID %in% keeps),]
# plot(tmp); plot(tmp2, add=T, col="cornflowerblue")


# rsf.strata =="Y"
#loop over rsf bins
#draw hip sample from each bin with n.samples proportional to inclusion prob
#join resulting dataframes
#add it if-else with clusters



# Calculate shortest route between hexagons
# "Traveling Salesman Problem"
hex.samp.cents <- SpatialPointsDataFrame(gCentroid(hex.samp, byid=T), hex.samp@data) #centroids of hex sample
dMat <- spDists(hex.samp.cents) #distance matrix 
tsp <- TSP(dMat,labels=hex.samp.cents@data$hexID,method='euclidean') #create TSP object
solution <- solve_TSP(tsp) #solve using TSP methods
route.length <- tour_length(solution)*0.001 #tour length in km


#nests in hexagons
nest.hex <- over(samp.thin, hex.samp)
nest.samp.tab <- nest.hex[!is.na(nest.hex$hexID),]
nest.samp.tab$geometryID <- rownames(nest.samp.tab)
nest.samp <- samp.thin[samp.thin$geometryID %in% nest.samp.tab$geometryID,]
nest.samp <- merge(nest.samp, nest.samp.tab, by="geometryID", all.x=T, all.y=F)

# plot(hex.samp)
# plot(nest.samp,col="blue", add=T)

## summarize
n.nests.samp <- nrow(nest.samp) #N nests in sample
n.nests.hex <- table(as.numeric(table(nest.hex$hexID))) #table of N nests per hex
n.nests.cluster <- table(as.numeric(table(nest.hex$cluster))) #table of N nests per hex
## save
# writeOGR(samp.thin, shpOut, "nestSample", driver="ESRI Shapefile")
# writeOGR(hex.samp, shpOut, "hexSample", driver="ESRI Shapefile")


# Select demographic data for nests in hex sample
samp.dat <- demo.dat[demo.dat$nestID %in% nest.samp$nestID,] #nests in hex sample
samp.dat.mrg <- merge(samp.dat, nest.samp@data, by="nestID", all.x=T, all.y=T) #add hex data to nest table
samp.dat.mrg <- samp.dat.mrg[order(samp.dat.mrg$year, samp.dat.mrg$nestID),] #order by year and sample id

# Make table of true occupancy for all hexagons in sample
hex.nest <- ddply(samp.dat.mrg, .(hexID,year), summarise, occ=sum(occ)) #hexagons with nests
hex.noNest <- hex.samp$hexID[!(hex.samp$hexID %in% hex.nest$hexID)] #hexagons without nests
hex.noNest.df <- data.frame(hexID = rep(hex.noNest, each=n.years),
                            year=rep(1:n.years, length(hex.noNest)),
                            occ=0)
hex.all <- rbind(hex.nest,hex.noNest.df) #combined
hex.all <- hex.all[order(hex.all$year,hex.all$hexID),]

#Mean true occupancy probability
occu.truth <- hex.all
occu.truth$occ[occu.truth$occ>1] <- 1
# ddply(hex.all, .(year), summarise, meanOccu=mean(occ))


#Simulate detection process for all hexagons in sample
# For hexagons with multiple occupied nests, I used the cumulative detection probability (p*), 
# based on the idea that the probability of detecting an occupied hexagon was the cumulative probability of detecting at least 1 occupied nest
# we could use the same approach for occupancy with reproductive success

# Simulate detection
det.hist <- data.frame(matrix(NA,nrow=nrow(hex.all),ncol=n.vis.occ))
  for(j in seq_along(1:n.vis.occ)){ #loop over number of occupancy visits
    det.hist[,j] <- rbinom(n.samples*n.years,1,
                           (1-((1-occ.det)^hex.all$occ))*occu.truth$occ) #p* for hexagons with multiple occupied nests
    }

#replace data following detections with NAs for removal design
if (removal=="Y"){
for(k in 1:(ncol(det.hist)-1)){
  det.hist[,(k+1)][det.hist[,k]>0 | is.na(det.hist[,k])] <- NA
}
}
colnames(det.hist) <- paste0(rep("vis",n.vis.occ),1:n.vis.occ)

# Unique histories
# table(apply(det.hist[,1:ncol(det.hist)],1,paste,collapse=","))


# Use data to run occupancy models

end_time <- Sys.time()
run.time <- end_time - start_time

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
                                                 hex.cluster=hex.cluster,
                                                 rsf.strata=rsf.strata),
                                     route.length=route.length,
                                     n.nests.samp=n.nests.samp,
                                     n.nests.hex=n.nests.hex,
                                     n.nests.cluster=n.nests.cluster,
                                     nest.samp=nest.samp,
                                     hex.samp=hex.samp),
                     occu.sim=list(occu.truth=occu.truth,
                                   removal=removal),
                     run.time=list(run.time=run.time)
                     )

} #end simulation loop

return(results)

} #end function

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

run1 <- occu.sim()

results

#summary

# figure out how to lapply to summarize key results

##parameters
##number of nests per township
##number of nests dropped

# min(nndist(samp.thin$x, samp.thin$y))>min.dist #check min dist



## It would be helpful to see these results before we go much further
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
