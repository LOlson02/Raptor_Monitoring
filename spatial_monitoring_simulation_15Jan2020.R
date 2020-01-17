#####################################################################################################################################################################
#
#  Ferruginous Hawk and Golden Eagle Monitoring Program Simulations
#
#####################################################################################################################################################################

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
library(viridis) #accessible color palettes

options(scipen=999)

# #setwd
# setwd("E:/Lucretia/Ferruginous_Hawk_Project/Monitoring")

#--- Zach's computer
#directories and paths
setwd("C:/Users/zwallac2/Dropbox (UW WYNDD)/Proj_GovESA_FEHA_GOEA_Monitoring/Data/Analyses/Raptor_Monitoring")
shpOut <- ("C:/Users/zwallac2/Dropbox (UW WYNDD)/Proj_GovESA_FEHA_GOEA_Monitoring/Data/Analyses/Raptor_Monitoring_shpOut")
#--- End Zach's computer

###############################################################################################################
#
#  Import and process files
#
###############################################################################################################

# source('E:/Lucretia/Ferruginous_Hawk_Project/Monitoring/Raptor_Monitoring/mode_func.R')
# source('E:/Lucretia/Ferruginous_Hawk_Project/Monitoring/Raptor_Monitoring/hex_grid_func.R')
#
# #import species RSFs
# hawk.rsf <- raster("largehawksGLM3_7bins.tif")
# 
# #import townships
# townships <- readOGR(dsn="E:/Lucretia/Ferruginous_Hawk_Project/WY_background_layers",layer="townships_all_possible")
# #import boundaries
# state <- readOGR(dsn="E:/Lucretia/Ferruginous_Hawk_Project/WY_background_layers",layer="State_Boundary")
# study.bound <- readOGR(dsn="E:/Lucretia/Ferruginous_Hawk_Project/WY_background_layers/Wyo_ecoregions_EPA",layer="study_areas")

#--- Zach's computer
source('C:/Users/zwallac2/Dropbox (UW WYNDD)/Proj_GovESA_FEHA_GOEA_Monitoring/Data/Analyses/Raptor_Monitoring/mode_func.R') #mode function for zonal stats

# Import species RSFs
hawk.rsf <- raster("C:/Users/zwallac2/Dropbox (UW WYNDD)/Proj_GovESA_FEHA_GOEA_Monitoring/GIS/RSF/largehawksGLM3_7bins.tif")

# Import sample strata RSF
strata.rsf <- hawk.rsf

# Import townships
townships <- readOGR(dsn="C:/Users/zwallac2/Dropbox (UW WYNDD)/Proj_BI_USFS_FEHA_Movement/GIS/Base_Data/townships",layer="townships_all_possible")

# Import boundaries
state <- readOGR(dsn="C:/Users/zwallac2/Dropbox (UW WYNDD)/Proj_GovESA_FEHA_GOEA_Monitoring/GIS/Base_Data/Administrative",layer="Wyoming_Boundary")
study.bound <- spTransform(readOGR(dsn="C:/Users/zwallac2/Dropbox (UW WYNDD)/Proj_GovESA_FEHA_GOEA_Monitoring/GIS/Base_Data/Study_area_shapefile",layer="study_areas"),proj4string(state))
#--- End Zach's computer

#check crs
identicalCRS(hawk.rsf, townships)
identicalCRS(state, townships)
identicalCRS(study.bound, townships)

#plot
# plot(hawk.rsf)
# plot(state, add=T)
# plot(study.bound, add=T)
# plot(townships, add=T)

#convert raster to points
hawk.rsf.pts <- rasterToPoints(hawk.rsf, spatial=TRUE)
colnames(hawk.rsf.pts@data) <- "rsf" #rename rsf bins column
hawk.rsf.pts <- hawk.rsf.pts[complete.cases(hawk.rsf.pts@data),] #drop NA rows, if necessary
#not necessary to clip points by study area polygon, but could be for goea
saveRDS(hawk.rsf.pts, file="hawkRsfPts") #just in case because rasterToPoints is slow
# hawk.rsf.pts <- readRDS("hawkRsfPts")


###############################################################################################################
#
#  HYPER PARAMETERS
#
###############################################################################################################

## Simulation reps
n.reps <- 50

## Effort
n.years <- 20


###############################################################################################################
#
#  POPULATION SIMULATIONS
#
###############################################################################################################

# Set parameters

## Demographic rates
occupancy.mean <- 0.6 #average reoccupancy rate from 2011-2013 study
occupancy.sd <- 0.02 #add inter-annual variation to occupancy
success <- 0.6
# success.sd <- 0.1
young <- 0.7
# young.sd <- 0.2

## Decline
decline <- 0.10 #decline over n.years as decimal (e.g. 0.20)

## RSF bins
rsf.bins <- 1:7 #rsf raster bins
nest.probs <- c(0,0.02,0.03,0.09,0.14,0.2,0.52) #proportion of nests within each bin from rsf results

## Number of nesting territories in population
n.nests.mean <- round(sum(study.bound@data$area_km2 / c(83,127.1,85.5)) / occupancy.mean, 0) #area of each ecoregion / estimated mean density from Olson et al. 2015
# N nests at yr0 is calculated as:
# (1) the sum of the area of each ecoregion
# (2) divided by the density of occupied nests (km2/nest) per region estimated from distance sampling in 2010-2011 
# (3) Divided by the average re-occupancy rate during the 2011-2013 study
# Step 3 is necessary because the density estimates from  distance sampling are for occupied nests only.
# If we started with the sample size of occupied nests at t0 and multipled that by our chosen re-occupancy rate, we would underestimate the number of nests at t1.
# Instead, we estimate a population size of nests that is 60% occupied at t0.
# Do you agree that this makes sense?
n.nests.sd <- 100 #stdev of n nests

## Inter-nest distance
min.dist.mean <- 1500 # minimum distance between nests (m)
min.dist.sd <- 100 #stdev of inter-nest distance


# Function to simulate population of nesting territories
pop.sim <- function(){

  pop.sim.out <- list()
  
  for (j in 1:n.reps){
    start_time <- Sys.time()
    
    # Draw parameter values for initial population
    n.nests <- round(rnorm(1,n.nests.mean,n.nests.sd),0)
    min.dist <- rnorm(1,min.dist.mean,min.dist.sd)
    size.samp <- round(nest.probs*n.nests,0) #sample size per bin
    occupancy <- rnorm(n.years,
                       seq(occupancy.mean, occupancy.mean-decline, length.out=n.years), #occupancy with declining trend
                       occupancy.sd) #sd of occupancy
    
    # Loop over rsf bins to allocate nests
    d.list <- list() #empty list
    for (i in 1:length(rsf.bins)){
      d.list[[i]] <- srs.point(hawk.rsf.pts[hawk.rsf.pts@data$rsf==rsf.bins[i],],n=size.samp[i]) #srs.point from SDraw package
    }
    d.list[sapply(d.list,is.null)] <- NULL #drop any bins with 0 nests
    samp <- do.call(rbind,d.list) #bind into single dataframe
    colnames(samp@data)[1] <- "nestID" #rename
    samp@data$nestID <- 1:nrow(samp@data) #renumber
    
    # #check
    # length(samp);n.nests
    # table(samp$rsf); size.samp
    # plot(samp, col=samp@data$rsf)

    # Drop one of any pair of nests that are closer than the min distance
    ## messy code, but works
    dmat <- spDists(samp,diagonal=F) #matrix of distances
    pairs <- data.frame(which(dmat<=min.dist & dmat>0, arr.in=TRUE)) #pairs of points with <min distance
    pairs.unique <- data.frame(unique(t(apply(pairs, 1, sort)))) #unique pairs only
    colnames(pairs.unique) <- colnames(pairs) #replace column ids lost with apply
    pairs.unique$rsf1 <- merge(pairs.unique,samp@data,by.x="row",by.y="nestID",all.x=TRUE,all.y=FALSE)[,4] #add rsf value for rows
    pairs.unique$rsf2 <- merge(pairs.unique,samp@data,by.x="col",by.y="nestID",all.x=TRUE,all.y=FALSE)[,5] #add rsf value for cols
    pairs.unique$diff <- pairs.unique$rsf1-pairs.unique$rsf2 #add difference in rsf values
    
    ## Apply function to choose one of each pair with highest rsf
    drops <- apply(pairs.unique, 1, function(x){
      if (x['diff']<0){x['row']
      } else if (x['diff']>0){x['col']
      } else if (x['diff']==0){sample(x[1:2],1)} #choose randomly if rsf values of pair are equal
    })
    
    samp.thin <- samp[!(samp$nestID %in% drops),] #thinned nest sample spatial data file
    # drops.sp <- samp[(samp$nestID %in% drops),] #dropped nests spatial data file
    
    # Check number of nests per township
    nests.twp <- over(samp.thin,townships)
    nests.per.twp <- table(as.numeric(table(nests.twp$joinfield)))
    
    # This approach works, but does not draw more nests to replace those dropped. This decreases the sample size  by ~3% for 1500 nests
    # We could consider adding code to replace nests and match the sample size exactly
    
    # Demographic data for simulated nests
    sim.nests <- samp.thin@data #make table
    sim.nests$X <- samp.thin$x #add coordinates
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

    end_time <- Sys.time()
    run.time <- end_time - start_time

# Output list        
pop.sim.out[[j]] <- list(params=list(n.reps=n.reps,
                                     n.years=n.years,
                                     
                                     occupancy.mean=occupancy.mean,
                                     occupancy.sd=occupancy.sd,
                                     occupancy.draw=occupancy,
                                     success=success,
                                     young=young,
                                     
                                     decline=decline,
                                     
                                     rsf.bins=rsf.bins,
                                     nest.probs=nest.probs,
                                     
                                     n.nests.mean=n.nests.mean,
                                     n.nests.sd=n.nests.sd,
                                     n.nests.draw=n.nests,
                                     min.dist.mean=min.dist.mean,
                                     min.dist.sd=min.dist.sd,
                                     min.dist.draw=min.dist,
                                     n.nests.bins=size.samp),
                         
                         n.drops=table(drops.sp@data$rsf),
                         nests.per.twp=nests.per.twp,
                         
                         sim.pop=sim.nests,
                         samp.thin.sp=samp.thin,
                         run.time=run.time)
  } #end simulation loop
  
  return(pop.sim.out)
  
} #end function
pop.sims <- pop.sim() #call

# Summarize results
## Total run time
sum(sapply(pop.sims, function(x){x$run.time}))/60 #minutes run time -- 6.13 min for 50 reps

# # Stochastic variables
# hist(sapply(pop.sims, function(x){x$params$n.nests.draw}),xlab="Nests (N)",main="N nests in population")
# hist(sapply(pop.sims, function(x){x$params$min.dist.draw}),xlab="Distance (m)",main="Minimum inter-nest distance")
# rowMeans(sapply(pop.sims,function(x){x$params$occupancy.draw}))
## Occupancy plot
plot(1:n.years, seq(occupancy.mean, occupancy.mean-decline, length.out=n.years), ylim=c(0,1)) #without variation
points(rep(1:n.years,n.reps),
       unlist(lapply(pop.sims,function(x){x$params$occupancy.draw})),
       col="blue", pch=20, cex=0.5) #with sd
# N nests per twp
# annoying to average ragged tables - consider summarizing in a different way



###############################################################################################################
#
#  HEX SAMPLING
#
###############################################################################################################

# Set parameters
n.samples <- 200 #number of hexagons sampled
hex.area <- 1200 #area of hexagons (ha)
hex.cluster <- "Y"
rsf.strata <- "Y" #stratify by RSF bins?
strata.rsf.bins <- 1:7 #in this case we are using the same rsf that we used for the population simulations, but that will not be the case once we combine hawks and eagles
inc.probs <- 1:7/sum(1:7)  #inclusion probabilities proportional to rsf bins, for now -- increasing order, must sum to 1


# Function to make hexagon grid
make.hex.grid <- function(study.bound,strata.rsf,proj){

  grid <- st_make_grid(st_as_sf(study.bound), #convert to sf type
                       cellsize = 2*sqrt(hex.area*1e+4/(2*sqrt(3))), #calculate short axis (m) from area (ha)
                       square = FALSE)
  ##clean up attributes
  hex.grid <- spTransform(as(grid, 'Spatial'), proj) #convert back to st and reproject
  hex.ids <- sapply(slot(hex.grid, "polygons"), function(x) slot(x, "ID")) #get IDs from hex.grid
  hex.df <- data.frame(ID=1:length(hex.grid), row.names=hex.ids) #make data frame of IDs
  hex.grid <- SpatialPolygonsDataFrame(hex.grid, hex.df) #join dataframe to spatial polygons
  colnames(hex.grid@data)[1] <- "hexID"
  
  ##summarize rsf within hex grid -- this will be used for stratified sampling
  #get average (mode) RSF value in each hex
  rsf.v <- velox(strata.rsf)
  hex.modes <- rsf.v$extract(sp=hex.grid,fun=getmode) #use getmode function ##added na.omit to mode_func.R to work for hexagons that partially overlap rsf
  hex.grid@data$rsfMode <- hex.modes[,1]
  saveRDS(hex.grid, file="hexGrid") #save
  writeOGR(hex.grid, shpOut, "hexGridTest", driver="ESRI Shapefile", overwrite_layer=TRUE)
  
  return(hex.grid)
}
  start_time <- Sys.time()
hex.grid <- make.hex.grid(study.bound,strata.rsf,proj4string(townships)) #call
  end_time <- Sys.time(); end_time - start_time #<1 min

hex.cents <- SpatialPointsDataFrame(gCentroid(hex.grid,byid=T), hex.grid@data) #hexagon centroids
n.hex <- length(hex.grid) #save hex grid sample size


# Function to sample hex grid
hex.sim <- function(pop.sims){

    hex.sim.out <- list()
    
  # Adjust sample size if clusters are used to select only centroids
  if(hex.cluster=="Y"){
    tot.samp = n.samples
    n.samp = n.samples/7}else{
      n.samp = n.samples}
  
  # Loop over number of simulation replicates
  for (j in 1:n.reps){
    
    start_time <- Sys.time()
  
    #Stratified sample by RSF
    if(rsf.strata=="Y"){
      hex.hip.rsf <- list() #list to store results
      for(b in 1:length(rsf.bins)){ #loop over rsf bins
        strat <- hex.cents[hex.cents$rsfMode==b & !is.na(hex.cents$rsfMode),] #subset by bin/stratum
        hex.hip <- hip.point(strat, round(n.samp*inc.probs[b],0)) #sample from stratum proportional to inclusion probability
        colnames(hex.hip@data)[1] <- "hexSampleOrder" #rename master sample order column
        hex.hip@data$rsfStratum <- b #add column indicating stratum
        hex.hip.rsf[[b]] <- hex.hip #add to list
      }
      hex.hip <- do.call(rbind,hex.hip.rsf) #bind list into one spdf
      hex.hip$sampleID <- 1:n.samp
    
      # Select hexagon sample using hip point and bind necessary data columns
      df <- cbind(over(hex.hip, hex.grid), hex.hip@data)
      df <- subset(df, select=-c(rsfMode)) #drop column to avoid duplicates later
      h <- merge(hex.grid, df, by="hexID")
      h$clusterID <- NA #add columns to match cluster sample df
      h$type <- "sample"
      hex.samp <- h[!is.na(h$sampleID),]
      hex.samp <- hex.samp[order(hex.samp$sampleID),]
    
    #Un-stratified sample  
    }else{
      hex.hip <- hip.point(hex.cents, n.samp)
      colnames(hex.hip@data)[1] <- "hexSampleOrder" #rename master sample order column
      hex.hip@data$rsfStratum <- NA
      hex.hip$sampleID <- 1:n.samp #redundant, but makes columns consistent with stratified sample df
      
      # Select hexagon sample using hip point and bind necessary data columns
      df <- cbind(over(hex.hip, hex.grid), hex.hip@data)
      df <- subset(df, select=-c(rsfMode)) #drop column to avoid duplicates later
      h <- merge(hex.grid, df, by="hexID")
      h$clusterID <- NA #add columns to match cluster sample df
      h$type <- "sample"
      hex.samp <- h[!is.na(h$sampleID),]
      hex.samp <- hex.samp[order(hex.samp$sampleID),]
      }
    
    #Cluster sample
    if(hex.cluster=="Y"){
      dMat <- gIntersects(hex.samp, hex.grid, byid=T) #matrix of neighboring hexagons
      pairs <- data.frame(which(dMat==TRUE, arr.in=TRUE)) #collapse to table of pairs
      sampHex <- data.frame(col=1:ncol(dMat),
                            sampHex=colnames(dMat))
      pairs <- pairs[order(pairs$row),] #reorder
      colnames(pairs) <- c("hexID", "clusterID")
      pairs <- pairs[!duplicated(pairs$hexID),] #drop overlapping hexagons from second cluster
      hex.clust <- hex.grid[hex.grid@data$hexID %in% pairs$hexID,] #select all hexagons in sample
      hex.clust@data$clusterID <- pairs$clusterID
      hex.clust$type <- rep("cluster",nrow(hex.clust@data)) #add column indicating if hex is first of second stage cluster sample
      hex.clust$type[hex.clust$hexID %in% hex.samp$hexID] <- "sample" #indicate first stage sample hexagons
      h <- merge(hex.clust, hex.samp@data[,c(1,3:5)], by="hexID") #merge, dropping duplicate columns
      h@data <- h@data[colnames(hex.samp@data)] #order columns to match unclustered df
      #assign centorids of overlapping hexagons back to their cluster 
      h@data$clusterID[!is.na(h@data$sampleID) & !(h@data$sampleID==h@data$clusterID)] <- h@data$sampleID[!is.na(h@data$sampleID) & !(h@data$sampleID==h@data$clusterID)]
      h <- h[order(h$clusterID,h$hexSampleOrder),]
      hex.samp <- h #assign final sample to same name as unclustered method above
    }
    
    # #plot unstratified
    # plot(study.bound, border="grey40", col="grey90")
    # plot(hex.samp, border='black', col="cornflowerblue", add=T)
    # 
    # #plot stratified, unclustered
    # v.pal <- rep(viridis(length(rsf.bins)), table(hex.samp@data$rsfStratum))
    # plot(study.bound, border="grey40", col="grey90")
    # plot(hex.samp, col=v.pal, add=T)
    # 
    # #plot stratified, clustered
    # v.pal <- rep(viridis(length(rsf.bins)), table(hex.samp@data$rsfStratum[hex.samp@data$type=="sample"]))
    # plot(study.bound, border="grey40", col="grey90")
    # plot(hex.samp[hex.samp@data$type=="sample",], col=v.pal, add=T)
    # plot(hex.samp[hex.samp@data$type=="cluster",], col="white", add=T)
    # 
    # # cluster codes correctly assigned?
    # plot(hex.samp, col=hex.samp@data$clusterID)
    
    # #cluster sample is slightly smaller because some cluster cells fall outside the study area
    # length(hex.samp); tot.samp
    
    
    # Select nests in sampled hexagons
    nests.in.hex <- function(){
      pop.sim <- pop.sims[[j]]$samp.thin.sp
      nest.hex <- over(pop.sim, hex.samp)
      nest.samp.tab <- nest.hex[!is.na(nest.hex$hexID),]
      nest.samp.tab$geometryID <- rownames(nest.samp.tab)
      nest.samp <- pop.sim[pop.sim@data$geometryID %in% nest.samp.tab$geometryID,]
      nest.samp <- merge(nest.samp, nest.samp.tab, by="geometryID", all.x=T, all.y=F)
      return(nest.samp)
    }
    nest.samp <- nests.in.hex() #call
    
    
    end_time <- Sys.time()
    run.time <- end_time - start_time
    
    # Output list        
    hex.sim.out[[j]] <- list(params=list(n.hex.frame=n.hex,
                                         n.samples=tot.samp,
                                         hex.area.ha=hex.area,
                                         hex.cluster=hex.cluster,
                                         rsf.strata=rsf.strata,
                                         strata.rsf.bins=strata.rsf.bins,
                                         inc.probs=inc.probs),
                             hex.sample=list(hex.n.target=tot.samp,
                                             hex.n.realized=length(hex.samp),
                                             hex.samp=hex.samp),
                             nest.sample=list(n.nests.samp=nrow(nest.samp),
                                              n.nests.hex=table(as.numeric(table(nest.hex$hexID))),
                                              n.nests.cluster=table(as.numeric(table(nest.hex$cluster))),
                                              nest.samp=nest.samp),
                             run.time=run.time)

    ## save
    # writeOGR(samp.thin, shpOut, "nestSample", driver="ESRI Shapefile")
    # writeOGR(hex.samp, shpOut, "hexSample", driver="ESRI Shapefile")
  
  } #end simulation loop
    
    return(hex.sim.out)

} #end function
hex.sims <- hex.sim(pop.sims) #call

# Summarize results
## Total run time
sum(sapply(hex.sims, function(x){x$run.time}))/60 #minutes run time -- 6.13 min for 50 reps = 20 hrs for 10,000 reps


###############################################################################################################
#
#  OCCUPANCY SURVEY AND MODELS
#
###############################################################################################################


#########################################
#### WORKING ON THIS SECTION ############
#########################################


# Set parameters

## Occupancy survey design
n.vis.occ <- 3
# n.vis.prod <- 2
removal <- "Y" #removal sampling?
samp.freq <- 1 #Number of years between sampling (1 = annual)

## Detection rates
occ.det <- 0.5 ### This is a tough one to estimate becuase our detection rates are for RE-OCCUPANCY
# occ.prod.det <- 0.80 #detection probability for occupancy with reporductive success
# occ.misclass <- 0.10 #probability of misclassifying a site with reproduction as not reproducing


# Function to simulate occupancy survey
occu.sim <- function(){
  
  occu.sim.out <- list()
  
  # Loop over number of simulation replicates
  for (j in 1:n.reps){
    
  start_time <- Sys.time()

  # Select demographic data for nests in hex sample
  demo.dat <- pop.sims[[j]]$
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
    det.hist[,j] <- rbinom(n.samp*n.years,1,
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
  occu.sim.out[[j]] <- list(nest.sim=list(params=list(n.nests.mean=n.nests.mean,
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
                                                   n.samples=n.samp,
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
                       run.time=list(run.time=run.time))

  } #end simulation loop

  return(occu.sim.out)

} #end function
occu.sims <- occu.sim() #call



###############################################################################################################
#
#  EFFORT AND COST
#
###############################################################################################################

# Set parameters
min.km2 <- 0.63 # Navajo Nation hexagon survey time: 0.63 min/km2
min.hex <- min.km2 * (hex.area*0.01) #7.56 min/hexagon
ferry.speed <- 200 #km/hr
cost.hr <- 400 #dollars -- rough estimate for plane, pilot, hotels, observer


# Function to simulate effort and cost
effort.sim <- function(hex.sims){
  
  effort.sim.out <- list()
  
  # Loop over number of simulation replicates
  for (j in 1:n.reps){
    
    start_time <- Sys.time()

    # Calculate shortest route between hexagons
    ## "Traveling Salesman Problem"
    traveling.salesman <- function(){
        hex.samp <- hex.sims[[j]]$hex.sample$hex.samp   ###Currently set to calculate route lenght for complete sample (occasion 1) -- need to use occupancy data to calculate effort for multiple occasions with removal design
        hex.samp.cents <- SpatialPointsDataFrame(gCentroid(hex.samp, byid=T), hex.samp@data) #centroids of hex sample
        dMat <- spDists(hex.samp.cents) #distance matrix 
        tsp <- TSP(dMat,labels=hex.samp.cents@data$hexID,method='euclidean') #create TSP object
        solution <- solve_TSP(tsp) #solve using TSP methods
        route.length <- tour_length(solution)*0.001 #tour length in km
        
        #Plotting code -- comment out for simulation runs
        path <- cut_tour(solution, paste0(hex.samp.cents@data$hexID[1])) #cut circular tour in arbitrary location
        path.order <- data.frame(hexID=as.numeric(names(path)), #make df
                                 rowID=path,
                                 flyOrder=1:length(path))
        tmp <- merge(hex.samp.cents, path.order, by="hexID") #merge with hex.cents
        tmp$line <- "a" #id to group all points into single polyline
        tmp <- tmp[order(tmp@data$flyOrder),] #sort by fly order
        
        x <- lapply(split(tmp, tmp$line), function(x) Lines(list(Line(coordinates(x))), x$line[1L]))
        lines <- SpatialLines(x)
        data <- data.frame(id = unique(tmp$line))
        rownames(data) <- data$id
        l <- SpatialLinesDataFrame(lines, data)
        proj4string(l) <- proj4string(tmp)
        plot <- plot(l, lty=2, col="coral"); plot(tmp, pch=20, add=T)
        return(list(route.length, plot))
        
        # return(route.length)
      }
    route.length <- traveling.salesman() #call
    
    # Estimate time and cost of route
    # Make this a function, too
    n.hex.samples <- hex.sims[[j]]$hex.sample$hex.n.realized
    hrs.round.1 <- (route.length[[1]] / ferry.speed) + (n.hex.samples * min.hex / 60) #Hours to complete one full round of survey
    cost.round.1 <- hrs.round.1 * cost.hr


        
    end_time <- Sys.time()
    run.time <- end_time - start_time
    
    # Output list        
    effort.sim.out[[j]] <- list(params=list(min.km2=min.km2,
                                            min.hex=min.hex,
                                            ferry.speed=ferry.speed,
                                            cost.hr=cost.hr),
                                hex.sample=list(n.hex.samples=n.hex.samples),
                                effort.results=list(route.length=route.length[[1]],
                                                 hrs.round.1=hrs.round.1,
                                                 cost.round.1=cost.round.1),
                             run.time=run.time)
    
  } #end simulation loop
  
  return(effort.sim.out)
  
} #end function
effort.sims <- effort.sim(hex.sims) #call



###############################################################################################################
#
#  RESULTS SUMMARY
#
###############################################################################################################

# Summarize results here



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Lucretia's old code

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
