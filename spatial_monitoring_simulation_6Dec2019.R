library(raster)
library(rgdal)
library(spatstat)
library(sp)
library(maptools)
library(spatialEco)
library(SDraw) #srs functions from SDraw package sample from point objects and retains attributes


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


#####Zach's computer
source('C:/Users/zwallac2/Dropbox (UW WYNDD)/Proj_GovESA_FEHA_GOEA_Monitoring/Data/Analyses/Raptor_Monitoring/mode_func.R')
source('C:/Users/zwallac2/Dropbox (UW WYNDD)/Proj_GovESA_FEHA_GOEA_Monitoring/Data/Analyses/Raptor_Monitoring/hex_grid_func.R')

#setwd
setwd("C:/Users/zwallac2/Dropbox (UW WYNDD)/Proj_GovESA_FEHA_GOEA_Monitoring/Data/Analyses/Raptor_Monitoring")

#import species RSFs
hawk.rsf <- raster("C:/Users/zwallac2/Dropbox (UW WYNDD)/Proj_GovESA_FEHA_GOEA_Monitoring/GIS/RSF/largehawksGLM3_7bins.tif")

#import townships
townships <- readOGR(dsn="C:/Users/zwallac2/Dropbox (UW WYNDD)/Proj_BI_USFS_FEHA_Movement/GIS/Base_Data/townships",layer="townships_all_possible")
#import boundaries
state <- readOGR(dsn="C:/Users/zwallac2/Dropbox (UW WYNDD)/Proj_GovESA_FEHA_GOEA_Monitoring/GIS/Base_Data/Administrative",layer="Wyoming_Boundary")
hawk.bound <- spTransform(readOGR(dsn="C:/Users/zwallac2/Dropbox (UW WYNDD)/Proj_GovESA_FEHA_GOEA_Monitoring/GIS/Base_Data/Study_area_shapefile",layer="study_areas"),proj4string(state))
#####

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


#set paramaters for simulation
n.nests <- 500 #number of nests
rsf.bins <- c(1,2,3,4,5,6,7) 
nest.probs <- c(0,.02,.03,.09,.14,.2,.52) #Proportion of nests within each RSF bin
min.dist <- 1500 # minimum distance between nests (m)

#Draw sample stratified by rsf bins
#sample sizes
size.samp <- round(nest.probs*n.nests,0) #sample size per bin


nDrops <- list()
for (j in 1:100){

#loop over rsf bins
d.list <- list() #empty list
for (i in 1:length(rsf.bins)){
  d.list[[i]] <- srs.point(hawk.rsf.pts[hawk.rsf.pts@data$rsf==rsf.bins[i],],n=size.list[i])
}
d.list[sapply(d.list,is.null)] <- NULL #drop any binds with 0 nests
samp <- do.call(rbind,d.list) #bind into single dataframe
samp@data$sampleID <- 1:nrow(samp@data) #renumber

#check
# length(samp);n.nests
# table(samp$rsf); size.list
# plot(samp, col=samp$rsf)

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

samp.thin <- samp[!(samp$sampleID %in% drops),] #thinned nest sample
drops.sp <- samp[(samp$sampleID %in% drops),] #dropped nests

length(drops) #how many nests dropped
data.frame(table(drops.sp@data$rsf)) #how many nests dropped from each rsf bin

nDrops[[j]] <- length(drops)
}

nDropped <- unlist(nDrops)
summary(nDropped/500) # Dropping 3.6% of 1500 nests, 10.9% of 500 nests


# #create random points inside hawk boundary
# # samp <- spsample(hawk.bound,n=10000,type="random")
# samp <- spsample(hawk.bound,n=5000,type="random") #Don't have enough memory to do matrix operation, so reduced dim
# samp$bin <- extract(hawk.rsf,samp)
# samp <- samp[complete.cases(samp@data$bin),]

#calculate distance matrix for all random points
dmat <- spDists(samp,diagonal=F)
min.dist <- 1500 #set minimum distance between nests (m)
dmat[dmat<=min.dist & dmat>0]<- NA  #make anything less than the minimum distance = NA

#for any pair of nests that are closer than the min distance, drop one
rownames(dmat) <- 1:nrow(dmat)
dmat2 <- na.omit(dmat) ###ZW: Wouldn't this drop both nests in pairs that are closer than min distance?
samp2 <- samp[as.numeric(rownames(dmat2)),]

#loop through each RSF bin and pick the specified number of nests, with a minimum dist between nests enforced
n.nests <- 500#1162
# nest.probs <- c(0,.02,.03,.09,.14,.2,.52) #probabilities based on the number of nests in each RSF bin
nest.probs <- rep(.25,3)
size.list <- nest.probs*n.nests

samp.comb <- list()
for (i in 1:3){
	temp.sub <- samp2[samp2@data$bin==i,]  #subset data by bins
	temp.samp <- temp.sub[sample(nrow(temp.sub),size=size.list[i]),]  #sample for specified size
	samp.comb[[i]] <- temp.samp #put results in lsit
}
new.samp <- do.call("rbind",samp.comb)  #combine list into single spatial data frame
new.samp@data$siteID <- seq(1,nrow(new.samp@data))

#verify nearest neighbor nest distances
near.neigh <- as.data.frame(nndist(X=new.samp$x,Y=new.samp$y))
names(near.neigh) <- "dist"
near.neigh$recno <- 1:nrow(near.neigh)
summary(near.neigh)
plot(new.samp)

#check number of nests per township
nests.town <- over(new.samp,townships)
nests.per <- as.numeric(table(nests.town$joinfield))
summary(nests.per)
hist(nests.per)
table(nests.per)

#demographic data for simulated nests
sim.nests <- new.samp@data
#sim.nests$site <- seq(1:nrow(sim.nests))
sim.nests$X <- new.samp$x
sim.nests$Y <- new.samp$y

#set simulation parameters
occupancy <- 0.6  ###ZW: We should define these parameters as distributions, with hyperparameters based on the SE of estimates
# e.g.
occupancy.se <- 0.1
success <- 0.6
young <- 0.7
occ.det <- 0.5
state.det1 <- 0.6 #number of states to be adjusted
state.det2 <- 0.7
state.det3 <- 0.8
n.years <- 10

#results <- matrix(nrow=nrow(sim.nests),ncol=8)
results <- matrix(nrow=n.years,ncol=8)
results2 <- NULL
for (j in 1:nrow(sim.nests)){
	site <- j
	for (i in 1:n.years){
		year <- i
		occ <- rbinom(1,1,rnorm(1,occupancy,occupancy.se))
		if (occ==0){
			succ=0;yng=0
		} else {
		succ <- rbinom(1,1,success)
		} 
		if (occ==0 & succ==0) {
			yng=0
		} else if (occ==1 & succ==1) {	
			yng <- rpois(1,young)
		}
		if (occ==0){
			vis1=0;vis2=0;vis3=0;vis4=0
		} else if (occ==1 & succ==0) {
			vis1=rbinom(1,1,occ.det)
			vis2=rbinom(1,1,occ.det)
			vis3=rbinom(1,1,occ.det)
			vis4=rbinom(1,1,occ.det)
		} else if (occ==1 & yng>=1) {
			vis1=rbinom(1,1,occ.det)
			vis2=rbinom(1,1,state.det1)
			vis3=rbinom(1,1,state.det2)
			vis4=rbinom(1,1,state.det3)
		}
	results[i,c(1:8)] <-cbind(year,occ,succ,yng,vis1,vis2,vis3,vis4)
	}
	results2 <- rbind(results2,cbind(results,j))
}

colnames(results2) <- c("year","occ","succ","yng","vis1","vis2","vis3","vis4","siteID")
head(results2)

data.sim <- merge(sim.nests,results2)

summary(data.sim$occ)
sd(data.sim$occ) #0.49
sd(data.sim$occ)/sqrt(length(data.sim$occ)) #0.008

#make hexagon grid MUST RUN FUNCTION "make_grid" AT END OF CODE FIRST
hex.grid2 <- make_grid(state,cell_diameter=25000)
hex.grid2 <- SpatialPolygonsDataFrame(hex.grid2,data.frame(ID=1:length(hex.grid2)))
#get average (mode) RSF value in each hex
###this code takes too long- have to find an alternative
hex.modes <- zonal.stats(x=hex.grid2,y=hawk.rsf,stat=getmode,plot=F) #function for getmode is at end
hex.grid2@data$rsfmode <- hex.modes

spplot(hex.grid2,zcol="rsfmode")


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












