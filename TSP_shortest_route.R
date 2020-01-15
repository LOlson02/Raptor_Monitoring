library(TSP) #Create efficient routes using 'traveling salesman problem'

hex.samp.cents <- SpatialPointsDataFrame(gCentroid(hex.samp, byid=T), hex.samp@data) #centroids of hex sample
dMat <- spDists(hex.samp.cents) #distance matrix 
tsp <- TSP(dMat,labels=hex.samp.cents@data$hexID,method='euclidean') #create tsp object
solution <- solve_TSP(tsp) #solve using TSP methods
tour_length(solution) #tour length in m


# create path from tour
# not necessary if all we want to do is compare lengths of tours for different sampling methods
path <- cut_tour(solution, "1") #cut circular tour in arbitrary location
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


plot(l, lty=2, col="coral")
plot(tmp, pch=20, add=T)





# compare optimization methods
methods <- c("nearest_insertion","farthest_insertion","cheapest_insertion","arbitrary_insertion","nn","repetitive_nn","two_opt")
tours <- sapply(methods, FUN = function(m) solve_TSP(tsp, method = m), simplify = FALSE)
dotchart(sort(c(sapply(tours, tour_length), optimal = 6273700)), xlab = "tour length", xlim = c(5500000, 8000000))

# tsp.ham <- insert_dummy(tsp, label="cut")
# solution.ham <- solve_TSP(tsp.ham, method="farthest_insertion")
# solution.ham