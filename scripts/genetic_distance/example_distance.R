load('Rda/Fst-Nei/fst_data.rda')


myd <- fst_data

# calculating distance and storing them 
grps <- as.character(unique(myd$pop))
out <- list()
for ( i in 1:length(grps)){
  grp = grps[i]
  mydi <- myd[myd$pop==grp,-1:-3]
  distance <- dist(mydi, method = "euclidean")
  distanceV <- as.matrix(distance)[lower.tri(as.matrix(distance))]
  out[[i]] <- data.frame(group=grp, distanceV=distanceV)
}
# single data for all groups 
resultdistance <- do.call("rbind", out)

# scaling with maximum distance 
resultdistance$distanceV_SCALED <- resultdistance$distanceV/max(resultdistance$distanceV)


resultdistance$group <- factor(resultdistance$group, levels(resultdistance$group)[c(3:6,1,2)])
levels(resultdistance$group)

# boxplot with added datapoints 
with(resultdistance, boxplot(distanceV_SCALED~group))

# Add data points
mylevels <- levels(resultdistance$group)
levelProportions <- summary(resultdistance$group)/nrow(resultdistance)

for(i in 1:length(mylevels)){
  thislevel <- mylevels[i]
  thisvalues <- resultdistance[resultdistance$group==thislevel, "distanceV_SCALED"]
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
  points(myjitter, thisvalues, pch=20, col=rgb(0,0.1,0.5,.9), cex=0.5) 
} 

# plotting with beanplot 
library("beanplot")
with(resultdistance, beanplot(distanceV_SCALED~group))



