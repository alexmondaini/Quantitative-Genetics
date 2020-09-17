#load packages
require(adegenet)
require(here)
require(gtools)
require(ape)

#load data
load(here('FST','fst_hierfstat.rda'))
fst_hierfstat[1:7,1:10]

# transform group vector into a factor with appropriate levels
fst_hierfstat$group <- factor(fst_hierfstat$group)
fst_hierfstat$group <- factor(fst_hierfstat$group,levels = mixedsort(levels(fst_hierfstat$group)))
levels(fst_hierfstat$group)

# extract chromosome information
chromo <- substr(names(fst_hierfstat)[6:length(fst_hierfstat)],start = 1,stop = 2)
unique(chromo)
chromo <- factor(chromo)

# create a genlight object to store data in an efficient format usign adegenet
obj <-new('genlight',gen=snp[,4:ncol(snp)],ploidy=(as.integer(rep(2,nrow(snp)))),ind.names=snp$Genotype)
pop(obj) <- fst_hierfstat$group
chromosome(obj) <- chromo
save(obj,file = here('adegenet_dapc','obj.rda'))

# we can remove the dataframe now from memory to save space, everything is stored into obj.
rm(fst_hierfstat)

# pca with adegenet
pca <- glPca(obj)

#pca plot
myCol <- colorplot(pca$scores,pca$scores,transp = T,cex=1)
abline(h=0,v=0,col='grey')
title('PCA of ESWYT-SAWYT lines')

#find clusters for dapc
grp <- find.clusters(obj,max.n.clust=50,glPca=pca) # 2000 pcas and 42 clusters chosen

#dapc
dapc1 <- dapc(obj,grp$grp,glPca = pca)
dapc2 <- dapc(obj,pop(obj),glPca = pca)

#plot dapc
scatter(dapc2,scree.da = F,cell=0,cstar=0,solid = .4,clab=0,leg=T,posi.leg = 'bottomleft',
        col = c("green","purple","orange","red","darkblue","lightblue"))
