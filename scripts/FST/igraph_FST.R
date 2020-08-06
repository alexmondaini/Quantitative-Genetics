df <- read.csv(file = 'files/csv/pairwise_nei_fst_data_frame2.csv',stringsAsFactors = F)
library(tidyr)
library(gtools)

#tidyr
df1 <- gather(df,'POP2','FST',2:23,na.rm = T)
rownames(df1) <- 1:nrow(df1)
str(df1)
# gtools
df1<-df1[mixedorder(df1$POP1),]
df1$POP2 <- gsub('\\.', '-', df1$POP2)
rownames(df1) <- 1:nrow(df1)
# preparing the new vector for network graph
pop1_unique <- unique(df1$POP1)
pop1_unique <- rep(pop1_unique,each=2)
pop1_unique <- pop1_unique[-1]
pop1_unique <- pop1_unique[-length(pop1_unique)]
# igraph
require(igraph)
pop1_unique2 <- c("ES1-3", "ES4-6",
                  "ES4-6", "ES6-8",
                  "ES6-8", "ES7-9",
                  "ES7-9", "ES10-12",
                  "ES10-12", "ES13-15",
                  "ES13-15", "ES16-18",
                  "ES16-18", "ES19-21",
                  "ES19-21", "ES22-24",
                  "ES22-24", "ES25-27",
                  "ES25-27", "ES28-30",
                  "ES28-30", "ES31-33",
                  "ES31-33", "ES34-36",
                  "ES34-36", "ES37-38",
                  "ES13-15", "SA1-3",
                  "SA1-3", "SA4-6", 
                  "SA4-6", "SA7-9",
                  "SA7-9", "SA10-12",
                  "SA10-12", "SA13-15", 
                  "SA13-15","SA16-18",
                  "SA16-18", "SA19-21",
                  "SA19-21","SA22-25")

index_of_Fst <- match(unique(df1$POP1),df1$POP1)
which(index_of_Fst=='196')
index_of_Fst  <- index_of_Fst[-14]

# insert function to add an index that I want 

insert.at <- function(a, pos, ...){
  dots <- list(...)
  stopifnot(length(dots)==length(pos))
  result <- vector("list",2*length(pos)+1)
  result[c(TRUE,FALSE)] <- split(a, cumsum(seq_along(a) %in% (pos+1)))
  result[c(FALSE,TRUE)] <- dots
  unlist(result)
}

index_of_Fst <- insert.at(index_of_Fst,13,104)

Fst <- df1[index_of_Fst,3]
Fst

# plotting igraph
g1 <- make_graph(edges=pop1_unique2)

# how many edges the graph has ?

E(g1)$weight = Fst
E(g1)$weight

# plotting 
ColMap = rep('green',vcount(g1))
ColMap[grep('^SA',V(g1)$name)] = 'red'

## Create an explicit layout of the vertices 
## that will separate the labels
LO = matrix(0, nrow=22, ncol=2)
LO[1:5,1]   = -5:-1
LO[7:14,1]  = (1:8)/1.7
LO[15:22,1] = (1:8)/1.7
LO[7:14,2]  = -(1:8)
LO[15:22,2] = 2:9


plot(g1, layout=LO,
     edge.label = E(g1)$weight,
     edge.label.cex = 0.5,
     edge.arrow.size=.4,
     vertex.color=ColMap,
     vertex.size=7, 
     vertex.frame.color="black",
     vertex.label.color="black", 
     vertex.label.cex=0.7,
     vertex.label.dist=c(rep(1.2,5), rep(2.2,17)),
     vertex.label.degree = c(rep(-pi/2, 5), 0, rep(-0.1,8), rep(0.1,8)),
     margin=-.2,
     vertex.shape='circle')

#### another plot with base population

df2 <- df1[df1$POP1=='ES1-3',]
df2

g2 <- graph_from_data_frame(d = df2[,c(1,2)],directed = TRUE)

ColMap = rep('seagreen3',vcount(g2))
ColMap[grep('^SA',V(g2)$name)] = 'plum3'

E(g2)$weight = df2$FST


LO = matrix(0, nrow=22, ncol=2)
LO <- layout.circle(g2)
LO[1,1]  <- -.5
LO[1,2]  <- 0
LO[2,1]  <- 1.09
LO[22,1] <- 1.09
LO[21,1] <- .94
LO[21,2] <- -0.67
LO[20,1] <- .65
LO[20,2] <- -.85

plot(g2,
     edge.arrow.size=0.5,
     vertex.shape="circle",
     layout = LO,
     vertex.color=ColMap,
     vertex.frame.color="white",
     vertex.label.color="black",
     vertex.label.cex = .4,
     vertex.label.dist= c(0,rep(2.2,7),0,0,0,0,0,0,rep(0,8)),
     edge.label = E(g2)$weight,
     edge.label.cex = 0.6,
     edge.curved=T)

layout.circle(g2)


##### get one more plot of the parallel populations

### g3 

df3 <- df1[c(119,134,148,154,161,173,184,194,203,104,109,index_of_Fst[c(-1:-5,-14)]),]
df3
df3<-df3[mixedorder(df3$POP1),]
rownames(df3) <- 1:nrow(df3)
df3

df3 <- df3[c(-2,-12),]
rownames(df3) <- 1:nrow(df3)
df3

require(igraph)
g3 <- graph_from_data_frame(d = df3[,c(1,2)],directed = F)

ColMap = rep('seagreen3',vcount(g3))
ColMap[grep('^SA',V(g3)$name)] = 'plum3'

E(g3)$weight = df3$FST
E(g3)$weight[c(4,6,8,9,12,14,16,18,19,20,21,22,23,24)] <- NA
E(g3)$weight
df3

set.seed(345)
plot(g3,
     vertex.color=ColMap,
     vertex.frame.color="white",
     vertex.label.color="black",
     vertex.label.cex = .6,
     edge.label.cex = 0.6,
     edge.label = E(g3)$weight,
     layout = L1,
     edge.curved=F,
     main='FST parallel')


L1 <- layout.auto(g3)
L1[1,1] <- 17

df3 <- structure(list(POP1 = c("ES13-15", "ES13-15", "ES16-18", "ES16-18", 
                              "ES19-21", "ES19-21", "ES22-24", "ES22-24", "ES25-27", "ES25-27", 
                              "ES28-30", "ES28-30", "ES31-33", "ES31-33", "ES34-36", "ES34-36", 
                              "ES37-38", "SA1-3", "SA4-6", "SA7-9", "SA10-12", "SA13-15", "SA16-18", 
                              "SA19-21"), POP2 = c("SA1-3", "ES16-18", "SA1-3", "ES19-21", 
                              "SA4-6", "ES22-24", "SA7-9", "ES25-27", "ES28-30", "SA10-12", 
                              "SA13-15", "ES31-33", "SA16-18", "ES34-36", "SA19-21", "ES37-38", 
                              "SA22-25", "SA4-6", "SA7-9", "SA10-12", "SA13-15", "SA16-18", 
                              "SA19-21", "SA22-25"), FST = c(0.0129, 0.0109, 0.0176, 0.012, 
                              0.0564, 0.0275, 0.0193, 0.0108, 0.029, 0.0416, 0.0398, 0.0061, 
                              0.0541, 0.0134, 0.0402, 0.0205, 0.024, 0.0081, 0.0035, 0.02, 
                              0.0094, 0.0095, 0.0161, 0.0199)), row.names = c(NA, 24L), class = "data.frame")
df3
