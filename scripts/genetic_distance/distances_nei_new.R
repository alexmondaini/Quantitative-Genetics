require(poppr)

nei_within <- nei.dist(obj,warning = T)
load(file = "Rda/nei_within.rda")
class(nei_within)

#-------------------------

nei_within_mat <- data.matrix(nei_within)
nei_within_mat[1:10,1:10]


load(file = 'Rda/Fst-Nei/fst_hierfstat.rda')

dt <- fst_hierfstat[,1:3]
dt[1:10,1:3]
row.names(nei_within_mat) <- dt[,2]

indices <- which(lower.tri(nei_within_mat,diag = F),arr.ind = T)
indices

nei_dist_df <- data.frame(GID= dimnames(nei_within_mat)[[2]][indices[,2]],
                          Trial = dimnames(nei_within_mat)[[1]][indices[,1]],
                          Nei_distances = nei_within_mat[indices]
                          )


load(file = "Rda/nei_dist_df.rda")

nei_dist_df[1:10,1:3]


nei_dist_df <- nei_dist_df[grep('SA',nei_dist_df$Trial),]
nei_dist_df <- nei_dist_df[grep('ES',nei_dist_df$Trial),]
nei_dist_df <- nei_dist_df[grep('ES',nei_dist_df$Trial),]
nei_dist_df <- nei_dist_df[grep('ES[1-9]$|ES1[1-3]$',nei_dist_df$Trial),]


#---------------- grouping together the trials

nei_dist_df$Trial <- gsub("ES[1-4]$", "ES1-4", nei_dist_df$Trial)
nei_dist_df$Trial <- gsub("ES[5-7]$", "ES5-7", nei_dist_df$Trial)
nei_dist_df$Trial <- gsub("ES[8-9]$|ES10", "ES8-10", nei_dist_df$Trial)
nei_dist_df$Trial <- gsub("ES1[1-3]$", "ES11-13", nei_dist_df$Trial)

nei_dist_df$Trial <- gsub("SA[1-3]$|ES1[4-6]$", "ESSA1-3", nei_dist_df$Trial)
nei_dist_df$Trial <- gsub("SA[4-6]$|ES1[7-9]$", "ESSA4-6", nei_dist_df$Trial)
nei_dist_df$Trial <- gsub("SA[7-9]$|ES2[0-2]$", "ESSA7-9", nei_dist_df$Trial)
nei_dist_df$Trial <- gsub("SA1[0-2]$|ES2[3-5]$", "ESSA10-12", nei_dist_df$Trial)
nei_dist_df$Trial <- gsub("SA1[3-5]$|ES2[6-8]$", "ESSA13-15", nei_dist_df$Trial)
nei_dist_df$Trial <- gsub("SA1[6-8]$|ES29|ES3[0-1]$", "ESSA16-18", nei_dist_df$Trial)
nei_dist_df$Trial <- gsub("SA19|SA2[0-1]$|ES3[2-4]$", "ESSA19-21", nei_dist_df$Trial)
nei_dist_df$Trial <- gsub("SA2[2-5]$|ES3[5-8]$", "ESSA22-25", nei_dist_df$Trial)

#-------------------------- grouping separetely the trials

nei_dist_df$Trial <- gsub("ES[1-4]$", "ES1_4", nei_dist_df$Trial)
nei_dist_df$Trial <- gsub("ES[5-7]$", "ES5_7", nei_dist_df$Trial)
nei_dist_df$Trial <- gsub("ES[8-9]$|ES10", "ES8_10", nei_dist_df$Trial)
nei_dist_df$Trial <- gsub("ES1[1-3]$", "ES11_13", nei_dist_df$Trial)

nei_dist_df$Trial <- gsub("SA[1-3]$", "SA1_3", nei_dist_df$Trial)
nei_dist_df$Trial <- gsub("ES1[4-6]$", "ES14_16", nei_dist_df$Trial)

nei_dist_df$Trial <- gsub("SA[4-6]$", "SA4_6", nei_dist_df$Trial)
nei_dist_df$Trial <- gsub("ES1[7-9]$", "ES17_19", nei_dist_df$Trial)

nei_dist_df$Trial <- gsub("SA[7-9]$", "SA7_9", nei_dist_df$Trial)
nei_dist_df$Trial <- gsub("ES2[0-2]$", "ES20_22", nei_dist_df$Trial)

nei_dist_df$Trial <- gsub("SA1[0-2]$", "SA10_12", nei_dist_df$Trial)
nei_dist_df$Trial <- gsub("ES2[3-5]$", "ES23_25", nei_dist_df$Trial)

nei_dist_df$Trial <- gsub("SA1[3-5]$", "SA13_15", nei_dist_df$Trial)
nei_dist_df$Trial <- gsub("ES2[6-8]$", "ES26_28", nei_dist_df$Trial)

nei_dist_df$Trial <- gsub("SA1[6-8]$", "SA16_18", nei_dist_df$Trial)
nei_dist_df$Trial <- gsub("ES29|ES3[0-1]$", "ES29_31", nei_dist_df$Trial)

nei_dist_df$Trial <- gsub("SA19|SA2[0-1]$", "SA19_21", nei_dist_df$Trial)
nei_dist_df$Trial <- gsub("ES3[2-4]$", "ES32_34", nei_dist_df$Trial)

nei_dist_df$Trial <- gsub("SA2[2-5]$", "SA22_25", nei_dist_df$Trial)
nei_dist_df$Trial <- gsub("ES3[5-8]$", "ES35_38", nei_dist_df$Trial)

#--------------------- tranbsforming the string into a factor and ordering the levels 

nei_dist_df$Trial <- factor(nei_dist_df$Trial)
levels(nei_dist_df$Trial)

nei_dist_df$Trial <- factor(nei_dist_df$Trial, levels(nei_dist_df$Trial)[c(1,3,4,2,5,11,12,6,7,8,9,10)]) # ordering for grouping together

nei_dist_df$Trial <- factor(nei_dist_df$Trial, 
                            rev(levels(nei_dist_df$Trial))[c(1,11,12,2,13,3,19,4,20,5,14,6,15,7,16,8,17,9,18,10)]) # ordering for grouping separetely

nei_dist_df$Trial <- factor(nei_dist_df$Trial, levels(nei_dist_df$Trial)[c(1,7,8,2,3,4,5,6)]) # ordering for SA only

nei_dist_df$Trial <- factor(nei_dist_df$Trial, levels(nei_dist_df$Trial)[c(1,11,12,2:10)]) # ordering for ES only

nei_dist_df$Trial <- factor(nei_dist_df$Trial, levels(nei_dist_df$Trial) [c(1,3,4,2)])


#-------------------- ploting boxplot

bp <- with(nei_dist_df, boxplot(Nei_distances~Trial))


bp <- with(nei_dist_df, boxplot(Nei_distances~Trial, width=proportion, col=c("orange","seagreen")) )

new_order <- with(nei_dist_df,reorder(Trial, Nei_distances, median, na.rm=T)) # order by increasing median
bp <- with(nei_dist_df, boxplot(Nei_distances~new_order, col=c("orange","seagreen")) )

require(ggplot2)

ggplot(nei_dist_df,aes(x=Trial,y=Nei_distances)) +  
  geom_boxplot()



#----------------------------- boxplot plus Tukey


model = lm(nei_dist_df$Nei_distances~nei_dist_df$Trial)
ANOVA = aov(model)

TUKEY <- TukeyHSD(x=ANOVA,'nei_dist_df$Trial',conf.level = 0.95)
TUKEY

par(mar=c(5.1,10,2,1))
plot(TUKEY,las=1,col="brown")

generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}

# Apply the function on my dataset
LABELS <- generate_label_df(TUKEY , "nei_dist_df$Trial")
LABELS

# A panel of colors to draw each group with the same color :
my_colors <- c( 
  rgb(143,199,74,maxColorValue = 255),
  rgb(139,69,19,maxColorValue = 255), 
  rgb(128,128,128,maxColorValue = 255),
  rgb(128,0,128,maxColorValue = 255),
  rgb(128,128,128,maxColorValue = 255),
  rgb(128,128,128,maxColorValue = 255),
  rgb(127,255,212,maxColorValue = 255),
  rgb(128,128,128,maxColorValue = 255),
  rgb(128,128,128,maxColorValue = 255),
  rgb(242,104,34,maxColorValue = 255),
  rgb(111,145,202,maxColorValue = 255),
  rgb(128,128,128,maxColorValue = 255),
  rgb(128,128,128,maxColorValue = 255),
  rgb(128,128,128,maxColorValue = 255)
)

par(cex.lab=1.5) # is for y-axis

par(cex.axis=1.1) # is for x-axis


a <- boxplot(nei_dist_df$Nei_distances ~ nei_dist_df$Trial , ylim=c(min(nei_dist_df$Nei_distances) , 1.1*max(nei_dist_df$Nei_distances)) 
             , col=my_colors[as.numeric(LABELS[,1])] , ylab="Nei_distances" , main="")
a

# I want to write the letter over each box. Over is how high I want to write it.
over <- 0.18*max( a$stats[nrow(a$stats),] )


#Add the labels
text( c(1:nlevels(nei_dist_df$Trial)) , a$stats[nrow(a$stats),]+over , LABELS[,1]  , col=my_colors[as.numeric(LABELS[,1])],cex=2 )

#---------------------- ridgline plot

nei_dist_df[1:10,1:3]

require(ggplot2)
require(ggridges)

ggplot(nei_dist_df,aes( x=Nei_distances,  y= factor(Trial,rev(levels(Trial))) , fill= ..x..) ) +
       geom_density_ridges_gradient(rel_min_height = 0.001)+
       scale_fill_viridis_c(name='Nei distances',option = 'C') + ylab('Trial')


ggplot(nei_dist_df,  aes (x=Nei_distances, y=factor(Trial,rev(levels(Trial))),fill=factor(..quantile..))) +
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,
    quantiles = 4, quantile_lines = TRUE
  )+
  scale_fill_viridis_d(name='Quartiles',option = 'C') + ylab('Trial')


ggplot(nei_dist_df, aes(x=Nei_distances, y=factor(Trial,rev(levels(Trial)))) ) +
  stat_density_ridges(quantile_lines=T ,quantiles = 4) + ylab('Trial')


#--------------------- facet_wrap() multiple plots laid one next to the other

require(ggplot2)
require(hrbrthemes)


ggplot(data = nei_dist_df, aes(x=Nei_distances,fill=Trial)) +
  geom_density()+
  theme_ipsum_rc() +
  facet_wrap(~Trial)+
  theme(
    legend.position = 'none',
    panel.spacing = unit(1,'lines'),
    axis.ticks.x = element_blank()
  )




