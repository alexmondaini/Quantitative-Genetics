require(dplyr)
require(ggplot2)

load(file = 'Rda/allele_frequencies/al_freq_dev.rda')
head(al_freq_dev)
levels(al_freq_dev$Trial)

load(file='Rda/allele_frequencies/env_selection.rda')
load(file='Rda/allele_frequencies/selection_time.rda')

levels(al_freq_dev$Trial)
levels(env_selection$Trial)
levels(selection_time$Trial)

################################################

q95_s <- quantile(al_freq_dev$MAF,.99)
q05_s <- quantile(al_freq_dev$MAF,.01)

al_freq_dev_1 <- al_freq_dev[al_freq_dev$MAF>q95_s|al_freq_dev$MAF<q05_s,]

# the deviation from the first to the last

q95_l <- quantile(selection_time$MAF,.99)
q05_l <- quantile(selection_time$MAF,.01)

selection_time_1 <- selection_time[selection_time$MAF>q95_l|selection_time$MAF<q05_l,]

# the deviation from the last two

q95_v <- quantile(env_selection$MAF,.99)
q05_v <- quantile(env_selection$MAF,.01)

env_selection_1 <- env_selection[env_selection$MAF>q95_v|env_selection$MAF<q05_v,]

which(env_selection$chr=='1D')
i <- env_selection[4285,]
which(env_selection$chr=='3D')
j <- env_selection[17027,]

env_selection_1 <- rbind(env_selection_1,i,j)
tes

which(tes$chr=='3D')
tes[793,]
####### plotting the combined dataframes MAF distribution

library(gridExtra)

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


p0 <- ggplot(al_freq_dev,aes(x=MAF)) +
      geom_density(lwd=1)+
      labs(title='a)   Temporal MAF deviation (ES1-7-all)',
      x='',
      y='')+
      xlim(-0.4,0.4)

archt0 <- ggplot_build(p0)$data[[1]]

p0 <- p0 + geom_area(data = subset(archt0,  x<q05_s ), aes(x=x, y=y,fill="lower 1%"))+
           geom_area(data = subset(archt0,  x>q95_s ), aes(x=x, y=y,fill="upper 1%"))+
           scale_fill_manual(values = c('#E54B2A','#2B41E7'),name='Percentiles')

legend_r <- get_legend(p0)
legend_r

p0 <- p0 + theme(legend.position = "none")
p0



p1 <- ggplot(env_selection,aes(x=MAF)) +
      theme(legend.position = "none") +
      geom_density(lwd=1)+
      labs(title="b)   Environmental MAF deviation (ES26-38-SA13-25)",
      x='MAF deviation',
      y = '') + 
      xlim(-0.4,0.4)


p1

archt1 <- ggplot_build(p1)$data[[1]]

p1 <- p1 + geom_area(data = subset(archt1,  x<q05_v ), aes(x=x, y=y,fill="lower 1%"))+
           geom_area(data = subset(archt1,  x>q95_v ), aes(x=x, y=y,fill="upper 1%"))+
           scale_fill_manual(values = c('#E54B2A','#2B41E7'),name='Percentiles')

p1

require(gridExtra)

grid.arrange(arrangeGrob(p0,p1),legend_r,ncol=2,widths=c(3.5,0.8) )
  
  
  
  
 # p0,p1,legend_r,nrow=2,ncol=3,widths = c(2.3,2.3,0.8)


##########################

save(env_selection_1,file='Rda/allele_frequencies/env_selection_1.rda')
save(al_freq_dev_1,file='Rda/allele_frequencies/al_freq_dev_1.rda')
save(selection_time_1,file = 'Rda/allele_frequencies/selection_time_1.rda')
