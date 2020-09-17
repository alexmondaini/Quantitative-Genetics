# load libraries
require(dplyr)
require(ggplot2)
require(gtools)
require(gridExtra)

# load data

# Defining the threshold quantiles of 1% from both sides of the distribution
# of MAF deviations

#over time
q95_time <- quantile(time_selection$MAF,.99)
q05_time <- quantile(time_selection$MAF,.01)

time_selection <- time_selection[time_selection$MAF>q95_time|time_selection$MAF<q05_time,]

#over environment

q95_env <- quantile(env_selection$MAF,.99)
q05_env <- quantile(env_selection$MAF,.01)

env_selection_1 <- env_selection[env_selection$MAF>q95_env|env_selection$MAF<q05_env,]

# note that chr 1D and 3D are missing let's just add one data point in each
# for plotting purposes
mixedsort(unique(env_selection_1$chr))
which(env_selection_1$chr=='1D')
which(env_selection_1$chr=='3D')

i <- env_selection[4285,]
j <- env_selection[17027,]

env_selection_1 <- rbind(env_selection_1,i,j)
env_selection <- env_selection_1
rm(env_selection_1)

# confirm the chr are back
mixedsort(unique(env_selection$chr))


# plotting the quantiles

# get legend function
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#over time
p0 <- ggplot(time_selection,aes(x=MAF)) +
      geom_density(lwd=1)+
      labs(title='a)   Temporal MAF deviation (ES1-7-all)',
      x='',
      y='')+
      scale_x_continuous(breaks = as.numeric(sprintf('%.1f',seq(-0.5,0.5,0.1))))

archt0 <- ggplot_build(p0)$data[[1]]

p0 <- p0 + geom_area(data = subset(archt0,  x<q05_time ), aes(x=x, y=y,fill="lower 1%"))+
           geom_area(data = subset(archt0,  x>q95_time ), aes(x=x, y=y,fill="upper 1%"))+
           scale_fill_manual(values = c('#E54B2A','#2B41E7'),name='Percentiles')

legend_r <- get_legend(p0)
p0 <- p0 + theme(legend.position = "none")
p0


#over environment
p1 <- ggplot(env_selection,aes(x=MAF)) +
      theme(legend.position = "none") +
      geom_density(lwd=1)+
      labs(title="b)   Environmental MAF deviation (ES26-38-SA13-25)",
      x='MAF deviation',
      y = '') + 
      scale_x_continuous(breaks = as.numeric(sprintf('%.1f',seq(-0.3,0.3,0.1))))


p1

archt1 <- ggplot_build(p1)$data[[1]]

p1 <- p1 + geom_area(data = subset(archt1,  x<q05_env ), aes(x=x, y=y,fill="lower 1%"))+
           geom_area(data = subset(archt1,  x>q95_env ), aes(x=x, y=y,fill="upper 1%"))+
           scale_fill_manual(values = c('#E54B2A','#2B41E7'),name='Percentiles')
p1



# combine plot time and environment
grid.arrange(arrangeGrob(p0,p1),legend_r,ncol=2,widths=c(3.5,0.8) )
  
