load(file = 'Rda/numeric_hapmap/gbs_data2.rda')
gbs_data2[1:6,1:6]

av_data <- gbs_data2[,3:4] 
av_data

require(dplyr)

av <- av_data %>% group_by(chrom) %>% mutate(dev= lead(pos,1,order_by = pos)- pos) %>% 
      mutate(dev = replace(dev,is.na(dev),0)) %>% mutate(mean_dist = round(mean(dev),0))
av

av <- av[!duplicated(av$mean_dist),]
av <- av[,c(1,4)]
av

write.table(av,file = 'mean_dist_pos.txt',sep = ',',quote = F,row.names = F)
