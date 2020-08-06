library(tidyverse)
load(file = 'Rda/allele_frequencies/al_freq.rda')

al_freq[1:6,1:7]
al_freq <- t(al_freq)
al_freq <- as.data.frame(al_freq,row.names = row.names(al_freq))

al_freq <- rownames_to_column(al_freq)
names(al_freq)[1] <- 'locus'
head(al_freq)

UN_out <- grep('UN-',al_freq$locus)
al_freq <- al_freq[-UN_out,]
head(al_freq)
al_freq$locus <- gsub('-','_',al_freq$locus)
al_freq$locus <- paste0('S',al_freq$locus)

head(al_freq)
str(al_freq)

save(al_freq,file = 'Rda/allele_frequencies/al_freq.rda')

##############################

al_freq_dev <- al_freq[,2]-al_freq[,3:ncol(al_freq)]
head(al_freq_dev)
al_freq_dev <- data.frame(al_freq[,1,drop=F],al_freq_dev,stringsAsFactors = F,check.names = F)
al_freq_dev[1:6,1:6]
str(al_freq_dev)
## grabbing chromosome
chr <- substr(al_freq_dev$locus,start = 2,stop = 3)
al_freq_dev$locus <- gsub('S[1-7][A-D]_','',al_freq_dev$locus)
al_freq_dev[1:6,1:6]
al_freq_dev <- add_column(al_freq_dev,chr,.before = 'locus')
head(al_freq_dev)

al_freq_dev <-  al_freq_dev %>%
                gather(Trial,MAF,c(-locus,-chr)) %>%
                mutate(Trial = factor(Trial))
al_freq_dev
al_freq_dev$Trial = factor(al_freq_dev$Trial,levels(al_freq_dev$Trial)[c(3,1,2,4,5)])
levels(al_freq_dev$Trial)
save(al_freq_dev,file = 'Rda/allele_frequencies/al_freq_dev.rda')

################################

selection_time  <-  al_freq_dev %>% filter(Trial == 'SA13-25'|Trial=='ES26-38')
selection_time$Trial <- factor(selection_time$Trial)
levels(selection_time$Trial)
save(selection_time,file = 'Rda/allele_frequencies/selection_time.rda')
head(selection_time)

################################
library(tibble)


env_selection <-   data.frame(al_freq[,1,drop=F],al_freq[,5,drop=F]-al_freq[,7,drop=F]
                              ,stringsAsFactors = F,check.names = F )         

colnames(env_selection)[2] <- "MAF"
head(env_selection)

chr <- substr(env_selection$locus,start = 2,stop = 3)
env_selection$locus <- gsub('S[1-7][A-D]_','',env_selection$locus)
head(env_selection)
env_selection <- add_column(env_selection,chr,.before = 'locus')
head(env_selection)

Trial <- factor(rep("ES26_38-SA13_25",nrow(env_selection)))

env_selection <- add_column(env_selection,Trial,.after = 'locus')
head(env_selection)
str(env_selection)

save(env_selection,file = 'Rda/allele_frequencies/env_selection.rda')
#######################
