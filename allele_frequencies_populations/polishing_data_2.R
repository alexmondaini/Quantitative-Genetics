library(tidyverse)

# load the data
load(here('allele_frequencies_populations','al_freq.rda'))

# transpose the data
al_freq[1:6,1:6]
al_freq <- t(al_freq)
al_freq <- as.data.frame(al_freq,row.names = row.names(al_freq))

# rownames to column
al_freq <- rownames_to_column(al_freq)
names(al_freq)[1] <- 'locus'
head(al_freq)

# removing UN (unidentified columns)
UN_out <- grep('UN-',al_freq$locus)
al_freq <- al_freq[-UN_out,]
head(al_freq)
al_freq$locus <- gsub('-','_',al_freq$locus)
al_freq$locus <- paste0('S',al_freq$locus)

# view of data and its structure
head(al_freq)
str(al_freq)

## grabbing chromosome
chr <- substr(al_freq$locus,start = 2,stop = 3)
al_freq$locus <- gsub('S[1-7][A-D]_','',al_freq$locus)
al_freq[1:6,1:6]
al_freq <- add_column(al_freq,chr,.before = 'locus')
head(al_freq)


# we need the deviation from base population ES1-7 for time selection
time_selection <- al_freq[,3]-al_freq[,4:ncol(al_freq)]
head(time_selection)
time_selection <- data.frame(al_freq[,1:2,drop=F],time_selection,stringsAsFactors = F,check.names = F)
time_selection[1:6,1:6]
str(time_selection)

## arrange the data back long format
time_selection <-  time_selection %>%
                gather(Trial,MAF,c(-locus,-chr)) %>%
                mutate(Trial = factor(Trial))
head(time_selection)

## factor trial by by its levels in temporal fashion
time_selection$Trial = factor(time_selection$Trial,levels(time_selection$Trial)[c(3,1,2,4,5)])
levels(time_selection$Trial)


# create an environment selection data frame

env_selection <-   data.frame(al_freq[,1:2,drop=F],al_freq[,6,drop=F]-al_freq[,8,drop=F]
                              ,stringsAsFactors = F,check.names = F )       

# create a trial factor column with one level
Trial <- factor(rep("ES26_38-SA13_25",nrow(env_selection)))
env_selection <- add_column(env_selection,Trial,.after = 'locus')
names(env_selection)[4] <- 'MAF'
head(env_selection)
str(env_selection)

