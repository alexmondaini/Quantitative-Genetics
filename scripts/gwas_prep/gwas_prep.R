data1 <- read.csv(file = 'excel_text_files/gwas/ES_GWAS.csv',stringsAsFactors = F)
data2 <- read.csv(file = 'excel_text_files/csv/excel_to_run.csv',stringsAsFactors = F)
str(data1)
str(data2)
head(data1)
head(data2)

data1 <- data1[!duplicated(data1[c(9,10)]),]
require(dplyr)
data3 <- inner_join(data1,data2)
data3 <- data3[,c(-18,-20)]

# calculating the BLUE of the BLUE's
require(lme4)
head(data3)
str(data3)
data3$GID <- as.factor(data3$GID)
data3$Occ <- as.factor(data3$Occ)
nlevels(data3$Occ)
nlevels(data3$GID)

which(is.na(data3),arr.ind = T)

blue_model <- lmer( Value ~ 0 + GID + (1|Occ) , data = data3)

blues <- fixef(blue_model)
GID <- names(blues)

df <- data.frame(GID,blues)
row.names(df) <- 1:nrow(df)
df
df$GID<-gsub("^GID","",df$GID)
df
sub("GID","",df$GID)
write.csv(data3,file = "blueee.csv")

require(dplyr)

data3 <- inner_join(data3,df)
head(data3)

data3 <- data3 %>% select(GID,blues)
write.csv(data3,file='excel_text_files/gwas/data3.csv',row.names = F)


####
data3 <- read.csv(file = 'excel_text_files/gwas/data3.csv',stringsAsFactors = F)
head(data3)
#### read membership prob file
membership_prob <- read.csv(file = 'excel_text_files/gwas/membership_prob.csv',stringsAsFactors = F)
head(membership_prob)
colnames(membership_prob)[1] <- 'GID'
data3 <- inner_join(data3,membership_prob)
head(data3)
colnames(data3)[2] <- 'BLUE'
head(data3)
#colnames(data3)[3:8] <- gsub('\\.','-',colnames(data3)[3:8])
#head(data3)
write.csv(data3,'excel_text_files/gwas/data3.csv',row.names = F)

########

load(file = 'Rda/numeric_hapmap/gbs_data2.rda')
info <- gbs_data2[-1:-4]
info <-  info[colnames(info)%in%data3$GID]
info <-  info[match(data3$GID,colnames(info))]
info[1:6,1:6]
gbs_data2 <- cbind(gbs_data2[1:4],info)
gbs_data2[1:6,1:6]
class(gbs_data2)

#######

write.csv(gbs_data2,'excel_text_files/gwas/genodata.csv',row.names = F)
