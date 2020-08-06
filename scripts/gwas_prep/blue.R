require(lme4)
head(data3)
str(data3)
data3$GID <- as.factor(data3$GID)
data3$Occ <- as.factor(data3$Occ)
nlevels(data3$Occ)
nlevels(data3$GID)

write.csv(data3,file = "blueee.csv")


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
?grep
?gsub