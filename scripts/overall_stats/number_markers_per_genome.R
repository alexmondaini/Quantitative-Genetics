load(file = 'Rda/Fst-Nei/fst_hierfstat.rda')

a <- colnames(fst_hierfstat)[-1:-3]
a <- substr(a,1,2)
table(a)


A<-sum(table(a)[c(1,4,7,10,13,16,19)])
B<-sum(table(a)[c(2,5,8,11,14,17,20)])
D<-sum(table(a)[c(3,6,9,12,15,18,21)])
UN <- sum(table(a)[22])


require(dplyr)

fst_hierfstat %>% group_by(pop) %>% count()

table(fst_hierfstat$pop)

cfo <- read.csv(file = 'excel_text_files/cfo.csv',stringsAsFactors = F)
cfo

colnames(cfo)[1] <- 'chr'
colnames(cfo)



