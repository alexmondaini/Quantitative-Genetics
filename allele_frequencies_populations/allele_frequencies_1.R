require(here)
require(gtools)
require(hierfstat)

load(here('FST','fst_hierfstat.rda'))
fst_hierfstat[1:6,1:9]

##### basic stats from hierfstat to calculate several statistics

basic_stats <- basic.stats(fst_hierfstat[,c(-1,-2,-4,-5)],diploid = T,digits = 3)

# check the class of different objects produced by the function

class(basic_stats$n.ind.samp)
class(basic_stats$perloc)
class(basic_stats$pop.freq)
class(basic_stats$pop.freq$X1A.1145442)

# we want to extract population allelic frequencies
al_freq <- basic_stats$pop.freq
al_freq[1]

al_freq <- sapply(al_freq, function(x) apply(x,2,min))
al_freq[,1:3]

al_freq <- al_freq[mixedorder(rownames(al_freq)),]
al_freq[,1:3]

colnames(al_freq) <- sub('X','',colnames(al_freq))
colnames(al_freq) <- sub('.','-',colnames(al_freq),fixed = T)

al_freq[,1:3]

save(al_freq,file=here('allele_frequencies_populations','al_freq.rda'))
