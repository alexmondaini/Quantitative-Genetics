require(gtools)
require(hierfstat)

load('Rda/Fst-Nei/fst_hierfstat.rda')
fst_hierfstat[1:6,1:9]

##### basic stats form hierfstat

basic_stats <- basic.stats(fst_hierfstat[,-1:-2],diploid = T,digits = 3)

class(basic_stats$n.ind.samp)
class(basic_stats$perloc)
class(basic_stats$pop.freq)
class(basic_stats$pop.freq$X1A.1145442)

save(basic_stats,file = 'Rda/Fst-Nei/basic_stats.rda')

# the allele frequencies are stored in a list of table objects 
load('Rda/Fst-Nei/basic_stats.rda')

al_freq <- basic_stats$pop.freq
al_freq[1:3]

al_freq <- sapply(al_freq, function(x) apply(x,2,min))
al_freq[,1:3]

al_freq <- al_freq[mixedorder(rownames(al_freq)),]
al_freq[,1:3]

colnames(al_freq) <- sub('X','',colnames(al_freq))
colnames(al_freq) <- sub('.','-',colnames(al_freq),fixed = T)

al_freq[,1:3]

save(al_freq,file = 'Rda/allele_frequencies/al_freq.rda')
