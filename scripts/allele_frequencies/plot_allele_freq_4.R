###########
centromere <- read.csv2(file = 'excel_text_files/chromosome_features/centromere.csv',
                        header = T,stringsAsFactors = F)
head(centromere)

load(file = 'Rda/allele_frequencies/al_freq_dev.rda')
load(file = 'Rda/allele_frequencies/al_freq_dev_1.rda')
load(file = 'Rda/allele_frequencies/env_selection_1.rda')
load(file = 'Rda/allele_frequencies/selection_time_1.rda')


al_freq_dev$locus   <- as.numeric(al_freq_dev$locus)
al_freq_dev_1$locus <- as.numeric(al_freq_dev_1$locus)
env_selection_1$locus <- as.numeric(env_selection_1$locus)
selection_time_1$locus <- as.numeric(selection_time_1$locus)

# get from kernel script the data to plot

require(lattice)
require(latticeExtra)
options(scipen=999)

Mytheme <- standard.theme()
Mytheme$background$col <- 'white'
Mytheme$panel.background$col <- 'white'
Mytheme$superpose.symbol$alpha <- 0.5
Mytheme$superpose.symbol$pch <- 19
Mytheme$superpose.symbol$col <- c('red','lightblue')
Mytheme$superpose.symbol$col <- c("purple","orange","red","darkblue","lightblue")
Mytheme$superpose.symbol$col <- '#399F42'
leg.txt<-c("ES1_7 - ES26_38","ES1_7 - SA13_25")
leg.txt2 <- c("ES26_38 - SA13_25")

A1_x <- rep(508725917,2)
A1_y <- c(0,0.1)
A1_text <- c('Glu-A1','*')

B1_x <- c(rep(555767104,2),rep(670233646,2))   
B1_y <- rep(c(0,0.1),2)
B1_text <- c('Glu-B1','*','Lr46','*')

D1_x <- rep(412162629,2)
D1_y <- c(0,0.1)
D1_text <- c('Glu-D1','*')

A2_x <- c(rep(3965176,2),rep(36934661,2))
A2_y <- c(-.07,-.2,0.15,0.2)
A2_text <- c('2NS','*','Ppd-A1','*')

B2_x <- rep(56238084,2)
B2_y <- c(0,0.1)
B2_text <- c('Ppd-B1','*')

D2_x <- rep(33955649,2)
D2_y <- c(0,0.1)
D2_text <- c('Ppd-D1','*')

B3_x <- rep(8848496,2)
B3_y <- c(0,0.1)
B3_text <- c('Sr2/Fhb1','*')

A4_x <- rep(688099006,2)
A4_y <- c(0,0.1)
A4_text <- c('Wx-B1','*')

B4_x <- rep(30861349,2)
B4_y <- c(0,0.1)
B4_text <- c('Rht1','*')

D4_x <- c(rep(18781069,2),rep(221906243,2))
D4_y <- rep(c(0,0.1),2)
D4_text <- c('Rht2','*','Lr67','*')

A5_x <- rep(588550218,2)
A5_y <- c(0,0.1)
A5_text <- c('Vrn-A1','*')

B5_x <- rep(573815063,2)
B5_y <- c(0,0.1)
B5_text <- c('Vrn-B1','*')

D5_x <- c(rep(3590617,2),rep(3609650,2),rep(467183397,2))
D5_y <- c(-.07,-.2,0.15,0.2,0,0.1)
D5_text <- c('Pinb-D1','*','Pina-D1','*','Vrn-D1','*')

B7_x <- rep(740039246,2)
B7_y <- c(0,0.1)
B7_text <- c('Lr68','*')

D7_x <- rep(47409577,2)
D7_y <- c(0,0.1)
D7_text <- c('Lr34','*')
 
pdf('images/maf10.pdf')
pl <-  xyplot(MAF ~ locus | factor(chr), data = env_selection_1,
       layout = c(3,7),
       as.table=T,
       grid=F,
       groups = Trial,
       type = c('p'),
       cex = 0.1,
       alpha = 0.3,
       xlab=list(label='position (Mb)',cex=1),
       ylab=list(label='minor allele frequency deviation',cex=1), 
       par.settings = Mytheme,
       scales = list(x=list(at=c(5000000,100000000,200000000,300000000,400000000,500000000,
                                 600000000,700000000,800000000),
                            labels=c('50','100','200','300','400','500','600','700','800'),
                            cex=0.5),
                     y=list(#at=c('-0.2','-0.1','0','0.1','0.2','0.2'),
                            #labels=c('-0.2','-0.1','0','0.1','0.2','0.2'),
                            cex=0.5)),
       auto.key = list(title='',space='top',columns=1,lines=F,points=T,cex=0.9,leg.txt2),
       panel = function(x,y,...){
       panel.xyplot(x,y,...)
       panel.text(A1_x[panel.number()==1],A1_y[panel.number()==1],
                  label=A1_text[panel.number()==1],cex=.3)       
       panel.text(B1_x[panel.number()==2],B1_y[panel.number()==2],
                  label=B1_text[panel.number()==2],cex=.3)
       panel.text(D1_x[panel.number()==3],D1_y[panel.number()==3],
                  label=D1_text[panel.number()==3],cex=.3) 
       panel.text(A2_x[panel.number()==4],A2_y[panel.number()==4],
                  label=A2_text[panel.number()==4],cex=.3)  
       panel.text(B2_x[panel.number()==5],B2_y[panel.number()==5],
                  label=B2_text[panel.number()==5],cex=.3)
       panel.text(D2_x[panel.number()==6],D2_y[panel.number()==6],
                  label=D2_text[panel.number()==6],cex=.3)
       panel.text(B3_x[panel.number()==8],B3_y[panel.number()==8],
                  label=B3_text[panel.number()==8],cex=.3)
       panel.text(A4_x[panel.number()==10],A4_y[panel.number()==10],
                  label=A4_text[panel.number()==10],cex=.3)
       panel.text(B4_x[panel.number()==11],B4_y[panel.number()==11],
                  label=B4_text[panel.number()==11],cex=.3)
       panel.text(D4_x[panel.number()==12],D4_y[panel.number()==12],
                  label=D4_text[panel.number()==12],cex=.3)
       panel.text(A5_x[panel.number()==13],A5_y[panel.number()==13],
                  label=A5_text[panel.number()==13],cex=.3)
       panel.text(B5_x[panel.number()==14],B5_y[panel.number()==14],
                  label=B5_text[panel.number()==14],cex=.3)
       panel.text(D5_x[panel.number()==15],D5_y[panel.number()==15],
                  label=D5_text[panel.number()==15],cex=.3)
       panel.text(B7_x[panel.number()==20],B7_y[panel.number()==20],
                  label=B7_text[panel.number()==20],cex=.3)
       panel.text(D7_x[panel.number()==21],D7_y[panel.number()==21],
                  label=D7_text[panel.number()==21],cex=.3)
       A1 <- 210200000:215800000
       B1 <- 237700000:243500000
       D1 <- 166200000:173800000
       A2 <- 326300000:327000000
       B2 <- 344400000:351300000
       D2 <- 264400000:272500000
       A3 <- 316900000:319900000
       B3 <- 345800000:347000000
       D3 <- 237100000:243200000
       A4 <- 264100000:267900000
       B4 <- 303900000:304400000
       D4 <- 182300000:188200000
       A5 <- 108900000:109100000
       B5 <- 198900000:202500000
       D5 <- 185600000:188700000
       A6 <- 283300000:288700000
       B6 <- 323000000:327500000
       D6 <- 211900000:217400000
       A7 <- 360200000:363800000
       B7 <- 288200000:288300000
       D7 <- 336300000:341700000
       panel.xblocks(A1[panel.number()==1],A1,col = 'lightblue',alpha=0.4)
       panel.xblocks(B1[panel.number()==2],B1,col = 'lightblue',alpha=0.4)
       panel.xblocks(D1[panel.number()==3],D1,col = 'lightblue',alpha=0.4)
       panel.xblocks(A2[panel.number()==4],A2,col = 'lightblue',alpha=0.4)
       panel.xblocks(B2[panel.number()==5],B2,col = 'lightblue',alpha=0.4)
       panel.xblocks(D2[panel.number()==6],D2,col = 'lightblue',alpha=0.4)
       panel.xblocks(A3[panel.number()==7],A3,col = 'lightblue',alpha=0.4)
       panel.xblocks(B3[panel.number()==8],B3,col = 'lightblue',alpha=0.4)
       panel.xblocks(D3[panel.number()==9],D3,col = 'lightblue',alpha=0.4)
       panel.xblocks(A4[panel.number()==10],A4,col = 'lightblue',alpha=0.4)
       panel.xblocks(B4[panel.number()==11],B4,col = 'lightblue',alpha=0.4)
       panel.xblocks(D4[panel.number()==12],D4,col = 'lightblue',alpha=0.4)
       panel.xblocks(A5[panel.number()==13],A5,col = 'lightblue',alpha=0.4)
       panel.xblocks(B5[panel.number()==14],B5,col = 'lightblue',alpha=0.4)
       panel.xblocks(D5[panel.number()==15],D5,col = 'lightblue',alpha=0.4)
       panel.xblocks(A6[panel.number()==16],A6,col = 'lightblue',alpha=0.4)
       panel.xblocks(B6[panel.number()==17],B6,col = 'lightblue',alpha=0.4)
       panel.xblocks(D6[panel.number()==18],D6,col = 'lightblue',alpha=0.4)
       panel.xblocks(A7[panel.number()==19],A7,col = 'lightblue',alpha=0.4)
       panel.xblocks(B7[panel.number()==20],B7,col = 'lightblue',alpha=0.4)
       panel.xblocks(D7[panel.number()==21],D7,col = 'lightblue',alpha=0.4)
       }) 
       
print(pl)
dev.off()



