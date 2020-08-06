require(tidyverse)

deep  <- read.csv(file = 'excel_text_files/csv_deepmala/deep_mala_markers.csv',stringsAsFactors = F)
head(deep)
str(deep)
colnames(deep)[2] <- 'rs.'


load(file = 'Rda/numeric_hapmap/gbs_data2.rda')
gbs_data2[1:6,1:6]

deep <- inner_join(deep,gbs_data2)
deep[1:6,1:9]

deep <- deep %>% gather(GID,marker,-1:-7) 
head(deep)

exce <- read.csv(file = 'excel_text_files/csv/excel_to_run.csv',stringsAsFactors = F)
head(exce)

deep <- deep %>% inner_join(select(exce,GID,trial))
head(deep)

deep$trial <- gsub("ES[1-7]$",          "ES1-7",    deep$trial)
deep$trial <- gsub("ES[8-9]$|ES1[0-3]$",  "ES8-13",  deep$trial)
deep$trial <- gsub("ES1[4-9]$|ES2[0-5]$", "ES14-25",  deep$trial) 
deep$trial <- gsub("ES2[6-9]$|ES3[0-8]$", "ES26-38",  deep$trial)
deep$trial <- gsub("SA[1-9]$|SA1[0-2]$",  "SA1-12",   deep$trial)
deep$trial <- gsub("SA1[3-9]$|SA2[0-5]$", "SA13-25",  deep$trial)

head(deep)

oi <- deep %>% group_by(trial, marker, rs.) %>% tally()
oi
oi <- oi %>% group_by(trial, rs.) %>% mutate(prop = prop.table(n))
oi
oi <-oi %>%
     group_by(trial,rs.) %>%
     mutate(prop2 = ifelse(marker != 1, prop + prop[marker == 1]/2, prop))

oi <- oi %>% mutate(prop3 = ifelse(is.na(prop2),prop,prop2)) 
oi <- oi %>% select(-c(prop,prop2))
oi <- oi %>% rename(prop=prop3)


deep  <- read.csv(file = 'excel_text_files/csv_deepmala/deep_mala_markers.csv',stringsAsFactors = F)
head(deep)
colnames(deep)[2] <- 'rs.'
deep <- inner_join(deep,gbs_data2)
deep[1:6,1:9]

oi2 <- oi %>% inner_join(select(deep,BASE,rs.,alleles))
oi2

write.csv(oi2,file = 'excel_text_files/csv_deepmala/oi2.csv',row.names = F)
oi2 <- read.csv(file = 'excel_text_files/csv_deepmala/oi2.csv',stringsAsFactors = F)

oi2 <- oi2 %>%
              group_by(rs.) %>%
                                filter(
                                       case_when(BASE==substr(alleles,1,1)~marker==2,
                                       BASE==substr(alleles,3,3)~marker==0)
                                       )
oi2
View(oi2)
# checking the variables 
str(oi2)
oi2$trial <- factor(oi2$trial)
oi2$trial <- factor(oi2$trial,levels(oi2$trial)[c(1,4,2,3,5,6)])
levels(oi2$trial)
oi2 <- oi2 %>% ungroup()

#### plotting

oi3 <- oi2 %>%

   # create appropriate x-axis values based on trial values
   mutate(x = case_when(trial == "ES1-7" ~ 1,
                        trial == "ES8-13" ~ 2,
                        trial %in% c("ES14-25", "SA1-12") ~ 3,
                        trial %in% c("ES26-38", "SA13-25") ~ 4,
                        TRUE ~ 0)) %>%
   
   # expand data frame by repeating the last point before divergence
   # for each rs. facet
   group_by(rs.) %>%
   mutate(last.point.before.divergence = x == max(x[x <= 2])) %>%
   ungroup() %>% 
   slice(c(1:n(),
           which(last.point.before.divergence))) %>% 
   
   # create group for line
   group_by(rs., x) %>%
   arrange(trial) %>% 
   mutate(group = seq(1, n())) %>%
   ungroup()
oi3



############### plotting

cust_theme <- theme_bw() + theme( panel.spacing =  unit(0, "lines"), 
                                  panel.border = element_rect(size = 0.25, color = "black"), 
                                  panel.grid = element_blank(), 
                                  axis.text.x = element_text(angle = 45,  hjust = 1),
                                  legend.title = element_blank(),
                                  axis.title.x=element_blank())



reto <-ggplot(data=oi3, aes(x = x, y = prop,
                     linetype = factor(group,labels = c('ESWYT','SAWYT')),
                     color=factor(group),
                     group = group)) +
       guides(color=FALSE)+
       scale_color_hue(l=40, c=35)+
       geom_point(size=1) + 
       geom_line(size=1) + 
       cust_theme + 
       facet_wrap( .~ rs., nrow =7)+
       scale_x_continuous(breaks = seq(1, 4),
                          labels = 
                          c("ES1-7","ES8-13",
                           "ES14-25 SA1-12",
                           "ES26-38 SA13-25")) +
       labs(y='Frequency')+
       expand_limits(y=c(0,1))


pdf('ggplot2.pdf',width=12,height = 10)
print(reto)
dev.off()




c <-   ggplot(oi3,aes(x = x, y = prop,linetype = factor(group,labels = c('ESWYT','SAWYT')),
                      color=factor(group)
                      )) +
       guides(color=FALSE)+
       scale_color_hue(l=40, c=35)+
       geom_point() +
       geom_line ( data = . %>%
                   group_by(rs., group) %>%
                   summarise(x1 = list(spline(x, prop, n = 50, method = "natural")[["x"]]),
                             y1 = list(spline(x, prop, n = 50, method = "natural")[["y"]])) %>%
                   
                   tidyr::unnest(),
                   aes(x = x1, y = y1)
                  ) +
       facet_wrap(.~rs.,nrow = 7) +
       scale_x_continuous(breaks = seq(1, 4),
                          labels = 
                          c("ES1-7","ES8-13",
                           "ES14-25 SA1-12",
                           "ES26-38 SA13-25")) +
       cust_theme +
       labs(y='Frequency')+
       expand_limits(y=c(0,1))
   

c

pdf('ggplot.pdf',width=12,height = 10)
print(c)
dev.off()






################# old plots


### function for x axis labels

addline_format <- function(x,...){
   gsub('\\s','\n',x)
}


old_version <-   ggplot(oi3[grep('S5',oi3$rs.),],aes(x = x, y = prop, color = rs., 
                                           linetype = factor(group,labels = c('ESWYT','SAWYT')))) +
   geom_point() +
   geom_line(data = . %>%
                group_by(rs., group) %>%
                summarise(x1 = list(spline(x, prop, n = 50, method = "natural")[["x"]]),
                          y1 = list(spline(x, prop, n = 50, method = "natural")[["y"]])) %>%
                
                tidyr::unnest(),
             aes(x = x1, y = y1)) +
   facet_grid(.~rs.) +
   scale_x_continuous(breaks = seq(1, 4),
                      labels = addline_format(c("ES1-7", "ES8-13",
                                                "ES14-25 SA1-13",
                                                "ES26-38 SA14-25"))) +
   labs(subtitle="Favorable allele over time", 
        y="allele frequency", 
        x="Groups", 
        title="Yield QTL markers",
        col='markers',
        linetype='Environment')+
   theme(axis.text.x = element_text(hjust = 0.7,size = 7.5,angle = 45),panel.spacing.x = unit(6,'mm')
         ,legend.position = 'none')


old_version_2 <- ggplot(oi3[grep('S1',oi3$rs.),]
                 ,aes(x = x, y = prop,
                      linetype = factor(group,labels = c('ESWYT','SAWYT')))) +
   geom_point() +
   geom_line(   data = . %>%
                   group_by(rs., group) %>%
                   summarise(x1 = list(spline(x, prop, n = 50, method = "natural")[["x"]]),
                             y1 = list(spline(x, prop, n = 50, method = "natural")[["y"]])) %>%
                   
                   tidyr::unnest(),
                aes(x = x1, y = y1)) +
   facet_grid(.~rs.) +
   scale_x_continuous(breaks = seq(1, 4),
                      labels = addline_format(
                         c("ES1-7", "ES8-13",
                           "ES14-25 SA1-13",
                           "ES26-38 SA14-25"))) +
   theme(legend.position = 'none',
         axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         panel.spacing = unit(0,'lines') 
   )
