# load packages
require(here)
require(dplyr)
require(gtools)

#load objects
load(here('PCA','df.rda'))
load(here('PCA','trial_df.rda'))

# Glimpse of the data
df[1:7,1:7]
head(trial_df)

# We need to transform df columns into rows
snp <- data.frame(Genotype=names(df)[-1:-4],t(df[-1:-4]),row.names = NULL)
dim(snp)

# rename the columns with the locus information
names(snp)[-1] <- paste0(df$chrom,'-',df$pos)

# bind with the trial_df group column
snp <- cbind(group=trial_df$group,trial=trial_df$trial,snp)

# Now we need to transform into hierfstat format
transform <- snp[,-1:-3]
transform[1:7,1:7]
# coerce to character and duplicate string
transform <- sapply(transform, strrep,2)
transform[1:7,1:7]
# replace 11 for 12
transform <- replace(transform,transform=='11','12')
transform[1:7,1:7]
# replace 00 for 11
transform <- replace(transform,transform=='00','11')
transform[1:7,1:7]
# coerce to numeric
mode(transform) <- 'numeric'
str(transform)

# let's bind back the columns
fst_hierfstat <- cbind(snp[,1:3],transform)
str(fst_hierfstat)
fst_hierfstat[1:7,1:7]

# Create smaller populations to identify FST trends
fst_hierfstat <- fst_hierfstat %>% 
  mutate(pop = case_when(
    trial == "ES1" | trial == "ES2" | trial == "ES3" | trial == "ES4" ~ "ES1-4",
    trial == "ES5" | trial == "ES6" | trial == "ES7" ~ "ES5-7",
    trial == "ES8" | trial == "ES9" | trial == "ES10" ~ "ES8-10",
    trial == "ES11" | trial == "ES12" | trial == "ES13" ~ "ES11-13",
    trial == "ES14" | trial == "ES15" | trial == "ES16" ~ "ES14-16",
    trial == "ES17" | trial == "ES18" | trial == "ES19" ~ "ES17-19",
    trial == "ES20" | trial == "ES21" | trial == "ES22" ~ "ES20-22",
    trial == "ES23" | trial == "ES24" | trial == "ES25" ~ "ES23-25",
    trial == "ES26" | trial == "ES27" | trial == "ES28" ~ "ES26-28",
    trial == "ES29" | trial == "ES30" | trial == "ES31" ~ "ES29-31",
    trial == "ES32" | trial == "ES33" | trial == "ES34" ~ "ES32-34",
    trial == "ES35" | trial == "ES36" | trial == "ES37"  | trial == "ES38" ~ "ES35-38",
    trial == "SA1" | trial == "SA2" | trial == "SA3" ~ "SA1-3",
    trial == "SA4" | trial == "SA5" | trial == "SA6" ~ "SA4-6",
    trial == "SA7" | trial == "SA8" | trial == "SA9" ~ "SA7-9",
    trial == "SA10" | trial == "SA11" | trial == "SA12" ~ "SA10-12",
    trial == "SA13" | trial == "SA14" | trial == "SA15" ~ "SA13-15",
    trial == "SA16" | trial == "SA17" | trial == "SA18" ~ "SA16-18",
    trial == "SA19" | trial == "SA20" | trial == "SA21" ~ "SA19-21",
    trial == "SA22" | trial == "SA23" | trial == "SA24" | trial == "SA25" ~ "SA22-25")) %>%
  select(Genotype,trial,group,pop,everything())

# Create a column with pop as factor and sort them in temporal order, this will be useful for plots        
fst_hierfstat$pop_factor <- factor(fst_hierfstat$pop)
fst_hierfstat$pop_factor <- factor(fst_hierfstat$pop,levels = mixedsort(levels(fst_hierfstat$pop_factor)))
levels(fst_hierfstat$pop_factor)

# Rearrange columns
fst_hierfstat <- fst_hierfstat %>% select(Genotype,trial,group,pop,pop_factor,everything())
fst_hierfstat[1:7,1:7]

# And save the object
save(fst_hierfstat,file = here('FST','fst_hierfstat.rda'))
save(snp,file = here('adegenet_dapc','snp.rda'))

