require(here)

df_gbs = data.frame(read.delim(here('hapmap','Beagle.hmp'),header=T,stringsAsFactors = F))
class(df_gbs)
str(df_gbs)

# first eleven rows are not needed for the numerical conversion, the first 5 columns we will need later.
df_col = df_gbs[,1:4]
df_num = df_gbs[,-1:-11]

df_col[1:5,1:4]
df_num[1:5,1:5]

# Numeric encoding to transform  hapmap format into a SNP allele count matrix
# 1 to all heterozygotic sites and NA for missing data 

# for missing value
df_num[df_num=='N'] = NA
df_num[df_num=='.'] = NA
df_num[df_num=='-'] = NA

# for tri-allelic sites
df_num[df_num=='B'] = NA
df_num[df_num=='D'] = NA
df_num[df_num=='H'] = NA
df_num[df_num=='V'] = NA

# for heterozygous sites
df_num[df_num=='R'] = 1
df_num[df_num=='Y'] = 1
df_num[df_num=='S'] = 1
df_num[df_num=='W'] = 1
df_num[df_num=='K'] = 1
df_num[df_num=='M'] = 1


# Let's save the reference and alternative allele form column alleles in df_col into different objects

ref  <- substr(df_col$alleles,1,1)
alt  <- substr(df_col$alleles,3,3)

# Create a function to evaluate if the allele is the reference or the alternative and run sapply

f = function(x){
  x = ifelse(x==ref,2,x)
  x = ifelse(x==alt,0,x)
  return(x)
}

df_num <- sapply(df_num,f,simplify = TRUE)

# transform the the characters into numbers
mode(df_num) = 'numeric'
df_num[1:5,1:5]

# checking for na values
df_num[is.na(df_num)]

df_num[is.nan(df_num)]
df_num[is.null(df_num)]


# Calculate minor allele-frequency

MAF <- apply(df_num, 1, function(x) 1-mean(x)/2)
length(MAF) == nrow(df_num)

# Proportion of MAF >0.05 = 62%
table(MAF>0.05)
prop.table(table(MAF>0.05))

# Subset df_num rows by the 62% boolean True values which are greater than 0.05,
# and also subset df_col the remaining dataframe so as to have the complete df.

dim(df_num)
df_num <- df_num[MAF>0.05,]
dim(df_num)

dim(df_col)
df_col <- df_col[MAF>0.05,]
dim(df_col)

df <- cbind(df_col,df_num)
df[1:5,1:10]

save(df,file = here('hapmap','df.rda'))

