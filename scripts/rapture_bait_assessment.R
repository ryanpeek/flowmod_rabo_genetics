# RAPTURE bait assessment

# Thu Jul  6 16:18:23 2017 ------------------------------

library(tidyverse)
library(magrittr)
#source("https://bioconductor.org/biocLite.R") 
#biocLite("Biostrings")
suppressPackageStartupMessages(library(Biostrings))

# GET DATA ----------------------------------------------------------------

## USE SFTP TO GRAB FILE:

# cd Documents/github/rabo_genomics/data/RAPTURE
# sftp -P 2022 rapeek@farm.cse.ucdavis.edu
# cd projects/rana_rapture/fastq 
# get bamlist.flt.mafs.gz

# angsd data on all reads
dat <- read_delim("data/bamlist.flt.mafs.gz", delim = "\t")
head(dat)
dim(dat)

# filter to position 50
dat50 <- dat %>% select(chromo, position, nInd) %>% filter(position==50) 

# baits
baits<-fasta.index(filepath = "data/baits_filt5_120.fasta")
head(baits) ## can use gunzip or gzip to get file directly from .gz format.

baits <- baits %>% select(desc) %>% dplyr::rename("chromo"=desc)
head(baits)
dim(baits)

library(readr)
write_delim(baits, path = "data_output/bait_list")

# FILTER TO BAITS ONLY  ----------------------------------------------------

# use an inner_join to match to baits

df <- inner_join(dat50, baits, by="chromo")

mean(df$nInd) # 661
median(df$nInd) # 720

# random sample of 100 baits around median number:

#df_med <- df %>% filter(nInd>720, nInd<770) %>% sample_n(100) %T>% glimpse

# df_med %>%
#   ggplot(data=., aes(x=nInd)) + 
#   geom_histogram(binwidth=5, color="gray60") + 
#   geom_vline(xintercept=median(df$nInd), col="red", lty=2, lwd=1, alpha=0.6)+
#   labs(title="Random Sample of Baits Near Median: BAITS (n=100)",
#        x="# Individuals Sharing Bait")

# All baits at pos 50
df %>%
  ggplot(data=., aes(x=nInd)) + 
  geom_histogram(binwidth=25, color="gray60") + 
  geom_vline(xintercept=median(df$nInd), col="red", lty=2, lwd=1, alpha=0.6)+
  geom_vline(xintercept=500, col="blue", lty=2, lwd=1, alpha=0.6)+
  geom_vline(xintercept=1000, col="blue", lty=2, lwd=1, alpha=0.6)+
  labs(title="Number of Individuals at Position 50: BAITS (n=9,377)",
       x="# Individuals Sharing Bait")

# all non-baits
dfNB <- anti_join(dat50, baits, by="chromo")

dfNB %>% ggplot(data=., aes(x=nInd)) + 
  geom_histogram(binwidth=25, color="gray60") + 
  geom_vline(xintercept=mean(dfNB$nInd), col="red", lty=2, lwd=1, alpha=0.6)+
  labs(title="Number of Individuals at Position 50: NON-BAITS (n=23,907)",
       x="# Individuals Sharing Bait")
  
# all chromos at position 50
dat50 %>% 
  ggplot(data=., aes(x=nInd)) + 
  geom_histogram(binwidth=25, color="gray60") + 
  geom_vline(xintercept=mean(dat50$nInd), col="red", lty=2, lwd=1, alpha=0.6)+
  labs(title="Number of Individuals at Position 50: ALL CHROMOS (n=33,284)", 
       x="# Individuals Sharing Bait")

#ggsave(filename = "figs/rapture/mean_reads_by_collection_sampletype.png", width = 11, height=8.5, units = "in", dpi=200)



# SAMPLE BAITS ------------------------------------------------------------

# 1. filter to >0.1 minor allele freq (unknownEM) and < 0.3, and position should be < 84;
# 2. Then find baits that are close to "median" range (500 to 1000)

head(dat)
dat_filt<- dat %>% filter(unknownEM>0.1, unknownEM<0.3, position<84) %>% inner_join(., baits, by="chromo") %>% filter(nInd>500, nInd<1000)

dat_filt %>% ggplot(.) + geom_histogram(aes(x=nInd),binwidth = 10, color="gray70")

# randomly sample one SNP per bait for plots:
bait_snp_filt <- dat_filt %>% group_by(chromo) %>% sample_n(1) %>% as.data.frame() %>% 
  select(chromo, position)

#write_delim(bait_snp_filt, "data_output/RAPTURE/bait_snp_filt",col_names = F)

# need to index on cluster
# angsd sites index bait_snp_filt

# sbatch -t 48:00:00 --mem=60G 02_angsd_sites_cov.sh bamlist.flt bait_sub_cov.flt

# count how many NN (uncalled genos in each column)
# count how many reads (x-axis)
# plot NN/total genos (total number of rows in geno file). 
# then split up by plates (each col will be each line in bamlist), plot for each plate

# BAIT GENOS --------------------------------------------------------------

# bamlist.names
bams <- read_tsv("data/bamlist.names", col_names=F)
head(bams)
bams$plateID<-substr(bams$X1, 9,14)
bams$wellID<-substr(bams$X1, 21,28)
bams$sampID<-paste0(bams$plateID,"_",bams$wellID)
length(bams$sampID)

# add col names
bamnames<- c("chromo","pos",bams$sampID)
bamnames

# barcode_wells
#bams <- read_tsv("data/RAPTURE/bamlist.names", col_names=F)

# angsd data on reads
#dat2 <- read_delim("data/RAPTURE/bait_sub_cov.geno.gz", col_names = F, delim = "\t", col_types = c("c","i","l"))
dat2 <- read.delim(file = "data/bait_sub_cov.geno.gz", header = F, sep = "\t")
head(dat2)
dat2 <- dat2[,-1635]
colnames(dat2)<-bamnames
str(dat2)

# count NN's per col

nntot<-gather(dat2, key, value, -pos) %>%
  group_by(key, value) %>%
  tally %>% 
  spread(value, n, fill = 0) %>% 
  select(1:14) %>% 
  as.data.frame()


# JOIN WITH METADATA ------------------------------------------------------

# plate metadata
metadata <- read_csv(file = "data/Plates_RADseq-RAPTURE.csv") %>% as.data.frame()
head(metadata)

# get reads and mapped counts
counts <- read.table("data/rapture_run.txt", header = F, sep = "\t",
                     col.names = c("Seq","reads","mapped", "paired", "unique"))
head(counts)

counts$plateID<-substr(counts$Seq, 9,14)
counts$wellID<-substr(counts$Seq, 21,28)
counts$sampID<-paste0(counts$plateID,"_",counts$wellID)

# add reads
nntot_counts <- left_join(nntot, counts, by=c("key"="sampID")) %>% 
  select(key, NN:wellID)
names(nntot_counts)

# add prop called:
nntot_counts <- nntot_counts %>% 
  mutate(prop_called= (1-(NN/2512))) %>% 
  filter(!is.na(Seq))

summary(nntot_counts)

# add plateMetadata
nntot_counts <- left_join(nntot_counts, metadata, by=c("plateID"="PlateBarcode"))
names(nntot_counts)

# plot reads by prop_called to see how each plate performed
ggplot() + geom_point(data=nntot_counts, aes(x=reads, y=prop_called, color=Collection, shape=SampleType), alpha=0.7) + geom_vline(xintercept = 500000, col="red", lty=2)+
  geom_smooth(data=nntot_counts, aes(x=reads, y=prop_called, color=Collection),
              se=FALSE, method="gam", formula=y~s(log(x)))

# as above but for a single plate
ggplot() + geom_point(data=nntot_counts[nntot_counts$Collection=="WSU1",], aes(x=reads, y=prop_called, color=Collection), alpha=0.7) + 
  geom_vline(xintercept = 200000, col="red", lty=2) + 
  geom_smooth(data=nntot_counts[nntot_counts$Collection=="WSU1",], aes(x=reads, y=prop_called, color=Collection),
              se=FALSE, method="gam", formula=y~s(log(x)))

# proportion of reads currently above 500,000
nntot_500k <- nntot_counts %>% 
  filter(reads>500000) 

# plot reads by prop_called to see how each plate performed
ggplot() + geom_point(data=nntot_500k, aes(x=reads, y=prop_called, color=Collection, shape=SampleType), alpha=0.7) + geom_vline(xintercept = 500000, col="red", lty=2)

# Number of samples above the 500k read threshold
nntot_500k %>% group_by(plateID, RAD_ID) %>% 
  tally() %>% 
  ggplot(.) + geom_point(aes(x=RAD_ID, y=n, color=RAD_ID), size=4)

# proportion of reads currently above 200,000
nntot_200k <- nntot_counts %>% 
  filter(reads>200000) 

head(nntot_200k)

# write to a new bamlist to run PCA with
wsuRABO<-filter(nntot_200k, RAD_ID=="RAD-220" | RAD_ID=="RAD-221" | RAD_ID=="RAD-222") %>% select(Seq)

wsuRABO_b <- str_replace(string = wsuRABO$Seq, pattern = "_RA_", replacement = "_RB_")
head(wsuRABO_b)
wsuRABO$b<-wsuRABO_b
head(wsuRABO)
write_delim(wsuRABO, path="data_output/bamlist_wsu_rabo.txt", col_names = F)


# plot reads by prop_called to see how each plate performed
ggplot() + geom_point(data=nntot_200k, aes(x=reads, y=prop_called, color=Collection, shape=SampleType), alpha=0.7) + geom_vline(xintercept = 200000, col="red", lty=2)

# Number of samples above the 200k read threshold
nntot_200k %>% group_by(plateID, RAD_ID) %>% 
  tally() %>% 
  ggplot(.) + geom_point(aes(x=RAD_ID, y=n, color=RAD_ID), size=4)

# summarize the reads, to see how it performed
nntot_counts %>% filter(plateID=="CGGAAT") %>% 
  summarize(reads_sum=sum(reads)/2) # all R1 reads...this number *17 plates is = ~200,000,000 which is half a lane, so makes sense

# now let's calc how much additional seq will be req to meet 95% threshold
test2 <- nntot_counts %>% group_by(plateID) %>% 
  mutate(readsW50mill = ((reads / sum(reads))*50000000) + reads, 
         n_cutoff = if_else(readsW50mill > 500000, 1, 0))
head(test2)

ggplot() + geom_histogram(data=test2, aes(readsW50mill, fill=n_cutoff), binwidth = 25000) + #ylim(0,20)+
  geom_vline(xintercept = 500000, col="red", lty=2) + 
  facet_wrap(nrow = 5, facets = "plateID", scales = "free_y")

# now summarize number of individuals per plate that would work above the 500k threshold using our calc with 500m reads
plateThresh <- test2 %>% group_by(plateID) %>% 
  summarize(nInds=sum(n_cutoff)) %>% 
  left_join(., metadata, by=c("plateID"="PlateBarcode"))

library(ggrepel)
library(viridis)

ggplot() + 
  geom_point(data=plateThresh, aes(x=nInds, y=RAD_ID, fill=nInds/Samples), size=4, pch=21) + 
  geom_label_repel(data=plateThresh, aes(x=nInds, y=RAD_ID, label=Samples)) +
  scale_fill_viridis()

# write out nntot
write_delim(nntot_counts, path = "data_output/nntot_counts.txt")

# for each plate, how many more reads will be required to reach 500000?