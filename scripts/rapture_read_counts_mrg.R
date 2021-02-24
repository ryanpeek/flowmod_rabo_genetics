# quantify RAPTURE seqs

# Thu Sep  7 00:14:15 2017 ------------------------------

# LOAD DATA & LIBRARIES ---------------------------------------------------

library(tidyverse)

# READ IN DATA ------------------------------------------------------------

counts <- read.table("data/fastq_counts_mrg.txt", header = F, sep = "\t",
                    col.names = c("Seq","reads","mapped", "paired", "unique"))
head(counts)

metadata <- read_csv(file = "data/Plates_RADseq-RAPTURE.csv") %>% as.data.frame()
head(metadata)

# SPLIT OUT BARCODES FROM PLATES ------------------------------------------

# get plate ID
counts$plateID<-substr(counts$Seq, 9, 14 )
head(counts)
counts$wellID<-substr(counts$Seq, 21,28)
counts$sampID<-paste0(counts$plateID,"_",counts$wellID)

# get plate ID meta
metadata$plateID<-substr(metadata$Plates, 12, 17 )
head(metadata)

# JOIN Plates ------------------------------------------------------------

df <- left_join(counts, metadata[,c(2:13)], by="plateID")

# join with metadata
head(df)

# rename a few cols
#save(df, file = "data/RAPTURE/plate_metadata.rda")

# SUMMARIZE ---------------------------------------------------------------

# need to calc CV by plate
dfCV <- df %>% select(plateID, reads:unique) %>% 
  group_by(plateID) %>% 
  summarize_all(funs(mean, sd)) %>% 
  mutate(readsCV=reads_mean/reads_sd,
         mappedCV=mapped_mean/mapped_sd,
         pairedCV=paired_mean/paired_sd,
         uniqueCV=unique_mean/unique_sd) %>% 
  left_join(., select(metadata, plateID, RAD_ID, Collection, SampleType, Samples), by="plateID")


# PLOTS -------------------------------------------------------------------

library(ggrepel)
library(viridis)

#options("scipen"=8, "digits"=4) # set scientific notation lower

dfCV %>% mutate(propmean=unique_mean/reads_mean) %>% ggplot(.) + geom_point(aes(x=RAD_ID, y=propmean, color=Collection, shape=SampleType), size=4) + theme(axis.text.x = element_text(angle = 90, vjust = .5, size=12))

# READS -------------------------------------------------------------------

# mean with +/- SD
ggplot() + 
  geom_point(data=dfCV, aes(x=RAD_ID, y=reads_mean, color=Collection, shape=SampleType), size=4.5) +
  geom_pointrange(data=dfCV, aes(x=RAD_ID, y=reads_mean, ymin=reads_mean, ymax=reads_mean + reads_sd, color=Collection, shape=SampleType)) + ylim(0, 3e6) +
  #scale_y_continuous(breaks=c(seq(1e5, 3e5, 25000)))+
  guides(shape=guide_legend(title="Sample Types"))+
  theme(axis.text.x = element_text(angle = 90, vjust = .5, size=12),
        legend.key=element_blank())+ 
  labs(x="", y="Mean Reads", title="Mean Reads")

ggsave(filename = "figs/rapture/mean_reads_by_collection_sampletype.png", width = 11, height=8.5, units = "in", dpi=200)

# CV
ggplot() + 
  geom_point(data=dfCV, aes(x=RAD_ID, y= readsCV*100, color=Collection, shape=SampleType), size=4.5) +
  guides(shape=guide_legend(title="Sample Types"))+
  theme(axis.text.x = element_text(angle = 90, vjust = .5, size=12),
        legend.key=element_blank())+
  geom_hline(yintercept = 100, color="red", lty=2)+
  ylim(c(20,160))+
  labs(x="", y="%CV of Total Reads", title="%CV of Total Reads")

ggsave(filename = "figs/rapture/CV_mean_reads_by_collection_sampletype.png", width = 11, height=8.5, units = "in", dpi=200)


# MAPPED ------------------------------------------------------------------

# mean
ggplot() + 
  geom_point(data=dfCV, aes(x=RAD_ID, y=mapped_mean, color=Collection, shape=SampleType), size=4.5) +
  #scale_y_continuous(breaks=c(seq(1e5, 3e5, 25000)))+
  guides(shape=guide_legend(title="Sample Types"))+
  theme(axis.text.x = element_text(angle = 90, vjust = .5, size=12),
        legend.key=element_blank())+ 
  labs(x="", y="Mapped Reads", title="Mapped Reads")

ggsave(filename = "figs/rapture/mapped_reads_by_collection_sampletype.png", width = 11, height=8.5, units = "in", dpi=200)

# CV
ggplot() + 
  geom_point(data=dfCV, aes(x=RAD_ID, y= mappedCV*100, color=Collection, shape=SampleType), size=4.5) +
  guides(shape=guide_legend(title="Sample Types"))+
  theme(axis.text.x = element_text(angle = 90, vjust = .5, size=12),
        legend.key=element_blank())+
  geom_hline(yintercept = 100, color="red", lty=2)+
  ylim(c(20,160))+
  labs(x="", y="%CV of Mapped Reads", title="%CV of Mapped Reads")

ggsave(filename = "figs/rapture/CV_mapped_reads_by_collection_sampletype.png", width = 11, height=8.5, units = "in", dpi=200)


# PAIRED ------------------------------------------------------------------

# mean
ggplot() + 
  geom_point(data=dfCV, aes(x=RAD_ID, y=paired_mean, color=Collection, shape=SampleType), size=4.5) +
  #scale_y_continuous(breaks=c(seq(1e5, 3e5, 25000)))+
  guides(shape=guide_legend(title="Sample Types"))+
  theme(axis.text.x = element_text(angle = 90, vjust = .5, size=12),
        legend.key=element_blank())+ 
  labs(x="", y="Prop Paired", title="Paired Reads")
ggsave(filename = "figs/rapture/prop_paired_by_collection_sampletype.png", width = 11, height=8.5, units = "in", dpi=200)

# CV
ggplot() + 
  geom_point(data=dfCV, aes(x=RAD_ID, y= pairedCV*100, color=Collection, shape=SampleType), size=4.5) +
  guides(shape=guide_legend(title="Sample Types"))+
  theme(axis.text.x = element_text(angle = 90, vjust = .5, size=12),
        legend.key=element_blank())+
  geom_hline(yintercept = 100, color="red", lty=2)+
  ylim(c(20,160))+
  labs(x="", y="%CV of Prop Paired Reads", title="%CV of Paired Reads")

ggsave(filename = "figs/rapture/CV_paired_reads_by_collection_sampletype.png", width = 11, height=8.5, units = "in", dpi=200)


# UNIQUE ------------------------------------------------------------------

# mean
ggplot() + 
  geom_point(data=dfCV, aes(x=RAD_ID, y=unique_mean, color=Collection, shape=SampleType), size=4.5) +
  #scale_y_continuous(breaks=c(seq(1e5, 3e5, 25000)))+
  guides(shape=guide_legend(title="Sample Types"))+
  theme(axis.text.x = element_text(angle = 90, vjust = .5, size=12),
        legend.key=element_blank())+ 
  labs(x="", y="Prop Paired", title="Paired Reads")

ggsave(filename = "figs/rapture/prop_paired_by_collection_sampletype.png", width = 11, height=8.5, units = "in", dpi=200)

# CV
ggplot() + 
  geom_point(data=dfCV, aes(x=RAD_ID, y= uniqueCV*100, color=Collection, shape=SampleType), size=4.5) +
  guides(shape=guide_legend(title="Sample Types"))+
  theme(axis.text.x = element_text(angle = 90, vjust = .5, size=12),
        legend.key=element_blank())+
  geom_hline(yintercept = 100, color="red", lty=2)+
  ylim(c(20,160))+
  labs(x="", y="%CV of Unique Reads", title="%CV of Unique Reads")

ggsave(filename = "figs/rapture/CV_unique_reads_by_collection_sampletype.png", width = 11, height=8.5, units = "in", dpi=200)



# RATIOS ------------------------------------------------------------------

# mapped/reads

ggplot() + 
  geom_point(data=dfCV, aes(x=RAD_ID, y=mapped_mean/reads_mean, color=Collection, shape=SampleType), size=4.5) +
  #scale_y_continuous(breaks=c(seq(1e5, 3e5, 25000)))+
  guides(shape=guide_legend(title="Sample Types"))+
  theme(axis.text.x = element_text(angle = 90, vjust = .5, size=12),
        legend.key=element_blank())+ 
  labs(x="", y="Mapped/Reads", title="Mapped/Reads")

ggsave(filename = "figs/rapture/ratio_mapped-to-reads_by_collection_sampletype.png", width = 11, height=8.5, units = "in", dpi=200)

# pp/mapped

ggplot() + 
  geom_point(data=dfCV, aes(x=RAD_ID, y=paired_mean/mapped_mean, color=Collection, shape=SampleType), size=4.5) +
  #scale_y_continuous(breaks=c(seq(1e5, 3e5, 25000)))+
  guides(shape=guide_legend(title="Sample Types"))+
  theme(axis.text.x = element_text(angle = 90, vjust = .5, size=12),
        legend.key=element_blank())+ 
  labs(x="", y="Paired/Mapped", title="Paired/Mapped")

ggsave(filename = "figs/rapture/ratio_pairs-to-mapped_by_collection_sampletype.png", width = 11, height=8.5, units = "in", dpi=200)

# uniq/pp

ggplot() + 
  geom_point(data=dfCV, aes(x=RAD_ID, y=unique_mean/paired_mean, color=Collection, shape=SampleType), size=4.5) +
  guides(shape=guide_legend(title="Sample Types"))+
  theme(axis.text.x = element_text(angle = 90, vjust = .5, size=12),
        legend.key=element_blank())+ 
  labs(x="", y="unique_mean/Paired", title="unique_mean/Paired")
ggsave(filename = "figs/rapture/ratio_uniq-to-pairs_by_collection_sampletype.png", width = 11, height=8.5, units = "in", dpi=200)

