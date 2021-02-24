# read counts for reg paper

counts <- read.table("data_output/sierra_samples_v3_fq_counts.txt")


head(counts)

colnames(counts) <- c("Seq", "tot_reads", "paired_reads", "prop_pairs")

summary(counts)

