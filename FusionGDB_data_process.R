library(dplyr)
library(tidyr)

dd <- read.delim("combinedFGDB2genes_information_genesymbol_arranged_real.txt", header = F)
#unique(dd$V1)
colnames(dd) <- c("database", "cohort", "ID", "gene1", "chr1", "bp1", "g/l1", "gene2", "chr2", "bp2", "g/l2")
dd <- dd[! is.na(dd$chr1),]
dd <- dd[! is.na(dd$bp1),]

dd <- dd[! is.na(dd$chr2),]
dd <- dd[! is.na(dd$bp2),]

data1 <- data.frame(chr= noquote(dd$chr1), st1= dd$bp1, st2= dd$bp1)
#data2 <- data.frame(chr= noquote(dd$chr2), st1= dd$bp2, st2= dd$bp2)

write.table(data1, file = "FusionGDB_conversion1.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#write.table(data2, file = "FusionGDB_conversion2.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#sd <- read.table("FusionGDB_conversion.BED")

##### load converted data
faulty1 <- read.table("failure1.txt")
faulty1 <- faulty1[,1:2]

#faulty2 <- read.table("failure2.txt")
#faulty2 <- faulty2[,1:2]


##### filter out failures
filtered_df <- dd %>%
  anti_join(faulty1, by = c("chr1" = "V1", "bp1" = "V2"))

filtered1 <- data.frame(chr= noquote(filtered_df$chr2), st1= filtered_df$bp2, st2= filtered_df$bp2)

write.table(filtered1, file= "FusionGDB_conversion2.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

####second breakpoints filtered

faulty2 <- read.table("failure2.txt")
faulty2 <- faulty2[,1:2]

##### filter out failures
converted_filtered_df <- filtered_df %>%
  anti_join(faulty2, by = c("chr2" = "V1", "bp2" = "V2"))

converted_filtered_df$chr1 <-  gsub("chr","",converted_filtered_df$chr1)
converted_filtered_df$chr2 <-  gsub("chr","",converted_filtered_df$chr2)

chrs <- c(1:22, "X", "Y")

converted_filtered_df <- converted_filtered_df[converted_filtered_df$chr1 %in% chrs, ]
converted_filtered_df <- converted_filtered_df[converted_filtered_df$chr2 %in% chrs, ]

write.table(converted_filtered_df, "fusionGDB_filtered_final.txt", row.names = F)
#################
chitar <- converted_filtered_df[converted_filtered_df$database == "ChiTaRS5.0",]
non_chitar <- converted_filtered_df[converted_filtered_df$database != "ChiTaRS5.0",]

#### summary chr pairs

summary_chr <- matrix(0, nrow = 24, ncol = 24)
colnames(summary_chr) <- chrs
rownames(summary_chr) <- chrs

for(i in 1:nrow(converted_filtered_df)){
  ind1 <- which(colnames(summary_chr) == converted_filtered_df$chr1[i])
  ind2 <- which(colnames(summary_chr) == converted_filtered_df$chr2[i])
  if(ind1 > ind2){
  summary_chr[ind1,ind2] <- summary_chr[ind1,ind2] + 1 
  } else {
    summary_chr[ind2,ind1] <- summary_chr[ind2,ind1] + 1 
  }
}

write.csv(summary_chr, "FusionGDB_chr-matrix_summary.csv")
##### plot whole data for chrs
# Convert matrix to data frame
df <- as.data.frame(as.table(summary_chr))

# Load the ggplot2 package
library(ggplot2)

# Plot heatmap using ggplot2
# Define the order of X and Y axes
order_levels <- c("Y", "X", 22:1)

# Reorder factors
df$Var1 <- factor(df$Var1, levels = order_levels)
df$Var2 <- factor(df$Var2, levels = order_levels)

# Plot heatmap using ggplot2 with reordered axes
ggplot(df, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("#FFFFFF", "#D6E6F5", "#A3CDE7", "#6F9FD8", "#3C86CA", "#0057B8"), na.value = "white") +
  labs(title = "Heatmap Example",
       x = "X-axis Label",
       y = "Y-axis Label") +
  theme_minimal()

######until here
# rotate heatmap. count the cis and trans fusions and report them