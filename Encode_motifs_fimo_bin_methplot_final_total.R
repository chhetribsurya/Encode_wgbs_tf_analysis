library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(gtools)
library(gplots)
library(reshape)
library(data.table)
library(grid)

# weights <- ifelse(pcaOutput2$PC1 < -5 & abs(pcaOutput2$PC2) > 10, 2, 1)

####################

#HEATMAP.2 plots with clustering:

####################
output_dir <- "/home/surya/Dropbox/encode_3/fimo_motif_based_analysis/motifs_bin_methplot_total/heatmaps_motif_updownstream"
### Load the TF annoation file:
anno_df <- fread("~/Dropbox/TFs_Annotation_file.txt", sep="\t")

### 2kb up and downstream heatmap:
input_file <- "~/Dropbox/encode_3/motifs_methplot_total/all_tf_files/all_tfs_motif_bin_methplot_heatmap_data.bed"
read_file <- fread(input_file, sep="\t", header=TRUE)
read_file_merged <- merge(read_file, anno_df, by.x="tf_name", by.y="Target", all.x=TRUE)
#read_file_merged[,c("tf_name", "Category")]
read_file_merged$tf_name_edit <- paste(read_file_merged$tf_name, read_file_merged$annotation)

read_df <- as.data.frame(read_file_merged)
data <- read_df[,2:(ncol(read_df)-4)] ## Subtract last 4 columns
mat_data <- as.matrix(data)
rownames(mat_data) <- read_df$tf_name_edit

colnames(mat_data)
rownames(mat_data)


output_file_name <- paste0(output_dir, "/", "final_motifs_binned_methplot_2kb_up&downstream.pdf")       
pdf(output_file_name)

#custom_col <- c(green_col(1), red_col1(1))
#custom_col <- c(grey_col(1), red_col(1))
custom_col <- colorRampPalette(c("blue", "yellow", "red"))(10)

summary(read_df) #max_val = 325
#improved_col <- c("#808080",greenred(10) )
par(cex.main=0.80)
test <- heatmap.2(mat_data,
  # dendrogram = "none",
  # Rowv = TRUE, 
  Colv = FALSE,
  #scale="none",
  main = "Methylation distribution from the motif centres", 
  xlab = "Distance(bp) from the centre of motifs",
  ylab = "Transcription Factors",
  #col=greenred(10), 
  col=custom_col, 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.15,
  cexCol = 0.35,
  na.color="grey",
  RowSideColors=read_df$label,
  key.xlab="Methylation fraction"
  #symm=F,symkey=T,symbreaks=T,
  #breaks=c(-1,0.8,1.2,3)

  #breaks=c(c(0,0.9),seq(1,325,(325-1)/9))
  #scale="column", tracecol="#303030"
  )   
 
legend("topright", legend=c("Prom assoc","Enh assoc","CTCF assoc","Multiple assoc","Others"),
fill=c("red","orange","blueviolet","burlywood","green"), border=TRUE, bty="n", y.intersp = 0.9, cex=0.7) 

test
dev.off()





### 1kb up and downstream heatmap:

# input_file <- "~/Dropbox/encode_3/motifs_methplot_total/all_tf_files/all_tfs_motif_bin_methplot_heatmap_data.bed"
# read_file <- fread(input_file, sep="\t", header=TRUE)

# read_df <- as.data.frame(read_file)
# data <- read_df[,2:ncol(read_df)]
# mat_data <- as.matrix(data)
# column_ranges <- as.numeric(colnames(mat_data))
# column_select <- column_ranges[column_ranges >= -1000 & column_ranges <= 1000] 
# mat_data <- mat_data[,colnames(mat_data) %in% column_select]

# rownames(mat_data) <- read_df[,1]

# colnames(mat_data)
# rownames(mat_data)



### 1kb up and downstream heatmap:

### Load the TF annoation file:
anno_df <- fread("~/Dropbox/TFs_Annotation_file.txt", sep="\t")

### 1kb up and downstream heatmap:
input_file <- "~/Dropbox/encode_3/motifs_methplot_total/all_tf_files/all_tfs_motif_bin_methplot_heatmap_data.bed"
read_file <- fread(input_file, sep="\t", header=TRUE)
read_file_merged <- merge(read_file, anno_df, by.x="tf_name", by.y="Target", all.x=TRUE)
#read_file_merged[,c("tf_name", "Category")]
read_file_merged$tf_name_edit <- paste(read_file_merged$tf_name, read_file_merged$annotation)

read_df <- as.data.frame(read_file_merged)
data <- read_df[,2:(ncol(read_df)-4)] ## Subtract last 4 columns
mat_data <- as.matrix(data)

column_ranges <- as.numeric(colnames(mat_data))
column_select <- column_ranges[column_ranges >= -1000 & column_ranges <= 1000] 
mat_data <- mat_data[,colnames(mat_data) %in% column_select]

rownames(mat_data) <- read_df$tf_name_edit

colnames(mat_data)
rownames(mat_data)



output_file_name <- paste0(output_dir, "/", "final_motifs_binned_methplot_1kb_up&downstream.pdf")       
pdf(output_file_name)

#custom_col <- c(green_col(1), red_col1(1))
#custom_col <- c(grey_col(1), red_col(1))
custom_col <- colorRampPalette(c("blue", "yellow", "red"))(10)

summary(read_df) #max_val = 325
#improved_col <- c("#808080",greenred(10) )
par(cex.main=0.80)
test <- heatmap.2(mat_data,
  # dendrogram = "none",
  # Rowv = TRUE, 
  Colv = FALSE,
  #scale="none",
  main = "Methylation distribution from the motif centres", 
  xlab = "Distance(bp) from the centre of motifs",
  ylab = "Transcription Factors",
  #col=greenred(10), 
  col=custom_col, 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.15,
  cexCol = 0.35,
  na.color="grey",
  RowSideColors=read_df$label,
  key.xlab="Methylation fraction"
  #symm=F,symkey=T,symbreaks=T,
  #breaks=c(-1,0.8,1.2,3)

  #breaks=c(c(0,0.9),seq(1,325,(325-1)/9))
  #scale="column", tracecol="#303030"
  )    

 
legend("topright", legend=c("Prom assoc","Enh assoc","CTCF assoc","Multiple assoc","Others"),
fill=c("red","orange","blueviolet","burlywood","green"), border=TRUE, bty="n", y.intersp = 0.9, cex=0.7) 

test
dev.off()




### 0.6kb up and downstream heatmap:

# input_file <- "~/Dropbox/encode_3/motifs_methplot_total/all_tf_files/all_tfs_motif_bin_methplot_heatmap_data.bed"
# read_file <- fread(input_file, sep="\t", header=TRUE)

# read_df <- as.data.frame(read_file)
# data <- read_df[,2:ncol(read_df)]
# mat_data <- as.matrix(data)
# column_ranges <- as.numeric(colnames(mat_data))
# column_select <- column_ranges[column_ranges >= -600 & column_ranges <= 600] 
# mat_data <- mat_data[,colnames(mat_data) %in% column_select]

# rownames(mat_data) <- read_df[,1]

# colnames(mat_data)
# rownames(mat_data)


### 0.6kb up and downstream heatmap:

### Load the TF annoation file:
anno_df <- fread("~/Dropbox/TFs_Annotation_file.txt", sep="\t")

### 0.6kb up and downstream heatmap:
input_file <- "~/Dropbox/encode_3/motifs_methplot_total/all_tf_files/all_tfs_motif_bin_methplot_heatmap_data.bed"
read_file <- fread(input_file, sep="\t", header=TRUE)
read_file_merged <- merge(read_file, anno_df, by.x="tf_name", by.y="Target", all.x=TRUE)
#read_file_merged[,c("tf_name", "Category")]
read_file_merged$tf_name_edit <- paste(read_file_merged$tf_name, read_file_merged$annotation)

read_df <- as.data.frame(read_file_merged)
data <- read_df[,2:(ncol(read_df)-4)] ## Subtract last 4 columns
mat_data <- as.matrix(data)


column_ranges <- as.numeric(colnames(mat_data))
column_select <- column_ranges[column_ranges >= -600 & column_ranges <= 600] 
mat_data <- mat_data[,colnames(mat_data) %in% column_select]
rownames(mat_data) <- read_df$tf_name_edit

colnames(mat_data)
rownames(mat_data)



output_file_name <- paste0(output_dir, "/", "final_motifs_binned_methplot_0.6kb_up&downstream.pdf")       
pdf(output_file_name)

#custom_col <- c(green_col(1), red_col1(1))
#custom_col <- c(grey_col(1), red_col(1))
custom_col <- colorRampPalette(c("blue", "yellow", "red"))(10)

summary(read_df) #max_val = 325
#improved_col <- c("#808080",greenred(10) )
par(cex.main=0.80)
test <- heatmap.2(mat_data,
  # dendrogram = "none",
  # Rowv = TRUE, 
  Colv = FALSE,
  #scale="none",
  main = "Methylation distribution from the motif centres", 
  xlab = "Distance(bp) from the centre of motifs",
  ylab = "Transcription Factors",
  #col=greenred(10), 
  col=custom_col, 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.15,
  cexCol = 0.35,
  na.color="grey",
  RowSideColors=read_df$label,
  key.xlab="Methylation fraction"
  #symm=F,symkey=T,symbreaks=T,
  #breaks=c(-1,0.8,1.2,3)

  #breaks=c(c(0,0.9),seq(1,325,(325-1)/9))
  #scale="column", tracecol="#303030"
  )    

 
legend("topright", legend=c("Prom assoc","Enh assoc","CTCF assoc","Multiple assoc","Others"),
fill=c("red","orange","blueviolet","burlywood","green"), border=TRUE, bty="n", y.intersp = 0.9, cex=0.7) 

test
dev.off()



##################

#Same analysis but distal enhancer associated motifs:

##################


### Load the TF annoation file:
anno_df <- fread("~/Dropbox/TFs_Annotation_file.txt", sep="\t")

### 2kb up and downstream heatmap:
input_file <- "~/Dropbox/encode_3/motifs_methplot_total/all_tf_files/all_tfs_motif_bin_methplot_heatmap_data_for_dist_enh_assoc.bed"
read_file <- fread(input_file, sep="\t", header=TRUE)
read_file_merged <- merge(read_file, anno_df, by.x="tf_name", by.y="Target", all.x=TRUE)
#read_file_merged[,c("tf_name", "Category")]
read_file_merged$tf_name_edit <- paste(read_file_merged$tf_name, read_file_merged$annotation)

read_df <- as.data.frame(read_file_merged)
data <- read_df[,2:(ncol(read_df)-4)] ## Subtract last 4 columns
mat_data <- as.matrix(data)
rownames(mat_data) <- read_df$tf_name_edit

colnames(mat_data)
rownames(mat_data)


output_file_name <- paste0(output_dir, "/", "final_motifs_binned_methplot_2kb_up&downstream_distal_enhancer.pdf")       
pdf(output_file_name)

#custom_col <- c(green_col(1), red_col1(1))
#custom_col <- c(grey_col(1), red_col(1))
custom_col <- colorRampPalette(c("blue", "yellow", "red"))(10)

summary(read_df) #max_val = 325
#improved_col <- c("#808080",greenred(10) )
par(cex.main=0.80)
test <- heatmap.2(mat_data,
  # dendrogram = "none",
  # Rowv = TRUE, 
  Colv = FALSE,
  #scale="none",
  main = "Methylation distribution from the motif centres", 
  xlab = "Distance(bp) from the centre of motifs",
  ylab = "Transcription Factors",
  #col=greenred(10), 
  col=custom_col, 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.15,
  cexCol = 0.35,
  na.color="grey",
  RowSideColors=read_df$label,
  key.xlab="Methylation fraction"
  #symm=F,symkey=T,symbreaks=T,
  #breaks=c(-1,0.8,1.2,3)

  #breaks=c(c(0,0.9),seq(1,325,(325-1)/9))
  #scale="column", tracecol="#303030"
  )    
 
legend("topright", legend=c("Prom assoc","Enh assoc","CTCF assoc","Multiple assoc","Others"),
fill=c("red","orange","blueviolet","burlywood","green"), border=TRUE, bty="n", y.intersp = 0.9, cex=0.7) 

test
dev.off()



### 1kb up and downstream heatmap:

### Load the TF annoation file:
anno_df <- fread("~/Dropbox/TFs_Annotation_file.txt", sep="\t")

### 1kb up and downstream heatmap:
input_file <- "~/Dropbox/encode_3/motifs_methplot_total/all_tf_files/all_tfs_motif_bin_methplot_heatmap_data_for_dist_enh_assoc.bed"
read_file <- fread(input_file, sep="\t", header=TRUE)
read_file_merged <- merge(read_file, anno_df, by.x="tf_name", by.y="Target", all.x=TRUE)
#read_file_merged[,c("tf_name", "Category")]
read_file_merged$tf_name_edit <- paste(read_file_merged$tf_name, read_file_merged$annotation)

read_df <- as.data.frame(read_file_merged)
data <- read_df[,2:(ncol(read_df)-4)] ## Subtract last 4 columns
mat_data <- as.matrix(data)

column_ranges <- as.numeric(colnames(mat_data))
column_select <- column_ranges[column_ranges >= -1000 & column_ranges <= 1000] 
mat_data <- mat_data[,colnames(mat_data) %in% column_select]

rownames(mat_data) <- read_df$tf_name_edit

colnames(mat_data)
rownames(mat_data)



output_file_name <- paste0(output_dir, "/", "final_motifs_binned_methplot_1kb_up&downstream_distal_enhancer.pdf")       
pdf(output_file_name)

#custom_col <- c(green_col(1), red_col1(1))
#custom_col <- c(grey_col(1), red_col(1))
custom_col <- colorRampPalette(c("blue", "yellow", "red"))(10)

summary(read_df) #max_val = 325
#improved_col <- c("#808080",greenred(10) )
par(cex.main=0.80)
test <- heatmap.2(mat_data,
  # dendrogram = "none",
  # Rowv = TRUE, 
  Colv = FALSE,
  #scale="none",
  main = "Methylation distribution from the motif centres", 
  xlab = "Distance(bp) from the centre of motifs",
  ylab = "Transcription Factors",
  #col=greenred(10), 
  col=custom_col, 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.15,
  cexCol = 0.35,
  na.color="grey",
  RowSideColors=read_df$label,
  key.xlab="Methylation fraction"
  #symm=F,symkey=T,symbreaks=T,
  #breaks=c(-1,0.8,1.2,3)

  #breaks=c(c(0,0.9),seq(1,325,(325-1)/9))
  #scale="column", tracecol="#303030"
  )    

 
legend("topright", legend=c("Prom assoc","Enh assoc","CTCF assoc","Multiple assoc","Others"),
fill=c("red","orange","blueviolet","burlywood","green"), border=TRUE, bty="n", y.intersp = 0.9, cex=0.7) 

test
dev.off()



### 0.6kb up and downstream heatmap:

### Load the TF annoation file:
anno_df <- fread("~/Dropbox/TFs_Annotation_file.txt", sep="\t")

### 0.6kb up and downstream heatmap:
input_file <- "~/Dropbox/encode_3/motifs_methplot_total/all_tf_files/all_tfs_motif_bin_methplot_heatmap_data_for_dist_enh_assoc.bed"
read_file <- fread(input_file, sep="\t", header=TRUE)
read_file_merged <- merge(read_file, anno_df, by.x="tf_name", by.y="Target", all.x=TRUE)
#read_file_merged[,c("tf_name", "Category")]
read_file_merged$tf_name_edit <- paste(read_file_merged$tf_name, read_file_merged$annotation)

read_df <- as.data.frame(read_file_merged)
data <- read_df[,2:(ncol(read_df)-4)] ## Subtract last 4 columns
mat_data <- as.matrix(data)


column_ranges <- as.numeric(colnames(mat_data))
column_select <- column_ranges[column_ranges >= -600 & column_ranges <= 600] 
mat_data <- mat_data[,colnames(mat_data) %in% column_select]
rownames(mat_data) <- read_df$tf_name_edit

colnames(mat_data)
rownames(mat_data)



output_file_name <- paste0(output_dir, "/", "final_motifs_binned_methplot_0.6kb_up&downstream_distal_enhancer.pdf")       
pdf(output_file_name)

#custom_col <- c(green_col(1), red_col1(1))
#custom_col <- c(grey_col(1), red_col(1))
custom_col <- colorRampPalette(c("blue", "yellow", "red"))(10)

summary(read_df) #max_val = 325
#improved_col <- c("#808080",greenred(10) )
par(cex.main=0.80)
test <- heatmap.2(mat_data,
  # dendrogram = "none",
  # Rowv = TRUE, 
  Colv = FALSE,
  #scale="none",
  main = "Methylation distribution from the motif centres", 
  xlab = "Distance(bp) from the centre of motifs",
  ylab = "Transcription Factors",
  #col=greenred(10), 
  col=custom_col, 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.15,
  cexCol = 0.35,
  na.color="grey",
  RowSideColors=read_df$label,
  key.xlab="Methylation fraction"
  #symm=F,symkey=T,symbreaks=T,
  #breaks=c(-1,0.8,1.2,3)

  #breaks=c(c(0,0.9),seq(1,325,(325-1)/9))
  #scale="column", tracecol="#303030"
  )    
test

 
legend("topright", legend=c("Prom assoc","Enh assoc","CTCF assoc","Multiple assoc","Others"),
fill=c("red","orange","blueviolet","burlywood","green"), border=TRUE, bty="n", y.intersp = 0.9, cex=0.7) 

dev.off()










#######################

## GGplot heatmap, but not necessary for the Encode3 publication:

#######################


### Heatmap:

input_file <- "~/Dropbox/encode_3/motifs_methplot_total/all_tf_files/all_tfs_motif_bin_methplot_heatmap_data.bed"
read_file <- fread(input_file, sep="\t", header=TRUE)
#read_df[is.na(read_df)] <- 0
melted_stacked_df <- melt(read_df, id.vars="tf_name")

output_file_name <- paste0(output_dir, "/", "final_motifs_bin_methplot_2kb.pdf")       
pdf(output_file_name)

p <- ggplot(melted_stacked_df, aes(variable, tf_name)) +
  geom_tile(aes(fill = value)) + 
  #geom_tile(aes(fill = value), colour = "black") + 
  #scale_fill_gradient(low = "white", high = "steelblue") +
  #scale_fill_gradient2(low=muted("blue"), high=muted("red"))+
  #scale_fill_gradient(low="red", high="yellow") +
  scale_fill_gradientn(name="Methylation %", breaks=c(0.1, 0.3, 0.5, 0.7), colours=colorRampPalette(c('blue', 'yellow', "red"))(50)) + 
  labs(x="Distance(bp) from the centre of motifs", y="Transcription Factors", title="Methylation distribution 2kb up & downstream from the motif centres") +
  #geom_text(aes(label = round(value, 1)), size=0.6) +
  theme(axis.text.x = element_text(size = 3.5, angle = -90, hjust = 1)) +
  theme(axis.text.y = element_text(size = 2.2, colour = "black")) +
  theme(legend.text = element_text(size = 6, face="bold", colour = "black")) +
  theme(legend.title = element_text(size = 10, face="bold", hjust= 0.5 )) +
  theme(plot.title = element_text(size = 10, face="bold", hjust= 0.5 )) +
  theme(axis.title=element_text(size=10,face="bold"))
  #coord_flip() 
  # coord_equal()
  #theme_dendro()
  #facet_wrap(~tf_name)
p
#print(p, vp=viewport(angle=180))

dev.off()


### Flag based heatmap:
input_file <- "~/Dropbox/encode_3/motifs_methplot_total/all_tf_files/all_tfs_motif_bin_methplot_heatmap_data.bed"
read_file <- fread(input_file, sep="\t", header=TRUE)
read_file <- read_file[grepl("FLAG", read_file$tf_name)]
read_df <- as.data.frame(read_file)
#read_df[is.na(read_df)] <- 0
melted_stacked_df <- melt(read_df, id.vars="tf_name")

output_file_name <- paste0(output_dir, "/", "flag_based_final_motifs_bin_methplot_2kb.pdf")       
pdf(output_file_name)

p <- ggplot(melted_stacked_df, aes(variable, tf_name)) +
  geom_tile(aes(fill = value)) + 
  #geom_tile(aes(fill = value), colour = "black") + 
  #scale_fill_gradient(low = "white", high = "steelblue") +
  #scale_fill_gradient2(low=muted("blue"), high=muted("red"))+
  #scale_fill_gradient(low="red", high="yellow") +
  scale_fill_gradientn(name="Methylation %", breaks=c(0.1, 0.3, 0.5, 0.7), colours=colorRampPalette(c('blue', 'yellow', "red"))(50)) + 
  labs(x="Distance(bp) from the centre of motifs", y="Transcription Factors", title="Methylation distribution 2kb up & downstream from the motif centres") +
  #geom_text(aes(label = round(value, 1)), size=0.6) +
  theme(axis.text.x = element_text(size = 3.5, angle = -90, hjust = 1)) +
  theme(axis.text.y = element_text(size = 2.2, colour = "black")) +
  theme(legend.text = element_text(size = 6, face="bold", colour = "black")) +
  theme(legend.title = element_text(size = 10, face="bold", hjust= 0.5 )) +
  theme(plot.title = element_text(size = 10, face="bold", hjust= 0.5 )) +
  theme(axis.title=element_text(size=10,face="bold"))
  #coord_flip() 
  # coord_equal()
  #theme_dendro()
  #facet_wrap(~tf_name)
p
#print(p, vp=viewport(angle=180))

dev.off()


### Non-flag_based heatmap:
input_file <- "~/Dropbox/encode_3/motifs_methplot_total/all_tf_files/all_tfs_motif_bin_methplot_heatmap_data.bed"
read_file <- fread(input_file, sep="\t", header=TRUE)
read_file <- read_file[!grepl("FLAG", read_file$tf_name)]
read_df <- as.data.frame(read_file)
#read_df[is.na(read_df)] <- 0
melted_stacked_df <- melt(read_df, id.vars="tf_name")

output_file_name <- paste0(output_dir, "/", "non_flag_based_final_motifs_bin_methplot_2kb.pdf")       
pdf(output_file_name)

p <- ggplot(melted_stacked_df, aes(variable, tf_name)) +
  geom_tile(aes(fill = value)) + 
  #geom_tile(aes(fill = value), colour = "black") + 
  #scale_fill_gradient(low = "white", high = "steelblue") +
  #scale_fill_gradient2(low=muted("blue"), high=muted("red"))+
  #scale_fill_gradient(low="red", high="yellow") +
  scale_fill_gradientn(name="Methylation %", breaks=c(0.1, 0.3, 0.5, 0.7), colours=colorRampPalette(c('blue', 'yellow', "red"))(50)) + 
  labs(x="Distance(bp) from the centre of motifs", y="Transcription Factors", title="Methylation distribution 2kb up & downstream from the motif centres") +
  #geom_text(aes(label = round(value, 1)), size=0.6) +
  theme(axis.text.x = element_text(size = 3.5, angle = -90, hjust = 1)) +
  theme(axis.text.y = element_text(size = 2.2, colour = "black")) +
  theme(legend.text = element_text(size = 6, face="bold", colour = "black")) +
  theme(legend.title = element_text(size = 10, face="bold", hjust= 0.5 )) +
  theme(plot.title = element_text(size = 10, face="bold", hjust= 0.5 )) +
  theme(axis.title=element_text(size=10,face="bold"))
  #coord_flip() 
  # coord_equal()
  #theme_dendro()
  #facet_wrap(~tf_name)
p
#print(p, vp=viewport(angle=180))

dev.off()

