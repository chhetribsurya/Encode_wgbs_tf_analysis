library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(gtools)
library(gplots)

# weights <- ifelse(pcaOutput2$PC1 < -5 & abs(pcaOutput2$PC2) > 10, 2, 1)

### Heatmap:

### With genebody:
input_file <- "~/Dropbox/encode_3/motifs_methplot_total/all_tf_files/final_wgbs_tf_motif_ideas_intersect_heatmap_data_with_cpg_count20_and_genebody.bed"
read_file <- fread(input_file, sep="\t", header=TRUE)

read_file <- read_file[grepl("FLAG", read_file$tf_name)]
read_file <- read_file[!grepl("FLAG", read_file$tf_name)]

### Also, can be done using stringr:
#read_file[str_detect(read_file$tf_name, "FLAG")] #str_detect(chars, value)
# > library(stringr)
# > chars <- "test"
# > value <- "es"
# > str_detect(chars, value)
# [1] TRUE
# > value <- c("es", "l", "est", "a", "test")
# > str_detect(chars, value)
# [1]  TRUE FALSE  TRUE FALSE  TRUE

read_df <- as.data.frame(read_file)
read_df[is.na(read_df)] <- -0.05
#read_df[is.na(read_df)] <- 0
rnames <- read_df[,1]
data <- read_df[,2:ncol(read_df)]
mat_data <- as.matrix(data)
rownames(mat_data) <- rnames
rownames(mat_data)
colnames(mat_data)


output_file_name <- paste0("~/Dropbox", "/", "final_motif_ideas_intersect_heatmap_data_with_cpg_count20_w_genebody.pdf")       
pdf(output_file_name)

#green_col <-  colorRampPalette("green")
#red_col <-  colorRampPalette("red")

# grey_col <-  colorRampPalette("grey")
# red_col <- colorRampPalette(c("indianred1"))
# wes_palette("Zissou")
# my_palette <- wes_palette("Zissou", 100, type = "continuous")
# custom_col <- c(green_col(1), red_col(1))
# custom_col <- c(grey_col(1), red_col(1))

summary(read_df) #max_val = 325
improved_col <- c("#808080",greenred(50) )
par(cex.main=0.55)
test <- heatmap.2(mat_data,
  #dendrogram = "none",
  #Rowv = "none", 
  #Colv = "none",
  #scale="none",
  main = "Avg. methylation distribution across IDEAS states", 
  xlab = "IDEAS states",
  ylab = "Transcription Factors",
  col=improved_col,
  #col=greenred(200), 
  #col=custom_col, 
  #col=my_palette, 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.15,
  cexCol = 0.45,
  na.color="grey",
  #symm=F,symkey=T,symbreaks=T,
  breaks=c(c(-0.1,-0.02),seq(0,1,(1-0)/49))
  #breaks=c(c(0,0.9),seq(1,325,(325-1)/9))
  #scale="column", tracecol="#303030"
  )    
test
dev.off()


### For the flag based data:
input_file <- "~/Dropbox/encode_3/motifs_methplot_total/all_tf_files/final_wgbs_tf_motif_ideas_intersect_heatmap_data_with_cpg_count20_and_genebody.bed"
read_file <- fread(input_file, sep="\t", header=TRUE)

read_file <- read_file[grepl("FLAG", read_file$tf_name)]
read_df <- as.data.frame(read_file)
read_df[is.na(read_df)] <- -0.05
#read_df[is.na(read_df)] <- 0
rnames <- read_df[,1]
data <- read_df[,2:ncol(read_df)]
mat_data <- as.matrix(data)
rownames(mat_data) <- rnames
rownames(mat_data)
colnames(mat_data)


output_file_name <- paste0("~/Dropbox", "/", "flag_based_final_motif_ideas_intersect_heatmap_data_w_genebody.pdf")       
pdf(output_file_name)

summary(read_df) #max_val = 325
improved_col <- c("#808080",greenred(50) )
par(cex.main=0.55)
test <- heatmap.2(mat_data,
  #dendrogram = "none",
  #Rowv = "none", 
  #Colv = "none",
  #scale="none",
  main = "Avg. methylation distribution across IDEAS states", 
  xlab = "IDEAS states",
  ylab = "Transcription Factors",
  col=improved_col,
  #col=greenred(200), 
  #col=custom_col, 
  #col=my_palette, 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.15,
  cexCol = 0.45,
  na.color="grey",
  #symm=F,symkey=T,symbreaks=T,
  breaks=c(c(-0.1,-0.02),seq(0,1,(1-0)/49))
  #breaks=c(c(0,0.9),seq(1,325,(325-1)/9))
  #scale="column", tracecol="#303030"
  )    
test
dev.off()



### For the non-flag based data:
input_file <- "~/Dropbox/encode_3/motifs_methplot_total/all_tf_files/final_wgbs_tf_motif_ideas_intersect_heatmap_data_with_cpg_count20_and_genebody.bed"
read_file <- fread(input_file, sep="\t", header=TRUE)

read_file <- read_file[!grepl("FLAG", read_file$tf_name)]
read_df <- as.data.frame(read_file)
read_df[is.na(read_df)] <- -0.05
#read_df[is.na(read_df)] <- 0
rnames <- read_df[,1]
data <- read_df[,2:ncol(read_df)]
mat_data <- as.matrix(data)
rownames(mat_data) <- rnames
rownames(mat_data)
colnames(mat_data)


output_file_name <- paste0("~/Dropbox", "/", "non_flag_based_final_motif_ideas_intersect_heatmap_data_w_genebody.pdf")       
pdf(output_file_name)

summary(read_df) #max_val = 325
improved_col <- c("#808080",greenred(50) )
par(cex.main=0.55)
test <- heatmap.2(mat_data,
  #dendrogram = "none",
  #Rowv = "none", 
  #Colv = "none",
  #scale="none",
  main = "Avg. methylation distribution across IDEAS states", 
  xlab = "IDEAS states",
  ylab = "Transcription Factors",
  col=improved_col,
  #col=greenred(200), 
  #col=custom_col, 
  #col=my_palette, 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.15,
  cexCol = 0.45,
  na.color="grey",
  #symm=F,symkey=T,symbreaks=T,
  breaks=c(c(-0.1,-0.02),seq(0,1,(1-0)/49))
  #breaks=c(c(0,0.9),seq(1,325,(325-1)/9))
  #scale="column", tracecol="#303030"
  )    
test
dev.off()







### Without genebody:
input_file <- "~/Dropbox/encode_3/motifs_methplot_total/all_tf_files/final_wgbs_tf_motif_ideas_intersect_heatmap_data_with_cpg_count20.bed"
read_file <- fread(input_file, sep="\t", header=TRUE)
#read_file_flag <- read_file[grepl("FLAG", read_file$tf_name)]
#read_file_nonflag <- read_file[!grepl("FLAG", read_file$tf_name)]

read_df <- as.data.frame(read_file)
read_df[is.na(read_df)] <- -0.05
#read_df[is.na(read_df)] <- 0
rnames <- read_df[,1]
data <- read_df[,2:ncol(read_df)]
mat_data <- as.matrix(data)
rownames(mat_data) <- rnames
rownames(mat_data)
colnames(mat_data)


output_file_name <- paste0("~/Dropbox", "/", "final_motif_ideas_intersect_heatmap_data_with_cpg_count20.pdf")       
pdf(output_file_name)

green_col <-  colorRampPalette("green")
red_col <-  colorRampPalette("red")

# grey_col <-  colorRampPalette("grey")
# red_col <- colorRampPalette(c("indianred1"))
# wes_palette("Zissou")
# my_palette <- wes_palette("Zissou", 100, type = "continuous")
# custom_col <- c(green_col(1), red_col(1))
#custom_col <- c(grey_col(1), red_col(1))

summary(read_df) #max_val = 325
improved_col <- c("#808080",greenred(50) )
par(cex.main=0.55)
test <- heatmap.2(mat_data,
  #dendrogram = "none",
  #Rowv = "none", 
  #Colv = "none",
  #scale="none",
  main = "Avg. methylation distribution across IDEAS states", 
  xlab = "IDEAS states",
  ylab = "Transcription Factors",
  col=improved_col,
  #col=greenred(200), 
  #col=custom_col, 
  #col=my_palette, 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.15,
  cexCol = 0.45,
  na.color="grey",
  #symm=F,symkey=T,symbreaks=T,
  breaks=c(c(-0.1,-0.02),seq(0,1,(1-0)/49))
  #breaks=c(c(0,0.9),seq(1,325,(325-1)/9))
  #scale="column", tracecol="#303030"
  )    
test
dev.off()




#### If needed, here's the ggplot heatmap, but not needed for the above case:
# Plot the strip heatmap of RPKM for genes in 
# the same order as the main/central heatmap:
order_of_heatmap <- mat_data[rev(test$rowInd), test$colInd]
gene_names <- colnames(order_of_heatmap)
tf_names <- rownames(order_of_heatmap)

heatmap_ordered_df <- tf_count_sp_gene_df[match(gene_names, tf_count_sp_gene_df$gene_name),]
melted_stacked_df <- melt(heatmap_ordered_df, id.vars="gene_name", measure.vars =c("HepG2_RNAseq", "Liver_RNAseq") )
melted_stacked_df$re_gene_name <- factor(melted_stacked_df$gene_name, levels=melted_stacked_df$gene_name)
casted_unstacked_df <- cast(melted_stacked_df, variable~re_gene_name) # but, not used for these analysis

output_file_name <- paste0("~/Dropbox", "/", "Liver_hepg2_promoters_heatmap_TPM_chart_1.pdf")       
pdf(output_file_name)

p <- ggplot(melted_stacked_df, aes(re_gene_name, variable)) +
  geom_tile(aes(fill = value), colour = "black") + 
  scale_fill_gradient(low = "white", high = "steelblue") +
  #scale_fill_gradient2(low=muted("blue"), high=muted("red"))+
  #scale_fill_gradient(low="red", high="yellow") +
  #scale_fill_gradient(low = "white", high = "red") + 
  labs(x="Genes", y="", title="") +
  geom_text(aes(label = round(value, 1)), size=0.6) +
  theme(axis.text.x = element_text(size = 2.5, angle = 90, hjust = 1)) +
  theme(axis.text.y = element_text(size = 4, colour = "black")) +
  #coord_flip() 
  coord_fixed(ratio = 1.5)
  #theme_dendro()
  #facet_wrap(~variable)

p

dev.off()





# library(ggplot2)
# #Reading in the data
# chicagoMVT <- read.csv('motor_vehicle_theft.csv', stringsAsFactors = FALSE)
# #Converting the date to a recognizable format
# chicagoMVT$Date <- strptime(chicagoMVT$Date, format = '%m/%d/%Y %I:%M:%S %p')
# #Getting the day and hour of each crime
# chicagoMVT$Day <- weekdays(chicagoMVT$Date)
# chicagoMVT$Hour <- chicagoMVT$Date$hour
# #Sorting the weekdays
# dailyCrimes <- as.data.frame(table(chicagoMVT$Day, chicagoMVT$Hour))
# names(dailyCrimes) <- c('Day', 'Hour', 'Freq')
# dailyCrimes$Hour <- as.numeric(as.character(dailyCrimes$Hour))
# dailyCrimes$Day <- factor(dailyCrimes$Day, ordered = TRUE, 
#                          levels = c('Sunday', 'Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday'))
# #Plotting the number of crimes each day (line graph)
# #ggplot(dailyCrimes, aes(x = Hour, y = Freq)) + geom_line(aes(group = Day, color = Day)) + xlab('Hour') + ylab('Number of thefts') + ggtitle('Daily number of Motor Vehicle Thefts')

# ggplot(dailyCrimes, aes(x = Hour, y = Day)) + geom_tile(aes(fill = Freq)) + scale_fill_gradient(name = 'Total Motor Vehicle Thefts', low = 'white', high = 'red') + theme(axis.title.y = element_blank())


# library(ggplot2)
# library(reshape2)      # for melt(...)
# library(grid)          # for unit(...)

# set.seed(1)            # for reproducible example
# df <- data.frame(matrix(rnorm(100*10), nr=10))
# df.melt <- melt(cbind(x=1:nrow(df),df),id="x")
# ggplot(df.melt,aes(x=factor(x),y=variable,fill=value)) +
#   geom_tile() +
#   labs(x="",y="")+
#   scale_x_discrete(expand=c(0,0))+
#   scale_fill_gradientn(name="", limits=c(-3,3),
#                        colours=colorRampPalette(c('blue', 'yellow'))(12))+
#   theme(legend.position="top", 
#         legend.key.width=unit(.1,"npc"),legend.key.height=unit(.05,"npc"),
#         axis.text=element_blank(),axis.ticks=element_blank())


