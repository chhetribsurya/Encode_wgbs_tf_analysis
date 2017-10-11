library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
library(gtools)
library(gplots)
library(xlsx)
library("cluster")

args <-  commandArgs(TRUE)
input_file <- args[1]
tf_name <- args[2]
out_dir_name <- args[3]

read_file <- fread("/Users/suryachhetri/Dropbox/encode_3/wgbs_tf_intersect_total/files/final_wgbs_tf_ideas_intersect_heatmap_data_1.bed", sep="\t", header=TRUE)
read_file <- fread("/Users/suryachhetri/Dropbox/encode_3/wgbs_tf_intersect_total/files/final_wgbs_tf_ideas_intersect_heatmap_data_2.bed", sep="\t", header=TRUE)
read_df <- as.data.frame(read_file)

#read_df <- read_df[,select_cols]
read_df <- as.data.frame(read_df)


read_df[is.na(read_df)] <- -0.05

data <- read_df[,2:ncol(read_df)]
mat_data <- as.matrix(data)
rownames(mat_data) <- rnames
rownames(mat_data)
colnames(mat_data) 


output_file_name <- paste0("~/Dropbox", "/", "wgbs_tf_ideas_state_heatmap_all_sites.pdf")        
pdf(output_file_name)

summary(read_df) #max_val = 325
#pair_break =  c(c(-1,0.099), seq(0,1,length=9))
improved_col <- c("#808080",greenred(10) )
par(cex.main=0.6)
heatmap.2(mat_data,
  main = "Methylation Distribution of TFs over the associated IDEAS states", 
  xlab = "IDEAS States",
  ylab = "Transcription Factors",
  col=improved_col,  
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.15,
  cexCol = 0.6,
  #symbreaks = min(mat_data, na.rm=TRUE),
  na.color="grey",
  breaks=c(c(-0.1,-0.02),seq(0,1,(1-0)/9))
  )    

dev.off()

# allmisscols <- apply(read_df,1, function(x)all(is.na(x)));  
#  colswithallmiss <-names(allmisscols[allmisscols>0]);    
#  print("the columns with all values missing");    
#  print(colswithallmiss);

read_file <- fread("~/Dropbox/encode_3/wgbs_tf_intersect_total/files/final_wgbs_tf_chromHMM_intersect_heatmap_data.bed", sep="\t", header=TRUE)
read_df <- as.data.frame(read_file)
read_df[is.na(read_df)] <- 0
rnames <- read_df[,1]
data <- read_df[,2:ncol(read_df)]
mat_data <- as.matrix(data)
rownames(mat_data) <- rnames
rownames(mat_data)
colnames(mat_data)

output_file_name <- paste0("~/Dropbox", "/", "wgbs_tf_ideas_state_heatmap_total.pdf")        
pdf(output_file_name)

par(cex.main=0.6)
heatmap.2(mat_data,
  main = "Methylation Distribution of TFs over the associated IDEAS states", 
  xlab = "IDEAS States",
  ylab = "Transcription Factors",
  col=greenred(10),  
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.15,
  cexCol = 0.6,
  na.color="grey"
  )    

dev.off()


#############################

#### For only the regulatory_regions of the IDEAS states:

#read_file <- fread("/Users/suryachhetri/Dropbox/encode_3/wgbs_tf_intersect_total/files/final_wgbs_tf_ideas_intersect_heatmap_data_1.bed", sep="\t", header=TRUE)
#read_file <- fread("/Users/suryachhetri/Dropbox/encode_3/wgbs_tf_intersect_total/files/final_wgbs_tf_ideas_intersect_heatmap_data_2.bed", sep="\t", header=TRUE)
read_file <- fread("/Users/suryachhetri/Dropbox/encode_3/wgbs_tf_intersect_total/files/final_heatmap_files/final_wgbs_tf_ideas_intersect_heatmap_data_20.bed", sep="\t", header=TRUE)
read_file <- fread("/Users/suryachhetri/Dropbox/encode_3/wgbs_tf_intersect_total/files/final_heatmap_files/final_wgbs_tf_ideas_intersect_heatmap_data_30.bed", sep="\t", header=TRUE)


ideas_mnemonics_df = fread("~/Dropbox/ideas_table.txt", sep="\t")
regulatory_regions <- c("tf_name", ideas_mnemonics_df$Mnemonics[0:20])

select_cols <- colnames(read_file) %in% regulatory_regions
read_select_file <- read_file[,select_cols, with=FALSE]

read_df <- as.data.frame(read_select_file)

read_df$ElonW <- NULL
read_df$Pol2 <- NULL
read_df$EnhWF3 <- NULL

read_df[is.na(read_df)] <- -0.05

rnames <- read_df[,1]
data <- read_df[,2:ncol(read_df)]
#data <- apply(data, 2, as.numeric)
mat_data <- as.matrix(data)
rownames(mat_data) <- rnames
rownames(mat_data)
colnames(mat_data) 


output_file_name <- paste0("~/Dropbox", "/", "wgbs_tf_ideas_state_heatmap_20_cpgs_cutoff_confirm_final_final.pdf")        
output_file_name <- paste0("~/Dropbox", "/", "wgbs_tf_ideas_state_heatmap_30_cpgs_cutoff_confirm_final.pdf")        
pdf(output_file_name)

summary(read_df) #max_val = 325
#pair_break =  c(c(-1,0.099), seq(0,1,length=9))
improved_col <- c("#808080",greenred(10) )
par(cex.main=0.6)
heatmap.2(mat_data,
  main = "Methylation Distribution of TFs over the associated IDEAS states", 
  xlab = "IDEAS States",
  ylab = "Transcription Factors",
  col=improved_col,  
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.15,
  cexCol = 0.6,
  #symbreaks = min(mat_data, na.rm=TRUE),
  na.color="grey",
  breaks=c(c(-0.1,-0.02),seq(0,1,(1-0)/9))
  )    

dev.off()




###############################################
# For chromHMM models:


#### For only the regulatory_regions of the IDEAS states:

read_file <- fread("/Users/suryachhetri/Dropbox/encode_3/wgbs_tf_intersect_total/files/final_heatmap_files/final_wgbs_tf_chromHMM_intersect_heatmap_data.bed", sep="\t", header=TRUE)
read_file <- fread("/Users/suryachhetri/Dropbox/encode_3/wgbs_tf_intersect_total/files/final_heatmap_files/final_wgbs_tf_chromHMM_intersect_heatmap_data_20.bed", sep="\t", header=TRUE)
read_file <- fread("/Users/suryachhetri/Dropbox/encode_3/wgbs_tf_intersect_total/files/final_heatmap_files/final_wgbs_tf_chromHMM_intersect_heatmap_data_30.bed", sep="\t", header=TRUE)


#ideas_mnemonics_df = fread("~/Dropbox/ideas_table.txt", sep="\t")
#regulatory_regions <- c("tf_name", ideas_mnemonics_df$Mnemonics[0:20])

#select_cols <- colnames(read_file) %in% regulatory_regions
#read_select_file <- read_file[,select_cols, with=FALSE]

read_df <- as.data.frame(read_file)

#read_df$ElonW <- NULL
#read_df$Pol2 <- NULL
#read_df$EnhWF3 <- NULL

read_df[is.na(read_df)] <- -0.05

rnames <- read_df[,1]
data <- read_df[,2:ncol(read_df)]
#data <- apply(data, 2, as.numeric)
mat_data <- as.matrix(data)
rownames(mat_data) <- rnames
rownames(mat_data)
colnames(mat_data) 


output_file_name <- paste0("~/Dropbox", "/", "wgbs_tf_chromHMM_state_heatmap_20_cpgs_cutoff.pdf")        
output_file_name <- paste0("~/Dropbox", "/", "wgbs_tf_chromHMM_state_heatmap_all_cpgs.pdf")        
pdf(output_file_name)

summary(read_df) #max_val = 325
#pair_break =  c(c(-1,0.099), seq(0,1,length=9))
improved_col <- c("#808080",greenred(10) )
par(cex.main=0.6)
heatmap.2(mat_data,
  main = "Methylation Distribution of TFs over the associated IDEAS states", 
  xlab = "IDEAS States",
  ylab = "Transcription Factors",
  col=improved_col,  
  #col=greenred(10),  
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.15,
  cexCol = 0.6,
  #symbreaks = min(mat_data, na.rm=TRUE),
  breaks=c(c(-0.1,-0.02),seq(0,1,(1-0)/9)),
  na.color="grey"
  )    

dev.off()













