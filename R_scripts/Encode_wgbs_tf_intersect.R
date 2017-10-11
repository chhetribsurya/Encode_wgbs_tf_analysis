library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
library(gtools)
library(gplots)
library(xlsx)

args <-  commandArgs(TRUE)
input_file <- args[1]
tf_name <- args[2]
out_dir_name <- args[3]

#avg_df <- fread("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/wgbs_tf_intersect/final_wgbs_tf_intersect_average_meth.bed", sep="\t", header=TRUE)
xls_df <- read.xlsx2("~/Dropbox/Encode_full_hepg2_datasets_DBF_CR.xls", sheetIndex=4)
avg_df <- fread("~/Dropbox/encode_3/wgbs_tf_intersect_total/files/final_wgbs_tf_intersect_average_meth.bed", sep="\t", header=TRUE)
meth_tf_df <- fread("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/wgbs_tf_intersect/final_wgbs_tf_intersect.bed", sep="\t", header=TRUE)

names(avg_df) <- c("TF_name", "meth_percent")
names(meth_tf_df) <- c("TF_name", "index", "peak_chrom", "peak_start", "peak_end", "meth_percent")

target_vector <- avg_df$TF_name
xls_df_ordered <- xls_df[match(target_vector, xls_df$Target),]

# check if the order b/w 2 columns are equal. 
if(all(avg_df$TF_name == xls_df_ordered$Target))
{
   print("Column A(tf_name) and B(target) are identical")
}

xls_df_ordered <- xls_df_ordered %>% select(Target, Lab, Category)
avg_df$tf_category <-  xls_df_ordered$Category

cr_df <- avg_df %>% as.data.frame %>% filter(tf_category == "CR/CF")
dbf_df <- avg_df %>% as.data.frame %>% filter(tf_category == "DBF")
all_tf_df <- avg_df %>% as.data.frame


# or, do the following matching by ordering:
df1 <- xls_df %>%
		select(Target, Lab, Category) %>%
		arrange(Target) %>% as.data.table()

df2 <- avg_df %>% 
		arrange(TF_name) %>% as.data.table()


output_file_name <- paste0("~/Dropbox", "/", "all_tfs_methylation_distribution.pdf")				
pdf(output_file_name)

barplot <-ggplot(all_tf_df, aes(x=tf_category, y=meth_percent, fill=tf_category)) + 
	#geom_point() + 
    geom_boxplot(width=0.15, fill = "grey80", outlier.colour = "red", color="blue") + # outlier.shape = NA to remove the outliers 
	geom_dotplot(data=subset(all_tf_df, meth_percent < 0.15), binaxis='y', stackdir='centerwhole', dotsize=0.2, fill="black") +
	stat_boxplot( aes(tf_category, meth_percent), 
    geom='errorbar', linetype=2, width=0.15, color="blue") +  
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    xlab("Transcription Factor Category") + ylab("Average methylation percent") + theme_bw()+ 
    ggtitle("Average Methylation Distribution of 208 TFs") + 
    scale_y_continuous(limits=c(0,1)) +
    theme(
    axis.text.x = element_text(size=10, face="bold"),
	plot.title=element_text(size=14, face="bold", hjust = 0.6)
	)

barplot <- barplot + geom_text(data = subset(all_tf_df, meth_percent > 0.15), aes(x = tf_category, y =meth_percent, label = TF_name), 
	hjust= -0.15, size = 2, check_overlap = TRUE)

print(barplot)

dev.off()


# check if the order b/w 2 columns are equal. 
# if(all(df1$Target == df2$TF_name))
# {
#    print("Column A and B are identical")
# }

# xls_ordered <- xls_df[mixedorder(xls_df$Target),] %>% as.data.table()
# avg_df_ordered <- avg_df[mixedorder(avg_df$TF_name),]


# barplot <-ggplot(all_tf_df, aes(x=" ", y=meth_percent, fill=tf_category)) + 
# 	geom_point() + 
#     geom_boxplot(width=0.15, fill = "grey80", color="blue") + 
# 	geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2) +
# 	stat_boxplot( aes(" ", meth_percent), 
#     geom='errorbar', linetype=2, width=0.15) +  
#     geom_jitter(shape=16, position=position_jitter(0.004)) + 
#     xlab("Transcription Factors") + ylab("Average methylation percent") + theme_bw()+ 
#     ggtitle("Average Methylation Distribution of 208 TFs") + 
#     scale_y_continuous(limits=c(0,1)) +
#     theme(
# 	plot.title=element_text(size=14, face="bold", hjust = 0.6)
# 	)

# # dot plot with mean points
# p + stat_summary(fun.y=mean, geom="point", shape=18,
#                  size=3, color="red")
# # dot plot with median points
# p + stat_summary(fun.y=median, geom="point", shape=18,
#                  size=3, color="red")


# filter(msleep, sleep_total >= 16)
# filter(msleep, sleep_total >= 16, bodywt >= 1) # multiple conditions
# head(select(msleep, -name)) # drop column "name"

### Add new column to the dataframe:
# msleep %>% 
#     mutate(rem_proportion = sleep_rem / sleep_total) %>%
#     head

### add many new columns using comma(,)
# msleep %>% 
#     mutate(rem_proportion = sleep_rem / sleep_total, 
#            bodywt_grams = bodywt * 1000) %>%
#     head

# msleep %>% 
#     group_by(order) %>%
#     summarise(avg_sleep = mean(sleep_total), 
#               min_sleep = min(sleep_total), 
#               max_sleep = max(sleep_total),
#               total = n())

# msleep %>% 
#     select(name, order, sleep_total) %>%
#     arrange(order, sleep_total) %>% 
#     filter(sleep_total >= 16)


#####################
#####################


avg_df$meth_proportion = ""
less_than_5_indices <- (avg_df$meth_percent < 0.05)
less_than_10_indices <- (avg_df$meth_percent >= 0.05 & avg_df$meth_percent <0.10 )
less_than_15_indices <- (avg_df$meth_percent >= 0.10 & avg_df$meth_percent < 0.15 )
more_than_15_indices <- (avg_df$meth_percent >= 0.15)

avg_df[less_than_5_indices, ]$meth_proportion <- "<0.05"
avg_df[less_than_10_indices, ]$meth_proportion <- "0.05-0.10"
avg_df[less_than_15_indices, ]$meth_proportion <- "0.10-0.15"
avg_df[more_than_15_indices, ]$meth_proportion <- ">=0.15"

factor_order = c("<0.05", "0.05-0.10", "0.10-0.15", ">=0.15")


### Maintaining the order of factor:
factor(avg_df$meth_proportion)
avg_df$meth_proportion <- factor(avg_df$meth_proportion, levels=factor_order)

output_file_name <- paste0("~/Dropbox", "/", "categorised_tf_barplot.pdf")				
pdf(output_file_name)

barplot <-ggplot(avg_df, aes(x=meth_proportion, y=meth_percent, color=meth_proportion)) + 
    geom_boxplot(width= 0.9, fill = "grey80",outlier.colour = "red", outlier.shape = 1) + 
    geom_point(colour = "blue", position = "jitter", alpha=0.5) +
    stat_boxplot( aes(meth_proportion, meth_percent), 
    geom='errorbar', linetype=1, width=0.5)+  #whiskers
	#geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25, color="black") + 
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    #stat_summary(fun.y=mean, geom="point", size=2) + 
    #stat_summary(fun.data = mean_se, geom = "errorbar") +

    xlab("Transcription Factors Category") + ylab("Average methylation percent") + 
    theme_bw() + labs(color='% Methylation Level') +
    ggtitle("Average Methylation Distribution of 102 TFs") + 
    theme(
	plot.title=element_text(size=14, face="bold", hjust = 0.6),
	legend.title=element_text(face="bold")
	)

print(barplot)

dev.off()


### Split the set2 into 3 groups:
##############
set2_df_split <- avg_df[less_than_10_indices, ]
set2_df_ordered <- set2_df_split[order(set2_df_split$meth_percent),]

### Match the TF with group(avg_df) to extract the TF from the main file that is meth_tf_df:
target_vector_1 =  avg_df[(avg_df$meth_proportion == "<0.05"),]$TF_name
target_vector_2 =  avg_df[(avg_df$meth_proportion == "0.05-0.10"),]$TF_name
target_vector_3 =  avg_df[(avg_df$meth_proportion == "0.10-0.15"),]$TF_name
target_vector_2_1 <- set2_df_ordered[1:17,]$TF_name
target_vector_2_2 <- set2_df_ordered[18:34,]$TF_name
target_vector_2_3 <- set2_df_ordered[35:50,]$TF_name
target_vector_4 =  avg_df[(avg_df$meth_proportion == ">=0.15"),]$TF_name

# meth_df_1 <- meth_tf_df[match(target_vector_1, meth_tf_df$TF_name),]
# meth_df_2 <- meth_tf_df[match(target_vector_2, meth_tf_df$TF_name),]
# meth_df_3 <- meth_tf_df[match(target_vector_3, meth_tf_df$TF_name),]
# meth_df_4 <- meth_tf_df[match(target_vector_4, meth_tf_df$TF_name),]

meth_df_list_1 <- lapply(target_vector_1, function(x){ 
	match_indices <- which(x == meth_tf_df$TF_name); 
	meth_tf_df[match_indices,]}
	)

meth_df_list_2 <- lapply(target_vector_2, function(x){ 
	match_indices <- which(x == meth_tf_df$TF_name); 
	meth_tf_df[match_indices,]}
	)


meth_df_list_3 <- lapply(target_vector_3, function(x){ 
	match_indices <- which(x == meth_tf_df$TF_name); 
	meth_tf_df[match_indices,]}
	)


meth_df_list_2_1 <- lapply(target_vector_2_1, function(x){ 
	match_indices <- which(x == meth_tf_df$TF_name); 
	meth_tf_df[match_indices,]}
	)


meth_df_list_2_2 <- lapply(target_vector_2_2, function(x){ 
	match_indices <- which(x == meth_tf_df$TF_name); 
	meth_tf_df[match_indices,]}
	)


meth_df_list_2_3 <- lapply(target_vector_2_3, function(x){ 
	match_indices <- which(x == meth_tf_df$TF_name); 
	meth_tf_df[match_indices,]}
	)


meth_df_list_4 <- lapply(target_vector_4, function(x){ 
	match_indices <- which(x == meth_tf_df$TF_name); 
	meth_tf_df[match_indices,]}
	)


### Combining the list of lists with rbind:
meth_df_1 <- do.call(rbind, meth_df_list_1)
meth_df_2 <- do.call(rbind, meth_df_list_2)
meth_df_3 <- do.call(rbind, meth_df_list_3)
meth_df_2_1 <- do.call(rbind, meth_df_list_2_1)
meth_df_2_2 <- do.call(rbind, meth_df_list_2_2)
meth_df_2_3 <- do.call(rbind, meth_df_list_2_3)
meth_df_4 <- do.call(rbind, meth_df_list_4)

### Sanity check if the target vector matches the new list:
unique(meth_df_1$TF_name) == target_vector_1
unique(meth_df_2$TF_name) == target_vector_2
unique(meth_df_3$TF_name) == target_vector_3
unique(meth_df_2_1$TF_name) == target_vector_2_1
unique(meth_df_2_2$TF_name) == target_vector_2_2
unique(meth_df_2_3$TF_name) == target_vector_2_3
unique(meth_df_4$TF_name) == target_vector_4


### Check the levels:
factor(meth_df_2_1$TF_name)
factor(meth_df_2_2$TF_name)
factor(meth_df_2_3$TF_name)

### Maintaining the order of factor/Correct the factors, or maintain the order
meth_df_2_1$TF_name <- factor(meth_df_2_1$TF_name, levels= target_vector_2_1)
meth_df_2_2$TF_name <- factor(meth_df_2_2$TF_name, levels= target_vector_2_2)
meth_df_2_3$TF_name <- factor(meth_df_2_3$TF_name, levels= target_vector_2_3)

### Crosscheck the factor, if its reordered:
factor(meth_df_2_1$TF_name)
factor(meth_df_2_2$TF_name)
factor(meth_df_2_3$TF_name)
################


output_file_name <- paste0("~/Dropbox", "/", "category_1_tf_barplot.pdf")				
pdf(output_file_name)

barplot <- ggplot(meth_df_1, aes(x=TF_name, y=meth_percent, color=TF_name)) + 
    geom_boxplot(outlier.shape = NA) + geom_jitter(size=0.01, color="black", alpha=0.01) +
    #geom_point(colour = "blue", position = "jitter", alpha=0.5) +
    stat_boxplot( aes(TF_name, meth_percent), 
    geom='errorbar', linetype=1, width=0.5)+  #whiskers
	#geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25, color="black") + 
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    #stat_summary(fun.y=mean, geom="point", size=2) + 
    #stat_summary(fun.data = mean_se, geom = "errorbar") +

    xlab("Transcription Factors") + ylab("Average methylation percent") + 
    theme_bw() + labs(color='Transcription Factors(<0.5)') +
    ggtitle("Avg. Methylation across TF binding sites") + 
    theme(
	plot.title=element_text(size=14, face="bold", hjust = 0.6),
	legend.title=element_text(face="bold")
	) + coord_flip() + scale_y_continuous(limits = quantile(meth_df_1$meth_percent, c(0.1, 0.8)))

print(barplot)

dev.off()


output_file_name <- paste0("~/Dropbox", "/", "category_1_tf_violinplot.pdf")				
pdf(output_file_name)

violinplot <- ggplot(meth_df_1, aes(x=TF_name, y=meth_percent, color=TF_name)) +
    geom_violin() + 
    geom_point(size = 0.05, position = "jitter", alpha=0.1) +
    #stat_boxplot( aes(TF_name, meth_percent), 
    #geom='errorbar', linetype=1, width=0.5)+  #whiskers
	#geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25, color="black") + 
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    #stat_summary(fun.y=mean, geom="point", size=2) + 
    #stat_summary(fun.data = mean_se, geom = "errorbar") +

    xlab("Transcription Factors") + ylab("Average methylation percent") + 
    theme_bw() + labs(color='Transcription Factors(<0.5)') +
    ggtitle("Avg. Methylation across TF binding sites") + 
    theme(
	plot.title=element_text(size=14, face="bold", hjust = 0.6),
	legend.title=element_text(face="bold")
	) +coord_flip() + scale_y_continuous(limits = quantile(meth_df_1$meth_percent, c(0.1, 1.0)))


print(violinplot)


dev.off()


##########
##########

output_file_name <- paste0("~/Dropbox", "/", "category_2_tf_barplot.pdf")				
pdf(output_file_name)

barplot <- ggplot(meth_df_2, aes(x=TF_name, y=meth_percent, color=TF_name)) + 
    geom_boxplot(outlier.shape = NA) + geom_jitter(size=0.01, color="black", alpha=0.02) +
    #geom_point(colour = "blue", position = "jitter", alpha=0.5) +
    stat_boxplot( aes(TF_name, meth_percent), 
    geom='errorbar', linetype=1, width=0.5)+  #whiskers
	#geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25, color="black") + 
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    #stat_summary(fun.y=mean, geom="point", size=2) + 
    #stat_summary(fun.data = mean_se, geom = "errorbar") +

    xlab("Transcription Factors ") + ylab("Average methylation percent") + 
    theme_bw() + labs(color='Transcription Factors(0.05-0.10)') +
    ggtitle("Avg. Methylation across TF binding sites") + 
    theme(
	plot.title=element_text(size=14, face="bold", hjust = 0.6),
	legend.title=element_text(face="bold")
	) + coord_flip() + scale_y_continuous(limits = c(0.00, 0.4))

print(barplot)

dev.off()


output_file_name <- paste0("~/Dropbox", "/", "category_2_tf_violinplot.pdf")				
pdf(output_file_name)

violinplot <- ggplot(meth_df_2, aes(x=TF_name, y=meth_percent, color=TF_name)) +
    geom_violin() + 
    geom_point(size = 0.05, position = "jitter", alpha=0.1) +
    #stat_boxplot( aes(TF_name, meth_percent), 
    #geom='errorbar', linetype=1, width=0.5)+  #whiskers
	#geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25, color="black") + 
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    #stat_summary(fun.y=mean, geom="point", size=2) + 
    #stat_summary(fun.data = mean_se, geom = "errorbar") +

    xlab("Transcription Factors") + ylab("Average methylation percent") + 
    theme_bw() + labs(color='Transcription Factors(0.05-0.10)') +
    ggtitle("Avg. Methylation across TF binding sites") + 
    theme(
	plot.title=element_text(size=14, face="bold", hjust = 0.6),
	legend.title=element_text(face="bold")
	) +coord_flip() + scale_y_continuous(limits = quantile(meth_df_1$meth_percent, c(0.1, 1.0)))


print(violinplot)


dev.off()


output_file_name <- paste0("~/Dropbox", "/", "category_3_tf_barplot.pdf")				
pdf(output_file_name)

barplot <- ggplot(meth_df_3, aes(x=TF_name, y=meth_percent, color=TF_name)) + 
    geom_boxplot(outlier.shape = NA) + geom_jitter(size=0.01, color="black", alpha=0.02) +
    #geom_point(colour = "blue", position = "jitter", alpha=0.5) +
    stat_boxplot( aes(TF_name, meth_percent), 
    geom='errorbar', linetype=1, width=0.5)+  #whiskers
	#geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25, color="black") + 
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    #stat_summary(fun.y=mean, geom="point", size=2) + 
    #stat_summary(fun.data = mean_se, geom = "errorbar") +

    xlab("Transcription Factors") + ylab("Average methylation percent") + 
    theme_bw() + labs(color='Transcription Factors(0.10-0.15)') +
    ggtitle("Avg. Methylation across TF binding sites") + 
    theme(
	plot.title=element_text(size=14, face="bold", hjust = 0.6),
	legend.title=element_text(face="bold")
	) + coord_flip() + scale_y_continuous(limits = c(0.0, 0.5))

print(barplot)

dev.off()


output_file_name <- paste0("~/Dropbox", "/", "category_3_tf_violinplot.pdf")				
pdf(output_file_name)

violinplot <- ggplot(meth_df_3, aes(x=TF_name, y=meth_percent, color=TF_name)) +
    geom_violin() + 
    geom_point(size = 0.05, position = "jitter", alpha=0.01) +
    #stat_boxplot( aes(TF_name, meth_percent), 
    #geom='errorbar', linetype=1, width=0.5)+  #whiskers
	#geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25, color="black") + 
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    #stat_summary(fun.y=mean, geom="point", size=2) + 
    #stat_summary(fun.data = mean_se, geom = "errorbar") +

    xlab("Transcription Factors") + ylab("Average methylation percent") + 
    theme_bw() + labs(color='Transcription Factors(0.10-0.15)') +
    ggtitle("Avg. Methylation across TF binding sites") + 
    theme(
	plot.title=element_text(size=14, face="bold", hjust = 0.6),
	legend.title=element_text(face="bold")
	) +coord_flip() + scale_y_continuous(limits = quantile(meth_df_1$meth_percent, c(0.1, 1.0)))


print(violinplot)


dev.off()



output_file_name <- paste0("~/Dropbox", "/", "category_4_tf_barplot.pdf")				
pdf(output_file_name)

barplot <- ggplot(meth_df_4, aes(x=TF_name, y=meth_percent, color=TF_name)) + 
    geom_boxplot(width = 0.3, outlier.shape = NA) + geom_jitter(size=0.01, color="black", alpha=0.02) +
    #geom_point(colour = "blue", position = "jitter", alpha=0.5) +
    stat_boxplot( aes(TF_name, meth_percent), 
    geom='errorbar', linetype=1, width=0.5)+  #whiskers
	#geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25, color="black") + 
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    #stat_summary(fun.y=mean, geom="point", size=2) + 
    #stat_summary(fun.data = mean_se, geom = "errorbar") +

    xlab("Transcription Factors") + ylab("Average methylation percent") + 
    theme_bw() + labs(color='Transcription Factors(>=0.15)') +
    ggtitle("Avg. Methylation across TF binding sites") + 
    theme(
	plot.title=element_text(size=14, face="bold", hjust = 0.6),
	legend.title=element_text(face="bold")
	) + coord_flip() + scale_y_continuous(limits = quantile(meth_df_1$meth_percent, c(0.1, 1.0)))

print(barplot)

dev.off()


output_file_name <- paste0("~/Dropbox", "/", "category_4_tf_violinplot.pdf")				
pdf(output_file_name)

violinplot <- ggplot(meth_df_4, aes(x=TF_name, y=meth_percent, color=TF_name)) +
    geom_violin() + 
    geom_point(size = 0.05, position = "jitter", alpha=0.1) +
    #stat_boxplot( aes(TF_name, meth_percent), 
    #geom='errorbar', linetype=1, width=0.5)+  #whiskers
	#geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25, color="black") + 
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    #stat_summary(fun.y=mean, geom="point", size=2) + 
    #stat_summary(fun.data = mean_se, geom = "errorbar") +

    xlab("Transcription Factors") + ylab("Average methylation percent") + 
    theme_bw() + labs(color='Transcription Factors(>=0.15)') +
    ggtitle("Avg. Methylation across TF binding sites") + 
    theme(
	plot.title=element_text(size=14, face="bold", hjust = 0.6),
	legend.title=element_text(face="bold")
	) +coord_flip() + scale_y_continuous(limits = quantile(meth_df_1$meth_percent, c(0.1, 1.0)))


print(violinplot)


dev.off()





### Split the set2 into 3 groups:
####################
####################


output_file_name <- paste0("~/Dropbox", "/", "category_2_tf_barplot_group1.pdf")				
pdf(output_file_name)

barplot <- ggplot(meth_df_2_1, aes(x=TF_name, y=meth_percent, color=TF_name)) + 
    geom_boxplot(outlier.shape = NA) + geom_jitter(size=0.01, color="black", alpha=0.02) +
    #geom_point(colour = "blue", position = "jitter", alpha=0.5) +
    stat_boxplot( aes(TF_name, meth_percent), 
    geom='errorbar', linetype=1, width=0.5)+  #whiskers
	#geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25, color="black") + 
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    #stat_summary(fun.y=mean, geom="point", size=2) + 
    #stat_summary(fun.data = mean_se, geom = "errorbar") +

    xlab("Transcription Factors ") + ylab("Average methylation percent") + 
    theme_bw() + labs(color='Transcription Factors(0.05-0.10)') +
    ggtitle("Avg. Methylation across TF binding sites") + 
    theme(
	plot.title=element_text(size=14, face="bold", hjust = 0.6),
	legend.title=element_text(face="bold")
	) + coord_flip() 
print(barplot)

dev.off()


output_file_name <- paste0("~/Dropbox", "/", "category_2_tf_violinplot_group1.pdf")				
pdf(output_file_name)

violinplot <- ggplot(meth_df_2_1, aes(x=TF_name, y=meth_percent, color=TF_name)) +
    geom_violin() + 
    geom_point(size = 0.05, position = "jitter", alpha=0.1) +
    #stat_boxplot( aes(TF_name, meth_percent), 
    #geom='errorbar', linetype=1, width=0.5)+  #whiskers
	#geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25, color="black") + 
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    #stat_summary(fun.y=mean, geom="point", size=2) + 
    #stat_summary(fun.data = mean_se, geom = "errorbar") +

    xlab("Transcription Factors") + ylab("Average methylation percent") + 
    theme_bw() + labs(color='Transcription Factors(0.05-0.10)') +
    ggtitle("Avg. Methylation across TF binding sites") + 
    theme(
	plot.title=element_text(size=14, face="bold", hjust = 0.6),
	legend.title=element_text(face="bold")
	) +coord_flip() + scale_y_continuous(limits = quantile(meth_df_1$meth_percent, c(0.1, 1.0)))


print(violinplot)


dev.off()





output_file_name <- paste0("~/Dropbox", "/", "category_2_tf_barplot_group2.pdf")				
pdf(output_file_name)

barplot <- ggplot(meth_df_2_2, aes(x=TF_name, y=meth_percent, color=TF_name)) + 
    geom_boxplot(outlier.shape = NA) + geom_jitter(size=0.01, color="black", alpha=0.02) +
    #geom_point(colour = "blue", position = "jitter", alpha=0.5) +
    stat_boxplot( aes(TF_name, meth_percent), 
    geom='errorbar', linetype=1, width=0.5)+  #whiskers
	#geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25, color="black") + 
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    #stat_summary(fun.y=mean, geom="point", size=2) + 
    #stat_summary(fun.data = mean_se, geom = "errorbar") +

    xlab("Transcription Factors ") + ylab("Average methylation percent") + 
    theme_bw() + labs(color='Transcription Factors(0.05-0.10)') +
    ggtitle("Avg. Methylation across TF binding sites") + 
    theme(
	plot.title=element_text(size=14, face="bold", hjust = 0.6),
	legend.title=element_text(face="bold")
	) + coord_flip() 
print(barplot)

dev.off()


output_file_name <- paste0("~/Dropbox", "/", "category_2_tf_violinplot_group2.pdf")				
pdf(output_file_name)

violinplot <- ggplot(meth_df_2_2, aes(x=TF_name, y=meth_percent, color=TF_name)) +
    geom_violin() + 
    geom_point(size = 0.05, position = "jitter", alpha=0.1) +
    #stat_boxplot( aes(TF_name, meth_percent), 
    #geom='errorbar', linetype=1, width=0.5)+  #whiskers
	#geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25, color="black") + 
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    #stat_summary(fun.y=mean, geom="point", size=2) + 
    #stat_summary(fun.data = mean_se, geom = "errorbar") +

    xlab("Transcription Factors") + ylab("Average methylation percent") + 
    theme_bw() + labs(color='Transcription Factors(0.05-0.10)') +
    ggtitle("Avg. Methylation across TF binding sites") + 
    theme(
	plot.title=element_text(size=14, face="bold", hjust = 0.6),
	legend.title=element_text(face="bold")
	) +coord_flip() + scale_y_continuous(limits = quantile(meth_df_1$meth_percent, c(0.1, 1.0)))


print(violinplot)


dev.off()




output_file_name <- paste0("~/Dropbox", "/", "category_2_tf_barplot_group3.pdf")				
pdf(output_file_name)

barplot <- ggplot(meth_df_2_3, aes(x=TF_name, y=meth_percent, color=TF_name)) + 
    geom_boxplot(outlier.shape = NA) + geom_jitter(size=0.01, color="black", alpha=0.02) +
    #geom_point(colour = "blue", position = "jitter", alpha=0.5) +
    stat_boxplot( aes(TF_name, meth_percent), 
    geom='errorbar', linetype=1, width=0.5)+  #whiskers
	#geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25, color="black") + 
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    #stat_summary(fun.y=mean, geom="point", size=2) + 
    #stat_summary(fun.data = mean_se, geom = "errorbar") +

    xlab("Transcription Factors ") + ylab("Average methylation percent") + 
    theme_bw() + labs(color='Transcription Factors(0.05-0.10)') +
    ggtitle("Avg. Methylation across TF binding sites") + 
    theme(
	plot.title=element_text(size=14, face="bold", hjust = 0.6),
	legend.title=element_text(face="bold")
	) + coord_flip() 
print(barplot)

dev.off()


output_file_name <- paste0("~/Dropbox", "/", "category_2_tf_violinplot_group3.pdf")				
pdf(output_file_name)

violinplot <- ggplot(meth_df_2_3, aes(x=TF_name, y=meth_percent, color=TF_name)) +
    geom_violin() + 
    geom_point(size = 0.05, position = "jitter", alpha=0.1) +
    #stat_boxplot( aes(TF_name, meth_percent), 
    #geom='errorbar', linetype=1, width=0.5)+  #whiskers
	#geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25, color="black") + 
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    #stat_summary(fun.y=mean, geom="point", size=2) + 
    #stat_summary(fun.data = mean_se, geom = "errorbar") +

    xlab("Transcription Factors") + ylab("Average methylation percent") + 
    theme_bw() + labs(color='Transcription Factors(0.05-0.10)') +
    ggtitle("Avg. Methylation across TF binding sites") + 
    theme(
	plot.title=element_text(size=14, face="bold", hjust = 0.6),
	legend.title=element_text(face="bold")
	) +coord_flip() + scale_y_continuous(limits = quantile(meth_df_1$meth_percent, c(0.1, 1.0)))


print(violinplot)


dev.off()




library(gplots)
read_df <- read.table("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/wgbs_tf_intersect/final_wgbs_tf_ideas_intersect_heatmap_data.bed", sep="\t", header=TRUE)
#heatmap_df <- as.matrix(read_df)
rnames <- read_df[,1]
mat_data <- data.matrix(read_df[,2:ncol(read_df)])
rownames(mat_data) <- rnames
rownames(mat_data)
colnames(mat_data)

output_file_name <- paste0("~/Dropbox", "/", "wgbs_tf_ideas_state_heatmap.pdf")				
pdf(output_file_name)

heatmap.2(mat_data,
  main = "Methylation Distribution", 
  xlab = "IDEAS States",
  ylab = "Transcription Factors",
  col=greenred(10), 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.3,
  cexCol = 0.6,
  na.color="grey"
  #scale="column", tracecol="#303030"
  )    

dev.off()



read_file <- fread("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/wgbs_tf_intersect/final_wgbs_tf_chromHMM_intersect_heatmap_data.bed", sep="\t", header=TRUE)
#heatmap_df <- as.matrix(read_df)
read_df <- as.data.frame(read_file)
rnames <- read_df[,1]
data <- read_df[,2:ncol(read_df)]
mat_data <- as.matrix(data)
rownames(mat_data) <- rnames
rownames(mat_data)
colnames(mat_data)


output_file_name <- paste0("~/Dropbox", "/", "wgbs_tf_chromHMM_state_heatmap.pdf")				
pdf(output_file_name)

heatmap.2(mat_data,
  main = "Methylation Distribution", 
  xlab = "ChromHMM States",
  ylab = "Transcription Factors",
  col=greenred(10), 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.3,
  cexCol = 0.6,
  na.color="grey"
  #scale="column", tracecol="#303030"
  )    

dev.off()


#### tf-tf cobinding analyis:
library(gplots)
#read_file <- fread("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/tf_cobinding_analyis/final_cobinding_peak_count.bed", sep="\t", header=TRUE)
read_file <- fread("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/tf_cobinding_analyis/final_cobinding_peak_count_heatmap.bed", sep="\t", header=TRUE)
#heatmap_df <- as.matrix(read_df)
read_df <- as.data.frame(read_file)
rnames <- read_df[,1]
data <- read_df[,2:ncol(read_df)]
mat_data <- as.matrix(data)
rownames(mat_data) <- rnames
rownames(mat_data)
colnames(mat_data)


output_file_name <- paste0("~/Dropbox", "/", "Tf_cobinding_heatmap.pdf")				
pdf(output_file_name)

heatmap.2(mat_data,
  main = "TF Cobind Distribution", 
  xlab = "TF cobind",
  ylab = "TF name",
  col=greenred(10), 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.35,
  cexCol = 0.35,
  na.color="grey"
  #scale="column", tracecol="#303030"
  )    
dev.off()




#### Unique hepg2 ideas tf analysis:


#### tf-tf cobinding analyis:
library(gplots)
#read_file <- fread("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/tf_cobinding_analyis/final_cobinding_peak_count.bed", sep="\t", header=TRUE)
read_file <- fread("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/unique_hepg2_analysis/final_tf_ideas_peak_count_heatmap.bed", sep="\t", header=TRUE)
#heatmap_df <- as.matrix(read_df)
read_df <- as.data.frame(read_file)
rnames <- read_df[,1]
data <- read_df[,2:ncol(read_df)]
mat_data <- as.matrix(data)
rownames(mat_data) <- rnames
rownames(mat_data)
colnames(mat_data)


output_file_name <- paste0("~/Dropbox", "/", "hepg2_specific_tf_distribution_heatmap.pdf")				
pdf(output_file_name)

heatmap.2(mat_data,
  main = "HepG2 specific TF Distribution", 
  xlab = "TF name",
  ylab = "IDEAS state",
  col=greenred(10), 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.35,
  cexCol = 0.35,
  na.color="grey"
  #scale="column", tracecol="#303030"
  )    
dev.off()


#### tf-tf cobinding analyis:
library(gplots)
#read_file <- fread("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/tf_cobinding_analyis/final_cobinding_peak_count.bed", sep="\t", header=TRUE)
read_file <- fread("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/unique_hepg2_analysis/final_tf_ideas_peak_count_heatmap_qval.bed", sep="\t", header=TRUE)
#heatmap_df <- as.matrix(read_df)
read_df <- as.data.frame(read_file)
rnames <- read_df[,1]
data <- read_df[,2:ncol(read_df)]
mat_data <- as.matrix(data)
rownames(mat_data) <- rnames
rownames(mat_data)
colnames(mat_data)


output_file_name <- paste0("~/Dropbox", "/", "hepg2_specific_tf_distribution_heatmap_qval.pdf")				
pdf(output_file_name)

heatmap.2(mat_data,
  main = "HepG2 specific TF Distribution", 
  xlab = "IDEAS state",
  ylab = "TF name enriched with qvalue: < .001",
  col=greenred(10), 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.35,
  cexCol = 0.35,
  na.color="grey"
  #scale="column", tracecol="#303030"
  )    

dev.off()


#### tf-tf cobinding analyis:
library(gplots)
#read_file <- fread("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/tf_cobinding_analyis/final_cobinding_peak_count.bed", sep="\t", header=TRUE)
read_file <- fread("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/unique_hepg2_analysis/final_tf_ideas_peak_count_heatmap_2.bed", sep="\t", header=TRUE)
#heatmap_df <- as.matrix(read_df)
read_df <- as.data.frame(read_file)
rnames <- read_df[,1]
data <- read_df[,2:ncol(read_df)]
mat_data <- as.matrix(data)
rownames(mat_data) <- rnames
rownames(mat_data)
colnames(mat_data)


output_file_name <- paste0("~/Dropbox", "/", "hepg2_specific_tf_distribution_heatmap.pdf")				
pdf(output_file_name)

heatmap.2(mat_data,
  main = "HepG2 specific TF Distribution", 
  xlab = "TF name",
  ylab = "IDEAS state",
  col=greenred(10), 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.35,
  cexCol = 0.35,
  na.color="grey"
  #scale="column", tracecol="#303030"
  )    
dev.off()









### Sanity check of the match:

# > do.call(rbind, meth_df_1)
#            TF_name index peak_chrom peak_start  peak_end meth_percent
#      1: ATF1[FLAG]     0       chr1     713825    714385  0.001302083
#      2: ATF1[FLAG]     1       chr1     936060    936620  0.013035382
#      3: ATF1[FLAG]     2       chr1     941647    942207  0.047511312
#      4: ATF1[FLAG]     3       chr1     941850    942410  0.230179028
#      5: ATF1[FLAG]     4       chr1     948541    949101  0.002164502
#     ---                                                              
# 288025:       CREM 17473       chrX  153626602 153627132  0.008830022
# 288026:       CREM 17474       chrX  153770003 153770533  0.000000000
# 288027:       CREM 17475       chrX  154027974 154028504  0.023809524
# 288028:       CREM 17476       chrX  154216170 154216700  0.000000000
# 288029:       CREM 17477       chrX  155110654 155111184  0.241086587


# > unique(do.call(rbind, meth_df_1)$TF_name)
#  [1] "ATF1[FLAG]"      "TCF25[FLAG]"     "HBP1[FLAG]"      "KDM3A[FLAG]"    
#  [5] "RERE[FLAG]"      "KMT2B[FLAG]"     "ZNF792[FLAG]"    "MLX[FLAG]"      
#  [9] "KLF9[FLAG]"      "DRAP1[FLAG]"     "ZSCAN9[FLAG]"    "KLF6_v2[FLAG]"  
# [13] "KAT7[FLAG]"      "MXD3_v1[FLAG]"   "HMGXB4[FLAG]"    "DMAP1[FLAG]"    
# [17] "MXD4[FLAG]"      "SP5[FLAG]"       "TFDP1[FLAG]"     "RFXANK[FLAG]"   
# [21] "HMG20B_v2[FLAG]" "SOX13"           "GATA4"           "CBX1"           
# [25] "KLF10"           "TGIF2[FLAG]"     "GMEB2[FLAG]"     "RUVBL1"         
# [29] "ETV4"            "CREM"           
# > unique(do.call(rbind, meth_df_1)$TF_name) == target_vector_1
#  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
# [16] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
# > target_vector_1
#  [1] "ATF1[FLAG]"      "TCF25[FLAG]"     "HBP1[FLAG]"      "KDM3A[FLAG]"    
#  [5] "RERE[FLAG]"      "KMT2B[FLAG]"     "ZNF792[FLAG]"    "MLX[FLAG]"      
#  [9] "KLF9[FLAG]"      "DRAP1[FLAG]"     "ZSCAN9[FLAG]"    "KLF6_v2[FLAG]"  
# [13] "KAT7[FLAG]"      "MXD3_v1[FLAG]"   "HMGXB4[FLAG]"    "DMAP1[FLAG]"    
# [17] "MXD4[FLAG]"      "SP5[FLAG]"       "TFDP1[FLAG]"     "RFXANK[FLAG]"   
# [21] "HMG20B_v2[FLAG]" "SOX13"           "GATA4"           "CBX1"           
# [25] "KLF10"           "TGIF2[FLAG]"     "GMEB2[FLAG]"     "RUVBL1"         
# [29] "ETV4"            "CREM"    


# In [284]: df_test
# Out[284]: 
#    name   Tf  value
# 0   Foo  fox     10
# 1  Baar  Sp1     15
# 2   Foo  FOX     10
# 3  Baar  sp1     25


# In [276]: df_test.pivot_table(index="name", columns='Tf', values='value')
# Out[276]: 
# Tf     FOX   Sp1   fox   sp1
# name                        
# Baar   NaN  15.0   NaN  25.0
# Foo   10.0   NaN  10.0   NaN

# In [277]: df_test.pivot(index="name", columns='Tf', values='value')
# Out[277]: 
# Tf     FOX   Sp1   fox   sp1
# name                        
# Baar   NaN  15.0   NaN  25.0
# Foo   10.0   NaN  10.0   NaN

# In [278]: df_test.set_index(["name","Tf"]).stack().unstack()
# Out[278]: 
#           value
# name Tf        
# Baar Sp1     15
#      sp1     25
# Foo  FOX     10
#      fox     10

# In [279]: df_test.set_index(["name","Tf"]).unstack()
# Out[279]: 
#      value                  
# Tf     FOX   Sp1   fox   sp1
# name                        
# Baar   NaN  15.0   NaN  25.0
# Foo   10.0   NaN  10.0   NaN

# In [280]: df_test.set_index(["name","Tf"], append=True).stack()
# Out[280]: 
#    name  Tf        
# 0  Foo   fox  value    10
# 1  Baar  Sp1  value    15
# 2  Foo   FOX  value    10
# 3  Baar  sp1  value    25
# dtype: int64

# In [281]: df_test.set_index(["name","Tf"], append=True).unstack()
# Out[281]: 
#        value                  
# Tf       FOX   Sp1   fox   sp1
#   name                        
# 0 Foo    NaN   NaN  10.0   NaN
# 1 Baar   NaN  15.0   NaN   NaN
# 2 Foo   10.0   NaN   NaN   NaN
# 3 Baar   NaN   NaN   NaN  25.0

# MATPLOTLIB
# fig, ax = plt.subplots(1, 1, figsize=(7.5, 7.5))
 
# def scatter(group):
#     plt.plot(group['petalLength'],
#              group['petalWidth'],
#              'o', label=group.name)
 
# df.groupby('species').apply(scatter)
 
# ax.set(xlabel='Petal Length',
#        ylabel='Petal Width',
#        title='Petal Width v. Length -- by Species')
 
# ax.legend(loc=2)


