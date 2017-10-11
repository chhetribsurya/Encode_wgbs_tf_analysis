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



######################################################################
#### Average methylation distribution of all unique TFs (208 factors):
######################################################################

#avg_df <- fread("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/wgbs_tf_intersect/final_wgbs_tf_intersect_average_meth.bed", sep="\t", header=TRUE)
#avg_df <- fread("~/Dropbox/encode_3/wgbs_tf_intersect_total/files/avg_meth_distribution/final_wgbs_tf_intersect_average_meth_unique_TFs.bed", sep="\t", header=TRUE)
#avg_df <- fread("~/Dropbox/encode_3/wgbs_tf_intersect_total/files/avg_meth_distribution/final_wgbs_tf_intersect_average_meth_all_TFs.bed", sep="\t", header=TRUE)

xls_df <- read.xlsx2("~/Dropbox/for_chris/Encode_full_hepg2_datasets_DBF_CR.xlsx", sheetIndex=4)
avg_df <- fread("~/Dropbox/encode_3/fimo_motif_based_analysis/wgbs_tf_intersect_avg_of_motifs/files/final_wgbs_tf_intersect_average_meth_unique_TFs.txt", sep="\t", header=TRUE)

names(avg_df) <- c("TF_name", "cpg_density(cpgs per motif)", "meth_percent", "total_cpg_sites", "total_motif_sites")
target_vector <- avg_df$TF_name
xls_df_ordered <- xls_df[match(target_vector, xls_df$Target),]


# check if the order b/w 2 columns are equal: 
if(all(avg_df$TF_name == xls_df_ordered$Target))
{
   print("Column A(tf_name) and B(target) are identical")
}

xls_df_ordered <- xls_df_ordered %>% select(Target, Lab, Category)
avg_df$tf_category <-  xls_df_ordered$Category

cr_df <- avg_df %>% as.data.frame %>% filter(tf_category == "CR/CF")
dbf_df <- avg_df %>% as.data.frame %>% filter(tf_category == "DBF")
all_tf_df <- avg_df %>% as.data.frame

### eliminate BMI1 and ZNF274, the ones with 301 and 198 low peak counts; else just comment the line to include all the tfs. 
# all_tf_df <- all_tf_df %>% filter(!(TF_name == "BMI1_human"| TF_name == "ZNF274_human"))


#output_file_name <- paste0("~/Dropbox", "/", "all_tfs_methylation_distribution_BMI1_znf274_filtered.pdf")				
output_file_name <- paste0("~/Dropbox", "/", "average_methylation_distribution_for_unique_tfs.pdf")				
pdf(output_file_name)

barplot <-ggplot(all_tf_df, aes(x=tf_category, y=meth_percent, fill=tf_category)) + 
	#geom_point() + 
    geom_boxplot(width=0.15, fill = "grey80", outlier.colour = "red", color="blue") + # outlier.shape = NA to remove the outliers 
	geom_dotplot(data=subset(all_tf_df, meth_percent < 0.20), binaxis='y', stackdir='centerwhole', dotsize=0.2, fill="black") +
	stat_boxplot( aes(tf_category, meth_percent), 
    geom='errorbar', linetype=2, width=0.15, color="blue") +  
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    xlab("Transcription Factor Category") + ylab("Average methylation percent") + theme_bw()+ 
    ggtitle("Average methylation distribution across motifs of 208 TFs") + 
    scale_y_continuous(limits=c(0,1)) +
    theme(
    axis.text.x = element_text(size=10, face="bold"),
	plot.title=element_text(size=14, face="bold", hjust = 0.6)
	)

barplot <- barplot + geom_text(data = subset(all_tf_df, meth_percent > 0.15), aes(x = tf_category, y =meth_percent, label = TF_name), 
	hjust= -0.15, size = 1.5, check_overlap = TRUE)

print(barplot)

dev.off()



######################################################
### Average Methylation distribution of dup and unique with total of 239 TFs :
######################################################

xls_df1 <- read.xlsx2("~/Dropbox/for_chris/Encode_full_hepg2_datasets_DBF_CR.xlsx", sheetIndex=4)
xls_df2 <- read.xlsx2("~/Dropbox/for_chris/Encode_full_hepg2_datasets_DBF_CR.xlsx", sheetIndex=5)
xls_df <-  rbind(xls_df1, xls_df2)

#avg_df <- fread("~/Dropbox/encode_3/wgbs_tf_intersect_total/files/avg_meth_distribution/final_wgbs_tf_intersect_average_meth_all_TFs.bed", sep="\t", header=TRUE)
avg_df <- fread("~/Dropbox/encode_3/fimo_motif_based_analysis/wgbs_tf_intersect_avg_of_motifs/files/final_wgbs_tf_intersect_average_meth_total.txt", sep="\t", header=TRUE)

names(avg_df) <- c("TF_name", "meth_percent")

target_vector <- avg_df$TF_name
xls_df_ordered <- xls_df[match(target_vector, xls_df$Target),]

# check if the order b/w 2 columns are equal: 
if(all(avg_df$TF_name == xls_df_ordered$Target))
{
   print("Column A(tf_name) and B(target) are identical")
}

xls_df_ordered <- xls_df_ordered %>% select(Target, Lab, Category)
avg_df$tf_category <-  xls_df_ordered$Category

cr_df <- avg_df %>% as.data.frame %>% filter(tf_category == "CR/CF")
dbf_df <- avg_df %>% as.data.frame %>% filter(tf_category == "DBF")
all_tf_df <- avg_df %>% as.data.frame


#output_file_name <- paste0("~/Dropbox", "/", "all_tfs_methylation_distribution_BMI1_znf274_filtered.pdf")              
output_file_name <- paste0("~/Dropbox", "/", "final_methylation_distribution_all_dup_tfs_239.pdf")              
pdf(output_file_name)

barplot <-ggplot(all_tf_df, aes(x=tf_category, y=meth_percent, fill=tf_category)) + 
    #geom_point() + 
    geom_boxplot(width=0.15, fill = "grey80", outlier.colour = "red", color="blue") + # outlier.shape = NA to remove the outliers 
    geom_dotplot(data=subset(all_tf_df, meth_percent < 0.15), binaxis='y', stackdir='centerwhole', dotsize=0.2, fill="black") +
    stat_boxplot( aes(tf_category, meth_percent), 
    geom='errorbar', linetype=2, width=0.15, color="blue") +  
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    xlab("Transcription Factor Category") + ylab("Average methylation percent") + theme_bw()+ 
    ggtitle("Average Methylation Distribution of 239 TFs") + 
    scale_y_continuous(limits=c(0,1)) +
    theme(
    axis.text.x = element_text(size=10, face="bold"),
    plot.title=element_text(size=14, face="bold", hjust = 0.6)
    )

barplot <- barplot + geom_text(data = subset(all_tf_df, meth_percent > 0.15), aes(x = tf_category, y =meth_percent, label = TF_name), 
    hjust= -0.15, size = 2, check_overlap = TRUE)

print(barplot)

dev.off()




###################
###################
###################

######################################################
### Sorted detailed peak Methylation distribution for 208 factors:
######################################################

### Peak methylation distribution for unique TFs 208 factors:
### Peak methylation distribution for 3 categories: DBF + CR/CF, DBF, CR 
### Both DBF and CR/CF combined methylation distribution for each peaks:

#xls_df <- read.xlsx2("~/Dropbox/for_chris/Encode_full_hepg2_datasets_DBF_CR.xlsx", sheetIndex=4)
#meth_tf_df <- fread("~/Dropbox/encode_3/wgbs_tf_intersect_total/files/avg_meth_distribution/final_wgbs_tf_intersect_unique_TFs.bed", sep="\t", header=TRUE)
avg_df <- fread("~/Dropbox/encode_3/fimo_motif_based_analysis/wgbs_tf_intersect_avg_of_motifs/files/final_wgbs_tf_intersect_average_meth_unique_TFs.txt", sep="\t", header=TRUE)
meth_tf_df <- fread("~/Dropbox/encode_3/fimo_motif_based_analysis/wgbs_tf_intersect_avg_of_motifs/files/final_wgbs_tf_motif_intersect_unique_TFs.txt", sep="\t", header=TRUE)
names(meth_tf_df) <- c("TF_name", "index", "peak_chrom", "peak_start", "peak_end", "meth_percent")

target_vector <- meth_tf_df$TF_name
xls_df_ordered_peaks <- xls_df[match(target_vector, xls_df$Target),]

# check if the order b/w 2 columns are equal: 
if(all(meth_tf_df$TF_name == xls_df_ordered_peaks$Target))
{
   print("Column A(tf_name) and B(target) are identical")
}

xls_df_ordered_peaks <- xls_df_ordered_peaks %>% select(Target, Lab, Category)
meth_tf_df$tf_category <-  xls_df_ordered_peaks$Category

## Now, sort the order of TF wrt methylation value:
sorted_avg_df <- avg_df[order(avg_df$meth_percent),]
factor(sorted_avg_df$TF_name) %>% levels
sorted_avg_df$TF_name <- factor(sorted_avg_df$TF_name, levels=sorted_avg_df$TF_name)

cr_df <- meth_tf_df %>% as.data.frame %>% filter(tf_category == "CR/CF")
dbf_df <- meth_tf_df %>% as.data.frame %>% filter(tf_category == "DBF")
all_tf_df <- meth_tf_df %>% as.data.frame

# change the order of factor in TF Category:
meth_tf_df$tf_category <- factor(meth_tf_df$tf_category, levels = c("DBF", "CR/CF" ))
meth_tf_df$TF_name <- factor(meth_tf_df$TF_name, levels=sorted_avg_df$TF_name)

dbf_df$TF_name <- factor(dbf_df$TF_name, levels=sorted_avg_df$TF_name)
cr_df$TF_name <- factor(cr_df$TF_name, levels=sorted_avg_df$TF_name)

output_file_name <- paste0("~/Dropbox", "/", "motif_based_dbf_cr_meth_boxplot_for_208_TFs.pdf")				
pdf(output_file_name)

# z <- c(0.15, 0.15)
test <- data.frame(X = c("DBF", "CR/CF" ), Z = c(0.15, 0.15))

barplot <- ggplot(meth_tf_df, aes(x=TF_name, y=meth_percent, fill=tf_category)) + 
    geom_boxplot(aes(fill=tf_category),outlier.shape = NA) +
    geom_hline(data= test, aes(yintercept = Z), colour="red", linetype = "longdash" ) +
    #geom_jitter(size=0.01, color="black", alpha=0.01) +
    #geom_point(colour = "blue", position = "jitter", aeslpha=0.5) +
    #stat_boxplot( aes(TF_name, meth_percent), 
    #geom='errorbar', linetype=1, width=0.5) +  #whiskers
	#geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25, color="black") + 
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    #stat_summary(fun.y=mean, geom="point", size=2) + 
    #stat_summary(fun.data = mean_se, geom = "errorbar") +

    xlab("Transcription Factors") + ylab("Average methylation percent") + 
    theme_bw() + 
    ggtitle("Avg. methylation distribution across TF binding sites") + 
    theme(
    axis.text.y = element_text(size=2.2),
	plot.title=element_text(size=14, face="bold", hjust = 0.6),
	legend.title=element_text(face="bold")
	) + 
	coord_flip() + 
	guides(fill=guide_legend(title="TF Category")) + 
	scale_y_continuous(limits = quantile(meth_tf_df$meth_percent, c(0.0, 1.0))) +
	facet_wrap(~tf_category)

print(barplot)

dev.off()


### the same plot, instead of facet wrap plotted individually:
### Only DNA binding factors peak methylation distribution:

cr_df <- meth_tf_df %>% as.data.frame %>% filter(tf_category == "CR/CF")
dbf_df <- meth_tf_df %>% as.data.frame %>% filter(tf_category == "DBF")
all_tf_df <- meth_tf_df %>% as.data.frame

output_file_name <- paste0("~/Dropbox", "/", "motif_based_dbf_meth_boxplot_for_208_TFs.pdf")				
pdf(output_file_name)

barplot <- ggplot(dbf_df, aes(x=TF_name, y=meth_percent, fill=tf_category)) + 
    geom_boxplot(aes(fill=tf_category),outlier.shape = NA) +
    geom_hline(aes(yintercept = 0.15), colour="red", linetype = "longdash" ) +
    #geom_jitter(size=0.01, color="black", alpha=0.01) +
    #geom_point(colour = "blue", position = "jitter", aeslpha=0.5) +
    #stat_boxplot( aes(TF_name, meth_percent), 
    #geom='errorbar', linetype=1, width=0.5) +  #whiskers
	#geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25, color="black") + 
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    #stat_summary(fun.y=mean, geom="point", size=2) + 
    #stat_summary(fun.data = mean_se, geom = "errorbar") +

    xlab("Transcription Factors") + ylab("Average methylation percent") + 
    theme_bw() + 
    ggtitle("Avg. methylation distribution across TF binding sites") + 
    theme(
    axis.text.y = element_text(size=2.3),
	plot.title=element_text(size=14, face="bold", hjust = 0.6),
	legend.title=element_text(face="bold")
	) + 
	coord_flip() + 
	guides(fill=guide_legend(title="TF Category")) + 
	scale_y_continuous(limits = quantile(dbf_df$meth_percent, c(0.0, 1.0))) 

print(barplot)

dev.off()


### Only Chromatin Regulators peak methylation distribution:

cr_df <- meth_tf_df %>% as.data.frame %>% filter(tf_category == "CR/CF")
dbf_df <- meth_tf_df %>% as.data.frame %>% filter(tf_category == "DBF")
all_tf_df <- meth_tf_df %>% as.data.frame

output_file_name <- paste0("~/Dropbox", "/", "motif_based_cr_meth_boxplot_for_208_TFs.pdf")				
pdf(output_file_name)

barplot <- ggplot(cr_df, aes(x=TF_name, y=meth_percent, fill=tf_category)) + 
    geom_boxplot(aes(fill=tf_category),outlier.shape = NA) +
    geom_hline(aes(yintercept = 0.15), colour="red", linetype = "longdash" ) +
    #geom_text(aes(x=5, label="0.15", y=0.20)) + 
    #geom_jitter(size=0.01, color="black", alpha=0.01) +
    #geom_point(colour = "blue", position = "jitter", aeslpha=0.5) +
    #stat_boxplot( aes(TF_name, meth_percent), 
    #geom='errorbar', linetype=1, width=0.5) +  #whiskers
	#geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25, color="black") + 
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    #stat_summary(fun.y=mean, geom="point", size=2) + 
    #stat_summary(fun.data = mean_se, geom = "errorbar") +

    xlab("Transcription Factors") + ylab("Average methylation percent") + 
    theme_bw() + 
    ggtitle("Avg. methylation distribution across TF binding sites") + 
    theme(
    axis.text.y = element_text(size=5.5),
	plot.title=element_text(size=14, face="bold", hjust = 0.6),
	legend.title=element_text(face="bold")
	) + 
	coord_flip() + 
	guides(fill=guide_legend(title="TF Category")) + 
	scale_y_continuous(limits = quantile(cr_df$meth_percent, c(0.0, 1.0))) 

print(barplot)

dev.off()



##########################################################
### Sorted detailed Peak methylation distribution for all TFs 239 factors:
##########################################################

### Peak methylation distribution for 3 categories: DBF + CR/CF, DBF, CR 
### Both DBF and CR/CF combined methylation distribution for each peaks:

xls_df1 <- read.xlsx2("~/Dropbox/for_chris/Encode_full_hepg2_datasets_DBF_CR.xlsx", sheetIndex=4)
xls_df2 <- read.xlsx2("~/Dropbox/for_chris/Encode_full_hepg2_datasets_DBF_CR.xlsx", sheetIndex=5)
xls_df <-  rbind(xls_df1, xls_df2)

#avg_df <- fread("~/Dropbox/encode_3/wgbs_tf_intersect_total/files/avg_meth_distribution/final_wgbs_tf_intersect_average_meth_all_TFs.bed", sep="\t", header=TRUE)
#meth_tf_df <- fread("~/Dropbox/encode_3/wgbs_tf_intersect_total/files/avg_meth_distribution/final_wgbs_tf_intersect_all_TFs.bed", sep="\t", header=TRUE)
avg_df <- fread("~/Dropbox/encode_3/fimo_motif_based_analysis/wgbs_tf_intersect_avg_of_motifs/files/final_wgbs_tf_intersect_average_meth_total.txt", sep="\t", header=TRUE)
meth_tf_df <- fread("~/Dropbox/encode_3/fimo_motif_based_analysis/wgbs_tf_intersect_avg_of_motifs/files/final_wgbs_tf_motif_intersect_total_TFs.txt", sep="\t", header=TRUE)
names(meth_tf_df) <- c("TF_name", "index", "peak_chrom", "peak_start", "peak_end", "meth_percent")

target_vector <- meth_tf_df$TF_name
xls_df_ordered_peaks <- xls_df[match(target_vector, xls_df$Target),]

# check if the order b/w 2 columns are equal: 
if(all(meth_tf_df$TF_name == xls_df_ordered_peaks$Target))
{
   print("Column A(tf_name) and B(target) are identical")
}

xls_df_ordered_peaks <- xls_df_ordered_peaks %>% select(Target, Lab, Category)
meth_tf_df$tf_category <-  xls_df_ordered_peaks$Category

## Now, sort the order of TF wrt methylation value:
sorted_avg_df <- avg_df[order(avg_df$meth_percent),]
factor(sorted_avg_df$TF_name) %>% levels
sorted_avg_df$TF_name <- factor(sorted_avg_df$TF_name, levels=sorted_avg_df$TF_name)

cr_df <- meth_tf_df %>% as.data.frame %>% filter(tf_category == "CR/CF")
dbf_df <- meth_tf_df %>% as.data.frame %>% filter(tf_category == "DBF")
all_tf_df <- meth_tf_df %>% as.data.frame

# change the order of factor in TF Category:
meth_tf_df$tf_category <- factor(meth_tf_df$tf_category, levels = c("DBF", "CR/CF" ))
meth_tf_df$TF_name <- factor(meth_tf_df$TF_name, levels=sorted_avg_df$TF_name)

dbf_df$TF_name <- factor(dbf_df$TF_name, levels=sorted_avg_df$TF_name)
cr_df$TF_name <- factor(cr_df$TF_name, levels=sorted_avg_df$TF_name)

output_file_name <- paste0("~/Dropbox", "/", "motif_based_dbf_cr_meth_boxplot_for_239_TFs.pdf")               
pdf(output_file_name)

# z <- c(0.15, 0.15)
test <- data.frame(X = c("DBF", "CR/CF" ), Z = c(0.15, 0.15))

barplot <- ggplot(meth_tf_df, aes(x=TF_name, y=meth_percent, fill=tf_category)) + 
    geom_boxplot(aes(fill=tf_category),outlier.shape = NA) +
    geom_hline(data= test, aes(yintercept = Z), colour="red", linetype = "longdash" ) +
    #geom_jitter(size=0.01, color="black", alpha=0.01) +
    #geom_point(colour = "blue", position = "jitter", aeslpha=0.5) +
    #stat_boxplot( aes(TF_name, meth_percent), 
    #geom='errorbar', linetype=1, width=0.5) +  #whiskers
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25, color="black") + 
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    #stat_summary(fun.y=mean, geom="point", size=2) + 
    #stat_summary(fun.data = mean_se, geom = "errorbar") +

    xlab("Transcription Factors") + ylab("Average methylation percent") + 
    theme_bw() + 
    ggtitle("Avg. methylation distribution across TF binding sites") + 
    theme(
    axis.text.y = element_text(size=2.0),
    plot.title=element_text(size=14, face="bold", hjust = 0.6),
    legend.title=element_text(face="bold")
    ) + 
    coord_flip() + 
    guides(fill=guide_legend(title="TF Category")) + 
    scale_y_continuous(limits = quantile(meth_tf_df$meth_percent, c(0.0, 1.0))) +
    facet_wrap(~tf_category)

print(barplot)

dev.off()


### the same plot, instead of facet wrap plotted individually:
### Only DNA binding factors peak methylation distribution:

cr_df <- meth_tf_df %>% as.data.frame %>% filter(tf_category == "CR/CF")
dbf_df <- meth_tf_df %>% as.data.frame %>% filter(tf_category == "DBF")
all_tf_df <- meth_tf_df %>% as.data.frame

output_file_name <- paste0("~/Dropbox", "/", "motif_based_dbf_meth_boxplot_for_239_TFs.pdf")              
pdf(output_file_name)

barplot <- ggplot(dbf_df, aes(x=TF_name, y=meth_percent, fill=tf_category)) + 
    geom_boxplot(aes(fill=tf_category),outlier.shape = NA) +
    geom_hline(aes(yintercept = 0.15), colour="red", linetype = "longdash" ) +
    #geom_jitter(size=0.01, color="black", alpha=0.01) +
    #geom_point(colour = "blue", position = "jitter", aeslpha=0.5) +
    #stat_boxplot( aes(TF_name, meth_percent), 
    #geom='errorbar', linetype=1, width=0.5) +  #whiskers
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25, color="black") + 
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    #stat_summary(fun.y=mean, geom="point", size=2) + 
    #stat_summary(fun.data = mean_se, geom = "errorbar") +

    xlab("Transcription Factors") + ylab("Average methylation percent") + 
    theme_bw() + 
    ggtitle("Avg. methylation distribution across TF binding sites") + 
    theme(
    axis.text.y = element_text(size=2.2),
    plot.title=element_text(size=14, face="bold", hjust = 0.6),
    legend.title=element_text(face="bold")
    ) + 
    coord_flip() + 
    guides(fill=guide_legend(title="TF Category")) + 
    scale_y_continuous(limits = quantile(dbf_df$meth_percent, c(0.0, 1.0))) 

print(barplot)

dev.off()


### Only Chromatin Regulators peak methylation distribution:

cr_df <- meth_tf_df %>% as.data.frame %>% filter(tf_category == "CR/CF")
dbf_df <- meth_tf_df %>% as.data.frame %>% filter(tf_category == "DBF")
all_tf_df <- meth_tf_df %>% as.data.frame

output_file_name <- paste0("~/Dropbox", "/", "motif_based_cr_meth_boxplot_for_239_TFs.pdf")               
pdf(output_file_name)

barplot <- ggplot(cr_df, aes(x=TF_name, y=meth_percent, fill=tf_category)) + 
    geom_boxplot(aes(fill=tf_category),outlier.shape = NA) +
    geom_hline(aes(yintercept = 0.15), colour="red", linetype = "longdash" ) +
    #geom_text(aes(x=5, label="0.15", y=0.20)) + 
    #geom_jitter(size=0.01, color="black", alpha=0.01) +
    #geom_point(colour = "blue", position = "jitter", aeslpha=0.5) +
    #stat_boxplot( aes(TF_name, meth_percent), 
    #geom='errorbar', linetype=1, width=0.5) +  #whiskers
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25, color="black") + 
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    #stat_summary(fun.y=mean, geom="point", size=2) + 
    #stat_summary(fun.data = mean_se, geom = "errorbar") +

    xlab("Transcription Factors") + ylab("Average methylation percent") + 
    theme_bw() + 
    ggtitle("Avg. methylation distribution across TF binding sites") + 
    theme(
    axis.text.y = element_text(size=5.5),
    plot.title=element_text(size=14, face="bold", hjust = 0.6),
    legend.title=element_text(face="bold")
    ) + 
    coord_flip() + 
    guides(fill=guide_legend(title="TF Category")) + 
    scale_y_continuous(limits = quantile(cr_df$meth_percent, c(0.0, 1.0))) 

print(barplot)

dev.off()




#######################################################
### Violinplot for the highest methylation distribution:
#######################################################

xls_df <- read.xlsx2("~/Dropbox/for_chris/Encode_full_hepg2_datasets_DBF_CR.xlsx", sheetIndex=4)
avg_df <- fread("~/Dropbox/encode_3/wgbs_tf_intersect_total/files/avg_meth_distribution/final_wgbs_tf_intersect_average_meth_unique_TFs.bed", sep="\t", header=TRUE)
meth_tf_df <- fread("~/Dropbox/encode_3/wgbs_tf_intersect_total/files/avg_meth_distribution/final_wgbs_tf_intersect_unique_TFs.bed", sep="\t", header=TRUE)

names(avg_df) <- c("TF_name", "meth_percent")
names(meth_tf_df) <- c("TF_name", "index", "peak_chrom", "peak_start", "peak_end", "meth_percent")

target_vector <- meth_tf_df$TF_name
xls_df_ordered_peaks <- xls_df[match(target_vector, xls_df$Target),]

# check if the order b/w 2 columns are equal: 
if(all(meth_tf_df$TF_name == xls_df_ordered_peaks$Target))
{
   print("Column A(tf_name) and B(target) are identical")
}

xls_df_ordered_peaks <- xls_df_ordered_peaks %>% select(Target, Lab, Category)
meth_tf_df$tf_category <-  xls_df_ordered_peaks$Category

## Now, sort the order of TF wrt methylation value:
sorted_avg_df <- avg_df[order(avg_df$meth_percent),]
factor(sorted_avg_df$TF_name) %>% levels
sorted_avg_df$TF_name <- factor(sorted_avg_df$TF_name, levels=sorted_avg_df$TF_name)

cr_df <- meth_tf_df %>% as.data.frame %>% filter(tf_category == "CR/CF")
dbf_df <- meth_tf_df %>% as.data.frame %>% filter(tf_category == "DBF")
all_tf_df <- meth_tf_df %>% as.data.frame

# change the order of factor in TF Category:
meth_tf_df$tf_category <- factor(meth_tf_df$tf_category, levels = c("DBF", "CR/CF" ))
meth_tf_df$TF_name <- factor(meth_tf_df$TF_name, levels=sorted_avg_df$TF_name)

dbf_df$TF_name <- factor(dbf_df$TF_name, levels=sorted_avg_df$TF_name)
cr_df$TF_name <- factor(cr_df$TF_name, levels=sorted_avg_df$TF_name)


################
high_meth_dbf_df <- meth_tf_df

high_meth_df <- sorted_avg_df[sorted_avg_df$meth_percent >= 0.20,]
factor(high_meth_df$TF_name) %>% levels
#high_meth_df$TF_name <-  factor(high_meth_df$TF_name, levels=high_meth_df$TF_name )
high_meth_dbf_df$TF_name <- factor(high_meth_dbf_df$TF_name, levels=high_meth_df$TF_name)
factor(high_meth_dbf_df$TF_name) %>% levels

### The factor level for high_meth_df$TF_name is 13, whereas factor(high_meth_dbf_df$TF_name) level is 208 factors
## thus, while forcing the factor level of just 13 factors, rest of the 195 factors will turn to NA. 
## So, make sure that you eliminate the rows containing null values in TF_name
unique(high_meth_dbf_df$TF_name)
high_meth_dbf_df <- high_meth_dbf_df %>% na.omit()
### filter BMI1_human which had 301 peaks, and ZNF274_human with 198 peaks
high_meth_dbf_df <- high_meth_dbf_df %>% filter(!(TF_name == "BMI1_human"| TF_name == "ZNF274_human"))
unique(high_meth_dbf_df$TF_name)

#################

output_file_name <- paste0("~/Dropbox", "/", "final_dbf_cr_meth_dist_violinplot_filtered_without_scatter.pdf")              
output_file_name <- paste0("~/Dropbox", "/", "final_dbf_cr_meth_dist_violinplot_filtered_with_scatter.pdf")              
pdf(output_file_name)

violinplot <- ggplot(high_meth_dbf_df, aes(x=TF_name, y=meth_percent, fill=tf_category)) +
    geom_violin(fill="white", aes(colour=high_meth_dbf_df$tf_category), trim=FALSE) + 
    geom_point(size = 0.01, position = "jitter", alpha=0.01) +
    #stat_boxplot( aes(TF_name, meth_percent), 
    #geom='errorbar', linetype=1, width=0.5)+  #whiskers
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25, color="black") + 
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    #stat_summary(fun.y=mean, geom="point", size=2) + 
    #stat_summary(fun.data = mean_se, geom = "errorbar") +

    xlab("Transcription Factors") + ylab("Average methylation percent") + 
    theme_bw() + 
    ggtitle("Avg. methylation distribution across TF binding sites") + 
    theme(
    axis.text.y = element_text(size=7.5, face="bold"),
    plot.title=element_text(size=14, face="bold", hjust = 0.6),
    legend.title=element_text(face="bold")
    ) +
    coord_flip() + 
    guides(colour=guide_legend(title="TF Category")) + 
    scale_y_continuous(limits = quantile(high_meth_dbf_df$meth_percent, c(0.0, 1.0))) 

print(violinplot)

dev.off()



######################################################
### Motif based High, Intermediate and Low Methylation distribution for 208 factors:
######################################################

### Peak methylation distribution for unique TFs 208 factors:
### Peak methylation distribution for 3 categories: DBF + CR/CF, DBF, CR 
### Both DBF and CR/CF combined methylation distribution for each peaks:

#xls_df <- read.xlsx2("~/Dropbox/for_chris/Encode_full_hepg2_datasets_DBF_CR.xlsx", sheetIndex=4)
#meth_tf_df <- fread("~/Dropbox/encode_3/wgbs_tf_intersect_total/files/avg_meth_distribution/final_wgbs_tf_intersect_unique_TFs.bed", sep="\t", header=TRUE)
avg_df <- fread("~/Dropbox/encode_3/fimo_motif_based_analysis/wgbs_tf_intersect_avg_of_motifs/files/final_wgbs_tf_intersect_average_meth_unique_TFs.txt", sep="\t", header=TRUE)
meth_tf_df <- fread("~/Dropbox/encode_3/fimo_motif_based_analysis/ideas_cpg_motif_intersect_total/files/final_wgbs_tf_motif_ideas_intersect_HIL_barplot_data.bed", sep="\t", header=TRUE)
#names(meth_tf_df) <- c("TF_name", "index", "peak_chrom", "peak_start", "peak_end", "meth_percent")

target_vector <- meth_tf_df$tf_name
#xls_df_ordered_peaks <- xls_df[match(target_vector, xls_df$Target),]

# check if the order b/w 2 columns are equal: 
#if(all(meth_tf_df$TF_name == xls_df_ordered_peaks$Target))
#{
#   print("Column A(tf_name) and B(target) are identical")
#}

#xls_df_ordered_peaks <- xls_df_ordered_peaks %>% select(Target, Lab, Category)
#meth_tf_df$tf_category <-  xls_df_ordered_peaks$Category

## Now, sort the order of TF wrt methylation value:
sorted_avg_df <- avg_df[order(avg_df$meth_percent),]
factor(sorted_avg_df$TF_name) %>% levels
sorted_avg_df$TF_name <- factor(sorted_avg_df$TF_name, levels=sorted_avg_df$TF_name)

meth_tf_df$tf_name <- factor(meth_tf_df$tf_name, levels=sorted_avg_df$TF_name)
#meth_tf_df$meth_group <- factor(meth_tf_df$meth_group, levels=c("High", "Low", "Intermediate"))
#cr_df <- meth_tf_df %>% as.data.frame %>% filter(tf_category == "CR/CF")
#dbf_df <- meth_tf_df %>% as.data.frame %>% filter(tf_category == "DBF")
all_tf_df <- meth_tf_df %>% as.data.frame

# change the order of factor in TF Category:
#meth_tf_df$tf_category <- factor(meth_tf_df$tf_category, levels = c("DBF", "CR/CF" ))
#meth_tf_df$TF_name <- factor(meth_tf_df$TF_name, levels=sorted_avg_df$TF_name)

#dbf_df$TF_name <- factor(dbf_df$TF_name, levels=sorted_avg_df$TF_name)
#cr_df$TF_name <- factor(cr_df$TF_name, levels=sorted_avg_df$TF_name)

output_file_name <- paste0("~/Dropbox", "/", "motif_based_dbf_meth_group_barplot_for_208_TFs.pdf")             
pdf(output_file_name)

# z <- c(0.15, 0.15)
#test <- data.frame(X = c("DBF", "CR/CF" ), Z = c(0.15, 0.15))

barplot <- ggplot(meth_tf_df, aes(x=tf_name, y=meth_group_percent, fill=meth_group)) + 
    geom_bar(stat="identity") +
    #geom_jitter(size=0.01, color="black", alpha=0.01) +
    #geom_point(colour = "blue", position = "jitter", aeslpha=0.5) +
    #stat_boxplot( aes(TF_name, meth_percent), 
    #geom='errorbar', linetype=1, width=0.5) +  #whiskers
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25, color="black") + 
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    #stat_summary(fun.y=mean, geom="point", size=2) + 
    #stat_summary(fun.data = mean_se, geom = "errorbar") +

    xlab("Transcription Factors") + ylab("Fraction of motifs in different methylation group bins") + 
    theme_bw() + 
    ggtitle("Fraction of motifs in different methylation bins") + 
    theme(
    axis.text.y = element_text(face="bold", size=2.2),
    plot.title=element_text(size=14, face="bold", hjust = 0.6),
    legend.title=element_text(face="bold")
    ) + 
    coord_flip() + 
    guides(fill=guide_legend(title="Category")) +
    scale_fill_manual(values=c("#F8766D", "#619CFF", "#00BA38"))+
    theme(axis.text=element_text(size=5),
        #axis.title=element_text(size=10,face="bold")) + coord_flip() + scale_fill_gradient(low="light green", high="dark red")
        axis.title=element_text(size=10,face="bold"))
    #scale_y_continuous(limits = quantile(meth_tf_df$meth_percent, c(0.0, 1.0))) +
    #facet_wrap(~tf_category)

print(barplot)

dev.off()


avg_df <- fread("~/Dropbox/encode_3/fimo_motif_based_analysis/wgbs_tf_intersect_avg_of_motifs/files/final_wgbs_tf_intersect_average_meth_unique_TFs.txt", sep="\t", header=TRUE)
meth_tf_df <- fread("~/Dropbox/encode_3/fimo_motif_based_analysis/ideas_cpg_motif_intersect_total/files/final_wgbs_tf_motif_ideas_intersect_facetwrap_barplot_data.bed", sep="\t", header=TRUE)
#names(meth_tf_df) <- c("TF_name", "index", "peak_chrom", "peak_start", "peak_end", "meth_percent")

target_vector <- meth_tf_df$tf_name

## Now, sort the order of TF wrt methylation value:
sorted_avg_df <- avg_df[order(avg_df$meth_percent),]
factor(sorted_avg_df$TF_name) %>% levels
sorted_avg_df$TF_name <- factor(sorted_avg_df$TF_name, levels=sorted_avg_df$TF_name)

meth_tf_df$tf_name <- factor(meth_tf_df$tf_name, levels=sorted_avg_df$TF_name)
all_tf_df <- meth_tf_df %>% as.data.frame

# change the order of factor in TF Category:
#meth_tf_df$tf_category <- factor(meth_tf_df$tf_category, levels = c("DBF", "CR/CF" ))
#meth_tf_df$TF_name <- factor(meth_tf_df$TF_name, levels=sorted_avg_df$TF_name)

#dbf_df$TF_name <- factor(dbf_df$TF_name, levels=sorted_avg_df$TF_name)
#cr_df$TF_name <- factor(cr_df$TF_name, levels=sorted_avg_df$TF_name)

output_file_name <- paste0("~/Dropbox", "/", "motif_based_dbf_meth_group_ideas_barplot_for_208_TFs.pdf")             
pdf(output_file_name)

barplot <- ggplot(subset(meth_tf_df, meth_group != "Intermediate"), aes(x=tf_name, y=meth_group_percent, fill=ideas_state_map)) + 
    geom_bar(stat="identity") +
    #geom_jitter(size=0.01, color="black", alpha=0.01) +
    #geom_point(colour = "blue", position = "jitter", aeslpha=0.5) +
    #stat_boxplot( aes(TF_name, meth_percent), 
    #geom='errorbar', linetype=1, width=0.5) +  #whiskers
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25, color="black") + 
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    #stat_summary(fun.y=mean, geom="point", size=2) + 
    #stat_summary(fun.data = mean_se, geom = "errorbar") +

    xlab("Transcription Factors") + ylab("Fraction of motifs in different methylation group bins") + 
    theme_bw() + 
    ggtitle("Fraction of motifs in different methylation bins") + 
    theme(
    axis.text.y = element_text(face="bold", size=2.2),
    plot.title=element_text(size=14, face="bold", hjust = 0.6),
    legend.title=element_text(face="bold")
    ) + 
    coord_flip() + 
    guides(fill=guide_legend(title="Category")) +
    #scale_fill_manual(values=c("#F8766D", "#619CFF", "#00BA38"))+
    theme(axis.text=element_text(size=5),
        #axis.title=element_text(size=10,face="bold")) + coord_flip() + scale_fill_gradient(low="light green", high="dark red")
        axis.title=element_text(size=10,face="bold")) +
    #scale_y_continuous(limits = quantile(meth_tf_df$meth_percent, c(0.0, 1.0))) +
    facet_wrap(~meth_group)

print(barplot)

dev.off()





