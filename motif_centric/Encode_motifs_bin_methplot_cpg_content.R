library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(gtools)
library(gplots)
library(RColorBrewer)


### For motif based methylation plots:
tf_name = "GATA"
input_file <- "~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/motifs_methplot/final_wgbs_motifs_intersect_2kb_GATA.bed"
binned_perc_meth_table <- fread(input_file, sep="\t", header= TRUE)

names(binned_perc_meth_table) <- c("index", "annotation", "bin_start", "bin_end", "meth_percent")
binned_perc_meth_table$bin_mid <- with(binned_perc_meth_table, (bin_start + bin_end)/2)


output_file_name <- paste0("~/Dropbox", "/", "motif_composite_plot.pdf")       
pdf(output_file_name, width=4.5, height=5)

this_plot <- ggplot2::ggplot(binned_perc_meth_table, aes(bin_mid, meth_percent)) +
ggplot2::geom_line(aes(color=annotation)) +
ggplot2::ylab("Percent Methylation") +
ggplot2::xlab("Distance from Center (bp)") +
ggplot2::ggtitle(paste(tf_name,"motifs methylation profile"))+
#ggplot2::ggsave(output_file_name, width = 4, height = 10, dpi = 200)+
ggplot2::scale_x_continuous() +
ggplot2::scale_y_continuous(limits=c(0,1)) +
ggplot2::theme( legend.key = element_blank(),
      #legend.title=element_text(size=8, face="bold"),
      #legend.text=element_text(size=5),
      legend.key.size = unit(0.2, "cm"),
      axis.line = element_line(),
      panel.background = element_blank(),
      axis.text.x = element_text(color="black", size=12),
      axis.title.x = element_text(color="black", size=14, vjust=-1.5),
      axis.text.y = element_text(color="black", size=12),
      axis.title.y = element_text(color="black", size=14, vjust=3),
      plot.title=element_text(size=14, face="bold"),
      plot.margin = grid::unit(c(1,1,1,1), "cm"),
      panel.border=element_blank(),
      axis.ticks=element_line(size=0.6, color="black")) + theme_bw()

this_plot

dev.off()



### For heatmaps associated to motifs:
# CG, CHG, CHH containing motif Distribution:
#read_motif_file = fread("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/motif_analysis/final_bed_files/final_cpg_chg_chh_motif_percent.bed")
read_motif_file = fread("~/Dropbox/encode_3/fimo_motif_based_analysis/cpg_chg_chh_motif_percent_content/files/final_cpg_chg_chh_motif_percent.txt")
read_df = as.data.frame(read_motif_file)
read_df = read_df[,c(1,2,3,4)]

rnames <- read_df[,1]
data <- read_df[,2:ncol(read_df)]
mat_data <- as.matrix(data)
rownames(mat_data) <- rnames
rownames(mat_data)
colnames(mat_data)


output_file_name <- paste0("~/Dropbox", "/", "final_cg_cpg_chg_motif_distribution_heatmap.pdf")				
pdf(output_file_name)

heatmap.2(mat_data,
  main = "Motif Sites Distribution", 
  xlab = "",
  ylab = "208 Transcription Factors",
  col= brewer.pal(9,"Blues"), #greenred(10), 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(8,10),      # widens margins around plot)
  cexRow = 0.15,
  cexCol = 0.6,
  lhei = c(2, 9), #sepwidth=c(2,2),
  Colv = FALSE,
  na.color="grey"
  #scale="column", tracecol="#303030"
  )    

dev.off()


### Not necessary, same as above just non-C added up i.e CG, CHG, CHH, non C motif Distribution:
#read_motif_file = fread("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/motif_analysis/final_bed_files/final_cpg_chg_chh_motif_percent.bed")
read_motif_file = fread("~/Dropbox/encode_3/fimo_motif_based_analysis/cpg_chg_chh_motif_percent_content/files/final_cpg_chg_chh_motif_percent.txt")
read_df = as.data.frame(read_motif_file)
read_df = read_df[,c(1,2,3,4,6)]

rnames <- read_df[,1]
data <- read_df[,2:ncol(read_df)]
mat_data <- as.matrix(data)
rownames(mat_data) <- rnames
rownames(mat_data)
colnames(mat_data)


output_file_name <- paste0("~/Dropbox", "/", "final_cg_cpg_chg_nonC_motif_distribution_heatmap.pdf")				
pdf(output_file_name)

heatmap.2(mat_data,
  main = "Motif Sites Distribution", 
  xlab = "",
  ylab = "208 Transcription Factors",
  col= brewer.pal(9,"Blues"), #greenred(10), 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(8,10),      # widens margins around plot)
  cexRow = 0.15,
  cexCol = 0.6,
  lhei = c(2, 9), #sepwidth=c(2,2),
  #Colv = FALSE,
  na.color="grey"
  #scale="column", tracecol="#303030"
  )
dev.off()



### Barplot Distribution of the CG containing motifs
read_motif_file = fread("~/Dropbox/encode_3/fimo_motif_based_analysis/cpg_chg_chh_motif_percent_content/files/final_cpg_chg_chh_motif_percent.txt")
#read_motif_file = fread("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/motif_analysis/final_bed_files/final_cpg_chg_chh_motif_percent.txt")

read_df = as.data.frame(read_motif_file)
#read_df$"" <- 0
read_df = read_df[,c(1,2)]
names(read_df) <- c("tf_name","CG_motif_percent")
read_df <- read_df[order(read_df$CG_motif_percent),]

factor_order <- read_df$tf_name
read_df$tf_name <- factor(read_df$tf_name, levels=factor_order)
factor(read_df$tf_name)


output_file_name <- paste0("~/Dropbox", "/", "final_CG_sites_distribution_barplot_1.pdf")				
pdf(output_file_name)

ggplot(read_df, aes(x=tf_name, y=CG_motif_percent, fill=CG_motif_percent)) + geom_bar(stat="identity") + theme_bw() +
      xlab("Transcription Factors") + ylab("CG motifs percent") + 
      ggtitle("CG containing motifs distribution") + 
      guides(fill=guide_legend(title="CG Percent(%)")) +
      theme(
      axis.text.y = element_text(size=2.2, face="bold" ),
      plot.title=element_text(size=12, face="bold", hjust = 0.6),
      legend.title=element_text(face="bold")
      ) +
      theme(axis.text=element_text(size=5),
        #axis.title=element_text(size=10,face="bold")) + coord_flip() + scale_fill_gradient(low="light green", high="dark red")
        axis.title=element_text(size=10,face="bold")) + coord_flip() + scale_fill_gradient(low="light green", high="dark red")

dev.off()


# ggplot(read_df, aes(x=factor(tf_name, levels = read_df$tf_name ), y=CG_motif_percent, fill=CG_motif_percent)) + geom_bar(stat="identity") + 
# theme(axis.text=element_text(size=5),
#         axis.title=element_text(size=10,face="bold")) + coord_flip() 


# bar <- ggplot(data = df, aes(x = factor(dat, levels = month.abb), y = val, 
#               fill=val)) +
#        geom_bar(stat = 'identity') + 
#        scale_fill_gradient2(low=LtoM(100), mid='snow3', 
#        high=MtoH(100), space='Lab')
