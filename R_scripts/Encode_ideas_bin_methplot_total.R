library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
library(gtools)

args <-  commandArgs(TRUE)
ideas_cpg_tf_intersect_file <- args[1]
tf_name <- args[2]
out_dir_name <- args[3]
bin_size <-  args[4]
bin_size <- as.numeric(50)

print(paste(ideas_cpg_tf_intersect_file, "filename...."))
print(paste(tf_name, "tf name...."))
print(paste(out_dir_name, "output dir name...."))

### For ideas distinct regulatory elements plot:
binned_perc_meth_table <- fread(ideas_cpg_tf_intersect_file, sep="\t", header= TRUE)

names(binned_perc_meth_table) <- c("index", "annotation", "bin_start", "bin_end", "state", "meth_percent")
binned_perc_meth_table$bin_mid <- with(binned_perc_meth_table, (bin_start + bin_end)/2)

### Maintain the order of the list as per the targetl list:
ordered_df <- binned_perc_meth_table[mixedorder(binned_perc_meth_table$state),]
binned_perc_meth_table$state <- factor(binned_perc_meth_table$state, levels= ordered_df$state)
factor(binned_perc_meth_table$state)

#bin_size <-  100
min_size <- min(binned_perc_meth_table$bin_mid) #original upstream
max_size <- max(binned_perc_meth_table$bin_mid) #original downstream

label_min_size <- min(binned_perc_meth_table$bin_mid) + (-bin_size/2) #Just to label upstream (but not original upstream)
label_max_size <- max(binned_perc_meth_table$bin_mid) + (bin_size/2) #Just to label downstream (but not original downstream)

# second_max_size <- sort(unique(binned_perc_meth_table$bin_mid))[2]
# second_min_size <- sort(unique(binned_perc_meth_table$bin_mid), decreasing=TRUE)[2]

output_file_name <- paste0(out_dir_name, "/", tf_name, "_ideas_composite_plot.pdf")
pdf(output_file_name)

this_plot <-ggplot2::ggplot(binned_perc_meth_table, aes(bin_mid, meth_percent)) +
            ggplot2::geom_line(aes(color=annotation)) +
            ggplot2::ylab("Percent Methylation") +
            ggplot2::xlab("Distance from Center (bp)") +
            ggplot2::scale_x_continuous(breaks=c(min_size,0,max_size), labels=c(label_min_size/1000, "0", label_max_size/1000)) +
            #ggplot2::scale_x_continuous(breaks=c(second_min_size,0,second_max_size), labels=c(paste0(second_min_size/1000), "0", paste0(second_max_size/1000))) +
            ggplot2::scale_y_continuous(limits=c(0,1)) +
            ggplot2::theme( legend.key = element_blank(),
                  axis.line = element_line(),
                  panel.background = element_blank(),
                  axis.text.x = element_text(color="black", size=12),
                  axis.title.x = element_text(color="black", size=14, vjust=-1.5),
                  axis.text.y = element_text(color="black", size=12),
                  axis.title.y = element_text(color="black", size=14, vjust=3),
                  plot.margin = grid::unit(c(1,1,1,1), "cm"),
                  panel.border=element_blank(),
                  axis.ticks=element_line(size=0.6, color="black")) + theme_bw() +
            ggplot2::facet_wrap(~state) 

this_plot

dev.off()




