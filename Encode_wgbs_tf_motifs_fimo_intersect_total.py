#!/usr/bin/env python

import pandas as pd
import pybedtools
import os
import re
from glob import glob
from os.path import basename
from os.path import join
from os.path import splitext


def meth_percent_calc(x):
    x_meth = x["cpg_meth"].sum()
    x_unmeth = x["cpg_meth"].sum() + x["cpg_unmeth"].sum()
    x_perc = x_meth/float(x_unmeth)
    return(x_perc)

### Input and output file and dir paths:
output_dir = os.path.expanduser("~/for_chris/batch_I/fimo_motifs_total/motif_meth_analysis/cpg_motif_intersect_total")
cpg_bed_file = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/SL83768_SL83771_300x_merged_CpG_context_deduplicated.bismark.cov"
#peak_bed_list = glob("/gpfs/gpfs1/home/schhetri/for_chris/batch_I/fimo_motifs_total/idr_passed_motifs_total/unique_TFs/SL*fimo_motif*")
peak_bed_list = glob("/gpfs/gpfs1/home/schhetri/for_chris/batch_I/fimo_motifs_total/idr_passed_motifs_total/SL*fimo_motif*")

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

""" \nProcessing and creating a Cpg bed file with a strand cols from raw cpg file 
followed by pybedtool object formation...\n """
print " \nProcessing Cpg bed file\n "
os.environ["output_dir"] = output_dir
os.environ["cpg_bed_file"] = cpg_bed_file
cpg_bed_edited = splitext(basename(cpg_bed_file))[0] + "_edited" + splitext(cpg_bed_file)[1]
print cpg_bed_edited
os.environ["cpg_bed_edited"] = cpg_bed_edited

if not os.path.exists(join(output_dir,cpg_bed_edited)):
    # os.system(CMD) #subprocess.call(COMMAND, shell=True) # import subprocess
    os.system("ln -fs $cpg_bed_file $output_dir/$cpg_bed_edited")
else:
    print "CpG file already converted to bed format..."

print "Generating pybedtools object..."
Cpg_bed_file = pybedtools.BedTool(join(output_dir,cpg_bed_edited))


master_tf_dict = {}
for each_file in peak_bed_list:
    print "processing", each_file
    """ \nProcessing and creating a sorted bed file for the cpg intersect
    followed by pybedtools object formation...\n """
    os.environ["each_file"] = each_file

    selected_motif_file = splitext(basename(each_file))[0] + (".motif_1_2" + ".txt")
    os.environ["selected_motif_file"] = selected_motif_file

    if not os.path.exists(join(output_dir,selected_motif_file)):
        #os.system('''egrep "^1\\b|^2\\b" $each_file | cut -f2-4 | awk 'BEGIN{OFS="\t"} {chr=$1;start=$2;end=$3; print chr, start+1+(-2), end+1+(2+1)}' > $output_dir/$selected_motif_file''')
        os.system('''egrep "^1\\b|^2\\b" $each_file | cut -f2-4 | awk 'BEGIN{OFS="\t"} {chr=$1;start=$2;end=$3; print chr, start+1, end+1+(1)}' > $output_dir/$selected_motif_file''')
    else:
        print "%s : Motif file exists"%(basename(each_file))

    sorted_peak_file = splitext(basename(each_file))[0] + "_sorted" + splitext(basename(each_file))[1]
    os.environ["sorted_peak_file"] = sorted_peak_file
    os.system('sort -k1,1 -k2,2n $output_dir/$selected_motif_file > $output_dir/$sorted_peak_file')
    TF_name_list = re.findall(r'.*fimo_motif_(.*).txt$', basename(each_file))
    TF_name =  TF_name_list[0]
    peak_bed_file = pd.read_csv(join(output_dir,sorted_peak_file), sep="\t", header=None)
    peak_count = peak_bed_file.drop_duplicates().shape[0]
    peak_bed = pybedtools.BedTool(join(output_dir,sorted_peak_file))

    """ \nPybedtool intersection of Cpg_bed_file and peak_bed_file...\n """
    print " Processing the intersection for", TF_name
    pybed_outfile = join(output_dir, (TF_name + "_sorted_cpg_intersect.bed"))
    pybed_outfile_v = join(output_dir, (TF_name + "_sorted_cpg_not_intersect.bed"))
    
    if not os.path.exists(pybed_outfile):
        Cpg_bed_file.intersect(peak_bed, wa = True, wb = True, output=pybed_outfile)
        # Cpg_bed_file.intersect(peak_bed, wa = True, wb = True, v = True, output=pybed_outfile_v)
    else:
        print "Pybedtools object already present"

    if os.stat(pybed_outfile).st_size == 0:
        if TF_name not in master_tf_dict:
            master_tf_dict[TF_name] = None
    else:
        ### Working with the dataframes; reading the output file of pybedtool intersect:
        df_file = pd.read_csv(pybed_outfile, sep = "\t", header = None)
        df_ordered = df_file.iloc[:, [6, 7, 8, 0, 1, 2, 3, 4]]
        df_ordered.columns = ["peak_chrom", "peak_start", "peak_end", "cpg_chrom", "cpg_start", "cpg_end", "cpg_meth", "cpg_unmeth" ]
        df_ordered = df_ordered.drop_duplicates(subset=["peak_chrom", "peak_start", "peak_end"])
        df_ordered[["cpg_meth", "cpg_unmeth"]] = df_ordered[["cpg_meth", "cpg_unmeth"]].astype(int)
        df_grouped =  df_ordered.groupby(["peak_chrom", "peak_start", "peak_end"]).apply(lambda x : x["cpg_meth"].sum()/float(x["cpg_meth"].sum() + x["cpg_unmeth"].sum()))
        #peak_grouped = df_ordered.groupby(["peak_chrom", "peak_start", "peak_end"]).apply(lambda x: x.sum())  ### For quick test of the grouped items:
        #peak_grouped = df_ordered.groupby(["peak_chrom", "peak_start", "peak_end"]).apply(meth_percent_calc)  
        print "Dimension of current TF is", df_grouped.reset_index().shape
        if TF_name not in master_tf_dict:
            master_tf_dict[TF_name] = df_grouped.reset_index()

        print "\nIntersection of %s completed!!!...\n\n" %(TF_name)

### Combine all the dataframes (i.e all values of master_tf_dict) representing the tf intersects with CpG:
tf_combined_df = pd.concat(master_tf_dict).reset_index()
#tf_combined_df.to_csv(join(output_dir, "final_wgbs_tf_motif_intersect_unique_TFs.txt"), sep ="\t", header = True, index = False)
tf_combined_df.to_csv(join(output_dir, "final_wgbs_tf_motif_intersect_total_TFs.txt"), sep ="\t", header = True, index = False)

#tf_combined_grouped  = tf_combined_df.groupby("level_0").apply(lambda x : x["0"].sum()/float(len(x["0"])))
print "Job completed successfully!!!"
