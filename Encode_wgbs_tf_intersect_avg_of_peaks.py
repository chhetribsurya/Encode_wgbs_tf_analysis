#!/usr/bin/env python

import pandas as pd
import pybedtools
import glob
import os
import re
from os.path import basename
from os.path import join
from os.path import splitext


def meth_percent_calc(x):
    x_meth = x["cpg_meth"].sum()
    x_unmeth = x["cpg_meth"].sum() + x["cpg_unmeth"].sum()
    x_perc = x_meth/float(x_unmeth)
    return(x_perc)

### Input and output file and dir paths:
output_dir = os.path.expanduser("~/for_chris/batch_I/wgbs_tf_intersect_total")
cpg_bed_file = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/SL83768_SL83771_300x_merged_CpG_context_deduplicated.bismark.cov"
#cpg_bed_file = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/SL83768_CpG_context_deduplicated.bismark.cov_head"
#peak_bed_list = glob.glob("/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/SL*narrowPeak*")
peak_bed_list = glob.glob("/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL*narrowPeak*")

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
    #CMD = '''awk 'BEGIN{FS=" "; OFS="\t"} { print $1,$2,$3+1,$5,$6,"." }' $cpg_bed_file | sort -k1,1 -k2,2n | grep -v "^*" > $output_dir/$cpg_bed_edited'''
    #os.system(CMD) #subprocess.call(COMMAND, shell=True) # import subprocess
    os.system("ln -fs $cpg_bed_file $output_dir/$cpg_bed_edited")
else:
    print "CpG file already converted to bed format..."

print "Generating pybedtools object..."
Cpg_bed_file = pybedtools.BedTool(join(output_dir,cpg_bed_edited))



header = ["TF_name", "meth_percent" ]
master_tf_values = [[],[]]

for each_file in peak_bed_list:
    print "processing", each_file
    """ \nProcessing and creating a sorted bed file for the cpg intersect
    followed by pybedtools object formation...\n """
    os.environ["each_file"] = each_file
    sorted_peak_file = splitext(basename(each_file))[0] + "_sorted" + splitext(basename(each_file))[1]
    os.environ["sorted_peak_file"] = sorted_peak_file
    os.system('sort -k1,1 -k2,2n $each_file > $output_dir/$sorted_peak_file')
    #TF_name = basename(each_file).split("_")[-1] 
    TF_name_list = re.findall(r'.*narrowPeak_(.*)$', basename(each_file))
    TF_name =  TF_name_list[0]
    peak_bed = pybedtools.BedTool(join(output_dir,sorted_peak_file))
    #peak_bed = pybedtools.BedTool(each_file)

    """ \nPybedtool intersection of Cpg_bed_file and peak_bed_file...\n """
    print " Processing the intersection for", TF_name
    pybed_outfile = join(output_dir, (TF_name + "_sorted_cpg_intersect.bed"))
    pybed_outfile_v = join(output_dir, (TF_name + "_sorted_cpg_not_intersect.bed"))

    if not os.path.exists(pybed_outfile):
        Cpg_bed_file.intersect(peak_bed, wa = True, wb = True, output=pybed_outfile)
		#Cpg_bed_file.intersect(peak_bed, wa = True, wb = True, v = True, output=pybed_outfile_v)
    else:
        print "Pybedtools object already present"


    ### Working with the dataframes; reading the output file of pybedtool intersect:
    df_file = pd.read_csv(pybed_outfile, sep = "\t", header = None)
    df_ordered = df_file.iloc[:, [6, 7, 8, 0, 1, 2, 3, 4]]
    df_ordered.columns = ["peak_chrom", "peak_start", "peak_end", "cpg_chrom", "cpg_start", "cpg_end", "cpg_meth", "cpg_unmeth" ]
    df_ordered[["cpg_meth", "cpg_unmeth"]] = df_ordered[["cpg_meth", "cpg_unmeth"]].astype(int)
    avg_meth_percent =  df_ordered["cpg_meth"].sum()/float(df_ordered["cpg_meth"].sum() + df_ordered["cpg_unmeth"].sum())
    master_tf_values[0].append(TF_name)
    master_tf_values[1].append(avg_meth_percent)
    #peak_grouped = df_ordered.groupby(["peak_chrom", "peak_start", "peak_end"]).apply(lambda x: x.sum())  ### For quick test of the grouped items:
    #peak_grouped = df_ordered.groupby(["peak_chrom", "peak_start", "peak_end"]).apply(meth_percent_calc)  
    print "\nIntersection of %s completed!!!...\n\n" %(TF_name)


### Generate master_tf_dict, and combine all the values of master_tf_dict to form dataframe:
master_tf_dict = dict(zip(header, master_tf_values))
master_tf_df = pd.DataFrame(master_tf_dict)
master_tf_df.to_csv(join(output_dir, "final_wgbs_tf_intersect_average_meth.bed"), sep ="\t", header = True, index = False)
print master_tf_df

#df_grouped = df.groupby("level_0").apply(lambda x : x["0"].sum()/float(len(x["0"])))
print "Job completed successfully!!!"

