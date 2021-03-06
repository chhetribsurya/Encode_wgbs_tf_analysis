#!/usr/bin/env python
import pandas as pd
import pybedtools
import glob
import os
import re
import pickle
from os.path import basename
from os.path import join
from os.path import splitext


def meth_percent_calc(x):
    x_meth = x["cpg_meth"].sum()
    x_unmeth = x["cpg_meth"].sum() + x["cpg_unmeth"].sum()
    x_perc = x_meth/float(x_unmeth)
    x_cpg_count = x["cpg_count"].sum()
    return(x_perc)

### Input the hepg2 ideas segmenation:
ideas_hepg2_file = os.path.expanduser("~/for_chris/batch_I/hepg2_ideas_36_dense.bed")  

### Input and output file and dir paths:
output_dir = os.path.expanduser("~/for_chris/batch_I/wgbs_tf_intersect_total")
cpg_bed_file = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/SL83768_SL83771_300x_merged_CpG_context_deduplicated.bismark.cov"
# cpg_bed_file = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/SL83768_CpG_context_deduplicated.bismark.cov_head"
peak_bed_list = glob.glob("/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL*narrowPeak*")

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

os.environ["output_dir"] = output_dir

wgbs_tf_file_for_ideas = join(output_dir, "final_wgbs_tf_intersect_preping_for_ideas.bed")

if not os.path.exists(wgbs_tf_file_for_ideas):
    """ \nProcessing and creating a Cpg bed file with a strand cols from raw cpg file 
    followed by pybedtool object formation...\n """
    print " \nProcessing Cpg bed file\n "
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

    print "Generating bedtools object..."
    Cpg_bed_file = pybedtools.BedTool(join(output_dir,cpg_bed_edited))



    master_tf_dict = {}

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

        """ \nPybedtool intersection of Cpg_bed_file and peak_bed_file...\n """
        print " Processing the intersection for", TF_name
        pybed_outfile = join(output_dir, (TF_name + "_sorted_cpg_intersect.bed"))
        #pybed_outfile_v = join(output_dir, (TF_name + "_cpg_not_intersect.bed"))

        if not os.path.exists(pybed_outfile):
            Cpg_bed_file.intersect(peak_bed, wa = True, wb = True, output=pybed_outfile)
        #   Cpg_bed_file.intersect(peak_bed, wa = True, wb = True, v = True, output=pybed_outfile_v)
        else:
            print "Pybedtools object already present"

        ### Working with the dataframes; reading the output file of pybedtool intersect:
        df_file = pd.read_csv(pybed_outfile, sep = "\t", header = None)
        df_ordered = df_file.iloc[:, [0, 1, 2, 3, 4, 6, 7, 8]]
        df_ordered.columns = ["cpg_chrom", "cpg_start", "cpg_end", "cpg_meth", "cpg_unmeth", "peak_chrom", "peak_start", "peak_end" ]
        df_ordered[["cpg_meth", "cpg_unmeth"]] = df_ordered[["cpg_meth", "cpg_unmeth"]].astype(int)
        df_ordered["strand"] = "."
        ### Making sure its in proper bed format with strand in 6th column:
        df_ordered = df_ordered.loc[:,["cpg_chrom", "cpg_start", "cpg_end", "cpg_meth", "cpg_unmeth", "strand", "peak_chrom", "peak_start", "peak_end"]]
        print "Dimension of current TF is", df_ordered.shape

        if TF_name not in master_tf_dict:
          master_tf_dict[TF_name] = df_ordered    

    with open(join(output_dir,"final_wgbs_tf_intersect_preping_for_ideas.pkl"), "w") as outfile:
        pickle.dump(master_tf_dict, outfile)

    with open(join(output_dir,"final_wgbs_tf_intersect_preping_for_ideas.pkl")) as infile:
        master_tf_dict = pickle.load(infile)

    ### Combine all the dataframes (i.e all values of master_tf_dict) representing the tf intersects with CpG:
    tf_combined_df = pd.concat(master_tf_dict).reset_index()
    tf_combined_df.to_csv(join(output_dir, "final_wgbs_tf_intersect_preping_for_ideas_1.bed"), sep ="\t", header = False, index = False)


#########################################
#                                       #
# Main starting point for the           #
# instersection b/w final wgbs_tf       #
# intersect file and ideas hepg2        #
# file                                  #
#                                       #
#########################################


### Intersection with ideas file, after sorting the wgbs_tf_intersect bed:
print " \nProcessing final wgbs_tf_intersect file for strand and cols rearrangement\n "
wgbs_tf_to_edit_file = join(output_dir, "final_wgbs_tf_intersect_preping_for_ideas_1.bed")
os.environ["final_wgbs_tf_intersect_file"] = wgbs_tf_to_edit_file
wgbs_tf_bed_edited = splitext(basename(wgbs_tf_to_edit_file))[0] + "_edited" + splitext(wgbs_tf_to_edit_file)[1]
print wgbs_tf_bed_edited
os.environ["wgbs_tf_bed_edited"] = wgbs_tf_bed_edited
os.environ["output_dir"] = output_dir

if not os.path.exists(join(output_dir,wgbs_tf_bed_edited)):
    CMD = '''awk 'BEGIN{FS= " "; OFS="\t"} { print $3,$4,$5,$6,$7,$8,$9,$10,$11,$1 }' $final_wgbs_tf_intersect_file | sort -k1,1 -k2,2n  > $output_dir/$wgbs_tf_bed_edited'''
    os.system(CMD) #subprocess.call(COMMAND, shell=True) # import subprocess
else:
    print "\nwgbs_tf_bed_edited file already present\n"

#wgbs_tf_intersect_bed = join(output_dir, "final_wgbs_tf_intersect_preping_for_ideas.bed")
#sorted_wgbs_tf_file = splitext(basename(wgbs_tf_intersect_bed))[0] + "_sorted" + splitext(basename(wgbs_tf_intersect_bed))[1]
#os.environ["wgbs_tf_file"] = wgbs_tf_intersect_bed
#os.environ["sorted_wgbs_tf_file"] = sorted_wgbs_tf_file
#print "Sorting of the wgbs_tf_intersect file..."
#if not os.path.exists(join(output_dir,sorted_wgbs_tf_file)):
#    os.system('sort -k1,1 -k2,2n $wgbs_tf_file > $output_dir/$sorted_wgbs_tf_file')
#else:
#    print "wgbs_tf sorted file already exists"


sorted_ideas_file = splitext(basename(ideas_hepg2_file))[0] + "_sorted" + splitext(basename(ideas_hepg2_file))[1]
os.environ["ideas_hepg2_file"] = ideas_hepg2_file
os.environ["sorted_ideas_file"] = sorted_ideas_file
print "\nSorting of the hepg2 ideas file...\n"
if not os.path.exists(join(output_dir,sorted_ideas_file)):
    os.system('sort -k1,1 -k2,2n $ideas_hepg2_file > $output_dir/$sorted_ideas_file')
else:
    print "\nHepG2 sorted file already exists\n"


print "\nStarting the intersection b/w wgbs_tf_file and ideas segmentation...\n"
wgbs_tf_bed = pybedtools.BedTool(join(output_dir,wgbs_tf_bed_edited))
ideas_hepg2_bed =  pybedtools.BedTool(join(output_dir,sorted_ideas_file))
pybed_ideas_outfile = join(output_dir, "final_wgbs_tf_ideas_intersect.bed") 
if not os.path.exists(pybed_ideas_outfile):
        wgbs_tf_bed.intersect(ideas_hepg2_bed, wa = True, wb = True, output=pybed_ideas_outfile)
else:
   print "\nPybedtools object present for wgbs_tf ideas intersection\n"

print "\n Processing the final_wgbs_tf_ideas intersection file for heatmap data generation...\n"
if not os.path.exists(join(output_dir, "final_wgbs_tf_ideas_intersect_heatmap_data.bed")):
    read_df = pd.read_csv(join(output_dir, "final_wgbs_tf_ideas_intersect.bed"), sep="\t", header=None)
    df_select = read_df.iloc[:,[0,1,2,3,4,9,13]]
    df_select.columns = ["cpg_chrom", "cpg_start", "cpg_end", "cpg_meth", "cpg_unmeth","tf_name", "ideas_state"]
    
    # To find unique cpg binding to each TF and state:
    #df_unique_cpg_df = df_select.drop_duplicates(["cpg_chrom", "cpg_start", "cpg_end"])
    #df_unique_cpg_df["cpg_count"] = 1
    #df_grouped_cpg_df = df_unique_cpg_df.groupby(["tf_name", "ideas_state"]).apply(lambda x : (x["cpg_meth"].sum()/float(x["cpg_meth"].sum() + x["cpg_unmeth"].sum()), x["cpg_count"]))
    #df_grouped_cpg_df_reset = df_grouped_cpg_df.reset_index()
    #df_grouped_cpg_df1_reset[["meth_percent","cpg_count"]] = df_grouped_cpg_df_reset.loc[:,0].apply(pd.Series) 
    #df_select["cpg_count"] = .groupby(["tf_name", "ideas_state"]).apply(lambda x : (x["cpg_meth"].sum()/float(x["cpg_meth"].sum() + x["cpg_unmeth"].sum()), x["cpg_count"]))
    
    print "\nGrouping the data for heatmap...\n"
    df_grouped = df_select.groupby(["tf_name", "ideas_state"]).apply(lambda x : x["cpg_meth"].sum()/float(x["cpg_meth"].sum() + x["cpg_unmeth"].sum()))
    print "Dimension after grouping", df_grouped.shape
    final_df_unstacked = df_grouped.unstack()
    print "Dimension after stacking", final_df_unstacked.shape
    final_df_unstacked.to_csv(join(output_dir, "final_wgbs_tf_ideas_intersect_heatmap_data.bed"), sep = "\t", header = True, index = True )

#tf_combined_grouped  = tf_combined_df.groupby("level_0").apply(lambda x : x["0"].sum()/float(len(x["0"])))
print "Job completed successfully!!!"
