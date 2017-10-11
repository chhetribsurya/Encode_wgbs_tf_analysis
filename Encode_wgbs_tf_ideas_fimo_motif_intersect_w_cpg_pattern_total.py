#!/usr/bin/env python

import pandas as pd
import pybedtools
import os
import re
import errno
from glob import glob
from os.path import basename
from os.path import join
from os.path import splitext
import pickle


def meth_percent_calc(x):
    x_meth = x["cpg_meth"].sum()
    x_unmeth = x["cpg_meth"].sum() + x["cpg_unmeth"].sum()
    x_perc = x_meth/float(x_unmeth)
    return(x_perc)


### Input the hepg2 ideas segmenation:
ideas_hepg2_file = os.path.expanduser("~/for_chris/batch_I/hepg2_ideas_36_dense.bed")  

### Input and output file and dir paths:
output_dir = os.path.expanduser("~/for_chris/batch_I/fimo_motifs_total/motif_meth_analysis/ideas_cpg_motif_intersect_total_corrected")
cpg_bed_file = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/SL83768_SL83771_300x_merged_CpG_context_deduplicated.bismark.cov"
# fimo_list = glob("/gpfs/gpfs1/home/schhetri/for_chris/batch_I/fimo_motifs_total/idr_passed_motifs_total/SL*fimo_motif*")
fimo_list = glob("/gpfs/gpfs1/home/schhetri/for_chris/batch_I/fimo_motifs_total/idr_passed_motifs_total/unique_TFs/SL*fimo_motif*")
# if not os.path.exists(output_dir):
#     os.makedirs(output_dir)

try:
    os.makedirs(output_dir)
except OSError as exc:
    if exc.errno != errno.EEXIST:
        raise
    pass


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

motif_count_dict = {"tf_name":[], "motif_count_nodups":[]}
master_tf_dict = {}
for each_file in fimo_list:
    print "processing", each_file
    """ \nProcessing and creating a sorted bed file for the cpg intersect
    followed by pybedtools object formation...\n """
    os.environ["each_file"] = each_file

    selected_motif_file = splitext(basename(each_file))[0] + (".motif_1_2" + ".txt")
    os.environ["selected_motif_file"] = selected_motif_file

    if not os.path.exists(join(output_dir,selected_motif_file)):
        os.system('''egrep "^1\\b|^2\\b" $each_file | cut -f2-4 | awk 'BEGIN{OFS="\t"} {chr=$1;start=$2;end=$3; print chr, start+1+(-2), end+1+(2+1)}' > $output_dir/$selected_motif_file''')
        # os.system('''egrep "^1\\b|^2\\b" $each_file | cut -f2-4 | awk 'BEGIN{OFS="\t"} {chr=$1;start=$2;end=$3; print chr, start+1, end+1+(1)}' > $output_dir/$selected_motif_file''')
    else:
       print "%s : Motif file exists"%(basename(each_file))

    sorted_peak_file = splitext(basename(each_file))[0] + "_sorted" + splitext(basename(each_file))[1]
    os.environ["sorted_peak_file"] = sorted_peak_file
    os.system('sort -k1,1 -k2,2n $output_dir/$selected_motif_file > $output_dir/$sorted_peak_file')
    TF_name_list = re.findall(r'.*fimo_motif_(.*).txt$', basename(each_file))
    TF_name =  TF_name_list[0]
    peak_bed_file = pd.read_csv(join(output_dir,sorted_peak_file), sep="\t", header=None)
    peak_count = peak_bed_file.drop_duplicates().shape[0]
    peak_bed_file = peak_bed_file.drop_duplicates()
    motif_count_dict["motif_count_nodups"].append(peak_bed_file.shape[0])
    motif_count_dict["tf_name"].append("TF_name")
       
    peak_bed = pybedtools.BedTool.from_dataframe(peak_bed_file)
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

    # if os.stat(pybed_outfile.fn).st_size == 0:
    #     master_tf_values[0].append(TF_name)
    #     master_tf_values[1].append(None)

    # else:
    ### Working with the dataframes; reading the output file of pybedtool intersect:
    df_file = pd.read_csv(pybed_outfile, sep = "\t", header = None)
    df_ordered = df_file.iloc[:, [6, 7, 8, 0, 1, 2, 3, 4]]
    df_ordered.columns = ["motif_chrom", "motif_start", "motif_end", "cpg_chrom", "cpg_start", "cpg_end", "cpg_meth", "cpg_unmeth" ]
    # df_ordered = df_ordered.drop_duplicates(subset=["motif_chrom", "motif_start", "motif_end"])
    df_ordered[["cpg_meth", "cpg_unmeth"]] = df_ordered[["cpg_meth", "cpg_unmeth"]].astype(int)
    df_ordered["strand"] = "."

    ### Making sure its in proper bed format with strand in 6th column:
    df_ordered = df_ordered.loc[:,["cpg_chrom", "cpg_start", "cpg_end", "cpg_meth", "cpg_unmeth", "strand", "motif_chrom", "motif_start", "motif_end"]]
    print "Dimension of current TF is", df_ordered.shape

    if TF_name not in master_tf_dict:
        master_tf_dict[TF_name] = df_ordered    

    # with open(join(output_dir,"final_wgbs_tf_motif_intersect_preping_for_ideas.pkl"), "w") as outfile:
    #    pickle.dump(master_tf_dict, outfile)

    # with open(join(output_dir,"final_wgbs_tf_motif_intersect_preping_for_ideas.pkl")) as infile:
    #    master_tf_dict = pickle.load(infile)

    ### Combine all the dataframes (i.e all values of master_tf_dict) representing the tf intersects with CpG:
    tf_combined_df = pd.concat(master_tf_dict).reset_index()
    tf_combined_df.to_csv(join(output_dir, "final_wgbs_tf_motif_intersect_preping_for_ideas.bed"), sep ="\t", header = False, index = False)

motif_count_df = pd.DataFrame(motif_count_dict)
motif_count_df.to_csv(join(output_dir, "final_all_tfs_non_duplicated_motif_count"), sep ="\t", header = True, index = False)

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
wgbs_tf_to_edit_file = join(output_dir, "final_wgbs_tf_motif_intersect_preping_for_ideas.bed")
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


#### Data plots info generated from hereon for motif based analysis and distribution wrt ideas regulatory regions:
print "\nGrouping the data for motif based analysis and distribution wrt ideas regulatory regions...\n"
read_ideas_df = pd.read_csv(join(output_dir, "final_wgbs_tf_ideas_intersect.bed"), sep="\t", header=None)
ideas_df_select = read_ideas_df.iloc[:,[0,1,2,3,4,6,7,8,9,13]]
ideas_df_select.columns = ["cpg_chrom", "cpg_start", "cpg_end", "cpg_meth", "cpg_unmeth", "motif_chrom", "motif_start", "motif_end", "tf_name", "ideas_state"]

map_ideas = {"Tss" : "Strong Prom", "TssF" :  "Strong Prom", "TssW" : "Weak Prom", "PromP": "Weak Prom", 
"CtcfO" : "Ctcf assoc", "Ctcf" : "Ctcf assoc", "Enh" : "Strong Enh", "EnhF" : "Strong Enh", "EnhWF1" : "Weak Enh", "EnhWF2" : "Weak Enh", "EnhWF3" : "Weak Enh", "EnhW" : "Weak Enh" } 
ideas_df_select["ideas_state_map"] =  ideas_df_select["ideas_state"].map(map_ideas)

### Selecting only the cpg motif associated to the regulatory regions:
ideas_df_select = ideas_df_select.dropna(subset=["ideas_state_map"])
# df_grouped = ideas_df_select.groupby(["tf_name", "motif_chrom", "motif_start", "motif_end", "ideas_state_map"]).apply(lambda x : x["cpg_meth"].sum()/float(x["cpg_meth"].sum() + x["cpg_unmeth"].sum()))
ideas_motif_group = ideas_df_select.groupby(["tf_name", "motif_chrom", "motif_start", "motif_end", "ideas_state_map"]).apply(lambda x : x["cpg_meth"].sum()/float(x["cpg_meth"].sum() + x["cpg_unmeth"].sum()))
ideas_motif_group_df = ideas_motif_group.reset_index(name="meth_percent")

ideas_motif_group_df["meth_group"] = " "
cond_1_low = (ideas_motif_group_df["meth_percent"] < 0.25)
cond_2_int = (ideas_motif_group_df["meth_percent"] >= 0.25) & (ideas_motif_group_df["meth_percent"] <= 0.75)
cond_3_high = (ideas_motif_group_df["meth_percent"] > 0.75)
ideas_motif_group_df.loc[cond_1_low, "meth_group"] = "Low"
ideas_motif_group_df.loc[cond_2_int, "meth_group"] = "Intermediate"
ideas_motif_group_df.loc[cond_3_high, "meth_group"] = "High"

ideas_motif_meth_group = ideas_motif_group_df.groupby(["tf_name", "meth_group", "ideas_state_map"]).size()
ideas_motif_meth_group_df = ideas_motif_meth_group.reset_index(name="meth_group_size")
ideas_motif_meth_group_df["meth_group_sizetotal"] = ideas_motif_meth_group_df.groupby(["tf_name", "meth_group"])["meth_group_size"].transform(lambda x: x.sum())   
ideas_motif_meth_group_df["meth_group_percent"] = (ideas_motif_meth_group_df["meth_group_size"]/ideas_motif_meth_group_df["meth_group_sizetotal"]*100)  
ideas_motif_meth_group_df.to_csv(join(output_dir, "final_wgbs_tf_motif_ideas_intersect_facetwrap_barplot_data.bed"), sep = "\t", header = True, index = True )


#### Data plots info generated from hereon for motif based analysis and high, low, intermediate annotation for the CpG status:
print "\nGrouping the data for motif based analysis and high, low, intermediate annotation for the CpG status\n"
motif_df = pd.read_csv(join(output_dir, "final_wgbs_tf_motif_intersect_preping_for_ideas.bed"), sep="\t", header=None)
motif_df_select = motif_df.iloc[:,[2,3,4,5,6,8,9,10,0]]
motif_df_select.columns = ["cpg_chrom", "cpg_start", "cpg_end", "cpg_meth", "cpg_unmeth", "motif_chrom", "motif_start", "motif_end", "tf_name"]

# df_grouped = motif_df_select.groupby(["tf_name", "motif_chrom", "motif_start", "motif_end", "ideas_state_map"]).apply(lambda x : x["cpg_meth"].sum()/float(x["cpg_meth"].sum() + x["cpg_unmeth"].sum()))
motif_group = motif_df_select.groupby(["tf_name", "motif_chrom", "motif_start", "motif_end"]).apply(lambda x : x["cpg_meth"].sum()/float(x["cpg_meth"].sum() + x["cpg_unmeth"].sum()))
motif_group_df = motif_group.reset_index(name="meth_percent")

motif_group_df["meth_group"] = " "
cond_1_low = (motif_group_df["meth_percent"] < 0.25)
cond_2_int = (motif_group_df["meth_percent"] >= 0.25) & (motif_group_df["meth_percent"] <= 0.75)
cond_3_high = (motif_group_df["meth_percent"] > 0.75)
motif_group_df.loc[cond_1_low, "meth_group"] = "Low"
motif_group_df.loc[cond_2_int, "meth_group"] = "Intermediate"
motif_group_df.loc[cond_3_high, "meth_group"] = "High"

motif_meth_group = motif_group_df.groupby(["tf_name", "meth_group"]).size()
motif_meth_group_df = motif_meth_group.reset_index(name="meth_group_size")
motif_meth_group_df["meth_group_sizetotal"] = motif_meth_group_df.groupby(["tf_name"])["meth_group_size"].transform(lambda x: x.sum())   
motif_meth_group_df["meth_group_percent"] = (motif_meth_group_df["meth_group_size"]/motif_meth_group_df["meth_group_sizetotal"]*100)  
motif_meth_group_df.to_csv(join(output_dir, "final_wgbs_tf_motif_ideas_intersect_HIL_barplot_data.bed"), sep = "\t", header = True, index = True )


### Data plots for the pattern of CpG in each motifs of each factor:
cpg_pattern_df = pd.read_csv(join(output_dir, "final_wgbs_tf_motif_intersect_preping_for_ideas.bed"), sep="\t", header=None)
cpg_pattern_df_select = cpg_pattern_df.iloc[:,[2,3,4,5,6,8,9,10,0]]
cpg_pattern_df_select.columns = ["cpg_chrom", "cpg_start", "cpg_end", "cpg_meth", "cpg_unmeth", "motif_chrom", "motif_start", "motif_end", "tf_name"]

cpg_pattern_df_select["total_meth_unmeth"] = cpg_pattern_df_select["cpg_meth"] + cpg_pattern_df_select["cpg_unmeth"]
cpg_pattern_df_select["meth_percent"] = cpg_pattern_df_select["cpg_meth"]/cpg_pattern_df_select["total_meth_unmeth"]
motif_cg_meth_df = cpg_pattern_df_select.loc[:, ["cpg_chrom", "cpg_start", "cpg_end", "meth_percent", "motif_chrom", "motif_start", "motif_end", "tf_name"]]

motif_cg_meth_df["meth_group"] = " "
cond_1_low = (motif_cg_meth_df["meth_percent"] < 0.25)
cond_2_int = (motif_cg_meth_df["meth_percent"] >= 0.25) & (motif_cg_meth_df["meth_percent"] <= 0.75)
cond_3_high = (motif_cg_meth_df["meth_percent"] > 0.75)
motif_cg_meth_df.loc[cond_1_low, "meth_group"] = "L"
motif_cg_meth_df.loc[cond_2_int, "meth_group"] = "I"
motif_cg_meth_df.loc[cond_3_high, "meth_group"] = "H"

motif_cg_meth_df["cumcount"] = motif_cg_meth_df.groupby(["tf_name","motif_chrom", "motif_start", "motif_end"]).cumcount() + 1
# motif_cg_pattern_for_list = motif_cg_meth_df.groupby(["tf_name", "motif_chrom", "motif_start", "motif_end"])["meth_group"].apply(list)

motif_cg_pattern = motif_cg_meth_df.groupby(["tf_name", "motif_chrom", "motif_start", "motif_end"]).apply(lambda x: str(x["cumcount"].max())+"_"+ x["meth_group"].str.cat(sep='_'))
motif_cg_pattern_df = motif_cg_pattern.reset_index(name="motif_meth_pattern")
motif_cg_pattern_df.to_csv(join(output_dir, "final_wgbs_tf_motif_ideas_intersect_cg_meth_pattern_with_coordinates.bed"), sep = "\t", header = True, index = True )

motif_cg_pattern_group = motif_cg_pattern_df.groupby(["tf_name", "motif_meth_pattern"]).size()
motif_cg_pattern_group_df = motif_cg_pattern_group.reset_index(name="meth_pattern_size")

motif_cg_pattern_group_df["meth_pattern_sizetotal"] = motif_cg_pattern_group_df.groupby(["tf_name"])["meth_pattern_size"].transform(lambda x: x.sum())   
motif_cg_pattern_group_df["meth_pattern_percent"] = (motif_cg_pattern_group_df["meth_pattern_size"]/motif_cg_pattern_group_df["meth_pattern_sizetotal"]*100)  
motif_cg_pattern_group_df_sorted = motif_cg_pattern_group_df.groupby(["tf_name"]).apply(lambda x: x.sort_values(["meth_pattern_percent"], ascending=False))
motif_cg_pattern_group_df_top_10 = motif_cg_pattern_group_df.groupby(["tf_name"])["meth_pattern_percent"].nlargest(10)
motif_cg_pattern_group_df.to_csv(join(output_dir, "final_wgbs_tf_motif_ideas_intersect_cg_meth_pattern_barplot_data.bed"), sep = "\t", header = True, index = True )
motif_cg_pattern_group_df_top_10.to_csv(join(output_dir, "final_wgbs_tf_motif_ideas_intersect_cg_meth_pattern_barplot_data_top10.bed"), sep = "\t", header = True, index = True )

#### Data plots info generated from hereon for heatmap generation with tf moitf and ideas regulatory regions:
print "\nGrouping the data for heatmap...\n"
read_ideas_df = pd.read_csv(join(output_dir, "final_wgbs_tf_ideas_intersect.bed"), sep="\t", header=None)
ideas_df_select = read_ideas_df.iloc[:,[0,1,2,3,4,6,7,8,9,13]]
ideas_df_select.columns = ["cpg_chrom", "cpg_start", "cpg_end", "cpg_meth", "cpg_unmeth", "motif_chrom", "motif_start", "motif_end", "tf_name", "ideas_state"]


"""/gpfs/gpfs1/home/schhetri/for_chris/batch_I/fimo_motifs_total/motif_meth_analysis/ideas_cpg_motif_intersect_total_corrected/final_wgbs_tf_ideas_intersect.bed"""


### Includes gene body region, uncomment below to include:
# map_ideas = {"Tss" : "Strong Prom", "TssF" :  "Strong Prom", "TssW" : "Weak Prom", "PromP": "Weak Prom", "Gen3": "Gene_body", "Gen5" : "Gene_body",
# "CtcfO" : "Ctcf assoc", "Ctcf" : "Ctcf assoc", "Enh" : "Strong Enh", "EnhF" : "Strong Enh", "EnhWF1" : "Weak Enh", "EnhWF2" : "Weak Enh", "EnhWF3" : "Weak Enh", "EnhW" : "Weak Enh" } 
# ideas_df_select["ideas_state_map"] =  ideas_df_select["ideas_state"].map(map_ideas)

### Excludes gene_body region:
map_ideas = {"Tss" : "Strong Prom", "TssF" :  "Strong Prom", "TssW" : "Weak Prom", "PromP": "Weak Prom",
"CtcfO" : "Ctcf assoc", "Ctcf" : "Ctcf assoc", "Enh" : "Strong Enh", "EnhF" : "Strong Enh", "EnhWF1" : "Weak Enh", "EnhWF2" : "Weak Enh", "EnhWF3" : "Weak Enh", "EnhW" : "Weak Enh" } 
ideas_df_select["ideas_state_map"] =  ideas_df_select["ideas_state"].map(map_ideas)

### Selecting only the cpg motif associated to the regulatory regions:
ideas_df_select = ideas_df_select.dropna(subset=["ideas_state_map"])

ideas_df_select["cpg_count"] = 1
df_grouped_with_cpg_count = ideas_df_select.groupby(["tf_name", "ideas_state_map"]).apply(lambda x : (x["cpg_meth"].sum()/float(x["cpg_meth"].sum() + x["cpg_unmeth"].sum()), x["cpg_count"].sum()))
df_grouped_with_cpg_count = df_grouped_with_cpg_count.reset_index(name="meth_perc_info")
df_grouped_with_cpg_count[["meth_percent", "cpg_count"]] = df_grouped_with_cpg_count["meth_perc_info"].apply(pd.Series)
df_grouped_with_cpg_count_20 = df_grouped_with_cpg_count.loc[df_grouped_with_cpg_count["cpg_count"] >=20, :]

final_df_unstacked_with_cpg_count_20 = df_grouped_with_cpg_count_20.pivot_table(index="tf_name", columns= "ideas_state_map", values ="meth_percent")
#final_df_unstacked_with_cpg_count_20 = df_with_cpg_count_10.set_index(["tf_name", "ideas_state_map"]).unstack()
#final_df_unstacked_with_cpg_count_20.columns = final_df_unstacked_with_cpg_count_20.columns.droplevel(0).reset_index()

print "Dimension after grouping", df_grouped.shape
df_grouped = ideas_df_select.groupby(["tf_name", "ideas_state_map"]).apply(lambda x : x["cpg_meth"].sum()/float(x["cpg_meth"].sum() + x["cpg_unmeth"].sum()))
final_df_unstacked = df_grouped.unstack()

print "Dimension after stacking", final_df_unstacked.shape
final_df_unstacked.to_csv(join(output_dir, "final_wgbs_tf_motif_ideas_intersect_heatmap_data.bed"), sep = "\t", header = True, index = True )
final_df_unstacked_with_cpg_count_20.to_csv(join(output_dir, "final_wgbs_tf_motif_ideas_intersect_heatmap_data_with_cpg_count20.bed"), sep = "\t", header = True, index = True)
df_grouped_with_cpg_count.to_csv(join(output_dir, "final_wgbs_tf_motif_ideas_intersect_heatmap_raw_data_with_cpg_count.bed"), sep = "\t", header = True, index = True )

### Only ZNF274_human not present in this list:
# df_test2 = pd.read_csv("final_wgbs_tf_motif_intersect_preping_for_ideas.bed", sep="\t")
# pd.Series(df_test2.iloc[:,0].unique())[~pd.Series(df_test2.iloc[:,0].unique()).isin(df_grouped["tf_name"])]

print "Job completed successfully!!!"




