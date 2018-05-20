import pybedtools
import os
import re
import glob
import pandas as pd
import numpy as np

#Cpg_file reading via pandas_dataframe for adding the strands for standard pybed usage, and editing the coordinates:
def load_stringf_for_cpg_after_edits(input_file):

	df = pd.read_csv(input_file, sep= "\t", header= None)
	df.columns = ["chrom", "start", "end", "perc_meth", "meth", "unmeth"]
	df.loc[:,"end"] = df["end"] + 1
	df.loc[:,"strand"] = "."
	Cpg_file = df.loc[:,["chrom","start", "end", "meth", "unmeth", "strand"]]

	Cpg_list = Cpg_file.values.tolist()
	Cpg_string_list = [ "\t".join(map(str, line)) for line in Cpg_list]
	Cpg_string_bed_format = "\n".join(Cpg_string_list)

	Cpg_bed_file = pybedtools.BedTool(Cpg_string_bed_format,from_string=True)

	return(Cpg_bed_file)


#Cpg_file reading via pandas_dataframe for adding the strands for standard pybed usage, and editing the coordinates:
def load_pybedObject_for_cpg_after_edits(input_file):

	file_basename = os.path.splitext(basename(input_file))[0]
	file_ext = os.path.splitext(basename(input_file))[1]

	df = pd.read_csv(input_file, sep= "\t", header= None)
	df.columns = ["chrom", "start", "end", "perc_meth", "meth", "unmeth"]
	df.loc[:,"end"] = df["end"] + 1
	df.loc[:,"strand"] = "."
	Cpg_file = df.loc[:,["chrom","start", "end", "meth", "unmeth", "strand"]]
	
	file_name = output_dir + "/" + file_basename + "_" + "cpg_edit" + file_ext
	#Make sure that header=False else, pybedtools won't recognise the headers above:
	Cpg_file.to_csv(file_name, sep ="\t", header = False, index = False)
	Cpg_bed_file = pybedtools.BedTool(file_name)

	#Cpg_list = Cpg_bed_filter.values.tolist()
	# Cpg_list = Cpg_file.values.tolist()
	# Cpg_string_list = [ "\t".join(map(str, line)) for line in Cpg_list]
	# Cpg_string_bed_format = "\n".join(Cpg_string_list)
	#Cpg_bed_file = pybedtools.BedTool(Cpg_string_bed_format,from_string=True)

	return(Cpg_bed_file)

def meth_percent_calc(x):
    x_meth = x["cpg_meth"].sum()
    x_unmeth = x["cpg_meth"].sum() + x["cpg_unmeth"].sum()
    x_perc = x_meth/float(x_unmeth)
    return(x_perc)

#cpg_bed_list = glob.glob("SL*.cov")
peak_bed_list = glob.glob("SL*narrowPeak")
cpg_bed_file = " "
output_dir = os.path.expand("~/for_chris/batch_I/wgbs_tf_intersect")

if not os.path.exists(output_dir):
    os.makedirs(output_dir)


""" \nProcessing and creating a Cpg bed file with a strand cols from raw cpg file 
followed by pybedtool object formation...\n """
os.environ["output_dir"] = output_dir
os.environ["cpg_bed_file"] = cpg_bed_file
cpg_bed_edited = splitext(basename(cpg_bed_file))[0] + "_edited" + splitext(cpg_bed_file)[1]
os.environ["cpg_bed_edited "] = cpg_bed_edited
COMMAND = '''awk 'BEGIN{FS=" "; OFS="\t"} { print $1,$2,$3+1,$5,$6,"." }' $cpg_bed_file | sort -k1,1 -k2,2n  > $output_dir/$cpg_bed_edited'''
os.system(COMMAND) #subprocess.call(COMMAND, shell=True) # import subprocess
Cpg_bed_file = pybedtools.BedTool(cpg_bed_edited)



#df_intersect = Cpg_bed_file.intersect(peak_bed_file, wa = True, wb = True, output="intersect_cpg_bed.txt")
master_tf_dict = {}

for each_file in peak_bed_list:
    """ \nProcessing and creating a sorted bed file for the cpg intersect
    followed by pybedtools object formation...\n """
    os.environ["each_file"] = each_file
    sorted_peak_file = splitext(basename(each_file))[0] + "_sorted" + splitext(basename(each_file))[1]
    os.environ["sorted_peak_file"] = sorted_peak_file    
	  os.system('sort -k1,1, -k2,2n $each_file > $output_dir/$sorted_peak_file')
    TF_name = basename(each_file).split("_")[-1] 
    peak_bed = pybedtools.BedTool(each_file)

    """ \nPybedtool intersection of Cpg_bed_file and peak_bed_file...\n """
    pybed_outfile = join(output_dir, TF_name + "_cpg_intersect.bed")
    Cpg_bed_file.intersect(peak_bed, wa = True, wb = True, output=pybed_outfile)

    ### Working with the dataframes; reading the output file of pybedtool intersect:
    df_file = pd.read_csv(pybed_outfile, sep = "\t", header = None)
    df_ordered = df_file.iloc[:, [6, 7, 8, 0, 1, 2, 3, 4]]
    df_ordered.columns = ["peak_chrom", "peak_start", "peak_end", "cpg_chrom", "cpg_start", "cpg_end", "cpg_meth", "cpg_unmeth" ]
    df_ordered[["cpg_meth", "cpg_unmeth"]] = df_ordered[["cpg_meth", "cpg_unmeth"]].astype(int)
    df_grouped =  df_ordered.groupby(["peak_chrom", "peak_start", "peak_end"]).apply(lambda x : x["cpg_meth"].sum()/float(x["cpg_meth"].sum() + x["cpg_unmeth"].sum()))  
    #peak_grouped = df_ordered.groupby(["peak_chrom", "peak_start", "peak_end"]).apply(lambda x: x.sum())  ### For quick test of the grouped items:
    #peak_grouped = df_ordered.groupby(["peak_chrom", "peak_start", "peak_end"]).apply(meth_percent_calc)  

    if TF_name not in tf_dict:
      master_tf_dict[TF_name] = df_grouped.reset_index()

### Combine all the dataframes (i.e all values of master_tf_dict) representing the tf intersects with CpG:
tf_combined_df = pd.concat(master_dict).reset_index()
tf_combined_df.to_csv(join(output_dir, "wgbs_tf_intersect.bed", sep ="\t", header = True, index = False)
		      
