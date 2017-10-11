import numpy as np
import pandas as pd
import pybedtools
import re
import os
import shutil
import time
from os.path import basename
from os.path import join
from os.path import splitext


start_time = time.time()


input_file_1 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/pol2/MACS_SL68321_SL68325_15_peaks.bed")
input_file_2 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/pol2/MACS_SL68321_SL68325_15_peaks.bed")
input_file_3 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/gata/MACS_SL60583_SL60585_10_peaks.bed")
input_file_4 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/gata/MACS_SL63653_SL63655_10_peaks.bed")

### Can be chromHMM file or IDEAS segmentation file:
chromHMM_hepg2_file = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/chrom_impute_chromHMM_data/corrected/E118_25_imputed12marks_dense.bed")
#peak_input_file = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/pol2/MACS_SL68321_SL68325_15_peaks.bed") 
#refgene_file = os.path.expanduser("/Users/suryachhetri/Dropbox/local_miscellaneous_data/hg_19_misc/refGene_hg19")
#refgene_file = os.path.expanduser("~/Dropbox/local_miscellaneous_data/hg_19_misc/refGene_head.txt")
#fasta_file = os.path.expanduser("~/Dropbox/local_miscellaneous_data/hg_19_misc/hg19-male.fa")
#cpg_file_1 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/gata/wgbs_GATA2.cov")
#cpg_file_2 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/gata/SL60583_GATA2_rep1.cov")
#cpg_file_3 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/gata/SL63653_GATA2_rep2.cov")

cpg_file_1 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/pol2/wgbs_Pol2.cov")
cpg_file_2 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/pol2/SL68321_Pol2_rep1.cov")
cpg_file_3 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/pol2/SL88935_Pol2_rep2.cov")
final_wgbs_tf_chromHMM_intersect_file_rep1 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/wgbs_tf_intersect/final_wgbs_tf_intersect_preping_for_chromHMM_edited.bed")


#shutil.rmtree(os.path.expanduser("~/Dropbox/local_miscellaneous_data/pol2/chromHMM_methplot"))
output_dir = os.path.expanduser("~/Dropbox/local_miscellaneous_data/pol2/chromHMM_methplot")
if not os.path.exists(output_dir):
	os.makedirs(output_dir)
	  
final_output_file = "final_wgbs_tf_chromHMM_intersect.bed"


def intersect_all(*beds, **kwargs):

	#convert args tuple to a list, so as to have the flexibility, or mutability advantage of modifying the list later on, if needed:
	bedlist = list(beds)
	x = pybedtools.BedTool(bedlist[0])

	for bed in bedlist[1:]:
		x = x.intersect(bed, **kwargs)
	return(x)


#f(x):cleaned_peak_file_1 = intersect_all(input_file_1, input_file_2, input_file_3, input_file_4, wa = True, u = True)
#cleaned_peak_file_1.saveas('/home/surya/surya/for_gata_parse/pybed_overlap.bed')


#Merging all the cpg files with the common key (chrom, start, end):
def merge_all(*cpgs, **kwargs):

	#convert args tuple to a list, so as to have the flexibility, or mutability advantage of modifying the list later on, if needed:
	cpg_file_list = list(cpgs)
	x = pd.read_csv(cpg_file_list[0], sep="\t", names=["chrom","start","end", "perc_meth", "meth", "unmeth"])

	for cpg_file in cpg_file_list[1:]:
		current_file = pd.read_csv(cpg_file, sep="\t", names=["chrom","start","end", "perc_meth", "meth", "unmeth"])
		x = pd.merge(x, current_file, **kwargs)
	return(x)


#f(x):merged_df = merge_all(cpg_file_1, cpg_file_2, cpg_file_3, how='inner', on=['chrom', 'start', 'end'] )
#merged_df.to_csv(output_file, sep="\t",header=True, index= False)


def final_peaks_model(peak_input_file):
	input_file = peak_input_file
	peak_df = pd.read_csv(input_file, sep="\t", header=None)
	#peak_df = pd.read_csv(input_file, sep="\t", skiprows=[0], header=None)
	peak_select_df = peak_df.iloc[:, [0,1,2,3]]
	peak_select_df.columns = ["chrom", "start", "end", "state"]
	print "\nCurrent dimension of the peak model:\n", peak_select_df.shape
 	peak_select_df = peak_select_df.drop_duplicates()
 	print "Dropping duplicates if any, current dimension of the peak model:\n", peak_select_df.shape
	return(peak_select_df)


def generate_peaks_binned_coords(peaks_coordinates_info, upstream_range, downstream_range, bin_size):
	#upstream = 1000
	#downstream = 1000
	#bin_size=100
	peaks_df =  peaks_coordinates_info.sort_values(["chrom", "start", "end"])
	#peaks_df =  peaks_coordinates_info.sort_values(["chrom", "state", "start"])
	upstream = upstream_range
	downstream = downstream_range
	bin_size = bin_size
	nrows =  peaks_df.shape[0]

	bins = range(-upstream, (downstream), bin_size)
	bin_len = len(bins)
	peaks_concat_df = pd.concat([peaks_df]*bin_len, ignore_index="TRUE")
	#peaks_sorted_df = peaks_concat_df
	peaks_sorted_df = peaks_concat_df.sort_values(["chrom","start","end"])

	### Copy the bin list that is deep copy:
	bin_start_list = bins[:]
	bin_end_list = []
	for each in bin_start_list:
		bin_end_list.append(each+bin_size)

	bin_df = pd.DataFrame({"bin_start":bin_start_list, "bin_end":bin_end_list})
	bin_df = bin_df.iloc[:,[1,0]] # Arrange bin_start as first col and bin_start as 2nd col.
	bin_concat_df = pd.concat([bin_df]*nrows, ignore_index="TRUE")

	### Combine the peaks df and bin df by cbind or column wise:
	temp_peaks_df = pd.concat([peaks_sorted_df.reset_index(), bin_concat_df], axis = 1)
	temp_peaks_df["peaks_midpoint"] = (temp_peaks_df["start"] + temp_peaks_df["end"])/2
	temp_peaks_df["peaks_midpoint"] = temp_peaks_df["peaks_midpoint"].round().astype(int)
	final_peaks_df = temp_peaks_df.loc[:,["chrom", "peaks_midpoint", "bin_start", "bin_end", "start", "end", "state"]]

	""" 
	chrom_start = tss_midpt + (bin_start); 
	chrom_end = tss_midpt + (bin_end) """
	final_peaks_df["chrom_start"] = final_peaks_df["peaks_midpoint"] + final_peaks_df["bin_start"]
	final_peaks_df["chrom_end"] = final_peaks_df["peaks_midpoint"] + final_peaks_df["bin_end"]

	select_cols = [u'chrom', u'chrom_start', u'chrom_end', u'bin_start', u'bin_end', u'peaks_midpoint', u'state']
	final_peaks_df = final_peaks_df.loc[:,select_cols]
	final_peaks_df.to_csv(join(output_dir, "chromHMM_coordinate_info.bed"), sep="\t", index=False, header=False)

	return(final_peaks_df)



def load_peaks_coord_pybedtool_object(file_name_with_full_path): 
	each_file = file_name_with_full_path
	os.environ["output_dir"] = output_dir
	os.environ["each_file"] = each_file
    sorted_peak_file = splitext(basename(each_file))[0] + "_sorted" + splitext(basename(each_file))[1]
    os.environ["sorted_peak_file"] = sorted_peak_file

    #if not os.path.exists(join(output_dir, sorted_peak_file)):
    CMD = 'sort -k1,1 -k2,2n $each_file > $output_dir/$sorted_peak_file'
    os.system(CMD)   
  
    print "Generating bedtools object..."
    peak_bed = pybedtools.BedTool(join(output_dir,sorted_peak_file))
 
    return(peak_bed)



def load_wgbs_tf_intersect_pybedtool_object(file_name_with_full_path):
	print " \nProcessing Cpg tf intersected bed file\n "
	cpg_tf_bed_file = file_name_with_full_path
	os.environ["output_dir"] = output_dir
	os.environ["cpg_bed_file"] = cpg_bed_file
	cpg_bed_edited = splitext(basename(cpg_bed_file))[0] + "_bedEdited" + splitext(cpg_bed_file)[1]
	os.environ["cpg_bed_edited"] = cpg_bed_edited

	if not os.path.exists(join(output_dir,cpg_bed_edited)):
	    CMD = '''awk 'BEGIN{FS=" "; OFS="\t"} { print $1,$2,$3,$4,$5,$6,$10 }' $cpg_bed_file > $output_dir/$cpg_bed_edited'''
	    os.system(CMD) #subprocess.call(COMMAND, shell=True) # import subprocess
	else:
	    print "CpG file already converted to bed format..."

	print "Generating bedtools object..."
	Cpg_bed_file = pybedtools.BedTool(join(output_dir,cpg_bed_edited))

	return(Cpg_bed_file)



def generate_wgbs_tf_binned_perc_meth_with_ideas_chromHMM(peaks_final_bedfile, wgbs_tf_intersect_meth_file_list, **kwargs):
	#file_name =  kwargs["files_basename"]
	print "kwargs: ", kwargs
	peaks_bedfile = peaks_final_bedfile

	master_dict = {}

	for idx, each_file in enumerate(wgbs_tf_intersect_meth_file_list):
		file_basename = splitext(basename(each_file))[0]	
		cpg_bed_file = load_wgbs_tf_intersect_pybedtool_object(each_file)

	    """Pybedtool intersection of Cpg_bed_file and peak_bed_file"""
	    print " Processing the intersection for", file_basename
	    pybed_outfile = join(output_dir, (file_basename + "_cpg_chromHMM_intersect.bed"))
	    pybed_outfile_v = join(output_dir, (file_basename + "_cpg_chromHMM_not_intersect.bed"))

	    #if not os.path.exists(pybed_outfile):
	    cpg_bed_intersect = cpg_bed_file.intersect(peaks_bedfile, wa = True, wb = True, output=pybed_outfile)	    
	    print(cpg_bed_intersect.head())  
	       #Cpg_bed_file.intersect(peak_bed, wa = True, wb = True, v = True, output=pybed_outfile_v)
	    #else:
	    #	print "Pybedtools object already present"

	    ### Working with the dataframes; reading the output file of pybedtool intersect:
	    df_file = pd.read_csv(pybed_outfile, sep = "\t", header = None)
	    df_ordered = df_file.iloc[:, [0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13]]
	    df_ordered.columns = ["cpg_chrom", "cpg_start", "cpg_end", "cpg_meth", "cpg_unmeth", "tf_name", "peak_chrom", "peak_start", "peak_end", "bin_start", "bin_end", "peaks_midpoint", "state"]
	    df_ordered[["cpg_meth", "cpg_unmeth"]] = df_ordered[["cpg_meth", "cpg_unmeth"]].astype(int)
	    df_grouped =  df_ordered.groupby(["bin_start", "bin_end", "tf_name", "state"]).apply(lambda x : x["cpg_meth"].sum()/float(x["cpg_meth"].sum() + x["cpg_unmeth"].sum()))
	    print "Dimension of currently intersected peak file is", df_grouped.reset_index().shape

		if file_basename not in master_dict:
	      master_dict[file_basename] = df_grouped

		print "\nIntersection of %s completed!!!...\n\n" %(file_basename)

	### Combine all the dataframes (i.e all values of master_tf_dict) representing the tf intersects with CpG:
	cpg_tf_chromHMM_intersect_combined_df = pd.concat(master_dict).reset_index()
	cpg_tf_chromHMM_intersect_combined_df.to_csv(join(output_dir, final_output_file), sep ="\t", header = True, index = True)

	return(cpg_tf_chromHMM_intersect_combined_df)



def capture_binned_meth_perc_list(meth_file_list, motif_peak_file):	
	frames = []
	for idx, each_file in enumerate(meth_file_list):
		file_basename1 = os.path.splitext(basename(each_file))[0]	
		#cpg_bed_file = load_stringf_for_cpg_after_edits(each_file)
		cpg_bed_file = load_pybedObject_for_cpg_after_edits(each_file)
		meth_data = generate_motifs_binned_perc_meth(motif_peak_file, cpg_bed_file, file_basename=file_basename1)
		meth_data["file_name"] = os.path.splitext(basename(each_file))[0]
		frames.append(meth_data)
	appended_data = pd.concat(frames, ignore_index=True)
	appended_data.to_csv(output_dir + "/" + ggplot_outfile_suffix, sep ="\t", header = True, index = False)

	return(appended_data)


def main():
	peaks_coord_df = final_peaks_model(chromHMM_hepg2_file)
	peaks_midpoint_coord_df= generate_peaks_binned_coords(peaks_coord_df, 200, 200, 100)

	#join(output_dir, "peaks_coordinate_info.bed") : this file is a output of generate_peaks_binned_coords(), last function:
	peaks_bed_file = load_peaks_coord_pybedtool_object(join(output_dir, "chromHMM_coordinate_info.bed"))
	
	list_of_cpg_files = [final_wgbs_tf_chromHMM_intersect_file_rep1]
	plot_data_set = generate_wgbs_tf_binned_perc_meth_with_ideas_chromHMM(peaks_bed_file, list_of_cpg_files)
	#selective_gg_plottng1 = plot_data_set[plot_data_set["state"]=="10_TxEnh5'"]
	#selective_gg_plottin2 = plot_data_set[(plot_data_set["state"]=="10_TxEnh5'") | \
	#					   (plot_data_set["state"]=="13_EnhA1") ].sort_values(["state"])



if __name__ == '__main__':
	main()

else:
	print "Functions Imported from other module"

print "Time for analysis = ", time.time()-start_time
print "Task completed"



# In [169]: indices = plot_data_set[plot_data_set["state"]=="10_TxEnh5'"].index.tolist()

# In [170]: plot_data_set.iloc[indices]
# Out[170]: """
#                level_0  bin_start  bin_end       state         0
# 0    SL68321_Pol2_rep1       -200     -100  10_TxEnh5'  0.596359
# 25   SL68321_Pol2_rep1       -100        0  10_TxEnh5'  0.600159
# 50   SL68321_Pol2_rep1          0      100  10_TxEnh5'  0.601366
# 75   SL68321_Pol2_rep1        100      200  10_TxEnh5'  0.606350
# 100  SL88935_Pol2_rep2       -200     -100  10_TxEnh5'  0.590904
# 125  SL88935_Pol2_rep2       -100        0  10_TxEnh5'  0.590441
# 150  SL88935_Pol2_rep2          0      100  10_TxEnh5'  0.600995
# 175  SL88935_Pol2_rep2        100      200  10_TxEnh5'  0.589709
# 200          wgbs_Pol2       -200     -100  10_TxEnh5'  0.541527
# 225          wgbs_Pol2       -100        0  10_TxEnh5'  0.557262
# 250          wgbs_Pol2          0      100  10_TxEnh5'  0.557284
# 275          wgbs_Pol2        100      200  10_TxEnh5'  0.545581
# """
