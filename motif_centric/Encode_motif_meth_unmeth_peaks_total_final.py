#!/usr/bin/env python

# Follow on the output of script wgbs_tf_intersect.py ("wgbs_tf_intersect/final_wgbs_tf_intersect.bed")
import pandas as pd
import pybedtools
import glob
import os
import re
from os.path import splitext
from os.path import basename
from os.path import join

#input_file = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/wgbs_tf_intersect/final_wgbs_tf_intersect.bed"
input_dir = os.path.expanduser("~/for_chris/batch_I/wgbs_tf_intersect_total")
output_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/redo_motif_top20_meth_unmeth_peaks"

###################################################################################################################

low_meth_range = 0.20
high_meth_range = 0.80
sub_dir_name = "25_75_meth_motif"
#bsub_options = "bsub -We 10:00 -q priority -n 8 -R span[hosts=1]"
#bsub_options = "bsub -We 10:00 -q c7normal -n 8 -R span[hosts=1]"

TFs_of_interest_top20 = ["ZNF274_human", "PAF1_v1[FLAG]", "POLR2AphosphoS2_human", "CEBPA[FLAG]", "ZNF7_iso2[FLAG]", "CEBPA", "MAFF_human", "MAFK_human",\
"CEBPG[FLAG]", "ZMYM3", "ZNF12[FLAG]", "NFE2L2_human", "CEBPB_human", "MAFK_human_2", "BACH1_human", "ZC3H4[FLAG]", "HLF[FLAG]", "SSRP1[FLAG]", "ASH2L_human", "ZNF189"]

TFs_of_interest_bottom20 = ["SP1", "BRD4_human", "SP1[FLAG]", "BHLHE40", "KMT2B[FLAG]", "DRAP1[FLAG]", "SREBF1_human", "SUZ12_human", "GMEB2[FLAG]", "IRF3_human",\
"POLR2A_human_2", "MAZ_human", "NRF1_human", "BRCA1_human", "SP2", "CBX1", "ATF3", "MYC_human", "SIN3B_human", "NR3C1_human"]
##################################################################################################################

sub_output_dir = join(output_dir, sub_dir_name)

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

if not os.path.exists(sub_output_dir):
    os.makedirs(sub_output_dir)


master_tf_dict = {}
for each_tf in TFs_of_interest_top20:
    TF_name = each_tf 
    print "\nCurrently Processing %s \n" %(TF_name)

    tf_output_dir = join(sub_output_dir,TF_name)
    if not os.path.exists(tf_output_dir):
        os.makedirs(tf_output_dir)

    # lowly methylated peaks dir preparation:
    low_meth_dir = join(tf_output_dir, "low_meth")
    if not os.path.exists(low_meth_dir):
        os.makedirs(low_meth_dir)

    low_meth_temp_dir = join(low_meth_dir, "temp")
    if not os.path.exists(low_meth_temp_dir):
        os.makedirs(low_meth_temp_dir)

    low_meth_meme_dir = join(low_meth_dir, "MeMe_chip_customised")
    if not os.path.exists(low_meth_meme_dir):
        os.makedirs(low_meth_meme_dir)

    # highly methylated peaks dir preparation:
    high_meth_dir = join(tf_output_dir, "high_meth")
    if not os.path.exists(high_meth_dir):
        os.makedirs(high_meth_dir)

    high_meth_temp_dir = join(high_meth_dir, "temp")
    if not os.path.exists(high_meth_temp_dir):
        os.makedirs(high_meth_temp_dir)

    high_meth_meme_dir = join(high_meth_dir, "MeMe_chip_customised")
    if not os.path.exists(high_meth_meme_dir):
        os.makedirs(high_meth_meme_dir)

    # output of the pybed intersect, intersection b/w wgbs and TF in second step w/ wgbs_tf_intersect.py: 
    pybed_outfile = join(input_dir, (TF_name + "_sorted_cpg_intersect.bed"))

    # working with the dataframes; reading the output file of pybedtool intersect:
    df_file = pd.read_csv(pybed_outfile, sep = "\t", header = None)
    df_ordered = df_file.iloc[:, [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4]]
    df_ordered.columns = ["peak_chrom", "peak_start", "peak_end", "col_4", "col_5", "col_6", "col_7", "col_8", 
                            "col_9", "col_10", "cpg_chrom", "cpg_start", "cpg_end", "cpg_meth", "cpg_unmeth" ]
    df_ordered[["cpg_meth", "cpg_unmeth"]] = df_ordered[["cpg_meth", "cpg_unmeth"]].astype(int)

    # grouped by col_4, col_5, col_6, since motif pipeline is designed for typical SPP original peaks, so needs all the information of the peaks.
    df_grouped =  df_ordered.groupby(["peak_chrom", "peak_start", "peak_end", "col_4", "col_5", "col_6", "col_7", "col_8", 
                                    "col_9", "col_10"]).apply(lambda x : x["cpg_meth"].sum()/float(x["cpg_meth"].sum() + x["cpg_unmeth"].sum()))
    # peak_grouped = df_ordered.groupby(["peak_chrom", "peak_start", "peak_end"]).apply(meth_percent_calc)
    meth_peaks_df =  df_grouped.reset_index()
    meth_peaks_df = meth_peaks_df.rename(columns={0 : "meth_percent"})
    
    # select the peaks with low meth and high meth range:
    low_meth_peaks = meth_peaks_df[meth_peaks_df["meth_percent"] <= low_meth_range]   # (meth_peaks_df["meth_percent"] <= 0.25)
    high_meth_peaks = meth_peaks_df[meth_peaks_df["meth_percent"] >= high_meth_range] # (meth_peaks_df["meth_percent"] >= 0.75)
 
    # write the output to file:
    low_meth_peaks.to_csv(join(low_meth_dir, TF_name + "_wgbs_tf_intersect_low_meth.bed"), sep ="\t", header = False, index = False)
    low_meth_peaks.drop("meth_percent", axis=1).to_csv(join(low_meth_dir, TF_name + "_wgbs_tf_intersect_low_meth_spp_format.bed"), sep ="\t", header = False, index = False)

    high_meth_peaks.to_csv(join(high_meth_dir, TF_name + "_wgbs_tf_intersect_high_meth.bed"), sep ="\t", header = False, index = False)
    high_meth_peaks.drop("meth_percent", axis=1).to_csv(join(high_meth_dir, TF_name + "_wgbs_tf_intersect_high_meth_spp_format.bed"), sep ="\t", header = False, index = False)

    print "Dimension of current TF grouped by peaks is", df_grouped.reset_index().shape
    if TF_name not in master_tf_dict:
      master_tf_dict[TF_name] = meth_peaks_df

    print "\n Motif calling preparation for", TF_name

    # general environment variables, preparing for pipeline run:   
    #os.environ["bsub_options"] = bsub_options
    os.environ["log_name"] = "motif_calling.out"
    os.environ["reference_dir"] = "/gpfs/gpfs1/home/schhetri/for_encode/spp/reference"
    os.environ["genome"] = "hg19-male"
    os.environ["library_name"] = TF_name # Just for naming of file in temp dir
    os.environ["motif_calling_script"] = "/gpfs/gpfs1/home/schhetri/for_encode/spp/meme_chip_scripts/final_pipeline_main/python_customised_meme_scripts/call_motif_analysis_final.sh"
    
    # env. variables specific to lowly methylated region:  
    #os.environ["bjob_name"] = "Motif calling for lowly methylated peaks from " + TF_name.split("[")[0] # Because of flag suffix with TFs; HLF[FLAG].split("[")[0]
    os.environ["bjob_name"] = "Motif calling for lowly methylated peaks from " + TF_name.replace("[FLAG]", "_FLAG") # Because of flag suffix with TFs; HLF[FLAG].split("[")[0]
    os.environ["output_dir"] = low_meth_meme_dir 
    os.environ["temp_dir"] =  low_meth_temp_dir
    os.environ["peak_file"] = join(low_meth_dir, TF_name + "_wgbs_tf_intersect_low_meth.bed")
    os.system("bash motif_calling_script $peak_file $output_dir $library_name")

    # env. variables specific to highly methylated region:  
    os.environ["bjob_name"] = "Motif calling for highly methylated peaks from " + TF_name.replace("[FLAG]", "_FLAG") # Because of flag suffix with TFs; HLF[FLAG].split("[")[0]
    os.environ["log_dir"] = high_meth_dir
    os.environ["output_dir"] = high_meth_meme_dir 
    os.environ["temp_dir"] = high_meth_temp_dir
    os.environ["peak_file"] = join(high_meth_dir, TF_name + "_wgbs_tf_intersect_high_meth.bed")
    os.system("bash motif_calling_script $peak_file $output_dir $library_name")

### Combine all the dataframes (i.e all values of master_tf_dict) representing the tf intersects with CpG:
tf_combined_df = pd.concat(master_tf_dict).reset_index()
tf_combined_df.to_csv(join(sub_output_dir, "final_wgbs_tf_intersect.bed"), sep ="\t", header = True, index = False)

print "All the Jobs submitted to morgan!!!"


### For bottom 20 peaks:
#############################
#############################

input_dir = os.path.expanduser("~/for_chris/batch_I/wgbs_tf_intersect_total")
output_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/redo_motif_top20_meth_unmeth_peaks"

###################################################################################################################

low_meth_range = 0.20
high_meth_range = 0.80
sub_dir_name = "25_75_meth_motif"
#bsub_options = "bsub -We 10:00 -q priority -n 8 -R span[hosts=1]"
#bsub_options = "bsub -We 10:00 -q c7normal -n 8 -R span[hosts=1]"

TFs_of_interest_bottom20 = ["SP1", "BRD4_human", "SP1[FLAG]", "BHLHE40", "KMT2B[FLAG]", "DRAP1[FLAG]", "SREBF1_human", "SUZ12_human", "GMEB2[FLAG]", "IRF3_human",\
"POLR2A_human_2", "MAZ_human", "NRF1_human", "BRCA1_human", "SP2", "CBX1", "ATF3", "MYC_human", "SIN3B_human", "NR3C1_human"]

##################################################################################################################

sub_output_dir = join(output_dir, sub_dir_name)

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

if not os.path.exists(sub_output_dir):
    os.makedirs(sub_output_dir)

master_tf_dict = {}
for each_tf in TFs_of_interest_bottom20:
    TF_name = each_tf 
    print "\nCurrently Processing %s \n" %(TF_name)

    tf_output_dir = join(sub_output_dir,TF_name)
    if not os.path.exists(tf_output_dir):
        os.makedirs(tf_output_dir)

    # lowly methylated peaks dir preparation:
    low_meth_dir = join(tf_output_dir, "low_meth")
    if not os.path.exists(low_meth_dir):
        os.makedirs(low_meth_dir)

    low_meth_temp_dir = join(low_meth_dir, "temp")
    if not os.path.exists(low_meth_temp_dir):
        os.makedirs(low_meth_temp_dir)

    low_meth_meme_dir = join(low_meth_dir, "MeMe_chip_meth_unmeth")
    if not os.path.exists(low_meth_meme_dir):
        os.makedirs(low_meth_meme_dir)

    # highly methylated peaks dir preparation:
    high_meth_dir = join(tf_output_dir, "high_meth")
    if not os.path.exists(high_meth_dir):
        os.makedirs(high_meth_dir)

    high_meth_temp_dir = join(high_meth_dir, "temp")
    if not os.path.exists(high_meth_temp_dir):
        os.makedirs(high_meth_temp_dir)

    high_meth_meme_dir = join(high_meth_dir, "MeMe_chip_meth_unmeth")
    if not os.path.exists(high_meth_meme_dir):
        os.makedirs(high_meth_meme_dir)

    # output of the pybed intersect, intersection b/w wgbs and TF in second step w/ wgbs_tf_intersect.py: 
    pybed_outfile = join(input_dir, (TF_name + "_sorted_cpg_intersect.bed"))

    # working with the dataframes; reading the output file of pybedtool intersect:
    df_file = pd.read_csv(pybed_outfile, sep = "\t", header = None)
    df_ordered = df_file.iloc[:, [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4]]
    df_ordered.columns = ["peak_chrom", "peak_start", "peak_end", "col_4", "col_5", "col_6", "col_7", "col_8", 
                            "col_9", "col_10", "cpg_chrom", "cpg_start", "cpg_end", "cpg_meth", "cpg_unmeth" ]
    df_ordered[["cpg_meth", "cpg_unmeth"]] = df_ordered[["cpg_meth", "cpg_unmeth"]].astype(int)

    # grouped by col_4, col_5, col_6, since motif pipeline is designed for typical SPP original peaks, so needs all the information of the peaks.
    df_grouped =  df_ordered.groupby(["peak_chrom", "peak_start", "peak_end", "col_4", "col_5", "col_6", "col_7", "col_8", 
                                    "col_9", "col_10"]).apply(lambda x : x["cpg_meth"].sum()/float(x["cpg_meth"].sum() + x["cpg_unmeth"].sum()))
    # peak_grouped = df_ordered.groupby(["peak_chrom", "peak_start", "peak_end"]).apply(meth_percent_calc)
    meth_peaks_df =  df_grouped.reset_index()
    meth_peaks_df = meth_peaks_df.rename(columns={0 : "meth_percent"})
    
    # select the peaks with low meth and high meth range:
    low_meth_peaks = meth_peaks_df[meth_peaks_df["meth_percent"] <= low_meth_range]   # (meth_peaks_df["meth_percent"] <= 0.25)
    high_meth_peaks = meth_peaks_df[meth_peaks_df["meth_percent"] >= high_meth_range] # (meth_peaks_df["meth_percent"] >= 0.75)
 
    # write the output to file:
    low_meth_peaks.to_csv(join(low_meth_dir, TF_name + "_wgbs_tf_intersect_low_meth.bed"), sep ="\t", header = False, index = False)
    low_meth_peaks.drop("meth_percent", axis=1).to_csv(join(low_meth_dir, TF_name + "_wgbs_tf_intersect_low_meth_spp_format.bed"), sep ="\t", header = False, index = False)

    high_meth_peaks.to_csv(join(high_meth_dir, TF_name + "_wgbs_tf_intersect_high_meth.bed"), sep ="\t", header = False, index = False)
    high_meth_peaks.drop("meth_percent", axis=1).to_csv(join(high_meth_dir, TF_name + "_wgbs_tf_intersect_high_meth_spp_format.bed"), sep ="\t", header = False, index = False)

    print "Dimension of current TF grouped by peaks is", df_grouped.reset_index().shape
    if TF_name not in master_tf_dict:
      master_tf_dict[TF_name] = meth_peaks_df

    print "\n Motif calling preparation for", TF_name

    # general environment variables, preparing for pipeline run:   
    #os.environ["bsub_options"] = bsub_options
    os.environ["log_name"] = "motif_calling.out"
    os.environ["reference_dir"] = "/gpfs/gpfs1/home/schhetri/for_encode/spp/reference"
    os.environ["genome"] = "hg19-male"
    os.environ["library_name"] = TF_name # Just for naming of file in temp dir
    os.environ["motif_calling_script"] = "/gpfs/gpfs1/home/schhetri/for_encode/spp/meme_chip_scripts/final_pipeline_main/python_customised_meme_scripts/call_motif_analysis_final.sh"
    
    # env. variables specific to lowly methylated region:  
    #os.environ["bjob_name"] = "Motif calling for lowly methylated peaks from " + TF_name.split("[")[0] # Because of flag suffix with TFs; HLF[FLAG].split("[")[0]
    os.environ["bjob_name"] = "Motif calling for lowly methylated peaks from " + TF_name.replace("[FLAG]", "_FLAG") # Because of flag suffix with TFs; HLF[FLAG].split("[")[0]
    os.environ["output_dir"] = low_meth_meme_dir 
    os.environ["temp_dir"] =  low_meth_temp_dir
    os.environ["peak_file"] = join(low_meth_dir, TF_name + "_wgbs_tf_intersect_low_meth.bed")
    os.system("bash motif_calling_script $peak_file $output_dir $library_name")

    # env. variables specific to highly methylated region:  
    os.environ["bjob_name"] = "Motif calling for highly methylated peaks from " + TF_name.replace("[FLAG]", "_FLAG") # Because of flag suffix with TFs; HLF[FLAG].split("[")[0]
    os.environ["log_dir"] = high_meth_dir
    os.environ["output_dir"] = high_meth_meme_dir 
    os.environ["temp_dir"] = high_meth_temp_dir
    os.environ["peak_file"] = join(high_meth_dir, TF_name + "_wgbs_tf_intersect_high_meth.bed")
    os.system("bash motif_calling_script $peak_file $output_dir $library_name")

### Combine all the dataframes (i.e all values of master_tf_dict) representing the tf intersects with CpG:
tf_combined_df = pd.concat(master_tf_dict).reset_index()
tf_combined_df.to_csv(join(sub_output_dir, "final_wgbs_tf_intersect.bed"), sep ="\t", header = True, index = False)

print "All the Jobs submitted to morgan!!!"





