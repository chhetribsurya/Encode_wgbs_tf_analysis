#!/usr/bin/env python
import numpy as np
import pandas as pd
import pybedtools
import re
import os
import shutil
import time
import sys
import errno
from glob import glob
from pyfasta import Fasta
from os.path import basename
from os.path import join
from os.path import splitext

start_time = time.time()

#fimo_file = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/fimo/fimo.txt")
fimo_file = os.path.expanduser("~/for_chris/batch_I/fimo_motifs_total/idr_passed_motifs_total/unique_TFs/SL154016_SE_VS_SL154017_fimo_motif_CEBPA[FLAG].txt")
#fimo_file = os.path.expanduser("~/for_chris/batch_I/fimo_motifs_total/motif_meth_analysis/ideas_cpg_motif_intersect_total_corrected/SL154016_SE_VS_SL154017_fimo_motif_CEBPA[FLAG]_sorted.txt")

fasta_file = os.path.expanduser("~/for_chris/hg19-male.fa")
fasta_idx = Fasta(fasta_file)

### Can be chromHMM file or IDEAS segmentation file:
ideas_hepg2_file = os.path.expanduser("~/for_chris/batch_I/hepg2_ideas_36_dense.bed")
#chromHMM_hepg2_file = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/chrom_impute_chromHMM_data/corrected/E118_25_imputed12marks_dense.bed")
#refgene_file = os.path.expanduser("/Users/suryachhetri/Dropbox/local_miscellaneous_data/hg_19_misc/refGene_hg19")

#cpg_file_1 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/pol2/wgbs_Pol2.cov")
#cpg_file_2 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/pol2/SL68321_Pol2_rep1.cov")
#pg_file_3 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/pol2/SL88935_Pol2_rep2.cov")

#original_motif_cpg_intersect_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/fimo_motifs_total/motif_meth_analysis/cpg_motif_intersect_total"
cpg_file_1 = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/SL83768_SL83771_300x_merged_CpG_context_deduplicated.bismark.cov" 

#output_dir = os.path.expanduser("~/for_chris/batch_I/fimo_motifs_total/motif_meth_analysis/motifs_bin_methplot_total/motif_5")
#output_dir = os.path.expanduser("~/for_chris/batch_I/fimo_motifs_total/motif_meth_analysis/motifs_bin_methplot_total/unique_TFs")
output_dir = os.path.expanduser("~/for_chris/batch_I/fimo_motifs_total/motif_meth_analysis/motifs_bin_methplot_total/distinct_reg_regions")
output_sub_dir_1 = join(output_dir, "IDEAS_state_assigned_files")
output_sub_dir_2 = join(output_dir, "master_df_files")
final_output_file = "final_wgbs_motifs_w_IDEAS_state_intersect_2kb.bed"
final_output_file_2 = "final_wgbs_motifs_w_selective_IDEAS_state_intersect_2kb.bed"


### Needed for multiple threads running with a race condition to create the dir:
try:
    os.makedirs(output_dir)
except OSError as exc:
    if exc.errno != errno.EEXIST:
        raise
    pass

try:
    os.makedirs(output_sub_dir_1)
except OSError as exc:
    if exc.errno != errno.EEXIST:
        raise
    pass

try:
    os.makedirs(output_sub_dir_2)
except OSError as exc:
    if exc.errno != errno.EEXIST:
        raise
    pass


ideas_mnemonics_file = os.path.expanduser("~/for_encode/spp/analysis_scripts/ideas_table.txt")

""" Select the state number, to be analysed on """
select_state_num = range(1,10) + [15,16] + range(17,20+1) # Add +1; range() is exclusive:

""" Read the mnemonics file, and create a dictionary of states with state number ID"""
read_file = pd.read_csv(ideas_mnemonics_file, sep="\t")
state_list = read_file["Mnemonics"]
state_num_list = [ i for i in range(1,len(state_list)+1) ]
ideas_state_dict = dict(zip(state_num_list,state_list))

Target_state = [ideas_state_dict[each] for each in select_state_num] #Depends on the select_state_num variable input
vals_to_replace_dict = {value:str(key)+"_"+value for key,value in ideas_state_dict.iteritems()}



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


def parse_fimo_motif_coords(fimo_file):

    df = pd.read_csv(fimo_file, sep="\t")
    master_dict = {}
    for index, row in df.iterrows():                                       
        motif_name = str(row["#pattern name"])
        chrom_name = row["sequence name"]
        start = row["start"]
        end = row["stop"]
        strand = row["strand"]
        fimo_seq = row["matched sequence"]

        if strand == "+":          
            start = start + 1
            end = end + 1
            seq_pyfasta = fasta_idx.sequence({"chr": chrom_name, "start" : start, "stop" : end, "strand" : strand}) # one based coords
                
        elif strand == "-":
            start = start + 1
            end = end + 1
            seq_pyfasta = fasta_idx.sequence({"chr": chrom_name, "start" : start, "stop" : end, "strand" : strand}) # one based coords

        # print motif_name, chrom_name, start, end
        if motif_name not in master_dict:
            motif_name = motif_name
            master_dict[motif_name]={"chrom":[chrom_name],"start":[start], "end":[end], "strand": [strand], "fimo_seq": [fimo_seq], "pyfasta_seq":[seq_pyfasta]}
        else:
            master_dict[motif_name]["chrom"].append(chrom_name)
            master_dict[motif_name]["start"].append(start)
            master_dict[motif_name]["end"].append(end)
            master_dict[motif_name]["strand"].append(strand)
            master_dict[motif_name]["fimo_seq"].append(fimo_seq)
            master_dict[motif_name]["pyfasta_seq"].append(seq_pyfasta)
    return(master_dict)


#Master_motif_dict = parse_fimo_motif_coords(fimo_file)



def combine_motif_coords_after_parsing_fimo_motifs(master_motif_dict, tf_name, **kwargs):
    print "motifs to be combined..."
    combine_motif_dict = {}
    count = 0
    for key, value in kwargs.iteritems():
        motif_value = str(value)  
        print key, ":", motif_value
        motif_df = pd.DataFrame(master_motif_dict.get(motif_value))     
        count +=1           
        if motif_value not in combine_motif_dict:
            combine_motif_dict[motif_value] = motif_df
        
    concat_motif_df = pd.concat(combine_motif_dict)
    combined_motif_df = concat_motif_df.reset_index()
    final_combined_df = combined_motif_df.iloc[:,[2,6,3,7,4,5,0]]
    final_combined_df.to_csv(join(output_dir, tf_name + "_motif_coordinate_info.bed"), sep="\t", index=False, header=False)
    print "\nTotal count of motifs combined:", count
    return(final_combined_df)

#tf_name = "CEBPA_FLAG"
#final_motif_df = combine_motif_coords_after_parsing_fimo_motifs(Master_motif_dict, tf_name, motif_1= 1, motif_2= 2)


def CG_containing_motifs_coord_df(final_motif_df_info, tf_name):
    
    final_motif_df = final_motif_df_info
    ### For finding CG containing motif dataframes, ussing boolean style of slicing:
    CG_containing_motif_df = final_motif_df[final_motif_df["fimo_seq"].str.contains("CG")]
    C_lacking_motif_df = final_motif_df[~final_motif_df["fimo_seq"].str.contains("C")]

    #CG_regex_match_df = final_motif_df["fimo_seq"].str.extractall("(C([G]))").shape
    CG_regex_match_df = final_motif_df["fimo_seq"].str.extractall("((CG))")
    CG_regex_grouped_df = CG_regex_match_df.reset_index().groupby(["level_0"]).apply( lambda x: x["match"].max())

    final_CG_motif_joined_df = pd.concat([CG_containing_motif_df, CG_regex_grouped_df], axis=1)
    final_CG_motif_joined_df = final_CG_motif_joined_df.rename(columns={0:"CG_count"}) ## since 0 based indices so add 1 to count a match:
    final_CG_motif_joined_df["CG_count"] = final_CG_motif_joined_df["CG_count"]+1 ## since 0 based indices so add 1 to count a match:
    final_CG_motif_joined_df.to_csv(join(output_dir, tf_name + "_CG_containing_motif_coordinate_info.bed"), sep="\t", index=False, header=False)

    total_CG_motif_summary_dict = {}
    total_CG_motifs = CG_containing_motif_df.shape[0]
    total_motifs = final_motif_df.shape[0]
    total_C_lacking_motifs = C_lacking_motif_df.shape[0]

    total_CG_motif_summary_dict["total_CG_motifs"] = total_CG_motifs
    total_CG_motif_summary_dict["total_motifs"] = total_motifs
    total_CG_motif_summary_dict["CG_motif_percent"] = (total_CG_motifs/float(total_motifs)*100)
    total_CG_motif_summary_dict["total_CG_sites_from_total_CG_motifs"] = final_CG_motif_joined_df["CG_count"].sum() #or, CG_regex_match_df.shape[0]
    total_CG_motif_summary_dict["total_C_lacking_motifs"] = total_C_lacking_motifs
    total_CG_motif_summary_dict["non_C_motif_percent"] = (total_C_lacking_motifs/float(total_motifs)*100)

    return(final_CG_motif_joined_df, total_CG_motif_summary_dict)


#final_CG_containing_motif_df, final_CG_motif_summary_dict = CG_containing_motifs_coord_df(final_motif_df, tf_name)
#final_CG_containing_motif_df["fimo_seq"].str.replace("CG","0")


"""CHG (where H is A, C or T)"""
def CHG_containing_motifs_coord_df(final_motif_df_info, tf_name):

    final_motif_df = final_motif_df_info
    ### For finding CHG containing motif dataframes, ussing boolean style of slicing:
    CHG_containing_motif_df = final_motif_df[final_motif_df["fimo_seq"].str.contains("C[ATC]G")]
    C_lacking_motif_df = final_motif_df[~final_motif_df["fimo_seq"].str.contains("C")]

    # CHG_regex_match_df = final_motif_df["fimo_seq"].str.extractall("(C([ATC]G))")
    CHG_regex_match_df = final_motif_df["fimo_seq"].str.extractall("((C[ATC]G))")
    CHG_regex_grouped_df = CHG_regex_match_df.reset_index().groupby(["level_0"]).apply( lambda x: x["match"].max())

    final_CHG_motif_joined_df = pd.concat([CHG_containing_motif_df, CHG_regex_grouped_df], axis=1)
    final_CHG_motif_joined_df = final_CHG_motif_joined_df.rename(columns={0:"CHG_count"}) ## since 0 based indices so add 1 to count a match:
    final_CHG_motif_joined_df["CHG_count"] = final_CHG_motif_joined_df["CHG_count"]+1 ## since 0 based indices so add 1 to count a match:
    final_CHG_motif_joined_df.to_csv(join(output_dir, tf_name + "_CHG_containing_motif_coordinate_info.bed"), sep="\t", index=False, header=False)

    total_CHG_motif_summary_dict = {}
    total_CHG_motifs = CHG_containing_motif_df.shape[0]
    total_motifs = final_motif_df.shape[0]
    total_C_lacking_motifs = C_lacking_motif_df.shape[0]

    total_CHG_motif_summary_dict["total_CHG_motifs"] = total_CHG_motifs
    total_CHG_motif_summary_dict["total_motifs"] = total_motifs
    total_CHG_motif_summary_dict["CHG_motif_percent"] = (total_CHG_motifs/float(total_motifs)*100)
    total_CHG_motif_summary_dict["total_CHG_sites_from_total_CHG_motifs"] = final_CHG_motif_joined_df["CHG_count"].sum() #or, CHG_regex_match_df.shape[0]
    total_CHG_motif_summary_dict["total_C_lacking_motifs"] = total_C_lacking_motifs
    total_CHG_motif_summary_dict["non_C_motif_percent"] = (total_C_lacking_motifs/float(total_motifs)*100)

    return(final_CHG_motif_joined_df, total_CHG_motif_summary_dict)

#final_CHG_containing_motif_df, final_CHG_motif_summary_dict = CHG_containing_motifs_coord_df(final_motif_df, tf_name)
#final_CHG_containing_motif_df["fimo_seq"].str.replace("C[ATC}G","0")


"""CHH (where H is A, C or T)"""
def CHH_containing_motifs_coord_df(final_motif_df_info, tf_name):

    final_motif_df = final_motif_df_info
    ### For finding CHH containing motif dataframes, ussing boolean style of slicing:
    CHH_containing_motif_df = final_motif_df[final_motif_df["fimo_seq"].str.contains("C[ATC][ATC]")]
    C_lacking_motif_df = final_motif_df[~final_motif_df["fimo_seq"].str.contains("C")]
    C_lacking_motif_df.to_csv(join(output_dir, tf_name+"C_lacking_motif_coordinate_info.bed"), sep="\t", index=False, header=False)

    #CHH_regex_match_df = final_motif_df["fimo_seq"].str.extractall("(C([ATC][ATC]))")
    CHH_regex_match_df = final_motif_df["fimo_seq"].str.extractall("((C[ATC][ATC]))")
    CHH_regex_grouped_df = CHH_regex_match_df.reset_index().groupby(["level_0"]).apply( lambda x: x["match"].max())

    final_CHH_motif_joined_df = pd.concat([CHH_containing_motif_df, CHH_regex_grouped_df], axis=1)
    final_CHH_motif_joined_df = final_CHH_motif_joined_df.rename(columns={0:"CHH_count"}) ## since 0 based indices so add 1 to count a match:
    final_CHH_motif_joined_df["CHH_count"] = final_CHH_motif_joined_df["CHH_count"]+1 ## since 0 based indices so add 1 to count a match:
    final_CHH_motif_joined_df.to_csv(join(output_dir, tf_name + "_CHH_containing_motif_coordinate_info.bed"), sep="\t", index=False, header=False)
    total_C_lacking_motifs = C_lacking_motif_df.shape[0]

    total_CHH_motif_summary_dict = {}
    total_CHH_motifs = CHH_containing_motif_df.shape[0]
    total_motifs = final_motif_df.shape[0]
    total_CHH_motif_summary_dict["total_CHH_motifs"] = total_CHH_motifs
    total_CHH_motif_summary_dict["total_motifs"] = total_motifs
    total_CHH_motif_summary_dict["CHH_motif_percent"] = (total_CHH_motifs/float(total_motifs)*100)
    total_CHH_motif_summary_dict["total_CHH_sites_from_total_CHH_motifs"] = final_CHH_motif_joined_df["CHH_count"].sum() #or, CHH_regex_match_df.shape[0]
    total_CHH_motif_summary_dict["total_C_lacking_motifs"] = total_C_lacking_motifs
    total_CHH_motif_summary_dict["non_C_motif_percent"] = (total_C_lacking_motifs/float(total_motifs)*100)

    return(final_CHH_motif_joined_df, total_CHH_motif_summary_dict)


#final_CHH_containing_motif_df, final_CHH_motif_summary_dict = CHH_containing_motifs_coord_df(final_motif_df, tf_name)
#final_CHH_containing_motif_df["fimo_seq"].str.replace("C[ATC][ATC]","0")


def final_motifs_model(motif_input_df):

    #motif_df = motif_input_df
    motif_df = motif_input_df

    select_cols = ["chrom", "start", "end", "fimo_seq", "level_0"]
    motif_select_df = motif_df.loc[:,select_cols]
    motif_select_df = motif_select_df.rename(columns={"level_0":"motif_id"})
    motif_select_df["motif_id"] = map( lambda x: "_".join(x.split()), motif_select_df["motif_id"] )
    print "\nCurrent dimension of the motif model:\n", motif_select_df.shape
    motif_select_df = motif_select_df.drop_duplicates()
    print "Dropping duplicates if any, current dimension of the motif model:\n", motif_select_df.shape

    return(motif_select_df)


# final_motif_df = combine_motif_coords_after_parsing_fimo_motifs(Master_motif_dict, tf_name, motif_1= 1, motif_2= 2)
# motifs_coord_df = final_motifs_model(final_motif_df)


def assign_IDEAS_State(peaks_motif_df, cols_to_retain, tf_name):

    #peaks_motif_df = motifs_coord_df
    #cols_to_retain = [0, 1, 2, 3, 4]
    #tf_name = "CEBPA_FLAG"
    ideas_hepg2_file = os.path.expanduser("/gpfs/gpfs1/home/schhetri/encode_datasets/ideas_state/HepG2_ideas_whole_genome")
    ideas_mnemonics_file = os.path.expanduser("~/for_encode/spp/analysis_scripts/ideas_table.txt")

    #peaks_motif_df = pd.read_csv(peaks_motif_bed_file, sep = "\t", header = None)
    col_select = cols_to_retain
    peaks_motif_df = peaks_motif_df.iloc[:, col_select ]
    #peaks_motif_df.columns = ["chrom", "start", "end", "fimo_seq", "motif_id"]

    ### Intersect Hepg2 bed file and Filter the rows with max base overlap size for duplicated rows:
    hepg2_ideas_df = pd.read_csv(ideas_hepg2_file, sep="\t")
    hepg2_ideas_df = hepg2_ideas_df.iloc[:, 1:5]
    hepg2_ideas_df = hepg2_ideas_df.sort_values(["chrom", "chromStart"])

    hepg2_bed_file = pybedtools.BedTool.from_dataframe(hepg2_ideas_df)
    peaks_motif_bed_file = pybedtools.BedTool.from_dataframe(peaks_motif_df)

    motif_ideas_intersect = peaks_motif_bed_file.intersect(hepg2_bed_file, wao=True)

    ### remove the feature with 0 overlap size(last column)
    motif_ideas_intersect_df = pd.read_csv(motif_ideas_intersect.fn, sep="\t", header=None)
    intersect_df =  motif_ideas_intersect_df[motif_ideas_intersect_df.iloc[:, -1] > 0] 

    ### Filter the rows with max base overlap size for duplicated rows:(column 8 here represents the size for base overlap)
    duplicated_df = intersect_df[intersect_df.duplicated([0,1,2], keep=False)]
    last_col_ID = duplicated_df.columns.tolist()[-1]
    duplicated_filtered_df = duplicated_df.loc[duplicated_df.groupby([0,1,2])[last_col_ID].idxmax()]
    non_duplicated_df = intersect_df[~intersect_df.duplicated([0,1,2], keep=False)] # non-duplicated df
    if (intersect_df.shape[0] == duplicated_df.shape[0] + non_duplicated_df.shape[0]): # 47743 = 34578 + 13165 
        print "Duplicated Rows filtered after groupby and merged successfully"

    motif_uniq_df = pd.concat([duplicated_filtered_df, non_duplicated_df], ignore_index=True)
    second_last_col_ID = duplicated_df.columns.tolist()[-2]
    final_col_select = col_select + [second_last_col_ID]
    motif_ideas_final_df = motif_uniq_df.iloc[:,final_col_select]
    motif_ideas_final_df.columns = peaks_motif_df.columns.tolist() + ["ideas_state"]
    motif_ideas_final_df = motif_ideas_final_df.sort_values(["chrom", "start"])
    motif_ideas_final_df = motif_ideas_final_df.reset_index(drop=True) # reordering of the index
    motif_ideas_distribution = motif_ideas_final_df["ideas_state"].value_counts()

    motif_ideas_distribution.to_csv(join(output_sub_dir_1, tf_name + "_peaks_motifs_wholegenome_piechart_distribution_data.txt"), header=True, index=True, sep="\t")
    motif_ideas_final_df.to_csv(join(output_sub_dir_1, tf_name + "_peaks_motifs_annotation_with_ideas_states.bed"), header=True, index=True, sep="\t")

    return(motif_ideas_final_df)


#motifs_ideas_coord_df = assign_IDEAS_State(motifs_coord_df, [0,1,2,3,4], tf_name )


def generate_motifs_binned_coords(motifs_coordinates_info, tf_name, upstream_range, downstream_range, bin_size):
    #upstream_range = 2000
    #downstream_range = 2000
    #bin_size=50
    #motifs_coordinates_info = motifs_ideas_coord_df

    motifs_df =  motifs_coordinates_info.sort_values(["chrom","start","end"])
    upstream = upstream_range
    downstream = downstream_range
    bin_size = bin_size
    nrows =  motifs_df.shape[0]

    bins = range(-upstream, (downstream), bin_size)
    bin_len = len(bins)
    motifs_concat_df = pd.concat([motifs_df]*bin_len, ignore_index="TRUE")
    motifs_sorted_df = motifs_concat_df.sort_values(["chrom","start","end"])

    ### Copy the bin list that is deep copy:
    bin_start_list = bins[:]
    bin_end_list = []
    for each in bin_start_list:
        bin_end_list.append(each+bin_size)

    bin_df = pd.DataFrame({"bin_start":bin_start_list, "bin_end":bin_end_list})
    bin_df = bin_df.iloc[:,[1,0]] # Arrange bin_start as first col and bin_start as 2nd col.
    bin_concat_df = pd.concat([bin_df]*nrows, ignore_index="TRUE")

    ### Combine the motifs df and bin df by cbind or column wise:
    temp_motifs_df = pd.concat([motifs_sorted_df.reset_index(), bin_concat_df], axis = 1)
    temp_motifs_df["motifs_midpoint"] = (temp_motifs_df["start"] + temp_motifs_df["end"])/2
    temp_motifs_df["motifs_midpoint"] = temp_motifs_df["motifs_midpoint"].round().astype(int)
    #final_motifs_df = temp_motifs_df.loc[:,["chrom", "motifs_midpoint", "bin_start", "bin_end", "start", "end", "motif_id", "fimo_seq"]]
    final_motifs_df = temp_motifs_df.loc[:,["chrom", "motifs_midpoint", "bin_start", "bin_end", "start", "end", "motif_id", "fimo_seq", "ideas_state"]]

    """ 
    chrom_start = tss_midpt + (bin_start); 
    chrom_end = tss_midpt + (bin_end) """
    final_motifs_df["chrom_start"] = final_motifs_df["motifs_midpoint"] + final_motifs_df["bin_start"]
    final_motifs_df["chrom_end"] = final_motifs_df["motifs_midpoint"] + final_motifs_df["bin_end"]
    select_cols = [u'chrom', u'chrom_start', u'chrom_end', u'bin_start', u'bin_end', u'motifs_midpoint', u"motif_id", u"fimo_seq", "ideas_state"]
    final_motifs_df = final_motifs_df.loc[:,select_cols]

    ### Handle the case where your coordinates(mid point) are smaller than the upstream range (-2000) indicated; else start coordinate < 0(i.e -ve)
    final_motifs_df = final_motifs_df.loc[final_motifs_df["chrom_start"] > 0, :] 
    final_motifs_df.to_csv(join(output_dir, tf_name + "_binned_motifs_coordinate_info_w_IDEAS.bed"), sep="\t", index=False, header=False)

    return(final_motifs_df)


# motifs_midpoint_coord_df = generate_motifs_binned_coords(motifs_ideas_coord_df, tf_name, 2000, 2000, 50)


def load_motifs_coord_pybedtool_object(file_name_with_full_path): 
    each_file = file_name_with_full_path
    os.environ["output_dir"] = output_dir
    os.environ["each_file"] = each_file
    sorted_motif_file = splitext(basename(each_file))[0] + "_sorted" + splitext(basename(each_file))[1]
    os.environ["sorted_motif_file"] = sorted_motif_file

    #if not os.path.exists(join(output_dir, sorted_motif_file)):
    CMD = '''awk 'BEGIN{FS=" "; OFS="\t"} { print $1,$2,$3,$4,$5,".",$6,$7,$8,$9 }' $each_file | sort -k1,1 -k2,2n  > $output_dir/$sorted_motif_file'''
    os.system(CMD)   

    print "Generating bedtools object..."
    motif_bed = pybedtools.BedTool(join(output_dir,sorted_motif_file)) 
    return(motif_bed)

# motifs_bed_file = load_motifs_coord_pybedtool_object(join(output_dir, tf_name + "_binned_motifs_coordinate_info_w_IDEAS.bed"))
    

def load_cpg_pybedtool_object(file_name_with_full_path):
    print " \nProcessing Cpg bed file\n "
    cpg_bed_file = file_name_with_full_path
    os.environ["output_dir"] = output_dir
    os.environ["cpg_bed_file"] = cpg_bed_file
    cpg_bed_edited = splitext(basename(cpg_bed_file))[0] + "_edited" + splitext(cpg_bed_file)[1]
    os.environ["cpg_bed_edited"] = cpg_bed_edited

    if not os.path.exists(join(output_dir,cpg_bed_edited)):
        #CMD = '''awk 'BEGIN{FS=" "; OFS="\t"} { print $1,$2,$3+1,$5,$6,"." }' $cpg_bed_file | sort -k1,1 -k2,2n | grep -v "^*" > $output_dir/$cpg_bed_edited'''
        #os.system(CMD) #subprocess.call(COMMAND, shell=True) # import subprocess
        os.system("ln -fs $cpg_bed_file $output_dir/$cpg_bed_edited")
    else:
        print "CpG file already converted to bed format..."

    print "Generating bedtools object..."
    Cpg_bed_file = pybedtools.BedTool(join(output_dir,cpg_bed_edited))

    return(Cpg_bed_file)


def generate_motifs_binned_perc_meth(tf_name, motifs_final_bedfile, meth_file_list, **kwargs):
    #file_name =  kwargs["files_basename"]
    #tf_name = "CEBPA_FLAG"
    #motifs_final_bedfile = motifs_bed_file
    #meth_file_list = [cpg_file_1]

    print "kwargs: ", kwargs
    motifs_bedfile = motifs_final_bedfile

    master_dict = {}

    for idx, each_file in enumerate(meth_file_list):
        file_basename = splitext(basename(each_file))[0]    
        cpg_bed_file = load_cpg_pybedtool_object(each_file)

        """Pybedtool intersection of Cpg_bed_file and peak_bed_file"""
        print " Processing the intersection for", file_basename
        pybed_outfile = join(output_dir, (tf_name + "_sorted_cpg_intersect.bed" ))
        pybed_outfile_v = join(output_dir, (tf_name +  "_sorted_cpg_not_intersect.bed" ))

        if not os.path.exists(pybed_outfile): 
            cpg_bed_intersect = cpg_bed_file.intersect(motifs_bedfile, wa = True, wb = True, output=pybed_outfile)      
            cpg_bed_intersect.head()
            print "redo file:", pybed_outfile 
            #Cpg_bed_file.intersect(peak_bed, wa = True, wb = True, v = True, output=pybed_outfile_v)
        else:
            print "Pybedtools object already present"

        if os.stat(pybed_outfile).st_size == 0: # handles the empty file or interrupted and killed file on prior jobs
            cpg_bed_intersect = cpg_bed_file.intersect(motifs_bedfile, wa = True, wb = True, output=pybed_outfile)      
            cpg_bed_intersect.head()
        else:
            print "Prior incomplete or undone intersection completed"


        ### Working with the dataframes; reading the output file of pybedtool intersect:
        df_file = pd.read_csv(pybed_outfile, sep = "\t", header = None)
        df_ordered = df_file.iloc[:, [0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 13, 14, 15]]
        df_ordered.columns = ["cpg_chrom", "cpg_start", "cpg_end", "cpg_meth", "cpg_unmeth", "peak_chrom", "peak_start", "peak_end", "bin_start", "bin_end", "motif_id", "fimo_seq", "ideas_state" ]
        df_ordered[["cpg_meth", "cpg_unmeth"]] = df_ordered[["cpg_meth", "cpg_unmeth"]].astype(int)
        df_grouped =  df_ordered.groupby(["bin_start", "bin_end", "ideas_state"]).apply(lambda x : x["cpg_meth"].sum()/float(x["cpg_meth"].sum() + x["cpg_unmeth"].sum()))
        print "Dimension of currently intersected peak file is", df_grouped.reset_index().shape

        if file_basename not in master_dict:
          master_dict[file_basename] = df_grouped

        print "\nIntersection of %s completed!!!...\n\n" %(file_basename)

    ### Combine all the dataframes (i.e all values of master_tf_dict) representing the tf intersects with CpG:
    cpg_intersect_combined_df = pd.concat(master_dict).reset_index()
    cpg_intersect_combined_df.columns = ["annotation", "bin_start", "bin_end", "ideas_state",  "meth_percent"]
    cpg_intersect_combined_df["annotation"] = tf_name + "_motif_wgbs_profile"
    cpg_intersect_combined_df.to_csv(join(output_dir, final_output_file + "_" + tf_name + ".bed"), sep ="\t", header = True, index = True)

    return(cpg_intersect_combined_df)

# list_of_cpg_files = [cpg_file_2, cpg_file_3, cpg_file_1]
# list_of_cpg_files = [cpg_file_1]
# plot_data_set = generate_motifs_binned_perc_meth(tf_name, motifs_bed_file, list_of_cpg_files)


def main():
    tf_name = "CEBPA_FLAG"
    Master_motif_dict = parse_fimo_motif_coords(fimo_file)
    final_motif_df = combine_motif_coords_after_parsing_fimo_motifs(Master_motif_dict, tf_name, motif_1= 1, motif_2= 2)
    motifs_coord_df = final_motifs_model(final_motif_df)

    motifs_ideas_coord_df = assign_IDEAS_State(motifs_coord_df, [0,1,2,3,4], tf_name )
    motifs_midpoint_coord_df = generate_motifs_binned_coords(motifs_ideas_coord_df, tf_name, 2000, 2000, 100)

    #motifs_midpoint_coord_df= generate_motifs_binned_coords(motifs_coord_df, tf_name, 2000, 2000, 50)
    # final_CG_containing_motif_df, final_CG_motif_summary_dict = CG_containing_motifs_coord_df(final_motif_df, tf_name)
    # final_CHG_containing_motif_df, final_CHG_motif_summary_dict = CHG_containing_motifs_coord_df(final_motif_df, tf_name)
    # final_CHH_containing_motif_df, final_CHH_motif_summary_dict = CHH_containing_motifs_coord_df(final_motif_df, tf_name)

    motifs_bed_file = load_motifs_coord_pybedtool_object(join(output_dir, tf_name + "_binned_motifs_coordinate_info_w_IDEAS.bed"))
    #list_of_cpg_files = [cpg_file_2, cpg_file_3, cpg_file_1]
    list_of_cpg_files = [cpg_file_1]
    plot_data_set = generate_motifs_binned_perc_meth(tf_name, motifs_bed_file, list_of_cpg_files)

    plot_data_set = plot_data_set[plot_data_set["ideas_state"].isin(Target_state)]
    plot_data_set = plot_data_set.replace({"ideas_state" : vals_to_replace_dict})
    plot_data_set.to_csv(join(output_dir, final_output_file_2 + "_" + tf_name + ".bed"), sep ="\t", header = True, index = True)


# if __name__ == '__main__': main(); else: print "Functions Imported from other module"
# print "Time for analysis = ", time.time()-start_time


############################
############################
############################


""" Motif methylation analysis of all TFs, since above codes was designed for single TF """

#fimo_file_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/fimo_motifs_total/idr_passed_motifs_total/unique_TFs/test_analysis/SL*fimo_motif*"
#fimo_file_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/fimo_motifs_total/idr_passed_motifs_total/unique_TFs/SL*fimo_motif*"
#fimo_file_list = glob(fimo_file_path)

### This is useful to do run all the files parallell in the cluster at a time:
### Usage: bsub -We -n 1 -o "./fimo_motif_bin_methplot.out $RUNPATH/Encode_motifs_fimo_bin_methplot_final_total.py fimo file
fimo_file_list = [sys.argv[1]]
#fimo_file_list= ["/gpfs/gpfs1/home/schhetri/for_chris/batch_I/fimo_motifs_total/idr_passed_motifs_total/unique_TFs/test_analysis/SL984_SE_VS_SL1054_fimo_motif_RXRA.txt"]

master_tf_dict = {}
for fimo_file in fimo_file_list:
        tf_name = re.findall(r'.*fimo_motif_(.*).txt$', basename(fimo_file))[0] #[('SL151597', 'SL151598', 'KLF6_v2[FLAG]')]
        #SL_rep1 = re.compile(r".*IDR_(.*?)_(.*?)_(.*?)(?=\/fimo)").findall(fimo_file)[0][0]
        print "\nCurrently processing %s tf\n" %(tf_name)

        Master_motif_dict = parse_fimo_motif_coords(fimo_file)
        final_motif_df = combine_motif_coords_after_parsing_fimo_motifs(Master_motif_dict, tf_name, motif_1= 1, motif_2= 2)
        #final_motif_df = combine_motif_coords_after_parsing_fimo_motifs(Master_motif_dict, tf_name, motif_1= 5)
        
        """ if interested in CG containing motifs only, then generate final_CG_containing_motif_df and process via final_motifs_model(final_CG_containing_motif_df) fn """
        # final_CG_containing_motif_df, final_CG_motif_summary_dict = CG_containing_motifs_coord_df(final_motif_df, tf_name)
        # final_CHG_containing_motif_df, final_CHG_motif_summary_dict = CHG_containing_motifs_coord_df(final_motif_df, tf_name)
        # final_CHH_containing_motif_df, final_CHH_motif_summary_dict = CHH_containing_motifs_coord_df(final_motif_df, tf_name)
        # df1 = trim_and_parse_single_fimo_motif_coords(fimo_file, 1, 11, 0)
        # df2 = trim_and_parse_single_fimo_motif_coords(fimo_file, 2, 0, 0)
        # final_trimmed_motif_df = combine_trimmed_motif_df(df1,df2)

        motifs_coord_df = final_motifs_model(final_motif_df)
        motifs_ideas_coord_df = assign_IDEAS_State(motifs_coord_df, [0,1,2,3,4], tf_name ) ## Assign the ideas state and select the col to retain
        motifs_midpoint_coord_df = generate_motifs_binned_coords(motifs_ideas_coord_df, tf_name, 2000, 2000, 100)
        
        motifs_bed_file = load_motifs_coord_pybedtool_object(join(output_dir, tf_name + "_binned_motifs_coordinate_info_w_IDEAS.bed"))
        list_of_cpg_files = [cpg_file_1]
        plot_data_set = generate_motifs_binned_perc_meth(tf_name, motifs_bed_file, list_of_cpg_files)

        plot_data_set = plot_data_set[plot_data_set["ideas_state"].isin(Target_state)] ## Selective IDEAS states to be analyzed
        plot_data_set = plot_data_set.replace({"ideas_state" : vals_to_replace_dict})
        plot_data_set.to_csv(join(output_dir, final_output_file_2 + "_" + tf_name + ".bed"), sep ="\t", header = True, index = True)

        if tf_name not in master_tf_dict:
                master_tf_dict[tf_name] = plot_data_set

        plot_output_dir = join(output_dir,"composite_methplots_w_IDEAS_states") # os.path.expanduser("~/for_chris/batch_I/motifs_methplot/composite_methplots")
        if not os.path.exists(plot_output_dir):
            os.makedirs(plot_output_dir)

        wgbs_motif_intersect_file = os.path.join(output_dir, final_output_file_2 + "_" + tf_name + ".bed")
        os.environ["wgbs_motif_file"] =  wgbs_motif_intersect_file
        os.environ["tf_name"] = tf_name
        os.environ["plot_output_dir"] = plot_output_dir
        os.environ["bin_size"] = "100"
        os.environ["script_path"] = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I"
        os.system('Rscript $script_path/Encode_motifs_fimo_bin_methplot_final_total_for_distinct_regulatory_regions.R $wgbs_motif_file $tf_name $plot_output_dir $bin_size')
        ## subprocess.call("Rscript " + "Encode_motifs_bin_methplot_final.R " +  wgbs_motif_file + tf_name + output_dir], shell=True)
        ## subprocess.call("Rscript Encode_motifs_bin_methplot_final.R --args wgbs_motif_file tf_name output_dir", shell=True)
        print "\nRunning of Rscript for plotting completed!!!....\n"
        print "\nCheck your plots in %s dir\n" %(plot_output_dir)


master_tf_df = pd.concat(master_tf_dict)

""" Final containing information of all the TFs motif methylation profile """
#master_tf_df.to_csv(join(output_dir, "all_tfs_" + final_output_file), sep ="\t", header = True, index = True)

### Comment above line code and uncomment the below, if you running parallel in the cluster since master_tf_df dataframe is output of the motifs fimo bin methplot script result:
master_tf_df.to_csv(join(output_sub_dir_2, (tf_name + "_" + final_output_file)), sep ="\t", header = True, index = True)

print "Task completed!!!"
print "\n\nTime taken to complete the analayis", time.time()-start_time

