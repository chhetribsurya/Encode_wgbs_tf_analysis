import numpy as np
import pandas as pd
import pybedtools
import re
import os
import shutil
import time
import glob
from pyfasta import Fasta
from os.path import basename
from os.path import join
from os.path import splitext

start_time = time.time()

meme_file = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/MeMe/meme.txt")
fasta_file = os.path.expanduser("~/Dropbox/local_miscellaneous_data/hg_19_misc/hg19-male.fa")
fasta_idx = Fasta(fasta_file)
# query_motif = "Motif 2"

### Can be chromHMM file or IDEAS segmentation file:
chromHMM_hepg2_file = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/chrom_impute_chromHMM_data/corrected/E118_25_imputed12marks_dense.bed")
# refgene_file = os.path.expanduser("/Users/suryachhetri/Dropbox/local_miscellaneous_data/hg_19_misc/refGene_hg19")

# cpg_file_1 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/pol2/wgbs_Pol2.cov")
# cpg_file_2 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/pol2/SL68321_Pol2_rep1.cov")
# pg_file_3 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/pol2/SL88935_Pol2_rep2.cov")

cpg_file_1 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/gata/wgbs_GATA2.cov")
cpg_file_2 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/gata/SL60583_GATA2_rep1.cov")
cpg_file_3 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/gata/SL63653_GATA2_rep2.cov")


output_dir = os.path.expanduser("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/motifs_methplot")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# motif_output_file = "final_motifs_coordinate.bed"
final_output_file = "final_wgbs_motifs_intersect_2kb"
tf_name = "GATA"

"""Returns the list of motif_regex/motif_seq"""
def parse_meme_motif_regex(MeMe_file, motif_name=False):    
    MeMe_file = meme_file
    with open(MeMe_file, "r") as motif_file:
        meme_content = motif_file.read()
        motif_name_pattern = "Motif\s\d+\sregular\sexpression"
        motif_name_list = re.compile(motif_name_pattern).findall(meme_content)

        motif_seq_pattern = "Motif\s\d+\sregular\sexpression\n-*\n(.*?\n)"
        motif_seq_list = re.compile(motif_seq_pattern).findall(meme_content)

        motif_dict = {}
        for i in range(len(motif_name_list)):
            motif_key, motif_value = " ".join(motif_name_list[i].split()[0:2]), motif_seq_list[i].strip("\n")      
            if motif_key in motif_dict:
                motif_dict[motif_key].append(motif_seq_list)
            else:
                motif_dict[motif_key] = [motif_value]

    return(motif_dict)

motif_patterns = parse_meme_motif_regex(meme_file)
# gata_regex = motif_patterns[query_motif]


# meme_content = open(meme_file, "r").read()
def parse_meme_motif_coords(MeMe_file):
    with open(meme_file, "r") as inputfile:
        meme_content = inputfile.read()
        str_pattern  = r"Motif \d+ sites sorted by position p-value.*\n.*\n.*\n.*\n"
        motif_markers_matched_list = re.compile(str_pattern).findall(meme_content)
        header = ["chrom","start","end","strand","meme_seq", "pyfasta_seq"]
        motif_header = "\t".join(header)

        Master_dict = {}
        for each_match in motif_markers_matched_list:
            motif_key = each_match[0:8].strip()
            motif_data = []
            for i in range(len(header)):
                motif_data.append([])

            motif_match_len = len(each_match)
            start_index = meme_content.find(each_match) + motif_match_len
            end_index = meme_content.find("---", start_index)
            required_each_meme_content = meme_content[start_index: end_index]
            each_motif_lines = required_each_meme_content.strip().split("\n")
            
            for motif_line in each_motif_lines:     
                splitted = motif_line.split()
                chrom1 = splitted[0].split(":")[0]
                sequence1 = splitted[5]
                strand1 = splitted[1]
                
                if strand1 == "+":          
                    start1 = (int(splitted[0].split(":")[1]) + int(splitted[2])) - 1
                    end1 = start1 + len(sequence1)
                    seq_pyfasta = fasta_idx.sequence({"chr": chrom1, "start" : start1, "stop" : end1, "strand" : strand1}, one_based = False).upper()
                    
                if strand1 == "-":
                    start1 = (int(splitted[0].split(":")[1]) + int(splitted[2])) -1 
                    end1 = start1 + len(sequence1)
                    seq_pyfasta = fasta_idx.sequence({"chr": chrom1, "start" : start1, "stop" : end1, "strand" : strand1}, one_based = False).upper()

                req = [chrom1, start1, end1, strand1, sequence1, str(seq_pyfasta)]
                for i,item in enumerate(req):
                    motif_data[i].append(item)

            motif_zip = zip(header,motif_data)
            motif_dict = dict(motif_zip)

            if motif_key not in Master_dict:
                Master_dict[motif_key] = motif_dict
        
    return(Master_dict)

Master_motif_dict = parse_meme_motif_coords(meme_file)


def combine_motif_coords_after_parsing_meme_motifs(master_motif_dict, tf_name, **kwargs):
    print "motifs to be combined..."
    combine_motif_dict = {}
    count = 0
    for key, value in kwargs.iteritems():
        motif_value = "Motif "+ str(value)  
        print key, ":", motif_value
        motif_df = pd.DataFrame(master_motif_dict[motif_value])     
        count +=1           
        if motif_value not in combine_motif_dict:
            combine_motif_dict[motif_value] = motif_df
        
    concat_motif_df = pd.concat(combine_motif_dict)
    combined_motif_df = concat_motif_df.reset_index()
    final_combined_df = combined_motif_df.iloc[:,[2,6,3,7,4,5,0]]
    final_combined_df.to_csv(join(output_dir, tf_name + "_motif_coordinate_info.bed"), sep="\t", index=False, header=False)
    print "\nTotal count of motifs combined:", count
    return(final_combined_df)

final_motif_df = combine_motif_coords_after_parsing_meme_motifs(Master_motif_dict, tf_name, motif_1= 1, motif_2= 2)


def alternative_combine_motif_coords_after_parsing_meme_motifs(master_motif_dict, **kwargs):
    print "motifs to be combined:\n"
    combine_motif_list = []
    count = 0
    for key, value in kwargs.iteritems():
        motif_value = "Motif "+ str(value)  
        print key, ":", original_value
        motif_df = pd.DataFrame(master_motif_dict[motif_value])
        combine_motif_list.append(motif_df.head())
        count +=1
    combined_motif_df = pd.concat(combine_motif_list, ignore_index=True)
    print "Total count of motifs combined:\n", count
    return(combined_motif_df)


### For trimming purpose:
def trim_and_parse_single_meme_motif_coords(MeMe_file, motif_number, trim_up, trim_down, trim_all_motifs=False):
    # Trim_up = 11
    # Trim_down = 0 
    if trim_all_motifs:
        Motif_num = "\d+"
    else:
        Motif_num = motif_number
    motif_value = "Motif " + str(Motif_num)
    meme_file = MeMe_file
    Trim_up = trim_up
    Trim_down = trim_down
    with open(meme_file, "r") as inputfile:
        meme_content = inputfile.read()
        str_pattern  = r"Motif " + str(Motif_num) + " sites sorted by position p-value.*\n.*\n.*\n.*\n"
        motif_markers_matched_list = re.compile(str_pattern).findall(meme_content)
        header = ["chrom","start","end","strand","meme_seq", "pyfasta_seq"]
        motif_header = "\t".join(header)

        Master_dict = {}
        for each_match in motif_markers_matched_list:
            motif_key = each_match[0:8].strip()
            motif_data = []
            for i in range(len(header)):
                motif_data.append([])

            motif_match_len = len(each_match)
            start_index = meme_content.find(each_match) + motif_match_len
            end_index = meme_content.find("---", start_index)
            required_each_meme_content = meme_content[start_index: end_index]
            each_motif_lines = required_each_meme_content.strip().split("\n")
            
            for motif_line in each_motif_lines:     
                splitted = motif_line.split()
                chrom1 = splitted[0].split(":")[0]
                strand1 = splitted[1]
                sequence1 = splitted[5]     
                Effective_len = (len(sequence1) - Trim_up)  

                if strand1 == "+":
                    # fasta_idx.sequence({"chr": "chr18", "start" : (3603155 + 22 -1) + Trim_up, "stop" : (3603155 + 22 -1) + (len_3 - Trim_down), "strand" : "+"}, one_based = False).upper() # u'GATAAAG'
                    start1_temp = (int(splitted[0].split(":")[1]) + int(splitted[2]) - 1)
                    start1 = start1_temp + Trim_up
                    end1 = start1_temp + (len(sequence1) - Trim_down)
                    seq_pyfasta = fasta_idx.sequence({"chr": chrom1, "start" : start1, "stop" : end1, "strand" : strand1}, one_based = False).upper()

                if strand1 == "-":
                    #fasta_idx.sequence({"chr": "chr2", "start" : (112383865 + 59 -1) + Trim_down, "stop" : (112383865 + 59 -1) + (len_3 - Trim_up), "strand" : "-"}, one_based = False).upper() # u'GATAATC'
                    start1_temp = ((int(splitted[0].split(":")[1]) + int(splitted[2])) -1)  
                    start1 =start1_temp + Trim_down
                    end1 = start1_temp + (len(sequence1) - Trim_up)
                    seq_pyfasta = fasta_idx.sequence({"chr": chrom1, "start" : start1, "stop" : end1, "strand" : strand1}, one_based = False).upper()

                req = [chrom1, start1, end1, strand1, sequence1, str(seq_pyfasta)]
                for i,item in enumerate(req):
                    motif_data[i].append(item)

            motif_zip = zip(header,motif_data)
            motif_dict = dict(motif_zip)
            motif_df = pd.DataFrame(motif_dict)
            if motif_key not in Master_dict:
                print motif_key
                Master_dict[motif_key] = motif_df
    # final_trimmed_df = pd.DataFrame(Master_dict[motif_value])
    final_trimmed_df = pd.concat(Master_dict).reset_index()
    return(final_trimmed_df)


df1 = trim_and_parse_single_meme_motif_coords(meme_file, 1, 11, 0)
df2 = trim_and_parse_single_meme_motif_coords(meme_file, 2, 0, 0)


def combine_trimmed_motif_df(tf_name, *args):
    combine_df_list = []
    for each_df in args:
        combine_df_list.append(each_df)

    combined_trimmed_df = pd.concat(combine_df_list, ignore_index=True)
    final_trimmed_motif_df = combined_trimmed_df.iloc[:,[2,6,3,7,4,5,0]]
    final_trimmed_motif_df.to_csv(join(output_dir, tf_name + "_trim_based_motif_coordinate_info.bed"), sep="\t", index=False, header=False)
    return(final_trimmed_motif_df)

final_trimmed_motif_df = combine_trimmed_motif_df(tf_name, df1,df2)

### Using boolean style of slicing:
# final_trimmed_motif_df[final_trimmed_motif_df["pyfasta_seq"].str.contains("CG")]
# final_trimmed_motif_df[final_trimmed_motif_df["pyfasta_seq"].str.contains("(CG)")]

### For regex compatible:
# [m.span() for m in re.compile("GA").finditer(str_test)]
# for each in [m.span() for m in re.compile("GA").finditer(str_test)]:
#   print each[0], each[1]
# for each_str in final_trimmed_motif_df["pyfasta_seq"].tolist():
#   [m.span() for m in re.compile("GA").finditer(each_str)]


def CG_containing_motifs_coord_df(final_motif_df_info, tf_name):
    
    final_motif_df = final_motif_df_info
    ### For finding CG containing motif dataframes, ussing boolean style of slicing:
    CG_containing_motif_df = final_motif_df[final_motif_df["meme_seq"].str.contains("CG")]
    C_lacking_motif_df = final_motif_df[~final_motif_df["meme_seq"].str.contains("C")]

    #CG_regex_match_df = final_motif_df["meme_seq"].str.extractall("(C([G]))").shape
    CG_regex_match_df = final_motif_df["meme_seq"].str.extractall("((CG))")
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


final_CG_containing_motif_df, final_CG_motif_summary_dict = CG_containing_motifs_coord_df(final_motif_df, tf_name)
# CG_containing_motif_df["meme_seq"].str.replace("CG","0")


"""CHG (where H is A, C or T)"""
def CHG_containing_motifs_coord_df(final_motif_df_info, tf_name):

    final_motif_df = final_motif_df_info
    ### For finding CHG containing motif dataframes, ussing boolean style of slicing:
    CHG_containing_motif_df = final_motif_df[final_motif_df["meme_seq"].str.contains("C[ATC]G")]
    C_lacking_motif_df = final_motif_df[~final_motif_df["meme_seq"].str.contains("C")]

    # CHG_regex_match_df = final_motif_df["meme_seq"].str.extractall("(C([ATC]G))")
    CHG_regex_match_df = final_motif_df["meme_seq"].str.extractall("((C[ATC]G))")
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

final_CHG_containing_motif_df, final_CHG_motif_summary_dict = CHG_containing_motifs_coord_df(final_motif_df, tf_name)
# CHG_containing_motif_df["meme_seq"].str.replace("C[ATC}G","0")


"""CHH (where H is A, C or T)"""
def CHH_containing_motifs_coord_df(final_motif_df_info, tf_name):

    final_motif_df = final_motif_df_info
    ### For finding CHH containing motif dataframes, ussing boolean style of slicing:
    CHH_containing_motif_df = final_motif_df[final_motif_df["meme_seq"].str.contains("C[ATC][ATC]")]
    C_lacking_motif_df = final_motif_df[~final_motif_df["meme_seq"].str.contains("C")]
    C_lacking_motif_df.to_csv(join(output_dir, tf_name+"C_lacking_motif_coordinate_info.bed"), sep="\t", index=False, header=False)

    #CHH_regex_match_df = final_motif_df["meme_seq"].str.extractall("(C([ATC][ATC]))")
    CHH_regex_match_df = final_motif_df["meme_seq"].str.extractall("((C[ATC][ATC]))")
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


final_CHH_containing_motif_df, final_CHH_motif_summary_dict = CHH_containing_motifs_coord_df(final_motif_df, tf_name)
# CHH_containing_motif_df["meme_seq"].str.replace("C[ATC][ATC]","0")


def final_motifs_model(motif_input_df):
    motif_df = motif_input_df
    #peak_df = pd.read_csv(input_file, sep="\t", header=None)
    #peak_df = pd.read_csv(input_file, sep="\t", skiprows=[0], header=None)
    select_cols = ["chrom", "start", "end", "meme_seq", "level_0"]
    motif_select_df = motif_df.loc[:,select_cols]
    motif_select_df = motif_select_df.rename(columns={"level_0":"motif_id"})
    motif_select_df["motif_id"] = map( lambda x: "_".join(x.split()), motif_select_df["motif_id"] )
    print "\nCurrent dimension of the motif model:\n", motif_select_df.shape
    motif_select_df = motif_select_df.drop_duplicates()
    print "Dropping duplicates if any, current dimension of the motif model:\n", motif_select_df.shape
    return(motif_select_df)

#final_motif_df = combine_motif_coords_after_parsing_meme_motifs(Master_motif_dict, tf_name, motif_1= 1, motif_2= 2)

motifs_coord_df = final_motifs_model(final_motif_df)


def generate_motifs_binned_coords(motifs_coordinates_info, tf_name, upstream_range, downstream_range, bin_size):
    #upstream = 1000
    #downstream = 1000
    #bin_size=100
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
    final_motifs_df = temp_motifs_df.loc[:,["chrom", "motifs_midpoint", "bin_start", "bin_end", "start", "end", "motif_id", "meme_seq"]]

    """ 
    chrom_start = tss_midpt + (bin_start); 
    chrom_end = tss_midpt + (bin_end) """
    final_motifs_df["chrom_start"] = final_motifs_df["motifs_midpoint"] + final_motifs_df["bin_start"]
    final_motifs_df["chrom_end"] = final_motifs_df["motifs_midpoint"] + final_motifs_df["bin_end"]

    select_cols = [u'chrom', u'chrom_start', u'chrom_end', u'bin_start', u'bin_end', u'motifs_midpoint', u"motif_id", u"meme_seq"]
    final_motifs_df = final_motifs_df.loc[:,select_cols]
    final_motifs_df.to_csv(join(output_dir, tf_name + "_motifs_coordinate_info.bed"), sep="\t", index=False, header=False)

    return(final_motifs_df)

motifs_midpoint_coord_df= generate_motifs_binned_coords(motifs_coord_df, tf_name, 2000, 2000, 50)


def load_motifs_coord_pybedtool_object(file_name_with_full_path): 
    each_file = file_name_with_full_path
    os.environ["output_dir"] = output_dir
    os.environ["each_file"] = each_file
    sorted_motif_file = splitext(basename(each_file))[0] + "_sorted" + splitext(basename(each_file))[1]
    os.environ["sorted_motif_file"] = sorted_motif_file

    #if not os.path.exists(join(output_dir, sorted_motif_file)):
    CMD = '''awk 'BEGIN{FS=" "; OFS="\t"} { print $1,$2,$3,$4,$5,".",$6,$7,$8 }' $each_file | sort -k1,1 -k2,2n | grep -v "^*" > $output_dir/$sorted_motif_file'''
    os.system(CMD)   
  
    print "Generating bedtools object..."
    motif_bed = pybedtools.BedTool(join(output_dir,sorted_motif_file))
 
    return(motif_bed)

motifs_bed_file = load_motifs_coord_pybedtool_object(join(output_dir, tf_name + "_motifs_coordinate_info.bed"))
    

def load_cpg_pybedtool_object(file_name_with_full_path):
    print " \nProcessing Cpg bed file\n "
    cpg_bed_file = file_name_with_full_path
    os.environ["output_dir"] = output_dir
    os.environ["cpg_bed_file"] = cpg_bed_file
    cpg_bed_edited = splitext(basename(cpg_bed_file))[0] + "_bedEdited" + splitext(cpg_bed_file)[1]
    os.environ["cpg_bed_edited"] = cpg_bed_edited

    if not os.path.exists(join(output_dir,cpg_bed_edited)):
        CMD = '''awk 'BEGIN{FS=" "; OFS="\t"} { print $1,$2,$3+1,$5,$6,"." }' $cpg_bed_file | sort -k1,1 -k2,2n | grep -v "^*" > $output_dir/$cpg_bed_edited'''
        os.system(CMD) #subprocess.call(COMMAND, shell=True) # import subprocess
    else:
        print "CpG file already converted to bed format..."

    print "Generating bedtools object..."
    Cpg_bed_file = pybedtools.BedTool(join(output_dir,cpg_bed_edited))

    return(Cpg_bed_file)


def generate_motifs_binned_perc_meth(tf_name, motifs_final_bedfile, meth_file_list, **kwargs):
    #file_name =  kwargs["files_basename"]
    print "kwargs: ", kwargs
    motifs_bedfile = motifs_final_bedfile

    master_dict = {}

    for idx, each_file in enumerate(meth_file_list):
        file_basename = splitext(basename(each_file))[0]    
        cpg_bed_file = load_cpg_pybedtool_object(each_file)

        """Pybedtool intersection of Cpg_bed_file and peak_bed_file"""
        print " Processing the intersection for", file_basename
        pybed_outfile = join(output_dir, (file_basename + tf_name + "_cpg_motifs_intersect.bed" ))
        pybed_outfile_v = join(output_dir, (file_basename + tf_name +  "_cpg_motifs_not_intersect.bed" ))

        #if not os.path.exists(pybed_outfile):
        cpg_bed_intersect = cpg_bed_file.intersect(motifs_bedfile, wa = True, wb = True, output=pybed_outfile)      
        print(cpg_bed_intersect.head())  
            #Cpg_bed_file.intersect(peak_bed, wa = True, wb = True, v = True, output=pybed_outfile_v)
        #else:
        #   print "Pybedtools object already present"

        ### Working with the dataframes; reading the output file of pybedtool intersect:
        df_file = pd.read_csv(pybed_outfile, sep = "\t", header = None)
        df_ordered = df_file.iloc[:, [0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 13, 14]]
        df_ordered.columns = ["cpg_chrom", "cpg_start", "cpg_end", "cpg_meth", "cpg_unmeth", "peak_chrom", "peak_start", "peak_end", "bin_start", "bin_end", "motif_id", "meme_seq" ]
        df_ordered[["cpg_meth", "cpg_unmeth"]] = df_ordered[["cpg_meth", "cpg_unmeth"]].astype(int)
        df_grouped =  df_ordered.groupby(["bin_start", "bin_end"]).apply(lambda x : x["cpg_meth"].sum()/float(x["cpg_meth"].sum() + x["cpg_unmeth"].sum()))
        print "Dimension of currently intersected peak file is", df_grouped.reset_index().shape

        if file_basename not in master_dict:
          master_dict[file_basename] = df_grouped

        print "\nIntersection of %s completed!!!...\n\n" %(file_basename)

    ### Combine all the dataframes (i.e all values of master_tf_dict) representing the tf intersects with CpG:
    cpg_intersect_combined_df = pd.concat(master_dict).reset_index()
    cpg_intersect_combined_df.columns = ["annotation", "bin_start", "bin_end", "meth_percent"]
    cpg_intersect_combined_df["annotation"] = tf_name + "_motif_wgbs_profile"
    cpg_intersect_combined_df.to_csv(join(output_dir, final_output_file + "_" + tf_name + ".bed"), sep ="\t", header = True, index = True)

    return(cpg_intersect_combined_df)



list_of_cpg_files = [cpg_file_2, cpg_file_3, cpg_file_1]
plot_data_set = generate_motifs_binned_perc_meth(tf_name, motifs_bed_file, list_of_cpg_files)


def main():
    tf_name = "GATA"

    motif_patterns = parse_meme_motif_regex(meme_file)
    Master_motif_dict = parse_meme_motif_coords(meme_file)
    final_motif_df = combine_motif_coords_after_parsing_meme_motifs(Master_motif_dict, tf_name, motif_1= 1, motif_2= 2)
    motifs_coord_df = final_motifs_model(final_motif_df)
    motifs_midpoint_coord_df= generate_motifs_binned_coords(motifs_coord_df, tf_name, 2000, 2000, 50)

    # final_CG_containing_motif_df, final_CG_motif_summary_dict = CG_containing_motifs_coord_df(final_motif_df, tf_name)
    # final_CHG_containing_motif_df, final_CHG_motif_summary_dict = CHG_containing_motifs_coord_df(final_motif_df, tf_name)
    # final_CHH_containing_motif_df, final_CHH_motif_summary_dict = CHH_containing_motifs_coord_df(final_motif_df, tf_name)

    motifs_bed_file = load_motifs_coord_pybedtool_object(join(output_dir, tf_name + "_motifs_coordinate_info.bed"))
    list_of_cpg_files = [cpg_file_2, cpg_file_3, cpg_file_1]
    plot_data_set = generate_motifs_binned_perc_meth(tf_name, motifs_bed_file, list_of_cpg_files)


# if __name__ == '__main__':
#   main()

# else:
#   print "Functions Imported from other module"

# print "Time for analysis = ", time.time()-start_time


############################
############################
############################


""" For motif methylation analysis of all TFs, since above codes was designed for single TF """

#meme_file_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/chip_analysis/IDR*/MeMe/*.txt"
meme_file_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/chip_analysis/test_analysis/IDR*/MeMe/*.txt"
meme_file_list = glob.glob(meme_file_path)

master_tf_dict = {}
for meme_file in meme_file_list:
        # meme_file = '/gpfs/gpfs1/home/schhetri/for_chris/batch_I/chip_analysis/IDR_SL151597_SL151598_KLF6_v2[FLAG]/MeMe/meme.txt' #string instance for regex
        tf_name = re.compile(r".*IDR_(.*?)_(.*?)_(.*?)(?=\/MeMe)").findall(meme_file)[0][-1] #[('SL151597', 'SL151598', 'KLF6_v2[FLAG]')]
        SL_rep1 = re.compile(r".*IDR_(.*?)_(.*?)_(.*?)(?=\/MeMe)").findall(meme_file)[0][0]
        SL_rep2 = re.compile(r".*IDR_(.*?)_(.*?)_(.*?)(?=\/MeMe)").findall(meme_file)[0][1]
        print "\nCurrently processing %s tf\n" %(tf_name)
        # motif_patterns = parse_meme_motif_regex(meme_file)
        motif_patterns = parse_meme_motif_regex(meme_file)
        Master_motif_dict = parse_meme_motif_coords(meme_file)
        final_motif_df = combine_motif_coords_after_parsing_meme_motifs(Master_motif_dict, tf_name, motif_1= 1, motif_2= 2)
        
        """ if interested in CG containing motifs only, then generate final_CG_containing_motif_df and process via final_motifs_model(final_CG_containing_motif_df) fn """
        final_CG_containing_motif_df, final_CG_motif_summary_dict = CG_containing_motifs_coord_df(final_motif_df, tf_name)
        final_CHG_containing_motif_df, final_CHG_motif_summary_dict = CHG_containing_motifs_coord_df(final_motif_df, tf_name)
        final_CHH_containing_motif_df, final_CHH_motif_summary_dict = CHH_containing_motifs_coord_df(final_motif_df, tf_name)
        # df1 = trim_and_parse_single_meme_motif_coords(meme_file, 1, 11, 0)
        # df2 = trim_and_parse_single_meme_motif_coords(meme_file, 2, 0, 0)
        # final_trimmed_motif_df = combine_trimmed_motif_df(df1,df2)

        motifs_coord_df = final_motifs_model(final_motif_df)
        print "Generating the binned coords..."
        motifs_midpoint_coord_df= generate_motifs_binned_coords(motifs_coord_df, tf_name, 2000, 2000, 50)
        
        motifs_bed_file = load_motifs_coord_pybedtool_object(join(output_dir, tf_name + "_motifs_coordinate_info.bed"))
        list_of_cpg_files = [cpg_file_1]
        plot_data_set = generate_motifs_binned_perc_meth(tf_name, motifs_bed_file, list_of_cpg_files)

        if tf_name not in master_tf_dict:
                master_tf_dict[tf_name] = plot_data_set

        plot_output_dir = join(output_dir,"composite_methplots") # os.path.expanduser("~/for_chris/batch_I/motifs_methplot/composite_methplots")
        if not os.path.exists(plot_output_dir):
            os.makedirs(plot_output_dir)

        wgbs_motif_intersect_file = os.path.join(output_dir, final_output_file + "_" + tf_name + ".bed")
        os.environ["wgbs_motif_file"] =  wgbs_motif_intersect_file
        os.environ["tf_name"] = tf_name
        os.environ["plot_output_dir"] = plot_output_dir
        os.system('Rscript ./Encode_motifs_bin_methplot_final.R $wgbs_motif_file $tf_name $plot_output_dir')
        #subprocess.call("Rscript " + "Encode_motifs_bin_methplot_final.R " +  wgbs_motif_file + tf_name + output_dir], shell=True)
        #subprocess.call("Rscript Encode_motifs_bin_methplot_final.R --args wgbs_motif_file tf_name output_dir", shell=True)
        print "\nRunning of Rscript for plotting completed!!!....\n"
        print "\nCheck your plots in %s dir\n" %(plot_output_dir)



master_tf_df = pd.concat(master_tf_dict, ignore_index=True)
#final_master_tf_df = master_tf_df.T
""" Final containing information of all the TFs motif methylation profile """
master_tf_df.to_csv(join(output_dir, "all_tfs_" + final_output_file), sep ="\t", header = True, index = True)

print "Task completed!!!"
print "\n\nTime taken to complete the analayis", time.time()-start_time




### For methylation profile 50bp upstream and downstream:
### For motif analysis using fastinterval and bedtools logic, so this is zero-based so the end element is exclusive:
# bin_size = 1
# upstream = 50
# downstream = 50
# bins = range(-upstream, (downstream+1), 1)

# test_data = [] 
# analysis_dict = {}


# final_data = list()
# final_header = ["chrom", "start", "end", "bin_start", "mid_point", "bin_end", "sequence", "motif_start", "motif_end"]

# for i in range(len(final_header)):
#   final_data.append(list())

# coordinate_file = os.path.expanduser("~/Desktop/pub_data/python_output_files/for_gata_parse/gata_motif.bed")
# with open(coordinate_file, "w") as outfile:
#   print "chrom\tmid_point\tstart\tend\tbin_start\tbin_end\tsequence\tmotif_start\tmotif_end"
#   for i in range(len(chrom)):
#       chrome = chrom[i]
#       print chrome
#       motif_start = int(start[i])
#       motif_end = int(end[i])
#       mid_point = (motif_start + motif_end) / 2
#       sequence = sequence_motif_pyfasta[i]

#       for each in bins :
#           bin_start = each
#           bin_end = (each + bin_size)
#           #For fastinterval and bedtools case:
#           #bin_end = (each + bin_size)
#           chrom_start = mid_point + bin_start
#           chrom_end = mid_point + bin_end

#           line_coords = [chrome, str(chrom_start), str(chrom_end), str(bin_start), str(mid_point), str(bin_end), sequence, str(motif_start), str(motif_end)]
#           test_data.append("\t".join(line_coords))

#           if chrome in analysis_dict:
#               analysis_dict[chrome].append((chrome, int(chrom_start), int(chrom_end), "\t".join(line_coords)))
#           else:
#               analysis_dict[chrome] = [(chrome, int(chrom_start), int(chrom_end), "\t".join(line_coords))]

                
#           for i,item in enumerate(line_coords):
#               final_data[i].append(item)

#           print_coords = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(chrome, chrom_start, chrom_end, bin_start, mid_point, bin_end, sequence, motif_start, motif_end) 
#           print print_coords
#           outfile.write(print_coords + "\n")

#   final_zip = zip(final_header, final_data)
#   final_dict = dict(final_zip)




