#!/usr/bin/env python

import numpy as np
import pandas as pd
import pybedtools
import re
import os
import shutil
import time
import errno
from glob import glob
from pyfasta import Fasta
from os.path import basename
from os.path import join
from os.path import splitext

start_time = time.time()

fasta_file = os.path.expanduser("~/for_chris/hg19-male.fa")
fasta_idx = Fasta(fasta_file)

output_dir = os.path.expanduser("~/for_chris/batch_I/fimo_motifs_total/motif_meth_analysis/motifs_cpg_chg_chh_content/unique_TFs_no_fanking_bases")
final_output_file = "final_cpg_chg_chh_motif_percent.txt"

### Needed for multiple threads running with a race condition to create the dir:
#if not os.path.exists(output_dir):
#   os.makedirs(output_dir)

try:
    os.makedirs(output_dir)
except OSError as exc:
    if exc.errno != errno.EEXIST:
        raise
    pass


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
            seq_pyfasta = fasta_idx.sequence({"chr": chrom_name, "start" : start, "stop" : end, "strand" : strand}) # one based
                
        elif strand == "-":
            start = start + 1
            end = end + 1 
            seq_pyfasta = fasta_idx.sequence({"chr": chrom_name, "start" : start, "stop" : end, "strand" : strand}) # one based

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
    #final_combined_df.to_csv(join(output_dir, tf_name + "_motif_coordinate_info.bed"), sep="\t", index=False, header=False)
    print "\nTotal count of motifs combined:", count
    return(final_combined_df)

#final_motif_df = combine_motif_coords_after_parsing_fimo_motifs(Master_motif_dict, tf_name, motif_1= 1, motif_2= 2)


def CG_containing_motifs_coord_df(final_motif_df_info, tf_name):
    
    final_motif_df = final_motif_df_info
    ### For finding CG containing motif dataframes, ussing boolean style of slicing:
    CG_containing_motif_df = final_motif_df[final_motif_df["pyfasta_seq"].str.contains("[Cc][Gg]")]
    C_lacking_motif_df = final_motif_df[~final_motif_df["pyfasta_seq"].str.contains("[Cc]")]
    total_C_lacking_motifs = C_lacking_motif_df.shape[0]

    #CG_regex_match_df = final_motif_df["pyfasta_seq"].str.extractall("(C([G]))").shape
    CG_regex_match_df = final_motif_df["pyfasta_seq"].str.extractall("(([Cc][Gg]))")
    total_CG_motif_summary_dict = {}

    if CG_regex_match_df.empty:
        total_CG_motif_summary_dict["total_CG_motifs"] = 0
        total_CG_motif_summary_dict["total_motifs"] = total_motifs
        total_CG_motif_summary_dict["CG_motif_percent"] = 0
        total_CG_motif_summary_dict["total_CG_sites_from_total_CG_motifs"] = 0 #or, CG_regex_match_df.shape[0]
        total_CG_motif_summary_dict["C_lacking_motifs"] = total_C_lacking_motifs
        total_CG_motif_summary_dict["non_C_motif_percent"] = (total_C_lacking_motifs/float(total_motifs)*100)


    else:
        CG_regex_grouped_df = CG_regex_match_df.reset_index().groupby(["level_0"]).apply( lambda x: x["match"].max())
        final_CG_motif_joined_df = pd.concat([CG_containing_motif_df, CG_regex_grouped_df], axis=1)
        final_CG_motif_joined_df = final_CG_motif_joined_df.rename(columns={0:"CG_count"}) ## since 0 based indices so add 1 to count a match:
        final_CG_motif_joined_df["CG_count"] = final_CG_motif_joined_df["CG_count"]+1 ## since 0 based indices so add 1 to count a match:
        final_CG_motif_joined_df.to_csv(join(output_dir, tf_name+"_CG_containing_motif_coordinate_info.bed"), sep="\t", index=False, header=False)

        total_CG_motifs = CG_containing_motif_df.shape[0]
        total_motifs = final_motif_df.shape[0]
        total_CG_motif_summary_dict["total_CG_motifs"] = total_CG_motifs
        total_CG_motif_summary_dict["total_motifs"] = total_motifs
        total_CG_motif_summary_dict["CG_motif_percent"] = (total_CG_motifs/float(total_motifs)*100)
        total_CG_motif_summary_dict["total_CG_sites_from_total_CG_motifs"] = final_CG_motif_joined_df["CG_count"].sum() #or, CG_regex_match_df.shape[0]
        total_CG_motif_summary_dict["C_lacking_motifs"] = total_C_lacking_motifs
        total_CG_motif_summary_dict["non_C_motif_percent"] = (total_C_lacking_motifs/float(total_motifs)*100)

    #return(final_CG_motif_joined_df, total_CG_motif_summary_dict)
    return(total_CG_motif_summary_dict)


#final_CG_containing_motif_df, final_CG_motif_summary_dict = CG_containing_motifs_coord_df(final_motif_df)
# CG_containing_motif_df["pyfasta_seq"].str.replace("CG","0")


"""CHG (where H is A, C or T)"""
def CHG_containing_motifs_coord_df(final_motif_df_info, tf_name):

    final_motif_df = final_motif_df_info
    ### For finding CHG containing motif dataframes, ussing boolean style of slicing:
    CHG_containing_motif_df = final_motif_df[final_motif_df["pyfasta_seq"].str.contains("[Cc][AaTtCc][Gg]")]
    C_lacking_motif_df = final_motif_df[~final_motif_df["pyfasta_seq"].str.contains("[Cc]")]
    total_C_lacking_motifs = C_lacking_motif_df.shape[0]

    #CHG_regex_match_df = final_motif_df["pyfasta_seq"].str.extractall("(C([ATC]G))")
    CHG_regex_match_df = final_motif_df["pyfasta_seq"].str.extractall("(([Cc][AaTtCc][Gg]))")
    total_CHG_motif_summary_dict = {}

    if CHG_regex_match_df.empty:
        total_CHG_motif_summary_dict["total_CHG_motifs"] = 0
        total_CHG_motif_summary_dict["total_motifs"] = total_motifs
        total_CHG_motif_summary_dict["CHG_motif_percent"] = 0
        total_CHG_motif_summary_dict["total_CHG_sites_from_total_CHG_motifs"] = 0 #or, CHG_regex_match_df.shape[0]
        total_CHG_motif_summary_dict["C_lacking_motifs"] = total_C_lacking_motifs
        total_CHG_motif_summary_dict["non_C_motif_percent"] = (total_C_lacking_motifs/float(total_motifs)*100)

    else:
        CHG_regex_grouped_df = CHG_regex_match_df.reset_index().groupby(["level_0"]).apply( lambda x: x["match"].max())
        final_CHG_motif_joined_df = pd.concat([CHG_containing_motif_df, CHG_regex_grouped_df], axis=1)
        final_CHG_motif_joined_df = final_CHG_motif_joined_df.rename(columns={0:"CHG_count"}) ## since 0 based indices so add 1 to count a match:
        final_CHG_motif_joined_df["CHG_count"] = final_CHG_motif_joined_df["CHG_count"]+1 ## since 0 based indices so add 1 to count a match:
        final_CHG_motif_joined_df.to_csv(join(output_dir, tf_name+"_CHG_containing_motif_coordinate_info.bed"), sep="\t", index=False, header=False)

        total_CHG_motifs = CHG_containing_motif_df.shape[0]
        total_motifs = final_motif_df.shape[0]
        total_CHG_motif_summary_dict["total_CHG_motifs"] = total_CHG_motifs
        total_CHG_motif_summary_dict["total_motifs"] = total_motifs
        total_CHG_motif_summary_dict["CHG_motif_percent"] = (total_CHG_motifs/float(total_motifs)*100)
        total_CHG_motif_summary_dict["total_CHG_sites_from_total_CHG_motifs"] = final_CHG_motif_joined_df["CHG_count"].sum() #or, CHG_regex_match_df.shape[0]
        total_CHG_motif_summary_dict["C_lacking_motifs"] = total_C_lacking_motifs
        total_CHG_motif_summary_dict["non_C_motif_percent"] = (total_C_lacking_motifs/float(total_motifs)*100)

    #return(final_CHG_motif_joined_df, total_CHG_motif_summary_dict)
    return(total_CHG_motif_summary_dict)

#final_CHG_containing_motif_df, final_CHG_motif_summary_dict = CHG_containing_motifs_coord_df(final_motif_df)
# CHG_containing_motif_df["pyfasta_seq"].str.replace("C[ATC}G","0")


"""CHH (where H is A, C or T)"""
def CHH_containing_motifs_coord_df(final_motif_df_info, tf_name):

    final_motif_df = final_motif_df_info
    ### For finding CHH containing motif dataframes, ussing boolean style of slicing:
    CHH_containing_motif_df = final_motif_df[final_motif_df["pyfasta_seq"].str.contains("[Cc][AaTtCc][AaTtCc]")]
    C_lacking_motif_df = final_motif_df[~final_motif_df["pyfasta_seq"].str.contains("[Cc]")]
    C_lacking_motif_df.to_csv(join(output_dir, tf_name+"C_lacking_motif_coordinate_info.bed"), sep="\t", index=False, header=False)
    total_C_lacking_motifs = C_lacking_motif_df.shape[0]

    #CHH_regex_match_df = final_motif_df["pyfasta_seq"].str.extractall("(C([ATC][ATC]))")
    CHH_regex_match_df = final_motif_df["pyfasta_seq"].str.extractall("(([Cc][AaTtCc][AaTtCc]))")
    total_CHH_motif_summary_dict = {}

    if  CHH_regex_match_df.empty:
        total_CHH_motif_summary_dict = {}
        total_CHH_motif_summary_dict["total_CHH_motifs"] = 0
        total_CHH_motif_summary_dict["total_motifs"] = total_motifs
        total_CHH_motif_summary_dict["CHH_motif_percent"] = 0
        total_CHH_motif_summary_dict["total_CHH_sites_from_total_CHH_motifs"] = 0 #or, CHH_regex_match_df.shape[0]
        total_CHH_motif_summary_dict["C_lacking_motifs"] = total_C_lacking_motifs
        total_CHH_motif_summary_dict["non_C_motif_percent"] = (total_C_lacking_motifs/float(total_motifs)*100)

    else:
        CHH_regex_grouped_df = CHH_regex_match_df.reset_index().groupby(["level_0"]).apply( lambda x: x["match"].max())
        final_CHH_motif_joined_df = pd.concat([CHH_containing_motif_df, CHH_regex_grouped_df], axis=1)
        final_CHH_motif_joined_df = final_CHH_motif_joined_df.rename(columns={0:"CHH_count"}) ## since 0 based indices so add 1 to count a match:
        final_CHH_motif_joined_df["CHH_count"] = final_CHH_motif_joined_df["CHH_count"]+1 ## since 0 based indices so add 1 to count a match:
        final_CHH_motif_joined_df.to_csv(join(output_dir, tf_name+"_CHH_containing_motif_coordinate_info.bed"), sep="\t", index=False, header=False)
        

        total_CHH_motifs = CHH_containing_motif_df.shape[0]
        total_motifs = final_motif_df.shape[0]
        total_CHH_motif_summary_dict["total_CHH_motifs"] = total_CHH_motifs
        total_CHH_motif_summary_dict["total_motifs"] = total_motifs
        total_CHH_motif_summary_dict["CHH_motif_percent"] = (total_CHH_motifs/float(total_motifs)*100)
        total_CHH_motif_summary_dict["total_CHH_sites_from_total_CHH_motifs"] = final_CHH_motif_joined_df["CHH_count"].sum() #or, CHH_regex_match_df.shape[0]
        total_CHH_motif_summary_dict["C_lacking_motifs"] = total_C_lacking_motifs
        total_CHH_motif_summary_dict["non_C_motif_percent"] = (total_C_lacking_motifs/float(total_motifs)*100)

    #return(final_CHH_motif_joined_df, total_CHH_motif_summary_dict)
    return(total_CHH_motif_summary_dict)


#final_CHH_containing_motif_df, final_CHH_motif_summary_dict = CHH_containing_motifs_coord_df(final_motif_df)
# CHH_containing_motif_df["pyfasta_seq"].str.replace("C[ATC][ATC]","0")

fimo_file_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/fimo_motifs_total/idr_passed_motifs_total/unique_TFs/SL*fimo_motif*"
fimo_file_list = glob(fimo_file_path)

master_tf_dict = {}
for fimo_file in fimo_file_list:
    tf_name = re.findall(r'.*fimo_motif_(.*).txt$', basename(fimo_file))[0] #[('SL151597', 'SL151598', 'KLF6_v2[FLAG]')]
    #SL_rep1 = re.compile(r".*IDR_(.*?)_(.*?)_(.*?)(?=\/fimo)").findall(fimo_file)[0][0]
    print "\nCurrently processing %s tf\n" %(tf_name)

    Master_motif_dict = parse_fimo_motif_coords(fimo_file)
    final_motif_df = combine_motif_coords_after_parsing_fimo_motifs(Master_motif_dict, tf_name, motif_1= 1, motif_2= 2)
    #df1 = trim_and_parse_single_fimo_motif_coords(fimo_file, 1, 11, 0)
    #df2 = trim_and_parse_single_fimo_motif_coords(fimo_file, 2, 0, 0)
    #final_trimmed_motif_df = combine_trimmed_motif_df(df1,df2)

    # final_CG_containing_motif_df, final_CG_motif_summary_dict = CG_containing_motifs_coord_df(final_motif_df, tf_name)
    # final_CHG_containing_motif_df, final_CHG_motif_summary_dict = CHG_containing_motifs_coord_df(final_motif_df, tf_name)
    # final_CHH_containing_motif_df, final_CHH_motif_summary_dict = CHH_containing_motifs_coord_df(final_motif_df, tf_name)

    final_CG_motif_summary_dict = CG_containing_motifs_coord_df(final_motif_df, tf_name)
    final_CHG_motif_summary_dict = CHG_containing_motifs_coord_df(final_motif_df, tf_name)
    final_CHH_motif_summary_dict = CHH_containing_motifs_coord_df(final_motif_df, tf_name)

    if tf_name not in master_tf_dict:
            master_tf_dict[tf_name] = final_CG_motif_summary_dict
            master_tf_dict[tf_name].update(final_CHG_motif_summary_dict)
            master_tf_dict[tf_name].update(final_CHH_motif_summary_dict)

master_tf_df = pd.DataFrame(master_tf_dict)
final_master_tf_df = master_tf_df.T
final_master_tf_df.to_csv(join(output_dir, final_output_file), sep ="\t", header = True, index = True)

print "Task completed!!!"
print "\n\nTime taken to complete the analayis", time.time()-start_time


# fimo_file="/gpfs/gpfs1/home/schhetri/for_chris/batch_I/fimo_motifs_total/motif_meth_analysis/cpg_motif_intersect_total/MBD4_sorted_cpg_intersect.bed" 
# df=pd.read_csv(fimo_file, sep="\t")                                                                                                                  
# df = df.iloc[:,:3]
# df.columns = ["chrom", "start", "end"]

# for idx, rows in df.iterrows():   
#     chrom_name = rows["chrom"]   
#     start = rows["start"]    
#     end = rows["end"]
#     print fasta_idx.sequence({"chr": chrom_name, "start" : start, "stop" : end, "strand" : "+"})
#     if idx == 100:
#         break


