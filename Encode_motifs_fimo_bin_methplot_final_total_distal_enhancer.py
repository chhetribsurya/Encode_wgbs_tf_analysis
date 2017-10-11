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

#fimo_file = os.path.expanduser("~/for_chris/batch_I/fimo_motifs_total/idr_passed_motifs_total/SL146880_SE_VS_SL146881_fimo_motif_CEBPA.txt")
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

output_dir = os.path.expanduser("~/for_chris/batch_I/fimo_motifs_total/motif_meth_analysis/motifs_bin_methplot_total/dist_enhancer_assoc_motifs_bin_methplot")
final_output_file = "final_wgbs_motifs_intersect_2kb.bed"

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


def filter_motif_coords_associated_to_distal_enhancer(hepg2_ideas_state_file, final_motif_df, tf_name):
    """ Read HepG2 ideas state file """
    hepg2_ideas_state_file = ideas_hepg2_file
    read_ideas_df = pd.read_csv(hepg2_ideas_state_file, sep="\t", header=None)
    ideas_df_select = read_ideas_df.iloc[:,0:4]
    ideas_df_select.columns = ["ideas_chrom", "ideas_start", "ideas_end", "ideas_state"]

    ### Map genebody, to select those regions only:
    ideas_df_genebody = ideas_df_select.copy()
    map_ideas_genebody = {"Gen5" : "Genebody", "Gen3" :  "Genebody", "Gen3Ctcf" : "Genebody", "Elon": "Genebody", "ElonW" : "Genebody"} 
    ideas_df_genebody["ideas_state_map"] =  ideas_df_genebody["ideas_state"].map(map_ideas_genebody)
    ideas_df_genebody = ideas_df_genebody.dropna(subset=["ideas_state_map"]) #rest of the state except for genebody gets eliminated   
    """Flanking the genebody regions by 2kb up and downstream"""
    ideas_df_genebody["ideas_start"] = ideas_df_genebody["ideas_start"] - 3000
    ideas_df_genebody["ideas_end"] = ideas_df_genebody["ideas_end"] + 3000

    ### Map enhancers, to select those regions only:
    ideas_df_enhancer = ideas_df_select.copy()
    map_ideas_enhancer = {"Enh" : "Enhancer", "EnhF" : "Enhancer", "EnhWF1" : "Enhancer", "EnhWF2" : "Enhancer", "EnhWF3" : "Enhancer", "EnhW" : "Enhancer"} 
    ideas_df_enhancer["ideas_state_map"] =  ideas_df_enhancer["ideas_state"].map(map_ideas_enhancer)
    ideas_df_enhancer = ideas_df_enhancer.dropna(subset=["ideas_state_map"]) #rest of the state except for genebody gets eliminated

    ### Filter the genebody regions so as to eliminate any intronic or genic enhancers:
    enhancer_bed_file = pybedtools.BedTool.from_dataframe(ideas_df_enhancer)
    genebody_bed_file = pybedtools.BedTool.from_dataframe(ideas_df_genebody)

    enhancer_genebody_outersect = enhancer_bed_file.intersect(genebody_bed_file, wa = True, wb = True, v=True)
    #enhancer_genebody_intersect = enhancer_bed_file.intersect(genebody_bed_file, wa = True, wb = True, u=True)
    
    ### Read the bedfile from outersect and save:
    distal_enhancer_coords = pd.read_csv(enhancer_genebody_outersect.fn, sep="\t", header = None)
    distal_enhancer_coords.columns = ["ideas_chrom", "ideas_start", "ideas_end", "ideas_state", "ideas_anno"]
    distal_enhancer_coords.to_csv(join(output_dir, "final_distal_enhancers_coordinates.bed"), sep="\t", header=True, index=False)

    """ Filter out only those motifs that are associated to the distal enhancers """
    motif_bed_file = pybedtools.BedTool.from_dataframe(final_motif_df)
    enhancer_final_motif_df = motif_bed_file.intersect(enhancer_genebody_outersect, wa = True, wb = True, u=True)
    filtered_final_motif_df = pd.read_csv(enhancer_final_motif_df.fn, sep="\t", header=None)
    filtered_final_motif_df.columns = final_motif_df.columns
    filtered_final_motif_df.to_csv(join(output_dir, tf_name + "_motifs_coordinate_info.bed"), sep="\t", index=False, header=False)

    return(filtered_final_motif_df)

#enh_assoc_final_motif_df = filter_motif_coords_associated_to_distal_enhancer(ideas_hepg2_file, final_motif_df, tf_name)


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
# CG_containing_motif_df["fimo_seq"].str.replace("CG","0")


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
# CHG_containing_motif_df["fimo_seq"].str.replace("C[ATC}G","0")


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
# CHH_containing_motif_df["fimo_seq"].str.replace("C[ATC][ATC]","0")


def final_motifs_model(motif_input_df):
    motif_df = motif_input_df
    select_cols = ["chrom", "start", "end", "fimo_seq", "level_0"]
    motif_select_df = motif_df.loc[:,select_cols]
    motif_select_df = motif_select_df.rename(columns={"level_0":"motif_id"})
    motif_select_df["motif_id"] = motif_select_df["motif_id"].astype(str)
    motif_select_df["motif_id"] = map( lambda x: "_".join(x.split()), motif_select_df["motif_id"] )
    print "\nCurrent dimension of the motif model:\n", motif_select_df.shape
    motif_select_df = motif_select_df.drop_duplicates()
    print "Dropping duplicates if any, current dimension of the motif model:\n", motif_select_df.shape
    return(motif_select_df)

#final_motif_df = combine_motif_coords_after_parsing_fimo_motifs(Master_motif_dict, tf_name, motif_1= 1, motif_2= 2)
#motifs_coord_df = final_motifs_model(final_motif_df)


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
    final_motifs_df = temp_motifs_df.loc[:,["chrom", "motifs_midpoint", "bin_start", "bin_end", "start", "end", "motif_id", "fimo_seq"]]

    """ 
    chrom_start = tss_midpt + (bin_start); 
    chrom_end = tss_midpt + (bin_end) """
    final_motifs_df["chrom_start"] = final_motifs_df["motifs_midpoint"] + final_motifs_df["bin_start"]
    final_motifs_df["chrom_end"] = final_motifs_df["motifs_midpoint"] + final_motifs_df["bin_end"]

    select_cols = [u'chrom', u'chrom_start', u'chrom_end', u'bin_start', u'bin_end', u'motifs_midpoint', u"motif_id", u"fimo_seq"]
    final_motifs_df = final_motifs_df.loc[:,select_cols]

    ### Handle the case where your coordinates(mid point) are smaller than the upstream range (-2000) indicated; else start coordinate < 0(i.e -ve)
    final_motifs_df = final_motifs_df.loc[final_motifs_df["chrom_start"] > 0, :] 
    final_motifs_df.to_csv(join(output_dir, tf_name + "_motifs_coordinate_info.bed"), sep="\t", index=False, header=False)

    return(final_motifs_df)

#motifs_midpoint_coord_df= generate_motifs_binned_coords(motifs_coord_df, tf_name, 2000, 2000, 50)


def load_motifs_coord_pybedtool_object(file_name_with_full_path): 
    each_file = file_name_with_full_path
    os.environ["output_dir"] = output_dir
    os.environ["each_file"] = each_file
    sorted_motif_file = splitext(basename(each_file))[0] + "_sorted" + splitext(basename(each_file))[1]
    os.environ["sorted_motif_file"] = sorted_motif_file

    #if not os.path.exists(join(output_dir, sorted_motif_file)):
    CMD = '''awk 'BEGIN{FS=" "; OFS="\t"} { print $1,$2,$3,$4,$5,".",$6,$7,$8 }' $each_file | sort -k1,1 -k2,2n  > $output_dir/$sorted_motif_file'''
    os.system(CMD)   
  
    print "Generating bedtools object..."
    motif_bed = pybedtools.BedTool(join(output_dir,sorted_motif_file))
 
    return(motif_bed)

#motifs_bed_file = load_motifs_coord_pybedtool_object(join(output_dir, tf_name + "_motifs_coordinate_info.bed"))


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
            #os.environ["output_dir"] = output_dir
            #os.environ["motif_cpg_bed_file"] = original_motif_cpg_intersect_dir/basename(pybed_outfile)
            #os.system("ln -fs $motif_cpg_bed_file $output_dir/$cpg_bed_file")
            cpg_bed_intersect = cpg_bed_file.intersect(motifs_bedfile, wa = True, wb = True, output=pybed_outfile)      
            cpg_bed_intersect.head() 
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
        df_ordered = df_file.iloc[:, [0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 13, 14]]
        df_ordered.columns = ["cpg_chrom", "cpg_start", "cpg_end", "cpg_meth", "cpg_unmeth", "peak_chrom", "peak_start", "peak_end", "bin_start", "bin_end", "motif_id", "fimo_seq" ]
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


# list_of_cpg_files = [cpg_file_2, cpg_file_3, cpg_file_1]
# plot_data_set = generate_motifs_binned_perc_meth(tf_name, motifs_bed_file, list_of_cpg_files)


def main():
    #tf_name = "GATA"
    Master_motif_dict = parse_fimo_motif_coords(fimo_file)
    final_motif_df = combine_motif_coords_after_parsing_fimo_motifs(Master_motif_dict, tf_name, motif_1= 1, motif_2= 2)
    motifs_coord_df = final_motifs_model(final_motif_df)
    motifs_midpoint_coord_df= generate_motifs_binned_coords(motifs_coord_df, tf_name, 2000, 2000, 50)

    # final_CG_containing_motif_df, final_CG_motif_summary_dict = CG_containing_motifs_coord_df(final_motif_df, tf_name)
    # final_CHG_containing_motif_df, final_CHG_motif_summary_dict = CHG_containing_motifs_coord_df(final_motif_df, tf_name)
    # final_CHH_containing_motif_df, final_CHH_motif_summary_dict = CHH_containing_motifs_coord_df(final_motif_df, tf_name)

    motifs_bed_file = load_motifs_coord_pybedtool_object(join(output_dir, tf_name + "_motifs_coordinate_info.bed"))
    list_of_cpg_files = [cpg_file_2, cpg_file_3, cpg_file_1]
    plot_data_set = generate_motifs_binned_perc_meth(tf_name, motifs_bed_file, list_of_cpg_files)


# if __name__ == '__main__': main(); else: print "Functions Imported from other module"
# print "Time for analysis = ", time.time()-start_time


############################
############################
############################


""" Motif methylation analysis of all TFs, since above codes was designed for single TF """

fimo_file_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/fimo_motifs_total/idr_passed_motifs_total/unique_TFs/test_analysis/SL*fimo_motif*"
fimo_file_list = glob(fimo_file_path)

### This is useful to do run all the files parallell in the cluster at a time:
### Usage: bsub -We -n 1 -o "./fimo_motif_bin_methplot.out $RUNPATH/Encode_motifs_fimo_bin_methplot_final_total.py fimo file
#fimo_file_list = [sys.argv[1]]

master_tf_dict = {}
for fimo_file in fimo_file_list:
        tf_name = re.findall(r'.*fimo_motif_(.*).txt$', basename(fimo_file))[0] #[('SL151597', 'SL151598', 'KLF6_v2[FLAG]')]
        #SL_rep1 = re.compile(r".*IDR_(.*?)_(.*?)_(.*?)(?=\/fimo)").findall(fimo_file)[0][0]
        print "\nCurrently processing %s tf\n" %(tf_name)

        Master_motif_dict = parse_fimo_motif_coords(fimo_file)
        final_motif_df = combine_motif_coords_after_parsing_fimo_motifs(Master_motif_dict, tf_name, motif_1= 1, motif_2= 2)
        enh_assoc_final_motif_df = filter_motif_coords_associated_to_distal_enhancer(ideas_hepg2_file, final_motif_df, tf_name)

        """ if interested in CG containing motifs only, then generate final_CG_containing_motif_df and process via final_motifs_model(final_CG_containing_motif_df) fn """
        # final_CG_containing_motif_df, final_CG_motif_summary_dict = CG_containing_motifs_coord_df(final_motif_df, tf_name)
        # final_CHG_containing_motif_df, final_CHG_motif_summary_dict = CHG_containing_motifs_coord_df(final_motif_df, tf_name)
        # final_CHH_containing_motif_df, final_CHH_motif_summary_dict = CHH_containing_motifs_coord_df(final_motif_df, tf_name)
        # df1 = trim_and_parse_single_fimo_motif_coords(fimo_file, 1, 11, 0)
        # df2 = trim_and_parse_single_fimo_motif_coords(fimo_file, 2, 0, 0)
        # final_trimmed_motif_df = combine_trimmed_motif_df(df1,df2)

        motifs_coord_df = final_motifs_model(enh_assoc_final_motif_df)
        motifs_midpoint_coord_df= generate_motifs_binned_coords(motifs_coord_df, tf_name, 2000, 2000, 50)
        #motifs_midpoint_coord_df= generate_motifs_binned_coords(motifs_coord_df, tf_name, 100, 100, 1)
        
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
        ## subprocess.call("Rscript " + "Encode_motifs_bin_methplot_final.R " +  wgbs_motif_file + tf_name + output_dir], shell=True)
        ## subprocess.call("Rscript Encode_motifs_bin_methplot_final.R --args wgbs_motif_file tf_name output_dir", shell=True)
        print "\nRunning of Rscript for plotting completed!!!....\n"
        print "\nCheck your plots in %s dir\n" %(plot_output_dir)


master_tf_df = pd.concat(master_tf_dict)

""" Final containing information of all the TFs motif methylation profile """
master_tf_df.to_csv(join(output_dir, "all_tfs_" + final_output_file), sep ="\t", header = True, index = True)

### Comment above line code and uncomment the below, if you running parallel in the cluster since master_tf_df dataframe is output of the motifs fimo bin methplot script result:
#master_tf_df.to_csv(join(output_dir, (tf_name + "_" + final_output_file)), sep ="\t", header = True, index = True)

print "Task completed!!!"
print "\n\nTime taken to complete the analayis", time.time()-start_time



### Concatenate all the motifs +2kb -2kb associated stuffs: 
#output_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/fimo_motifs_total/motif_meth_analysis/motifs_bin_methplot_total/motif_5"
output_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/fimo_motifs_total/motif_meth_analysis/motifs_bin_methplot_total/dist_enhancer_assoc_motifs_bin_methplot"
motif_file_list = glob(join(output_dir, "final_wgbs_motifs_intersect*"))

motif_file_list_df = []
for each_file in motif_file_list:
    df = pd.read_csv(each_file, sep="\t")
    motif_file_list_df.append(df)

df_combined = pd.concat(motif_file_list_df, ignore_index=True)
df_combined["tf_name"] = df_combined["annotation"].apply(lambda x: x.split("_motif")[0])
motif_bin_methplot_df = df_combined.pivot(index="tf_name", columns="bin_start", values="meth_percent")
motif_bin_methplot_df.to_csv(join(output_dir, "all_tfs_motif_bin_methplot_heatmap_data_for_dist_enh_assoc.bed"), sep ="\t", header = True, index = True)
#motif_bin_methplot_df.to_csv(join(output_dir, "all_tfs_motif_bin_methplot_heatmap_data_for_motif5.bed"), sep ="\t", header = True, index = True)
