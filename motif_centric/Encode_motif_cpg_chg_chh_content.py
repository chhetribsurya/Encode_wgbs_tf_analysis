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

output_dir = os.path.expanduser("~/Dropbox/local_miscellaneous_data/motifs_methplot")
if not os.path.exists(output_dir):
	os.makedirs(output_dir)

# motif_output_file = "final_motifs_coordinate.bed"
final_output_file = "final_wgbs_motifs_intersect.bed"


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


def combine_motif_coords_after_parsing_meme_motifs(master_motif_dict, **kwargs):
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
	final_combined_df.to_csv(join(output_dir, "motif_coordinate_info.bed"), sep="\t", index=False, header=False)
	print "\nTotal count of motifs combined:", count
	return(final_combined_df)

final_motif_df = combine_motif_coords_after_parsing_meme_motifs(Master_motif_dict, motif_1= 1, motif_2= 2)


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


def combine_trimmed_motif_df(*args):
	combine_df_list = []
	for each_df in args:
		combine_df_list.append(each_df)

	combined_trimmed_df = pd.concat(combine_df_list, ignore_index=True)
	final_trimmed_motif_df = combined_trimmed_df.iloc[:,[2,6,3,7,4,5,0]]
	final_trimmed_motif_df.to_csv(join(output_dir, "trim_based_motif_coordinate_info.bed"), sep="\t", index=False, header=False)
	return(final_trimmed_motif_df)

final_trimmed_motif_df = combine_trimmed_motif_df(df1,df2)

### Using boolean style of slicing:
# final_trimmed_motif_df[final_trimmed_motif_df["pyfasta_seq"].str.contains("CG")]
# final_trimmed_motif_df[final_trimmed_motif_df["pyfasta_seq"].str.contains("(CG)")]

### For regex compatible:
# [m.span() for m in re.compile("GA").finditer(str_test)]
# for each in [m.span() for m in re.compile("GA").finditer(str_test)]:
# 	print each[0], each[1]
# for each_str in final_trimmed_motif_df["pyfasta_seq"].tolist():
# 	[m.span() for m in re.compile("GA").finditer(each_str)]

# df_match = final_trimmed_motif_df[~final_trimmed_motif_df["pyfasta_seq"].str.extract("(CG)").isnull()]
# df_match = final_trimmed_motif_df["pyfasta_seq"].str.extractall("(C([CAT]))")
# df_match = final_trimmed_motif_df["pyfasta_seq"].str.extractall("(C([\w+]G))").head(60)
# df_match = final_trimmed_motif_df["pyfasta_seq"].str.contains("CG").astype(int).head(60).sum()
# df_match.reset_index().groupby(["level_0"]).apply( lambda x: x["match"].max())
# df_match = final_trimmed_motif_df["pyfasta_seq"].str.extractall("(C([A-Z]))")
# select_indices = df_match.index.get_level_values(0).tolist() ## includes repetition of index
# final_trimmed_motif_df["pyfasta_seq"].str.extractall("(C(G))")
# select_indices = CG_match_df.index.levels[0].tolist() ## includes unique indexes
# motif_seq_df = final_motif_df.iloc[select_indices]


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
    final_CG_motif_joined_df.to_csv(join(output_dir, tf_name+"CG_containing_motif_coordinate_info.bed"), sep="\t", index=False, header=False)

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


final_CG_containing_motif_df, final_CG_motif_summary_dict = CG_containing_motifs_coord_df(final_motif_df)
# CG_containing_motif_df["meme_seq"].str.replace("CG","0")


"""CHG (where H is A, C or T)"""
def CHG_containing_motifs_coord_df(final_motif_df_info, tf_name):

    final_motif_df = final_motif_df_info
    ### For finding CHG containing motif dataframes, ussing boolean style of slicing:
    CHG_containing_motif_df = final_motif_df[final_motif_df["meme_seq"].str.contains("C[ATC]G")]
    C_lacking_motif_df = final_motif_df[~final_motif_df["meme_seq"].str.contains("C")]

    #CHG_regex_match_df = final_motif_df["meme_seq"].str.extractall("(C([ATC]G))")
    CHG_regex_match_df = final_motif_df["meme_seq"].str.extractall("((C[ATC]G))")
    CHG_regex_grouped_df = CHG_regex_match_df.reset_index().groupby(["level_0"]).apply( lambda x: x["match"].max())

    final_CHG_motif_joined_df = pd.concat([CHG_containing_motif_df, CHG_regex_grouped_df], axis=1)
    final_CHG_motif_joined_df = final_CHG_motif_joined_df.rename(columns={0:"CHG_count"}) ## since 0 based indices so add 1 to count a match:
    final_CHG_motif_joined_df["CHG_count"] = final_CHG_motif_joined_df["CHG_count"]+1 ## since 0 based indices so add 1 to count a match:
    final_CHG_motif_joined_df.to_csv(join(output_dir, tf_name+"CHG_containing_motif_coordinate_info.bed"), sep="\t", index=False, header=False)

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

final_CHG_containing_motif_df, final_CHG_motif_summary_dict = CHG_containing_motifs_coord_df(final_motif_df)
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
    final_CHH_motif_joined_df.to_csv(join(output_dir, tf_name+"CHH_containing_motif_coordinate_info.bed"), sep="\t", index=False, header=False)
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


final_CHH_containing_motif_df, final_CHH_motif_summary_dict = CHH_containing_motifs_coord_df(final_motif_df)
# CHH_containing_motif_df["meme_seq"].str.replace("C[ATC][ATC]","0")


meme_file_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/chip_analysis/IDR*/MeMe/*.txt"
meme_file_list = glob.glob(meme_file_path)

master_tf_dict = {}
for meme_file in meme_file_list:
        #meme_file = '/gpfs/gpfs1/home/schhetri/for_chris/batch_I/chip_analysis/IDR_SL151597_SL151598_KLF6_v2[FLAG]/MeMe/meme.txt' #string instance for regex
        tf_name = re.compile(r".*IDR_(.*?)_(.*?)_(.*?)(?=\/MeMe)").findall(meme_file)[0][-1] #[('SL151597', 'SL151598', 'KLF6_v2[FLAG]')]
        SL_rep1 = re.compile(r".*IDR_(.*?)_(.*?)_(.*?)(?=\/MeMe)").findall(meme_file)[0][0]
        SL_rep2 = re.compile(r".*IDR_(.*?)_(.*?)_(.*?)(?=\/MeMe)").findall(meme_file)[0][1]
        print "\nCurrently processing %s tf\n" %(tf_name)
        #motif_patterns = parse_meme_motif_regex(meme_file)
        Master_motif_dict = parse_meme_motif_coords(meme_file)
        final_motif_df = combine_motif_coords_after_parsing_meme_motifs(Master_motif_dict, motif_1= 1, motif_2= 2)
        #df1 = trim_and_parse_single_meme_motif_coords(meme_file, 1, 11, 0)
        #df2 = trim_and_parse_single_meme_motif_coords(meme_file, 2, 0, 0)
        #final_trimmed_motif_df = combine_trimmed_motif_df(df1,df2)

        final_CG_containing_motif_df, final_CG_motif_summary_dict = CG_containing_motifs_coord_df(final_motif_df, tf_name)
        final_CHG_containing_motif_df, final_CHG_motif_summary_dict = CHG_containing_motifs_coord_df(final_motif_df, tf_name)
        final_CHH_containing_motif_df, final_CHH_motif_summary_dict = CHH_containing_motifs_coord_df(final_motif_df, tf_name)

        if tf_name not in master_tf_dict:
                master_tf_dict[tf_name] = final_CG_motif_summary_dict
                master_tf_dict[tf_name].update(final_CHG_motif_summary_dict)
                master_tf_dict[tf_name].update(final_CHH_motif_summary_dict)

master_tf_df = pd.DataFrame(master_tf_dict)
final_master_tf_df = master_tf_df.T
final_master_tf_df.to_csv(join(output_dir, final_output_file), sep ="\t", header = True, index = True)

print "Task completed!!!"
print "\n\nTime taken to complete the analayis", time.time()-start_time

final_motif_df["fimo_seq"].str.extractall("(([a-zA-Z]+C))")

	                  0        1
     match                  
1    0      AGTAAAC  AGTAAAC
2    0      AGTAAAC  AGTAAAC
3    0      AGTAAAC  AGTAAAC
4    0      AGTAAAC  AGTAAAC
5    0      AGTAAAC  AGTAAAC
6    0      AGTAAAC  AGTAAAC
7    0      AGTAAAC  AGTAAAC
8    0      AGTAAAC  AGTAAAC
9    0      AGTAAAC  AGTAAAC
10   0      agtaaaC  agtaaaC
11   0      AGTAAAC  AGTAAAC
14   0      AGTAAAC  AGTAAAC
15   0      AGTAAAC  AGTAAAC
18   0      AGTAAAC  AGTAAAC
19   0      agtaaaC  agtaaaC
20   0      AGTAAAC  AGTAAAC
21   0      AGTAAAC  AGTAAAC
23   0      AGTAAAC  AGTAAAC
24   0      AGTAAAC  AGTAAAC
25   0      AGTAAAC  AGTAAAC
26   0      AGTAAAC  AGTAAAC
27   0      AGTAAAC  AGTAAAC
28   0      AGTAAAC  AGTAAAC
30   0      AGTAAAC  AGTAAAC
31   0      AGTAAAC  AGTAAAC
33   0      AGTAAAC  AGTAAAC
34   0      AGTAAAC  AGTAAAC
35   0      AGTAAAC  AGTAAAC
36   0      AGTAAAC  AGTAAAC
37   0      AGTAAAC  AGTAAAC
...             ...      ...
7443 0      GTGAGTC  GTGAGTC
7445 0      GTGAGTC  GTGAGTC
7446 0      GTGAGTC  GTGAGTC
7447 0      GTGAGTC  GTGAGTC
7448 0      GTGAGTC  GTGAGTC
7449 0      GTGAGTC  GTGAGTC
7451 0      GTGAGTC  GTGAGTC
7452 0      GTGAGTC  GTGAGTC
7453 0      GTGAGTC  GTGAGTC
7455 0      GTGAGTC  GTGAGTC
7456 0      GTGAGTC  GTGAGTC
7458 0      GTGAGTC  GTGAGTC
7460 0      GTGAGTC  GTGAGTC
7461 0      GTGAGTC  GTGAGTC
7463 0      GTGAGTC  GTGAGTC
7464 0      GTGAGTC  GTGAGTC
7465 0      GTGAGTC  GTGAGTC
7466 0      GTGAGTC  GTGAGTC
7471 0      GTGAGTC  GTGAGTC
7473 0      GTGAGTC  GTGAGTC
7474 0      GTGAGTC  GTGAGTC
7475 0      GTGAGTC  GTGAGTC
7476 0      GTGAGTC  GTGAGTC
7478 0      GTGAGTC  GTGAGTC
7479 0      GTGAGTC  GTGAGTC
7480 0      GTGAGTC  GTGAGTC
7481 0      GTGAGTC  GTGAGTC
7482 0      GTGAGTC  GTGAGTC
7483 0      GTGAGTC  GTGAGTC
7484 0      GTGAGTC  GTGAGTC

[5385 rows x 2 columns]


