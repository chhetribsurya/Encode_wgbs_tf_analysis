import pandas as pd
import numpy as np
from pyfasta import Fasta
import json
import os
import time
import re 

# output_dir_path = os.path.expanduser("~/Desktop/pub_data/python_output_files")
# out_directory = os.path.join(output_dir_path, "for_gata_parse")
# if not os.path.exists(out_directory):
#     os.makedirs(out_directory)


output_dir = os.path.expanduser("~/Desktop/pub_data/python_output_files/for_gata_parse")
if not os.path.exists(output_dir):
	os.makedirs(output_dir)


fasta_file = os.path.expanduser("~/Dropbox/local_miscellaneous_data/hg_19_misc/hg19-male.fa")
f = Fasta(fasta_file)
meme_file = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/MeMe/meme.txt")
query_motif = "Motif 2"

# Returns the list of motif_regex/motif_seq
def find_motif_regex(MeMe_file, motif_name=False):
	
	MeMe_file = meme_file

    with open(MeMe_file, "r") as motif_file:
        data= motif_file.read()

        motif_name_pattern = "Motif\s\d+\sregular\sexpression"
        regex = re.compile(motif_name_pattern)
        motif_name_list = regex.findall(data)

        motif_seq_pattern = "Motif\s\d+\sregular\sexpression\n-*\n(.*?\n)"
        regex_next = re.compile(motif_seq_pattern)
        motif_seq_list = regex_next.findall(data)

        motif_dict = {}

        for i in range(len(motif_name_list)):
            motif_key, motif_value = " ".join(motif_name_list[i].split()[0:2]), motif_seq_list[i].strip("\n")
         
            if motif_key in motif_dict:
                motif_dict[motif_key].append(motif_seq_list)
            else:
                motif_dict[motif_key] = [motif_value]

    #return(motif_dict[query_motif])
    return(motif_dict)

motif_patterns = find_motif_regex(meme_file)
gata_regex = motif_pattern[query_motif]



content = open(os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/MeMe/meme.txt"), "r").read()
#pattern  = r"Motif 2 sites sorted by position p-value.*\n.*\n.*\n.*\n"
pattern  = r"Motif \d+ sites sorted by position p-value.*\n.*\n.*\n.*\n"
regex = re.compile(pattern)
matched_list = regex.findall(content)
header = ["chrom","start","end","strand","meme_seq", "pyfasta_seq"]
motif_header = "\t".join(header)


Master_dict = dict()

# print matched_list
for each in matched_list:
	M_key = each[0:8].strip()
	print M_key

	motif_data = []
	for i in range(len(header)):
		motif_data.append([])

	el = len(each)
	start_index = content.find(each) + el
	end_index = content.find("-----", start_index)
	required_each_content = content[start_index: end_index]
	each_motif_lines = required_each_content.strip().split("\n")
	
	for motif_line in each_motif_lines:		
		splitted = motif_line.split()
		chrom1 = splitted[0].split(":")[0]
		sequence1 = splitted[5]
		strand1 = splitted[1]
		
		if strand1 == "+":			
			#f.sequence({"chr": "chr18", "start" : (3603155 + 22 -1), "stop" : (3603155 + 22 -1) + len_3, "strand" : "+"}, one_based = False).upper()
			# u'TCAGGTCACCAGATAAAG'
			start1 = (int(splitted[0].split(":")[1]) + int(splitted[2])) - 1
			end1 = start1 + len(sequence1)
			seq_pyfasta = f.sequence({"chr": chrom1, "start" : start1, "stop" : end1, "strand" : strand1}, one_based = False).upper()
			
		if strand1 == "-":
			#f.sequence({"chr": "chr2", "start" : (112383865 + 59 -1), "stop" : (112383865 + 59 -1) + len_3, "strand" : "-"}, one_based = False).upper()
			# u'AATCTTTGTCAGATAATC'	
			start1 = (int(splitted[0].split(":")[1]) + int(splitted[2])) -1 
			end1 = start1 + len(sequence1)
			seq_pyfasta = f.sequence({"chr": chrom1, "start" : start1, "stop" : end1, "strand" : strand1}, one_based = False).upper()
			#convert unicode string to python string
			#str(seq_pyfasta)
			#seq_pyfasta.encode("ascii", "replace")
			#seq_pyfasta.encode("ascii", "ignore")

		req = [chrom1, start1, end1, strand1, sequence1, str(seq_pyfasta)]

		for i,item in enumerate(req):
			motif_data[i].append(item)
		#for i in range(len(header)):
		# motif_data[0].append(chrom)
		# motif_data[1].append(start)

	motif_zip = zip(header,motif_data)
	motif_dict = dict(motif_zip)
	#$$$
	Master_dict[M_key] = motif_dict
	#$$$

### For quick overview of the dictionary:
df_view = pd.DataFrame(Master_dict).unstack()
motif_1_df = pd.DataFrame(Master_dict["Motif 1"]).head(10)


dict_1 = Master_dict["Motif 1"]
dict_1.keys()
chrom = dict_1["chrom"]
start = dict_1["start"]
end = dict_1["end"]
strand = dict_1["strand"]
sequence_motif_meme = dict_1["meme_seq"]
sequence_motif_pyfasta = dict_1["pyfasta_seq"]



print("Final dict_1 data starts here...\n\n")
count = 0
#for i in range(len(chrom)):
for i in range(10):
	motif_coordinates = "%s\t%s\t%s\t%s\t%s\t%s" %(chrom[i], start[i], end[i], strand[i], sequence_motif_meme[i], sequence_motif_pyfasta[i])
	print motif_coordinates
	count +=1
print "Total site :", count






### For trimming purpose:



Master_dict_1 = dict()
Trim_up = 11
Trim_down = 0

# print matched_list
for each in matched_list:
	M_key = each[0:8].strip()
	print M_key

	motif_data = []
	for i in range(len(header)):
		motif_data.append([])

	el = len(each)
	start_index = content.find(each) + el
	end_index = content.find("-----", start_index)
	required_each_content = content[start_index: end_index]
	each_motif_lines = required_each_content.strip().split("\n")
	
	for motif_line in each_motif_lines:		
		splitted = motif_line.split()
		chrom1 = splitted[0].split(":")[0]
		strand1 = splitted[1]
		sequence1 = splitted[5]
		
		Effective_len = (len(sequence1) - Trim_up)
		
		if strand1 == "+":
			# f.sequence({"chr": "chr18", "start" : (3603155 + 22 -1) + Trim_up, "stop" : (3603155 + 22 -1) + (len_3 - Trim_down), "strand" : "+"}, one_based = False).upper()
			# u'GATAAAG'
			start1_temp = (int(splitted[0].split(":")[1]) + int(splitted[2]) - 1)
			start1 = start1_temp + Trim_up
			end1 = start1_temp + (len(sequence1) - Trim_down)
			seq_pyfasta = f.sequence({"chr": chrom1, "start" : start1, "stop" : end1, "strand" : strand1}, one_based = False).upper()
			#sequence1 = splitted[5]

		if strand1 == "-":
			#f.sequence({"chr": "chr2", "start" : (112383865 + 59 -1) + Trim_down, "stop" : (112383865 + 59 -1) + (len_3 - Trim_up), "strand" : "-"}, one_based = False).upper()
			# u'GATAATC'
			start1_temp = ((int(splitted[0].split(":")[1]) + int(splitted[2])) -1)  
			start1 =start1_temp + Trim_down
			end1 = start1_temp + (len(sequence1) - Trim_up)
			seq_pyfasta = f.sequence({"chr": chrom1, "start" : start1, "stop" : end1, "strand" : strand1}, one_based = False).upper()

		req = [chrom1, start1, end1, strand1, sequence1, str(seq_pyfasta)]

		for i,item in enumerate(req):
			motif_data[i].append(item)
		#for i in range(len(header)):
		# motif_data[0].append(chrom)
		# motif_data[1].append(start)

	motif_zip = zip(header,motif_data)
	motif_dict = dict(motif_zip)
	#$$$
	Master_dict_1[M_key] = motif_dict
	#$$$


dict_2 = Master_dict_1["Motif 1"]
dict_2.keys()
chrom = dict_2["chrom"]
start = dict_2["start"]
end = dict_2["end"]
strand = dict_2["strand"]
sequence_motif_meme = dict_2["meme_seq"]
sequence_motif_pyfasta = dict_2["pyfasta_seq"]

print("Final dict_2 data_2 starts here...\n\n")
count = 0
#for i in range(len(chrom)):
for i in range(10):
	motif_coordinates = "%s\t%s\t%s\t%s\t%s\t%s" %(chrom[i], start[i], end[i], strand[i], sequence_motif_meme[i], sequence_motif_pyfasta[i])
	print motif_coordinates
	count +=1
print "Total site :", count



#Merging the dictionaries, also includes the concatenation.
dict_items = [dict_1, dict_2]

#Concatenating the dictionaries.
#pd.concat([dict_1, dict_2], ignore_index=True)

merged_dict = {}

for each_dict in dict_items:
    for key, value in each_dict.iteritems():
        if key in merged_dict:
            merged_dict[key].extend(value)
        else:
            merged_dict[key] = []
            merged_dict[key].extend(value)

df_test = pd.DataFrame(merged_dict)

chrom = merged_dict["chrom"]
start = merged_dict["start"]
end = merged_dict["end"]
strand = merged_dict["strand"]
sequence_motif_meme = merged_dict["meme_seq"]
sequence_motif_pyfasta = merged_dict["pyfasta_seq"]





### For motif analysis using fastinterval and bedtools logic, so this is zero-based so the end element is exclusive:
bin_size = 1
upstream = 50
downstream = 50
bins = range(-upstream, (downstream+1), 1)


test_data = [] 
analysis_dict = {}

#for i in range(5):
#    print analysis_dict["chr1"][i][3]

final_data = list()
final_header = ["chrom", "start", "end", "bin_start", "mid_point", "bin_end", "sequence", "motif_start", "motif_end"]

for i in range(len(final_header)):
	final_data.append(list())

coordinate_file = os.path.expanduser("~/Desktop/pub_data/python_output_files/for_gata_parse/gata_motif.bed")
with open(coordinate_file, "w") as outfile:
	print "chrom\tmid_point\tstart\tend\tbin_start\tbin_end\tsequence\tmotif_start\tmotif_end"
	for i in range(len(chrom)):
		chrome = chrom[i]
		print chrome
		motif_start = int(start[i])
		motif_end = int(end[i])
		mid_point = (motif_start + motif_end) / 2
		sequence = sequence_motif_pyfasta[i]

		for each in bins :
			bin_start = each
			bin_end = (each + bin_size)
			#For fastinterval and bedtools case:
			#bin_end = (each + bin_size)
			chrom_start = mid_point + bin_start
			chrom_end = mid_point + bin_end

			line_coords = [chrome, str(chrom_start), str(chrom_end), str(bin_start), str(mid_point), str(bin_end), sequence, str(motif_start), str(motif_end)]
			test_data.append("\t".join(line_coords))

			if chrome in analysis_dict:
				analysis_dict[chrome].append((chrome, int(chrom_start), int(chrom_end), "\t".join(line_coords)))
			else:
				analysis_dict[chrome] = [(chrome, int(chrom_start), int(chrom_end), "\t".join(line_coords))]

				
			for i,item in enumerate(line_coords):
				final_data[i].append(item)

			print_coords = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(chrome, chrom_start, chrom_end, bin_start, mid_point, bin_end, sequence, motif_start, motif_end) 
			print print_coords
			outfile.write(print_coords + "\n")

	final_zip = zip(final_header, final_data)
	final_dict = dict(final_zip)




###
#awk '{print $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t"$6}' SL60583_GATA2_rep1.cov > SL60583_GATA2_rep1_0_based.cov
import pybedtools

#Cpg_file = os.path.expanduser("~/surya/for_gata_parse/head_cpg.bed")
Cpg_file = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/for_gata_parse/SL60583_GATA2_rep1_0_based_1.cov")
coordinate_file = os.path.expanduser("~/Desktop/pub_data/python_output_files/for_gata_parse/gata_motif.bed")
c = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/for_gata_parse/fastinterval_overlaps.bed")

cpg_bed = pybedtools.example_bedtool(Cpg_file)
gata_bed = pybedtools.example_bedtool(coordinate_file)
c = pybedtools.example_bedtool(c)
a_with_b = cpg_bed.intersect(gata_bed, wa = True, wb = True)
# a_with_b = a.intersect(a, u = True)
# a_with_b = cpg_bed.intersect(gata_bed, wo= True)
# a_with_b.groupby(g=[1, 4, 5], c=10, o =[sum])
a_with_b_grouped = a_with_b.groupby(g=[11, 15, 16], c=[5,7])
c = a_with_b_grouped.saveas(os.path.expanduser('~/Desktop/pub_data/python_output_files/for_gata_parse/grouped_gata_motif.bed'))


a_with_b.groupby(g=[11, 15, 16], c=[5,7]).head()
a_with_b.head()

a_with_b.groupby(g=[11, 15, 16], c=[5,7]).tail()
a_with_b.tail()
# df = pandas.read_table(cpg_bed.fn, names=['chrom', 'start', 'stop', 'perc', 'meth', 'strand', 'unmeth'])
# print len(a_with_b)
# print(a_with_b)
# print(open(a_with_b.fn).read())
d= a_with_b.groupby(g=[1, 4, 5], c=[5,10], o=["sum","max"])

# cpg_intersect = a.intersect( 
# 								[b.fn, c.fn],
# 								names = ["b", "c"],
# 								filenames = True,
# 								wa = True,
# 								u = True

# 							)
# cpg_intersect.head()



# line_splitted = open(a_with_b.fn).read().split("\n")
# for each in line_splitted:
#     tab_splitted = each.split()
#     chrom = tab_splitted[0]
#     start = tab_splitted[1]
#     stop = tab_splitted[2]

# count = 0
# for each in line_splitted:
#     if each is not "" or None:
#         print each.split()
#         count +=1
# print count

# c= a_with_b
# d= c.groupby(g=[1, 4, 5], c=10, ops=['sum'])
# print d

# c = a_with_b.saveas('/home/surya/surya/for_gata_parse/pybed_overlap.bed', trackline='track name="a and b"')
# c = a_with_b_grouped.saveas('/home/surya/Desktop/pub_data/python_output_files/for_gata_parse/grouped_gata_motif.bed')