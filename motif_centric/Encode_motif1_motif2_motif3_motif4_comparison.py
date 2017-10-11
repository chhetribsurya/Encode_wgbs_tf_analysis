from os.path import basename
from os.path import join
from glob import glob
import pandas as pd



motif1_tf = []

for idx, line in enumerate(open("motif1_file_list", "r")): 
	tf_name = re.compile("final_wgbs_motifs_intersect_2kb.bed_(.*).bed").findall(basename(line))[0]   
	print idx, tf_name 
	motif1_tf.append(tf_name)

df_1 = pd.Series(motif1_tf)
df_1[~df_1.isin(motif2_tf)] 



motif2_tf = []

for idx, line in enumerate(open("motif2_file_list", "r")): 
	tf_name = re.compile("final_wgbs_motifs_intersect_2kb.bed_(.*).bed").findall(basename(line))[0]   
	print idx, tf_name 
	motif2_tf.append(tf_name)


motif3_tf = []

for idx, line in enumerate(open("motif3_file_list", "r")): 
	tf_name = re.compile("final_wgbs_motifs_intersect_2kb.bed_(.*).bed").findall(basename(line))[0] 
	print idx, tf_name 
	motif3_tf.append(tf_name)



motif4_tf = []

for idx, line in enumerate(open("motif4_file_list", "r")): 
	tf_name = re.compile("final_wgbs_motifs_intersect_2kb.bed_(.*).bed").findall(basename(line))[0]   
	print idx, tf_name 
	motif4_tf.append(tf_name)