import re
import pandas as pd
from glob import glob
import os
from os.path import join
from os.path import splitext
from os.path import basename

file_dir = "/gpfs/gpfs1/home/schhetri/encode_hepg2_data_analysis/TF_fimo_motifs"
output_dir = "/gpfs/gpfs1/home/schhetri/encode_hepg2_data_analysis/avg_gc_content_cg_count_for_motifs"
motif_file_list = glob(join(file_dir, "*fimo*.txt"))

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

master_dict = {}
for each_file in motif_file_list: 
    file_name = "Avg_GC_content_CG_count_for_" + basename(each_file)
    tf_name = re.compile("fimo_motif_*(.*).txt").findall(basename(each_file))[0]
    print "processing %s"%(tf_name)
    df = pd.read_csv(each_file, sep="\t")
    ### (count("G") + count("C"))/(len(sequence)) i.e "G" or "C" counts per total length
    df["GC_content"] = df["matched sequence"].str.count("[GC]", re.IGNORECASE)/df["matched sequence"].str.len()
    df["motif_length"] = df["matched sequence"].str.len()
    df["CpG_count"] = df["matched sequence"].str.count("CG", re.IGNORECASE)*2
    df_grouped = df.groupby("#pattern name")[["GC_content","motif_length", "CpG_count"]].apply(lambda x : x.mean()).reset_index()
    df_grouped.to_csv(join(output_dir, file_name), sep="\t", header=True, index=False)

    if tf_name not in master_dict:
        master_dict[tf_name] = df_grouped

df_combined = pd.concat(master_dict).reset_index().drop("level_1",1)
df_combined.rename(columns={"level_0" : "tf_name"}, inplace=True)
df_combined.to_csv(join(output_dir, "all_tfs_combined_for_avg_GC_content_CG_count.txt"), sep="\t", header=True, index=False)

print "Task completed!!!"


