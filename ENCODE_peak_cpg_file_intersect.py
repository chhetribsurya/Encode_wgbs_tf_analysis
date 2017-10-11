import pybedtools
import os
import re
import glob

#Cpg_file reading via pandas_dataframe for adding the strands for standard pybed usage, and editing the coordinates:
def load_stringf_for_cpg_after_edits(input_file):

	df = pd.read_csv(input_file, sep= "\t", header= None)
	df.columns = ["chrom", "start", "end", "perc_meth", "meth", "unmeth"]
	df.loc[:,"end"] = df["end"] + 1
	df.loc[:,"strand"] = "."
	Cpg_file = df.loc[:,["chrom","start", "end", "meth", "unmeth", "strand"]]

	Cpg_list = Cpg_file.values.tolist()
	Cpg_string_list = [ "\t".join(map(str, line)) for line in Cpg_list]
	Cpg_string_bed_format = "\n".join(Cpg_string_list)

	Cpg_bed_file = pybedtools.BedTool(Cpg_string_bed_format,from_string=True)

	return(Cpg_bed_file)


#Cpg_file reading via pandas_dataframe for adding the strands for standard pybed usage, and editing the coordinates:
def load_pybedObject_for_cpg_after_edits(input_file):

	file_basename = os.path.splitext(basename(input_file))[0]
	file_ext = os.path.splitext(basename(input_file))[1]

	df = pd.read_csv(input_file, sep= "\t", header= None)
	df.columns = ["chrom", "start", "end", "perc_meth", "meth", "unmeth"]
	df.loc[:,"end"] = df["end"] + 1
	df.loc[:,"strand"] = "."
	Cpg_file = df.loc[:,["chrom","start", "end", "meth", "unmeth", "strand"]]
	
	file_name = output_dir + "/" + file_basename + "_" + "cpg_edit" + file_ext
	#Make sure that header=False else, pybedtools won't recognise the headers above:
	Cpg_file.to_csv(file_name, sep ="\t", header = False, index = False)
	Cpg_bed_file = pybedtools.BedTool(file_name)

	#Cpg_list = Cpg_bed_filter.values.tolist()
	# Cpg_list = Cpg_file.values.tolist()
	# Cpg_string_list = [ "\t".join(map(str, line)) for line in Cpg_list]
	# Cpg_string_bed_format = "\n".join(Cpg_string_list)
	#Cpg_bed_file = pybedtools.BedTool(Cpg_string_bed_format,from_string=True)

	return(Cpg_bed_file)


import pybedtools
import glob
import os
from os.path import basename
from os.path import join
from os.path import splitext


def meth_percent_calc(x):
    x_meth = x["cpg_meth"].sum()
    x_unmeth = x["cpg_meth"].sum() + x["cpg_unmeth"].sum()
    x_perc = x_meth/float(x_unmeth)
    return(x_perc)

#cpg_bed_list = glob.glob("SL*.cov")
peak_bed_list = glob.glob("SL*narrowPeak")
cpg_bed_file = " "
output_dir = os.path.expand("~/for_chris/batch_I/wgbs_tf_intersect")

if not os.path.exists(output_dir):
    os.makedirs(output_dir)


""" \nProcessing and creating a Cpg bed file with a strand cols from raw cpg file 
followed by pybedtool object formation...\n """
os.environ["output_dir"] = output_dir
os.environ["cpg_bed_file"] = cpg_bed_file
cpg_bed_edited = splitext(basename(cpg_bed_file))[0] + "_edited" + splitext(cpg_bed_file)[1]
os.environ["cpg_bed_edited "] = cpg_bed_edited
COMMAND = '''awk 'BEGIN{FS=" "; OFS="\t"} { print $1,$2,$3+1,$5,$6,"." }' $cpg_bed_file | sort -k1,1 -k2,2n  > $output_dir/$cpg_bed_edited'''
os.system(COMMAND) #subprocess.call(COMMAND, shell=True) # import subprocess
Cpg_bed_file = pybedtools.BedTool(cpg_bed_edited)



#df_intersect = Cpg_bed_file.intersect(peak_bed_file, wa = True, wb = True, output="intersect_cpg_bed.txt")
master_tf_dict = {}

for each_file in peak_bed_list:
    """ \nProcessing and creating a sorted bed file for the cpg intersect
    followed by pybedtools object formation...\n """
    os.environ["each_file"] = each_file
    sorted_peak_file = splitext(basename(each_file))[0] + "_sorted" + splitext(basename(each_file))[1]
    os.environ["sorted_peak_file"] = sorted_peak_file    
	  os.system('sort -k1,1, -k2,2n $each_file > $output_dir/$sorted_peak_file')
    TF_name = basename(each_file).split("_")[-1] 
    peak_bed = pybedtools.BedTool(each_file)

    """ \nPybedtool intersection of Cpg_bed_file and peak_bed_file...\n """
    pybed_outfile = join(output_dir, TF_name + "_cpg_intersect.bed")
    Cpg_bed_file.intersect(peak_bed, wa = True, wb = True, output=pybed_outfile)

    ### Working with the dataframes; reading the output file of pybedtool intersect:
    df_file = pd.read_csv(pybed_outfile, sep = "\t", header = None)
    df_ordered = df_file.iloc[:, [6, 7, 8, 0, 1, 2, 3, 4]]
    df_ordered.columns = ["peak_chrom", "peak_start", "peak_end", "cpg_chrom", "cpg_start", "cpg_end", "cpg_meth", "cpg_unmeth" ]
    df_ordered[["cpg_meth", "cpg_unmeth"]] = df_ordered[["cpg_meth", "cpg_unmeth"]].astype(int)
    df_grouped =  df_ordered.groupby(["peak_chrom", "peak_start", "peak_end"]).apply(lambda x : x["cpg_meth"].sum()/float(x["cpg_meth"].sum() + x["cpg_unmeth"].sum()))  
    #peak_grouped = df_ordered.groupby(["peak_chrom", "peak_start", "peak_end"]).apply(lambda x: x.sum())  ### For quick test of the grouped items:
    #peak_grouped = df_ordered.groupby(["peak_chrom", "peak_start", "peak_end"]).apply(meth_percent_calc)  

    if TF_name not in tf_dict:
      master_tf_dict[TF_name] = df_grouped.reset_index()

### Combine all the dataframes (i.e all values of master_tf_dict) representing the tf intersects with CpG:
tf_combined_df = pd.concat(master_dict).reset_index()
tf_combined_df.to_csv(join(output_dir, "wgbs_tf_intersect.bed", sep ="\t", header = True, index = False)



Dimension of current TF is (7931,)

Intersection of FOXO1 completed!!!...


processing /gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks/test_analysis/SL151607.filt.nodup.srt.SE_VS_SL151608.filt.nodup.srt.SE.IDR0.02.filt.narrowPeak_FOXA3[FLAG]
 Processing the intersection for FOXA3[FLAG]
Dimension of current TF is (35827,)

Intersection of FOXA3[FLAG] completed!!!...


processing /gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks/test_analysis/SL88833.filt.nodup.srt.SE_VS_SL88834.filt.nodup.srt.SE.IDR0.02.filt.narrowPeak_RAD21[FLAG]
 Processing the intersection for RAD21[FLAG]
Dimension of current TF is (50486,)

Intersection of RAD21[FLAG] completed!!!...


 #    for key in df_grouped.groups:
	# 	  grouped_data = Cpg_grouped.get_group(key)
	# 	  print grouped_data


	# f = {"cpg_meth": np.sum, "cpg_unmeth": np.sum}
	# peak_agg = peak_grouped.agg(f)
	# peak_agg["Percent_meth"] = (peak_agg["cpg_meth"]) / (peak_agg["cpg_meth"] + peak_agg["cpg_unmeth"] )
	# #peak_agg["group_name"] = ggplot_group_name

	# peak_agg_1 = peak_agg.reset_index()
	# peak_agg_1["bin_start"] = peak_agg_1["bin_start"].astype(int)
	# peak_agg_sorted = peak_agg_1.sort_values(["bin_start"], ascending = [1])
	# #peak_agg_sorted.index = peak_agg_sorted.index.sort_values()

	# peak_agg_sorted.to_csv(directory + "/" + ggplot_outfile_suffix, sep ="\t", header = True, index = False)








df = pd.DataFrame({"name":["Foo", "Baar", "Foo", "Baar"], "score_1":[5,10,15,10], "score_2" :[10,15,10,25], "score_3" : [10,20,30,40]})


In [61]: df
Out[61]:
   name  score_1  score_2  score_3
0   Foo        5       10       10
1  Baar       10       15       20
2   Foo       15       10       30
3  Baar       10       25       40


In [52]: def f(x):
        x_meth = x["score_2"].sum()
        x_unmeth = x["score_2"].sum() + x["score_3"].sum()
        x_perc = x_meth/float(x_unmeth)
        return(x_perc)
   ....:



In [54]: df.groupby(["name", "score_1"]).apply(f)
Out[54]:
name  score_1
Baar  10         0.40
Foo   5          0.50
      15         0.25
dtype: float64


In [53]: df.groupby(["name", "score_1"]).apply(lambda x : x["score_2"].sum()/float(x["score_2"].sum() + x["score_3"].sum()))
Out[53]:
name  score_1
Baar  10         0.40
Foo   5          0.50
      15         0.25
dtype: float64

In [55]: df.groupby(["name", "score_1"]).apply(lambda x : x["score_2"].sum()/float(x["score_2"].sum() + 10))
Out[55]:
name  score_1
Baar  10         0.8
Foo   5          0.5
      15         0.5
dtype: float64

In [56]: df.groupby(["name", "score_1"]).apply(lambda x : x.sum())
Out[56]:
                  name  score_1  score_2  score_3
name score_1
Baar 10       BaarBaar       20       40       60
Foo  5             Foo        5       10       10
     15            Foo       15       10       30



# In [60]: df.groupby(["name", "score_1"])["score_2"].sum().apply(lambda x : x.sum())
# Out[60]:
# name  score_1
# Baar  10         40
# Foo   5          10
#       15         10
# Name: score_2, dtype: int64


In [74]: df.groupby(["name", "score_1"]).apply(lambda x : x["score_2"].sum())
Out[74]:
name  score_1
Baar  10         40
Foo   5          10
      15         10
dtype: int64

In [69]: df_test.reset_index()
Out[69]:
   name  score_1         0
0  Baar       10  0.666667
1   Foo        5  0.333333
2   Foo       15  0.333333

In [70]: df_set = df_test.reset_index()

In [71]: df_set.columns = ["start", "end", "perc"]

In [72]: df_set
Out[72]:
  start  end      perc
0  Baar   10  0.666667
1   Foo    5  0.333333
2   Foo   15  0.333333





In [88]: def f(x):
        x_meth = x["score_2"].sum()
        x_unmeth = x["score_2"].sum() + x["score_3"].sum()
        x_perc = x_meth/float(x_unmeth)
        return(x_perc,x_unmeth, x_meth)
   ....:

In [90]: df_test = df.groupby(["name", "score_1"]).apply(f)

In [91]: df_test
Out[91]:
name  score_1
Baar  10         (0.4, 100, 40)
Foo   5           (0.5, 20, 10)
      15         (0.25, 40, 10)
dtype: object



In [93]: df_test.reset_index()
Out[93]:
   name  score_1               0
0  Baar       10  (0.4, 100, 40)
1   Foo        5   (0.5, 20, 10)
2   Foo       15  (0.25, 40, 10)

In [94]: df_reset= df_test.reset_index()


### name of the column is 0 itself, so don't confuse it as a indexing 0:
In [120]: df_reset.columns
Out[120]: Index([u'name', u'score_1', 0], dtype='object')


### name of the column is 0 itself, so don't confuse it as a indexing 0:
In [95]: df_reset[0]
Out[95]:
0    (0.4, 100, 40)
1     (0.5, 20, 10)
2    (0.25, 40, 10)
Name: 0, dtype: object



In [100]: df_reset["first"] = df_reset[0].apply( lambda x: x[0])

In [101]: df_reset
Out[101]:
   name  score_1               0  first
0  Baar       10  (0.4, 100, 40)   0.40
1   Foo        5   (0.5, 20, 10)   0.50
2   Foo       15  (0.25, 40, 10)   0.25

In [102]: df_reset["second"] = df_reset[0].apply( lambda x: x[1])

In [103]: df_reset["third"] = df_reset[0].apply( lambda x: x[2])

In [104]: df_reset
Out[104]:
   name  score_1               0  first  second  third
0  Baar       10  (0.4, 100, 40)   0.40     100     40
1   Foo        5   (0.5, 20, 10)   0.50      20     10
2   Foo       15  (0.25, 40, 10)   0.25      40     10


### Drop the column or delete the column, 1 for col, and 0 for rows:
In [108]: df_reset.drop(0, 1)
Out[108]:
   name  score_1  first  second  third
0  Baar       10   0.40     100     40
1   Foo        5   0.50      20     10
2   Foo       15   0.25      40     10

In [109]: df_reset.drop(0, 1, inplace=True)

In [110]: df_reset
Out[110]:
   name  score_1  first  second  third
0  Baar       10   0.40     100     40
1   Foo        5   0.50      20     10
2   Foo       15   0.25      40     10



In [112]: df_reset.drop("first", axis=1 )
Out[112]:
   name  score_1  second  third
0  Baar       10     100     40
1   Foo        5      20     10
2   Foo       15      40     10

In [113]: df_reset
Out[113]:
   name  score_1  first  second  third
0  Baar       10   0.40     100     40
1   Foo        5   0.50      20     10
2   Foo       15   0.25      40     10

In [114]: df_reset.drop("first", axis=1, inplace = True )

In [115]: df_reset
Out[115]:
   name  score_1  second  third
0  Baar       10     100     40
1   Foo        5      20     10
2   Foo       15      40     10



In [781]: pd.concat(master_dict)
Out[781]: 
       name  score_1  score_2         0
fox 0  Baar       10       15  0.428571
    1  Baar       10       25  0.384615
    2   Foo        5       10  0.500000
    3   Foo       15       10  0.250000
pol 0  Baar       10       15  0.428571
    1  Baar       10       25  0.384615
    2   Foo        5       10  0.500000
    3   Foo       15       10  0.250000
sp1 0  Baar       10       15  0.428571
    1  Baar       10       25  0.384615
    2   Foo        5       10  0.500000
    3   Foo       15       10  0.250000
sp5 0  Baar       10       15  0.428571
    1  Baar       10       25  0.384615
    2   Foo        5       10  0.500000
    3   Foo       15       10  0.250000

In [782]: pd.concat(master_dict).reset_index()
Out[782]: 
   level_0  level_1  name  score_1  score_2         0
0      fox        0  Baar       10       15  0.428571
1      fox        1  Baar       10       25  0.384615
2      fox        2   Foo        5       10  0.500000
3      fox        3   Foo       15       10  0.250000
4      pol        0  Baar       10       15  0.428571
5      pol        1  Baar       10       25  0.384615
6      pol        2   Foo        5       10  0.500000
7      pol        3   Foo       15       10  0.250000
8      sp1        0  Baar       10       15  0.428571
9      sp1        1  Baar       10       25  0.384615
10     sp1        2   Foo        5       10  0.500000
11     sp1        3   Foo       15       10  0.250000
12     sp5        0  Baar       10       15  0.428571
13     sp5        1  Baar       10       25  0.384615
14     sp5        2   Foo        5       10  0.500000
15     sp5        3   Foo       15       10  0.250000

In [783]: df_final = pd.concat(master_dict).reset_index()

In [784]: df_final
Out[784]: 
   level_0  level_1  name  score_1  score_2         0
0      fox        0  Baar       10       15  0.428571
1      fox        1  Baar       10       25  0.384615
2      fox        2   Foo        5       10  0.500000
3      fox        3   Foo       15       10  0.250000
4      pol        0  Baar       10       15  0.428571
5      pol        1  Baar       10       25  0.384615
6      pol        2   Foo        5       10  0.500000
7      pol        3   Foo       15       10  0.250000
8      sp1        0  Baar       10       15  0.428571
9      sp1        1  Baar       10       25  0.384615
10     sp1        2   Foo        5       10  0.500000
11     sp1        3   Foo       15       10  0.250000
12     sp5        0  Baar       10       15  0.428571
13     sp5        1  Baar       10       25  0.384615
14     sp5        2   Foo        5       10  0.500000
15     sp5        3   Foo       15       10  0.250000

In [785]: df_final["level_0"]
Out[785]: 
0     fox
1     fox
2     fox
3     fox
4     pol
5     pol
6     pol
7     pol
8     sp1
9     sp1
10    sp1
11    sp1
12    sp5
13    sp5
14    sp5
15    sp5
Name: level_0, dtype: object

In [786]: type(df_final)
Out[786]: pandas.core.frame.DataFrame





