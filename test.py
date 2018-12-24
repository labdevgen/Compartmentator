from __future__ import print_function
import numpy as np
import pandas as pd
from shared import get_bin_size
from clustering2compartments import get_compartments

fname = "Aste_VARYA_comb_allValidPairs.hic.10000.oe.10MB"
data = pd.read_csv(
    fname,
    sep="\t",
    header=None,
    names=["chr","st","end","oe"],
    dtype={"chr":str,"st":np.uint32, "end":np.uint32, "oe":np.float32}
)

print(data.head())
assert np.all((data["end"]-data["st"]).values>=0)

#print (data.max(),data.min())
null_values_count = data.isnull().values.sum()
if null_values_count != 0:
    print ("Found null values")
    print (null_values_count," out of ",len(data))
    data = data.dropna()

assert not data.isnull().values.any()

assert not data.isna().values.any()
binsize = get_bin_size(data,fields=["st","end"])

data["end"] = data["end"] // binsize
data["st"] = data["st"] // binsize

groupped = data.groupby(by="chr")
boundaries = []

length = 1500000 // binsize
step = length // 5
# end = min(start+length+step*10,len(array) - length)

for chr,groupped_data in groupped:
    start = 0
    print ("Computing chr ",chr)
    result = get_compartments(groupped_data,binsize, start,length,step)
    #print (result.tail(5))
    boundaries.append(result)
    #break
boundaries = pd.concat(boundaries)
boundaries[["chr", "st", "end", "Strength"]].to_csv(
    fname + ".win" + str((length * binsize) // 1000) + "kb.boundaries.bedGraph",
    sep="\t", index=False, header=None)
