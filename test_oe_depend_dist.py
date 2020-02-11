from __future__ import print_function
import numpy as np
import pandas as pd
import seaborn
import matplotlib.pyplot as plt
import math
source_path = "D:/Lab Archive/Projects Archive/3Dpredictor/source"
import sys
sys.path.append(source_path)
from shared import get_bin_size
from clustering2compartments import get_compartments


fname = "AcolNg_V3.1000.hic.100000.oe"
data = pd.read_csv(
    fname,
    sep="\t",
    header=None,
    names=["chr","st","end","oe"],
    dtype={"chr":str,"st":np.uint64, "end":np.uint64, "oe":np.float32}
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

data["dist"] = ((data["end"] - data["st"]) / binsize).astype("int")
data["log_oe"] = data.oe.apply(math.log)

#groupped = data.groupby(by="chr")
step = 5

for chr in np.unique(data.chr.values):
    curr_chr = data.query("chr==@chr")
    #groupped = curr_chr.groupby("dist")
    distances = np.unique(curr_chr.dist.values)
    subset_size = int(len(distances) / 20)
    print (len(distances))
    subset = curr_chr.query("dist % @subset_size == 0")

    averages = subset.groupby("dist").aggregate(np.average)
    print (range(len(averages.index)),
           averages.log_oe)
    seaborn.set_style("whitegrid")
    v = seaborn.violinplot(x = subset.dist, y=subset.log_oe)
    v.set_xticklabels(v.get_xticklabels(), rotation=45)
    plt.axhline(y=0)
    plt.title(fname + chr)
    plt.title(fname + chr)
    plt.xlabel("Distance * "+str(binsize))
    plt.savefig(fname+chr+"oeViolin.png",dpi=300)
    plt.clf()