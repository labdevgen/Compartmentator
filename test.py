from __future__ import print_function
import numpy as np
import pandas as pd
from shared import get_bin_size
from sklearn.cluster import AgglomerativeClustering
import datetime

fname = "Aste_VARYA_comb_allValidPairs.2L.oe.10MB"
data = pd.read_csv(
    fname,
    sep="\t",
    header=None,
    names=["st","end","oe"],
    dtype={"st":np.uint32, "end":np.uint32, "oe":np.float32}
)

print(data.head())
assert np.all((data["end"]-data["st"]).values>=0)
data = data.dropna()
print (data.max(),data.min())
assert not data.isnull().values.any()
assert not data.isna().values.any()
binsize = get_bin_size(data,fields=["st","end"])

data["end"] = data["end"] // binsize
data["st"] = data["st"] // binsize

# Sparse to dense
# assume 'data' is pd.DataFrame
# with columns 'st', 'end' (coordinates in bins)
# and column 'oe' (values)
L = data["end"].max() - data["st"].min() + 1
array = np.zeros(shape=(L,L))
X = (data["st"]-data["st"].min()).values
Y = (data["end"]-data["st"].min()).values
array[X,Y] = data["oe"]
array[Y,X] = data["oe"]

# slow; use ones for debuging only
# global_min = data["st"].min()
# def check_func(series):
#     x = int(series["st"] - global_min)
#     y = int(series["end"] - global_min)
#     return array[x,y] == array[y,x] == series["oe"]
#
# data["check"] = data.apply(check_func,axis = "columns")
# assert np.all(data["check"].values)
# print ("passed check!")

# Now we start the Comparator
# So what shell we do:
# 1. Compute Pearson correlation within 5x5 MB matrices (step by 1MB)
# 2. Cluster
# 3. Draw distribution of distances

start = 0
length = 1500000 // binsize
step = length // 5
#end = min(start+length+step*10,len(array) - length)
end = len(array) - length
assert end > start + length

clustering = AgglomerativeClustering(n_clusters=2)

boundaries = pd.DataFrame()
boundaries["st"] = np.arange(data["st"].min()*binsize,
                             data["st"].min()*binsize+(data["end"].max()+1)*binsize,
                            binsize)
boundaries["chr"] = ["2L"]*len(boundaries["st"])
boundaries["end"] = boundaries["st"] + binsize
boundaries["Strength"] = np.zeros(len(boundaries["st"]))


while start <= end:
    submatrix = array[start:start+length,start:start+length]
    #print (submatrix,np.max(submatrix),np.min(submatrix))
    corr = np.corrcoef(submatrix)
    #print (corr,np.max(corr),np.min(corr))
    corr[np.where(~np.isfinite(corr))] = 0

    clusters = clustering.fit_predict(corr)
    #print (clusters)
    changes = np.equal(clusters[:-1],clusters[1:]) # compare i to i+1
    #print (clusters)
    changes_ids = np.array([ind for ind,val in enumerate(changes) if not val])
    #print (start + changes_ids)
    boundaries.iloc[start + changes_ids,-1] += 1
    start += step
    if (start % (step*10) == 0) and start // step >= 10:
        print (datetime.datetime.now(),"Step ",start // step)

boundaries[["chr","st","end","Strength"]].to_csv(
    fname+".win"+str((length*binsize)//1000)+"kb.boundaries.bedGraph",sep="\t",index=False,header=None)
