from __future__ import print_function
import datetime
import numpy as np
import pandas as pd
from shared import sparse2dense
from sklearn.cluster import AgglomerativeClustering


def get_compartments(data, binsize, start, length, step):
    array = sparse2dense(data, default_value = 0.)
    #print (len(array))
    #print (array[len(array)-10:len(array),len(array)-10:len(array)])
    #print (array[1430,1430])

    # Now we start the Comparator
    # So what shell we do:
    # 1. Compute Pearson correlation within 5x5 MB matrices (step by 1MB)
    # 2. Cluster
    # 3. Draw distribution of distances

    end = len(array) - length
    assert end > start + length

    clustering = AgglomerativeClustering(n_clusters=2)

    boundaries = pd.DataFrame()
    boundaries["st"] = np.arange(data["st"].min() * binsize,
                                 data["st"].min() * binsize + (data["end"].max() + 1) * binsize,
                                 binsize)
    print (data.loc[:,"chr"].iloc[0])
    boundaries["chr"] = [data.loc[:,"chr"].iloc[0]] * len(boundaries["st"])
    boundaries["end"] = boundaries["st"] + binsize
    boundaries["Strength"] = np.zeros(len(boundaries["st"]))

    while start <= end:
        submatrix = array[start:start + length, start:start + length]
        # print (submatrix,np.max(submatrix),np.min(submatrix))
        corr = np.corrcoef(submatrix)
        if not np.all(np.isfinite(corr)):
            print ("Error on point ",
                   data.loc[:,"chr"].iloc[0],start*binsize)
            #print (submatrix)
            #break
        #print (submatrix)
        # print (corr,np.max(corr),np.min(corr))
        corr[np.where(~np.isfinite(corr))] = 0

        clusters = clustering.fit_predict(corr)
        # print (clusters)
        changes = np.equal(clusters[:-1], clusters[1:])  # compare i to i+1
        # print (clusters)
        changes_ids = np.array([ind for ind, val in enumerate(changes) if not val])
        # print (start + changes_ids)
        boundaries.iloc[start + changes_ids, -1] += 1
        start += step
        if (start % (step * 10) == 0) and start // step >= 10:
            print(datetime.datetime.now(), "Step ", start // step)

    return boundaries
