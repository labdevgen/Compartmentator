import numpy as np
import matplotlib.pyplot as plt

input_file = "D:/Lab Archive/anopheles/Aatr_VARYA_comb_allValidPairs.2L.25000.pears.txt"
data = np.loadtxt(input_file,skiprows=2,dtype=np.float)
print (data[0:3,0:3])

#data2 = np.zeros(data)
#for i in range(2,len(data))
#    for j in range(i,len(data))
#        if

#err_count = 0
#for i in range(len(data)):
#    diag_ids1 = (np.arange(len(data)-i),np.arange(i,len(data)))
#    diag_ids2 = (np.arange(i, len(data)),np.arange(len(data) - i))
#    av = np.nanmedian(np.abs(data[diag_ids1].flat))
#    if (av == 0 or not np.isfinite(av)):
#        print (av)
#        av = 1
#        err_count += 1

#    data[diag_ids1] = data[diag_ids1] / av
#    data[diag_ids2] = data[diag_ids2] / av
#print (err_count)
#data[np.triu_indices(len(data))] = data[np.tril_indices(len(data))]

plt.imshow(data)
plt.colorbar()
plt.title(input_file)
plt.show()