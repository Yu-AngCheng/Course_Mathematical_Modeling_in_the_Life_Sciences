import os
import requests
import numpy as np

# Data retrieval
fname = []
for j in range(3):
  fname.append('steinmetz_part%d.npz'%j)
url = ["https://osf.io/agvxh/download"]
url.append("https://osf.io/uv3mw/download")
url.append("https://osf.io/ehmw2/download")

for j in range(len(url)):
  if not os.path.isfile(fname[j]):
    try:
      r = requests.get(url[j])
    except requests.ConnectionError:
      print("!!! Failed to download data !!!")
    else:
      if r.status_code != requests.codes.ok:
        print("!!! Failed to download data !!!")
      else:
        with open(fname[j], "wb") as fid:
          fid.write(r.content)

# Data loading
alldat = np.array([])
for j in range(len(fname)):
    alldat = np.hstack((alldat, np.load('steinmetz_part%d.npz' % j, allow_pickle=True)['dat']))

# sort the spikes by different brain groups
brain_groups = ["VISa", "VISam", "VISl", "VISp", "VISpm", "VISrl",  # visual cortex
                "CL", "LD", "LGd", "LH", "LP", "MD", "MG", "PO", "POL", "PT", "RT", "SPF", "TH", "VAL", "VPL", "VPM",
                # thalamus
                "CA", "CA1", "CA2", "CA3", "DG", "SUB", "POST",  # hippocampal
                "ACA", "AUD", "COA", "DP", "ILA", "MOp", "MOs", "OLF", "ORB", "ORBm", "PIR", "PL", "SSp", "SSs", "RSP",
                " TT",  # non-visual cortex
                "APN", "IC", "MB", "MRN", "NB", "PAG", "RN", "SCs", "SCm", "SCig", "SCsg", "ZI",  # midbrain
                "ACB", "CP", "GPe", "LS", "LSc", "LSr", "MS", "OT", "SNr", "SI",  # basal ganglia
                "BLA", "BMA", "EP", "EPd", "MEA"  # cortical subplate
                ]
nGroup = len(brain_groups)

for k in range(len(alldat)):

    dat = alldat[k]
    tmp = dat["spks"].shape  # Neuron * Trial * Time Bin
    nNeuron = tmp[0]
    nTrial = tmp[1]
    nBin = tmp[2]

    nareas = 6  # only the top nareas regions are taken into account
    NN = len(dat['brain_area'])  # number of neurons
    groupindex = nGroup * np.ones(NN, )  # last one is "other"
    for j in range(nareas):
        groupindex[np.isin(dat['brain_area'], brain_groups[j])] = j
    dat['groupindex'] = groupindex

    groupSpk = np.full((0, nBin), np.nan)
    groupindex_sorted = np.array([])
    sessionnumber = np.array([])
    trialnumber = np.array([])
    neuronnumber = np.array([])
    for i in range(nareas):
        number_neurons = np.sum(groupindex == i)
        if number_neurons != 0:
            temp = dat['spks'][groupindex == i, :, :]
            groupSpk = np.concatenate((groupSpk, temp.reshape(-1,nBin)), axis=0)
            groupindex_sorted = np.concatenate((groupindex_sorted, np.tile(np.array([i]),temp.shape[1]*temp.shape[0])))
            sessionnumber = np.concatenate((sessionnumber, np.tile(np.array([k]),temp.shape[1]*temp.shape[0])))
            trialnumber = np.concatenate((trialnumber, np.tile(np.arange(temp.shape[1]),temp.shape[0])))
            neuronnumber = np.concatenate((neuronnumber, np.repeat(np.arange(temp.shape[0]),temp.shape[1])))
            a = 1

    dat['groupSpk'] = groupSpk
    dat['groupindex_sorted'] = groupindex_sorted
    dat['trialnumber'] = trialnumber
    dat['session'] = sessionnumber
    dat['neuronnumber'] = neuronnumber

    alldat[k] = dat

np.save('Steinmetz_main_singleneuron.npy', alldat)
