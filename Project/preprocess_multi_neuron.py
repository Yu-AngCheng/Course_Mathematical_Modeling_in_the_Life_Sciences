import os
import numpy as np

# Data retrieval
fname = []
for j in range(3):
    fname.append('steinmetz_part%d.npz' % j)
    if not os.path.isfile(fname[j]):
        raise Exception("!!! No data !!!")

# Data loading
alldat = np.array([])
for j in range(len(fname)):
    alldat = np.hstack((alldat, np.load('steinmetz_part%d.npz' % j, allow_pickle=True)['dat']))

# calculate the Spikes by every brain group

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

    nareas = nGroup  # only the top nareas regions are taken into account
    NN = len(dat['brain_area'])  # number of neurons
    groupindex = nGroup * np.ones(NN, )  # last one is "other"
    for j in range(nareas):
        groupindex[np.isin(dat['brain_area'], brain_groups[j])] = j
    dat['groupindex'] = groupindex

    groupSpk = np.full((nareas, nTrial, nBin), np.nan)
    for i in range(nareas):
        if np.sum(groupindex == i) != 0:
            groupSpk[i, :, :] = np.sum(dat['spks'][groupindex == i, :, :], axis=0)
        else:
            groupSpk[i, :, :] = np.full((1, nTrial, nBin), np.nan)

    dat['groupSpk'] = groupSpk

    dat['session'] = np.tile(k, nTrial)

    alldat[k] = dat

np.save('Steinmetz_main_multineuron.npy', alldat)
