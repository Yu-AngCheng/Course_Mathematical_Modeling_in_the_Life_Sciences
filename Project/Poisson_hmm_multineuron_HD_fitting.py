import numpy as np
import scipy.io as io
import os
from datetime import datetime
import matplotlib.pyplot as plt
from PoissonHMM import PoissonHMM


def fitHMM(Y, nstate, nfeature, maxrun, fix_means):
    modellist = list()
    scorelist = list()
    for run in range(maxrun):
        tempmodel = PoissonHMM(n_components=nstate,
                               means_prior=0,
                               n_iter=10000, tol=1e-5,
                               params="st", init_params="st")
        tempmodel.means_ = fix_means
        tempmodel.fit(Y)
        modellist.append(tempmodel)
        scorelist.append(tempmodel.score(Y))
    index = scorelist.index(max(scorelist))
    finalmodel = modellist[index]
    return finalmodel


# load data
if os.path.isfile("Poisson_multineuron_HD_HMMFIT.npy"):
    raise Exception("!!'Poisson_multineuron_HD_HMMFIT.npy' already exists!!")
alldat = np.load('Steinmetz_main_multineuron.npy', allow_pickle=True)
PoissonHMMFIT = np.load('Poisson_multineuron_HMMFIT.npy', allow_pickle=True)
allSpk = np.concatenate([alldat[k]['groupSpk']
                         for k in range(len(alldat))], axis=1)
allSession = np.concatenate([alldat[k]['session']
                             for k in range(len(alldat))], axis=0)
allGroupIndex = np.concatenate([alldat[k]['groupindex']
                                for k in range(len(alldat))], axis=0)
del alldat

# Cut the first 50 Time Bins
BinCut = 50
allSpk = allSpk[:, :, BinCut:]

# fitting parameters
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
statenumberlist = [4]
maxrun = 10
nFeature = 2
nTrial = allSpk.shape[1]
nAreasConnection = 15
nBin = allSpk.shape[2]
Poisson_multineuron_HD_HMMFIT = np.array([])
count = 0
arealists = [[0, 1], [0, 2], [0, 3], [0, 4], [0, 5],
             [1, 2], [1, 3], [1, 4], [1, 5],
             [2, 3], [2, 4], [2, 5],
             [3, 4], [3, 5],
             [4, 5]]

for nState in statenumberlist:

    errorTrial = []
    statechain = np.full([nTrial, nAreasConnection, nBin], np.nan)
    statechain_probability = np.full([nTrial, nAreasConnection, nBin, nState], np.nan)
    converged = np.full([nTrial, nAreasConnection], np.nan)
    minusloglikelihood = np.full([nTrial, nAreasConnection], np.nan)
    means = np.full([nTrial, nAreasConnection, nState, nFeature], np.nan)
    transmat = np.full([nTrial, nAreasConnection, nState, nState], np.nan)

    for trial in range(nTrial):
        for arealist in arealists:

            tic = datetime.now()
            areaidx = arealists.index(arealist)
            print('Percent:  {0:1.3g}%'.format(count / (len(statenumberlist) * len(arealists) * nTrial) * 100))

            Y = allSpk[arealist, trial, :].transpose()
            count = count + 1
            if np.isnan(Y).any():
                continue
            else:
                try:
                    tempmeans = PoissonHMMFIT[1]['means'][trial, :, :]
                    B1_S1 = tempmeans[arealist[0], 0].item()
                    B1_S2 = tempmeans[arealist[0], 1].item()
                    if B1_S1 > B1_S2:
                        temp = B1_S1
                        B1_S1 = B1_S2
                        B1_S2 = temp
                    B2_S1 = tempmeans[arealist[1], 0].item()
                    B2_S2 = tempmeans[arealist[1], 1].item()
                    if B2_S1 > B2_S2:
                        temp = B2_S1
                        B2_S1 = B2_S2
                        B2_S2 = temp
                    fixmeans = np.array([[B1_S1, B2_S1], [B1_S1, B2_S2], [B1_S2, B2_S1], [B1_S2, B2_S2]])
                    finalmodel = fitHMM(Y, nState, nFeature, maxrun, fixmeans)
                    statechain[trial, areaidx, :] = finalmodel.predict(Y)
                    statechain_probability[trial, areaidx, :, :] = finalmodel.predict_proba(Y)
                    converged[trial, areaidx] = finalmodel.monitor_.converged
                    minusloglikelihood[trial, areaidx] = - finalmodel.score(Y)
                    means[trial, areaidx, :, :] = finalmodel.means_
                    transmat[trial, areaidx, :, :] = finalmodel.transmat_
                except:
                    errorTrial.append(trial)

            toc = datetime.now()
            print('Elapsed time: %f seconds' % (toc - tic).total_seconds())

    idx = statenumberlist.index(nState)
    temp = dict()
    temp['statechain'] = statechain
    temp['statechain_probability'] = statechain_probability
    temp['minusloglikelihood'] = minusloglikelihood
    temp['means'] = means
    temp['transmat'] = transmat
    temp['converged'] = converged
    temp['allSession'] = allSession
    temp['errorTrial'] = errorTrial
    Poisson_multineuron_HD_HMMFIT = np.concatenate((Poisson_multineuron_HD_HMMFIT, np.array([temp])))

np.save('Poisson_multineuron_HD_HMMFIT.npy', Poisson_multineuron_HD_HMMFIT)
