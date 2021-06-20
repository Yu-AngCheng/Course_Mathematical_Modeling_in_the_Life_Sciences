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
nAreas = 6
nBin = allSpk.shape[2]
Poisson_multineuron_HD_HMMFIT = np.array([])
count = 0
arealist = [1, 2]

for nState in statenumberlist:

    errorTrial = []
    statechain = np.full([nTrial, nAreas, nBin], np.nan)
    statechain_probability = np.full([nTrial, nAreas, nBin, nState], np.nan)
    converged = np.full([nTrial, nAreas], np.nan)
    minusloglikelihood = np.full([nTrial, nAreas], np.nan)
    means = np.full([nTrial, nAreas, nState, nFeature], np.nan)
    transmat = np.full([nTrial, nAreas, nState, nState], np.nan)

    for trial in range(nTrial):

        tic = datetime.now()

        print('Percent:  {0:1.3g}%'.format(count / (len(statenumberlist) * nTrial) * 100))

        Y = allSpk[arealist, trial, :].transpose()
        count = count + 1
        if np.isnan(Y).any():
            continue
        else:
            try:
                tempmeans = PoissonHMMFIT[1]['means'][trial, :, :]
                B1_S1 = tempmeans[arealist[0], 0].item()
                B1_S2 = tempmeans[arealist[0], 1].item()
                B2_S1 = tempmeans[arealist[1], 0].item()
                B2_S2 = tempmeans[arealist[1], 1].item()
                fixmeans = np.array([[B1_S1, B2_S1], [B1_S1, B2_S2], [B1_S2, B2_S1], [B1_S2, B2_S2]])
                finalmodel = fitHMM(Y, nState, nFeature, maxrun, fixmeans)
                statechain[trial, arealist, :] = finalmodel.predict(Y)
                statechain_probability[trial, arealist, :, :] = finalmodel.predict_proba(Y)
                converged[trial, arealist] = finalmodel.monitor_.converged
                minusloglikelihood[trial, arealist] = - finalmodel.score(Y)
                means[trial, arealist, :, :] = finalmodel.means_
                transmat[trial, arealist, :, :] = finalmodel.transmat_
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
