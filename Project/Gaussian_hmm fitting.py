import numpy as np
import hmmlearn.hmm as hmm
import scipy.io as io
import time as tm
import os
from datetime import datetime
import matplotlib.pyplot as plt


def fitHMM(Y, nstate, nfeature, maxrun):
    modellist = list()
    scorelist = list()
    for run in range(maxrun):
        tempmodel = hmm.GaussianHMM(n_components=nstate,
                                    covariance_type='diag',
                                    n_iter=10000, tol=1e-5)
        tempmodel.fit(Y)
        modellist.append(tempmodel)
        scorelist.append(tempmodel.score(Y))
    index = scorelist.index(max(scorelist))
    finalmodel = modellist[index]
    return finalmodel


# load data
if os.path.isfile("HMMFIT.npy"):
    raise Exception("!!'HMMFIT.npy' already exists!!")
alldat = np.load('Steinmetz_main_braingroup.npy', allow_pickle=True)
allFiringrate = np.concatenate([alldat[k]['groupFiringRate']
                                for k in range(len(alldat))], axis=1)
allSession = np.concatenate([alldat[k]['session']
                             for k in range(len(alldat))], axis=0)
del alldat

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
statenumberlist = [1, 2, 3, 4, 5]
maxrun = 10
nFeature = 1
nTrial = allFiringrate.shape[1]
nAreas = 6
nBin = allFiringrate.shape[2]
HMMFIT = np.array([])
count = 0

for nState in statenumberlist:

    statechain = np.full([nTrial, nAreas, nBin], np.nan)
    statechain_probability = np.full([nTrial, nAreas, nBin, nState], np.nan)
    converged = np.full([nTrial, nAreas], np.nan)
    minusloglikelihood = np.full([nTrial, nAreas], np.nan)
    means = np.full([nTrial, nAreas, nState, nFeature], np.nan)
    covars = np.full([nTrial, nAreas, nState, nFeature, nFeature], np.nan)
    transmat = np.full([nTrial, nAreas, nState, nState], np.nan)

    for trial in range(nTrial):
        tic = datetime.now()
        for area in range(nAreas):

            print('Percent:  {0:1.3g}%'.format(count / (len(statenumberlist) * nTrial * nAreas) * 100))

            Y = allFiringrate[area, trial, :].reshape(-1, 1)
            count = count + 1
            if np.isnan(Y).all():
                continue
            else:
                finalmodel = fitHMM(Y, nState, nFeature, maxrun)

            statechain[trial, area, :] = finalmodel.predict(Y)
            statechain_probability[trial, area, :, :] = finalmodel.predict_proba(Y)
            converged[trial, area] = finalmodel.monitor_.converged
            minusloglikelihood[trial, area] = - finalmodel.score(Y)
            means[trial, area, :, :] = finalmodel.means_
            covars[trial, area, :, :, :] = finalmodel.covars_
            transmat[trial, area, :, :] = finalmodel.transmat_

        toc = datetime.now()
        print('Elapsed time: %f seconds' % (toc - tic).total_seconds())

    idx = statenumberlist.index(nState)
    temp = dict()
    temp['statechain'] = statechain
    temp['statechain_probability'] = statechain_probability
    temp['minusloglikelihood'] = minusloglikelihood
    temp['means'] = means
    temp['covars'] = covars
    temp['transmat'] = transmat
    temp['converged'] = converged
    HMMFIT = np.concatenate((HMMFIT, np.array([temp])))

np.save('HMMFIT.npy', HMMFIT)
