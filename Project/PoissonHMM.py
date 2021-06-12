import numpy as np
from sklearn import cluster
from sklearn.utils import check_random_state
from scipy.stats import poisson
from hmmlearn.base import _BaseHMM


# __all__ = ["PoissonHMM"]

def _check_and_set_n_features(model, X):
    _, n_features = X.shape
    if hasattr(model, "n_features") and model.n_features != n_features:
        raise ValueError("Unexpected number of dimensions, got {} but "
                         "expected {}".format(n_features, model.n_features))
    model.n_features = n_features


def log_possion_mass(X, lams):
    n_samples, n_dim = X.shape
    logp = np.zeros((n_samples, 1))
    for dim in np.arange(n_dim):
        logp = logp + poisson.logpmf(X[:, dim], lams[dim])[:, None]
    return logp


class PoissonHMM(_BaseHMM):
    r"""Hidden Markov Model with Possion emissions.
    Parameters
    ----------
    n_components : int
        Number of states.
    startprob_prior : array, shape (n_components, ), optional
        Parameters of the Dirichlet prior distribution for
        :attr:`startprob_`.
    transmat_prior : array, shape (n_components, n_components), optional
        Parameters of the Dirichlet prior distribution for each row
        of the transition probabilities :attr:`transmat_`.
    means_prior, means_weight : array, shape (n_components, ), optional
        Mean and precision of the Normal prior distribtion for
        :attr:`means_`.
    algorithm : string, optional
        Decoder algorithm. Must be one of "viterbi" or`"map".
        Defaults to "viterbi".
    random_state: RandomState or an int seed, optional
        A random number generator instance.
    n_iter : int, optional
        Maximum number of iterations to perform.
    tol : float, optional
        Convergence threshold. EM will stop if the gain in log-likelihood
        is below this value.
    verbose : bool, optional
        When ``True`` per-iteration convergence reports are printed
        to :data:`sys.stderr`. You can diagnose convergence via the
        :attr:`monitor_` attribute.
    params : string, optional
        Controls which parameters are updated in the training
        process.  Can contain any combination of 's' for startprob,
        't' for transmat and 'm' for means. Defaults to all parameters.
    init_params : string, optional
        Controls which parameters are initialized prior to
        training.  Can contain any combination of 's' for
        startprob, 't' for transmat and 'm' for means.
        Defaults to all parameters.
    Attributes
    ----------
    n_features : int
        Dimensionality of the Possion emissions
        Features are assumed to be independent from each other.
    monitor\_ : ConvergenceMonitor
        Monitor object used to check the convergence of EM.
    startprob\_ : array, shape (n_components, )
        Initial state occupation distribution.
    transmat\_ : array, shape (n_components, n_components)
        Matrix of transition probabilities between states.
    means\_ : array, shape (n_components, n_features)
        Mean parameters for each state.
    """

    def __init__(self, n_components=1,
                 startprob_prior=1.0, transmat_prior=1.0,
                 means_prior=0, means_weight=0,
                 algorithm="viterbi", random_state=None,
                 n_iter=10, tol=1e-2, verbose=False,
                 params="stm", init_params="stm"):
        _BaseHMM.__init__(self, n_components,
                          startprob_prior=startprob_prior,
                          transmat_prior=transmat_prior, algorithm=algorithm,
                          random_state=random_state, n_iter=n_iter,
                          tol=tol, params=params, verbose=verbose,
                          init_params=init_params)

        self.means_prior = means_prior
        self.means_weight = means_weight

    def _get_n_fit_scalars_per_param(self):
        nc = self.n_components
        nf = self.n_features
        return {
            "s": nc - 1,
            "t": nc * (nc - 1),
            "m": nc * nf,
        }

    def _init(self, X, lengths=None):
        _check_and_set_n_features(self, X)
        super()._init(X, lengths=lengths)

        if self._needs_init("m", "means_"):
            kmeans = cluster.KMeans(n_clusters=self.n_components,
                                    random_state=self.random_state)
            kmeans.fit(X)
            self.means_ = kmeans.cluster_centers_

    def _check(self):
        super()._check()

        self.means_ = np.asarray(self.means_)
        self.n_features = self.means_.shape[1]

    def _compute_log_likelihood(self, X):
        logp = np.zeros_like(X)
        for lams in self.means_:
            logp = np.concatenate((logp, log_possion_mass(X, lams)), axis=1)
        return logp[:, 1:]

    def _generate_sample_from_state(self, state, random_state=None):
        random_state = check_random_state(random_state)
        sample = np.empty(self.n_features)
        for i in np.arange(self.n_features):
            sample[i] = random_state.poisson(self.means_[state, i])
        return sample

    def _initialize_sufficient_statistics(self):
        stats = super()._initialize_sufficient_statistics()
        stats['post'] = np.zeros(self.n_components)
        stats['obs'] = np.zeros((self.n_components, self.n_features))
        return stats

    def _accumulate_sufficient_statistics(self, stats, obs, framelogprob,
                                          posteriors, fwdlattice, bwdlattice):
        super()._accumulate_sufficient_statistics(
            stats, obs, framelogprob, posteriors, fwdlattice, bwdlattice)

        if 'm' in self.params:
            stats['post'] += posteriors.sum(axis=0)
            stats['obs'] += np.dot(posteriors.T, obs)

    def _do_mstep(self, stats):
        super()._do_mstep(stats)

        means_prior = self.means_prior
        means_weight = self.means_weight

        denom = stats['post'][:, None]
        if 'm' in self.params:
            self.means_ = ((means_weight * means_prior + stats['obs'])
                           / (means_weight + denom))
