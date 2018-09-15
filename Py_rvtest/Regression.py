"""Regression models"""
import scipy.stats
import numpy as np


class LinearRegression:

    def __init__(self, geno, pheno):
        self._geno = geno
        self._y = pheno
        self._size = len(self._y)

    @property
    def cal_p(self):
        if self._y.shape[1] == 1:
            beta = np.dot(self._geno, self._y) / self._size
            res = self._y.T - (self._geno * beta[:, None])
            sigma_sq = np.einsum('ij,ji->i', res, res.T)
            se = np.sqrt((sigma_sq / self._size) / self._size)
            z_scores = beta.T / se
            p_values = scipy.stats.norm.sf(abs(z_scores)) * 2
        return beta, se, z_scores, p_values
