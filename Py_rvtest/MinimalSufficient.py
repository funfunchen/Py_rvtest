"""Calculate and prepare the minimal sufficient statistics."""
from scipy import linalg
import numpy as np


class GenoMinSuffStat:

    def __init__(self, geno, trait, cov):
        self._geno = geno
        self._trait = trait
        self._cov = cov
        self._size = len(self._trait)

    @property
    def get_GmY(self):  # Gm: n*p, Y: p*1
        return self._geno.dot(self._trait)

    @property
    def get_GmGm(self): # Gm: n*p, GmGm: just need diagonal n
        return np.sum(np.square(self._geno), axis=1)

    @property
    def get_GmZc(self):
        return self._geno.dot(self._cov)


class PhenoMinSuffStat(object):
    def __init__(self, t, c):
        self._trait = t
        self._cov = c
        self._size = len(self._trait)

    @property
    def get_ZcYt(self):
        return self._cov.T.dot(self._trait)

    @property
    def get_ZcZc(self):
        return self._cov.T.dot(self._cov)

    @property
    def get_YtYt(self):
        return self._trait.T.dot(self._trait)

