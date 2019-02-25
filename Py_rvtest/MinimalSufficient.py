"""Calculate and prepare the minimal sufficient statistics."""
from scipy import linalg
import numpy as np
import os
os.environ["CUDA_VISIBLE_DEVICES"]="-1"
import tensorflow as tf


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

    def tensor_cal(self, batch_size=50000, z=None):
        sample_size = self._geno.get_sample_size
        #generator
        geno_gr = self._geno.data_generator(batch_size)
        #placeholders
        t_geno = tf.placeholder(dtype=tf.float64, name="geno")
        t_trait = tf.constant(self._trait, dtype=tf.float64, name="trait")
        t_trait = tf.reshape(t_trait, [-1,1])
        t_cov = tf.constant(self._cov, dtype=tf.float64, name="cov")
        geno_data = tf.data.Dataset.from_generator(geno_gr, tf.float64)
        iter_geno = geno_data.make_one_shot_iterator()
        next_geno = iter_geno.get_next()
        num_iter = self._geno.get_num_batches
        
        with tf.Session() as sess:



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

    def tensor_cal(self):
        t_trait = tf.constant(self._trait, dtype=tf.float64, name="trait")
        t_trait = tf.reshape(t_trait, [-1,1])
        t_cov = tf.constant(self._cov, dtype=tf.float64, name="cov")
        t_ZcYt = tf.matmul(tf.transpose(t_cov), t_trait, [-1,1])
        t_ZcZc = tf.matmul(tf.transpose(t_cov), t_cov)
        t_YtYt = tf.matmul(tf.transpose(t_trait), t_trait)

        with tf.Session() as sess:
            _t1, _t2, _t3 = sess.run(t_ZcYt, t_ZcZc, t_YtYt)

        return _t1, _t2, _t3
