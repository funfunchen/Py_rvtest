"""Prepare the matrix for model fitting

Todo:
    * read all formats of vcf-related files
    * clean the data to model-fitting ready matrix
"""
import allel
import h5py
import pandas as pd
import numpy.ma as ma
import numpy as np
from sklearn import preprocessing


class H5pyToMatrix:
    """read h5py file and clean it to the model-fitting ready matrix, produce generator if needed"""

    def __init__(self, data, batch_size=50000):
        try:
            self._call = h5py.File(data, mode='r')
        except OSError as e:
            print('Could not open the file because {}'.format(str(e)))
            raise
        self._gt = self._call['calldata/GT']
        self._samples = self._call['samples']
        self._nx, self._ny, self._nz = self._gt.shape # the length of each dim
        self._batch_size = batch_size
        self._num_batches = int(np.ceil(self._nx / batch_size))
    
    def __str__(self):
        return 'vcf size: ({}, {})'.format(self._nx, self._ny)

    def __repr__(self):
        return 'H5py file shape: ({}, {}, {})'.format(self._nx, self._ny, self._nz)

    def data_generator(self, z=None):
        """generate batchs of genetype data

	    :param z: the selected column #, should be a list or None as default
        """
        batch_size = self._batch_size
        genotype_batch_indexes = [[i*batch_size, (i+1)*batch_size] for i in range(self._num_batches)]
        for k, (x,y) in enumerate(genotype_batch_indexes, 1):
            if z is not None:
                batch_geno = allel.GenotypeArray(self._gt[x:y, z])
            else:
                batch_geno =  allel.GenotypeArray(self._gt[x:y])
                if k == genotype_batch_indexes:
                    batch_geno =  allel.GenotypeArray(self._gt[x:y]) # deal with the last batch
            batch_alt = batch_geno.to_n_alt(fill=0).astype('float64')   # missing is '0'
            yield batch_alt

    def get_geno(self, m=0, n=0, z=None):
        """return the subset or whole genotype data in the vcf files

        :param m: the beginning row #
        :param n: the ending row #
        :param z: the selected column #, should be a list or None as default
        :return: the genotype data, which fill the missing cells with average value
        """
        if m + n > 0 and z is not None:  # need to be more flexible
            gc = allel.GenotypeArray(self._gt[m:n, z])
        elif m + n > 0 and z is None:
            gc = allel.GenotypeArray(self._gt[m:n])
        else:
            gc = allel.GenotypeArray(self._gt[...])
        gc_alt = gc.to_n_alt(fill=-1).astype('float64')   # missing is '-1'
        gc_alt_ma = ma.masked_less(gc_alt, 0)
        ma_mean = gc_alt_ma.mean(axis=1)
        np.copyto(gc_alt_ma, ma_mean[..., None], where=gc_alt_ma.mask)
        return gc_alt_ma.data

    @property  # read-only
    def get_snp(self):
        return self._call['variants/ID']

    @property
    def get_sample_size(self):
        return self._ny
    
    @property
    def get_num_batches(self):
        return self._num_batches

    def check_sample_matched(self, y_sam):
        samples = self._samples[:].astype(y_sam.dtype)
        return np.all(samples == y_sam)

    def match_index(self, y_sam):
        samples = list(self._samples[:])
        matched_index = [samples.index(s) for s in y_sam
                         if s in samples]
        return matched_index


class PheToMatrix:
    """read ped or phenotype file and clean it to vector or matrix"""

    def __init__(self, data):
        self._ped = pd.read_csv(data, sep='\t', index_col=None)
        self._pheno = self._ped.iloc[:, 5]
        self._y = self._pheno.values.ravel()
        # self._yc = self._y - self._y.mean()  # centered the data

    def __str__(self):
        return 'Ped file shape sample size: {} '.format(self._y.shape)

    @property
    def get_id(self):
        return self._ped.iloc[:, 0].values

    @property
    def get_pheno(self):
        trait = self._ped.iloc[:, 5].values
        return preprocessing.scale(trait, with_mean=True, with_std=False)

    @property
    def get_cov(self):
        cov = self._ped.iloc[:, 6:].values
        return preprocessing.scale(cov, with_mean=True, with_std=False)


class VcfToMatrix:

    """only for small vcf files"""

    def __init__(self, data):
        self._call = allel.read_vcf(data)
        self._gt = self._call['calldata/GT']
        self._nx, self._ny, self._nz = self._gt.shape   # the length of each dim

def main():
    file = './data/chr3_1000.h5'
    s_vcf = H5pyToMatrix(file)
    s_ge = sub_vcf.get_geno(2, 6, [1, 2, 3])
    print(s_ge.shape)


if __name__ == '__main__':
    main()
