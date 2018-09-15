from Py_rvtest.DataConsolidator import *
import argparse
from tqdm import tqdm
import pandas as pd
import Py_rvtest.Regression as Reg
import Py_rvtest.Results as Res
import Py_rvtest.MinimalSufficient as MinS

# parse the arguments
parser = argparse.ArgumentParser()
parser.add_argument("--geno", "-g", help="genotype file, HDF5 format")
parser.add_argument("--out", "-o", type=str, help="output file name")
parser.add_argument("--chunk", type=int, help="the size of a chunk")
parser.add_argument("--ped", "-p", help="phenotype file, ped format")
args = parser.parse_args()

ge = H5pyToMatrix(args.geno)
pheno = PheToMatrix(args.ped)
y_samples = pheno.get_id

if ge.check_sample_matched(y_samples):
    mat_index = None
else:
    mat_index = ge.match_index(y_samples)  # match the index, y_sample should be a list


snps = ge.get_snp
snp_num = len(snps)
chunk = args.chunk
traits = pheno.get_pheno
cos = pheno.get_cov
_, cos_col = cos.shape
# cols = ['SNP', 'Beta', 'SE', 'Z-score', 'P']
# out_name = args.out + '.csv'
# the results df
# df = pd.DataFrame(columns=cols)

with h5py.File("temp_min.h5", "w") as f:
    pheno_min = MinS.PhenoMinSuffStat(traits, cos)
    pheno_min_ct = pheno_min.get_ZcYt
    pheno_min_tt = pheno_min.get_YtYt
    pheno_min_cc = pheno_min.get_ZcZc
    f.create_dataset("min_ct", data=pheno_min_ct, dtype=np.float64)
    f.create_dataset("min_tt", data=pheno_min_tt, dtype=np.float64)
    f.create_dataset("min_cc", data=pheno_min_cc, dtype=np.float64)
    pheno_min_tt = None

min_gy_store = Res.HDF5Store("temp_min.h5", "min_gy", d_shape=(1, 1))
min_gg_store = Res.HDF5Store("temp_min.h5", "min_gg", d_shape=(1, 1))
min_gz_store = Res.HDF5Store("temp_min.h5", "min_gz", d_shape=(1, cos_col))
for i in tqdm(range(0, snp_num, chunk)):
    # Todo
    # iterate the object to get matrix and compute
    # snp_id = snps[i, i+chunk]
    mat = ge.get_geno(i, i+chunk, mat_index)
    geno_min = MinS.GenoMinSuffStat(mat, traits, cos)
    min_gy_store.append(geno_min.get_GmY)
    min_gg_store.append(geno_min.get_GmGm)
    min_gz_store.append(geno_min.get_GmZc)

    # beta, se, z_score, p = Reg.LinearRegression(mat, ped_pheno).cal_p()
    # res = pd.DataFrame({'SNP': snp_id, 'Beta': beta, 'SE': se, 'Z-score': z_score, 'P': p}).round({'Beta': 6,
    #                                                                                                'SE': 5,
    #                                                                                                'Z-score': 6,
    #                                                                                                'P': 6})
    # df = df.append(res, ignore_index=True)




