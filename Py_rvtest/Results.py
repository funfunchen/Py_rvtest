import numpy as np
import h5py

class SaveData:

    def __init__(self, res, out):
        self._res = res
        self._out = out

    def save_csv(self):
        self._res.to_csv(self._out, sep='\t', index=False, header=True)

class HDF5Store:
    """
    Simple class to append > 2D value to a hdf5 file on disc

    Params:
        d_path: filepath of h5 file
        d_set: dataset name within the file
        d_shape: dataset coulumn dim(not counting main/batch axis)
        d_type: numpy dtype

    Usage:
        hdf5_store.shape = (20, 20, 3)
        hdf5_store = HDF5Store('/tmp/hdf5_store.h5','X', d_shape=(20,20,3))
        x = np.random.random(hdf5_store.shape)
        hdf5_store.append(x)
        hdf5_store.append(x)

    From https://gist.github.com/wassname/a0a75f133831eed1113d052c67cf8633
    """

    def __init__(self, d_path, d_set, d_shape, d_type=np.float64):
        self._d_path = d_path
        self._d_set = d_set
        self._d_col = d_shape  # d_shape should be a tuple
        self._row = 0
        self._d_type = d_type

        with h5py.File(self._d_path, mode='w') as h5f:
            self.dset = h5f.create_dataset(
                self._d_set,
                # shape=(0, ) + self._d_shape[1:],
                shape=(0, ) + self._d_col,
                maxshape=(None,) + self._d_col,
                dtype=self._d_type)

    def append(self, arr):
        with h5py.File(self._d_path, mode='a') as h5f:
            dset = h5f[self._d_set]
            dset.resize((self._row + arr.shape[0], ) + self._d_col)
            dset[self._row:, :] = arr
            self._row += arr.shape[0]

    def append_1d(self, arr):
        with h5py.File(self._d_path, mode='a') as h5f:
            dset = h5f[self._d_set]
            dset.resize((self._row + arr.shape[0], ) + self._d_col)
            dset[self._row:, ] = arr
            self._row += arr.shape[0]


def main():
    shape_1 = (20,)
    shape_2 = (20, 20, 3)
    hdf5_store_1 = HDF5Store('/tmp/hdf5_store.h5', 'X', d_shape=())
    for _ in range(10):
        hdf5_store_1.append_1d(np.random.random(20))
    hdf5_store_2 = HDF5Store('/tmp/hdf5_store.h5', 'Y', d_shape=(20, 3))
    for _ in range(10):
        hdf5_store_2.append(np.random.random(shape_2))


if __name__ == '__main__':
    main()