import numpy as np
import argparse
import pandas
import cv2
from mayavi.mlab import volume_slice
from mayavi import mlab
from scipy.sparse.linalg import eigs


def resize(array, size):
    temp = np.zeros((size[0], size[1], array.shape[2]))
    res = np.zeros(size)

    for idx in range(array.shape[2]):
        img = array[:, :, idx]
        img_sm = cv2.resize(
            img, (size[0], size[1]), interpolation=cv2.INTER_CUBIC)
        temp[:, :, idx] = img_sm

    for idx in range(size[0]):
        img = temp[idx, :, :]
        img_sm = cv2.resize(
            img, (size[1], size[2]), interpolation=cv2.INTER_CUBIC)
        res[idx, :, :] = img_sm

    return res


def stress_mat(csv, x, y, z):
    def val(csv, x, y, z, w):
        return csv.iat[x + args.size * y + args.size * args.size * z, w]

    return np.array([
        [val(csv, x, y, z, 0), val(csv, x, y, z, 3), val(csv, x, y, z, 5)],
        [val(csv, x, y, z, 3), val(csv, x, y, z, 1), val(csv, x, y, z, 4)],
        [val(csv, x, y, z, 5), val(csv, x, y, z, 4), val(csv, x, y, z, 2)]
    ])


parser = argparse.ArgumentParser()
parser.add_argument("input")
parser.add_argument("-s", "--size", type=int, default=64)

args = parser.parse_args()

csv = pandas.read_csv(args.input, header=None, dtype=float)
print(type(csv))

u = np.zeros((args.size, args.size, args.size))
v = np.zeros((args.size, args.size, args.size))
w = np.zeros((args.size, args.size, args.size))
data = np.zeros((args.size, args.size, args.size))
for x in range(args.size):
    for y in range(args.size):
        for z in range(args.size):
            val, vec = eigs(stress_mat(csv, x, y, z), 1)
            data[x, y, z] = val[0]
            u[x, y, z] = vec[0, 0]
            v[x, y, z] = vec[1, 0]
            w[x, y, z] = vec[2, 0]

# src = mlab.pipeline.vector_field(u, v, w)
# mlab.pipeline.vector_cut_plane(src, mask_points=2, scale_factor=3)

mlab.pipeline.volume(mlab.pipeline.scalar_field(data), vmin=0, vmax=0.8)
mlab.show()
