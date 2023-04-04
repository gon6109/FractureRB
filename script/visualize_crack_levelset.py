import matplotlib.pyplot as plt
import numpy as np
import argparse
import pandas
import cv2
from mayavi.mlab import volume_slice
from mayavi import mlab


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


parser = argparse.ArgumentParser()
parser.add_argument("input")
parser.add_argument("-s", "--size", type=int, default=256)

args = parser.parse_args()

csv = pandas.read_csv(args.input, header=None, dtype=float)
data = csv.to_numpy()

print(type(data))

data = np.split(data, args.size)
data = np.stack(data, axis=2)

fig = volume_slice(data)


@mlab.animate(delay=30, ui=False)
def anim():
    scene = mlab.gcf()
    scene.scene.camera.position = [
        872.539173801294, 786.7491893734637, 731.321600028871]
    scene.scene.camera.focal_point = [
        61.65661489633652, 111.77788010112808, 71.09760306504057]
    scene.scene.camera.view_angle = 30.0
    scene.scene.camera.view_up = [0.0, 0.0, 1.0]
    scene.scene.camera.clipping_range = [794.4074033069271, 1438.6129430793414]
    scene.scene.camera.compute_view_plane_normal()

    for i in range(200):
        yield

    for i in range(args.size):
        fig.ipw.origin = [i, 0.5,  0.5]
        fig.ipw.point1 = [i, 256.5,   0.5]
        fig.ipw.point2 = [i, 0.5, 256.5]
        fig.ipw.slice_index = i
        fig.ipw.slice_position = i
        yield


anim()
mlab.show()
