import argparse
import os
import subprocess
from sys import argv
from scipy.spatial import ConvexHull
import numpy as np
import vtk
from glob import glob
import re
import pandas as pd
from vtk.util import numpy_support

def start_index(args):
    res = 0
    for file in glob(f"{args.out}*_stress.*"):
        i = int(re.match(fr"{args.out}(\d+)_stress\..*", file.replace("\\", "/")).groups()[0])
        if i >= res:
            res = i + 1

    return res


def generate_3d_convex_hull_obj(path, points_num=30):
    points = np.random.uniform(-1, 1, size=(points_num, 3))

    # 凸包を計算
    hull = ConvexHull(points)

    # 凸包の面のインデックスを取得
    faces = hull.simplices.tolist()

    # objファイルとして出力
    with open(path, "w") as f:
        for v in points:
            f.write("v {} {} {}\n".format(v[0], v[1], v[2]))
        for face in faces:
            f.write("f {} {} {}\n".format(
                face[0] + 1, face[1] + 1, face[2] + 1))


def obj_to_voxel_and_sdf(path):
    # OBJファイルからメッシュを読み込む
    reader = vtk.vtkOBJReader()
    reader.SetFileName(path)

    # vtkDelaunay3Dを作成する
    delaunay = vtk.vtkDelaunay3D()
    delaunay.SetInputConnection(reader.GetOutputPort())

    # メッシュをボクセルに変換する
    voxelGrid = vtk.vtkVoxelModeller()
    voxelGrid.SetInputConnection(delaunay.GetOutputPort())
    voxelGrid.SetSampleDimensions(64, 64, 64)
    voxelGrid.SetModelBounds(-1, 1, -1, 1, -1, 1)
    voxelGrid.SetScalarTypeToFloat()

    distanceFilter = vtk.vtkImageEuclideanDistance()
    distanceFilter.SetInputConnection(voxelGrid.GetOutputPort())
    distanceFilter.Update()

    # Convert voxel to NumPy array
    voxelArray = numpy_support.vtk_to_numpy(voxelGrid.GetOutput().GetPointData().GetScalars())
    voxelArray = voxelArray.reshape((64, 64, 64))

    # Convert SDF to NumPy array
    sdfArray = numpy_support.vtk_to_numpy(distanceFilter.GetOutput().GetPointData().GetScalars())
    sdfArray = sdfArray.reshape((64, 64, 64))

    return voxelArray, sdfArray


def save_feather(feather_path, grid: np.ndarray):
    grid_2d = grid.reshape(
        (grid.shape[0] * grid.shape[2], grid.shape[1]))
    df = pd.DataFrame(data=grid_2d, columns=[
                      str(i) for i in range(grid_2d.shape[1])])
    df.to_feather(feather_path)


os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.chdir("../build/RELEASE/")

parser = argparse.ArgumentParser()
parser.add_argument("param", type=str)
parser.add_argument("--shape", type=int, default=64)
parser.add_argument("-i", "--iter", type=int, default=64)
parser.add_argument("-o", "--out", default="")
parser.add_argument("--start", default="")

parsed_args, _ = parser.parse_known_args()
args = argv

error = 0
i = 0
while start_index(parsed_args) < parsed_args.iter * parsed_args.shape:
    start = start_index(parsed_args) % parsed_args.iter
    i = start_index(parsed_args) // parsed_args.iter
    if start == 0:
        generate_3d_convex_hull_obj(
            f"{parsed_args.out}{i * parsed_args.iter}_input.obj")
        voxel, sdf = obj_to_voxel_and_sdf(
            f"{parsed_args.out}{i * parsed_args.iter}_input.obj")
        for l in range(parsed_args.iter):
            save_feather(f"{parsed_args.out}{i * parsed_args.iter + l}_voxel.fea", voxel)
            save_feather(f"{parsed_args.out}{i * parsed_args.iter + l}_sdf.fea", sdf)
    offset = i * parsed_args.iter
    while start_index(parsed_args) < parsed_args.iter * (i + 1):
        start = start_index(parsed_args) % parsed_args.iter
        print(f"start: {start}, offset:{offset}")
        subprocess.run([
            "GenerateFractureData.exe",
            f"{parsed_args.out}{i * parsed_args.iter}_input.obj",
            parsed_args.param,
            "-o",
            parsed_args.out,
            "-i",
            f"{parsed_args.iter}",
            "--start",
            f"{start}",
            "--offset",
            f"{offset}",
        ])
    
