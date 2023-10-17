import argparse
import json
import os
import subprocess
import sys
import socket
from datetime import datetime
from sys import argv
import requests
from scipy.spatial import ConvexHull
import numpy as np
import vtk
from glob import glob
import re
import pandas as pd
from vtk.util import numpy_support

sys.path.append(
    os.path.dirname(os.path.abspath(__file__)) + "/../thirdparty/neural-fracture/src"
)
from dataset import shape
from dataset import crack
from dataset import fragment


def start_index(args):
    res = 0
    out = args.out.replace("\\", "/")
    for file in glob(f"{out}*_contact.*"):
        i = int(
            re.match(rf"{out}(\d+)_contact\..*", file.replace("\\", "/")).groups()[
                0
            ]
        )
        if i >= res:
            res = i + 1

    return res


def compute_normal(triangle):
    normal = np.cross(triangle[1] - triangle[0], triangle[2] - triangle[0])
    normal /= np.linalg.norm(normal)
    return normal


def is_inside_convex_hull(point, hull):
    hull_center = np.mean(hull.points, axis=0)
    for simplex in hull.simplices:
        triangle = hull.points[simplex]
        normal = compute_normal(triangle)

        # 法線が外向きであることを保証
        if np.dot(normal, hull_center - triangle[0]) < 0:
            normal = -normal

        # ポイントが面の外側にある場合
        if np.dot(normal, point - triangle[0]) < 0:
            return False
    return True


def generate_3d_convex_hull_obj(path, points):
    # 凸包を計算
    hull = ConvexHull(points)

    # 凸包の面のインデックスを取得
    faces = hull.simplices.tolist()

    # objファイルとして出力
    with open(path, "w") as f:
        for v in points:
            f.write("v {} {} {}\n".format(v[0], v[1], v[2]))
        for face in faces:
            f.write("f {} {} {}\n".format(face[0] + 1, face[1] + 1, face[2] + 1))


def generate_3d_shape_obj(path, points):
    # 凸包を計算
    hull = ConvexHull(points)

    # 凸包の面のインデックスを取得
    faces = hull.simplices.tolist()

    # 凹ます中心を決める
    centroid = np.random.uniform(-1, 1, size=(3))
    while not is_inside_convex_hull(centroid, hull):
        centroid = np.random.uniform(-1, 1, size=(3))

    # 頂点を重心方向に凹ます
    for index in hull.vertices:
        points[index] += (centroid - points[index]) * min(
            np.random.exponential(0.25), 0.9
        )

    # objファイルとして出力
    with open(path, "w") as f:
        for v in points:
            f.write("v {} {} {}\n".format(v[0], v[1], v[2]))
        for face in faces:
            f.write("f {} {} {}\n".format(face[0] + 1, face[1] + 1, face[2] + 1))


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
    voxelArray = numpy_support.vtk_to_numpy(
        voxelGrid.GetOutput().GetPointData().GetScalars()
    )
    voxelArray = voxelArray.reshape((64, 64, 64))

    # Convert SDF to NumPy array
    sdfArray = numpy_support.vtk_to_numpy(
        distanceFilter.GetOutput().GetPointData().GetScalars()
    )
    sdfArray = sdfArray.reshape((64, 64, 64))

    return voxelArray, sdfArray


def convert_data(root):
    for file in glob(f"{root}*_crack.csv"):
        try:
            convert_path = os.path.splitext(file)[0] + f".npz"
            _crack = crack.read_crack(file)
            crack.save_crack_npz(convert_path, _crack)
            os.remove(file)
        except Exception as e:
            print(f"convert error: {e}({file})")


def log(log_type, log_data):
    url = "https://l2m7mz314k.execute-api.ap-northeast-1.amazonaws.com/api/logs"
    headers = {
        "Content-Type": "application/json"
    }
    data = {
        "device": socket.gethostname(),
        "task": "generate_data",
        "type": log_type,
        "data": log_data,
        "created_at": datetime.now().isoformat()
    }

    response = requests.post(url, headers=headers, data=json.dumps(data))

    if response.status_code == 200:
        print(response.json())
    else:
        print(f"Failed to post data. Status code: {response.status_code}")


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
    if start_index(parsed_args) == 0:
        log("start", {})

    start = start_index(parsed_args) % parsed_args.iter
    i = start_index(parsed_args) // parsed_args.iter
    if start == 0:
        points_num = int(np.random.uniform(10, 30, size=(1))[0])
        points = np.random.uniform(-1, 1, size=(points_num, 3))
        generate_3d_shape_obj(
            f"{parsed_args.out}{i * parsed_args.iter}_input.obj", points
        )
        voxel, sdf = obj_to_voxel_and_sdf(
            f"{parsed_args.out}{i * parsed_args.iter}_input.obj"
        )
        for l in range(parsed_args.iter):
            shape.save_shape_npz(
                f"{parsed_args.out}{i * parsed_args.iter + l}_voxel.npz", voxel
            )
            shape.save_shape_npz(
                f"{parsed_args.out}{i * parsed_args.iter + l}_sdf.npz", sdf
            )
    offset = i * parsed_args.iter
    retry = 0
    while start_index(parsed_args) < parsed_args.iter * (i + 1):
        print(f"iter:{offset}")

        if retry > 5:
                print("retry: 5")
                
                for l in range(i * parsed_args.iter, (i + 1) * parsed_args.iter):
                    for file in glob(f"{parsed_args.out}{l}_*"):
                        os.remove(file)
                        print(f"remove: {file}")

                break

        try:
            subprocess.run(
                [
                    "GenerateFractureData.exe",
                    f"{parsed_args.out}{i * parsed_args.iter}_input.obj",
                    parsed_args.param,
                    "-o",
                    parsed_args.out,
                    "-i",
                    f"1",
                    "--start",
                    f"0",
                    "--offset",
                    f"{offset}",
                ],
                timeout=parsed_args.iter * 20,
            )

            crack_path = f"{parsed_args.out}{offset}_crack.csv"
            convert_path = os.path.splitext(crack_path)[0] + f".npz"
            _crack = crack.read_crack(crack_path)
            crack.save_crack_npz(convert_path, _crack)
            os.remove(crack_path)

            subprocess.run(
                ["python", "../../thirdparty/neural-fracture/tool/crack_to_fragment.py", convert_path, "--crack-info", f"{parsed_args.out}crack_info.csv"]
            )

            fragment_path = f"{parsed_args.out}{offset}_fragment.npz"
            _fragment = fragment.read_fragment(fragment_path)
            if _fragment.max() != 2:
                print("number of fragments is not 2")
                retry += 1
                continue

            offset += 1
            retry = 0
            log("in_progress", {"index": offset})

        except Exception as e:
            print(f"{e}")
            retry += 1

    convert_data(parsed_args.out)

log("complete", {})
