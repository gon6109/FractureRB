from scipy.spatial import ConvexHull
import numpy as np

# ランダムな3D座標を生成
points = np.random.uniform(-1, 1, size=(30, 3))

# 凸包を計算
hull = ConvexHull(points)

# 凸包の頂点座標を取得
vertices = hull.vertices.tolist()

# 凸包の面のインデックスを取得
faces = hull.simplices.tolist()

# 出力ファイル名
output_file = "convex_hull.obj"

# objファイルとして出力
with open(output_file, "w") as f:
    for v in points:
        f.write("v {} {} {}\n".format(v[0], v[1], v[2]))
    for face in faces:
        f.write("f {} {} {}\n".format(face[0] + 1, face[1] + 1, face[2] + 1))
