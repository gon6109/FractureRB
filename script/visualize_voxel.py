import vtk
import numpy as np
from mayavi import mlab
from mayavi.mlab import volume_slice

# OBJファイルからメッシュを読み込む
reader = vtk.vtkOBJReader()
reader.SetFileName("E:\\fracture-data\\test\\0_input.obj")

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

# Convert SDF to NumPy array
sdfArray = vtk.util.numpy_support.vtk_to_numpy(distanceFilter.GetOutput().GetPointData().GetScalars())
sdfArray = sdfArray.reshape((64, 64, 64))

volume_slice(sdfArray)
mlab.show()