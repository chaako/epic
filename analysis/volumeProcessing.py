# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from vtk import *
from math import *
import numpy as np
import matplotlib
# Save as file and display instead of %pylab inline
#matplotlib.use('svg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as tri
import matplotlib.cm as cm
import IPython.core.display as IPdisp

# <codecell>

file_names = ["/home/chaako/meshtest/fromLoki/densityGradient18/densityGradient/spheres_iter07.vtk",
              "/home/chaako/meshtest/fromLoki/densityGradient31/densityGradient/spheres_iter01.vtk"]
vtk_readers = []
for file_name in file_names:
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(file_name)
    reader.Update()
    vtk_readers.append(reader.GetOutput())
vtk_mesh = vtk_readers[0]

# <codecell>

coordinates = []
number_of_points = vtk_mesh.GetNumberOfPoints()
vtk_points = vtk_mesh.GetPoints();
for i in range(0,number_of_points):
    coordinates.append(vtk_points.GetPoint(i))
coordinates_np = np.array(coordinates)

# <codecell>

def getValueArrayFromVtk(vtk_reader, array_name):
    vtk_array = vtk_reader.GetPointData().GetArray(array_name)
    array = []
    for i in range(0,vtk_reader.GetNumberOfPoints()):
        array.append(vtk_array.GetValue(i))
    return array

# <codecell>

def getTupleArrayFromVtk(vtk_reader, array_name):
    vtk_array = vtk_reader.GetPointData().GetArray(array_name)
    array = []
    for i in range(0,vtk_reader.GetNumberOfPoints()):
        array.append(vtk_array.GetTuple(i))
    return array

# <codecell>

def getTupleArrayFromVtkXYZ(vtk_reader, array_name):
    vtk_arrayX = getTupleArrayFromVtk(vtk_reader, array_name+"X");
    vtk_arrayY = getTupleArrayFromVtk(vtk_reader, array_name+"Y");
    vtk_arrayZ = getTupleArrayFromVtk(vtk_reader, array_name+"Z");
    array = []
    for i in range(0,vtk_reader.GetNumberOfPoints()):
        array.append((vtk_arrayX[i][0],vtk_arrayY[i][0],vtk_arrayZ[i][0]))
    return array

# <codecell>

ion_densities = []
potentials = []
ion_velocities = []
for vtk_reader in vtk_readers:
    ion_densities.append(getValueArrayFromVtk(vtk_reader, "ionDensity"))
    potentials.append(getValueArrayFromVtk(vtk_reader, "potential"))
    ion_velocities.append(getTupleArrayFromVtkXYZ(vtk_reader, "ionVelocity"))
ion_densities_np = np.array(ion_densities)
potentials_np = np.array(potentials)
ion_velocities_np = np.array(ion_velocities)

# <codecell>

coordinates_np_tranpose = np.transpose(coordinates_np)
#indices = range(0,number_of_points)
indices = np.logical_and(abs(potentials_np[1]-potentials_np[0])==0.0, potentials_np[0]!=0.)
#indices = potentials_np[0]==0.
x = coordinates_np_tranpose[0][indices]
y = coordinates_np_tranpose[1][indices]
z = coordinates_np_tranpose[2][indices]

# <codecell>

plt.clf()

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(24,12), subplot_kw=dict(projection='3d'))
#for ax in axes.flat:
#    ax.scatter(x,y,z)
sct = axes[0].scatter(x,y,z)
axes[0].view_init(0,0)

#plt.savefig("data.svg")
#IPdisp.SVG(filename="data.svg")
plt.savefig("data.png")
IPdisp.Image(filename="data.png")

# <codecell>


