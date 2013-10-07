# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from vtk import *
from math import *
import numpy as np
import matplotlib
import glob
# Save as file and display instead of %pylab inline
#matplotlib.use('svg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as tri
import matplotlib.cm as cm
import IPython.core.display as IPdisp

# <codecell>

#file_names = ["/home/chaako/epic/src/testCase00.vtu"]
#file_names = ["/home/chaako/epic/src/spheres_scan1000_narrow00.vtu",
#              "/home/chaako/epic/src/spheres_scan4000_narrow00.vtu",
#              "/home/chaako/epic/src/spheres_scan1000_00.vtu"]
#file_names = ["/home/chaako/epic/src/runsTesting/dG16_01/spheres5_scan4000_narrow00.vtu",
#              "/home/chaako/epic/src/spheres6_scan4000_narrow00.vtu",
#              "/home/chaako/epic/src/spheres_scan4000_narrow00.vtu"]
#file_names = ["/home/chaako/epic/src/runsTesting/dG22_10/spheres00_scan1000_narrow00.vtu",
#              "/home/chaako/epic/src/runsTesting/dG22_11/spheres01_scan1000_narrow00.vtu",
#             "/home/chaako/epic/src/runsTesting/dG22_12/spheres02_scan1000_narrow00.vtu",
#              "/home/chaako/epic/src/runsTesting/dG22_13/spheres03_scan1000_narrow00.vtu",
#              "/home/chaako/epic/src/runsTesting/dG22_14/spheres04_scan1000_narrow00.vtu",
#              "/home/chaako/epic/src/runsTesting/dG22_15/spheres05_scan1000_narrow00.vtu",
#              "/home/chaako/epic/src/runsTesting/dG22_16/spheres06_scan1000_narrow00.vtu"]
#file_names = ["/home/chaako/meshtest/fromLoki/densityGradient40/densityGradient/spheres_iter00.vtu"]
#file_names = ["/home/chaako/meshtest/fromLoki/iterationTest23/iterationTest/spheres_iter00.vtu"]
folder = "/home/chaako/meshtest/fromLoki/densityGradient97/densityGradient"
file_names = glob.glob(folder + "/*.vtu")
file_names.sort()
print file_names

# <codecell>

vtk_readers = []
for file_name in file_names:
    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(file_name)
    reader.Update()
    vtk_readers.append(reader.GetOutput())
vtk_mesh = vtk_readers[0]

# <codecell>

#target_point = np.array([-0.160938, 0.752753, 2.82129])
#target_point = np.array([-0.165899, 0.990274, 6.2703])
#target_point = np.array([-0.387776, 0.357946, -2.1117])
#target_point = np.array([0.252365, 0.769792, 5.80838])
#target_point = np.array([0.801777, 0.72545, 13.814])
#target_point = np.array([0.131742, 0.528943, -3.84111])
#target_point = np.array([0.406818, 1.12104, 14.7508])
#target_point = np.array([-0.109546, 0.794057, -4.72881])
#target_point = np.array([0.0501395, 1.3243, -7.43127])
target_point = np.array([-0.19836, 1.27436, 3.76091])
closeness_threshold = 0.001

# <codecell>

coordinates = []
target_point_index = -1
number_of_points = vtk_mesh.GetNumberOfPoints()
vtk_points = vtk_mesh.GetPoints();
for i in range(0,number_of_points):
    coordinates.append(vtk_points.GetPoint(i))
    #print np.linalg.norm(np.array(coordinates[i])-target_point)
    if (np.linalg.norm(np.array(coordinates[i])-target_point)<closeness_threshold):
        target_point_index = i
coordinates_np = np.array(coordinates)
print target_point_index, number_of_points

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
#ion_densities_scan = []
potentials = []
previousPotentials = []
#potentials_scan = []
ion_velocities = []
for vtk_reader in vtk_readers:
    ion_densities.append(getTupleArrayFromVtk(vtk_reader, "ionDensity"))
#    ion_densities_scan.append(getTupleArrayFromVtk(vtk_reader, "ionDensityScan"))
    potentials.append(getTupleArrayFromVtk(vtk_reader, "potential"))
    previousPotentials.append(getTupleArrayFromVtk(vtk_reader, "previousPotential"))
#    potentials_scan.append(getTupleArrayFromVtk(vtk_reader, "potentialScan"))
    ion_velocities.append(getTupleArrayFromVtk(vtk_reader, "ionVelocity"))
ion_densities_np = np.squeeze(np.array(ion_densities))
potentials_np = np.squeeze(np.array(potentials))
previousPotentials_np = np.squeeze(np.array(previousPotentials))
ion_velocities_np = np.squeeze(np.array(ion_velocities))

# <codecell>

plt.clf()

fig, axes = plt.subplots(nrows=1, ncols=7, figsize=(24,4))
i=0
j=0
startIter = 0
endIter = len(file_names)
for ax in axes:
    #x = range(0,len(file_names))
    x = range(startIter,endIter)
    y = potentials_np.T[target_point_index+j][startIter:endIter]
    sct = ax.scatter(x,y)
    j = j+100

#plt.savefig("data.svg")
#IPdisp.SVG(filename="data.svg")
plt.savefig("data.png")
IPdisp.Image(filename="data.png")

# <codecell>

plt.clf()

fig, axes = plt.subplots(nrows=1, ncols=7, figsize=(24,4))
i=0
j=0
for ax in axes:
    #x = range(0,len(file_names))
    x = range(startIter,endIter)
    y = ion_densities_np.T[target_point_index+j][startIter:endIter]
    sct = ax.scatter(x,y)
    #y = np.exp(previousPotentials_np.T[target_point_index+j][startIter:endIter])
    #sct = ax.scatter(x,y,color='yellow')
    y = np.exp(potentials_np.T[target_point_index+j][startIter:endIter])
    sct = ax.scatter(x,y,color='green')
    #y = np.exp(3./4.*previousPotentials_np.T[target_point_index+j][startIter:endIter] + 1./4.*np.log(ion_densities_np.T[target_point_index+j][startIter-1:endIter-1]))
    #sct = ax.scatter(x,y,color='red')
    j = j+100

#plt.savefig("data.svg")
#IPdisp.SVG(filename="data.svg")
plt.savefig("data.png")
IPdisp.Image(filename="data.png")

# <codecell>

plt.clf()

fig, axes = plt.subplots(nrows=1, ncols=7, figsize=(24,4))
i=0
j=0
for ax in axes:
    if (i<len(file_names)):
        print i, target_point_index+j
        x = potentials_scan[i][target_point_index+j]
        y = ion_densities_scan[i][target_point_index+j]
        sct = ax.scatter(x,y)
        #print y[0], y[1]
        y = np.exp(potentials_scan[i][target_point_index+j])
        sct = ax.scatter(x,y,color='red')
        x = potentials_scan[i][target_point_index+j][0]
        y = ion_densities[i][target_point_index+j]
        sct = ax.scatter(x,y,color='green')
        ax.set_xlim([-0.8,0.2])
        ax.set_ylim([0,1.2])
        #print y[0]
        #i = i+1
        j = j+100

#plt.savefig("data.svg")
#IPdisp.SVG(filename="data.svg")
plt.savefig("data.png")
IPdisp.Image(filename="data.png")

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


# <codecell>


