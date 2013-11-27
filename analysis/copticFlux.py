# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from vtk import *
from math import *
import numpy as np
import matplotlib
import glob
import scipy.stats as stats
import scipy.optimize as optimize
import scipy.signal as signal
import scipy as sp
# Save as file and display instead of %pylab inline
#matplotlib.use('svg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as tri
import matplotlib.cm as cm
import IPython.core.display as IPdisp

# <codecell>

#folder = "/home/chaako/meshtest/fromLoki/ctzhou/v0.1"
folder = "/scratch/chaako/fromLoki/unsteadyness013"
file_names = glob.glob(folder + "/*flx*.vtk")
file_names.sort()
print len(file_names), file_names[0]

# <codecell>

vtk_readers = []
for file_name in file_names:
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(file_name)
    reader.Update()
    vtk_readers.append(reader.GetOutput())
vtk_mesh = vtk_readers[0]

# <codecell>

#for reader in vtk_readers:
#    print reader.GetNumberOfPoints()
#print vtk_readers[0]

# <codecell>

coordinates = []
number_of_points = vtk_mesh.GetNumberOfPoints()
vtk_points = vtk_mesh.GetPoints();
for i in range(0,number_of_points):
    coordinates.append(vtk_points.GetPoint(i))
coordinates_np = np.array(coordinates)

# <codecell>

coordinates_np_tranpose = np.transpose(coordinates_np)
x = coordinates_np_tranpose[0]
y = coordinates_np_tranpose[1]
z = coordinates_np_tranpose[2]

# <codecell>

def getValueArrayFromVtk(vtk_reader, array_name):
    vtk_array = vtk_reader.GetPointData().GetArray(array_name)
    array = []
    for i in range(0,vtk_reader.GetNumberOfPoints()):
        array.append(vtk_array.GetValue(i))
    return array

# <codecell>

fluxes = []
for vtk_reader in vtk_readers:
    fluxes.append(getValueArrayFromVtk(vtk_reader, "flux"))
fluxes_np = np.squeeze(np.array(fluxes))

# <codecell>

plt.clf()
plt.plot(fluxes_np)
#plt.savefig("data.svg")
#IPdisp.SVG(filename="data.svg")
plt.savefig("data.png")
IPdisp.Image(filename="data.png")

# <codecell>

def gaussian_smooth(data_array, sigma):
    mu = 0
    bins = np.arange(-4*sigma, 4*sigma, 1, dtype=np.float32)
    gaussian = 1/(sigma * np.sqrt(2 * np.pi))*np.exp( - (bins - mu)**2 / (2 * sigma**2) )
    return signal.convolve(data_array,gaussian,mode='same')

# <codecell>

print fluxes_np.shape

# <codecell>

fluxes_np_t = fluxes_np.transpose()
fluxes_np_smoothed_t = np.zeros(fluxes_np_t.shape, dtype=np.float32)
for j in range(0,fluxes_np_t.shape[0]):
    fluxes_np_smoothed_t[j][:] = gaussian_smooth(fluxes_np_t[j],20)
fluxes_np_smoothed = fluxes_np_smoothed_t.transpose()

# <codecell>

plt.clf()
plt.plot(fluxes_np_smoothed)
#plt.savefig("data.svg")
#IPdisp.SVG(filename="data.svg")
plt.savefig("data.png")
IPdisp.Image(filename="data.png")

# <codecell>

print coordinates

# <codecell>

plt.clf()
fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(26,11), subplot_kw=dict(projection='3d'))
i=0
scts = []
for ax in axes.flat:
    if i<fluxes_np.shape[0]:
        sct = ax.scatter(x,y,z,c=fluxes[i],vmin=0,vmax=max(max(fluxes)))
        i=i+1
    else:
        sct = ax.scatter(x,y,z,vmin=0,c=fluxes[0],vmax=max(max(fluxes)))
        # Note that colorbar seems to screw up transparency/shading
        sct_ref = sct
        sct = ax.scatter(x,y,z,vmin=0,vmax=max(max(fluxes)))
    scts.append(sct)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    #ax.vew_init(0,0)
#sct_ref = scts[len(scts)-1]
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(sct_ref, cax=cbar_ax)

#plt.savefig("data.svg")
#IPdisp.SVG(filename="data.svg")
plt.savefig("data.png")
IPdisp.Image(filename="data.png")

# <codecell>

print fluxes_np.shape[1]

# <codecell>

print 'hello'

# <codecell>


