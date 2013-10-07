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
import matplotlib.tri as tri
import matplotlib.cm as cm
import IPython.core.display as IPdisp

# <codecell>

file_names = [
              "/Users/chaako/transferEPS/elSphG1.0EvalG1.0N4000_13_00.vtu",
              "/Users/chaako/transferEPS/elSphG0.0EvalG0.0N4000_13_00.vtu",
              "/Users/chaako/transferEPS/elSphG-1.0EvalG-1.0N4000_13_00.vtu",
              "/Users/chaako/transferEPS/elSphG1.0EvalG1.0N4000_13_00.vtu",
              "/Users/chaako/transferEPS/elSphG1.0EvalG0.0N4000_13_00.vtu",
              "/Users/chaako/transferEPS/elSphG-1.0EvalG-1.0N4000_13_00.vtu",
              "/Users/chaako/transferEPS/elSphG-1.0EvalG0.0N4000_13_00.vtu",
              "/Users/chaako/transferEPS/elSphG0.0EvalG1.0N4000_13_00.vtu",
              "/Users/chaako/transferEPS/elSphG0.0EvalG0.0N4000_13_00.vtu",
              "/Users/chaako/transferEPS/elSphG0.0EvalG-1.0N4000_13_00.vtu",
              "/Users/chaako/transferEPS/elSphG0.0EvalG0.0N4000_13_00.vtu"
              ]
vtk_readers = []
for file_name in file_names:
    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(file_name)
    reader.Update()
    vtk_readers.append(reader.GetOutput())
vtk_mesh = vtk_readers[0]

# <codecell>

for reader in vtk_readers:
    print reader.GetNumberOfPoints()

# <codecell>

cell_code_name = "cell_code";
cell_codes = vtk_mesh.GetCellData().GetArray(cell_code_name);

# <codecell>

coordinates = []
number_of_points = vtk_mesh.GetNumberOfPoints()
vtk_points = vtk_mesh.GetPoints();
for i in range(0,number_of_points):
    coordinates.append(vtk_points.GetPoint(i))
coordinates_np = np.array(coordinates)

# <codecell>

coordinates_np_tranpose = np.transpose(coordinates_np)
Xs = coordinates_np_tranpose[0]
Ys = coordinates_np_tranpose[1]
Zs = coordinates_np_tranpose[2]

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
        #array.append(vtk_array.GetTuple(i))
        # TODO: generalize this for N dimensions
        v1 = vtk_array.GetComponent(i,0)
        v2 = vtk_array.GetComponent(i,1)
        v3 = vtk_array.GetComponent(i,2)
        array.append((v1,v2,v3))
    return array

# <codecell>

ion_densities = []
reference_densities = []
ion_velocities = []
for vtk_reader in vtk_readers:
    ion_densities.append(getValueArrayFromVtk(vtk_reader, "ionSurfaceDensity"))
    reference_densities.append(getValueArrayFromVtk(vtk_reader, "surfaceReferenceDensity"))
    ion_velocities.append(getTupleArrayFromVtk(vtk_reader, "ionSurfaceVelocity"))

# <codecell>

v_E = 0.1
v_n = 0.1

#innerCoordinates = [];
zetas = []
xis = []
alphas = []
thetas = []
reference_fluxes = []
reference_fluxes_scaled = []
reference_fluxes_grad = []
reference_fluxes_grad_scaled = []
for i in range(0,vtk_mesh.GetNumberOfPoints()):
    c = coordinates[i]
    if (c[0]*c[0]+c[1]*c[1]+c[2]*c[2]<1.5):
        #innerCoordinates.append(coordinates[i])
        zeta = -atan2(c[0],c[2])
        zetas.append(zeta)
        xi = -asin(c[1])
        xis.append(xi)
        alpha = abs(asin(c[2]))
        alphas.append(alpha)
        theta = -pi/2.-atan2(-c[1],abs(c[2]))
        thetas.append(theta)
        exponent = -1.0+v_E/(-0.0000001+tan(theta))
        reference_fluxes.append(exp(min(10.,max(-10.,exponent)))*sin(alpha))
        reference_fluxes_scaled.append(exp(min(10.,max(-10.,exponent)))*sin(alpha)*reference_densities[0][i])
        exponent = -1.0+(v_E+v_n*(1.-(1.0-sin(alpha))/(1.0+sin(alpha))))/(0.0000001+tan(theta))
        reference_fluxes_grad.append(exp(min(10.,max(-10.,exponent)))*sin(alpha))
        reference_fluxes_grad_scaled.append(exp(min(10.,max(-10.,exponent)))*sin(alpha)*reference_densities[0][i])
    else:
        xis.append(0.)
        zetas.append(0.)
        alphas.append(0.)
        thetas.append(0.)
        reference_fluxes.append(0.)
        reference_fluxes_scaled.append(0.)
        reference_fluxes_grad.append(0.)
        reference_fluxes_grad_scaled.append(0.)
flux_arrays = [];
for ion_density, reference_density, ion_velocity in zip(ion_densities, reference_densities, ion_velocities):
    fluxes = []
    for i in range(0,vtk_mesh.GetNumberOfPoints()):
        c = coordinates[i]
        v = ion_velocity[i] #ies[j][i]
        normal_velocity = c[0]*v[0]+c[1]*v[1]+c[2]*v[2]
        #theta = thetas[i]
        #alpha = alphas[i]
        if (c[0]*c[0]+c[1]*c[1]+c[2]*c[2]<1.5):
            fluxes.append(-normal_velocity*ion_density[i])
            #exponent = -1.0+(v_E+v_n*(1.-(1.0-sin(alpha))/(1.0+sin(alpha))))/(0.0000001+tan(theta))
            #reference_fluxes.append(exp(min(10.,max(-10.,exponent)))*sin(alpha))
        else:
            fluxes.append(0.)
            #reference_fluxes.append(0.)
    flux_arrays.append(fluxes)

# <codecell>

triangles = []
for i in range(0,vtk_mesh.GetNumberOfCells()):
    if (vtk_mesh.GetCellType(i) == VTK_TRIANGLE):
        if (cell_codes.GetValue(i) == 4):
            vtk_ids = vtkIdList()
            vtk_mesh.GetCellPoints(i,vtk_ids)
            ids = (vtk_ids.GetId(0),vtk_ids.GetId(1),vtk_ids.GetId(2))
            skip_triangle = False
            for j in range(0,3):
                for jj in range(j,3):
                    if (abs(zetas[ids[j]]-zetas[ids[jj]])>4.):
                        skip_triangle=True
            if (not skip_triangle):
                triangles.append(ids)

# <codecell>

zs = []
for fluxes in flux_arrays:
    zs.append(np.array(fluxes).flatten())
x = np.array(zetas).flatten()
y = np.array(xis).flatten()
triang = tri.Triangulation(x, y, triangles)
# sqrt(2) from c_s difference (?)
zs.append(sqrt(2.)*np.array(reference_fluxes_grad_scaled).flatten())
zs.append(sqrt(2.)*np.array(reference_fluxes).flatten())
zs.append(sqrt(2.)*np.array(reference_fluxes_scaled).flatten())
zs.append(sqrt(2.)*np.array(reference_fluxes_grad).flatten())

# <codecell>

angles = []
angles.append(np.array(zetas).flatten())
angles.append(np.array(xis).flatten())
angles.append(np.array(thetas).flatten())
angles.append(np.array(alphas).flatten())
#indices = Xs*Xs+Ys*Ys+Zs*Zs>1.5
#alphas = np.zeros(number_of_points)
#thetas = np.zeros(number_of_points)
#print thetas, indices
#print thetas[indices]
#alphas[indices] = np.absolute(np.array(xis[indices]))
#thetas[indices] = -pi/2.-np.arctan2(-ys[indices],zs[indices])
#reference_fluxes = np.exp(-1.0+(v_E+v_n*(1.-(1.0-np.sin(alphas))/(1.0+np.sin(alphas))))/np.tan(thetas))*np.sin(alphas)

# <codecell>

plt.clf()

#collection = matplotlib.collections.TriMesh(triang)
#collection.set_cmap(plt.cm.rainbow)
#collection.set_array(z)
#collection.autoscale_None()

#plt.gca().set_aspect('equal')
#plt.tripcolor(triang, z, shading='gouraud', cmap=plt.cm.rainbow)
#plt.colorbar()
#plt.scatter(x,y)
#plt.title('tripcolor of Delaunay triangulation, gouraud shading')

#fig = plt.figure(1, figsize=(12,6))
#fig = plt.figure(1)
#fig.clf()
#ax = fig.add_subplot(111, projection = 'mollweide')
#ax = plt.subplot(111, projection = 'mollweide')
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(6,9))#, subplot_kw=dict(projection='mollweide'))
ims = []
for ax,z in zip(axes.flat,zs):
    im = ax.tripcolor(triang, z, vmin=0., vmax=0.8, shading='flat', cmap=plt.cm.rainbow)
    ims.append(im)
    cont = ax.tricontour(triang, z, vmin=0., vmax=0.8, colors='k')
    #plt.colorbar(im, ax=ax, orientation='horizontal')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(ims[0], cax=cbar_ax)

#fig, ax = plt.subplots()
#p = ax.pcolor(X/(2*pi), Y/(2*pi), Z, cmap=cm.RdBu, vmin=abs(Z).min(), vmax=abs(Z).max())
#cb = fig.colorbar(p)

#im = imshow(Z, cmap=cm.RdBu, vmin=abs(Z).min(), vmax=abs(Z).max(), extent=[0, 1, 0, 1])
#im.set_interpolation('bilinear')
#cb = fig.colorbar(im)

#fig = plt.figure()
#axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # lower, bottom, width, height (range 0 to 1)
#axes.plot(x, y, 'r')

#axes.set_xlabel('x')
#axes.set_ylabel('y')
#axes.set_title('title');

#plt.savefig("data.svg")
#IPdisp.SVG(filename="data.svg")
plt.savefig("fluxesFlat.png")
IPdisp.Image(filename="fluxesFlat.png")

# <codecell>

plt.clf()

fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12,6))#, subplot_kw=dict(projection='mollweide'))
ims = []
for ax,z in zip(axes.flat,angles):
    im = ax.tripcolor(triang, z, shading='faceted', cmap=plt.cm.rainbow)
    ims.append(im)
    cont = ax.tricontour(triang, z, colors='k')
    #plt.colorbar(im, ax=ax, orientation='horizontal')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(ims[0], cax=cbar_ax)

#plt.savefig("data.svg")
#IPdisp.SVG(filename="data.svg")
plt.savefig("data.png")
IPdisp.Image(filename="data.png")

# <codecell>


