# iMesh
-include ../fmdb/fmdbSerial/lib/iMesh-Defs.inc

# epic
bin_PROGRAMS=epic mesher solutionTransfer
CLASS_FILES = custom_mImportVTK.cc Mesh.cpp Orbit.cpp Field.cpp IntegrandContainer.cpp \
              intersect_RayTriangle.cpp distributionFunctionFromBoundary.cpp
#              VelocityAndAcceleration.cpp
epic_SOURCES = epic.cpp $(CLASS_FILES)
epic_CFLAGS = $(AM_CFLAGS) $(IMESH_INCLUDES) -I../fmdb/SCOREC_Software/SCUtil/siter/src/ \
              -I ../cubature/cubature-20101018/ -I ../cuba/Cuba-3.0/ -I ../superlu/superLU_4.3/SRC/
epic_CXXFLAGS = $(AM_CXXFLAGS) $(IMESH_INCLUDES) -I../fmdb/SCOREC_Software/SCUtil/siter/src/ \
              -I ../eigen/eigen-eigen-43d9075b23ef/ -I ../cubature/cubature-20101018/ \
              -I ../cuba/Cuba-3.0/ -I ../netgen/opt/netgen/include/  -I ../superlu/SuperLU_4.3/SRC/ \
              -I ../boost/boost_1_49_0/ -I ../odeint/headmyshoulder-odeint-v2-47dbdac/ \
              -I ../vtk/VTK/install/include/vtk-5.10/ -Wno-deprecated
epic_LDFLAGS = $(AM_LDFLAGS) $(IMESH_LIBS) -lMeshAdapt -lrt \
			-L ../vtk/VTK/install/lib/vtk-5.10/ \
			-lLSDyna -lMapReduceMPI -lvtkalglib -lvtkCharts -lvtkCommon -lvtkDICOMParser \
			-lvtkexoIIc -lvtkexpat -lvtkFiltering -lvtkfreetype -lvtkftgl -lvtkGenericFiltering \
			-lvtkGeovis -lvtkGraphics -lvtkhdf5 -lvtkhdf5_hl -lvtkHybrid -lvtkImaging \
			-lvtkInfovis -lvtkIO -lvtkjpeg -lvtklibxml2 -lvtkmetaio -lvtkNetCDF -lvtkNetCDF_cxx \
			-lvtkpng -lvtkproj4 -lvtkRendering -lvtksqlite -lvtksys -lvtktiff -lvtkverdict \
			-lvtkViews -lvtkVolumeRendering -lvtkWidgets -lvtkzlib \
			-L ../netgen/opt/netgen/lib/ -lnglib
if USE_MPI
epic_CXXFLAGS += -DHAVE_MPI
else
epic_LDFLAGS += -lmpistubs
endif

# Note that linking order matters: if A has a reference to B, then the order should be A B
epic_LDADD   = ../cubature/libcubature.a ../cuba/Cuba-3.0/libcuba.a
epic_LDADD   += ../superlu/SuperLU_4.3/lib/libsuperlu_4.3.a ../superlu/SuperLU_4.3/lib/libblas.a

#mesher
mesher_SOURCES = mesher.cpp $(CLASS_FILES) playground.cpp SurfaceMesh.cpp
mesher_CFLAGS = $(epic_CFLAGS) -DMESHER
mesher_CXXFLAGS = $(epic_CXXFLAGS) -DMESHER
mesher_LDFLAGS = $(epic_LDFLAGS)
mesher_LDADD = $(epic_LDADD)

#solutionTransfer
../sceptic3D/.dummy:
	cd .. && git submodule init && git submodule update

../sceptic3D/interptool: ../sceptic3D/.dummy ../sceptic3D/interpolateFieldsSceptic3D.f
	cd ../sceptic3D && make interptool

../sceptic3D/libsceptic3Dhdf.a: ../sceptic3D/interptool
	cd ../sceptic3D && make libsceptic3Dhdf.a

DIRHDF = $(PWD)/../sceptic3D/hdf5-1.8.4
LIBHDF = -L$(DIRHDF)/lib -lhdf5hl_fortran -lhdf5_hl \
          -lhdf5_fortran -lhdf5 -lz -lm -Wl,-rpath -Wl,$(DIRHDF)/lib

solutionTransfer_SOURCES = solutionTransfer.cpp $(CLASS_FILES)
solutionTransfer_CFLAGS = $(epic_CFLAGS)
solutionTransfer_CXXFLAGS = $(epic_CXXFLAGS)
solutionTransfer_LDFLAGS = $(epic_LDFLAGS) $(LIBHDF) -lgfortran
solutionTransfer_LDADD = $(epic_LDADD) ../sceptic3D/libsceptic3Dhdf.a
#interpolateFieldsSceptic3D.o
#		../sceptic3D/advancing.o ../sceptic3D/readhdf.o ../sceptic3D/chargefield.o
