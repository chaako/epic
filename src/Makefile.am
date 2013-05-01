# iMesh
-include ../fmdb/fmdbSerial/lib/iMesh-Defs.inc

# epic
bin_PROGRAMS=epic mesher solutionTransfer
CLASS_FILES = custom_mImportVTK.cc Mesh.cpp Orbit.cpp Field.cpp IntegrandContainer.cpp \
              intersect_RayTriangle.cpp distributionFunctionFromBoundary.cpp Solver.cpp \
              MinimumBasisFunction.cpp
#              VelocityAndAcceleration.cpp
epic_SOURCES = epic.cpp $(CLASS_FILES)
epic_CFLAGS = $(AM_CFLAGS) $(IMESH_INCLUDES) -I../fmdb/SCOREC_Software/SCUtil/siter/src/ \
              -I @abs_top_builddir@/include/
if USE_MPI
epic_CFLAGS += -I ../superlu/SuperLU_DIST_3.2/SRC/
else
epic_CFLAGS += -I ../superlu/SuperLU_4.3/SRC/
endif
epic_CXXFLAGS = $(AM_CXXFLAGS) $(IMESH_INCLUDES) -I../fmdb/SCOREC_Software/SCUtil/siter/src/ \
              -I @abs_top_builddir@/include/ -I @abs_top_builddir@/include/vtk-5.10/ -Wno-deprecated
if USE_MPI
epic_CXXFLAGS += -I ../superlu/SuperLU_DIST_3.2/SRC/
else
epic_CXXFLAGS += -I ../superlu/SuperLU_4.3/SRC/
endif
epic_LDFLAGS = $(AM_LDFLAGS) $(IMESH_LIBS) -lMeshAdapt -lrt \
			-L @abs_top_builddir@/lib/ -L @abs_top_builddir@/lib/vtk-5.10/ \
			-lLSDyna -lMapReduceMPI -lvtkalglib -lvtkCharts -lvtkCommon -lvtkDICOMParser \
			-lvtkexoIIc -lvtkexpat -lvtkFiltering -lvtkfreetype -lvtkftgl -lvtkGenericFiltering \
			-lvtkGeovis -lvtkGraphics -lvtkhdf5 -lvtkhdf5_hl -lvtkHybrid -lvtkImaging \
			-lvtkInfovis -lvtkIO -lvtkjpeg -lvtklibxml2 -lvtkmetaio -lvtkNetCDF -lvtkNetCDF_cxx \
			-lvtkpng -lvtkproj4 -lvtkRendering -lvtksqlite -lvtksys -lvtktiff -lvtkverdict \
			-lvtkViews -lvtkVolumeRendering -lvtkWidgets -lvtkzlib \
			-lnglib
if USE_MPI
epic_CXXFLAGS += -DHAVE_MPI
else
epic_LDFLAGS += -lmpistubs
endif

# Note that linking order matters: if A has a reference to B, then the order should be A B
epic_LDADD   = ../cubature/libcubature.a ../cuba/Cuba-3.0/libcuba.a
if USE_MPI
epic_LDADD   += ../superlu/SuperLU_DIST_3.2/lib/libsuperlu_dist_3.2.a ../superlu/SuperLU_DIST_3.2/lib/libblas.a
else
epic_LDADD   += ../superlu/SuperLU_4.3/lib/libsuperlu_4.3.a ../superlu/SuperLU_4.3/lib/libblas.a
endif

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
