# iMesh
-include ../fmdb/fmdbSerial/lib/iMesh-Defs.inc

# epic
bin_PROGRAMS=epic
epic_SOURCES = epic.cpp custom_mImportVTK.cc Mesh.cpp Orbit.cpp Field.cpp IntegrandContainer.cpp \
              intersect_RayTriangle.cpp distributionFunctionFromBoundary.cpp playground.cpp
#              VelocityAndAcceleration.cpp
epic_CFLAGS = $(AM_CFLAGS) $(IMESH_INCLUDES) -I../fmdb/SCOREC_Software/SCUtil/siter/src/ \
              -I ../cubature/cubature-20101018/ -I ../cuba/Cuba-3.0/
epic_CXXFLAGS = $(AM_CXXFLAGS) $(IMESH_INCLUDES) -I../fmdb/SCOREC_Software/SCUtil/siter/src/ \
              -I ../eigen/eigen-eigen-43d9075b23ef/ -I ../cubature/cubature-20101018/ \
              -I ../cuba/Cuba-3.0/ -I ../netgen/opt/netgen/include/ \
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

epic_LDADD   = ../cubature/libcubature.a ../cuba/Cuba-3.0/libcuba.a
