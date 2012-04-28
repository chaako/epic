# iMesh
-include ../fmdb/fmdbSerial/lib/iMesh-Defs.inc

# epic
bin_PROGRAMS=epic
epic_SOURCES = epic.cpp custom_mImportVTK.cc Mesh.cpp Orbit.cpp Field.cpp IntegrandContainer.cpp \
              intersect_RayTriangle.cpp distributionFunctionFromBoundary.cpp
epic_CFLAGS = $(AM_CFLAGS) $(IMESH_INCLUDES) -I../fmdb/SCOREC_Software/SCUtil/siter/src/ \
              -I ../eigen/eigen-eigen-13a11181fc5a/ -I ../cubature/cubature-20101018/
epic_CXXFLAGS = $(AM_CXXFLAGS) $(IMESH_INCLUDES) -I../fmdb/SCOREC_Software/SCUtil/siter/src/ \
              -I ../eigen/eigen-eigen-13a11181fc5a/ -I ../cubature/cubature-20101018/ \
              -I ../boost/boost_1_49_0/
if USE_MPI
epic_CXXFLAGS += -DHAVE_MPI
endif
epic_LDFLAGS = $(AM_LDFLAGS) $(IMESH_LIBS) -lrt
epic_LDADD   = ../cubature/libcubature.a
