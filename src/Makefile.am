# iMesh
include ../fmdb/fmdbSerial/lib/iMesh-Defs.inc

# epic
bin_PROGRAMS=epic
epic_SOURCES = epic.cpp custom_mImportVTK.cc Mesh.cpp Orbit.cpp Field.cpp IntegrandContainer.cpp\
              ../cubature/cubature-20101018/cubature.c
epic_CFLAGS = $(AM_CFLAGS) $(IMESH_INCLUDES) -I../fmdb/SCOREC_Software/SCUtil/siter/src/ \
              -I ../eigen/eigen-eigen-13a11181fc5a/ -I ../cubature/cubature-20101018/ \
               -I ../cuba/Cuba-3.0/
epic_CXXFLAGS = $(AM_CXXFLAGS) $(IMESH_INCLUDES) -I../fmdb/SCOREC_Software/SCUtil/siter/src/ \
              -I ../eigen/eigen-eigen-13a11181fc5a/ -I ../cubature/cubature-20101018/ \
               -I ../cuba/Cuba-3.0/
epic_LDFLAGS = $(AM_LDFLAGS) $(IMESH_LIBS) -L../cuba/Cuba-3.0 -lcuba #Vegas.o
