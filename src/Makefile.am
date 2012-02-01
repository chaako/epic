# iMesh
include ../fmdb/fmdbSerial/lib/iMesh-Defs.inc

# epic
bin_PROGRAMS=epic
epic_SOURCES=epic.cpp custom_mImportVTK.cc
epic_CFLAGS = $(AM_CFLAGS) $(IMESH_INCLUDES) -I../fmdb/SCOREC_Software/SCUtil/siter/src/ \
              -I ../eigen/eigen-eigen-13a11181fc5a/
epic_CXXFLAGS = $(AM_CXXFLAGS) $(IMESH_INCLUDES) -I../fmdb/SCOREC_Software/SCUtil/siter/src/ \
              -I ../eigen/eigen-eigen-13a11181fc5a/
epic_LDFLAGS = $(AM_LDFLAGS) $(IMESH_LIBS)
