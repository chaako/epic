/*
 * types.h
 *
 *  Created on: Jun 18, 2012
 *      Author: chaako
 */

#ifndef TYPES_H_
#define TYPES_H_

#ifdef MESHER
#undef HAVE_MPI
#endif

using namespace std;

//const double VOLUME_TOLERANCE=1.e-8;
const double VOLUME_TOLERANCE=0.;
const double LENGTH_TOLERANCE=1.e-10;
const double DELTA_LENGTH=1.e-5;
//const double DELTA_LENGTH=1.e-5;
const double NODE_DISTANCE_THRESHOLD=1.e-3;
const double SMALL_VELOCITY=1.e-1;
const double SMALL_TIME=1.e-8;
const double SMALL_DENOMINATOR=1.e-5;
const double SMALL_POTENTIAL_CHANGE=1.e-2;
const double SMALL_DENSITY_CHANGE=1.e-2;
const double SMALL_DENSITY_DERIVATIVE=1.e-1;
const double LARGE_POTENTIAL_CHANGE=5.e-1;

const int WORKTAG=999999998;
const int DIETAG=999999999;

const int N_BASES_QUADRATIC=6;
const int N_BASES_CUBIC=22;

enum {
	OUTSIDE_DOMAIN,
	FAILURE_GETTING_FIELD,
	FAILURE_GETTING_TAG,
};

#include <map>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <numeric>
#include <vector>
#include <boost/array.hpp>
#include <boost/ref.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_array.hpp>
#include <time.h>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/thread/thread.hpp> // For sleep function
#include <boost/numeric/odeint.hpp> // Not true boost library
//#include <boost/numeric/quadrature/adaptive.hpp> // Not true boost library
#include <boost/math/tools/roots.hpp>
#include <set>
//#include <type_traits> // Requires -std=c++0x compiler flag
#include <boost/type_traits.hpp>
//#include <sys/time.h>
#include <time.h>
//#include <boost/timer/timer.hpp>
#include <algorithm>
#include <boost/program_options.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>

#include "iMesh.h"
#define CHECK(a) if (iBase_SUCCESS != ierr) printf("%s\n", a), exit(ierr)

#include "Eigen/Dense"
#include "Eigen/Sparse"
// There appear to be issues with FMDB and SuperLU in the mesher
#ifndef MESHER
#ifndef HAVE_MPI
#include <Eigen/SuperLUSupport>
#endif
#endif
#include "cubature.h"
#include "cuba.h"

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkTetra.h>
#include <vtkGenericCell.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkCellType.h>
#include <vtkDataSetMapper.h>
#include <vtkCellTreeLocator.h>
#include <vtkCellLocator.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>

class CodeField;

const int NDIM=3;
const int INTERPOLATIONORDER=1;
// Matrix gives better pretty printing than Vector3d in gdb
typedef Eigen::Matrix<double,NDIM,1> vect3d;
typedef map<vect3d,int,bool(*)(vect3d,vect3d)> vect3dMap;
typedef iBase_EntityHandle entHandle;

#endif /* TYPES_H_ */
