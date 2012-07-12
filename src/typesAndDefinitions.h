/*
 * types.h
 *
 *  Created on: Jun 18, 2012
 *      Author: chaako
 */

#ifndef TYPES_H_
#define TYPES_H_

using namespace std;

const double VOLUME_TOLERANCE=1.e-8;
const double LENGTH_TOLERANCE=1.e-10;
const double DELTA_LENGTH=1.e-10;
//const double DELTA_LENGTH=1.e-5;

const int WORKTAG=999999998;
const int DIETAG=999999999;

const int N_BASES_QUADRATIC=6;
const int N_BASES_CUBIC=22;

enum {
	OUTSIDE_DOMAIN,
	FAILURE_GETTING_FIELD,
};

#include <map>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <numeric>
#include <vector>
#include <boost/array.hpp>
#include <boost/ref.hpp>
#include <boost/numeric/odeint.hpp> // Not true boost library
#include <set>
//#include <type_traits> // Requires -std=c++0x compiler flag
#include <boost/type_traits.hpp>
//#include <sys/time.h>
#include <time.h>
//#include <boost/timer/timer.hpp>
#include <algorithm>

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "iMesh.h"
#define CHECK(a) if (iBase_SUCCESS != ierr) printf("%s\n", a), exit(ierr)

#include "Eigen/Dense"
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

const int NDIM=3;
const int INTERPOLATIONORDER=1;
// Matrix gives better pretty printing than Vector3d in gdb
typedef Eigen::Matrix<double,NDIM,1> vect3d;
typedef iBase_EntityHandle entHandle;


#endif /* TYPES_H_ */
