/*
 * epic.h
 *
 *  Created on: Jan 20, 2012
 *      Author: chaako
 */

#ifndef EPIC_H_
#define EPIC_H_

#define VOLUME_TOLERANCE 1.e-8
#define LENGTH_TOLERANCE 1.e-10

#include <map>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "iMesh.h"
#include <math.h>

#include "Eigen/Dense"
#include "cubature.h"
//#include "cuba.h"

#include "mpi.h"

#include <numeric>
#include <vector>
//#include <type_traits> // Requires -std=c++0x compiler flag
//#include <sys/time.h>
#include <time.h>
//#include <boost/timer/timer.hpp>
#include <algorithm>

#include "Mesh.h"
#include "IntegrandContainer.h"
#include "Field.h"
#include "Orbit.h"

#define CHECK(a) if (iBase_SUCCESS != ierr) printf("%s\n", a), exit(ierr)

class mMesh;
int custom_importVTK(mMesh *, const char *);

void distributionFunctionFromBoundary(unsigned ndim, const double *x,
		void *integrandContainer_ptr, unsigned fdim, double *fval);

int intersect_RayTriangle(std::vector<Eigen::Vector3d> R,
		std::vector<Eigen::Vector3d> T, Eigen::Vector3d* I);


// External timing variables
extern clock_t extern_findTet, extern_checkIfInNewTet;


#endif /* EPIC_H_ */
