/*
 * epic.h
 *
 *  Created on: Jan 20, 2012
 *      Author: chaako
 */

#ifndef EPIC_H_
#define EPIC_H_

#include "typesAndDefinitions.h"

#include "Mesh.h"
#include "IntegrandContainer.h"
#include "Field.h"
#include "Stepper.h"
#include "VelocityAndAcceleration.h"
#include "Orbit.h"

class mMesh;
int custom_importVTK(mMesh *, const char *);

void distributionFunctionFromBoundary(unsigned ndim, const double *x,
		void *integrandContainer_ptr, unsigned fdim, double *fval);

int intersect_RayTriangle(std::vector<vect3d> R,
		std::vector<vect3d> T, vect3d* I);


// External timing variables
extern clock_t extern_findTet, extern_checkIfInNewTet;


#endif /* EPIC_H_ */
