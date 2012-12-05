/*
 * functions.h
 *
 *  Created on: Dec 5, 2012
 *      Author: chaako
 */

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include "typesAndDefinitions.h"

int custom_importVTK(mMesh *, const char *);

void distributionFunctionFromBoundary(unsigned ndim, const double *x,
		void *integrandContainer_ptr, unsigned fdim, double *fval);
int distributionFunctionFromBoundaryCuba(const int *ndim, const double x[],
  const int *ncomp, double f[], void *integrandContainer_ptr);

int intersect_RayTriangle(vector<vect3d> R,
		vector<vect3d> T, vect3d* I);


#endif /* FUNCTIONS_H_ */
