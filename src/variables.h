/*
 * variables.h
 *
 *  Created on: Dec 5, 2012
 *      Author: chaako
 */

#ifndef VARIABLES_H_
#define VARIABLES_H_

#include "typesAndDefinitions.h"

// External timing variables
extern clock_t extern_findTet, extern_checkIfInNewTet;

// TODO: find better way to distinguish orbits in output
extern int extern_orbitNumber;

// TODO: don't use external pointer to vector
extern vector<vect3d> *extern_evalPositions_ptr;
extern bool extern_onlyDoEvalPositions;
extern bool extern_computeGlobDeriv;
extern int extern_nodeToComputeDerivAt;
extern vect3dMap *extern_vtkIdOfSurfacePoint_ptr;

extern int extern_numberOfSurfaceEvalPoints;
extern bool extern_saveOrbits;

// TODO: sort out units of E and B (electrons vs. ions etc.)
extern vect3d extern_B;
extern vect3d extern_E;
extern vect3d extern_VEXB;

#endif /* VARIABLES_H_ */
