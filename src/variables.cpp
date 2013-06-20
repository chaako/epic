#include "variables.h"

clock_t extern_findTet=0, extern_checkIfInNewTet=0;

// TODO: find better way to distinguish orbits in output
int extern_orbitNumber = 0;

// TODO: don't use external pointer to vector
vector<vect3d> *extern_evalPositions_ptr=NULL;
map<vect3d,vtkIdType,bool(*)(vect3d,vect3d)> *extern_vtkIdOfSurfacePoint_ptr=NULL;

int extern_numberOfSurfaceEvalPoints = 0;
bool extern_saveOrbits = false;

// TODO: sort out units of E and B (electrons vs. ions etc.)
vect3d extern_B(0.,0.,0.);
vect3d extern_E(0.,0.,0.);
vect3d extern_VEXB(0.,0.,0.);
