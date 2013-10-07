/*
 * functions.h
 *
 *  Created on: Dec 5, 2012
 *      Author: chaako
 */

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include "typesAndDefinitions.h"
class mMesh;

int custom_importVTK(mMesh *, const char *);

void distributionFunctionFromBoundary(unsigned ndim, const double *x,
		void *integrandContainer_ptr, unsigned fdim, double *fval);
int distributionFunctionFromBoundaryCuba(const int *ndim, const double x[],
  const int *ncomp, double f[], void *integrandContainer_ptr);

int intersect_RayTriangle(vector<vect3d> R,
		vector<vect3d> T, vect3d* I);
int intersect_SegmentTriangle(vector<vect3d> R,
		vector<vect3d> T, vect3d* I );

bool vect3dLessThan(vect3d a, vect3d b);

void getIterationDataFromFile(boost::filesystem::path& inputPath,
		string& inputMeshFile,
		Eigen::VectorXd *diisPotential, Eigen::VectorXd *residual);

#ifdef HAVE_MPI
void requestIterationDataFromSlaves(int numberOfFiles, int numberOfNodes,
		vector<Eigen::VectorXd> *diisPotentials, vector<Eigen::VectorXd> *diisResiduals);
MPI::Status receiveIterationData(int numberOfNodes,
		vector<Eigen::VectorXd> *diisPotentials, vector<Eigen::VectorXd> *diisResiduals);
void processIterationDataRequests(int numberOfNodes, vector<boost::filesystem::path>& inputPaths,
		vector<string>& inputMeshFiles);
#endif

#endif /* FUNCTIONS_H_ */
