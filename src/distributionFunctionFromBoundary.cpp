#include "epic.h"

void distributionFunctionFromBoundary(unsigned ndim, const double *x,
		void *integrandContainer_ptr, unsigned fdim, double *fval) {
	assert(ndim==3);
	assert(fdim==1);
	Mesh *mesh_ptr = ((IntegrandContainer*)integrandContainer_ptr)->mesh_ptr;
	entHandle node =
			((IntegrandContainer*)integrandContainer_ptr)->node;
	ElectricField *electricField_ptr =
			((IntegrandContainer*)integrandContainer_ptr)->electricField_ptr;
	PotentialField *potentialField_ptr =
			((IntegrandContainer*)integrandContainer_ptr)->potentialField_ptr;
	Field<int> *faceTypeField_ptr =
			((IntegrandContainer*)integrandContainer_ptr)->faceTypeField_ptr;
	CodeField *vertexTypeField_ptr =
			((IntegrandContainer*)integrandContainer_ptr)->vertexTypeField_ptr;
	ShortestEdgeField *shortestEdgeField_ptr =
			((IntegrandContainer*)integrandContainer_ptr)->shortestEdgeField_ptr;
	double charge = ((IntegrandContainer*)integrandContainer_ptr)->charge;
	FILE *orbitOutFile = ((IntegrandContainer*)integrandContainer_ptr)->orbitOutFile;
	vect3d velocity;
	double v = sqrt(-4.*log((1.+x[0])/2.));
	double theta = acos(x[1]);
	double phi = (1.+x[2])*M_PI;
	velocity[0] = v*cos(phi)*sin(theta);
	velocity[1] = v*sin(phi)*sin(theta);
	velocity[2] = v*cos(theta);
//	// TODO: shouldn't hard-code taking advantage of symmetry, but rather
//	//       get the preferred z-axis (now rHat) from the mesh (i.e. specified elsewhere)
//	vect3d position = mesh_ptr->getCoordinates(node);
//	vect3d rHat = position/position.norm();
////	std::cout << (velocity-rHat).norm() << std::endl;
////	std::cout << "rHat = " << rHat[0] << ", " << rHat[1] << ", " << rHat[2] << std::endl;
//	double beta = acos(rHat[2]);
//	double alpha = atan2(rHat[0],-rHat[1]);
////	std::cout << "vel = " << velocity[0] << ", " << velocity[1] <<
////			", " << velocity[2] << std::endl;
//	velocity = Eigen::AngleAxisd(-alpha, vect3d::UnitZ()) *
//			velocity;
//	velocity = Eigen::AngleAxisd(-beta, vect3d::UnitX()) *
//			velocity;
////	std::cout << (velocity-vect3d::UnitZ()).norm() << std::endl << std::endl;
////	std::cout << "rotVel = " << velocity[0] << ", " << velocity[1] <<
////			", " << velocity[2] << std::endl;

	Orbit orbit(mesh_ptr,node,velocity,charge);
//	orbit.integrate(*electricField_ptr, *potentialField_ptr,
//			*faceTypeField_ptr, *vertexTypeField_ptr, orbitOutFile);
	orbit.integrate(*potentialField_ptr, *electricField_ptr,
			*faceTypeField_ptr, *vertexTypeField_ptr,
			*shortestEdgeField_ptr, orbitOutFile);
	*fval = 0.;
	// TODO: shouldn't hard-code domain here
	if ( orbit.finalFaceType==5
			&& !orbit.negativeEnergy) {
		*fval = 1./pow(2.*M_PI,3./2.);
		*fval *= exp(-pow(orbit.finalVelocity.norm(),2.)/2.);
		*fval /= exp(-pow(v,2.)/2.);
		*fval *= M_PI;
		*fval *= (1.+x[0])*sqrt(-log((1.+x[0])/2.));
		// TODO: if fix outer potential this may not give equal ion and electron dens.
		if (charge*orbit.finalPotential>0.)
			*fval *= exp(-charge*orbit.finalPotential);
	} else if (orbit.finalPosition.norm()>1.) {
//		std::cout << "orbit terminated with final position in domain:" <<
//				std::endl << orbit.finalPosition.transpose() << " r=" <<
//				orbit.finalPosition.norm() << " " << orbit.finalFaceType <<
//				" " << orbit.negativeEnergy << std::endl;
//		std::cout << ".";
	}
	if (((IntegrandContainer*)integrandContainer_ptr)->outFile) {
		fprintf(((IntegrandContainer*)integrandContainer_ptr)->outFile,
				"%f %f %f %f\n", x[0], x[1], x[2], *fval);
	}
}


int distributionFunctionFromBoundaryCuba(const int *ndim, const double x[],
  const int *ncomp, double f[], void *integrandContainer_ptr) {
	assert(*ndim==3);
	assert(*ncomp==1);
	double y[3];
	for (int i=0; i<*ndim; i++) {
		y[i] = 2.*x[i]-1;
	}
	distributionFunctionFromBoundary((unsigned)*ndim, y,
			integrandContainer_ptr, *ncomp, f);
	*f *= pow(2.,(double)*ndim);

	return 0;
}
