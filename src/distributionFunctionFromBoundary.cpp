#include "epic.h"

void distributionFunctionFromBoundary(unsigned ndim, const double *x,
		void *integrandContainer_ptr, unsigned fdim, double *fval) {
	if (ndim!=3)
		throw string("only handle ndim=3 in distributionFunctionFromBoundary");
//	assert(fdim==1);
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
	// TODO: find better way to distinguish orbits in output
	extern_orbitNumber++;
//	orbit.integrate(*potentialField_ptr, *electricField_ptr,
//			*faceTypeField_ptr, *vertexTypeField_ptr,
//			*shortestEdgeField_ptr, orbitOutFile, extern_orbitNumber);
	try {
		orbit.integrate(*potentialField_ptr, *electricField_ptr,
				*faceTypeField_ptr, *vertexTypeField_ptr,
				*shortestEdgeField_ptr, NULL, extern_orbitNumber);
	} catch (string& message) {
		cout << "Caught in distributionFunctionFromBoundary():" << message << endl;
	} catch (int code) {
		cout << "Caught in distributionFunctionFromBoundary(): " << code << endl;
	}
//	// TODO: more transparent handling of external ExB drift?
//	vect3d finalVelocity = orbit.finalVelocity - VEXB;
	vect3d finalVelocity = orbit.finalVelocity;
	fval[0] = 0.;
	// TODO: shouldn't hard-code domain here
	if ( orbit.finalFaceType==5
			&& !orbit.negativeEnergy) {
		fval[0] = 1./pow(2.*M_PI,3./2.);
		fval[0] *= exp(-(pow(finalVelocity.norm(),2.)-pow(v,2.))/2.);
		fval[0] *= M_PI;
		fval[0] *= (1.+x[0])*sqrt(-log((1.+x[0])/2.));
//		// TODO: if fix outer potential this may not give equal ion and electron dens.
//		// TODO: this doesn't work when including external E in potential
//		if (charge*orbit.finalPotential>0.)
//			fval[0] *= exp(-charge*orbit.finalPotential);
	} else if (orbit.finalPosition.norm()>1.) {
//		std::cout << "orbit terminated with final position in domain:" <<
//				std::endl << orbit.finalPosition.transpose() << " r=" <<
//				orbit.finalPosition.norm() << " " << orbit.finalFaceType <<
//				" " << orbit.negativeEnergy << std::endl;
//		std::cout << ".";
	}
	// TODO: don't hard-code moment order?
	// -VEXB since orbit.initialVelocity=-velocity
	vect3d driftingVelocity = velocity - VEXB;
	// Since density isn't available yet, divide by it later
	if (fdim>=4) {
		for (int i=0; i<NDIM; i++)
			fval[i+1] = fval[0]*driftingVelocity[i];
	}
	// Since average velocity isn't available yet, subtract off KE from T later
	if (fdim>=5)
		fval[4] = fval[0]*pow(driftingVelocity.norm(),2.)/3.;
	// TODO: better way than integrating orbit again for output?
	if (orbitOutFile) {
		Orbit orbitForOutput(mesh_ptr,node,velocity,charge);
		try {
			orbitForOutput.integrate(*potentialField_ptr, *electricField_ptr,
					*faceTypeField_ptr, *vertexTypeField_ptr,
					*shortestEdgeField_ptr, orbitOutFile, fval[0]);
//					*shortestEdgeField_ptr, orbitOutFile, exp(-pow(finalVelocity.norm(),2.)/2.));
		} catch (string& message) {
			cout << "Caught in distributionFunctionFromBoundary():" << message << endl;
		} catch (int code) {
			cout << "Caught in distributionFunctionFromBoundary(): " << code << endl;
		}
	}
	if (((IntegrandContainer*)integrandContainer_ptr)->outFile) {
		fprintf(((IntegrandContainer*)integrandContainer_ptr)->outFile,
				"%f %f %f %f\n", x[0], x[1], x[2], fval[0]);
	}
}


int distributionFunctionFromBoundaryCuba(const int *ndim, const double x[],
  const int *ncomp, double f[], void *integrandContainer_ptr) {
	// TODO: Can get weight and iteration number from Vegas as two additional
	//       optional arguments, which would allow storing dist. func. for
	//       computation of other moments after density integral is complete.
	if (*ndim!=3)
		throw string("only handle ndim=3 in distributionFunctionFromBoundaryCuba");
//	assert(*ncomp==1);
	double y[3];
	for (int i=0; i<*ndim; i++) {
		y[i] = 2.*x[i]-1;
	}
	distributionFunctionFromBoundary((unsigned)*ndim, y,
			integrandContainer_ptr, *ncomp, f);
	for (int i=0; i<*ncomp; i++)
		f[i] *= pow(2.,(double)*ndim);

	return 0;
}
