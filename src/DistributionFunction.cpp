/*
 * DistributionFunction.cpp
 *
 *  Created on: Jun 6, 2013
 *      Author: chaako
 */

#include "DistributionFunction.h"

DistributionFunction::DistributionFunction() {
	// TODO Auto-generated constructor stub

}

DistributionFunction::~DistributionFunction() {
	// TODO Auto-generated destructor stub
}

double DistributionFunction::operator()(vect3d position, vect3d velocity) const {
	double distFunc = exp(-pow(velocity.norm(),2.)/2.)/pow(2.*M_PI,3./2.);
//	// TODO: debugging
//	cout << distFunc << " " << velocity.transpose() << endl;
	return distFunc;
}

// TODO: Allow for gradients not aligned with x?
double DistributionFunction::operator()(double xi) const {
	// Integrated over x and z, and rescaled to interval (-1.,1.)
	double vy = xi/(1-xi*xi);
	double weightedDistFunc = (1+xi*xi)/pow(1-xi*xi,2.)*2.*M_PI*
			this->operator()(presetPosition, vect3d(0.,vy,0.));
//	// TODO: for debugging:
//	cout << weightedDistFunc << " " << xi << " " << vy << " " <<
//			presetPosition.transpose() << " " <<
//			this->operator()(presetPosition, vect3d(0.,vy,0.)) << endl;
	return weightedDistFunc;
}

void DistributionFunction::setPresetPosition(vect3d position) {
	presetPosition = position;
}

Maxwellian::Maxwellian(SpatialDependence& density, vect3d magneticAxis,
		SpatialDependence& parallelTemperature,
		SpatialDependence& perpendicularTemperature,
		SpatialDependence& parallelDrift,
		vect3d perpendicularDrift) :
				density_ptr(&density),
				magneticAxis(magneticAxis),
				parallelTemperature_ptr(&parallelTemperature),
				perpendicularTemperature_ptr(&perpendicularTemperature),
				parallelDrift_ptr(&parallelDrift),
				perpendicularDrift(perpendicularDrift) {
}

Maxwellian::~Maxwellian() {
}

double Maxwellian::operator()(vect3d position, vect3d velocity) const {
	double density = density_ptr->operator()(position);
	double parallelTemperature = parallelTemperature_ptr->operator()(position);
	double perpendicularTemperature = perpendicularTemperature_ptr->operator()(position);
	vect3d parallelDrift = parallelDrift_ptr->operator()(position)*magneticAxis;
	vect3d drift = perpendicularDrift + parallelDrift;
	vect3d relativeVelocity = velocity-drift;
	vect3d parallelRelativeVelocity = relativeVelocity.dot(magneticAxis)*magneticAxis;
	vect3d perpendicularRelativeVelocity = relativeVelocity-parallelRelativeVelocity;
	double distributionFunction =
			density*exp(-pow(parallelRelativeVelocity.norm(),2.)/2./parallelTemperature)*
			exp(-pow(perpendicularRelativeVelocity.norm(),2.)/2./perpendicularTemperature)/
			pow(2.*M_PI,3./2.)/sqrt(parallelTemperature)/perpendicularTemperature;
//	// TODO: debugging
//	cout << distributionFunction << " " << position.transpose() << " " <<
//			velocity.transpose() << endl;
	return distributionFunction;
}

// TODO: Allow for gradients not aligned with x?
double Maxwellian::operator()(double xi) const {
	// Integrated over x and z, and rescaled to interval (-1.,1.)
	double vy = xi/(1-xi*xi);
	// TODO: don't restrict gradients to be in x-dir?
	//       (for now only vy gives displacement in gradient dir)
	vect3d relevantPerpendicularVelocity(0.,vy,0.);
	// TODO: check sign of Larmor vector subtraction
	vect3d guidingCenter = presetPosition +
			magneticAxis.cross(relevantPerpendicularVelocity)/extern_B.norm();
	double parallelTemperature = parallelTemperature_ptr->operator()(guidingCenter);
	double perpendicularTemperature = perpendicularTemperature_ptr->operator()(guidingCenter);
	vect3d parallelDrift = parallelDrift_ptr->operator()(guidingCenter)*magneticAxis;
	vect3d drift = perpendicularDrift + parallelDrift;
	return (1+xi*xi)/pow(1-xi*xi,2.)*2.*M_PI*sqrt(parallelTemperature*perpendicularTemperature)*
			this->operator()(guidingCenter, relevantPerpendicularVelocity+drift);
}

