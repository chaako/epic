/*
 * SpatialDependence.cpp
 *
 *  Created on: Jun 6, 2013
 *      Author: chaako
 */

#include "SpatialDependence.h"

SpatialDependence::SpatialDependence() :
		normalization(1.),
		direction(1.,0.,0.) {
}

SpatialDependence::SpatialDependence(double value) :
		normalization(value),
		direction(1.,0.,0.)  {
}

SpatialDependence::~SpatialDependence() {
	// TODO Auto-generated destructor stub
}

double SpatialDependence::operator()(vect3d position) {
	return normalization;
}


ExponentialDependence::ExponentialDependence(double scaleLength) :
		scaleLength(scaleLength),
		referencePosition(0.,0.,0.),
		referenceValue(1.) {
}

ExponentialDependence::ExponentialDependence(double scaleLength,
		vect3d referencePosition, double referenceValue) :
		scaleLength(scaleLength),
		referencePosition(referencePosition),
		referenceValue(referenceValue) {
}

ExponentialDependence::~ExponentialDependence() {
	// TODO Auto-generated destructor stub
}

double ExponentialDependence::operator()(vect3d position) {
	vect3d relativePosition = position-referencePosition;
	return referenceValue*exp(-relativePosition.dot(direction)/scaleLength);
}


TanhDependence::TanhDependence(double scaleLength) :
		scaleLength(scaleLength),
		referencePosition(0.,0.,0.),
		referenceValue(1.) {
}

TanhDependence::TanhDependence(double scaleLength,
		vect3d referencePosition, double referenceValue,
		double normalization) :
		SpatialDependence(normalization),
		scaleLength(scaleLength),
		referencePosition(referencePosition),
		referenceValue(referenceValue) {
}

TanhDependence::~TanhDependence() {
	// TODO Auto-generated destructor stub
}

double TanhDependence::operator()(vect3d position) {
	vect3d relativePosition = position-referencePosition;
	double eps=0.;
	// TODO: deal with exp(-bigNumber)==0 some other way?
	if (referenceValue==1.)
		eps = LENGTH_TOLERANCE;
	return normalization*(referenceValue + eps +
			tanh(-relativePosition.dot(direction)/scaleLength));
}
