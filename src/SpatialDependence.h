/*
 * SpatialDependence.h
 *
 *  Created on: Jun 6, 2013
 *      Author: chaako
 */

#ifndef SPATIALDEPENDENCE_H_
#define SPATIALDEPENDENCE_H_

#include "typesAndDefinitions.h"

class SpatialDependence {
public:
	SpatialDependence();
	SpatialDependence(double value);
	virtual ~SpatialDependence();

	virtual double operator()(vect3d position);

	double normalization;
	vect3d direction;
};

class ExponentialDependence : public SpatialDependence {
public:
	ExponentialDependence(double scaleLength);
	ExponentialDependence(double scaleLength, vect3d referencePosition,
			double referenceValue);
	virtual ~ExponentialDependence();

	double operator()(vect3d position);

	double scaleLength; // Allow to be negative
	vect3d referencePosition;
	double referenceValue;
};

class TanhDependence : public SpatialDependence {
public:
	TanhDependence(double scaleLength);
	TanhDependence(double scaleLength, vect3d referencePosition,
			double referenceValue, double normalization);
	virtual ~TanhDependence();

	double operator()(vect3d position);

	double scaleLength; // Allow to be negative
	vect3d referencePosition;
	double referenceValue;
};

#endif /* SPATIALDEPENDENCE_H_ */
