/*
 * DistributionFunction.h
 *
 *  Created on: Jun 6, 2013
 *      Author: chaako
 */

#ifndef DISTRIBUTIONFUNCTION_H_
#define DISTRIBUTIONFUNCTION_H_

#include "typesAndDefinitions.h"
#include "SpatialDependence.h"

class DistributionFunction {
public:
	DistributionFunction();
	virtual ~DistributionFunction();

	virtual double operator()(vect3d position, vect3d velocity) const;
	virtual double operator()(double xi) const;

	void setPresetPosition(vect3d position);

	vect3d presetPosition;
};

class Maxwellian : public DistributionFunction {
public:
	Maxwellian(SpatialDependence& density, vect3d magneticAxis,
			SpatialDependence& parallelTemperature,
			SpatialDependence& perpendicularTemperature,
			SpatialDependence& parallelDrift,
			vect3d perpendicularDrift);
	virtual ~Maxwellian();

	double operator()(vect3d position, vect3d velocity) const;
	double operator()(double xi) const;

	SpatialDependence *density_ptr;
	vect3d magneticAxis;
	SpatialDependence *parallelTemperature_ptr;
	SpatialDependence *perpendicularTemperature_ptr;
	SpatialDependence *parallelDrift_ptr;
	vect3d perpendicularDrift;
};

//// TODO: Allow for gradients not aligned with x
//struct DistributionFunctionVy {
//	double operator()(double vy) const {
//		// TODO: if include drift in velocity, should subtract
//		return distFunc_ptr->operator()(pos, vect3d(0.,vy,0.));
//	}
//	DistributionFunction *distFunc_ptr;
//	vect3d pos;
//};

#endif /* DISTRIBUTIONFUNCTION_H_ */
