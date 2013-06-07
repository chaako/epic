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

	double operator()(vect3d position);

	double normalization;
};

#endif /* SPATIALDEPENDENCE_H_ */
