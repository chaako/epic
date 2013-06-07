/*
 * SpatialDependence.cpp
 *
 *  Created on: Jun 6, 2013
 *      Author: chaako
 */

#include "SpatialDependence.h"

SpatialDependence::SpatialDependence() {
	// TODO Auto-generated constructor stub
	normalization = 1;
}

SpatialDependence::SpatialDependence(double value) :
		normalization(value) {
	// TODO Auto-generated constructor stub
}

SpatialDependence::~SpatialDependence() {
	// TODO Auto-generated destructor stub
}

double SpatialDependence::operator()(vect3d position) {
	return normalization;
}
