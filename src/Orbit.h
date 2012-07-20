/*
 * Orbit.h
 *
 *  Created on: Feb 10, 2012
 *      Author: chaako
 */

#ifndef ORBIT_H_
#define ORBIT_H_

#include "typesAndDefinitions.h"

class ElectricField;

class Orbit {
public:
//	Orbit(vect3d inputPosition, vect3d inputVelocity,
//			entHandle inputNode=NULL);
	Orbit(Mesh *inputMesh, entHandle inputNode,
			vect3d inputVelocity, double inputCharge=1.);
	~Orbit();
	void integrate(ElectricField& electricField,
			PotentialField& potentialField,
			Field<int>& faceType, CodeField& vertexTypeField,
			FILE *outFile=NULL);
	void integrate(PotentialField& potentialField,
			ElectricField& electricField,
			Field<int>& faceType, CodeField& vertexTypeField,
			ShortestEdgeField shortestEdgeField, FILE *outFile=NULL);
	Mesh *mesh_ptr;
	entHandle initialNode;
	entHandle currentElement;
	vect3d initialPosition;
	vect3d initialVelocity;
	vect3d finalPosition;
	vect3d finalVelocity;
	bool negativeEnergy;
	int vertexType;
	int finalFaceType;
	double finalPotential;
	double charge;
};

#endif /* ORBIT_H_ */
