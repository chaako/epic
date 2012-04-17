/*
 * Orbit.h
 *
 *  Created on: Feb 10, 2012
 *      Author: chaako
 */

#ifndef ORBIT_H_
#define ORBIT_H_

class ElectricField;

class Orbit {
public:
//	Orbit(Eigen::Vector3d inputPosition, Eigen::Vector3d inputVelocity,
//			iBase_EntityHandle inputNode=NULL);
	Orbit(Mesh *inputMesh, iBase_EntityHandle inputNode,
			Eigen::Vector3d inputVelocity, double inputCharge=1.);
	~Orbit();
	void integrate(ElectricField& electricField,
			PotentialField& potentialField,
			Field<int>& faceType, CodeField& vertexTypeField,
			FILE *outFile=NULL);
	Mesh *mesh_ptr;
	iBase_EntityHandle initialNode;
	iBase_EntityHandle currentElement;
	Eigen::Vector3d initialPosition;
	Eigen::Vector3d initialVelocity;
	Eigen::Vector3d finalPosition;
	Eigen::Vector3d finalVelocity;
	bool negativeEnergy;
	int vertexType;
	int finalFaceType;
	double charge;
};

#endif /* ORBIT_H_ */
