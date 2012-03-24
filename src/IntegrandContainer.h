/*
 * IntegrandContainer.h
 *
 *  Created on: Feb 21, 2012
 *      Author: chaako
 */

#ifndef INTEGRANDCONTAINER_H_
#define INTEGRANDCONTAINER_H_

class Mesh;
template <class T> class Field;
class CodeField;
class ElectricField;
class PotentialField;

class IntegrandContainer {
public:
	IntegrandContainer();
	virtual ~IntegrandContainer();

	Mesh *mesh_ptr;
	iBase_EntityHandle node;
	ElectricField *electricField_ptr;
	PotentialField *potentialField_ptr;
	Field<int> *faceTypeField_ptr;
	CodeField *vertexTypeField_ptr;
	FILE *outFile;
	double charge;
};

#endif /* INTEGRANDCONTAINER_H_ */
