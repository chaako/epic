/*
 * Field.h
 *
 *  Created on: Feb 10, 2012
 *      Author: chaako
 */

#ifndef FIELD_H_
#define FIELD_H_

#include "iMesh.h"
#include "Eigen/Dense"
#include "Mesh.h"

template <class T>
class Field {
public:
	Field(Mesh *inputMesh_ptr, std::string inputName,
			iBase_TagHandle inputTag=0);
	virtual ~Field() {}

	T getField(Eigen::Vector3d position);
	T getField(iBase_EntityHandle node);

	virtual void calcField() {}

	Mesh *mesh_ptr;
	std::string name;
	iBase_TagHandle tag;

};

class ElectricField : public Field<Eigen::Vector3d> {
public:
	ElectricField(Mesh *inputMesh_ptr, std::string inputName,
			iBase_TagHandle inputTag=0);
	virtual ~ElectricField() {}

	virtual void calcField();

};

class DensityField : public Field<double> {
public:
	DensityField(Mesh *inputMesh_ptr, std::string inputName,
			iBase_TagHandle inputTag=0);
	virtual ~DensityField() {}

	virtual void calcField();

};

class PotentialField : public Field<double> {
public:
	PotentialField(Mesh *inputMesh_ptr, std::string inputName,
			iBase_TagHandle inputTag=0);
	virtual ~PotentialField() {}

	virtual void calcField();

};

#endif /* FIELD_H_ */


