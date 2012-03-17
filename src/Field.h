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

#define CHECK(a) if (iBase_SUCCESS != ierr) printf("%s\n", a), exit(ierr)

template <class T>
class Field {
public:
	Field(Mesh *inputMesh_ptr, std::string inputName,
			iBase_TagHandle inputTag=0);
	virtual ~Field() {}

	T getField(Eigen::Vector3d position);
	T getField(iBase_EntityHandle node);

	void setField(iBase_EntityHandle node, T field);

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
	// Charge density
public:
	DensityField(Mesh *inputMesh_ptr, std::string inputName,
			iBase_TagHandle inputTag=0);
	virtual ~DensityField() {}

	virtual void calcField();
	void calcField(DensityField ionDensity, DensityField electronDensity);
	void calcField(ElectricField electricField,
			PotentialField potentialField, double charge=1.);

};

class PotentialField : public Field<double> {
public:
	PotentialField(Mesh *inputMesh_ptr, std::string inputName,
			iBase_TagHandle inputTag=0);
	virtual ~PotentialField() {}

	virtual void calcField();
	void calcField(DensityField ionDensity, DensityField electronDensity);

};


// It would seem template function definitions need to be in the same file
template <class T>
T Field<T>::getField(Eigen::Vector3d position) {
	T field;
	return field;
}

template <class T>
T Field<T>::getField(iBase_EntityHandle node) {
	T field;
	T *field_ptr=&field;
	int field_alloc = sizeof(T);
	int field_size = sizeof(T);
	int ierr;
	iMesh_getData(mesh_ptr->meshInstance, node, tag, &field_ptr,
			&field_alloc, &field_size, &ierr);
	CHECK("Failure getting field");
	return field;
}

template <class T>
void Field<T>::setField(iBase_EntityHandle node, T field) {
	int ierr;
	int field_size = sizeof(T);
	iMesh_setData(mesh_ptr->meshInstance, node, tag, &field,
			field_size, &ierr);
}

#endif /* FIELD_H_ */


