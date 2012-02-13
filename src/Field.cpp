/*
 * Field.cpp
 *
 *  Created on: Feb 10, 2012
 *      Author: chaako
 */

#include "epic.h"
#include "Field.h"

template <class T>
Field<T>::Field(Mesh *inputMesh_ptr, std::string inputName,
		iBase_TagHandle inputTag) {
	mesh_ptr = inputMesh_ptr;
	name = inputName;
	if (inputTag) {
		tag = inputTag;
	} else {
		// no tag specified, so create
		int ierr;
		iMesh_createTag(mesh_ptr->meshInstance, name.c_str(),
				(int)sizeof(T),	iBase_BYTES, &tag, &ierr,
				(int)name.length());
		CHECK("Failure creating tag");
	}
}


template <>
double Field<double>::getField(Eigen::Vector3d position) {
	double field;

	return field;
}

template <>
double Field<double>::getField(iBase_EntityHandle node) {
	double field;

	return field;
}


template <>
int Field<int>::getField(Eigen::Vector3d position) {
	int field;

	return field;
}

template <>
int Field<int>::getField(iBase_EntityHandle node) {
	int field;

	return field;
}


template <>
Eigen::Vector3d Field<Eigen::Vector3d>::getField(Eigen::Vector3d position) {
	Eigen::Vector3d field;

	return field;
}

template <>
Eigen::Vector3d Field<Eigen::Vector3d>::getField(iBase_EntityHandle node) {
	Eigen::Vector3d field;

	return field;
}

ElectricField::ElectricField(Mesh *inputMesh_ptr, std::string inputName,
		iBase_TagHandle inputTag)
		: Field(inputMesh_ptr, inputName, inputTag) {
}

void ElectricField::calcField() {

}


PotentialField::PotentialField(Mesh *inputMesh_ptr, std::string inputName,
	iBase_TagHandle inputTag)
	: Field(inputMesh_ptr, inputName, inputTag) {
}

void PotentialField::calcField() {

}


DensityField::DensityField(Mesh *inputMesh_ptr, std::string inputName,
		iBase_TagHandle inputTag)
		: Field(inputMesh_ptr, inputName, inputTag) {
}

void DensityField::calcField() {

}

