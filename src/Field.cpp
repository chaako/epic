/*
 * Field.cpp
 *
 *  Created on: Feb 10, 2012
 *      Author: chaako
 */

#include "epic.h"
//#include "Field.h"

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

//double ScalarField::getField(iBase_EntityHandle node) {
//	double field;
//
//	return field;
//}

ScalarField::ScalarField(Mesh *inputMesh_ptr, std::string inputName,
		iBase_TagHandle inputTag)
		: Field(inputMesh_ptr, inputName, inputTag) {
}

VectorField::VectorField(Mesh *inputMesh_ptr, std::string inputName,
		iBase_TagHandle inputTag)
		: Field(inputMesh_ptr, inputName, inputTag) {
}

ElectricField::ElectricField(Mesh *inputMesh_ptr, std::string inputName,
		iBase_TagHandle inputTag)
		: VectorField(inputMesh_ptr, inputName, inputTag) {
}

PotentialField::PotentialField(Mesh *inputMesh_ptr, std::string inputName,
	iBase_TagHandle inputTag)
	: ScalarField(inputMesh_ptr, inputName, inputTag) {
}

DensityField::DensityField(Mesh *inputMesh_ptr, std::string inputName,
		iBase_TagHandle inputTag)
		: ScalarField(inputMesh_ptr, inputName, inputTag) {
}


//Eigen::Vector3d VectorField::getField(
//		Eigen::Vector3d postion) {
//	Eigen::Vector3d field;
//
//	return field;
//};

