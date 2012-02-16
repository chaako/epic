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
Field<double>::Field(Mesh *inputMesh_ptr, std::string inputName,
		iBase_TagHandle inputTag) {
	mesh_ptr = inputMesh_ptr;
	name = inputName;
	if (inputTag) {
		tag = inputTag;
	} else {
		// no tag specified, so create
		int ierr;
		iMesh_createTag(mesh_ptr->meshInstance, name.c_str(),
				1,	iBase_DOUBLE, &tag, &ierr,
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
	int ierr;
	iBase_EntityHandle *ents0d = NULL;
	int ents0d_alloc = 0, ents0d_size;
	iBase_TagHandle potential_tag;
	std::string potentialName = "potential";
	iMesh_getTagHandle(mesh_ptr->meshInstance, potentialName.c_str(),
			&potential_tag, &ierr, potentialName.length());
	CHECK("Failed to get potential tag");
	// Calculate electric field for all 0d elements
	iMesh_getEntities(mesh_ptr->meshInstance, mesh_ptr->rootEntitySet,
			iBase_VERTEX, iMesh_ALL_TOPOLOGIES,
			&ents0d, &ents0d_alloc, &ents0d_size, &ierr);
	CHECK("Couldn't get vertex entities");
	for (int i = 0; i < ents0d_size; i++) {
		std::vector<iBase_EntityHandle> superCellFaces =
				getSuperCellFaces(mesh_ptr->meshInstance, ents0d[i]);
		Eigen::Vector3d eField(0.,0.,0.);

		double x,y,z;
		iMesh_getVtxCoord(mesh_ptr->meshInstance, ents0d[i], &x, &y, &z,
				&ierr);
		CHECK("Failure getting vertex coordinates");
		Eigen::Vector3d point(x,y,z);
		double volume=0;

		for (int j=0; j<superCellFaces.size(); j++)  {
			Eigen::Vector3d surfaceVector =
					getSurfaceVector(mesh_ptr->meshInstance, point,
							superCellFaces[j]);
			double potential = getAverageDblData(mesh_ptr->meshInstance,
					superCellFaces[j], potential_tag);
			volume += getTetVolume(mesh_ptr->meshInstance, point, superCellFaces[j]);
			eField -= potential*surfaceVector;
		}
		eField /= volume;

		iMesh_setData(mesh_ptr->meshInstance, ents0d[i], tag, &eField,
				sizeof(Eigen::Vector3d), &ierr);
		CHECK("Failure setting eField tag");
	}
}


PotentialField::PotentialField(Mesh *inputMesh_ptr, std::string inputName,
	iBase_TagHandle inputTag)
	: Field(inputMesh_ptr, inputName, inputTag) {
}

void PotentialField::calcField() {
	int ierr;
	iBase_EntityHandle *ents0d = NULL;
	int ents0d_alloc = 0, ents0d_size;
	// set potential value for all 0d elements
	iMesh_getEntities(mesh_ptr->meshInstance, mesh_ptr->rootEntitySet,
			iBase_VERTEX, iMesh_ALL_TOPOLOGIES,
			&ents0d, &ents0d_alloc, &ents0d_size, &ierr);
	CHECK("Couldn't get vertex entities");
	for (int i = 0; i < ents0d_size; i++) {
		double potential, x, y, z;
		iMesh_getVtxCoord(mesh_ptr->meshInstance, ents0d[i], &x, &y, &z,
				&ierr);
		CHECK("Failure getting vertex coordinates");
		potential = -1./sqrt(x*x+y*y+z*z);
		iMesh_setDblData(mesh_ptr->meshInstance, ents0d[i], tag, potential,
				&ierr);
		CHECK("Failure setting potential tag");
	}
}


DensityField::DensityField(Mesh *inputMesh_ptr, std::string inputName,
		iBase_TagHandle inputTag)
		: Field(inputMesh_ptr, inputName, inputTag) {
}

void DensityField::calcField() {

}

