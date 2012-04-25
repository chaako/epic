/*
 * Field.h
 *
 *  Created on: Feb 10, 2012
 *      Author: chaako
 */

#ifndef FIELD_H_
#define FIELD_H_

#include <stdio.h>
//#include <type_traits> // Requires -std=c++0x compiler flag

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "iMesh.h"
#include "Eigen/Dense"
#include "Mesh.h"

#define CHECK(a) if (iBase_SUCCESS != ierr) printf("%s\n", a), exit(ierr)

template <class T>
class Field {
public:
	Field(Mesh *inputMesh_ptr, std::string inputName,
			int entityDimension);
	virtual ~Field() {}

	T getField(Eigen::Vector3d position);
	T getField(iBase_EntityHandle node);
	T getAverageField(iBase_EntityHandle element);
	void setField(iBase_EntityHandle node, T field);

	Mesh *mesh_ptr;
	std::string name;
	iBase_TagHandle tag;
	std::vector<iBase_EntityHandle> entities;
};

class ElectricField : public Field<Eigen::Vector3d> {
public:
	ElectricField(Mesh *inputMesh_ptr, std::string inputName);
	virtual ~ElectricField() {}

	void calcField(PotentialField potentialField);
};

class CodeField : public Field<int> {
public:
	CodeField(Mesh *inputMesh_ptr, std::string inputName,
			int elementType);
	virtual ~CodeField() {}

	void calcField(Field<int> faceTypeField);
};

class DensityField : public Field<double> {
	// Charge density
public:
	DensityField(Mesh *inputMesh_ptr, std::string inputName);
	virtual ~DensityField() {}

	void calcField();
	void calcField(DensityField ionDensity, DensityField electronDensity);
	void calcField(ElectricField electricField,
			PotentialField potentialField,
			Field<int> faceType, CodeField vertexType, double charge=1.,
			FILE *outFile=NULL);
	double calculateDensity(int node, ElectricField electricField,
			PotentialField potentialField,
			Field<int> faceType, CodeField vertexType, double charge,
			double *error);
#ifdef HAVE_MPI
	void requestDensityFromSlaves(ElectricField electricField,
			PotentialField potentialField,
			Field<int> faceType, CodeField vertexType, double charge,
			FILE *outFile);
	MPI::Status receiveDensity(FILE *outFile);
	void processDensityRequests(ElectricField electricField,
			PotentialField potentialField,
			Field<int> faceType, CodeField vertexType, double charge);
#endif

};

class PotentialField : public Field<double> {
public:
	PotentialField(Mesh *inputMesh_ptr, std::string inputName);
	virtual ~PotentialField() {}

	void calcField();
	void calcField(DensityField ionDensity, DensityField electronDensity,
			FILE *outFile);

};

// gcc doesn't implement the export keyword, so define template functions here
template <class T>
Field<T>::Field(Mesh *inputMesh_ptr, std::string inputName,
		int entityDimension) {
	mesh_ptr = inputMesh_ptr;
	name = inputName;
	int ierr;
	int size, type;
//	if (std::is_same<T,double>::value) {
//		size = 1;
//		type = iBase_DOUBLE;
//	} else if (std::is_same<T,int>::value) {
//		size = 1;
//		type = iBase_INTEGER;
//	} else {
		size = (int)sizeof(T);
		type = iBase_BYTES;
//	}
	tag = mesh_ptr->createTagHandle(name, size, type);
	// TODO: consider more transparent handling of exiting tag here
	tag = mesh_ptr->getTagHandle(name);
	entities = mesh_ptr->getEntities(entityDimension);
}

template <class T>
T Field<T>::getField(Eigen::Vector3d position) {
	// TODO: write this function
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
T Field<T>::getAverageField(iBase_EntityHandle element) {
	// TODO: this doesn't work for non-vertex fields
	std::vector<iBase_EntityHandle> vertices =
			mesh_ptr->getAdjacentEntities(element, iBase_VERTEX);
	T average=0;

	for (int i=0; i<vertices.size(); i++) {
		average += this->getField(vertices[i]);
	}
	average /= (double)vertices.size();

	return average;
}

template <class T>
void Field<T>::setField(iBase_EntityHandle node, T field) {
	int ierr;
	int field_size = sizeof(T);
	iMesh_setData(mesh_ptr->meshInstance, node, tag, &field,
			field_size, &ierr);
}

#endif /* FIELD_H_ */


