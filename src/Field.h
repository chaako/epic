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
#include <boost/type_traits.hpp>

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
	// TODO: Does an empty destructor cause a memory leak?
	virtual ~Field() {}

	T getField(Eigen::Vector3d position, iBase_EntityHandle *entity=NULL);
	T getField(iBase_EntityHandle node);
	T getAverageField(iBase_EntityHandle element);
	void setField(iBase_EntityHandle node, T field);

	Mesh *mesh_ptr;
	std::string name;
	iBase_TagHandle tag;
	std::vector<iBase_EntityHandle> entities;
	iBase_EntityHandle currentElement;
	std::vector<iBase_EntityHandle> currentVertices;
	std::vector<T> currentFields;
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
	PotentialField(PotentialField potential, std::string inputName);
	virtual ~PotentialField() {}

	void calcField();
	void calcField(DensityField ionDensity, DensityField electronDensity,
			CodeField vertexType, FILE *outFile);

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
	if (boost::is_same<T,double>::value) {
		size = 1;
		type = iBase_DOUBLE;
//	} else if (std::is_same<T,int>::value) {
	} else if (boost::is_same<T,int>::value) {
		size = 1;
		type = iBase_INTEGER;
	} else {
		size = (int)sizeof(T);
		type = iBase_BYTES;
	}
	tag = mesh_ptr->createTagHandle(name, size, type);
	// TODO: consider more transparent handling of exiting tag here
	tag = mesh_ptr->getTagHandle(name);
	entities = mesh_ptr->getEntities(entityDimension);
	currentElement = NULL;
	// TODO: nVertices should be standardized across functions
	int nVertices = 4;
	currentFields = std::vector<T>(nVertices);
}

template <class T>
T Field<T>::getField(Eigen::Vector3d position, iBase_EntityHandle *entity) {
	// TODO: this doesn't work for non-vertex fields
	T field;
	bool inElement=false;
	Eigen::Vector4d interpolationCoeffs;
	if (*entity==NULL) {
		// TODO: handle case with no entity hint
		std::vector<iBase_EntityHandle> adjacentEntities =
				mesh_ptr->getAdjacentEntities(entities[0],iBase_REGION);
		*entity = adjacentEntities[0];
	} else if (*entity==currentElement) {
//		int dimension=mesh_ptr->getEntityDimension(*entity);
//		if (dimension==iBase_REGION) {
		// TODO: should probably reset currentElement somehow at beginning
		//       of each orbit to prevent data from an old field being used
		//       (probably very unlikely since end tet of one orbit would have
		//       to be start tet of the next)
			inElement=true;
			interpolationCoeffs =
					mesh_ptr->getInterpolationCoeffs(position, currentElement);
			for (int i=0; i<interpolationCoeffs.rows(); i++)
				if (interpolationCoeffs[i]<0.)
					inElement=false;
//		} else {
//			std::cout << currentElement << ": " << dimension << " " <<
//					position.transpose() <<
//					" " << position.norm() << std::endl;
//		}
	}
	bool isElement;
	if (!inElement) {
		Eigen::Vector3d centroid(0.,0.,0.);
		int dimension=mesh_ptr->getEntityDimension(*entity);
		if (dimension!=iBase_REGION) {
			// TODO: find adjacent element
	//		*entity = entities[0];
			isElement = false;
		} else {
			isElement = true;
			std::vector<Eigen::Vector3d> vVs =
					mesh_ptr->getVertexVectors(*entity);
			centroid = (vVs[0]+vVs[1]+vVs[2]+vVs[3])/4.;
		}
		currentElement = *entity;
		bool foundTet=false;
		currentElement = mesh_ptr->findTet(centroid,
				position, currentElement, &foundTet, isElement);
		// TODO: failure to find tet doesn't really mean outside domain yet
		if (!foundTet)
			throw int(OUTSIDE_DOMAIN);
		// TODO: should ensure order of handles and vectors is same
		currentVertices = mesh_ptr->getVertices(currentElement);
		for (int i=0;i<currentVertices.size();i++) {
			currentFields[i] = this->getField(currentVertices[i]);
		}
		interpolationCoeffs =
				mesh_ptr->getInterpolationCoeffs(position, currentElement);
	}
	field *= 0.;
	assert(currentVertices.size()==interpolationCoeffs.rows());
	for (int i=0;i<interpolationCoeffs.rows();i++) {
		field += interpolationCoeffs[i]*currentFields[i];
	}
//	std::cout << field << std::endl;
	*entity = currentElement;
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


