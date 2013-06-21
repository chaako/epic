/*
 * SurfaceField.h
 *
 *  Created on: Jun 20, 2013
 *      Author: chaako
 */

#ifndef SURFACEFIELD_H_
#define SURFACEFIELD_H_

#include "typesAndDefinitions.h"
#include "variables.h"
//#include <stdio.h>
////#include <type_traits> // Requires -std=c++0x compiler flag
//#include <boost/type_traits.hpp>
//
//#ifdef HAVE_MPI
//#include "mpi.h"
//#endif
//
//#include "iMesh.h"
//#include "Eigen/Dense"
#include "SurfaceMesh.h"
#include "Field.h"

template <class T, class vtkT, int N, class copyT=T>
class SurfaceField {
public:
	SurfaceField(SurfaceMesh *mesh_ptr, string name,
			int entityDimension);
	// TODO: Does an empty destructor cause a memory leak?
	virtual ~SurfaceField() {}

	void getField(vtkIdType node, T* field);
	void setField(vtkIdType node, T* field);

	void zero();

	void copyFromField(Field<copyT>& volumeField);

	SurfaceMesh *mesh_ptr;
	string name;
	int numberOfEntities;

	vtkSmartPointer<vtkT> field_ptr;

};

// gcc doesn't implement the export keyword, so define template functions here
template <class T, class vtkT, int N, class copyT>
SurfaceField<T,vtkT,N,copyT>::SurfaceField(SurfaceMesh *mesh_ptr, string name,
		int entityDimension=0) :
		mesh_ptr(mesh_ptr), name(name) {
	// TODO: don't assume entityDimension==0
	if (entityDimension!=0)
		throw string("nonzero entity dimension in SurfaceField not implemented yet");
	numberOfEntities = mesh_ptr->vtkMesh->GetNumberOfPoints();
	if (numberOfEntities>0) {
		// TODO: add try/catch?
		field_ptr = vtkSmartPointer<vtkT>::New();
		field_ptr->SetName(name.c_str());
		field_ptr->SetNumberOfComponents(N);
		field_ptr->SetNumberOfTuples(numberOfEntities);
		// TODO: verify that this doens't make a copy, i.e. can use field_ptr to modify
		mesh_ptr->vtkMesh->GetPointData()->AddArray(field_ptr);
		this->zero();
	}
}

template <class T, class vtkT, int N, class copyT>
void SurfaceField<T,vtkT,N,copyT>::getField(vtkIdType entity, T *field) {
	field_ptr->GetTuple(entity, field);
}

template <class T, class vtkT, int N, class copyT>
void SurfaceField<T,vtkT,N,copyT>::setField(vtkIdType entity, T *field) {
	field_ptr->SetTuple(entity, field);
}

template <class T, class vtkT, int N, class copyT>
void SurfaceField<T,vtkT,N,copyT>::zero() {
	T field[N];
	for (int i=0; i<N; i++)
		field[i] = 0;
	for (vtkIdType entity=0; entity<numberOfEntities; ++entity)
		this->setField(entity, field);
}

template <class T, class vtkT, int N, class copyT>
void SurfaceField<T,vtkT,N,copyT>::copyFromField(Field<copyT>& volumeField) {
	if (N!=3 || !boost::is_same<double,T>::value)
		throw string("trying to copy vect3d field to incompatible surfaceField");
	for (int i=0; i<volumeField.entities.size(); i++) {
		vect3d coordinates =
				volumeField.mesh_ptr->getCoordinates(volumeField.entities[i]);
		copyT values = volumeField[i];
		T *values_ptr;
		// TODO: this doesn't appear to be needed
		if (boost::is_same<copyT,vect3d>::value) {
			values_ptr = (T*)values.data();
		} else {
			values_ptr = (T*)&values;
		}
		vect3dMap::iterator lowerBound =
				mesh_ptr->idOfVertex.lower_bound(coordinates);
		vect3dMap::iterator upperBound =
				mesh_ptr->idOfVertex.upper_bound(coordinates);
		upperBound++;
		// TODO: make this more efficient by not checking all surface points for each point
		for (vect3dMap::iterator mapEntry=mesh_ptr->idOfVertex.begin();
				mapEntry!=mesh_ptr->idOfVertex.end(); mapEntry++) {
//		for (vect3dMap::iterator mapEntry=lowerBound; mapEntry!=upperBound; mapEntry++) {
			if ((coordinates-mapEntry->first).norm()<NODE_DISTANCE_THRESHOLD) {
				vtkIdType id = mapEntry->second;
//				// TODO: debugging
//				cout << id << "...";
				this->setField(id, values_ptr);
			}
		}
	}
}

#endif /* SURFACEFIELD_H_ */
