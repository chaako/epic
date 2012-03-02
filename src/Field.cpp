/*
 * Field.cpp
 *
 *  Created on: Feb 10, 2012
 *      Author: chaako
 */

#include "epic.h"
#include "Field.h"

#include <boost/type_traits.hpp>

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
		int size, type;
		// TODO: is_same is in std from C++0x, so perhaps switch to
		//       avoid boost dependency
		if (boost::is_same<T,double>::value) {
			size = 1;
			type = iBase_DOUBLE;
		} else if (boost::is_same<T,int>::value) {
			size = 1;
			type = iBase_INTEGER;
		} else {
			size = (int)sizeof(T);
			type = iBase_BYTES;
		}
		iMesh_createTag(mesh_ptr->meshInstance, name.c_str(),
				size, type, &tag, &ierr, (int)name.length());
		CHECK("Failure creating tag");
	}
}


template <class T>
T Field<T>::getField(Eigen::Vector3d position) {
	T field;
	return field;
}


//template <>
//double Field<double>::getField(Eigen::Vector3d position) {
//	double field;
//
//	return field;
//}
//
//template <>
//double Field<double>::getField(iBase_EntityHandle node) {
//	double field;
//
//	return field;
//}
//
//
//template <>
//int Field<int>::getField(Eigen::Vector3d position) {
//	int field;
//
//	return field;
//}
//
//template <>
//int Field<int>::getField(iBase_EntityHandle node) {
//	int field;
//
//	return field;
//}
//
//
//template <>
//Eigen::Vector3d Field<Eigen::Vector3d>::getField(Eigen::Vector3d position) {
//	Eigen::Vector3d field;
//
//	return field;
//}
//
//template <>
//Eigen::Vector3d Field<Eigen::Vector3d>::getField(iBase_EntityHandle node) {
//	Eigen::Vector3d field;
//
//	return field;
//}

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

void DensityField::calcField() {}

void DensityField::calcField(ElectricField electricField) {
	int ierr;
	iBase_EntityHandle *ents0d = NULL;
	int ents0d_alloc = 0, ents0d_size;
	clock_t startClock = clock(); // timing
	// set potential value for all 0d elements
	iMesh_getEntities(mesh_ptr->meshInstance, mesh_ptr->rootEntitySet,
			iBase_VERTEX, iMesh_ALL_TOPOLOGIES,
			&ents0d, &ents0d_alloc, &ents0d_size, &ierr);
	CHECK("Couldn't get vertex entities");
	for (int i = 0; i < ents0d_size; i++) {
		double density=0.;
		Eigen::Vector3d nodePosition = mesh_ptr->getCoordinates(ents0d[i]);
		Eigen::Vector3d xzPosition = nodePosition;
		xzPosition[1] = 0;
		double r = xzPosition.norm();
		if ( r<0.2*nodePosition[1] && 0.<nodePosition[1] ) {
		IntegrandContainer integrandContainer;
		integrandContainer.mesh_ptr = mesh_ptr;
		integrandContainer.node = ents0d[i];
		integrandContainer.electricField_ptr = &electricField;
		int vdim=3;
		double xmin[vdim], xmax[vdim];
		for (int j=0; j<vdim; j++) {
			xmin[j] = -1.;
			xmax[j] = 1.;
		}
		double error=0.;
		// TODO: should make number of orbits adaptive
		adapt_integrate(1, &valueFromBoundary, (void*)&integrandContainer,
				vdim, xmin, xmax, 100, 1.e-5, 1.e-5, &density, &error);
		std::cout << "density[" << i << "] = " << density << ", error ="
				<< error << std::endl;
		}
		iMesh_setDblData(mesh_ptr->meshInstance, ents0d[i], tag, density,
				&ierr);
		CHECK("Failure setting potential tag");
	}
	clock_t endClock = clock(); // timing
	std::cout << "calcField total (s)= "
			<< (double)(endClock-startClock)/(double)CLOCKS_PER_SEC << std::endl; // timing
	std::cout << "findTet total (s)= "
			<< (double)extern_findTet/(double)CLOCKS_PER_SEC << std::endl; // timing
//	std::cout << "checkIfInNewTet (s)= "
//			<< (double)extern_checkIfInNewTet/(double)CLOCKS_PER_SEC << std::endl; // timing
}

