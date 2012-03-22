/*
 * Field.cpp
 *
 *  Created on: Feb 10, 2012
 *      Author: chaako
 */

#include "epic.h"
#include "Field.h"

ElectricField::ElectricField(Mesh *inputMesh_ptr, std::string inputName)
		: Field(inputMesh_ptr, inputName, iBase_VERTEX) {
}

void ElectricField::calcField(PotentialField potentialField) {
	for (int i=0; i<entities.size(); i++) {
		std::vector<iBase_EntityHandle> superCellFaces =
				getSuperCellFaces(mesh_ptr->meshInstance, entities[i]);
		Eigen::Vector3d eField(0.,0.,0.);
		Eigen::Vector3d point = mesh_ptr->getCoordinates(entities[i]);
		double volume=0;
		for (int j=0; j<superCellFaces.size(); j++)  {
			Eigen::Vector3d surfaceVector =
					getSurfaceVector(mesh_ptr->meshInstance, point,
							superCellFaces[j]);
			double potential = getAverageDblData(mesh_ptr->meshInstance,
					superCellFaces[j], potentialField.tag);
			volume += getTetVolume(mesh_ptr->meshInstance, point, superCellFaces[j]);
			eField -= potential*surfaceVector;
		}
		eField /= volume;
		this->setField(entities[i], eField);
	}
}


PotentialField::PotentialField(Mesh *inputMesh_ptr, std::string inputName)
	: Field(inputMesh_ptr, inputName, iBase_VERTEX) {
}

void PotentialField::calcField() {
	for (int i=0; i<entities.size(); i++) {
		Eigen::Vector3d point = mesh_ptr->getCoordinates(entities[i]);
		// TODO: Change this from just a test Coulomb field
		double potential = -1./point.norm();
		this->setField(entities[i], potential);
	}
}

void PotentialField::calcField(DensityField ionDensity,
		DensityField electronDensity) {
	assert(ionDensity.mesh_ptr == mesh_ptr);
	assert(electronDensity.mesh_ptr == mesh_ptr);
	for (int i=0; i<entities.size(); i++) {
		// TODO: remove this restriction on vertices
		Eigen::Vector3d nodePosition = mesh_ptr->getCoordinates(entities[i]);
		Eigen::Vector3d xzPosition = nodePosition;
		xzPosition[1] = 0;
		double r = xzPosition.norm();
		if ( r<0.2*nodePosition[1] && 0.<nodePosition[1] ) {

		double potential = this->getField(entities[i]);
		potential -= log(ionDensity.getField(entities[i])/
				electronDensity.getField(entities[i]));
		this->setField(entities[i], potential);

//		if (i==381 || i==2543 || i==2540 || i==1052 || i==1489 || i==1598 || i==1597 || i==3499)
			std::cout << "potential[" << i << "] = " << potential <<
					", 1/r= " << 1./nodePosition.norm() << std::endl;
		}
	}
}

DensityField::DensityField(Mesh *inputMesh_ptr, std::string inputName)
		: Field(inputMesh_ptr, inputName, iBase_VERTEX) {
}

void DensityField::calcField() {}

void DensityField::calcField(DensityField ionDensity,
		DensityField electronDensity) {
	for (int i=0; i<entities.size(); i++) {
		double density = ionDensity.getField(entities[i]) -
				electronDensity.getField(entities[i]);
		this->setField(entities[i], density);
	}
}

void DensityField::calcField(ElectricField electricField,
		PotentialField potentialField, double charge) {
	clock_t startClock = clock(); // timing
	for (int i=0; i<entities.size(); i++) {
		double density=0.;
		Eigen::Vector3d nodePosition = mesh_ptr->getCoordinates(entities[i]);
		Eigen::Vector3d xzPosition = nodePosition;
		xzPosition[1] = 0;
		double r = xzPosition.norm();
		if ( r<0.2*nodePosition[1] && 0.<nodePosition[1] ) {
		IntegrandContainer integrandContainer;
		integrandContainer.mesh_ptr = mesh_ptr;
		integrandContainer.node = entities[i];
		integrandContainer.electricField_ptr = &electricField;
		integrandContainer.potentialField_ptr = &potentialField;
		std::stringstream fileNameStream;
		fileNameStream << "distFunc/distributionFunction_r" << nodePosition.norm()
				<< "_vert" << i << ".p3d";
		integrandContainer.outFile = fopen(fileNameStream.str().c_str(), "w");
		integrandContainer.charge = charge;
		fprintf(integrandContainer.outFile, "x y z f\n");
		int vdim=3;
		double xmin[vdim], xmax[vdim];
		for (int j=0; j<vdim; j++) {
			xmin[j] = -1.;
			xmax[j] = 1.;
		}
		double error=0.;
		// TODO: should make number of orbits adaptive
		int numberOfOrbits=100;
//		if (i==381 || i==2543 || i==2540 || i==1052 || i==1489 || i==1598 || i==1597 || i==3499)
//			numberOfOrbits = 10000;
		adapt_integrate(1, &valueFromBoundary, (void*)&integrandContainer,
				vdim, xmin, xmax, numberOfOrbits, 1.e-5, 1.e-5, &density, &error);
//		std::cout << "density[" << i << "] = " << density << ", error ="
//				<< error << ", r = " << nodePosition.norm() << std::endl;
//		if (i==381 || i==2543 || i==2540 || i==1052 || i==1489 || i==1598 || i==1597 || i==3499) {
//			std::cout << "potential[" << i << "] = " << potentialField.getField(ents0d[i]) << std::endl;
//			std::cout << "eField[" << i << "] = " << std::endl;
//			std::cout << electricField.getField(ents0d[i]) << std::endl;
//		}
//		if (i==381 || i==2543 || i==2540 || i==1052 || i==1489 || i==1598 || i==1597 || i==3499)
			std::cout << nodePosition.norm() << " " << density << " " << error << std::endl;
		fclose(integrandContainer.outFile);
		}
		this->setField(entities[i], density);
	}
	clock_t endClock = clock(); // timing
	std::cout << "calcField total (s)= "
			<< (double)(endClock-startClock)/(double)CLOCKS_PER_SEC << std::endl; // timing
	std::cout << "findTet total (s)= "
			<< (double)extern_findTet/(double)CLOCKS_PER_SEC << std::endl; // timing
//	std::cout << "checkIfInNewTet (s)= "
//			<< (double)extern_checkIfInNewTet/(double)CLOCKS_PER_SEC << std::endl; // timing
}

