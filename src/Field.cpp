/*
 * Field.cpp
 *
 *  Created on: Feb 10, 2012
 *      Author: chaako
 */

#include "epic.h"
#include "Field.h"

ElectricField::ElectricField(Mesh *inputMesh_ptr, string inputName)
		: Field<vect3d>(inputMesh_ptr, inputName, iBase_VERTEX) {
}

void ElectricField::calcField(PotentialField potentialField) {
	for (int i=0; i<entities.size(); i++) {
		vector<entHandle> superCellFaces =
				mesh_ptr->getSuperCellFaces(entities[i]);
		vect3d eField(0.,0.,0.);
		vect3d point = mesh_ptr->getCoordinates(entities[i]);
		double volume=0;
		for (int j=0; j<superCellFaces.size(); j++)  {
			vect3d surfaceVector =
					mesh_ptr->getSurfaceVector(superCellFaces[j], point);
			double potential = potentialField.getAverageField(superCellFaces[j]);
			volume += mesh_ptr->getTetVolume(point, superCellFaces[j]);
			eField += potential*surfaceVector;
		}
		eField /= volume;
		this->setField(entities[i], eField);
	}
}


PotentialField::PotentialField(Mesh *inputMesh_ptr, string inputName)
	: Field<double>(inputMesh_ptr, inputName, iBase_VERTEX) {
}

PotentialField::PotentialField(PotentialField potential, string inputName)
	: Field<double>(potential.mesh_ptr, inputName, iBase_VERTEX) {
	for (int i=0; i<entities.size(); i++) {
		this->setField(entities[i], potential.getField(entities[i]));
	}
}


void PotentialField::calcField() {
	for (int i=0; i<entities.size(); i++) {
		vect3d point = mesh_ptr->getCoordinates(entities[i]);
		// TODO: Change this from just a test Coulomb field
		double potential = -1./point.norm();
		this->setField(entities[i], potential);
	}
}

void PotentialField::calcField(DensityField ionDensity,
		DensityField electronDensity, CodeField vertexType, FILE *outFile) {
	assert(ionDensity.mesh_ptr == mesh_ptr);
	assert(electronDensity.mesh_ptr == mesh_ptr);
	for (int i=0; i<entities.size(); i++) {
		// TODO: remove this restriction on vertices
		vect3d nodePosition = mesh_ptr->getCoordinates(entities[i]);
		vect3d xzPosition = nodePosition;

		double potential;
		// TODO: don't hard-code boundary type and quasi-neutral operation
		if (vertexType.getField(entities[i])==4) {
			// TODO: don't hard-code sheath potential
			potential = -1./2.;
		} else {
			potential = this->getField(entities[i]);
			potential += 1./2.*log(ionDensity.getField(entities[i])/
					electronDensity.getField(entities[i]));
		}
		this->setField(entities[i], potential);

//#ifdef HAVE_MPI
//		if (MPI::COMM_WORLD.Get_rank() == 0)
//#endif
//			cout << nodePosition.norm() << " " << potential << endl;
		if (outFile)
			fprintf(outFile, "%g %g\n", nodePosition.norm(), potential);
//		}
	}
	if (outFile)
		fprintf(outFile, "\n\n\n\n");
}

DensityField::DensityField(Mesh *inputMesh_ptr, string inputName)
		: Field<double>(inputMesh_ptr, inputName, iBase_VERTEX) {
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
		PotentialField potentialField,
		Field<int> faceType, CodeField vertexType, double charge,
		FILE *outFile) {
	int mpiId = 0;
#ifdef HAVE_MPI
	double *density = new double[entities.size()];
	// TODO: Am assuming here that fields are identical on master and slaves
	mpiId = MPI::COMM_WORLD.Get_rank();
	if (mpiId == 0) {
		DensityField::requestDensityFromSlaves(electricField,
				potentialField, faceType, vertexType, charge, outFile);
		for (int node=0; node<entities.size(); node++) {
			density[node] = this->getField(entities[node]);
		}
	} else {
		DensityField::processDensityRequests(electricField,
				potentialField, faceType, vertexType, charge);
	}
	MPI::COMM_WORLD.Bcast(density, entities.size(), MPI::DOUBLE, 0);
	for (int node=0; node<entities.size(); node++) {
		this->setField(entities[node], density[node]);
	}
#else
	extern_findTet=0;
	clock_t startClock = clock(); // timing
	for (int node=0; node<entities.size(); node++) {
		double error;
	    double density = this->calculateDensity(node, electricField,
				potentialField, faceType, vertexType, charge, &error);
		this->setField(entities[node], density);
		vect3d nodePosition = mesh_ptr->getCoordinates(entities[node]);
//		cout << nodePosition.norm() << " " << density << " " << error << endl;
		if (outFile)
			fprintf(outFile, "%g %g %g\n", nodePosition.norm(), density, error);
	}
	clock_t endClock = clock(); // timing
	cout << "calcField total (s)= "
			<< (double)(endClock-startClock)/(double)CLOCKS_PER_SEC << endl; // timing
	cout << "findTet total (s)= "
			<< (double)extern_findTet/(double)CLOCKS_PER_SEC << endl; // timing
//	cout << "checkIfInNewTet (s)= "
//			<< (double)extern_checkIfInNewTet/(double)CLOCKS_PER_SEC << endl; // timing
#endif
	if (outFile)
		fprintf(outFile, "\n\n\n\n");
}

double DensityField::calculateDensity(int node, ElectricField electricField,
		PotentialField potentialField,
		Field<int> faceType, CodeField vertexType, double charge, double *error) {
	double density=0.;
	vect3d nodePosition = mesh_ptr->getCoordinates(entities[node]);
	IntegrandContainer integrandContainer;
	integrandContainer.mesh_ptr = mesh_ptr;
	integrandContainer.node = entities[node];
	integrandContainer.electricField_ptr = &electricField;
	integrandContainer.potentialField_ptr = &potentialField;
	integrandContainer.faceTypeField_ptr = &faceType;
	integrandContainer.vertexTypeField_ptr = &vertexType;
	stringstream fileNameStream;
//	fileNameStream << "distFunc/distributionFunction_r" << nodePosition.norm()
//			<< "_vert" << node << ".p3d";
	fileNameStream << "distFunc/distributionFunction_q" << charge << "_r"
			<< nodePosition.norm() << "_vert" << node << ".p3d";
	integrandContainer.outFile = NULL;
	stringstream fileNameStreamOrbit;
	fileNameStreamOrbit << "orbits/orbits_q" << charge << "_r"
			<< nodePosition.norm() << "_vert" << node << ".p3d";
	integrandContainer.orbitOutFile = NULL;
//	// TODO: don't hard-code output nodes
//	if (node<5) {
////	if (node==0) {
////	if (node==5 || node==2540) {
//		integrandContainer.outFile = fopen(fileNameStream.str().c_str(), "w");
//		fprintf(integrandContainer.outFile, "x y z f\n");
//		integrandContainer.orbitOutFile =
//				fopen(fileNameStreamOrbit.str().c_str(), "w");
//		fprintf(integrandContainer.orbitOutFile, "x y z energy\n");
//	}
	integrandContainer.charge = charge;
	int vdim=3;
	double xmin[vdim], xmax[vdim];
	for (int j=0; j<vdim; j++) {
		xmin[j] = -1.;
		xmax[j] = 1.;
	}
	// TODO: should make number of orbits adaptive
	int numberOfOrbits=100;
//	if (node==5 || node==2540)
//	if (node%1000==42)
//	if (node==0)
//	if (node<5)
	adapt_integrate(1, &distributionFunctionFromBoundary, (void*)&integrandContainer,
			vdim, xmin, xmax, numberOfOrbits, 1.e-5, 1.e-5, &density, error);
	if (integrandContainer.outFile)
		fclose(integrandContainer.outFile);
	if (integrandContainer.orbitOutFile)
		fclose(integrandContainer.orbitOutFile);
	return density;
}

#ifdef HAVE_MPI
void DensityField::requestDensityFromSlaves(ElectricField electricField,
		PotentialField potentialField,
		Field<int> faceType, CodeField vertexType, double charge,
		FILE *outFile) {
	int nProcesses = MPI::COMM_WORLD.Get_size();
	MPI::Status status;
	int node=0;

	// Send one node to each process
	for (int rank=1; rank<nProcesses; ++rank) {
		if (node<entities.size()) {
			MPI::COMM_WORLD.Send(&node, 1, MPI::INT, rank, WORKTAG);
			node++;
		}
	}

	// Process incoming density and send new nodes until all done
	while (node<entities.size()) {
		MPI::Status status = this->receiveDensity(outFile);
		MPI::COMM_WORLD.Send(&node, 1, MPI::INT, status.Get_source(),
				WORKTAG);
		node++;
	}

	// Process any outstanding densities
	// TODO: could run into problems here if fewer nodes than processes?
	assert(nProcesses<entities.size());
	for (int rank=1; rank<nProcesses; ++rank) {
		this->receiveDensity(outFile);
	}

	// Send empty message with DIETAG to signal done with nodes
	for (int rank=1; rank<nProcesses; ++rank) {
		MPI::COMM_WORLD.Send(0, 0, MPI::INT, rank, DIETAG);
	}
}
#endif

#ifdef HAVE_MPI
MPI::Status DensityField::receiveDensity(FILE *outFile) {
	MPI::Status status;
	double density[2];
	MPI::COMM_WORLD.Recv(&density[0], 2, MPI::DOUBLE, MPI_ANY_SOURCE,
			MPI_ANY_TAG, status);
	int incomingNode = status.Get_tag();
	if (incomingNode>=0 && incomingNode<entities.size())
		this->setField(entities[incomingNode], density[0]);
	vect3d nodePosition = mesh_ptr->getCoordinates(entities[incomingNode]);
//	cout << nodePosition.norm() << " " << density[0] << " " <<
//			density[1] << endl;
	if (outFile)
		fprintf(outFile, "%g %g %g\n", nodePosition.norm(), density[0], density[1]);
	return status;
}
#endif

#ifdef HAVE_MPI
void DensityField::processDensityRequests(ElectricField electricField,
		PotentialField potentialField,
		Field<int> faceType, CodeField vertexType, double charge=1.) {
	MPI::Status status;
	int node;

	while (1) {
		MPI::COMM_WORLD.Recv(&node, 1, MPI::INT, 0, MPI_ANY_TAG,
				status);
		if (status.Get_tag() == DIETAG) {
			return;
		}
		double density[2];
		density[0] = this->calculateDensity(node, electricField,
				potentialField, faceType, vertexType, charge, &density[1]);
		MPI::COMM_WORLD.Send(&density[0], 2, MPI::DOUBLE, 0, node);
	}
}
#endif

CodeField::CodeField(Mesh *inputMesh_ptr, string inputName, int elementType)
		: Field<int>(inputMesh_ptr, inputName, elementType) {
}

void CodeField::calcField(Field<int> faceTypeField) {
	for (int i=0; i<entities.size(); i++) {
		// TODO: this does not tag regions with an edge (but not whole face)
		//       on the boundary...could do 2ndAdjacent through vertex
		vector<entHandle> faces =
				mesh_ptr->getAdjacentEntities(entities[i],iBase_FACE);
		int elementType=0;
		for (int j=0; j<faces.size(); j++) {
			int faceType = faceTypeField.getField(faces[j]);
			// TODO: handle case where two adjacent faces are different boundaries?
			if (faceType>0)
				elementType=faceType;
		}
		this->setField(entities[i],elementType);
	}
}

