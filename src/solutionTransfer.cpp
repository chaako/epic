/*
 * solutionTransfer.cpp
 *
 *  Created on: Jan 27, 2013
 *      Author: chaako
 */

#include "epic.h"

extern "C" {void interpolatefieldssceptic3d_(
		const char*,double*,double*,double*,int&,
		double*,double*,double*,double*,
		double*,const int&);}

int main(int argc, char *argv[]) {
	if (argc<3) {
		printf("usage: %s hdfin meshin\n", argv[0]);
		exit(1);
	}

	cout.precision(3);
	string inputMeshFile(argv[2]);
//	if (inputMeshFile.find(".vtk")!=string::npos) {
	if (inputMeshFile.find(".sms")!=string::npos) {
		Mesh mesh(inputMeshFile);
		int numberOfVertices=mesh.entitiesVectors[iBase_VERTEX].size();
		double *xs = new double[numberOfVertices];
		double *ys = new double[numberOfVertices];
		double *zs = new double[numberOfVertices];
		double *potentials = new double[numberOfVertices];
		double *axs = new double[numberOfVertices];
		double *ays = new double[numberOfVertices];
		double *azs = new double[numberOfVertices];
		double *densities = new double[numberOfVertices];
		for (int i=0; i<numberOfVertices; i++) {
			entHandle entity=mesh.entitiesVectors[iBase_VERTEX][i];
			vect3d position=mesh.getCoordinates(entity);
			// Move nodes slightly to keep inside sceptic3D mesh
			// TODO: this is a little hacky
			double maxShift=0.00005;
			double shift=0.;
			double minimumRadius=1.;
			double unchangedRadius=1.5;
			double distanceToUnchangedR =
					unchangedRadius-position.norm();
			if (distanceToUnchangedR>0) {
				shift = maxShift*distanceToUnchangedR/
						(unchangedRadius-minimumRadius);
			} else {
				shift = maxShift*distanceToUnchangedR/
						position.norm();
			}
			double scaleFactor = 1. + shift/position.norm();
			position *= scaleFactor;
			xs[i] = position[0];
			ys[i] = position[1];
			zs[i] = position[2];
//			struct timespec sleepTime;
//			struct timespec returnTime;
//			sleepTime.tv_sec = 0;
//			sleepTime.tv_nsec = 10000000;
//			nanosleep(&sleepTime, &returnTime);
////			boost::this_thread::sleep( boost::posix_time::milliseconds(10) );
//			cout << "x y z: " << position.transpose() << endl;
//			cout << "r: " << position.norm() << endl;
		}
		string inputHdfFile(argv[1]);
		int filenameLength=inputHdfFile.length();
		interpolatefieldssceptic3d_(inputHdfFile.c_str(),xs,ys,zs,
				numberOfVertices,potentials,axs,ays,azs,
				densities,filenameLength);
//		cout << "first fields: " << potentials[0] << " "
//				<< axs[0] << " " << ays[0] << " " << azs[0] << endl;
		PotentialField potential(&mesh,string("potential"));
		DensityField ionDensity(&mesh,string("ionDensity"));
		ElectricField eField(&mesh,string("eField"));
		for (int i=0; i<numberOfVertices; i++) {
			entHandle entity=mesh.entitiesVectors[iBase_VERTEX][i];
			vect3d position=mesh.getCoordinates(entity);
			potential.setField(entity,potentials[i]);
			ionDensity.setField(entity,densities[i]);
			eField.setField(entity,vect3d(axs[i],ays[i],azs[i]));
		}
		int periodLocation = inputMeshFile.rfind(".vtk");
		stringstream outputFile;
		outputFile << inputMeshFile.substr(0,periodLocation)
				<< "_solutionTransfer" << ".sms";
		mesh.save(outputFile.str().c_str());
		delete[] xs;
		delete[] ys;
		delete[] zs;
		delete[] potentials;
		delete[] axs;
		delete[] ays;
		delete[] azs;
		delete[] densities;
	}
	return 0;
}
