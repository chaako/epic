/*
 * solutionTransfer.cpp
 *
 *  Created on: Jan 27, 2013
 *      Author: chaako
 */

#include "epic.h"

extern "C" {void interpolatefieldssceptic3d_(
		const char*,double*,double*,double*,int&,
		double*,double*,double*,double*,const int&);}

int main(int argc, char *argv[]) {
	if (argc<3) {
		printf("usage: %s hdfin meshin\n", argv[0]);
		exit(1);
	}

	cout.precision(3);
	string inputMeshFile(argv[2]);
	if (inputMeshFile.find(".vtk")!=string::npos) {
		Mesh mesh(inputMeshFile);
		int numberOfVertices=mesh.entitiesVectors[iBase_VERTEX].size();
		double *xs = new double[numberOfVertices];
		double *ys = new double[numberOfVertices];
		double *zs = new double[numberOfVertices];
		double *potentials = new double[numberOfVertices];
		double *axs = new double[numberOfVertices];
		double *ays = new double[numberOfVertices];
		double *azs = new double[numberOfVertices];
		for (int i=0; i<numberOfVertices; i++) {
			entHandle entity=mesh.entitiesVectors[iBase_VERTEX][i];
			vect3d position=mesh.getCoordinates(entity);
			double scaleFactor=1.001;
			xs[i] = position[0]*scaleFactor;
			ys[i] = position[1]*scaleFactor;
			zs[i] = position[2]*scaleFactor;
//			cout << position.transpose() << endl;
		}
		string inputHdfFile(argv[1]);
		int filenameLength=inputHdfFile.length();
		interpolatefieldssceptic3d_(inputHdfFile.c_str(),xs,ys,zs,
				numberOfVertices,potentials,axs,ays,azs,filenameLength);
		cout << "first fields: " << potentials[0] << " "
				<< axs[0] << " " << ays[0] << " " << azs[0] << endl;
		PotentialField potential(&mesh,string("potential"));
		ElectricField eField(&mesh,string("eField"));
		for (int i=0; i<numberOfVertices; i++) {
			entHandle entity=mesh.entitiesVectors[iBase_VERTEX][i];
			vect3d position=mesh.getCoordinates(entity);
			potential.setField(entity,potentials[i]);
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
	}
	return 0;
}

clock_t extern_findTet=0, extern_checkIfInNewTet=0;

// TODO: find better way to distinguish orbits in output
int extern_orbitNumber = 0;
