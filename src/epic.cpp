/*
 ============================================================================
 Name        : epic.cpp
 Author      : Christian Bernt Haakonsen
 Version     :
 Copyright   : Copyright
 Description : Hello World in C++,
 ============================================================================
 */

#include "epic.h"

using namespace std;

int main(int argc, char *argv[]) {

	if (argc<3) {
		printf("usage: %s meshin meshout\n", argv[0]);
		exit(1);
	}

	int mpiId = 0;
#ifdef HAVE_MPI
	MPI::Init(argc, argv);
	mpiId = MPI::COMM_WORLD.Get_rank();

	if (mpiId == 0) {
		cout << "  The number of processes is "
				<< MPI::COMM_WORLD.Get_size() << "\n";
	}
	if (MPI::COMM_WORLD.Get_size()<2) {
		cout << "ERROR: MPI version currently requires minimum two processes" <<
				endl << "Run './configure --with-mpi=no' to compile serial version" <<
				endl;
		return 1;
	}
//	cout << "  Process " << mpiId << " says 'Hello, world!'\n";
#endif

	Mesh mesh(argv[1]);
	mesh.printElementNumbers();

	FILE *densityFile=NULL;
	FILE *density_electronsFile=NULL;
	FILE *potentialFile=NULL;

	if (mpiId == 0) {
		string fileName;
		fileName = "density.dat";
		densityFile = fopen(fileName.c_str(), "w");
		fprintf(densityFile, "# r n_i dn_i\n");
		fileName = "density_electrons.dat";
		density_electronsFile = fopen(fileName.c_str(), "w");
		fprintf(density_electronsFile, "# r n_e dn_e\n");
		fileName = "potential.dat";
		potentialFile = fopen(fileName.c_str(), "w");
		fprintf(potentialFile, "# r phi\n");
	}

	Field<int> faceType(&mesh,string("cell_code"),iBase_FACE);
	CodeField vertexType(&mesh,string("vertex_code"),iBase_VERTEX);
	PotentialField potential(&mesh,string("potential"));
	ElectricField eField(&mesh,string("eField"));
	DensityField density(&mesh,string("density"));
	DensityField ionDensity(&mesh,string("ionDensity"));
	DensityField electronDensity(&mesh,string("electronDensity"));
	DensityField ionDensityPositivePerturbation(&mesh,string("PPionDensity"));
	DensityField ionDensityNegativePerturbation(&mesh,string("NPionDensity"));
	DensityField electronDensityPositivePerturbation(&mesh,string("PPelectronDensity"));
	DensityField electronDensityNegativePerturbation(&mesh,string("NPelectronDensity"));
	ShortestEdgeField shortestEdge(&mesh,string("shortestEdge"));

	double noPotentialPerturbation = 0.;
	double positivePotentialPerturbation = 0.05;
	double negativePotentialPerturbation = -0.05;

	// TODO: add more robust detection and handling of existing fields
	if (!mesh.vtkInputMesh) {
		if (mpiId == 0)
			cout << endl << "Calculating electric field..." << endl;
		eField.calcField(potential);
		if (mpiId == 0)
			cout << endl << "Calculating electron density..." << endl;
		electronDensity.calcField(eField, potential, faceType, vertexType,
				shortestEdge, -1., noPotentialPerturbation,
				density_electronsFile);
		if (mpiId == 0)
			cout << endl << "Calculating ion charge-density..." << endl;
		ionDensity.calcField(eField, potential, faceType, vertexType,
				shortestEdge, 1., noPotentialPerturbation,
				densityFile);
		// TODO: shouldn't return before closing files etc...
		cout << "Not in main iteration loop...improve handling of existing fields." << endl;
		return 0;
	}

	if (mpiId == 0)
		cout << endl << "Setting vertex codes..." << endl;
	vertexType.calcField(faceType);
	if (mpiId == 0)
		cout << endl << "Setting potential..." << endl;
	potential.calcField(vertexType);
	if (mpiId == 0)
		cout << endl << "Calculating shortest edge of each region..." << endl;
	shortestEdge.calcField();

//	// Integrate a circular test orbit (need to deactivate trapped orbit rejection)
//	{
//		entHandle node = ionDensity.entities[2540];
//		double nodePotential = potential.getField(node);
//		vect3d position = mesh.getCoordinates(node);
//		vect3d velocity, zHat(0.,0.,1.);
//		velocity = position.cross(zHat);
//		velocity /= velocity.norm();
//		velocity *= sqrt(-nodePotential);
//		Orbit orbit(&mesh,node,velocity,1.);
//		const char *fName = "integratedOrbitTest.p3d";
//		FILE* outFile = fopen(fName, "w");
//		fprintf(outFile, "x y z energy\n");
//		orbit.integrate(potential, eField,
//				faceType, vertexType, shortestEdge, outFile);
//		fclose(outFile);
//
//		return 0;
//	}

	for (int i=0; i<1; i++) {
		if (mpiId == 0)
			cout << endl  << endl << "ITERATION " << i << endl;
		if (mpiId == 0)
			cout << endl << "Calculating electric field..." << endl;
		eField.calcField(potential);
//		if (mpiId == 0)
//			cout << endl << "Calculating electron density..." << endl;
//		electronDensity.calcField(eField, potential, faceType, vertexType,
//				shortestEdge, -1., noPotentialPerturbation,
//				density_electronsFile);
//		if (mpiId == 0)
//			cout << endl << "Calculating PP electron density..." << endl;
//		electronDensityPositivePerturbation.calcField(eField,
//				potential, faceType, vertexType,
//				shortestEdge, -1., positivePotentialPerturbation,
//				density_electronsFile);
//		if (mpiId == 0)
//			cout << endl << "Calculating NP electron density..." << endl;
//		electronDensityNegativePerturbation.calcField(eField,
//				potential, faceType, vertexType,
//				shortestEdge, -1., negativePotentialPerturbation,
//				density_electronsFile);
		if (mpiId == 0)
			cout << endl << "Calculating ion charge-density..." << endl;
		ionDensity.calcField(eField, potential, faceType, vertexType,
				shortestEdge, 1., noPotentialPerturbation,
				densityFile);
//		if (mpiId == 0)
//			cout << endl << "Calculating PP ion charge-density..." << endl;
//		ionDensityPositivePerturbation.calcField(eField, potential,
//				faceType, vertexType,
//				shortestEdge, 1., positivePotentialPerturbation,
//				densityFile);
//		if (mpiId == 0)
//			cout << endl << "Calculating NP ion charge-density..." << endl;
//		ionDensityNegativePerturbation.calcField(eField, potential,
//				faceType, vertexType,
//				shortestEdge, 1., negativePotentialPerturbation,
//				densityFile);
//		if (mpiId == 0)
//			cout << endl << "Calculating charge density..." << endl;
//		density.calcField(ionDensity, electronDensity);
		if (mpiId == 0)
			cout << endl << "Saving current potential..." << endl;
		stringstream potentialCopyName;
		potentialCopyName << "potIter" << setfill('0') << setw(2) << i;
		PotentialField potentialCopy(potential,potentialCopyName.str());
		if (mpiId == 0)
			cout << endl << "Calculating updated potential..." << endl;
		potential.calcField(ionDensity, vertexType, potentialFile);
//		potential.calcField(ionDensity, electronDensity, vertexType, potentialFile);
//		potential.calcField(ionDensity,
//				ionDensityPositivePerturbation, ionDensityNegativePerturbation,
//				electronDensity,
//				electronDensityPositivePerturbation, electronDensityNegativePerturbation,
//				vertexType, positivePotentialPerturbation,
//				negativePotentialPerturbation, potentialFile);
		if (mpiId == 0)
			cout << endl << endl << endl;
	}

	if (mpiId == 0)
		mesh.save(argv[2]);

	if (mpiId == 0) {
		if (densityFile)
			fclose(densityFile);
		if (density_electronsFile)
			fclose(density_electronsFile);
		if (potentialFile)
			fclose(potentialFile);
	}

#ifdef HAVE_MPI
	MPI::Finalize();
#endif

	return 0;
}

clock_t extern_findTet=0, extern_checkIfInNewTet=0;
