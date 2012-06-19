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
		std::cout << "  The number of processes is "
				<< MPI::COMM_WORLD.Get_size() << "\n";
	}
	if (MPI::COMM_WORLD.Get_size()<2) {
		std::cout << "ERROR: MPI version currently requires minimum two processes" <<
				std::endl << "Run './configure --with-mpi=no' to compile serial version" <<
				std::endl;
		return 1;
	}
//	std::cout << "  Process " << mpiId << " says 'Hello, world!'\n";
#endif

	Mesh mesh(argv[1]);
	mesh.printElementNumbers();

	FILE *densityFile=NULL;
	FILE *density_electronsFile=NULL;
	FILE *potentialFile=NULL;

	if (mpiId == 0) {
		std::string fileName;
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

	Field<int> faceType(&mesh,std::string("cell_code"),iBase_FACE);
	CodeField vertexType(&mesh,std::string("vertex_code"),iBase_VERTEX);
	PotentialField potential(&mesh,std::string("potential"));
	ElectricField eField(&mesh,std::string("eField"));
	DensityField density(&mesh,std::string("density"));
	DensityField ionDensity(&mesh,std::string("ionDensity"));
	DensityField electronDensity(&mesh,std::string("electronDensity"));

	// TODO: add more robust detection and handling of existing fields
	if (!mesh.vtkInputMesh) {
		if (mpiId == 0)
			std::cout << std::endl << "Calculating electric field..." << std::endl;
		eField.calcField(potential);
		if (mpiId == 0)
			std::cout << std::endl << "Calculating electron density..." << std::endl;
		electronDensity.calcField(eField, potential, faceType, vertexType, -1.,
				density_electronsFile);
		if (mpiId == 0)
			std::cout << std::endl << "Calculating ion charge-density..." << std::endl;
		ionDensity.calcField(eField, potential, faceType, vertexType, 1.,
				densityFile);
		// TODO: shouldn't return before closing files etc...
		return 0;
	}

	if (mpiId == 0)
		std::cout << std::endl << "Setting vertex codes..." << std::endl;
	vertexType.calcField(faceType);
	if (mpiId == 0)
		std::cout << std::endl << "Setting potential..." << std::endl;
	potential.calcField();

//	// Integrate a circular test orbit
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
//				faceType, vertexType, outFile);
//		fclose(outFile);
//
//		return 0;
//	}

	for (int i=0; i<1; i++) {
		if (mpiId == 0)
			std::cout << std::endl  << std::endl << "ITERATION " << i << std::endl;
		if (mpiId == 0)
			std::cout << std::endl << "Calculating electric field..." << std::endl;
		eField.calcField(potential);
		if (mpiId == 0)
			std::cout << std::endl << "Calculating electron density..." << std::endl;
		electronDensity.calcField(eField, potential, faceType, vertexType, -1.,
				density_electronsFile);
		if (mpiId == 0)
			std::cout << std::endl << "Calculating ion charge-density..." << std::endl;
		ionDensity.calcField(eField, potential, faceType, vertexType, 1.,
				densityFile);
		if (mpiId == 0)
			std::cout << std::endl << "Calculating charge density..." << std::endl;
		density.calcField(ionDensity, electronDensity);
		if (mpiId == 0)
			std::cout << std::endl << "Saving current potential..." << std::endl;
		std::stringstream potentialCopyName;
		potentialCopyName << "potIter" << std::setfill('0') << std::setw(2) << i;
		PotentialField potentialCopy(potential,potentialCopyName.str());
		if (mpiId == 0)
			std::cout << std::endl << "Calculating updated potential..." << std::endl;
		potential.calcField(ionDensity, electronDensity, vertexType, potentialFile);
		if (mpiId == 0)
			std::cout << std::endl << std::endl << std::endl;
	}

	if (mpiId == 0)
		mesh.save(argv[2]);

	if (mpiId == 0) {
		fclose(densityFile);
		fclose(density_electronsFile);
		fclose(potentialFile);
	}

#ifdef HAVE_MPI
	MPI::Finalize();
#endif

	return 0;
}

clock_t extern_findTet=0, extern_checkIfInNewTet=0;
