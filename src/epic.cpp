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
//	std::cout << "  Process " << mpiId << " says 'Hello, world!'\n";
#endif

	Mesh mesh(argv[1]);
	mesh.printElementNumbers();

	Field<int> faceType(&mesh,std::string("cell_code"),iBase_FACE);
	CodeField vertexType(&mesh,std::string("vertex_code"),iBase_VERTEX);
	PotentialField potential(&mesh,std::string("potential"));
	ElectricField eField(&mesh,std::string("eField"));
	DensityField density(&mesh,std::string("density"));
	DensityField ionDensity(&mesh,std::string("ionDensity"));
	DensityField electronDensity(&mesh,std::string("electronDensity"));

	if (mpiId == 0)
		std::cout << std::endl << "Setting vertex codes..." << std::endl;
	vertexType.calcField(faceType);
	if (mpiId == 0)
		std::cout << std::endl << "Setting potential..." << std::endl;
	potential.calcField();

	for (int i=0; i<1; i++) {
		if (mpiId == 0)
			std::cout << std::endl << "Calculating electric field..." << std::endl;
		eField.calcField(potential);
		if (mpiId == 0)
			std::cout << std::endl << "Calculating electron density..." << std::endl;
		electronDensity.calcField(eField, potential, faceType, vertexType, -1.);
		if (mpiId == 0)
			std::cout << std::endl << "Calculating ion charge-density..." << std::endl;
		ionDensity.calcField(eField, potential, faceType, vertexType, 1.);
		if (mpiId == 0)
			std::cout << std::endl << "Calculating charge density..." << std::endl;
		density.calcField(ionDensity,electronDensity);
		if (mpiId == 0)
			std::cout << std::endl << "Calculating updated potential..." << std::endl;
		potential.calcField(ionDensity,electronDensity);
		if (mpiId == 0)
			std::cout << std::endl << std::endl << std::endl;
	}

	if (mpiId == 0)
		mesh.save(argv[2]);

#ifdef HAVE_MPI
	MPI::Finalize();
#endif

	return 0;
}

clock_t extern_findTet=0, extern_checkIfInNewTet=0;
