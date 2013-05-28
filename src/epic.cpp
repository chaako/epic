/*
 ============================================================================
 Name        : epic.cpp
 Author      : Christian Bernt Haakonsen
 Version     :
 Copyright   : Copyright
 Description : Hello World in C++,
 ============================================================================
 */

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "epic.h"

int main(int argc, char *argv[]) {

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

	string inputMeshFile, outputFile;
	vector<vect3d> evaluationPositions;
	extern_evalPositions_ptr = &evaluationPositions;
	bool doLuDecomposition, stopAfterEvalPos, doPoissonTest;
	int secondsToSleepForDebugAttach;
	// TODO: make names more transparent?
	{
		namespace po = boost::program_options;
		try {
			vector<double> evalPosX, evalPosY, evalPosZ;
			po::options_description desc("Allowed options");
			desc.add_options()
					("help", "produce help message")
					("inputFile", po::value<string>(&inputMeshFile), "input file")
					("outputFile", po::value<string>(&outputFile), "output file")
					("evalPositionX,x", po::value< vector<double> >(&evalPosX), "evaluation position(s) x")
					("evalPositionY,y", po::value< vector<double> >(&evalPosY), "evaluation position(s) y")
					("evalPositionZ,z", po::value< vector<double> >(&evalPosZ), "evaluation position(s) z")
					("doLuDecomposition", po::value<bool>(&doLuDecomposition)->default_value(true),
							"whether to do LU decomposition (true/false)")
					("stopAfterEvalPos", po::value<bool>(&stopAfterEvalPos)->default_value(false),
							"whether to stop after evaluating at positions (true/false)")
					("secondsToSleepForDebugAttach",
							po::value<int>(&secondsToSleepForDebugAttach)->default_value(0),
							"seconds to sleep to allow debugger to attach")
					("doPoissonTest", po::value<bool>(&doPoissonTest)->default_value(false),
							"whether to run in Poisson-test mode (true/false)")
			;

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);
			po::notify(vm);

			if (vm.count("help")) {
				cout << desc << "\n";
				exit(1);
			}

			if (vm.count("inputFile")) {
				cout << "Input file: " << inputMeshFile << endl;
			} else {
				cout << "Error: --inputFile was not set" << endl;
				exit(1);
			}
			if (vm.count("outputFile")) {
				cout << "Output file: " << outputFile << endl;
			} else {
				cout << "Error: --outputFile was not set" << endl;
				exit(1);
			}

			int numberOfEvalPositions=min(min(evalPosX.size(),evalPosY.size()),evalPosZ.size());
			if (numberOfEvalPositions>0) {
				cout << "Number of evaluation positions requested: " <<
						numberOfEvalPositions << endl;
				for (int i=0; i<numberOfEvalPositions; i++)
					evaluationPositions.push_back(vect3d(evalPosX[i],evalPosY[i],evalPosZ[i]));
			}
		} catch(exception& e) {
			cerr << "error: " << e.what() << "\n";
			exit(1);
		} catch(...) {
			cerr << "Exception of unknown type!\n";
		}
	}

	Mesh mesh(inputMeshFile);
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

    if (secondsToSleepForDebugAttach>0) {
		char hostname[256];
		gethostname(hostname, sizeof(hostname));
		printf("mpiId %d: PID %d on %s ready for attach\n", mpiId, getpid(), hostname);
		fflush(stdout);
		sleep(secondsToSleepForDebugAttach);
	}

	Field<int> faceType(&mesh,string("cell_code"),iBase_FACE);
	CodeField vertexType(&mesh,string("vertex_code"),iBase_VERTEX);
	if (mpiId == 0)
		cout << endl << "Setting vertex codes..." << endl;
	vertexType.calcField(faceType);
	PotentialField potential(&mesh,string("potential"));
	ElectricField eField(&mesh,string("eField"),vertexType,doLuDecomposition);
	Field<vect3d> ionVelocity(&mesh,string("ionVelocity"),iBase_VERTEX);
	Field<double> ionTemperature(&mesh,string("ionTemperature"),iBase_VERTEX);
	// TODO: updating other moments through density not very clean/transparent
	DensityField density(&mesh,string("density"));
	DensityField ionDensity(&mesh,string("ionDensity"),&ionVelocity,&ionTemperature);
	DensityField electronDensity(&mesh,string("electronDensity"));
	DensityField ionDensityPositivePerturbation(&mesh,string("PPionDensity"));
	DensityField ionDensityNegativePerturbation(&mesh,string("NPionDensity"));
	DensityField electronDensityPositivePerturbation(&mesh,string("PPelectronDensity"));
	DensityField electronDensityNegativePerturbation(&mesh,string("NPelectronDensity"));
	ShortestEdgeField shortestEdge(&mesh,string("shortestEdge"));

	double noPotentialPerturbation = 0.;
	double positivePotentialPerturbation = 0.05;
	double negativePotentialPerturbation = -0.05;

	if (mpiId == 0)
		cout << endl << "Calculating shortest edge of each region..." << endl;
	shortestEdge.calcField();
	if (mpiId == 0)
		cout << endl << "Setting density..." << endl;
	ionDensity.calcField(vertexType, potential, 1.);


	// TODO: add more robust detection and handling of existing fields
	if (!mesh.vtkInputMesh && !doPoissonTest) {
		if (mpiId == 0)
			cout << endl << ".sms file loaded, so assuming existing fields" << endl;
//		if (mpiId == 0)
//			cout << endl << "Calculating electric field..." << endl;
//		eField.calcField(potential);
//		if (mpiId == 0)
//			cout << endl << "Calculating electron density..." << endl;
//		electronDensity.calcField(eField, potential, faceType, vertexType,
//				shortestEdge, -1., noPotentialPerturbation,
//				density_electronsFile);
		if (mpiId == 0)
			cout << endl << "Calculating ion charge-density..." << endl;
		ionDensity.calcField(eField, potential, faceType, vertexType,
				shortestEdge, 1., noPotentialPerturbation,
				densityFile);
		if (stopAfterEvalPos) {
			// TODO: shouldn't return before closing files etc...
			cout << "Not in main iteration loop...improve handling of existing fields." << endl;
			return 0;
		}
	} else {
		if (mpiId == 0)
			cout << endl << "Setting potential..." << endl;
		potential.calcField(vertexType);
	}

	if (doPoissonTest) {
		if (mpiId == 0)
			cout << endl << "Setting poissonCubeTest density..." << endl;
		ionDensity.poissonCubeTest();
		if (mpiId == 0)
			cout << endl << "Calculating electric field..." << endl;
		eField.calcField(&potential, vertexType, ionDensity);
	}


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

	if (mpiId == 0){
		int i = 0;
		stringstream iterMeshFileName;
		int periodLocation = outputFile.rfind(".");
		iterMeshFileName << outputFile.substr(0,periodLocation)
				<< setfill('0') << setw(2) << i << outputFile.substr(periodLocation);
//		// TODO: debugging
//		cout << iterMeshFileName.str() << endl;
		mesh.save(iterMeshFileName.str());
		// mesh.save() destroys eField tag, so update
		eField.updateTagHandle();
	}

	// Exit for Poisson test
	if (doPoissonTest) {
#ifdef HAVE_MPI
		MPI::Finalize();
#endif
		return(0);
	}

	for (int i=1; i<2; i++) {
		if (mpiId == 0)
			cout << endl  << endl << "ITERATION " << i << endl;
		if (mpiId == 0)
			cout << endl << "Calculating electric field..." << endl;
		eField.calcField(&potential, vertexType, ionDensity);
//		eField.calcField(potential);
//		eField.calcField_Gatsonis(potential);
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
//		if (mpiId == 0)
//			cout << endl << "Saving current potential..." << endl;
//		stringstream potentialCopyName;
//		potentialCopyName << "potIter" << setfill('0') << setw(2) << i;
//		PotentialField potentialCopy(potential,potentialCopyName.str());
//		if (mpiId == 0)
//			cout << endl << "Calculating updated potential..." << endl;
//		potential.calcField(ionDensity, vertexType, potentialFile);
//		potential.calcField(ionDensity, electronDensity, vertexType, potentialFile);
//		potential.calcField(ionDensity,
//				ionDensityPositivePerturbation, ionDensityNegativePerturbation,
//				electronDensity,
//				electronDensityPositivePerturbation, electronDensityNegativePerturbation,
//				vertexType, positivePotentialPerturbation,
//				negativePotentialPerturbation, potentialFile);
		if (mpiId == 0)
			cout << endl << endl << endl;
		if (mpiId == 0){
			stringstream iterMeshFileName;
			int periodLocation = outputFile.rfind(".");
			iterMeshFileName << outputFile.substr(0,periodLocation)
					<< setfill('0') << setw(2) << i << outputFile.substr(periodLocation);
			mesh.save(iterMeshFileName.str());
			// mesh.save() destroys eField tag, so update
			eField.updateTagHandle();
		}
	}

	if (mpiId == 0)
		mesh.save(outputFile);

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
