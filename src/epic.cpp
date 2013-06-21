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
	string inputSurfaceMeshFile, outputSurfaceMeshFile;
	string noSurfaceMeshFile("noSurfaceMeshFile.vtu");
	SurfaceMesh surfaceMesh;
	int evalSurfaceCode;
	vect3dMap vtkIdOfSurfacePoint(vect3dLessThan);
	extern_vtkIdOfSurfacePoint_ptr = &vtkIdOfSurfacePoint;
	vector<vect3d> evaluationPositions;
	extern_evalPositions_ptr = &evaluationPositions;
	bool doLuDecomposition, stopAfterEvalPos, doPoissonTest, fixSheathPotential;
	bool usePotentialFromInput;
	int secondsToSleepForDebugAttach, numberOfIterations;
	double debyeLength, boundaryPotential, surfacePotential, sheathPotential;
	double densityGradient, parallelDriftGradient;
	double parallelTemperatureGradient, perpendicularTemperatureGradient;
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
					("numberOfIterations", po::value<int>(&numberOfIterations)->default_value(2),
							"number of iterations")
					("magneticFieldStrength", po::value<double>(&extern_B[2])->default_value(0.),
							"magnetic field strength")
					("electricFieldStrength", po::value<double>(&extern_E[0])->default_value(0.),
							"electric field strength")
					("debyeLength", po::value<double>(&debyeLength)->default_value(0.),
							"electron Debye length")
					("boundaryPotential", po::value<double>(&boundaryPotential)->default_value(0.),
							"potential at outer boundary")
					("surfacePotential", po::value<double>(&surfacePotential)->default_value(-4.),
							"potential at object surface")
					("sheathPotential", po::value<double>(&sheathPotential)->default_value(-0.5),
							"potential at sheath entrance (for debyeLength=0.)")
					("fixSheathPotential", po::value<bool>(&fixSheathPotential)->default_value(false),
							"fix potential at sheath entrance (true/false; for debyeLength=0.)")
					("usePotentialFromInput", po::value<bool>(&usePotentialFromInput)->default_value(false),
							"use potential from input file (true/false)")
					("densityGradient", po::value<double>(&densityGradient)->default_value(0.),
							"density gradient")
					("parallelTemperatureGradient", po::value<double>(&parallelTemperatureGradient)->default_value(0.),
							"parallel temperature gradient")
					("perpendicularTemperatureGradient", po::value<double>(&perpendicularTemperatureGradient)->default_value(0.),
							"perpendicular temperature gradient")
					("parallelDriftGradient", po::value<double>(&parallelDriftGradient)->default_value(0.),
							"parallel drift gradient")
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
					("inputSurfaceMeshFile", po::value<string>(&inputSurfaceMeshFile)->default_value(string(noSurfaceMeshFile)),
							"input surface mesh file")
					("outputSurfaceMeshFile", po::value<string>(&outputSurfaceMeshFile)->default_value(string(noSurfaceMeshFile)),
							"output surface mesh file")
					("evalSurfaceCode", po::value<int>(&evalSurfaceCode)->default_value(4),
							"cell_code of surface to get evalPositions from if inputSurfaceMeshFile set")
					("saveOrbits", po::value<bool>(&extern_saveOrbits)->default_value(false),
							"whether to output orbits from evaluated nodes (true/false)")
			;

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);
			po::notify(vm);

			if (vm.count("help")) {
				if (mpiId==0)
					cout << desc << "\n";
				exit(1);
			}

			if (vm.count("inputFile")) {
				if (mpiId==0)
					cout << "Input file: " << inputMeshFile << endl;
			} else {
				if (mpiId==0)
					cout << "Error: --inputFile was not set" << endl;
				exit(1);
			}
			if (vm.count("outputFile")) {
				if (mpiId==0)
					cout << "Output file: " << outputFile << endl;
			} else {
				if (mpiId==0)
					cout << "Error: --outputFile was not set" << endl;
				exit(1);
			}

			if (inputSurfaceMeshFile != noSurfaceMeshFile) {
				surfaceMesh.load(inputSurfaceMeshFile);
				// TODO: what if surface of volume mesh has been refined compared to orig.?
				evaluationPositions = surfaceMesh.getPoints(&vtkIdOfSurfacePoint, evalSurfaceCode);
				extern_numberOfSurfaceEvalPoints = evaluationPositions.size();
				if (mpiId==0)
					cout << "Number of evaluation positions from surface mesh: " <<
							extern_numberOfSurfaceEvalPoints << endl;
			}

			int numberOfEvalPositions=min(min(evalPosX.size(),evalPosY.size()),evalPosZ.size());
			if (numberOfEvalPositions>0 && mpiId==0) {
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

	if (extern_B.norm()>0.)
		extern_VEXB = extern_E.cross(extern_B)/pow(extern_B.norm(),2.);

	try {
	if (secondsToSleepForDebugAttach>0) {
		char hostname[256];
		gethostname(hostname, sizeof(hostname));
		printf("mpiId %d: PID %d on %s ready for attach\n", mpiId, getpid(), hostname);
		// TODO: this flush doesn't appear to be working...optimization problem?
		fflush(stdout);
		sleep(secondsToSleepForDebugAttach);
	}

	Mesh mesh(inputMeshFile);
	mesh.printElementNumbers();

	// TODO: replace this with something more flexible/modular?
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
	if (mpiId == 0)
		cout << endl << "Setting vertex codes..." << endl;
	vertexType.calcField(faceType);
	// TODO: create error fields with pointer in corresponding fields
	PotentialField potential(&mesh,string("potential"));
	ElectricField eField(&mesh,string("eField"),vertexType,debyeLength,doLuDecomposition);
	Field<vect3d> ionVelocity(&mesh,string("ionVelocity"),iBase_VERTEX);
	Field<double> ionTemperature(&mesh,string("ionTemperature"),iBase_VERTEX);
	// TODO: updating other moments through density not very clean/transparent
	DensityField ionDensity(&mesh,string("ionDensity"),&ionVelocity,&ionTemperature);
	DensityField referenceElectronDensity(&mesh,string("referenceElectronDensity"));
//	DensityField electronDensity(&mesh,string("electronDensity"));
//	DensityField density(&mesh,string("density"));
//	DensityField ionDensityPositivePerturbation(&mesh,string("PPionDensity"));
//	DensityField ionDensityNegativePerturbation(&mesh,string("NPionDensity"));
//	DensityField electronDensityPositivePerturbation(&mesh,string("PPelectronDensity"));
//	DensityField electronDensityNegativePerturbation(&mesh,string("NPelectronDensity"));
	ShortestEdgeField shortestEdge(&mesh,string("shortestEdge"));

	SurfaceField<double,vtkDoubleArray,3,vect3d> ionSurfaceVelocity(&surfaceMesh, "ionSurfaceVelocity");

	double noPotentialPerturbation = 0.;
//	double positivePotentialPerturbation = 0.05;
//	double negativePotentialPerturbation = -0.05;

	if (mpiId == 0)
		cout << endl << "Calculating shortest edge of each region..." << endl;
	shortestEdge.calcField();
	if (mpiId == 0)
		cout << endl << "Setting density..." << endl;
//	// TODO: implement Boltzman density calculation;
//	electronDensity.calcField(potential, -1.);
//	density.calcField(ionDensity, electronDensity);
	boost::shared_ptr<SpatialDependence> densityProfile_ptr;
	boost::shared_ptr<SpatialDependence> parallelTemperatureProfile_ptr;
	boost::shared_ptr<SpatialDependence> perpendicularTemperatureProfile_ptr;
	boost::shared_ptr<SpatialDependence> parallelDriftProfile_ptr;
	if (densityGradient==0.) {
		densityProfile_ptr = boost::shared_ptr<SpatialDependence>(new SpatialDependence(1.));
	} else {
//		densityProfile_ptr =
//				boost::shared_ptr<SpatialDependence>(new ExponentialDependence(1./densityGradient));
		densityProfile_ptr =
				boost::shared_ptr<SpatialDependence>(new TanhDependence(1./densityGradient));
	}
	if (parallelTemperatureGradient==0.) {
		parallelTemperatureProfile_ptr =
				boost::shared_ptr<SpatialDependence>(new SpatialDependence(1.));
	} else {
//		parallelTemperatureProfile_ptr =
//				boost::shared_ptr<SpatialDependence>(
//						new ExponentialDependence(1./parallelTemperatureGradient));
		parallelTemperatureProfile_ptr =
				boost::shared_ptr<SpatialDependence>(new TanhDependence(1./parallelTemperatureGradient));
	}
	if (perpendicularTemperatureGradient==0.) {
		perpendicularTemperatureProfile_ptr =
				boost::shared_ptr<SpatialDependence>(new SpatialDependence(1.));
	} else {
//		perpendicularTemperatureProfile_ptr =
//				boost::shared_ptr<SpatialDependence>(
//						new ExponentialDependence(1./perpendicularTemperatureGradient));
		perpendicularTemperatureProfile_ptr =
				boost::shared_ptr<SpatialDependence>(new TanhDependence(1./perpendicularTemperatureGradient));
	}
	if (parallelDriftGradient==0.) {
		parallelDriftProfile_ptr =
				boost::shared_ptr<SpatialDependence>(new SpatialDependence(0.));
	} else {
//		parallelDriftProfile_ptr =
//				boost::shared_ptr<SpatialDependence>(
//						new ExponentialDependence(1./parallelDriftGradient));
		parallelDriftProfile_ptr =
				boost::shared_ptr<SpatialDependence>(
						new TanhDependence(1./parallelDriftGradient,vect3d(0.,0.,0.), 0., 1.));
	}
	vect3d magneticAxis = extern_B/extern_B.norm();
	// TODO: consider including drift in evaluated velocities
	vect3d perpendicularDrift(0.,0.,0.);
//	DistributionFunction distributionFunction;
	// TODO: pass smart pointers to maintain reference count?
	Maxwellian distributionFunction(*densityProfile_ptr.get(), magneticAxis,
			*parallelTemperatureProfile_ptr.get(), *perpendicularTemperatureProfile_ptr.get(),
			*parallelDriftProfile_ptr.get(), perpendicularDrift);
	referenceElectronDensity.calcField(vertexType, distributionFunction, 1.);
	potential.setReferenceElectronDensity(referenceElectronDensity);
	// TODO: Allow parallel electron temperature to differ from ions?
	potential.setReferenceElectronTemperature(*parallelTemperatureProfile_ptr.get());
	ionDensity.setDistributionFunction(distributionFunction);
	ionDensity.calcField(vertexType, potential, referenceElectronDensity,
			*parallelTemperatureProfile_ptr.get(), 1.);

	// TODO: add more robust detection and handling of existing fields
//	if (!mesh.vtkInputMesh && !doPoissonTest) {
	if (usePotentialFromInput) {
		if (mpiId == 0)
			cout << endl << "Using potential from input file." << endl;
//		if (mpiId == 0)
//			cout << endl << "Calculating electric field..." << endl;
//		eField.calcField(potential);
//		if (mpiId == 0)
//			cout << endl << "Calculating electron density..." << endl;
//		electronDensity.calcField(eField, potential, faceType, vertexType,
//				shortestEdge, -1., noPotentialPerturbation,
//				density_electronsFile);
	} else {
		if (mpiId == 0)
			cout << endl << "Setting potential..." << endl;
		potential.calcField(vertexType, debyeLength, boundaryPotential,
				surfacePotential, sheathPotential);
	}

	if (doPoissonTest) {
		if (mpiId == 0)
			cout << endl << "Setting poissonCubeTest density..." << endl;
		ionDensity.poissonCubeTest(debyeLength);
		if (mpiId == 0)
			cout << endl << "Calculating electric field..." << endl;
		eField.calcField(&potential, vertexType, ionDensity, debyeLength);
		if (mpiId == 0){
			int i = 0;
			stringstream iterMeshFileName;
			int periodLocation = outputFile.rfind(".");
			iterMeshFileName << outputFile.substr(0,periodLocation)
					<< setfill('0') << setw(2) << i << outputFile.substr(periodLocation);
//			// TODO: debugging
//			cout << iterMeshFileName.str() << endl;
			mesh.save(iterMeshFileName.str());
			// mesh.save() destroys vector tags, so update
			// TODO: do this automatically?
			eField.updateTagHandle();
			ionVelocity.updateTagHandle();
		}

		if (doPoissonTest) {
#ifdef HAVE_MPI
			MPI::Finalize();
#endif
			return(0);
		}

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

	for (int i=0; i<numberOfIterations; i++) {
		if (mpiId == 0)
			cout << endl  << endl << "ITERATION " << i << endl;
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
			cout << endl << "Calculating ion density..." << endl;
		ionDensity.calcField(eField, potential, referenceElectronDensity, faceType, vertexType,
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
		if (extern_numberOfSurfaceEvalPoints>0)
			ionSurfaceVelocity.copyFromField(ionVelocity);
		if (mpiId == 0)
			cout << endl << endl << endl;
		if (mpiId == 0){
			stringstream iterMeshFileName;
			int periodLocation = outputFile.rfind(".");
			iterMeshFileName << outputFile.substr(0,periodLocation)
					<< setfill('0') << setw(2) << i << outputFile.substr(periodLocation);
			mesh.save(iterMeshFileName.str());
			// mesh.save() destroys vector tags, so update
			// TODO: do this automatically?
			eField.updateTagHandle();
			ionVelocity.updateTagHandle();

			if (extern_numberOfSurfaceEvalPoints>0) {
				stringstream iterSurfaceMeshFileName;
				periodLocation = outputSurfaceMeshFile.rfind(".");
				iterSurfaceMeshFileName << outputSurfaceMeshFile.substr(0,periodLocation)
						<< setfill('0') << setw(2) << i << outputSurfaceMeshFile.substr(periodLocation);
				surfaceMesh.save(iterSurfaceMeshFileName.str());
			}
		}
		if (stopAfterEvalPos) {
			break;
		}
//		if (mpiId == 0)
//			cout << endl << "Saving current potential..." << endl;
//		stringstream potentialCopyName;
//		potentialCopyName << "potIter" << setfill('0') << setw(2) << i;
//		PotentialField potentialCopy(potential,potentialCopyName.str());
		if (debyeLength==0.) {
			if (mpiId == 0)
				cout << endl << "Calculating updated potential..." << endl;
			potential.calcField(ionDensity, vertexType, potentialFile, boundaryPotential,
					sheathPotential, fixSheathPotential);
//			potential.calcField(ionDensity, electronDensity, vertexType, potentialFile);
//			potential.calcField(ionDensity,
//					ionDensityPositivePerturbation, ionDensityNegativePerturbation,
//					electronDensity,
//					electronDensityPositivePerturbation, electronDensityNegativePerturbation,
//					vertexType, positivePotentialPerturbation,
//					negativePotentialPerturbation, potentialFile);
			if (mpiId == 0)
				cout << endl << "Calculating electric field from potential..." << endl;
//			// TODO: parallelize FEM eField-from-potential solve
//			eField.calcField(potential);
			eField.calcField_Gatsonis(potential);
		} else {
			if (mpiId == 0)
				cout << endl << "Calculating updated potential and electric field..." << endl;
			eField.calcField(&potential, vertexType, ionDensity, debyeLength);
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
	} catch(exception& e) {
		cerr << "error: " << e.what() << "\n";
	} catch(string& message) {
		cerr << "error: " << message << "\n";
	} catch(int code) {
		cerr << "error: " << code << "\n";
	} catch(...) {
		cerr << "Exception of unknown type!\n";
	}

#ifdef HAVE_MPI
	MPI::Finalize();
#endif

	return 0;
}
