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

	string inputMeshFile, inputDirectory, outputFile;
	vector<string> inputMeshFiles;
	vector<boost::filesystem::path> inputPaths;
	vector<int> inputIterations;
	string inputSurfaceMeshFile, outputSurfaceMeshFile;
	string noSurfaceMeshFile("noSurfaceMeshFile.vtu");
	SurfaceMesh surfaceMesh;
	int evalSurfaceCode;
	vect3dMap vtkIdOfSurfacePoint(vect3dLessThan);
	extern_vtkIdOfSurfacePoint_ptr = &vtkIdOfSurfacePoint;
	vector<vect3d> evaluationPositions;
	extern_evalPositions_ptr = &evaluationPositions;
	bool doLuDecomposition, stopAfterEvalPos, doPoissonTest, fixSheathPotential;
	bool usePotentialFromInput, useDensityFromInput;
	int secondsToSleepForDebugAttach, numberOfIterations, numberOfIterationsToAveragePotentialOver;
	int numberOfPotentialValues;
	double debyeLength, boundaryPotential, objectPotential, sheathPotential;
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
					("inputDirectory", po::value<string>(&inputDirectory), "input directory to continue run from")
					("outputFile", po::value<string>(&outputFile), "output file")
					("numberOfIterations", po::value<int>(&numberOfIterations)->default_value(2),
							"number of iterations")
					("numberOfIterationsToAveragePotentialOver",
							po::value<int>(&numberOfIterationsToAveragePotentialOver)->default_value(1),
							"number of iterations to average potential over")
					("magneticFieldStrength", po::value<double>(&extern_B[2])->default_value(0.),
							"magnetic field strength")
					("electricFieldStrength", po::value<double>(&extern_E[0])->default_value(0.),
							"electric field strength")
					("debyeLength", po::value<double>(&debyeLength)->default_value(0.),
							"electron Debye length")
					("boundaryPotential", po::value<double>(&boundaryPotential)->default_value(0.),
							"potential at outer boundary")
					("objectPotential", po::value<double>(&objectPotential)->default_value(-4.),
							"potential at object surface")
					("sheathPotential", po::value<double>(&sheathPotential)->default_value(-0.5),
							"potential at sheath entrance (for debyeLength=0.)")
					("fixSheathPotential", po::value<bool>(&fixSheathPotential)->default_value(false),
							"fix potential at sheath entrance (true/false; for debyeLength=0.)")
					("usePotentialFromInput", po::value<bool>(&usePotentialFromInput)->default_value(false),
							"use potential from input file (true/false)")
					("numberOfPotentialValues", po::value<int>(&numberOfPotentialValues)->default_value(1),
							"number of potential values to evaluate density at")
					("useDensityFromInput", po::value<bool>(&useDensityFromInput)->default_value(false),
							"use density from input file (true/false)")
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
			} else if (vm.count("inputDirectory")) {
				if (mpiId==0)
					cout << "Continuing run from directory: " << inputDirectory << endl;
				boost::filesystem::path p(inputDirectory);

				try {
					if (boost::filesystem::exists(p)) {
						if (boost::filesystem::is_directory(p)) {
							vector<boost::filesystem::path> v;
							copy(boost::filesystem::directory_iterator(p),
									boost::filesystem::directory_iterator(),
									back_inserter(v));
							// TODO: this only works if all .sms files have same bareFilename
							sort(v.begin(), v.end());
							for (vector<boost::filesystem::path>::const_iterator
									it(v.begin()); it!=v.end(); ++it) {
								string fileExtension = it->extension().string();
								if (fileExtension.compare(".sms")==0) {
									inputMeshFile = it->string();
									// TODO: message about this or allow input override?
									usePotentialFromInput = true;
									useDensityFromInput = true;
									string bareFilename = it->stem().string();
									boost::regex e("(\\D+)(\\d{2})");
									boost::smatch what;
									if(boost::regex_match(bareFilename, what, e)) {
										int iterationNumber;
										istringstream(what[2]) >> iterationNumber;
										inputMeshFiles.push_back(it->string());
										inputPaths.push_back(*it);
										inputIterations.push_back(iterationNumber);
										if (mpiId==0) {
											cout << what[1] << setfill('0') << setw(2) << iterationNumber <<
													setfill(' ') << setw(0) << fileExtension
													<< " taken as existing iteration" << endl;
										}
										outputFile = what[1] + fileExtension;
									} else {
										if (mpiId==0)
											cout << "No iteration number found in " <<
											it->filename() << endl;
									}
								}
							}
						}
					} else {
						if (mpiId==0)
							cout << p << " does not exist" << endl;
					}
				} catch (const boost::filesystem::filesystem_error& ex) {
					cout << ex.what() << endl;
				}
			} else {
				if (mpiId==0)
					cout << "Error: neither --inputFile nor --inputDirectory was set" << endl;
				exit(1);
			}
			if (vm.count("outputFile")) {
				if (mpiId==0)
					cout << "Output file: " << outputFile << endl;
			} else if (vm.count("inputDirectory")) {
				if (mpiId==0)
					cout << "Output file set by --inputDirectory: " << outputFile << endl;
			} else {
				if (mpiId==0)
					cout << "Error: neither --outputFile nor --inputDirectory was set" << endl;
				exit(1);
			}

			// TODO: should be able to output surface mesh during normal run (i.e. not with evalPositions)
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

	// TODO: set random seed as input or based on time(NULL)?
	srand(54326);

	vector<Eigen::VectorXd> diisResiduals;
	vector<Eigen::VectorXd> diisPotentials;

	int lastInputIteration=0;
	for (int i=0; i<inputIterations.size(); i++) {
		int numberOfNodes=-1;
		if (mpiId == 0) {
			boost::filesystem::path cwd(".");
			boost::filesystem::path cwdPath = boost::filesystem::system_complete(cwd);
			// TODO: better way to strip off "."/get full path?
			try {
				boost::filesystem::path copyTarget(cwdPath.parent_path().string()+"/"+inputPaths[i].filename().string());
				boost::filesystem::copy_file(inputMeshFiles[i],copyTarget);
				boost::filesystem::path copyTargetVtu(cwdPath.parent_path().string()+"/"+inputPaths[i].stem().string() + ".vtu");
//				boost::filesystem::path vtuExtension(".vtu");
				boost::filesystem::path inputFileVtu(inputMeshFiles[i]);
				inputFileVtu.replace_extension(".vtu");
				boost::filesystem::copy_file(inputFileVtu,copyTargetVtu);
				cout << "Copied " << copyTarget.filename() << " and " <<
						copyTargetVtu.filename() << " to current dir." << endl;
			} catch (const boost::filesystem::filesystem_error& ex) {
				cout << ex.what() << endl;
			}
			// TODO: if loading takes a long time could have nodes do different files,
			//       and/or only load as many as needed by DIIS
			Mesh mesh(inputMeshFiles[i]);
			PotentialField potential(&mesh,string("potential"),boundaryPotential,
					sheathPotential,fixSheathPotential);
			DensityField ionDensity(&mesh,string("ionDensity"));
			numberOfNodes = potential.entities.size();
			Eigen::VectorXd residual(numberOfNodes);
			Eigen::VectorXd diisPotential(numberOfNodes);
			for (int k=0; k<numberOfNodes; k++) {
				diisPotential[k] = potential.getField(ionDensity.entities[k]);
				// TODO: make residual calculation a function so don't risk inconsistency
				residual[k] = ( ionDensity.getField(ionDensity.entities[k]) -
						exp(diisPotential[k]) );
			}
			diisResiduals.push_back(residual);
			diisPotentials.push_back(diisPotential);
		}
#ifdef HAVE_MPI
		MPI::COMM_WORLD.Bcast(&numberOfNodes, 1, MPI::INTEGER, 0);
		Eigen::VectorXd residual(numberOfNodes);
		Eigen::VectorXd diisPotential(numberOfNodes);
		residual = Eigen::VectorXd::Zero(numberOfNodes);
		diisPotential = Eigen::VectorXd::Zero(numberOfNodes);
		if (mpiId == 0) {
			residual = diisResiduals[i];
			diisPotential = diisPotentials[i];
		}
		MPI::COMM_WORLD.Bcast(residual.data(), residual.size(), MPI::DOUBLE, 0);
		MPI::COMM_WORLD.Bcast(diisPotential.data(), diisPotential.size(), MPI::DOUBLE, 0);
		if (mpiId != 0) {
			diisResiduals.push_back(residual);
			diisPotentials.push_back(diisPotential);
		}
#endif
		lastInputIteration = inputIterations[i];
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
	PotentialField potential(&mesh,string("potential"),boundaryPotential,
			sheathPotential,fixSheathPotential);
	PotentialField previousPotential(&mesh,string("previousPotential"),
			boundaryPotential,sheathPotential,fixSheathPotential);
	PotentialField potentialScan(&mesh,string("potentialScan"),boundaryPotential,
			sheathPotential,fixSheathPotential,numberOfPotentialValues);
	// TODO: this doesn't work properly with restart/continuation
//	Field<double> potentialHistory(&mesh,string("potentialHistory"),numberOfIterations);
//	Field<double> ionDensityHistory(&mesh,string("ionDensityHistory"),numberOfIterations);
	ElectricField eField(&mesh,string("eField"),vertexType,debyeLength,doLuDecomposition);
	Field<vect3d> ionVelocity(&mesh,string("ionVelocity"),iBase_VERTEX);
	Field<double> ionTemperature(&mesh,string("ionTemperature"),iBase_VERTEX);
	Field<vect3d> ionVelocityScan(&mesh,string("ionVelocityScan"),iBase_VERTEX,numberOfPotentialValues);
	Field<double> ionTemperatureScan(&mesh,string("ionTemperatureScan"),iBase_VERTEX,numberOfPotentialValues);
//	Field<vect3d> ionVelocityNegativePerturbation(&mesh,string("NPionVelocity"),iBase_VERTEX);
//	Field<double> ionTemperatureNegativePerturbation(&mesh,string("NPionTemperature"),iBase_VERTEX);
	// TODO: updating other moments through density not very clean/transparent
	DensityField ionDensity(&mesh,string("ionDensity"),&ionVelocity,&ionTemperature);
//	DensityField previousIonDensity(&mesh,string("previousIonDensity"),&ionVelocity,&ionTemperature);
	DensityField ionDensityScan(&mesh,string("ionDensityScan"),&ionVelocityScan,&ionTemperatureScan,
			numberOfPotentialValues);
	DensityField referenceElectronDensity(&mesh,string("referenceElectronDensity"));
//	DensityField electronDensity(&mesh,string("electronDensity"));
//	DensityField density(&mesh,string("density"));
//	DensityField ionDensityPositivePerturbation(&mesh,string("PPionDensity"));
//	DensityField ionDensityNegativePerturbation(&mesh,string("NPionDensity"),
//			&ionVelocityNegativePerturbation,&ionTemperatureNegativePerturbation);
//	DensityField electronDensityPositivePerturbation(&mesh,string("PPelectronDensity"));
//	DensityField electronDensityNegativePerturbation(&mesh,string("NPelectronDensity"));
//	DensityDerivativeField ionDensityDerivative(&mesh,string("ionDensDeriv"),iBase_VERTEX);
	ShortestEdgeField shortestEdge(&mesh,string("shortestEdge"));

	SurfaceField<double,vtkDoubleArray> surfacePotential(&surfaceMesh, "surfacePotential");
	SurfaceField<double,vtkDoubleArray,3,vect3d> surfaceEField(&surfaceMesh, "surfaceEField");
	SurfaceField<double,vtkDoubleArray> ionSurfaceDensity(&surfaceMesh, "ionSurfaceDensity");
	SurfaceField<double,vtkDoubleArray,3,vect3d> ionSurfaceVelocity(&surfaceMesh, "ionSurfaceVelocity");
	SurfaceField<double,vtkDoubleArray> ionSurfaceTemperature(&surfaceMesh, "ionSurfaceTemperature");
	SurfaceField<double,vtkDoubleArray> surfaceReferenceDensity(&surfaceMesh, "surfaceReferenceDensity");

	double noPotentialPerturbation = 0.;
	double positivePotentialPerturbation = 0.15;
	double negativePotentialPerturbation = -0.15;

	if (mpiId == 0)
		cout << endl << "Calculating shortest edge of each region..." << endl;
	shortestEdge.calcField();
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
	// TODO: add detection of forgotten distribution function...or add to constructor?
	// TODO: Allow parallel electron temperature to differ from ions?
	potential.setReferenceElectronTemperature(*parallelTemperatureProfile_ptr.get());
	ionDensity.setDistributionFunction(distributionFunction);
	ionDensityScan.setDistributionFunction(distributionFunction);
//	ionDensityNegativePerturbation.setDistributionFunction(distributionFunction);

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
				objectPotential, sheathPotential);
		previousPotential.copyValues(potential);
	}
//	potentialHistory.copyValues(potential,0);

	if (!useDensityFromInput) {
		if (mpiId == 0)
			cout << endl << "Setting density..." << endl;
		ionDensity.calcField(vertexType, potential, referenceElectronDensity,
				*parallelTemperatureProfile_ptr.get(), 1.);
	}
//	ionDensityHistory.copyValues(ionDensity,0);

	if (doPoissonTest) {
		if (mpiId == 0)
			cout << endl << "Setting poissonCubeTest density..." << endl;
		ionDensity.poissonCubeTest(debyeLength);
		if (mpiId == 0)
			cout << endl << "Calculating electric field..." << endl;
		eField.calcField(&potential, vertexType, ionDensity, debyeLength);
//		potentialHistory.copyValues(potential,1);
//		ionDensityHistory.copyValues(ionDensity,1);
		if (mpiId == 0){
			int i = 0;
			stringstream iterMeshFileName;
			int periodLocation = outputFile.rfind(".");
			iterMeshFileName << outputFile.substr(0,periodLocation)
					<< setfill('0') << setw(2) << i << setfill(' ') << setw(0) <<
					outputFile.substr(periodLocation);
//			// TODO: debugging
//			cout << iterMeshFileName.str() << endl;
			mesh.save(iterMeshFileName.str());
			// mesh.save() destroys vector tags, so update
			// TODO: do this automatically?
//			eField.updateTagHandle();
//			ionVelocity.updateTagHandle();
//			ionVelocityNegativePerturbation.updateTagHandle();
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

	for (int i=lastInputIteration+1; i<numberOfIterations; i++) {
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
//		previousIonDensity.copyValues(ionDensity);
		if (useDensityFromInput && i==0) {
			if (mpiId == 0)
				cout << endl << "Using ion density from input file." << endl;
		} else {
			if (numberOfPotentialValues>1) {
				if (mpiId == 0)
					cout << endl << "Setting potential scan values..." << endl;
				// TODO: this only works because copyValues hasn't been generalized to multi-comp.
				potentialScan.copyValues(potential);
				potentialScan.computePerturbedPotentials(negativePotentialPerturbation,
						positivePotentialPerturbation);
				if (mpiId == 0)
					cout << endl << "Calculating ion density scan..." << endl;
				// TODO: use different density file
				ionDensityScan.calcField(eField, &potentialScan, referenceElectronDensity, faceType, vertexType,
						shortestEdge, 1., noPotentialPerturbation,
						densityFile);
			}
			if (mpiId == 0)
				cout << endl << "Calculating ion density..." << endl;
			// TODO: Make these inputs?
			bool doAllNodes=true;
//			int doAllNodesEveryNIterations=10;
//			if ((i%doAllNodesEveryNIterations)==0 || i==numberOfIterations-1) {
//				doAllNodes=true;
//			} else {
//				doAllNodes=false;
//			}
			double unconvergednessThreshold=0.03;
			ionDensity.calcField(eField, &potential, referenceElectronDensity, faceType, vertexType,
					shortestEdge, 1., noPotentialPerturbation,
					doAllNodes, unconvergednessThreshold, densityFile);
//			if (mpiId == 0)
//				cout << endl << "Calculating NP ion charge-density..." << endl;
//			ionDensityNegativePerturbation.calcField(eField, potential,
//					referenceElectronDensity, faceType, vertexType,
//					shortestEdge, 1., negativePotentialPerturbation,
//					densityFile);
//			previousIonDensity.copyValues(ionDensityNegativePerturbation);
		}
//		ionDensityHistory.copyValues(ionDensity,i+1);
//		if (mpiId == 0)
//			cout << endl << "Calculating PP ion charge-density..." << endl;
//		ionDensityPositivePerturbation.calcField(eField, potential,
//				faceType, vertexType,
//				shortestEdge, 1., positivePotentialPerturbation,
//				densityFile);
//		if (mpiId == 0)
//			cout << endl << "Calculating charge density..." << endl;
//		density.calcField(ionDensity, electronDensity);
		if (extern_numberOfSurfaceEvalPoints>0) {
			surfacePotential.copyFromField(potential);
			surfaceEField.copyFromField(eField);
			ionSurfaceDensity.copyFromField(ionDensity);
			ionSurfaceVelocity.copyFromField(ionVelocity);
			ionSurfaceTemperature.copyFromField(ionTemperature);
			surfaceReferenceDensity.copyFromField(referenceElectronDensity);
		}
		if (mpiId == 0)
			cout << endl << endl << endl;
		if (mpiId == 0){
			stringstream iterMeshFileName;
			int periodLocation = outputFile.rfind(".");
			iterMeshFileName << outputFile.substr(0,periodLocation)
					<< setfill('0') << setw(2) << i << setfill(' ') << setw(0) <<
					outputFile.substr(periodLocation);
			mesh.save(iterMeshFileName.str());
			// mesh.save() destroys vector tags, so update
			// TODO: do this automatically?
//			eField.updateTagHandle();
//			ionVelocity.updateTagHandle();
//			ionVelocityNegativePerturbation.updateTagHandle();

			if (extern_numberOfSurfaceEvalPoints>0) {
				stringstream iterSurfaceMeshFileName;
				periodLocation = outputSurfaceMeshFile.rfind(".");
				iterSurfaceMeshFileName << outputSurfaceMeshFile.substr(0,periodLocation)
						<< setfill('0') << setw(2) << i << setfill(' ') << setw(0) <<
						outputSurfaceMeshFile.substr(periodLocation);
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
		previousPotential.copyValues(potential);
//		previousPotential += negativePotentialPerturbation;
//		ionDensityDerivative.calcField(ionDensity,previousIonDensity,potential,previousPotential);
		if (debyeLength==0.) {
			if (mpiId == 0)
				cout << endl << "Storing quantities for DIIS..." << endl;
			int numberOfNodes=ionDensity.entities.size();
			Eigen::VectorXd residual(numberOfNodes);
			Eigen::VectorXd diisPotential(numberOfNodes);
			for (int k=0; k<numberOfNodes; k++) {
				diisPotential[k] = potential.getField(ionDensity.entities[k]);
				residual[k] = ( ionDensity.getField(ionDensity.entities[k]) -
						exp(diisPotential[k]) );
			}
			diisResiduals.push_back(residual);
			diisPotentials.push_back(diisPotential);
			// TODO: rename numberOfIterationsToAveragePotentialOver since now different
			if (i<numberOfIterationsToAveragePotentialOver+1) {
				if (mpiId == 0)
					cout << endl << "Calculating updated potential..." << endl;
//				if (i==0 && !usePotentialFromInput) {
					potential.calcField(ionDensity, ionVelocity, vertexType, potentialFile, boundaryPotential,
							sheathPotential, fixSheathPotential);
			} else {
				if (mpiId == 0)
					cout << endl << "Updating potential using DIIS..." << endl;
				int nIterDIIS = min(i+1,numberOfIterationsToAveragePotentialOver);
				Eigen::MatrixXd leastSquaresMatrix;
				Eigen::VectorXd leastSquaresRHS;
				Eigen::VectorXd leastSquaresSolution;
				leastSquaresMatrix = Eigen::MatrixXd::Zero(nIterDIIS+1, nIterDIIS+1);
				leastSquaresRHS = Eigen::VectorXd::Zero(nIterDIIS+1);
				leastSquaresSolution = Eigen::VectorXd::Zero(nIterDIIS+1);
				vector<Eigen::VectorXd> residuals;
				vector<Eigen::VectorXd> potentials;
				vector<int> residualIndexToIteration;
				for (int j=0; j<nIterDIIS; j++) {
					int correspondingIteration = i-(nIterDIIS-1)+j;
//					int correspondingIteration = i-j;
					// TODO: could use residualIndexToIteration rather than copying to vectors
					residualIndexToIteration.push_back(correspondingIteration);
					residuals.push_back(diisResiduals[correspondingIteration]);
					potentials.push_back(diisPotentials[correspondingIteration]);
				}
				for (int j=0; j<nIterDIIS; j++) {
					// TODO: could do this with Eigen (sub)matrix operations
					leastSquaresMatrix(j,nIterDIIS) = -1.;
					leastSquaresMatrix(nIterDIIS,j) = -1.;
					for (int k=0; k<nIterDIIS; k++) {
						leastSquaresMatrix(j,k) = residuals[j].dot(residuals[k]);
					}
				}
				leastSquaresRHS[nIterDIIS] = -1.;
				leastSquaresSolution = leastSquaresMatrix.inverse()*leastSquaresRHS;
				// TODO: debugging
				if (mpiId == 0) {
					cout << leastSquaresMatrix << endl << endl <<
							leastSquaresSolution << endl << endl;
//							leastSquaresRHS << endl;
				}

				Eigen::VectorXd newPotential;
				newPotential = Eigen::VectorXd::Zero(numberOfNodes);
				for (int j=0; j<nIterDIIS; j++) {
					newPotential += leastSquaresSolution[j]*(potentials[j]+residuals[j]);
				}
				for (int k=0; k<numberOfNodes; k++) {
					potential.setField(potential.entities[k],newPotential[k]);
				}
			}

//			potentialHistory.copyValues(potential,i+1);
//			} else {
//				potential.calcField(ionDensity, ionDensityDerivative, vertexType, potentialFile, boundaryPotential,
//						sheathPotential, fixSheathPotential);
//			}
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

//			if (mpiId == 0)
//				cout << endl << "Averaging potential over iterations..." << endl;
//			vector<double> weights(numberOfIterations);
//			for (int j=0; j<numberOfIterations; j++) {
//				if (j>(i+1-numberOfIterationsToAveragePotentialOver) && j<i+2) {
//					if (i+2<numberOfIterationsToAveragePotentialOver) {
//						weights[j] = 1./(i+2);
//					} else {
//						weights[j] = 1./numberOfIterationsToAveragePotentialOver;
//					}
//				} else {
//					weights[j] = 0.;
//				}
////				// TODO: debugging
////				cout << weights[j] << "...";
//			}
//			potentialHistory.computeWeightedAverage(&potential,weights);
		} else {
			if (mpiId == 0)
				cout << endl << "Calculating updated potential and electric field..." << endl;
			eField.calcField(&potential, vertexType, ionDensity, debyeLength);
//			potentialHistory.copyValues(potential,i+1);
		}
//		ionDensityDerivative.calcField(ionDensity,previousIonDensity,potential,previousPotential);
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
