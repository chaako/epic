#include "functions.h"
#include "classes.h"

bool vect3dLessThan(vect3d a, vect3d b) {
	bool aLessThanB = false;
	// TODO: sort differently than by components?
	// TODO: do components with recursive function?
	if (a[0] < b[0]) {
		aLessThanB = true;
	} else if (a[0] == b[0]) {
		if (a[1] < b[1]) {
			aLessThanB = true;
		} else if (a[1] == b[1]) {
			if (a[2] < b[2]) {
				aLessThanB = true;
			}
		}
	}
	return aLessThanB;
}

void getIterationDataFromFile(boost::filesystem::path& inputPath,
		string &inputMeshFile,
		Eigen::VectorXd *diisPotential, Eigen::VectorXd *residual) {
	int numberOfNodes=diisPotential->size();
	boost::filesystem::path cwd(".");
	boost::filesystem::path cwdPath = boost::filesystem::system_complete(cwd);
	// TODO: better way to strip off "."/get full path?
	try {
		boost::filesystem::path copyTarget(cwdPath.parent_path().string()+
				"/"+inputPath.filename().string());
		boost::filesystem::copy_file(inputMeshFile,copyTarget);
		boost::filesystem::path copyTargetVtu(cwdPath.parent_path().string()+
				"/"+inputPath.stem().string() + ".vtu");
//			boost::filesystem::path vtuExtension(".vtu");
		boost::filesystem::path inputFileVtu(inputMeshFile);
		inputFileVtu.replace_extension(".vtu");
		boost::filesystem::copy_file(inputFileVtu,copyTargetVtu);
		cout << "Copied " << copyTarget.filename() << " and " <<
				copyTargetVtu.filename() << " to current dir." << endl;
	} catch (const boost::filesystem::filesystem_error& ex) {
		cout << ex.what() << endl;
	}
	// TODO: if loading takes a long time could have nodes do different files,
	//       and/or only load as many as needed by DIIS
	Mesh mesh(inputMeshFile, false);
	PotentialField potential(&mesh,string("potential"));
	DensityField ionDensity(&mesh,string("ionDensity"));
	for (int k=0; k<numberOfNodes; k++) {
		diisPotential->operator[](k) = potential.getField(ionDensity.entities[k]);
		// TODO: make residual calculation a function so don't risk inconsistency
		residual->operator[](k) = ( ionDensity.getField(ionDensity.entities[k]) -
				exp(diisPotential->operator[](k)) );
	}
}

#ifdef HAVE_MPI
void requestIterationDataFromSlaves(int numberOfFiles, int numberOfNodes,
		vector<Eigen::VectorXd> *diisPotentials, vector<Eigen::VectorXd> *diisResiduals) {
	int nProcesses = MPI::COMM_WORLD.Get_size();
	MPI::Status status;
	int fileCounter=0;
	int fileNumber=-1;

	// Send one file to each process (some may not get one)
	for (int rank=1; rank<min(nProcesses,int(numberOfFiles+1)); ++rank) {
		if (fileCounter<numberOfFiles) {
			fileNumber = fileCounter;
			fileCounter++;
			MPI::COMM_WORLD.Send(&fileNumber, 1, MPI::INT, rank, WORKTAG);
		} else {
			throw string("fileCounter out of bounds in requestIterationDataFromSlaves");
		}
	}

//	Eigen::VectorXd residual(numberOfNodes);
//	Eigen::VectorXd diisPotential(numberOfNodes);
//	for (int i=0; i<numberOfFiles; i++) {
//		diisPotentials->push_back(diisPotential);
//		diisResiduals->push_back(residual);
//	}

	// Process incoming density and send new nodes until all done
	while (fileCounter<numberOfFiles) {
		fileNumber = fileCounter;
		fileCounter++;
		status = receiveIterationData(numberOfNodes, diisPotentials, diisResiduals);
		MPI::COMM_WORLD.Send(&fileNumber, 1, MPI::INT, status.Get_source(),
				WORKTAG);
	}

	// Process any outstanding iterations
	for (int rank=1; rank<min(nProcesses,int(numberOfFiles+1)); ++rank) {
		status = receiveIterationData(numberOfNodes, diisPotentials, diisResiduals);
	}

	// Send empty message with DIETAG to signal done with nodes
	for (int rank=1; rank<nProcesses; ++rank) {
		MPI::COMM_WORLD.Send(0, 0, MPI::INT, rank, DIETAG);
	}
}
#endif

#ifdef HAVE_MPI
MPI::Status receiveIterationData(int numberOfNodes,
		vector<Eigen::VectorXd> *diisPotentials, vector<Eigen::VectorXd> *diisResiduals) {
	MPI::Status status;
	Eigen::VectorXd residual(numberOfNodes);
	Eigen::VectorXd diisPotential(numberOfNodes);
	MPI::COMM_WORLD.Recv(diisPotential.data(), diisPotential.size(), MPI::DOUBLE, MPI_ANY_SOURCE,
			MPI_ANY_TAG, status);
	int fileNumber = status.Get_tag();
	int source = status.Get_source();
	MPI::COMM_WORLD.Recv(residual.data(), residual.size(), MPI::DOUBLE, source, fileNumber, status);
	if (fileNumber>=0 && fileNumber<diisPotentials->size()) {
		diisPotentials->operator[](fileNumber) = diisPotential;
		diisResiduals->operator[](fileNumber) = residual;
	}
	return status;
}
#endif

#ifdef HAVE_MPI
void processIterationDataRequests(int numberOfNodes, vector<boost::filesystem::path>& inputPaths,
		vector<string>& inputMeshFiles) {
	MPI::Status status;
	int fileNumber;
	Eigen::VectorXd residual(numberOfNodes);
	Eigen::VectorXd diisPotential(numberOfNodes);

	while (1) {
		MPI::COMM_WORLD.Recv(&fileNumber, 1, MPI::INT, 0, MPI_ANY_TAG,
				status);
		if (status.Get_tag() == DIETAG) {
			return;
		}
		getIterationDataFromFile(inputPaths[fileNumber], inputMeshFiles[fileNumber],
				&diisPotential, &residual);
		MPI::COMM_WORLD.Send(diisPotential.data(), diisPotential.size(), MPI::DOUBLE, 0, fileNumber);
		MPI::COMM_WORLD.Send(residual.data(), residual.size(), MPI::DOUBLE, 0, fileNumber);
	}
}
#endif
