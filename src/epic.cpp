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

// TODO: find better way to distinguish orbits in output
int extern_orbitNumber=0;

int setBSizeField(pMesh mesh, pSField field, void *)
{
  pVertex vt;
  double h[3], dirs[3][3], xyz[3], norm;
  VIter vit=M_vertexIter(mesh);
  while( vt=VIter_next(vit) ) {
    h[0] = 1.0;
    h[1] = .2;
    h[2] = 2.5;

    dirs[0][0]=1.;
    dirs[0][1]=0;
    dirs[0][2]=0;
    dirs[1][0]=0;
    dirs[1][1]=1.;
    dirs[1][2]=0;
    dirs[2][0]=0;
    dirs[2][1]=0;
    dirs[2][2]=1.;

    ((PWLsfield *)field)->setSize((pEntity)vt,dirs,h);
  }
  VIter_delete (vit);
  return 1;
}

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

	if (argc<3) {
		printf("usage: %s meshin meshout\n", argv[0]);
		exit(1);
	}

//	// TODO: remove this testing ground
//	{
////		int playground_meshAdapt(int argc, char* argv[]);
////		playground_meshAdapt(argc, argv);
//		int playground_netgen(int argc, char* argv[]);
//		playground_netgen(argc, argv);
//		exit(0);
//	}

	string inputMeshFile(argv[1]);
	stringstream volumeMeshFile;
	// TODO: distinguish between surface and volume mesh in less obscure way than format
	if (inputMeshFile.find(".vtu")!=string::npos) {
		SurfaceMesh surfaceMesh(inputMeshFile);
		double rotationAngle = M_PI/5;
		// TODO: don't hard code collector cell_code
		surfaceMesh.rotateSurface(4, vect3d(0.,0.,0.), vect3d::UnitY(), rotationAngle);
		int periodLocation = inputMeshFile.rfind(".vtu");
		stringstream rotatedSurfaceMeshFile;
		rotatedSurfaceMeshFile << inputMeshFile.substr(0,periodLocation)
						<< "_rotated_" << rotationAngle << ".vtu";
		surfaceMesh.save(rotatedSurfaceMeshFile.str());
		surfaceMesh.createVolumeMesh();
		volumeMeshFile << inputMeshFile.substr(0,periodLocation)
						<< "_meshed_" << rotationAngle << ".vtk";
		surfaceMesh.saveVolumeMesh(volumeMeshFile.str());
	}
	Mesh mesh(volumeMeshFile.str());
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

	mesh.classifyBoundariesForMeshRefinement(faceType);
	pMesh part;
	FMDB_Mesh_GetPart((mMesh*)mesh.meshInstance, 0, part);
	pSField field=new PWLsfield(part);
	meshAdapt rdr(part,field,0,0);  // snap off; do refinement only
//	rdr.run(2,1, setBSizeField);
//
//	char mesh_file[256];
//	char outmesh[256];
//	char without_extension[256];
//
//	snprintf(without_extension,strlen(argv[1])-3,"%s",argv[1]);
//	sprintf(mesh_file,"%s",argv[1]);
//	sprintf(outmesh,"%s-refined.sms",without_extension);
//	FMDB_Mesh_WriteToFile(part->getMesh(), outmesh, 0);
//
//	{
//		using namespace nglib;
//		Ng_Init();
//		Ng_Exit();
//	}
//	exit(0);

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
//		if (mpiId == 0)
//			cout << endl << "Saving current potential..." << endl;
//		stringstream potentialCopyName;
//		potentialCopyName << "potIter" << setfill('0') << setw(2) << i;
//		PotentialField potentialCopy(potential,potentialCopyName.str());
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
		// TODO: can't use mesh.save() here because of eField destruction
		if (mpiId == 0){
			stringstream iterMeshFileName;
			string outFile(argv[2]);
			int periodLocation = outFile.rfind(".");
			iterMeshFileName << outFile.substr(0,periodLocation)
					<< setfill('0') << setw(2) << i << outFile.substr(periodLocation);
			mesh.save(iterMeshFileName.str());
//			char *options = NULL;
//			int options_len = 0;
//			int ierr;
//			iMesh_save(mesh.meshInstance, mesh.rootEntitySet, iterMeshFileName.str().c_str(),
//					options, &ierr, iterMeshFileName.str().length(), options_len);
//			CHECK("Save failed");
		}
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
