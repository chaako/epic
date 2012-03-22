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

	char *options = NULL;
	int options_len = 0;
	int dim, num, ierr;
	iMesh_Instance mesh;
	iBase_EntitySetHandle root;
	iBase_EntityHandle *ents0d = NULL;
	iBase_EntityHandle *ents2d = NULL, *ents3d = NULL;
	int ents0d_alloc = 0, ents0d_size;
	int ents2d_alloc = 0, ents2d_size;
	int ents3d_alloc = 0, ents3d_size;
	iBase_TagHandle code_tag, guard_tag, potential_tag, eField_tag;
	iBase_TagHandle eFieldX_tag, eFieldY_tag, eFieldZ_tag;

	if (argc<3) {
		printf("usage: %s meshin meshout\n", argv[0]);
		exit(1);
	}

	Mesh mesh2(argv[1]);
	mesh2.printElementNumbers();

	PotentialField potential(&mesh2,std::string("potential"));
	ElectricField eField(&mesh2,std::string("eField"));
	DensityField density(&mesh2,std::string("density"));
	DensityField ionDensity(&mesh2,std::string("ionDensity"));
	DensityField electronDensity(&mesh2,std::string("electronDensity"));

	std::cout << "Calculating potential..." << std::endl;
	potential.calcField();
	std::cout << "Calculating electric field..." << std::endl;
	eField.calcField(potential);
	std::cout << "Calculating ion charge-density..." << std::endl;
	ionDensity.calcField(eField, potential,1.);
	std::cout << "Calculating electron density..." << std::endl;
	electronDensity.calcField(eField, potential,-1.);
	std::cout << "Calculating charge density..." << std::endl;
	density.calcField(ionDensity,electronDensity);
	std::cout << "Calculating updated potential..." << std::endl;
	potential.calcField(ionDensity,electronDensity);

	{
	const char *fName = "integratedOrbits.p3d";
	FILE* outFile = fopen(fName, "w");
	fprintf(outFile, "# x y z density\n");

//	iBase_EntityHandle node = mesh2.getRandomVertex();
//	IntegrandContainer integrandContainer;
//	integrandContainer.mesh_ptr = &mesh2;
//	integrandContainer.node = node;
//	integrandContainer.electricField_ptr = &eField;
//	int vdim=3;
//	double xmin[vdim], xmax[vdim];
//	for (int i=0; i<vdim; i++) {
//		xmin[i] = -1.;
//		xmax[i] = 1.;
//	}
//	double density=0.;
//	double error=0.;
//	adapt_integrate(1, &valueFromBoundary, (void*)&integrandContainer,
//			vdim, xmin, xmax, 100, 1.e-2, 1.e-2, &density, &error);
//	std::cout << "error = " << error << std::endl;
//	std::cout << "density = " << density << std::endl;

//	int neval;
//	int fail;
//	double prob;
//	Vegas(3, 1, &valueFromBoundaryCuba, (void*)&integrandContainer,
//	  1.e-8, 1.e-8, 0, 999,  100, 100, 10, 10, 10, 0, NULL,
//	  &neval, &fail, &density, &error, &prob);
//	std::cout << "error = " << error << std::endl;
//	std::cout << "density = " << density << std::endl;
//	std::cout << "neval = " << neval << std::endl;
//	std::cout << "fail = " << fail << std::endl;

	fclose(outFile);
	}

	mesh2.save(argv[2]);
	return 0;

	/* create the Mesh instance */
	iMesh_newMesh(options, &mesh, &ierr, options_len);
	CHECK("Problems instantiating interface.");
	iMesh_getRootSet(mesh, &root, &ierr);
	CHECK("Problems getting root set");

	/* load the mesh */
//	iMesh_load(mesh, root, argv[1], options, &ierr,
//			strlen(argv[1]), options_len);
	// FMDB's importVTK can't handle our tags, so use custom version
	ierr = custom_importVTK((mMesh *)mesh, argv[1]);
	CHECK("Load failed");

	/* report the number of elements of each dimension */
	for (dim = iBase_VERTEX; dim <= iBase_REGION; dim++) {
		iMesh_getNumOfType(mesh, root, dim, &num, &ierr);
		CHECK("Failure in getNumOfType");
		printf("Number of %d-dimensional elements = %d\n", dim, num);
	}

	// create potential tag
	iMesh_createTag(mesh, "potential", 1, iBase_DOUBLE, &potential_tag,
			&ierr, 9);
	CHECK("Failure creating potential tag");

	// set potential value for all 0d elements
	iMesh_getEntities(mesh, root, iBase_VERTEX, iMesh_ALL_TOPOLOGIES,
			&ents0d, &ents0d_alloc, &ents0d_size, &ierr);
	CHECK("Couldn't get vertex entities");
	for (int i = 0; i < ents0d_size; i++) {
		double potential=0., x=0., y=0., z=0.;
		iMesh_getVtxCoord(mesh, ents0d[i], &x, &y, &z, &ierr);
		CHECK("Failure getting vertex coordinates");
		potential = -1./sqrt(x*x+y*y+z*z);
		iMesh_setDblData(mesh, ents0d[i], potential_tag, potential,
				&ierr);
		CHECK("Failure setting potential tag");
	}

	// create eField tag
	iMesh_createTag(mesh, "eField", sizeof(Eigen::Vector3d), iBase_BYTES,
			&eField_tag, &ierr, 6);
	CHECK("Failure creating eField tag");
	iMesh_createTag(mesh, "eFieldX", 1, iBase_DOUBLE,
			&eFieldX_tag, &ierr, 7);
	CHECK("Failure creating eField tag");
	iMesh_createTag(mesh, "eFieldY", 1, iBase_DOUBLE,
			&eFieldY_tag, &ierr, 7);
	CHECK("Failure creating eField tag");
	iMesh_createTag(mesh, "eFieldZ", 1, iBase_DOUBLE,
			&eFieldZ_tag, &ierr, 7);
	CHECK("Failure creating eField tag");

	// Calculate electric field for all 0d elements
	for (int i = 0; i < ents0d_size; i++) {
		vector<iBase_EntityHandle> superCellFaces =
				getSuperCellFaces(mesh, ents0d[i]);
		Eigen::Vector3d eField(0,0,0);

		double x=0,y=0,z=0;
		iMesh_getVtxCoord(mesh, ents0d[i], &x, &y, &z, &ierr);
		CHECK("Failure getting vertex coordinates");
		Eigen::Vector3d point(x,y,z);
		double volume=0;

		for (int j=0; j<superCellFaces.size(); j++)  {
			Eigen::Vector3d surfaceVector =
					getSurfaceVector(mesh, point, superCellFaces[j]);
			double potential = getAverageDblData(mesh, superCellFaces[j],
					potential_tag);
			volume += getTetVolume(mesh, point, superCellFaces[j]);
			eField -= potential*surfaceVector;
		}
		eField /= volume;

		iMesh_setData(mesh, ents0d[i], eField_tag, &eField,
				sizeof(Eigen::Vector3d), &ierr);
		CHECK("Failure setting eField tag");
		iMesh_setDblData(mesh, ents0d[i], eFieldX_tag, (double)eField[0], &ierr);
		CHECK("Failure setting eFieldX tag");
		iMesh_setDblData(mesh, ents0d[i], eFieldY_tag, (double)eField[1], &ierr);
		CHECK("Failure setting eFieldY tag");
		iMesh_setDblData(mesh, ents0d[i], eFieldZ_tag, (double)eField[2], &ierr);
		CHECK("Failure setting eFieldZ tag");
	}

	// start an orbit at a random node
	{
	iBase_TagHandle visited_tag;
	iMesh_createTag(mesh, "visited", 1, iBase_INTEGER,
			&visited_tag, &ierr, 7);
	CHECK("Failure creating visited tag");
	// set default guard_cells value for all 3d elements
	iMesh_getEntities(mesh, root, iBase_REGION, iMesh_ALL_TOPOLOGIES,
			&ents3d, &ents3d_alloc, &ents3d_size, &ierr);
	CHECK("Couldn't get region entities");
	for (int i = 0; i < ents3d_size; i++) {
		iMesh_setIntData(mesh, ents3d[i], visited_tag, 0, &ierr);
		CHECK("Failure setting default visited tag");
	}
	if (ents3d) free(ents3d);
	ents3d_alloc = 0;
	srand(999);
	const char *fName = "integratedOrbitTest.p3d";
	FILE* outFile = fopen(fName, "w");
	fprintf(outFile, "# x y z density\n");
	int nOrbits=0;
	for (double multiplier=0.5; multiplier<3; multiplier*=2.) {
		nOrbits++;
		int iSelectedNode = rand() % ents0d_size;
		double x, y, z;
		double eFieldX, eFieldY, eFieldZ;
		iMesh_getVtxCoord(mesh, ents0d[iSelectedNode], &x, &y, &z, &ierr);
		CHECK("Failure getting vertex coordinates");
		iBase_EntityHandle currentTet=NULL;
		Eigen::Vector3d currentPosition(x,y,z);
		Eigen::Vector3d currentVelocity, zHat(0.,0.,1.);
		double potential;
		iMesh_getDblData(mesh, ents0d[iSelectedNode], potential_tag, &potential,
				&ierr);
		CHECK("Failure getting potential tag");
		currentVelocity = currentPosition.cross(zHat);
		currentVelocity /= currentVelocity.norm();
		currentVelocity *= sqrt(-potential);
		currentVelocity *= multiplier;
		Eigen::Vector3d eField(0.,0.,0.);
		int eField_alloc = sizeof(Eigen::Vector3d);
		int eField_size = sizeof(Eigen::Vector3d);
		vector<Eigen::Vector3d> eFields, vertexVectors;
		// TODO: initialize these in some better way
		int nVertices=4;
		eFields.reserve(nVertices);
		vertexVectors.reserve(nVertices);
		eFields.push_back(eField);
		eFields.push_back(eField);
		eFields.push_back(eField);
		eFields.push_back(eField);
		vertexVectors.push_back(eField);
		vertexVectors.push_back(eField);
		vertexVectors.push_back(eField);
		vertexVectors.push_back(eField);
		// make output a text file in Point3D format for VisIt
		double dt=0.01, tMax=100;
		currentPosition -= currentVelocity*dt/2.;
		bool inNewTet = true;
		int nSteps=0, nNewTet=0;
		std::cout << "Initial radius=" << currentPosition.norm() << endl;
		for (double t=0; t<tMax; t+=dt) {
			nSteps++;
			currentPosition += dt*currentVelocity;
			inNewTet = !checkIfInTet(currentPosition, vertexVectors);
			if (inNewTet) {
				nNewTet++;
				// determine which tet currentPosition is in
				iBase_EntityHandle *entities = NULL;
				int entities_alloc = 0, entities_size;

				if (currentTet) {
					iMesh_getEnt2ndAdj(mesh, currentTet, iBase_VERTEX,
							iBase_REGION, &entities, &entities_alloc,
							&entities_size, &ierr);
					CHECK("Getting regions adjacent to entity failed");
				} else {
					iMesh_getEntAdj(mesh, ents0d[iSelectedNode],
							iBase_REGION, &entities, &entities_alloc,
							&entities_size, &ierr);
					CHECK("Getting regions adjacent to entity failed");
				}
				for (int i=0; i<entities_size; i++) {
					if (checkIfInTet(currentPosition, mesh, entities[i])) {
						currentTet = entities[i];
						inNewTet = false;
						break;
					}
				}
//				assert(inNewTet == false);
				if(entities) free (entities);
				entities_alloc = 0;

				if (inNewTet==true) {
					std::cout << "Failed to identify current tet: currentPosition ="
							<< currentPosition << endl;
					std::cout << "nSteps=" << nSteps << ", nNewTet=" << nNewTet << endl;
					break;
				}

				iMesh_setIntData(mesh, currentTet, visited_tag, nOrbits, &ierr);
				CHECK("Failure setting visited tag");

				// get coordinates and field at vertices
				// TODO: unify with getVertexVectors
				iBase_EntityHandle *vertices = NULL;
				int vertices_alloc = 0, vertices_size;

				iMesh_getEntAdj(mesh, currentTet, iBase_VERTEX,  &vertices, &vertices_alloc,
						&vertices_size, &ierr);
				CHECK("Getting vertices adjacent to entity failed");

				assert(vertices_size==vertexVectors.size() && vertices_size==eFields.size());
				for (int i=0; i<vertices_size; i++) {
					iMesh_getVtxCoord(mesh, vertices[i], &x, &y, &z, &ierr);
					CHECK("Failure getting vertex coordinates");
					vertexVectors[i] << x, y, z;

					iMesh_getDblData(mesh, vertices[i], eFieldX_tag, &eFieldX, &ierr);
					CHECK("Failure getting double data");
					iMesh_getDblData(mesh, vertices[i], eFieldY_tag, &eFieldY, &ierr);
					CHECK("Failure getting double data");
					iMesh_getDblData(mesh, vertices[i], eFieldZ_tag, &eFieldZ, &ierr);
					CHECK("Failure getting double data");
					eField << eFieldX, eFieldY, eFieldZ;

//					iMesh_getData(mesh, vertices[i], eField_tag, &eField,
//							&eField_alloc, &eField_size, &ierr);
//					CHECK("Failure getting eField tag");
					eFields[i] = eField;
				}
				if(vertices) free (vertices);
				vertices_alloc = 0;
			}

			assert(vertexVectors.size()==nVertices);
			std::vector<double> vertexWeights = getVertexWeights(currentPosition,
					vertexVectors);
			Eigen::Vector3d currentAcceleration(0.,0.,0.);
			assert(eFields.size()==vertexWeights.size());
			for (int i=0; i<vertexWeights.size(); i++) {
				currentAcceleration += eFields[i]*vertexWeights[i];
			}
			currentAcceleration = -currentPosition/pow(currentPosition.norm(),3.);
			double eFieldR = currentAcceleration.dot(currentPosition)/
					currentPosition.norm();
			currentVelocity += dt*currentAcceleration;
			fprintf(outFile, "%f %f %f %f\n", currentPosition[0], currentPosition[1],
					currentPosition[2], eFieldR);
//			fprintf(outFile, "%f %f %f %d\n", currentPosition[0], currentPosition[1],
//					currentPosition[2], nNewTet);
		}
		std::cout << "Final radius=" << currentPosition.norm() << endl;
	}
	fclose(outFile);
	}

	if (ents0d) free(ents0d);
	ents0d_alloc = 0;
	// destroy eField tag since VisIt doesn't understand
	iMesh_destroyTag (mesh, eField_tag, 1, &ierr);
	CHECK("Failed to destroy eField tag");

	// create guard_cells tag
	iMesh_createTag(mesh, "guard_cells", 1, iBase_INTEGER, &guard_tag,
			&ierr, 11);
	CHECK("Failure creating guard_cells tag");

	// set default guard_cells value for all 3d elements
	iMesh_getEntities(mesh, root, iBase_REGION, iMesh_ALL_TOPOLOGIES,
			&ents3d, &ents3d_alloc, &ents3d_size, &ierr);
	CHECK("Couldn't get region entities");
	for (int i = 0; i < ents3d_size; i++) {
		iMesh_setIntData(mesh, ents3d[i], guard_tag, 0, &ierr);
		CHECK("Failure setting default guard_cells tag");
	}
	if (ents3d) free(ents3d);
	ents3d_alloc = 0;

	// find and tag guard cells for surface faces
	iMesh_getEntities(mesh, root, iBase_FACE, iMesh_ALL_TOPOLOGIES,
			&ents2d, &ents2d_alloc, &ents2d_size, &ierr);
	CHECK("Couldn't get face entities");
	iMesh_getTagHandle(mesh, "cell_code", &code_tag, &ierr, 9);
	CHECK("Failure getting cell_code handle");
	for (int i = 0; i < ents2d_size; i++) {
		int cell_code = 0;
		iMesh_getIntData(mesh, ents2d[i], code_tag, &cell_code, &ierr);
		CHECK("Failure getting cell_code value");
		if (cell_code>1) {
			iBase_EntityHandle *entsGuard = NULL;
			int entsGuard_alloc = 0, entsGuard_size;
			// TODO: perhaps remove redundant guard cell returns and
			//       improve performance by using array version here
			iMesh_getEnt2ndAdj(mesh, ents2d[i], iBase_VERTEX,
					iBase_REGION, &entsGuard, &entsGuard_alloc,
					&entsGuard_size, &ierr);
			CHECK("Failure in getEnt2ndAdj");
			for (int j = 0; j < entsGuard_size; j++) {
				// TODO: what about cells adjacent to two boundaries?
				iMesh_setIntData(mesh, entsGuard[j], guard_tag,
						cell_code, &ierr);
				CHECK("Failure setting guard_cells tag");
			}
			if (entsGuard) free(entsGuard);
			entsGuard_alloc = 0;
		}
	}
	if (ents2d) free(ents2d);
	ents2d_alloc = 0;

	/* save the mesh */
	iMesh_save(mesh, root, argv[2], options, &ierr,
			strlen(argv[2]), options_len);
	CHECK("Save failed");

	iMesh_dtor(mesh, &ierr);
	CHECK("Failed to destroy interface");

	// create a new mesh for a test orbit
	iMesh_newMesh(options, &mesh, &ierr, options_len);
	CHECK("Problems instantiating interface.");

	// TODO: setting the geometric dimension appears not to work,
	//       possibly because there are no 3d elements
	dim = 3;
	iMesh_setGeometricDimension (mesh, dim, &ierr);
	CHECK("Error setting geometric dimension");
	dim = 0;
	iMesh_getGeometricDimension (mesh, &dim, &ierr);
	CHECK("Error getting geometric dimension");
	printf("dim = %d\n",dim);

	iMesh_getRootSet(mesh, &root, &ierr);
	CHECK("Problems getting root set");

	// add vertices for orbit
	iBase_EntityHandle vertices[4];
	int iVertex=0;
	// make output a text file in Point3D format for VisIt
	const char *fName = "orbitTest.p3d";
	FILE* outFile = fopen(fName, "w");
	fprintf(outFile, "# x y z density\n");
	for (double ang=0; ang<2*M_PI; ang+=0.01) {
		iBase_EntityHandle newVertex, newEdge, newRegion;
		int creationStatus;
		double x,y,z;
		x = 2.*cos(ang);
		y = 2.*sin(ang);
		z = sqrt(3.*3.-x*x-0.5*y*y);
		fprintf(outFile, "%f %f %f %f\n", x, y, z, 1.);
		iMesh_createVtx(mesh, x, y, z, &newVertex, &ierr);
		CHECK("Failure creating new vertex");
		iMesh_setVtxCoord(mesh, newVertex, x, y, z, &ierr);
		CHECK("Failure setting vertex coordinates");
		// TODO: because the geometric dimension isn't correct,
		//       getVtxCoord doesn't return the z-coordinate,
		//       but it is stored with set (though not create)
		iMesh_getVtxCoord(mesh, newVertex, &x, &y, &z, &ierr);
		CHECK("Failure getting vertex coordinates");
		vertices[2] = vertices[3];
		vertices[3] = vertices[0];
		vertices[0] = vertices[1];
		vertices[1] = newVertex;
		if (iVertex>0) {
			// create edges between nodes
			iMesh_createEnt(mesh, iMesh_LINE_SEGMENT, vertices, 2,
					&newEdge, &creationStatus, &ierr);
			CHECK("Failure creating edge");
		}
		if (iVertex==3) {
			// add a tet to change geometric dim. to 3
			iMesh_createEnt(mesh, iMesh_TETRAHEDRON, vertices, 4,
					&newRegion, &creationStatus, &ierr);
			CHECK("Failure creating region");
		}
		iVertex++;
	}
	fclose(outFile);
	iMesh_getGeometricDimension (mesh, &dim, &ierr);
	CHECK("Error getting geometric dimension");
	printf("dim = %d\n",dim);

	// save the test orbit mesh
	iMesh_save(mesh, root, "orbitTest.sms", options, &ierr,
			13, options_len);
	CHECK("Save failed");
	iMesh_save(mesh, root, "orbitTest.vtk", options, &ierr,
			13, options_len);
	CHECK("Save failed");
	iMesh_save(mesh, root, "orbitTest.nc", options, &ierr,
			12, options_len);
	CHECK("Save failed");

	iMesh_dtor(mesh, &ierr);
	CHECK("Failed to destroy interface");

	return 0;
}

Eigen::Vector3d getSurfaceVector(iMesh_Instance mesh, Eigen::Vector3d point,
		iBase_EntityHandle face) {
	int ierr;
	std::vector<Eigen::Vector3d> vertexVectors, edgeVectors;
	Eigen::Vector3d surfaceVector, referenceVector;

	iBase_TagHandle code_tag;
	int cell_code = 0;

	iMesh_getTagHandle(mesh, "cell_code", &code_tag, &ierr, 9);
	CHECK("Failure getting cell_code handle");

	iMesh_getIntData(mesh, face, code_tag, &cell_code, &ierr);
	CHECK("Failure getting cell_code value");

	vertexVectors = getVertexVectors(mesh, face);

	assert(3 == vertexVectors.size());
	edgeVectors.push_back(vertexVectors[1]-vertexVectors[0]);
	edgeVectors.push_back(vertexVectors[2]-vertexVectors[0]);

	surfaceVector = edgeVectors[0].cross(edgeVectors[1])/2.;
	referenceVector = vertexVectors[0] - point;

	// TODO: need a way to check orientation for boundary faces
	if (cell_code<=1 && referenceVector.dot(surfaceVector)<0)
		surfaceVector *= -1;
	return surfaceVector;
}

vector<iBase_EntityHandle> getSuperCellFaces(iMesh_Instance mesh,
		iBase_EntityHandle vertex) {
	int ierr;
	iBase_EntityHandle *faces = NULL;
	int faces_alloc = 0, faces_size;
	std::vector<iBase_EntityHandle> superCellFaces;
	iBase_TagHandle code_tag;

	iMesh_getTagHandle(mesh, "cell_code", &code_tag, &ierr, 9);
	CHECK("Failure getting cell_code handle");

	iMesh_getEnt2ndAdj(mesh, vertex, iBase_REGION,
			iBase_FACE, &faces, &faces_alloc,
			&faces_size, &ierr);
	CHECK("Failure in getEnt2ndAdj");
	for (int i = 0; i < faces_size; i++) {
		iBase_EntityHandle *vertices = NULL;
		int vertices_alloc = 0, vertices_size;
		bool onSuperCellSurface = true;
		int cell_code = 0;

		iMesh_getIntData(mesh, faces[i], code_tag, &cell_code, &ierr);
		CHECK("Failure getting cell_code value");

		iMesh_getEntAdj(mesh, faces[i], iBase_VERTEX,  &vertices, &vertices_alloc,
				&vertices_size, &ierr);
		CHECK("Getting vertices adjacent to face failed");

		for (int j=0; j<vertices_size; j++) {
			if (cell_code<=1 && vertex==vertices[j])
				onSuperCellSurface = false;
		}

		if (onSuperCellSurface) {
			superCellFaces.push_back(faces[i]);
		}

		if(vertices) free (vertices);
		vertices_alloc = 0;
	}
	if (faces) free(faces);
	faces_alloc = 0;

	return superCellFaces;
}

double getAverageDblData(iMesh_Instance mesh, iBase_EntityHandle entity,
		iBase_TagHandle dblData_tag) {
	// TODO: write this as a template?
	int ierr;
	iBase_EntityHandle *vertices = NULL;
	int vertices_alloc = 0, vertices_size;
	double averageDblData=0.;

	iMesh_getEntAdj(mesh, entity, iBase_VERTEX,  &vertices, &vertices_alloc,
			&vertices_size, &ierr);
	CHECK("Getting vertices adjacent to entity failed");

	for (int i=0; i<vertices_size; i++) {
		double dblData = 0;

		iMesh_getDblData(mesh, vertices[i], dblData_tag, &dblData, &ierr);
		CHECK("Failure getting double data");

		averageDblData += dblData;
	}
	averageDblData /= (double)vertices_size;

	if(vertices) free (vertices);
	vertices_alloc = 0;

	return averageDblData;
}

double getTetVolume(iMesh_Instance mesh, Eigen::Vector3d point,
		iBase_EntityHandle face) {
	std::vector<Eigen::Vector3d> vertexVectors = getVertexVectors(mesh, face);
	vertexVectors.push_back(point);

	return getTetVolume(vertexVectors);
}

double getTetVolume(std::vector<Eigen::Vector3d> vertexVectors) {
	int nEdges=3;
	std::vector<Eigen::Vector3d> edgeVectors(nEdges);

	int nVertices=4;
	assert(vertexVectors.size() == nVertices);

	edgeVectors[0] = vertexVectors[1]-vertexVectors[0];
	edgeVectors[1] = vertexVectors[2]-vertexVectors[0];
	edgeVectors[2] = vertexVectors[3]-vertexVectors[0];

	return fabs(
			edgeVectors[2].dot( edgeVectors[0].cross(edgeVectors[1]) )
			)/6.;
}

bool checkIfInTet(Eigen::Vector3d currentPosition,
		std::vector<Eigen::Vector3d> vertexVectors) {
	double tetVolume = getTetVolume(vertexVectors);
	if (tetVolume<VOLUME_TOLERANCE)
		return false;
	std::vector<double> subVolumes = getTetSubVolumes(currentPosition,
			vertexVectors);
	double sumSubVolumes =
			std::accumulate(subVolumes.begin(),subVolumes.end(),0.);
	return (fabs(sumSubVolumes-tetVolume)<VOLUME_TOLERANCE);
}

std::vector<double> getTetSubVolumes(Eigen::Vector3d point,
		std::vector<Eigen::Vector3d> vertexVectors) {
	int nVertices=4;
	assert(vertexVectors.size()==nVertices);
	std::vector<double> subVolumes(nVertices);
	for (int i=0; i<vertexVectors.size(); i++) {
		Eigen::Vector3d tmpVertex = vertexVectors[i];
		vertexVectors[i] = point;
		double volume = getTetVolume(vertexVectors);
		subVolumes[i] = volume;
		vertexVectors[i] = tmpVertex;
	}
	return subVolumes;
}

bool checkIfInTet(Eigen::Vector3d currentPosition, iMesh_Instance mesh,
		iBase_EntityHandle element) {
	std::vector<Eigen::Vector3d> vertexVectors = getVertexVectors(mesh,
			element);

	return checkIfInTet(currentPosition, vertexVectors);
}

std::vector<Eigen::Vector3d> getVertexVectors(iMesh_Instance mesh,
		iBase_EntityHandle entity) {
	int ierr;
	iBase_EntityHandle *vertices = NULL;
	int vertices_alloc = 0, vertices_size;
	std::vector<Eigen::Vector3d> vertexVectors;

	iMesh_getEntAdj(mesh, entity, iBase_VERTEX, &vertices, &vertices_alloc,
			&vertices_size, &ierr);
	CHECK("Getting vertices adjacent to entity failed");

	vertexVectors.reserve(vertices_size);
	for (int i=0; i<vertices_size; i++) {
		double x, y, z;
		iMesh_getVtxCoord(mesh, vertices[i], &x, &y, &z, &ierr);
		CHECK("Failure getting vertex coordinates");
		Eigen::Vector3d vertexVector(x,y,z);
		vertexVectors.push_back(vertexVector);
	}
	if(vertices) free (vertices);
	vertices_alloc = 0;

	return vertexVectors;
}

std::vector<double> getVertexWeights(Eigen::Vector3d point,
		std::vector<Eigen::Vector3d> vertexVectors) {
	std::vector<double> subVolumes = getTetSubVolumes(point,
			vertexVectors);
	double totalVolume = getTetVolume(vertexVectors);
	for(int i=0; i<subVolumes.size(); i++)
	    subVolumes[i] /= totalVolume;

	return subVolumes;
}


int valueFromBoundaryCuba(const int *ndim, const double x[],
  const int *ncomp, double f[], void *integrandContainer_ptr) {
	assert(*ndim==3);
	assert(*ncomp==1);
	double y[3];
	for (int i=0; i<*ndim; i++) {
		y[i] = 2.*x[i]-1;
	}
	valueFromBoundary((unsigned)*ndim, y, integrandContainer_ptr, *ncomp, f);
	*f *= pow(2.,(double)*ndim);

	return 0;
}

void valueFromBoundary(unsigned ndim, const double *x,
		void *integrandContainer_ptr, unsigned fdim, double *fval) {
	assert(ndim==3);
	assert(fdim==1);
	Eigen::Vector3d velocity;
//	for (int i=0; i<ndim; i++) {
//		velocity[i] = x[i]/(1.-pow(x[i],2.));
//	}
	double v = sqrt(-4.*log((1.+x[0])/2.));
	double theta = acos(x[1]);
	double phi = (1.+x[2])*M_PI;
	velocity[0] = v*cos(phi)*sin(theta);
	velocity[1] = v*sin(phi)*sin(theta);
	velocity[2] = v*cos(theta);
	Mesh *mesh_ptr = ((IntegrandContainer*)integrandContainer_ptr)->mesh_ptr;
	iBase_EntityHandle node =
			((IntegrandContainer*)integrandContainer_ptr)->node;
	ElectricField *electricField_ptr =
			((IntegrandContainer*)integrandContainer_ptr)->electricField_ptr;
	PotentialField *potentialField_ptr =
			((IntegrandContainer*)integrandContainer_ptr)->potentialField_ptr;
	double charge = ((IntegrandContainer*)integrandContainer_ptr)->charge;
	Orbit orbit(mesh_ptr,node,velocity,charge);
	orbit.integrate(*electricField_ptr, *potentialField_ptr);
	*fval = 0.;
	// TODO: shouldn't hard-code domain here
	if (orbit.finalPosition.norm()>4.9 && !orbit.negativeEnergy) {
		*fval = 1./pow(2.*M_PI,3./2.);
		*fval *= exp(-pow(orbit.finalVelocity.norm(),2.)/2.);
		*fval /= exp(-pow(v,2.)/2.);
		*fval *= M_PI;
		*fval *= (1.+x[0])*sqrt(-log((1.+x[0])/2.));
//		for (int i=0; i<ndim; i++) {
//			double t = x[i];
//			*fval *= (1+pow(t,2.))/pow(1.-pow(t,2.),2.);
//		}

	} else if (orbit.finalPosition.norm()>1.) {
//		std::cout << "orbit terminated with final position in domain:" <<
//				std::endl << orbit.finalPosition << std::endl;
	}
	if (((IntegrandContainer*)integrandContainer_ptr)->outFile) {
		fprintf(((IntegrandContainer*)integrandContainer_ptr)->outFile,
				"%f %f %f %f\n", x[0], x[1], x[2], *fval);
	}
}

clock_t extern_findTet=0, extern_checkIfInNewTet=0;
