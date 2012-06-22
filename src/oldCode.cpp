{
	char *options = NULL;
	int options_len = 0;
	int dim, num, ierr;
	iMesh_Instance mesh;
	iBase_EntitySetHandle root;
	entHandle *ents0d = NULL;
	entHandle *ents2d = NULL, *ents3d = NULL;
	int ents0d_alloc = 0, ents0d_size;
	int ents2d_alloc = 0, ents2d_size;
	int ents3d_alloc = 0, ents3d_size;
	iBase_TagHandle code_tag, guard_tag, potential_tag, eField_tag;
	iBase_TagHandle eFieldX_tag, eFieldY_tag, eFieldZ_tag;

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
	entHandle vertices[4];
	int iVertex=0;
	// make output a text file in Point3D format for VisIt
	const char *fName = "orbitTest.p3d";
	FILE* outFile = fopen(fName, "w");
	fprintf(outFile, "# x y z density\n");
	for (double ang=0; ang<2*M_PI; ang+=0.01) {
		entHandle newVertex, newEdge, newRegion;
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



// start an orbit at a random node
{
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
		entHandle currentTet=NULL;
		vect3d currentPosition(x,y,z);
		vect3d currentVelocity, zHat(0.,0.,1.);
		double potential;
		iMesh_getDblData(mesh, ents0d[iSelectedNode], potential_tag, &potential,
				&ierr);
		CHECK("Failure getting potential tag");
		currentVelocity = currentPosition.cross(zHat);
		currentVelocity /= currentVelocity.norm();
		currentVelocity *= sqrt(-potential);
		currentVelocity *= multiplier;
		vect3d eField(0.,0.,0.);
		int eField_alloc = sizeof(vect3d);
		int eField_size = sizeof(vect3d);
		vector<vect3d> eFields, vertexVectors;
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
		cout << "Initial radius=" << currentPosition.norm() << endl;
		for (double t=0; t<tMax; t+=dt) {
			nSteps++;
			currentPosition += dt*currentVelocity;
			inNewTet = !checkIfInTet(currentPosition, vertexVectors);
			if (inNewTet) {
				nNewTet++;
				// determine which tet currentPosition is in
				entHandle *entities = NULL;
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
				// assert(inNewTet == false);
				if(entities) free (entities);
				entities_alloc = 0;

				if (inNewTet==true) {
					cout << "Failed to identify current tet: currentPosition ="
							<< currentPosition << endl;
					cout << "nSteps=" << nSteps << ", nNewTet=" << nNewTet << endl;
					break;
				}

				iMesh_setIntData(mesh, currentTet, visited_tag, nOrbits, &ierr);
				CHECK("Failure setting visited tag");

				// get coordinates and field at vertices
				// TODO: unify with getVertexVectors
				entHandle *vertices = NULL;
				int vertices_alloc = 0, vertices_size;

				iMesh_getEntAdj(mesh, currentTet, iBase_VERTEX, &vertices, &vertices_alloc,
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

					// iMesh_getData(mesh, vertices[i], eField_tag, &eField,
					// &eField_alloc, &eField_size, &ierr);
					// CHECK("Failure getting eField tag");
					eFields[i] = eField;
				}
				if(vertices) free (vertices);
				vertices_alloc = 0;
			}

			assert(vertexVectors.size()==nVertices);
			vector<double> vertexWeights = getVertexWeights(currentPosition,
					vertexVectors);
			vect3d currentAcceleration(0.,0.,0.);
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
			// fprintf(outFile, "%f %f %f %d\n", currentPosition[0], currentPosition[1],
			// currentPosition[2], nNewTet);
		}
		cout << "Final radius=" << currentPosition.norm() << endl;
	}
	fclose(outFile);
}
