/*
 ============================================================================
 Name        : epic.cpp
 Author      : Christian Bernt Haakonsen
 Version     :
 Copyright   : Copyright
 Description : Hello World in C++,
 ============================================================================
 */

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "iMesh.h"
#include <math.h>

#include "Eigen/Dense"

#include "epic.h"

#define CHECK(a) if (iBase_SUCCESS != ierr) printf("%s\n", a), exit(ierr)

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
		potential = 1./sqrt(x*x+y*y+z*z);
		iMesh_setDblData(mesh, ents0d[i], potential_tag, potential,
				&ierr);
		CHECK("Failure setting potential tag");
	}

	// create eField tag
//	iMesh_createTag(mesh, "eField", sizeof(Eigen::Vector3d), iBase_BYTES,
//			&eField_tag, &ierr, 6);
//	CHECK("Failure creating eField tag");
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
	// TODO: add handling of boundaries with non-zero E
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
			volume += getVolumeBetweenPointAndFace(mesh, point,
					superCellFaces[j]);
			eField -= potential*surfaceVector;
		}
		eField /= volume;

//		iMesh_setData(mesh, ents0d[i], eField_tag, &eField,
//				sizeof(Eigen::Vector3d), &ierr);
//		CHECK("Failure setting eField tag");
		iMesh_setDblData(mesh, ents0d[i], eFieldX_tag, (double)eField[0], &ierr);
		CHECK("Failure setting eFieldX tag");
		iMesh_setDblData(mesh, ents0d[i], eFieldY_tag, (double)eField[1], &ierr);
		CHECK("Failure setting eFieldY tag");
		iMesh_setDblData(mesh, ents0d[i], eFieldZ_tag, (double)eField[2], &ierr);
		CHECK("Failure setting eFieldZ tag");
	}
	if (ents0d) free(ents0d);
	ents0d_alloc = 0;

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

	Eigen::Vector3d v(1,2,3);
	Eigen::Vector3d w(0,1,2);

	std::cout << "Dot product: " << v.dot(w) << endl;
	double dp = v.adjoint()*w; // automatic conversion of the inner product to a scalar
	std::cout << "Dot product via a matrix product: " << dp << endl;
	std::cout << "Cross product:\n" << v.cross(w) << endl;
	std::cout << "Element:\n" << (double)v[1] << endl;

	return 0;
}

Eigen::Vector3d getSurfaceVector(iMesh_Instance mesh, Eigen::Vector3d point,
		iBase_EntityHandle face) {
	int ierr;
	iBase_EntityHandle *vertices = NULL;
	int vertices_alloc = 0, vertices_size;
	std::vector<Eigen::Vector3d> vertexVectors, edgeVectors;
	Eigen::Vector3d surfaceVector, referenceVector;

	iMesh_getEntAdj(mesh, face, iBase_VERTEX,  &vertices, &vertices_alloc,
			&vertices_size, &ierr);
	CHECK("Getting vertices adjacent to face failed");
	for (int i=0; i<vertices_size; i++) {
		double x=0., y=0., z=0.;
		iMesh_getVtxCoord(mesh, vertices[i], &x, &y, &z, &ierr);
		CHECK("Failure getting vertex coordinates");
		Eigen::Vector3d vertexVector(x,y,z);
		vertexVectors.push_back(vertexVector);
	}
	if(vertices) free (vertices);
	vertices_alloc = 0;

	assert(3 == vertexVectors.size());
	edgeVectors.push_back(vertexVectors[1]-vertexVectors[0]);
	edgeVectors.push_back(vertexVectors[2]-vertexVectors[1]);

	surfaceVector = edgeVectors[0].cross(edgeVectors[1]);
	referenceVector = vertexVectors[0] - point;

	if (referenceVector.dot(surfaceVector) < 0 )
		surfaceVector *= -1;
	return surfaceVector;
}

vector<iBase_EntityHandle> getSuperCellFaces(iMesh_Instance mesh,
		iBase_EntityHandle vertex) {
	int ierr;
	iBase_EntityHandle *faces = NULL;
	int faces_alloc = 0, faces_size;
	std::vector<iBase_EntityHandle> superCellFaces;

	iMesh_getEnt2ndAdj(mesh, vertex, iBase_REGION,
			iBase_FACE, &faces, &faces_alloc,
			&faces_size, &ierr);
	CHECK("Failure in getEnt2ndAdj");
	for (int i = 0; i < faces_size; i++) {
		iBase_EntityHandle *vertices = NULL;
		int vertices_alloc = 0, vertices_size;
		bool onSuperCellSurface = true;

		iMesh_getEntAdj(mesh, faces[i], iBase_VERTEX,  &vertices, &vertices_alloc,
				&vertices_size, &ierr);
		CHECK("Getting vertices adjacent to face failed");

		for (int j=0; j<vertices_size; j++) {
			if (vertices[j]==vertex)
				onSuperCellSurface = false;
			// TODO: if vertex is on boundary, should probably include
			//       the corresponding faces to enclose volume
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
		CHECK("Failure getting potential tag");

		averageDblData += dblData;
	}
	averageDblData /= (double)vertices_size;

	if(vertices) free (vertices);
	vertices_alloc = 0;

	return averageDblData;
}

double getVolumeBetweenPointAndFace(iMesh_Instance mesh, Eigen::Vector3d point,
		iBase_EntityHandle face) {
	int ierr;
	iBase_EntityHandle *vertices = NULL;
	int vertices_alloc = 0, vertices_size;
	std::vector<Eigen::Vector3d> vertexVectors, edgeVectors;

	vertexVectors.push_back(point);

	iMesh_getEntAdj(mesh, face, iBase_VERTEX,  &vertices, &vertices_alloc,
			&vertices_size, &ierr);
	CHECK("Getting vertices adjacent to face failed");

	for (int i=0; i<vertices_size; i++) {
		double x=0., y=0., z=0.;
		iMesh_getVtxCoord(mesh, vertices[i], &x, &y, &z, &ierr);
		CHECK("Failure getting vertex coordinates");
		Eigen::Vector3d vertexVector(x,y,z);
		vertexVectors.push_back(vertexVector);
	}
	if(vertices) free (vertices);
	vertices_alloc = 0;

	assert(4 == vertexVectors.size());
	edgeVectors.push_back(vertexVectors[1]-vertexVectors[0]);
	edgeVectors.push_back(vertexVectors[2]-vertexVectors[1]);
	edgeVectors.push_back(vertexVectors[3]-vertexVectors[2]);
	return fabs(
			vertexVectors[3].dot( vertexVectors[1].cross(vertexVectors[2]) )
			)/6.;
}
