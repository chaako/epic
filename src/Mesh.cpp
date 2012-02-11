/*
 * Mesh.cpp
 *
 *  Created on: Feb 10, 2012
 *      Author: chaako
 */

#include "epic.h"
//#include "Mesh.h"

Mesh::Mesh(std::string inputMeshFile) {
	char *options = NULL;
	int options_len = 0;
	int ierr;
	iMesh_Instance mesh;
	iBase_EntitySetHandle root;

	// create mesh instance
	iMesh_newMesh(options, &mesh, &ierr, options_len);
	CHECK("Problems instantiating interface.");
	meshInstance = mesh;

	// get handle to root set, which contains all elements
	iMesh_getRootSet(mesh, &root, &ierr);
	CHECK("Problems getting root set");
	rootEntitySet = root;

//	iMesh_load(mesh, root, inputMeshFile.c_str(), options, &ierr,
//			strlen(inputMeshFile.c_str()), options_len);
	// FMDB's importVTK can't handle our tags, so use custom version
	ierr = custom_importVTK((mMesh *)mesh, inputMeshFile.c_str());
	CHECK("Load failed");
}

Mesh::~Mesh() {
	int ierr;
	iMesh_dtor(meshInstance, &ierr);
	CHECK("Failed to destroy interface");
}


void Mesh::printElementNumbers() {
	int dim, num;
	int ierr;
	/* report the number of elements of each dimension */
	for (dim = iBase_VERTEX; dim <= iBase_REGION; dim++) {
		iMesh_getNumOfType(meshInstance, rootEntitySet, dim, &num, &ierr);
		CHECK("Failure in getNumOfType");
		printf("Number of %d-dimensional elements = %d\n", dim, num);
	}
}
