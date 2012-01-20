/*
 ============================================================================
 Name        : epic.cpp
 Author      : Christian Bernt Haakonsen
 Version     :
 Copyright   : Copyright
 Description : Hello World in C++,
 ============================================================================
 */

//#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "iMesh.h"

#include "epic.h"

#define CHECK(a) if (iBase_SUCCESS != ierr) printf("%s\n", a), exit(ierr)

using namespace std;

int main(int argc, char *argv[]) {
	char *options = NULL;
	int options_len = 0;
	int dim, num, ierr;
	iMesh_Instance mesh;
	iBase_EntitySetHandle root;

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

	/* save the mesh */
	iMesh_save(mesh, root, argv[2], options, &ierr,
			strlen(argv[2]), options_len);
	CHECK("Save failed");

	iMesh_dtor(mesh, &ierr);
	CHECK("Failed to destroy interface");

//	cout << "Hello World" << endl; /* prints Hello World */
	return 0;
}
