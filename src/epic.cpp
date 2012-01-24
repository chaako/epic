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
#include <math.h>

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
	iBase_TagHandle code_tag, guard_tag, potential_tag;

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
		potential = 1./sqrt(x*x+y*y+z*z);
		iMesh_setDblData(mesh, ents0d[i], potential_tag, potential,
				&ierr);
		CHECK("Failure setting potential tag");
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

//	cout << "Hello World" << endl; /* prints Hello World */
	return 0;
}
