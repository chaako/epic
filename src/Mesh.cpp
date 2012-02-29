/*
 * Mesh.cpp
 *
 *  Created on: Feb 10, 2012
 *      Author: chaako
 */

#include "epic.h"
#include "Mesh.h"

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

	// store adjacency info for fast access
	adjacentTetsMap = getAdjacentsMap(iBase_REGION, iBase_REGION,
			iBase_VERTEX);
	adjacentVertsMap = getAdjacentsMap(iBase_REGION, iBase_VERTEX);

	// store vertex coordinates for fast access
	std::map<iBase_EntityHandle,std::vector<iBase_EntityHandle> >::iterator iter;
	for (iter=adjacentVertsMap.begin(); iter!=adjacentVertsMap.end(); ++iter) {
		vertexVectorsMap[iter->first] = Mesh::getVertexVectors(iter->first, false);
		assert(vertexVectorsMap[iter->first].size()==4);
	}

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

void Mesh::save(std::string outputMeshFile) {
	char *options = NULL;
	int options_len = 0;
	int ierr;
	// destroy eField tag since VisIt doesn't understand
	iBase_EntityHandle *ents0d = NULL;
	int ents0d_alloc = 0, ents0d_size;
	iBase_TagHandle eField_tag;
	std::string tagName = "eField";
	iMesh_getTagHandle(meshInstance, tagName.c_str(),
			&eField_tag, &ierr, tagName.length());
	CHECK("Failed to get eField tag");
	iBase_TagHandle eFieldX_tag, eFieldY_tag, eFieldZ_tag;
	iMesh_createTag(meshInstance, "eFieldX", 1, iBase_DOUBLE,
			&eFieldX_tag, &ierr, 7);
	CHECK("Failure creating eField tag");
	iMesh_createTag(meshInstance, "eFieldY", 1, iBase_DOUBLE,
			&eFieldY_tag, &ierr, 7);
	CHECK("Failure creating eField tag");
	iMesh_createTag(meshInstance, "eFieldZ", 1, iBase_DOUBLE,
			&eFieldZ_tag, &ierr, 7);
	CHECK("Failure creating eField tag");
	iMesh_getEntities(meshInstance, rootEntitySet,
			iBase_VERTEX, iMesh_ALL_TOPOLOGIES,
			&ents0d, &ents0d_alloc, &ents0d_size, &ierr);
	CHECK("Couldn't get vertex entities");
	for (int i = 0; i < ents0d_size; i++) {
		std::vector<iBase_EntityHandle> superCellFaces =
				getSuperCellFaces(meshInstance, ents0d[i]);
		Eigen::Vector3d eField(0.,0.,0.);
		Eigen::Vector3d *eField_ptr = &eField;
		int eField_alloc = sizeof(Eigen::Vector3d);
		int eField_size = sizeof(Eigen::Vector3d);
		iMesh_getData(meshInstance, ents0d[i], eField_tag, &eField_ptr,
				&eField_alloc, &eField_size, &ierr);
		CHECK("Failure getting eField tag");
		iMesh_setDblData(meshInstance, ents0d[i], eFieldX_tag,
				(double)eField[0], &ierr);
		CHECK("Failure setting eFieldX tag");
		iMesh_setDblData(meshInstance, ents0d[i], eFieldY_tag,
				(double)eField[1], &ierr);
		CHECK("Failure setting eFieldY tag");
		iMesh_setDblData(meshInstance, ents0d[i], eFieldZ_tag,
				(double)eField[2], &ierr);
		CHECK("Failure setting eFieldZ tag");
	}
	iMesh_destroyTag(meshInstance, eField_tag, 1, &ierr);
	CHECK("Failed to destroy eField tag");
	/* save the mesh */
	iMesh_save(meshInstance, rootEntitySet, outputMeshFile.c_str(),
			options, &ierr, outputMeshFile.length(), options_len);
	CHECK("Save failed");
}

std::map<iBase_EntityHandle,std::vector<iBase_EntityHandle> >
Mesh::getAdjacentsMap(int keyEntityType, int valueEntityType,
		int bridgeEntityType) {
	std::map<iBase_EntityHandle,std::vector<iBase_EntityHandle> > adjacentsMap;
	iBase_EntityHandle *ents3d = NULL;
	int ents3d_alloc = 0, ents3d_size;
	int ierr;
	iMesh_getEntities(meshInstance, rootEntitySet, keyEntityType,
			iMesh_ALL_TOPOLOGIES,
			&ents3d, &ents3d_alloc, &ents3d_size, &ierr);
	CHECK("Couldn't get entities");
	for(int i=0; i<ents3d_size; i++) {
		iBase_EntityHandle *entities = NULL;
		int entities_alloc = 0, entities_size;
		if (bridgeEntityType>=0) {
			iMesh_getEnt2ndAdj(meshInstance, ents3d[i], bridgeEntityType,
					valueEntityType, &entities, &entities_alloc,
					&entities_size, &ierr);
			CHECK("Getting second adjacency failed");
		} else {
			iMesh_getEntAdj(meshInstance, ents3d[i],
					valueEntityType, &entities, &entities_alloc,
					&entities_size, &ierr);
			CHECK("Getting adjacency failed");
		}
		std::vector<iBase_EntityHandle> handles;
		for (int j=0; j<entities_size; j++) {
			handles.push_back(entities[j]);
		}
		adjacentsMap[ents3d[i]] = handles;
		if(entities) free (entities);
		entities_alloc = 0;
	}
	if(ents3d_alloc) free (ents3d);
	ents3d_alloc = 0;

	return adjacentsMap;
}

Eigen::Vector3d Mesh::getCoordinates(iBase_EntityHandle node, bool useMap) {
	Eigen::Vector3d coordinates(0.,0.,0.);
	double x,y,z;
	int ierr;
	iMesh_getVtxCoord(meshInstance, node, &x, &y, &z, &ierr);
	CHECK("Failure getting vertex coordinates");
	coordinates << x, y, z;
	return coordinates;
}

iBase_EntityHandle Mesh::findTet(Eigen::Vector3d position,
		iBase_EntityHandle adjacentTet, bool *tetFound, bool isTet) {
	iBase_EntityHandle tet;
	iBase_EntityHandle *entities = NULL;
	int entities_alloc=0, entities_size=0;
	int ierr;
	std::vector<iBase_EntityHandle> ents;

	if (isTet) {
//	iMesh_getEnt2ndAdj(meshInstance, adjacentTet, iBase_VERTEX,
//			iBase_REGION, &entities, &entities_alloc,
//			&entities_size, &ierr);
//	CHECK("Getting regions adjacent to entity failed");
		ents = adjacentTetsMap[adjacentTet];
	} else {
		iMesh_getEntAdj(meshInstance, adjacentTet,
				iBase_REGION, &entities, &entities_alloc,
				&entities_size, &ierr);
		CHECK("Getting regions adjacent to entity failed");
	}
	for (int i=0; i<entities_size; i++) {
		ents.push_back(entities[i]);
	}
	if(entities) free (entities);
	entities_alloc = 0;
	for (int i=0; i<ents.size(); i++) {
		if (Mesh::checkIfInTet(position, meshInstance, ents[i])) {
			tet = ents[i];
			// TODO: could throw error if tetFound is false rather than pass
			*tetFound = true;
			break;
		}
	}

	return tet;
}

std::vector<iBase_EntityHandle> Mesh::getVertices(
		iBase_EntityHandle element) {
	int ierr;
	iBase_EntityHandle *vertices = NULL;
	int vertices_alloc = 0, vertices_size;
	int nVertices=4;
	std::vector<iBase_EntityHandle> vertexHandles(nVertices);

	iMesh_getEntAdj(meshInstance, element, iBase_VERTEX,
			&vertices, &vertices_alloc, &vertices_size, &ierr);
	CHECK("Getting vertices adjacent to entity failed");

	for (int i=0; i<vertices_size; i++) {
		vertexHandles[i] = vertices[i];
	}
	if(vertices) free (vertices);
	vertices_alloc = 0;

	return vertexHandles;
}

iBase_EntityHandle Mesh::getRandomVertex() {
	iBase_EntityHandle *ents0d = NULL;
	int ents0d_alloc = 0, ents0d_size;
	int ierr;

	iMesh_getEntities(meshInstance, rootEntitySet, iBase_VERTEX,
			iMesh_ALL_TOPOLOGIES,
			&ents0d, &ents0d_alloc, &ents0d_size, &ierr);
	CHECK("Couldn't get vertex entities");

	srand(999);
	int iSelectedNode = rand() % ents0d_size;

	if(ents0d_alloc) free (ents0d);
	ents0d_alloc = 0;

	return ents0d[iSelectedNode];
}

std::vector<Eigen::Vector3d> Mesh::getVertexVectors(iBase_EntityHandle entity,
		bool useMap) {
	int ierr;
	int nVerts = 4;
	std::vector<Eigen::Vector3d> vertexVectors(nVerts);
	std::vector<iBase_EntityHandle> vertices = adjacentVertsMap[entity];
	assert(vertexVectors.size()==vertices.size());
	if (useMap) {
		vertexVectors = vertexVectorsMap[entity];
	} else {
		for (int i=0; i<vertices.size(); i++) {
	//		double x, y, z;
	//		iMesh_getVtxCoord(mesh, vertices[i], &x, &y, &z, &ierr);
	//		CHECK("Failure getting vertex coordinates");
			vertexVectors[i] = Mesh::getCoordinates(vertices[i]);

	//		Eigen::Vector3d vertexVector(x,y,z);
	//		vertexVectors[i] = vertexVector;
		}
	}

	return vertexVectors;
}

bool Mesh::checkIfInTet(Eigen::Vector3d currentPosition, iMesh_Instance mesh,
		iBase_EntityHandle element) {
	std::vector<Eigen::Vector3d> vertexVectors = Mesh::getVertexVectors(element);

	return Mesh::checkIfInTet(currentPosition, vertexVectors);
}

bool Mesh::checkIfInTet(Eigen::Vector3d currentPosition,
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

