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
	adjacentFacesMap = getAdjacentsMap(iBase_REGION, iBase_FACE);
	adjacentVertsMap = getAdjacentsMap(iBase_REGION, iBase_VERTEX);
	adjacentTetsToFaceMap = getAdjacentsMap(iBase_FACE, iBase_REGION);

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

iBase_EntityHandle Mesh::findTet(Eigen::Vector3d oldPosition,
		Eigen::Vector3d position,
		iBase_EntityHandle adjacentTet, bool *tetFound, bool isTet) {
	iBase_EntityHandle tet;
	iBase_EntityHandle *entities = NULL;
	int entities_alloc=0, entities_size=0;
	int ierr;
	std::vector<iBase_EntityHandle> ents;
	std::vector<iBase_EntityHandle> faces;

	if (isTet) {
		iBase_EntityHandle faceCrossed = this->findFaceCrossed(
				adjacentTet, oldPosition, position);
		if (faceCrossed) {
			ents = adjacentTetsToFaceMap[faceCrossed];
//			ents = this->getAdjacentEntities(faceCrossed, iBase_REGION);
			for (int i=0; i<ents.size(); i++) {
				if (this->checkIfInTet(position, meshInstance, ents[i])) {
					tet = ents[i];
					*tetFound = true;
					return tet;
				}
			}
		}
//	iMesh_getEnt2ndAdj(meshInstance, adjacentTet, iBase_VERTEX,
//			iBase_REGION, &entities, &entities_alloc,
//			&entities_size, &ierr);
//	CHECK("Getting regions adjacent to entity failed");
		ents = adjacentTetsMap[adjacentTet];
	} else {
		ents = this->getAdjacentEntities(adjacentTet, iBase_REGION);
//		iMesh_getEntAdj(meshInstance, adjacentTet,
//				iBase_REGION, &entities, &entities_alloc,
//				&entities_size, &ierr);
//		CHECK("Getting regions adjacent to entity failed");
//		for (int i=0; i<entities_size; i++) {
//			ents.push_back(entities[i]);
//		}
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
	if (!*tetFound)
		tet = adjacentTet;

	return tet;
}

std::vector<iBase_EntityHandle> Mesh::getEntities(int dimension) {
	int ierr;
	iBase_EntityHandle *ents = NULL;
	int ents_alloc = 0, ents_size;
	iMesh_getEntities(meshInstance, rootEntitySet,
			dimension, iMesh_ALL_TOPOLOGIES,
			&ents, &ents_alloc, &ents_size, &ierr);
	CHECK("Couldn't get vertex entities");
	std::vector<iBase_EntityHandle> elements(ents_size);
	for (int i = 0; i < ents_size; i++) {
		elements[i] = ents[i];
	}
	if (ents) free(ents);
	ents_alloc = 0;
	return elements;
}

std::vector<iBase_EntityHandle> Mesh::getVertices(
		iBase_EntityHandle element) {
	return this->getAdjacentEntities(element, iBase_VERTEX);
}

std::vector<iBase_EntityHandle> Mesh::getAdjacentEntities(
		iBase_EntityHandle element, int dimension) {
	int ierr;
	iBase_EntityHandle *elements = NULL;
	int elements_alloc = 0, elements_size;

	iMesh_getEntAdj(meshInstance, element, dimension,
			&elements, &elements_alloc, &elements_size, &ierr);
	CHECK("Getting vertices adjacent to entity failed");

	std::vector<iBase_EntityHandle> elementHandles(elements_size);
	for (int i=0; i<elements_size; i++) {
		elementHandles[i] = elements[i];
	}
	if(elements) free (elements);
	elements_alloc = 0;

	return elementHandles;
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
	int dimension;
	iMesh_getEntType(meshInstance, entity, &dimension, &ierr);
	CHECK("Failure getting entity type");
	int nVerts = dimension + 1;
	std::vector<Eigen::Vector3d> vertexVectors(nVerts);
	std::vector<iBase_EntityHandle> vertices(nVerts);
	if (dimension == 3) {
		vertices = adjacentVertsMap[entity];
	} else {
		vertices = this->getAdjacentEntities(entity, iBase_VERTEX);
	}
	assert(vertexVectors.size()==vertices.size());
	if (useMap && dimension==3) {
		vertexVectors = vertexVectorsMap[entity];
	} else {
		for (int i=0; i<vertices.size(); i++) {
			vertexVectors[i] = Mesh::getCoordinates(vertices[i]);
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

bool Mesh::checkIfIntersectsTriangle(Eigen::Vector3d previousPosition,
		Eigen::Vector3d currentPosition,
		std::vector<Eigen::Vector3d> vertexVectors) {
	std::vector<Eigen::Vector3d> ray(2);
	ray[0] = previousPosition;
	ray[1] = currentPosition;
	Eigen::Vector3d intersectionPoint;
	int intersectionStatus =
			intersect_RayTriangle(ray, vertexVectors, &intersectionPoint);
	if (intersectionStatus == 1) {
		return true;
	} else {
		return false;
	}
}

//std::vector<Eigen::Vector3d> Mesh::getEdgeVectors(
//		std::vector<Eigen::Vector3d> vertexVectors) {
//	int nVerts = vertexVectors.size();
//	int nEdges=nVerts-1;
//	std::vector<Eigen::Vector3d> edgeVectors(nEdges);
//	switch (nVerts) {
//	case 4:
//		edgeVectors[2] = vertexVectors[3]-vertexVectors[0];
//	case 3:
//		edgeVectors[1] = vertexVectors[2]-vertexVectors[0];
//	case 2:
//		edgeVectors[0] = vertexVectors[1]-vertexVectors[0];
//		break;
//	}
//
//	return edgeVectors;
//}


iBase_EntityHandle Mesh::findFaceCrossed(iBase_EntityHandle previousElement,
		Eigen::Vector3d previousPosition, Eigen::Vector3d currentPosition) {
	iBase_EntityHandle faceCrossed=NULL;
	std::vector<iBase_EntityHandle> adjacentFaces =
			adjacentFacesMap[previousElement];
//	std::vector<iBase_EntityHandle> adjacentFaces =
//			getAdjacentEntities(previousElement,iBase_FACE);

	for (int i=0; i<adjacentFaces.size(); i++) {
		std::vector<Eigen::Vector3d> vertexVectors =
				this->getVertexVectors(adjacentFaces[i], true);
		bool intersectsTriangle = this->checkIfIntersectsTriangle(previousPosition,
				currentPosition, vertexVectors);
		if (intersectsTriangle)
			faceCrossed = adjacentFaces[i];
	}

	return faceCrossed;
}

iBase_TagHandle Mesh::getTagHandle(std::string tagName) {
	int ierr;
	iBase_TagHandle tag=0;
	iMesh_getTagHandle(meshInstance, tagName.c_str(),
			&tag, &ierr, tagName.length());
	return tag;
}

iBase_TagHandle Mesh::createTagHandle(std::string tagName, int size, int type) {
	int ierr;
	iBase_TagHandle tag=0;
	iMesh_createTag(meshInstance, tagName.c_str(),
			size, type, &tag, &ierr, (int)tagName.length());
	return tag;
}

Eigen::Vector3d Mesh::getNormalVector(iBase_EntityHandle face,
		Eigen::Vector3d point) {
	std::vector<Eigen::Vector3d> vertexVectors =
			this->getVertexVectors(face);
	assert(3 == vertexVectors.size());
	std::vector<Eigen::Vector3d> edgeVectors(vertexVectors.size()-1);

	edgeVectors[0] = vertexVectors[1]-vertexVectors[0];
	edgeVectors[1] = vertexVectors[2]-vertexVectors[0];

	Eigen::Vector3d surfaceVector = edgeVectors[0].cross(edgeVectors[1])/2.;
	if (point==Eigen::Vector3d(0.,0.,0.)) {
		std::vector<iBase_EntityHandle> tets =
				this->getAdjacentEntities(face,iBase_REGION);
		// TODO: handle interior faces with two adjacent tets?
		std::vector<Eigen::Vector3d> tetVertexVectors =
				this->getVertexVectors(tets[0]);
		int nPoints=0;
		for (int i=0; i<tetVertexVectors.size(); i++) {
			bool noMatches=true;
			for (int j=0; j<vertexVectors.size(); j++) {
				if (tetVertexVectors[i]==vertexVectors[j])
					noMatches=false;
			}
			if (noMatches) {
				point = tetVertexVectors[i];
				nPoints++;
			}
		}
		assert(nPoints==1);
	}
	assert(point!=Eigen::Vector3d(0.,0.,0.));
	Eigen::Vector3d referenceVector = point - vertexVectors[0];

	if (referenceVector.dot(surfaceVector)<0)
		surfaceVector *= -1.;

	return surfaceVector/surfaceVector.norm();
}

Eigen::Vector3d Mesh::getVertexNormalVector(iBase_EntityHandle vertex,
		Field<int> faceTypeField) {
	std::vector<iBase_EntityHandle> faces =
			this->getAdjacentEntities(vertex, iBase_FACE);
	// TODO: could pre-allocate normalVectors
	std::vector<Eigen::Vector3d> normalVectors;
	std::vector<iBase_EntityHandle>::iterator faceIter = faces.begin();
	for (; faceIter!=faces.end(); ++faceIter) {
		int faceType = faceTypeField.getField(*faceIter);
		if (faceType>0) {
			Eigen::Vector3d normalVector =
					this->getNormalVector(*faceIter);
			normalVectors.push_back(normalVector);
		}
	}
	Eigen::Vector3d normalVector(0.,0.,0.);
	std::vector<Eigen::Vector3d>::iterator normalIter =
			normalVectors.begin();
	for (; normalIter!=normalVectors.end(); ++normalIter) {
		normalVector += *normalIter;
	}
	normalVector /= normalVector.norm();
	Eigen::Vector3d testCoord = this->getCoordinates(vertex);
	return normalVector;
}
