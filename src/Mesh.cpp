/*
 * Mesh.cpp
 *
 *  Created on: Feb 10, 2012
 *      Author: chaako
 */

#include "epic.h"
#include "Mesh.h"

Mesh::Mesh(string inputMeshFile) {
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

	if (inputMeshFile.find(".vtk")==string::npos) {
		iMesh_load(mesh, root, inputMeshFile.c_str(), options, &ierr,
				strlen(inputMeshFile.c_str()), options_len);
		vtkInputMesh = false;
		// recreate vector field tags and destroy component tags
		try {
			this->convertComponentTagsToVectorTag("eField");
			this->convertComponentTagsToVectorTag("ionVelocity");
		} catch (string& error) {
			cout << error << endl;
		}
	} else {
		// FMDB's importVTK can't handle our tags, so use custom version
		ierr = custom_importVTK((mMesh *)mesh, inputMeshFile.c_str());
		CHECK("Load failed");
		vtkInputMesh = true;
	}
//	// TODO: debugging
//	cout << inputMeshFile << " " << vtkInputMesh << endl;

	// store adjacency info for fast access
	adjacentTetsMap = getAdjacentsMap(iBase_REGION, iBase_REGION,
			iBase_VERTEX);
	adjacentFacesMap = getAdjacentsMap(iBase_REGION, iBase_FACE);
	adjacentVertsMap = getAdjacentsMap(iBase_REGION, iBase_VERTEX);
	surroundingVertsMap = getSurroundingVerticesMap();
	adjacentTetsToFaceMap = getAdjacentsMap(iBase_FACE, iBase_REGION);

	// store vertex coordinates for fast access
	map<entHandle,vector<entHandle> >::iterator iter;
	for (iter=adjacentVertsMap.begin(); iter!=adjacentVertsMap.end(); ++iter) {
		vertexVectorsMap[iter->first] = Mesh::getVertexVectors(iter->first, false);
		if (vertexVectorsMap[iter->first].size()!=4)
			throw string("Error constructing vertexVectorsMap in Mesh.cpp");
	}

	previousCoordsToBasisElement = NULL;

	for (int i=0; i<=NDIM; i++) {
		entitiesVectors.push_back(this->getEntities(i));
		for (int j=0; j<entitiesVectors[i].size(); j++) {
			indicesOfEntities[entitiesVectors[i][j]] = j;
			dimensionsOfEntities[entitiesVectors[i][j]] = i;
		}
	}

	// Can't combine this loop with above because getAdjacentEntitiesIndices
	// uses indicesOfEntities
	for (int i=0; i<=NDIM; i++) {
		adjacentEntitiesVectors.push_back(
				vector<vector<vector<int> > >(entitiesVectors[i].size()));
		for (int j=0; j<entitiesVectors[i].size(); j++) {
			for (int k=0; k<=NDIM; k++) {
				adjacentEntitiesVectors[i][j].push_back(
						this->getAdjacentEntitiesIndices(j,i,k));
			}
		}
	}

	for (int i=0; i<entitiesVectors[iBase_VERTEX].size(); i++) {
		vect3d position =
				this->getCoordinates(entitiesVectors[iBase_VERTEX][i]);
		verticesPositions.push_back(position);
	}

	for (int i=0; i<entitiesVectors[iBase_REGION].size(); i++) {
		entHandle regionHandle = entitiesVectors[iBase_REGION][i];
		vector<entHandle> &surroundingVerts =
				surroundingVertsMap[regionHandle];
		vector<int> surroundingVertices(surroundingVerts.size());
		for (int j=0; j<surroundingVerts.size(); j++) {
			surroundingVertices[j] = indicesOfEntities[surroundingVerts[j]];
		}
		verticesSurroundingRegions.push_back(surroundingVertices);

		vector<entHandle> &surroundingRegs =
				adjacentTetsMap[regionHandle];
		vector<int> surroundingRegions(surroundingRegs.size());
		for (int j=0; j<surroundingRegs.size(); j++) {
			surroundingRegions[j] = indicesOfEntities[surroundingRegs[j]];
		}
		regionsSurroundingRegions.push_back(surroundingRegions);

		vector<int> &regionVertices =
				adjacentEntitiesVectors[iBase_REGION][i][iBase_VERTEX];
		vector<int> oppositeRegions(regionVertices.size());
		vector<int> &regionFaces =
				adjacentEntitiesVectors[iBase_REGION][i][iBase_FACE];
		for (int k=0; k<regionVertices.size(); k++) {
			int oppositeFace = -1;
			for (int j=0; j<regionFaces.size(); j++) {
				vector<int> &faceVertices =
						adjacentEntitiesVectors[iBase_FACE][regionFaces[j]][iBase_VERTEX];
				bool correctFace=true;
				for (int l=0; l<faceVertices.size(); l++) {
					if (faceVertices[l]==regionVertices[k])
						correctFace=false;
				}
				if (correctFace)
					oppositeFace = regionFaces[j];
			}
			vector<int> &faceRegions =
					adjacentEntitiesVectors[iBase_FACE][oppositeFace][iBase_REGION];
			if (faceRegions.size()==2) {
				oppositeRegions[k] = (faceRegions[0]!=i) ? faceRegions[0] : faceRegions[1];
			} else {
				oppositeRegions[k] = faceRegions[0];
			}
		}
		regionsOppositeVertices.push_back(oppositeRegions);

		vector<vect3d> vVs = this->getVertexVectors(regionHandle);
		positionsToBases.push_back(
				this->calculatePositionToBasesMatrix(vVs));
	}

//	cout << entitiesVectors[3][adjacentEntitiesVectors[2][15][3][1]] << " " <<
//			adjacentEntitiesVectors[2][15][3][1] << endl;
//	cout << (this->getAdjacentEntities(entitiesVectors[2][15],3))[1] << endl;

//	allVertices = this->getEntities(iBase_VERTEX);
//	for(int i=0; i<allVertices.size(); i++) {
//		indexOfVertices[allVertices[i]] = i;
//		indexOfEntity[allVertices[i]] = i;
//	}
//	allFaces = this->getEntities(iBase_FACE);
//	for(int i=0; i<allFaces.size(); i++) {
//		indexOfFaces[allFaces[i]] = i;
//		indexOfEntity[allFaces[i]] = i;
//	}
//	allElements = this->getEntities(iBase_REGION);
//	for(int i=0; i<allElements.size(); i++) {
//		indexOfElements[allElements[i]] = i;
//		indexOfEntity[allElements[i]] = i;
//	}

	vtkMesh_ptr = this->createVtkMesh();
	vtkCellTree_ptr = vtkSmartPointer<vtkCellTreeLocator>::New();
	vtkCellTree_ptr->SetDataSet(vtkMesh_ptr);
	vtkCellTree_ptr->BuildLocator();


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
#ifdef HAVE_MPI
		if (MPI::COMM_WORLD.Get_rank() == 0)
#endif
			printf("Number of %d-dimensional elements = %d\n", dim, num);
	}
}

void Mesh::save(string outputMeshFile) {
	char *options = NULL;
	int options_len = 0;
	int ierr;
	// destroy vector tags since VisIt doesn't understand
	this->convertVectorTagToComponentTags("eField");
	this->convertVectorTagToComponentTags("ionVelocity");

	iMesh_save(meshInstance, rootEntitySet, outputMeshFile.c_str(),
			options, &ierr, outputMeshFile.length(), options_len);
	CHECK("Save failed");

	// recreate vector tags and destroy component tags
	this->convertComponentTagsToVectorTag("eField");
	this->convertComponentTagsToVectorTag("ionVelocity");
}

iBase_TagHandle Mesh::createTag(string tagName, int size, int type) {
	int ierr;
	iBase_TagHandle tag;

	iMesh_createTag(meshInstance, tagName.c_str(), size, type,
			&tag, &ierr, tagName.length());
	if (ierr != iBase_SUCCESS)
		throw string("Failure creating tag named " + tagName);
	return tag;
}

iBase_TagHandle Mesh::getTag(string tagName) {
	int ierr;
	iBase_TagHandle tag;
	iMesh_getTagHandle(meshInstance, tagName.c_str(),
			&tag, &ierr, tagName.length());
	if (ierr != iBase_SUCCESS)
		throw ("Failed to get " + tagName + " tag.");

	return tag;
}


void Mesh::destroyTag(string tagName) {
	int ierr;
	iBase_TagHandle tag=this->getTag(tagName);
	iMesh_destroyTag(meshInstance, tag, 1, &ierr);
	CHECK("Failed to destroy tag");
}

void Mesh::convertVectorTagToComponentTags(string vectorTagName) {
	entHandle *ents0d = NULL;
	int ents0d_alloc = 0, ents0d_size;
	int ierr;
	iBase_TagHandle vectorTag=this->getTag(vectorTagName);

	iBase_TagHandle componentTagX, componentTagY, componentTagZ;
	componentTagX = this->createTag(vectorTagName+"X", 1, iBase_DOUBLE);
	componentTagY = this->createTag(vectorTagName+"Y", 1, iBase_DOUBLE);
	componentTagZ = this->createTag(vectorTagName+"Z", 1, iBase_DOUBLE);

	// TODO: clean this up, separating out into functions hiding iMesh details
	iMesh_getEntities(meshInstance, rootEntitySet,
			iBase_VERTEX, iMesh_ALL_TOPOLOGIES,
			&ents0d, &ents0d_alloc, &ents0d_size, &ierr);
	CHECK("Couldn't get vertex entities");
	for (int i = 0; i < ents0d_size; i++) {
		vector<entHandle> superCellFaces =
				this->getSuperCellFaces(ents0d[i]);
		vect3d vectorData(0.,0.,0.);
		vect3d *vectorData_ptr = &vectorData;
		int alloc = sizeof(vect3d);
		int size = sizeof(vect3d);
		iMesh_getData(meshInstance, ents0d[i], vectorTag, &vectorData_ptr,
				&alloc, &size, &ierr);
		CHECK("Failure getting vector data");
		iMesh_setDblData(meshInstance, ents0d[i], componentTagX,
				(double)vectorData[0], &ierr);
		CHECK("Failure setting component data");
		iMesh_setDblData(meshInstance, ents0d[i], componentTagY,
				(double)vectorData[1], &ierr);
		CHECK("Failure setting component data");
		iMesh_setDblData(meshInstance, ents0d[i], componentTagZ,
				(double)vectorData[2], &ierr);
		CHECK("Failure setting component data");
	}
	this->destroyTag(vectorTagName);
}

void Mesh::convertComponentTagsToVectorTag(string vectorTagName) {
	entHandle *ents0d = NULL;
	int ents0d_alloc = 0, ents0d_size;
	int ierr;
	iBase_TagHandle vectorTag;
	vectorTag = this->createTag(vectorTagName, sizeof(vect3d), iBase_BYTES);

	iBase_TagHandle componentTagX=this->getTag(vectorTagName+"X");
	iBase_TagHandle componentTagY=this->getTag(vectorTagName+"Y");
	iBase_TagHandle componentTagZ=this->getTag(vectorTagName+"Z");

	// TODO: clean this up, separating out into functions hiding iMesh details
	for (int i = 0; i < ents0d_size; i++) {
		double componentX, componentY, componentZ;
		iMesh_getDblData(meshInstance, ents0d[i], componentTagX,
				&componentX, &ierr);
		CHECK("Failure getting component");
		iMesh_getDblData(meshInstance, ents0d[i], componentTagY,
				&componentY, &ierr);
		CHECK("Failure getting component");
		iMesh_getDblData(meshInstance, ents0d[i], componentTagZ,
				&componentZ, &ierr);
		CHECK("Failure getting component");
		vect3d vectorData(componentX,componentY,componentZ);
		iMesh_setData(meshInstance, ents0d[i], vectorTag, &vectorData,
				sizeof(vect3d), &ierr);
		CHECK("Failure setting vector data");
	}

	this->destroyTag(vectorTagName+"X");
	this->destroyTag(vectorTagName+"Y");
	this->destroyTag(vectorTagName+"Z");
}

#ifdef MESHER
#include "FMDB.h"
void Mesh::classifyBoundariesForMeshRefinement(Field<int> faceTypeField){
	// meshAdapt requires mesh elements on the boundary to be classified
	// in terms of geometric elements.

	// Need to use native FMDB pointers since no iMesh interface
	mPart *part;
	FMDB_Mesh_GetPart((mMesh*)this->meshInstance, 0, part);
	typedef GEntity* (*entityBy_FP)(SGModel*,int,int);
	entityBy_FP cl_ptr = GM_entityByTag;

	// Start by setting all elements to be part of a volume geometry element
	for (int i=0; i<entitiesVectors.size(); i++) {
		for (int j=0; j<entitiesVectors[i].size(); j++) {
			entHandle entity = this->entitiesVectors[i][j];
			((mEntity*)entity)->classify(part->getGEntity(0,3,cl_ptr));
		}
	}

	// Take all faces, edges, and vertices of original surface mesh to be fixed
	int i=iBase_FACE;
	for (int j=0; j<entitiesVectors[i].size(); j++) {
		entHandle face = this->entitiesVectors[i][j];
		int faceType = faceTypeField.getField(face);
		if (faceType>1) {
			((mEntity*)face)->classify(part->getGEntity(0,2,cl_ptr));
			vector<entHandle> edges = this->getAdjacentEntities(face,iBase_EDGE);
			for (int k=0; k<edges.size(); k++) {
				((mEntity*)edges[k])->classify(part->getGEntity(0,1,cl_ptr));
			}
			vector<entHandle> vertices = this->getAdjacentEntities(face,iBase_VERTEX);
			for (int k=0; k<vertices.size(); k++) {
				((mEntity*)vertices[k])->classify(part->getGEntity(0,0,cl_ptr));
			}
		}
	}
}
#endif

map<entHandle,vector<entHandle> >
Mesh::getSurroundingVerticesMap() {
	map<entHandle,vector<entHandle> >
	surroundingVerticesMap;
	// TODO: should handle case where adjacentTetsMap isn't available yet
	if (adjacentTetsMap.begin()==adjacentTetsMap.end())
		throw string("See TODO in getSurroundingVerticesMap()");
	for(map<entHandle,vector<entHandle> >::iterator
			adjacentTets = adjacentTetsMap.begin();
			adjacentTets != adjacentTetsMap.end();
			adjacentTets++) {
		entHandle element = adjacentTets->first;
		set<entHandle> surroundingVerticesSet;
		for (vector<entHandle>::iterator adjacentTet =
				adjacentTets->second.begin();
				adjacentTet!=adjacentTets->second.end(); ++adjacentTet) {
			vector<entHandle> adjacentVertices =
					getAdjacentEntities(*adjacentTet, iBase_VERTEX);
			for (vector<entHandle>::iterator vertex =
					adjacentVertices.begin(); vertex!=adjacentVertices.end();
					++vertex) {
				surroundingVerticesSet.insert(*vertex);
			}
		}
		// Remove the vertices of the tet itself
		vector<entHandle> vertices =
				getAdjacentEntities(element, iBase_VERTEX);
//		vector<entHandle>& vec = adjacentTets->second;
		int initialSize = surroundingVerticesSet.size();
		for (vector<entHandle>::iterator vertex =
				vertices.begin(); vertex!=vertices.end(); ++vertex) {
			surroundingVerticesSet.erase(*vertex);
//			// Compact elements!=*vertex to front of vector and remove remaining
//			vec.erase(remove(vec.begin(), vec.end(), *vertex), vec.end());
		}
		vector<entHandle> surroundingVertices(
				surroundingVerticesSet.begin(), surroundingVerticesSet.end());
//		cout << element << endl;
//		cout << surroundingVerticesSet.size() << " " <<
//				initialSize << endl;
		if (surroundingVertices.size()!=initialSize-4)
			throw string("problem in getSurroundingVerticesMap()");
		surroundingVerticesMap[element] = surroundingVertices;
	}
	return surroundingVerticesMap;
}

map<entHandle,vector<entHandle> >
Mesh::getAdjacentsMap(int keyEntityType, int valueEntityType,
		int bridgeEntityType) {
	map<entHandle,vector<entHandle> > adjacentsMap;
	entHandle *ents3d = NULL;
	int ents3d_alloc = 0, ents3d_size;
	int ierr;
	iMesh_getEntities(meshInstance, rootEntitySet, keyEntityType,
			iMesh_ALL_TOPOLOGIES,
			&ents3d, &ents3d_alloc, &ents3d_size, &ierr);
	CHECK("Couldn't get entities");
	for(int i=0; i<ents3d_size; i++) {
		entHandle *entities = NULL;
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
		vector<entHandle> handles;
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

vect3d Mesh::getCoordinates(entHandle node, bool useMap) {
	vect3d coordinates(0.,0.,0.);
	double x,y,z;
	int ierr;
	iMesh_getVtxCoord(meshInstance, node, &x, &y, &z, &ierr);
	CHECK("Failure getting vertex coordinates");
	coordinates << x, y, z;
	return coordinates;
}

// TODO: eliminate code duplication with index-version
entHandle Mesh::findTet(vect3d oldPosition,
		vect3d position,
		entHandle adjacentTet, bool *tetFound, bool isTet) {
	entHandle tet;
	entHandle *entities = NULL;
	int entities_alloc=0, entities_size=0;
	int ierr;
	vector<entHandle> ents;
	vector<entHandle> faces;

	if (isTet) {
		entHandle faceCrossed = this->findFaceCrossed(
				adjacentTet, oldPosition, position);
		if (faceCrossed) {
			ents = adjacentTetsToFaceMap[faceCrossed];
//			ents = this->getAdjacentEntities(faceCrossed, iBase_REGION);
			for (int i=0; i<ents.size(); i++) {
				if (this->checkIfInTet(position, ents[i])) {
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
		if (Mesh::checkIfInTet(position, ents[i])) {
			tet = ents[i];
			// TODO: could throw error if tetFound is false rather than pass
			*tetFound = true;
			break;
		}
	}
	if (!*tetFound) {
		double pos[3];
		pos[0] = position[0];
		pos[1] = position[1];
		pos[2] = position[2];
		double pcoords[3], weights[3];
		vtkSmartPointer<vtkGenericCell> cell =
				vtkSmartPointer<vtkGenericCell>::New();
		vtkIdType cellId = vtkCellTree_ptr->FindCell(pos,0,
				cell, pcoords, weights);
//		cout << "In cell " << cellId << endl;
		if (cellId>=0) {
			tet = vtkToIMesh[int(cellId)];
			if (this->checkIfInTet(position, tet)) {
				*tetFound = true;
			} else {
//				throw;
			}
		} else {
//			throw int(OUTSIDE_DOMAIN);
		}
	}
	if (!*tetFound) {
//		tet = adjacentTet;
		tet = ents[0];
	}

//	int dimension=this->getEntityDimension(tet);
//	cout << tet << " " << *tetFound << " " << isTet << " " << dimension << endl;
	return tet;
}

int Mesh::findTet(vect3d oldPosition,
		vect3d position,
		int adjacentTetIndex, bool *tetFound, bool isTet) {
	int tetIndex=-1;
	vector<int> ents;
	vector<int> faces;

	*tetFound=false;
	if (isTet) {
		int vertexWithNegativeWeight=-1;
		this->checkIfInTet(position, adjacentTetIndex,
				&vertexWithNegativeWeight);
		int possibleNewRegionIndex;
		if (vertexWithNegativeWeight>=0) {
			possibleNewRegionIndex =
					this->regionsOppositeVertices[adjacentTetIndex][vertexWithNegativeWeight];
		} else {
			possibleNewRegionIndex = adjacentTetIndex;
		}
		if (this->checkIfInTet(position, possibleNewRegionIndex)) {
			tetIndex = possibleNewRegionIndex;
			*tetFound = true;
			return tetIndex;
		}
//		int faceCrossedIndex = this->findFaceCrossed(
//				adjacentTetIndex, oldPosition, position);
//		if (faceCrossedIndex>=0) {
//			ents = adjacentEntitiesVectors[iBase_FACE][faceCrossedIndex][iBase_REGION];
//			for (int i=0; i<ents.size(); i++) {
//				if (this->checkIfInTet(position, ents[i])) {
//					tetIndex = ents[i];
//					*tetFound = true;
//					return tetIndex;
//				}
//			}
//		}
		ents = regionsSurroundingRegions[adjacentTetIndex];
	} else {
		// TODO: should no longer get here
		cout << "findTet was passed non-region entity: " << adjacentTetIndex << endl;
		ents = regionsSurroundingRegions[0];
	}

	for (int i=0; i<ents.size(); i++) {
		if (Mesh::checkIfInTet(position, ents[i])) {
			tetIndex = ents[i];
			// TODO: could throw error if tetFound is false rather than pass
			*tetFound = true;
			break;
		}
	}
	if (!*tetFound) {
		double pos[3];
		pos[0] = position[0];
		pos[1] = position[1];
		pos[2] = position[2];
		double pcoords[3], weights[3];
		vtkSmartPointer<vtkGenericCell> cell =
				vtkSmartPointer<vtkGenericCell>::New();
		vtkIdType cellId = vtkCellTree_ptr->FindCell(pos,0,
				cell, pcoords, weights);
//		cout << "In cell " << cellId << endl;
		if (cellId>=0) {
			// TODO: assuming vtkIdType is int and ordering is same as in entitiesVectors
			if (this->checkIfInTet(position, cellId)) {
				*tetFound = true;
				tetIndex = cellId;
			} else {
//				throw;
			}
		} else {
//			throw int(OUTSIDE_DOMAIN);
		}
	}
	if (!*tetFound) {
//		tet = adjacentTet;
		tetIndex = ents[0];
	}

	if (tetIndex<0 || tetIndex>=entitiesVectors[iBase_REGION].size())
		throw string("tetIndex out of bounds in findTet()");
	return tetIndex;
}

entHandle Mesh::findStartingTet(vect3d const &position,
		vect3d const &velocity, entHandle vertex) {
	entHandle startingTet=NULL;
	// TODO: pertubation should be large enough that registers as only within tolerance of one tet
	// TODO: EXB drift wrong for electrons since uits are different
//	vect3d perturbedPosition =
//			position + sqrt(DELTA_LENGTH)*(velocity+VEXB)/(velocity+VEXB).norm();
	// TODO: standardize this
	vect3d perturbedPosition = position + (velocity+extern_VEXB)*SMALL_TIME;
	vector<entHandle> adjacentElements=getAdjacentEntities(vertex, iBase_REGION);
	int numberOfRegionsWithinTolerance=0;
	for (int i=0; i<adjacentElements.size(); i++) {
		if (checkIfInTet(perturbedPosition, adjacentElements[i])) {
			startingTet = adjacentElements[i];
			numberOfRegionsWithinTolerance++;
		}
	}
	if (numberOfRegionsWithinTolerance!=1)
		throw numberOfRegionsWithinTolerance;
	return startingTet;
}

bool Mesh::vertexLessThan(entHandle a, entHandle b) {
	vect3d aPos = this->getCoordinates(a);
	vect3d bPos = this->getCoordinates(b);
	return vect3dLessThan(aPos, bPos);
}

vector<entHandle> Mesh::getEntities(int dimension) {
	int ierr;
	entHandle *ents = NULL;
	int ents_alloc = 0, ents_size;
	iMesh_getEntities(meshInstance, rootEntitySet,
			dimension, iMesh_ALL_TOPOLOGIES,
			&ents, &ents_alloc, &ents_size, &ierr);
	CHECK("Couldn't get vertex entities");
	vector<entHandle> elements(ents_size);
	for (int i = 0; i < ents_size; i++) {
		elements[i] = ents[i];
	}
	if (dimension==0) {
		sort(elements.begin(), elements.end(), boost::bind(&Mesh::vertexLessThan, this, _1, _2));
	}
	if (ents) free(ents);
	ents_alloc = 0;
	return elements;
}

vector<entHandle> Mesh::getVertices(
		entHandle element) {
	return this->getAdjacentEntities(element, iBase_VERTEX);
}

vector<entHandle> Mesh::getAdjacentEntities(
		entHandle element, int dimension) {
	int ierr;
	entHandle *elements = NULL;
	int elements_alloc = 0, elements_size;

	iMesh_getEntAdj(meshInstance, element, dimension,
			&elements, &elements_alloc, &elements_size, &ierr);
	CHECK("Getting vertices adjacent to entity failed");

	vector<entHandle> elementHandles(elements_size);
	for (int i=0; i<elements_size; i++) {
		elementHandles[i] = elements[i];
	}
	if(elements) free (elements);
	elements_alloc = 0;

	return elementHandles;
}

vector<int> Mesh::getAdjacentEntitiesIndices(int entityIndex, int dimension,
		int adjacentsDimension) {
	entHandle entityHandle = entitiesVectors[dimension][entityIndex];
	vector<entHandle> adjacentEntities = this->getAdjacentEntities(
			entityHandle, adjacentsDimension);
	vector<int> adjacentEntitiesIndices(adjacentEntities.size());
	for (int i=0; i<adjacentEntities.size(); i++) {
		adjacentEntitiesIndices[i] = indicesOfEntities[adjacentEntities[i]];
	}

	return adjacentEntitiesIndices;
}

entHandle Mesh::getRandomVertex() {
	entHandle *ents0d = NULL;
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

int Mesh::getEntityDimension(entHandle entity) {
	int ierr;
	int dimension;
	iMesh_getEntType(meshInstance, entity, &dimension, &ierr);
	CHECK("Failure getting entity dimension");
	return dimension;
}

vector<vect3d> Mesh::getVertexVectors(entHandle entity,
		bool useMap) {
	int dimension = this->getEntityDimension(entity);
	int nVerts = dimension + 1;
	vector<vect3d> vertexVectors(nVerts);
	vector<entHandle> vertices(nVerts);
	if (dimension == 3) {
		vertices = adjacentVertsMap[entity];
	} else {
		vertices = this->getAdjacentEntities(entity, iBase_VERTEX);
	}
	if (vertexVectors.size()!=vertices.size())
		throw string("problem in getVertexVectors");
	if (useMap && dimension==3) {
		vertexVectors = vertexVectorsMap[entity];
	} else {
		for (int i=0; i<vertices.size(); i++) {
			vertexVectors[i] = Mesh::getCoordinates(vertices[i]);
		}
	}

	return vertexVectors;
}

vector<vect3d> Mesh::getVertexVectors(int index,
		int dimension) {
	int nVerts = dimension + 1;
	vector<vect3d> vertexVectors(nVerts);
	this->getVertexVectors(index, dimension, &vertexVectors);
	return vertexVectors;
}

void Mesh::getVertexVectors(int index,
		int dimension, vector<vect3d> *vertexVectors) {
	int nVerts = dimension + 1;
	vector<int> &vertices =
			adjacentEntitiesVectors[dimension][index][iBase_VERTEX];
	if (vertexVectors->size()!=vertices.size())
		throw string("problem in getVertexVector()");
	for (int i=0; i<vertices.size(); i++) {
		(*vertexVectors)[i] = verticesPositions[vertices[i]];
	}
}

bool Mesh::checkIfInTet(vect3d currentPosition,
		int elementIndex, int *nodeWithNegativeWeight) {
	vector<vect3d> vertexVectors = this->getVertexVectors(elementIndex, iBase_REGION);

	return this->checkIfInTet(currentPosition, vertexVectors, nodeWithNegativeWeight);
}

bool Mesh::checkIfInTet(vect3d currentPosition,
		entHandle element) {
	vector<vect3d> vertexVectors = this->getVertexVectors(element);

	return this->checkIfInTet(currentPosition, vertexVectors);
}

bool Mesh::checkIfInTet(vect3d currentPosition,
		vector<vect3d> vertexVectors, int *nodeWithNegativeWeight) {
//	double tetVolume = this->getTetVolume(vertexVectors);
//	if (tetVolume<VOLUME_TOLERANCE)
//		return false;
//	vector<double> subVolumes = this->getTetSubVolumes(currentPosition,
//			vertexVectors);
//	double sumSubVolumes =
//			accumulate(subVolumes.begin(),subVolumes.end(),0.);
//	return (fabs(sumSubVolumes-tetVolume)<VOLUME_TOLERANCE);
	Eigen::Vector4d linearBasisFunctions = this->evaluateLinearBasisFunctions(
			currentPosition, vertexVectors);
	bool inElement=true;
	for (int i=0; i<linearBasisFunctions.rows(); i++) {
		if (linearBasisFunctions[i]<0.-VOLUME_TOLERANCE) {
			inElement=false;
			if (nodeWithNegativeWeight)
				*nodeWithNegativeWeight = i;
			break;
		}
	}
	return inElement;
}

bool Mesh::checkIfIntersectsTriangle(vect3d previousPosition,
		vect3d currentPosition,
		vector<vect3d> vertexVectors) {
	vect3d intersectionPoint;
	return this->checkIfIntersectsTriangle(previousPosition, currentPosition,
			vertexVectors, &intersectionPoint);
}

bool Mesh::checkIfIntersectsTriangle(vect3d previousPosition,
		vect3d currentPosition,
		vector<vect3d> vertexVectors, vect3d *intersectionPoint) {
	vector<vect3d> ray(2);
	ray[0] = previousPosition;
	ray[1] = currentPosition;
	int intersectionStatus =
			intersect_SegmentTriangle(ray, vertexVectors, intersectionPoint);
//	int intersectionStatus =
//			intersect_RayTriangle(ray, vertexVectors, &intersectionPoint);
	if (intersectionStatus == 1) {
//		vect3d relativePosition = currentPosition - previousPosition;
//		// Check whether intersection point is between input points
//		if (relativePosition.dot(intersectionPoint-currentPosition)*
//				relativePosition.dot(intersectionPoint-previousPosition) <= 0.) {
			return true;
//		}
	}
	return false;
}

bool Mesh::checkIfRayIntersectsTriangle(vect3d previousPosition,
		vect3d currentPosition,
		vector<vect3d> vertexVectors) {
	vector<vect3d> ray(2);
	ray[0] = previousPosition;
	ray[1] = currentPosition;
	vect3d intersectionPoint;
	int intersectionStatus =
			intersect_RayTriangle(ray, vertexVectors, &intersectionPoint);
	if (intersectionStatus == 1) {
			return true;
	}
	return false;
}

bool Mesh::checkIfLineIntersectsTriangle(vect3d previousPosition,
		vect3d currentPosition,
		vector<vect3d> vertexVectors) {
	bool forwardIntersection = this->checkIfRayIntersectsTriangle(
			previousPosition, currentPosition, vertexVectors);
	bool backwardIntersection = this->checkIfRayIntersectsTriangle(
			currentPosition, previousPosition, vertexVectors);
	return (forwardIntersection || backwardIntersection);
}

//vector<vect3d> Mesh::getEdgeVectors(
//		vector<vect3d> vertexVectors) {
//	int nVerts = vertexVectors.size();
//	int nEdges=nVerts-1;
//	vector<vect3d> edgeVectors(nEdges);
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


entHandle Mesh::findFaceCrossed(entHandle previousElement,
		vect3d previousPosition, vect3d currentPosition) {
	entHandle faceCrossed=NULL;
	vector<entHandle> adjacentFaces =
			adjacentFacesMap[previousElement];
//	vector<entHandle> adjacentFaces =
//			getAdjacentEntities(previousElement,iBase_FACE);

	for (int i=0; i<adjacentFaces.size(); i++) {
		vector<vect3d> vertexVectors =
				this->getVertexVectors(adjacentFaces[i], true);
		bool intersectsTriangle = this->checkIfIntersectsTriangle(previousPosition,
				currentPosition, vertexVectors);
		if (intersectsTriangle)
			faceCrossed = adjacentFaces[i];
	}

	return faceCrossed;
}

int Mesh::findFaceCrossed(int previousElementIndex,
		vect3d previousPosition, vect3d currentPosition) {
	// TODO: which linearBasisFunction is negative tells you which face was
	//       crossed (face opposite vertex with negative weight)
	int faceCrossedIndex=-1;
	vector<int> &adjacentFaces =
			adjacentEntitiesVectors[iBase_REGION][previousElementIndex][iBase_FACE];
	for (int i=0; i<adjacentFaces.size(); i++) {
		vector<vect3d> vertexVectors =
				this->getVertexVectors(adjacentFaces[i], iBase_FACE);
		bool intersectsTriangle = this->checkIfIntersectsTriangle(previousPosition,
				currentPosition, vertexVectors);
		if (intersectsTriangle)
			faceCrossedIndex = adjacentFaces[i];
	}

	return faceCrossedIndex;
}

int Mesh::findBoundaryFaceCrossed(int previousElementIndex,
		vect3d previousPosition, vect3d currentPosition,
		Field<int>& faceTypeField, CodeField& vertexTypeField) {
	vect3d centroid(0.,0.,0.);
	vector<vect3d> vVs =
			this->getVertexVectors(previousElementIndex,iBase_REGION);
	centroid = (vVs[0]+vVs[1]+vVs[2]+vVs[3])/4.;
	// TODO: do something other than this hack to bring previousPosition away from node
	previousPosition += 0.00001*centroid;
	if (!this->checkIfInTet(previousPosition, previousElementIndex)) {
		bool tetFound=false;
		int actualPreviousElementIndex =
				this->findTet(previousPosition, previousPosition, previousElementIndex, &tetFound);
		if (tetFound) {
//			cout << "Element index corrected: " << previousElementIndex << " -> " <<
//					actualPreviousElementIndex << endl;
			previousElementIndex = actualPreviousElementIndex;
		}
	}
	int faceCrossedIndex=-1;
	int numberVisited=0;
	int previousFaceIndex=-1;
	int faceIndex=-1;
	// First try to traverse elements to find boundary face crossed
	int elementIndex = previousElementIndex;
	while (faceCrossedIndex<0 && numberVisited<100) {
		vector<int> &faces =
				adjacentEntitiesVectors[iBase_REGION][elementIndex][iBase_FACE];
		previousFaceIndex = faceIndex;
		for (int j=0; j<faces.size(); j++) {
			if (faces[j]!=previousFaceIndex) {
				vector<vect3d> vertexVectors =
						this->getVertexVectors(faces[j], iBase_FACE);
				bool intersectsTriangle = this->checkIfIntersectsTriangle(previousPosition,
						currentPosition, vertexVectors);
				// TODO: find out why reverse intersection not always same
				//       Probably when previousPosition is a vertex of the tet
				bool reverseIntersectsTriangle = this->checkIfIntersectsTriangle(
						currentPosition, previousPosition, vertexVectors);
//				if (intersectsTriangle!=reverseIntersectsTriangle)
//					cout << "Reversal different: " << previousElementIndex << " " <<
//					intersectsTriangle << " " <<
//						reverseIntersectsTriangle << " " << faces[j] << " " <<
//						currentPosition.transpose() << " " << previousPosition.transpose() << endl;
//				vector<vect3d> vVs =
//						this->getVertexVectors(elementIndex,iBase_REGION);
//				centroid = (vVs[0]+vVs[1]+vVs[2]+vVs[3])/4.;
//				if (intersectsTriangle || reverseIntersectsTriangle)
//					cout << "traversed face: " << previousElementIndex << " " <<
//					elementIndex << " " << centroid.transpose() << endl;
				if (intersectsTriangle || reverseIntersectsTriangle) {
					faceIndex = faces[j];
					if (faceTypeField[faces[j]]>0) {
						faceCrossedIndex = faces[j];
					} else {
						vector<int> &elements =
								adjacentEntitiesVectors[iBase_FACE][faces[j]][iBase_REGION];
						for (int i=0; i<elements.size(); i++) {
							// TODO: do this more transparently?
							if (elements[i]!=elementIndex) {
								elementIndex = elements[i];
								break;
							}
						}
					}
				}
			}
		}
		numberVisited++;
	}
	vVs = this->getVertexVectors(previousElementIndex,iBase_REGION);
	centroid = (vVs[0]+vVs[1]+vVs[2]+vVs[3])/4.;
	if (faceCrossedIndex<0) {
//		cout << "Boundary face fall-back for: " << previousElementIndex << " " <<
//			currentPosition.transpose() << " " << previousPosition.transpose() << endl;
		vector<int> facesCrossed;
		for (int j=0; j<entitiesVectors[iBase_FACE].size(); j++) {
			if (faceTypeField[j]>0) {
				vector<vect3d> vertexVectors =
						this->getVertexVectors(j, iBase_FACE);
				bool intersectsTriangle = this->checkIfLineIntersectsTriangle(centroid,
						currentPosition, vertexVectors);
				if (intersectsTriangle)
					facesCrossed.push_back(j);
			}
		}
		// Find closest boundary face intersected by line
		// TODO: don't hard-code large number
		double smallestDistanceBetweenCentroids = 1e10;
		for (int j=0; j<facesCrossed.size(); j++) {
			vVs = this->getVertexVectors(facesCrossed[j],iBase_FACE);
			vect3d faceCentroid = (vVs[0]+vVs[1]+vVs[2])/3.;
			double distanceBetweenCentroids = (faceCentroid-centroid).norm();
			if (distanceBetweenCentroids < smallestDistanceBetweenCentroids) {
				smallestDistanceBetweenCentroids = distanceBetweenCentroids;
				faceCrossedIndex = facesCrossed[j];
			}
		}
	}
	if (faceCrossedIndex<0) {
		cout << "Boundary face failure for: " << previousElementIndex << " " <<
			currentPosition.transpose() << " " << previousPosition.transpose() << endl;
	}
	return faceCrossedIndex;
}

iBase_TagHandle Mesh::getTagHandle(string tagName) {
	int ierr;
	iBase_TagHandle tag=0;
	iMesh_getTagHandle(meshInstance, tagName.c_str(),
			&tag, &ierr, tagName.length());
	return tag;
}

iBase_TagHandle Mesh::createTagHandle(string tagName, int size, int type) {
	int ierr;
	iBase_TagHandle tag=0;
	iMesh_createTag(meshInstance, tagName.c_str(),
			size, type, &tag, &ierr, (int)tagName.length());
	return tag;
}

vect3d Mesh::getSurfaceVector(entHandle face,
		vect3d point) {
	vector<vect3d> vertexVectors =
			this->getVertexVectors(face);
	if (vertexVectors.size()!=3)
		throw string("problem in getSurfaceVector");
	vector<vect3d> edgeVectors(vertexVectors.size()-1);

	edgeVectors[0] = vertexVectors[1]-vertexVectors[0];
	edgeVectors[1] = vertexVectors[2]-vertexVectors[0];

	vect3d surfaceVector = edgeVectors[0].cross(edgeVectors[1])/2.;
	vect3d referenceVector = point - vertexVectors[0];
	if (point==vect3d(0.,0.,0.) ||
			fabs(referenceVector.dot(surfaceVector))<LENGTH_TOLERANCE) {
		vector<entHandle> tets =
				this->getAdjacentEntities(face,iBase_REGION);
		// TODO: handle interior faces with two adjacent tets?
		vector<vect3d> tetVertexVectors =
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
		if (nPoints!=1)
			throw string("problem in getSurfaceVector()");
	}
	if (point==vect3d(0.,0.,0.))
		throw string("problem in getSurfaceVector().");
	referenceVector = point - vertexVectors[0];

	if (referenceVector.dot(surfaceVector)<0)
		surfaceVector *= -1.;

	return surfaceVector;
}

vect3d Mesh::getNormalVector(entHandle face,
		vect3d point) {
	vect3d surfaceVector = this->getSurfaceVector(face, point);
	return surfaceVector/surfaceVector.norm();
}

vect3d Mesh::getVertexNormalVector(entHandle vertex,
		Field<int> faceTypeField) {
	vector<entHandle> faces =
			this->getAdjacentEntities(vertex, iBase_FACE);
	// TODO: could pre-allocate normalVectors
	vector<vect3d> normalVectors;
	vector<entHandle>::iterator faceIter = faces.begin();
	for (; faceIter!=faces.end(); ++faceIter) {
		int faceType = faceTypeField.getField(*faceIter);
		if (faceType>0) {
			vect3d normalVector =
					this->getNormalVector(*faceIter);
			normalVectors.push_back(normalVector);
		}
	}
	vect3d normalVector(0.,0.,0.);
	vector<vect3d>::iterator normalIter =
			normalVectors.begin();
	for (; normalIter!=normalVectors.end(); ++normalIter) {
		normalVector += *normalIter;
	}
	normalVector /= normalVector.norm();
	vect3d testCoord = this->getCoordinates(vertex);
	return normalVector;
}

vector<entHandle> Mesh::getSuperCellFaces(
		entHandle vertex) {
	int ierr;
	entHandle *faces = NULL;
	int faces_alloc = 0, faces_size;
	vector<entHandle> superCellFaces;
	iBase_TagHandle code_tag;

	// TODO: cleaner way than hard-coding cell_code here?
	iMesh_getTagHandle(meshInstance, "cell_code", &code_tag,
			&ierr, 9);
	CHECK("Failure getting cell_code handle");

	iMesh_getEnt2ndAdj(meshInstance, vertex, iBase_REGION,
			iBase_FACE, &faces, &faces_alloc,
			&faces_size, &ierr);
	CHECK("Failure in getEnt2ndAdj");
	for (int i = 0; i < faces_size; i++) {
		entHandle *vertices = NULL;
		int vertices_alloc = 0, vertices_size;
		bool onSuperCellSurface = true;
		int cell_code = 0;

		iMesh_getIntData(meshInstance, faces[i], code_tag,
				&cell_code, &ierr);
		CHECK("Failure getting cell_code value");

		iMesh_getEntAdj(meshInstance, faces[i], iBase_VERTEX,
				&vertices, &vertices_alloc,
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

double Mesh::getTetVolume(vect3d point, entHandle face) {
	vector<vect3d> vertexVectors =
			this->getVertexVectors(face);
	vertexVectors.push_back(point);

	return getTetVolume(vertexVectors);
}

double Mesh::getTetVolume(vector<vect3d> vertexVectors) {
	int nEdges=3;
	vector<vect3d> edgeVectors(nEdges);

	int nVertices=4;
	if (vertexVectors.size()!=nVertices)
		throw string("problem in getTetVolume");

	edgeVectors[0] = vertexVectors[1]-vertexVectors[0];
	edgeVectors[1] = vertexVectors[2]-vertexVectors[0];
	edgeVectors[2] = vertexVectors[3]-vertexVectors[0];

	return fabs(
			edgeVectors[2].dot( edgeVectors[0].cross(edgeVectors[1]) )
			)/6.;
}

vector<double> Mesh::getTetSubVolumes(vect3d point,
		vector<vect3d> vertexVectors) {
	int nVertices=4;
	if (vertexVectors.size()!=nVertices)
		throw string("problem in getTetSubVolumes");
	vector<double> subVolumes(nVertices);
	for (int i=0; i<vertexVectors.size(); i++) {
		vect3d tmpVertex = vertexVectors[i];
		vertexVectors[i] = point;
		double volume = this->getTetVolume(vertexVectors);
		subVolumes[i] = volume;
		vertexVectors[i] = tmpVertex;
	}
	return subVolumes;
}

vector<double> Mesh::getVertexWeights(vect3d point,
		vector<vect3d> vertexVectors) {
	vector<double> subVolumes = this->getTetSubVolumes(point,
			vertexVectors);
	double totalVolume = this->getTetVolume(vertexVectors);
	for(int i=0; i<subVolumes.size(); i++)
	    subVolumes[i] /= totalVolume;

	return subVolumes;
}

Eigen::Vector4d Mesh::evaluateLinearBasisFunctions(vect3d position,
		entHandle element) {
	Eigen::Vector4d basisFunctions;
	if (element==previousCoordsToBasisElement && element!=NULL) {
		Eigen::Vector4d paddedPosition(1.,position[0],position[1],position[2]);
		basisFunctions = previousCoordsToBasis*paddedPosition;
	} else {
		vector<vect3d> vVs = this->getVertexVectors(element);
		basisFunctions = this->evaluateLinearBasisFunctions(position,vVs);
		// TODO: need to set this after evalLinBasFuncs call since it
		//       sets it to NULL, but not clean code
		previousCoordsToBasisElement = element;
	}
	return basisFunctions;
}

Eigen::Vector4d Mesh::evaluateLinearBasisFunctions(vect3d position,
		vector<vect3d> vVs) {
	// TODO: should replace 4 here with unified number across functions
	if (vVs.size()!=4)
		throw string("problem in evaluateLinearBasisFunctions");
	Eigen::Vector4d basisFunctions;
	Eigen::Matrix4d coordsToBasis = this->calculatePositionToBasesMatrix(vVs);
	Eigen::Vector4d paddedPosition(1.,position[0],position[1],position[2]);
	basisFunctions = coordsToBasis*paddedPosition;
	// TODO: have to set handle to NULL to prevent changing of matrix
	//       without changing handle, but not clean code...
	previousCoordsToBasisElement = NULL;
	previousCoordsToBasis = coordsToBasis;

	return basisFunctions;
}

Eigen::Vector4d Mesh::evaluateLinearBasisFunctions(vect3d position,
		int regionIndex) {
	Eigen::Vector4d basisFunctions;
	this->evaluateLinearBasisFunctions(position, regionIndex,
			&basisFunctions);
	return basisFunctions;
}

void Mesh::evaluateLinearBasisFunctions(const vect3d &position,
		int regionIndex, Eigen::Vector4d *basisFunctions) {
	Eigen::Vector4d paddedPosition(1.,position[0],position[1],position[2]);
	*basisFunctions = positionsToBases[regionIndex]*paddedPosition;
}

void Mesh::evaluateLinearBasisFunctionDerivatives(const vect3d &position,
		int regionIndex,
		Eigen::Matrix<double,NDIM,NDIM+1> *basisFunctionDerivatives) {
	Eigen::Vector4d basisFunctions;
	this->evaluateLinearBasisFunctions(position, regionIndex, &basisFunctions);
	for (int j=0; j<NDIM; j++) {
		vect3d perturbedPosition = position +
				vect3d::Unit(j)*DELTA_LENGTH;
		Eigen::Vector4d perturbedBasisFunctions;
		this->evaluateLinearBasisFunctions(perturbedPosition,
				regionIndex, &perturbedBasisFunctions);
		for (int i=0;i<basisFunctions.rows();i++) {
			(*basisFunctionDerivatives)(j,i) = (perturbedBasisFunctions[i] -
					basisFunctions[i])/DELTA_LENGTH;
		}
	}
}

Eigen::Matrix4d Mesh::calculatePositionToBasesMatrix(vector<vect3d> vVs) {
	// TODO: should replace 4 here with unified number across functions
	if (vVs.size()!=4)
		throw string("problem in calculatePositionToBasesMatrix");
	Eigen::Matrix4d basisToCoords;
	basisToCoords <<
			1.,			1.,			1.,			1.,
			vVs[0][0],	vVs[1][0],	vVs[2][0],	vVs[3][0],
			vVs[0][1],	vVs[1][1],	vVs[2][1],	vVs[3][1],
			vVs[0][2],	vVs[1][2],	vVs[2][2],	vVs[3][2];
	Eigen::Matrix4d coordsToBasis = basisToCoords.inverse();
	return coordsToBasis;
}

Eigen::VectorXd Mesh::evaluateQuadraticErrorBases(
		Eigen::Vector4d linearBasisFunctions) {
	if (linearBasisFunctions.rows()!=4)
		throw string("problem in evaluateQuadraticErrorBases");
	Eigen::VectorXd quadraticBasisFunctions(6);
	this->evaluateQuadraticErrorBases(linearBasisFunctions,
			&quadraticBasisFunctions);
	return quadraticBasisFunctions;
}

void Mesh::evaluateQuadraticErrorBases(
		Eigen::Vector4d linearBasisFunctions,
		Eigen::VectorXd *quadraticBasisFunctions) {
	if (linearBasisFunctions.rows()!=4)
		throw string("problem in evaluateQuadraticErrorBases");
	(*quadraticBasisFunctions)[0] =
			linearBasisFunctions[0]*linearBasisFunctions[1];
	(*quadraticBasisFunctions)[1] =
			linearBasisFunctions[0]*linearBasisFunctions[2];
	(*quadraticBasisFunctions)[2] =
			linearBasisFunctions[0]*linearBasisFunctions[3];
	(*quadraticBasisFunctions)[3] =
			linearBasisFunctions[1]*linearBasisFunctions[2];
	(*quadraticBasisFunctions)[4] =
			linearBasisFunctions[1]*linearBasisFunctions[3];
	(*quadraticBasisFunctions)[5] =
			linearBasisFunctions[2]*linearBasisFunctions[3];
}

Eigen::VectorXd Mesh::evaluateCubicErrorBases(
		Eigen::Vector4d linearBasisFunctions) {
	Eigen::VectorXd cubicBasisFunctions(16);
	this->evaluateCubicErrorBases(linearBasisFunctions,
			&cubicBasisFunctions);
	return cubicBasisFunctions;
}

void Mesh::evaluateCubicErrorBases(
		Eigen::Vector4d linearBasisFunctions,
		Eigen::VectorXd *cubicBasisFunctions) {
//	Eigen::VectorXd cubicBasisFunctions(22);
	// TODO: replace below with some sort of loop?
	(*cubicBasisFunctions)[0] = linearBasisFunctions[0]*
			linearBasisFunctions[0]*linearBasisFunctions[1];
	(*cubicBasisFunctions)[1] = linearBasisFunctions[0]*
			linearBasisFunctions[0]*linearBasisFunctions[2];
	(*cubicBasisFunctions)[2] = linearBasisFunctions[0]*
			linearBasisFunctions[0]*linearBasisFunctions[3];
	(*cubicBasisFunctions)[3] = linearBasisFunctions[0]*
			linearBasisFunctions[1]*linearBasisFunctions[2];
	(*cubicBasisFunctions)[4] = linearBasisFunctions[0]*
			linearBasisFunctions[1]*linearBasisFunctions[3];
	(*cubicBasisFunctions)[5] = linearBasisFunctions[0]*
			linearBasisFunctions[2]*linearBasisFunctions[3];

	(*cubicBasisFunctions)[6] = linearBasisFunctions[1]*
			linearBasisFunctions[0]*linearBasisFunctions[1];
	(*cubicBasisFunctions)[7] = linearBasisFunctions[1]*
			linearBasisFunctions[1]*linearBasisFunctions[2];
	(*cubicBasisFunctions)[8] = linearBasisFunctions[1]*
			linearBasisFunctions[1]*linearBasisFunctions[3];
	(*cubicBasisFunctions)[9] = linearBasisFunctions[1]*
			linearBasisFunctions[2]*linearBasisFunctions[3];

	(*cubicBasisFunctions)[10] = linearBasisFunctions[2]*
			linearBasisFunctions[0]*linearBasisFunctions[2];
	(*cubicBasisFunctions)[11] = linearBasisFunctions[2]*
			linearBasisFunctions[1]*linearBasisFunctions[2];
	(*cubicBasisFunctions)[12] = linearBasisFunctions[2]*
			linearBasisFunctions[2]*linearBasisFunctions[3];

	(*cubicBasisFunctions)[13] = linearBasisFunctions[3]*
			linearBasisFunctions[0]*linearBasisFunctions[3];
	(*cubicBasisFunctions)[14] = linearBasisFunctions[3]*
			linearBasisFunctions[1]*linearBasisFunctions[3];
	(*cubicBasisFunctions)[15] = linearBasisFunctions[3]*
			linearBasisFunctions[2]*linearBasisFunctions[3];

//	cubicBasisFunctions[16] =
//			linearBasisFunctions[0]*linearBasisFunctions[1];
//	cubicBasisFunctions[17] =
//			linearBasisFunctions[0]*linearBasisFunctions[2];
//	cubicBasisFunctions[18] =
//			linearBasisFunctions[0]*linearBasisFunctions[3];
//	cubicBasisFunctions[19] =
//			linearBasisFunctions[1]*linearBasisFunctions[2];
//	cubicBasisFunctions[20] =
//			linearBasisFunctions[1]*linearBasisFunctions[3];
//	cubicBasisFunctions[21] =
//			linearBasisFunctions[2]*linearBasisFunctions[3];

}

double Mesh::minimumBasisFunction(const vect3d &position, int &regionIndex) {
	Eigen::Vector4d bF(0.,0.,0.,0.);
	this->evaluateLinearBasisFunctions(position, regionIndex, &bF);
	return min(min(min(bF[0],bF[1]),bF[2]),bF[3]);
}

vtkSmartPointer<vtkUnstructuredGrid> Mesh::createVtkMesh() {
	vtkSmartPointer<vtkUnstructuredGrid> mesh =
			vtkSmartPointer<vtkUnstructuredGrid>::New();

	vtkSmartPointer<vtkPoints> points =
			vtkSmartPointer<vtkPoints>::New();
	vector<entHandle> vertices =
			this->getEntities(iBase_VERTEX);
	points->SetNumberOfPoints(vertices.size());
	map<entHandle,int> indexOfVertex;
	for(int i=0; i<vertices.size(); i++) {
		vect3d vV = this->getCoordinates(vertices[i]);
		points->InsertPoint(i, vV[0], vV[1], vV[2]);
		indexOfVertex[vertices[i]] = i;
	}
	mesh->SetPoints(points);

	vtkSmartPointer<vtkCellArray> cellArray =
			vtkSmartPointer<vtkCellArray>::New();
	vtkToIMesh = this->getEntities(iBase_REGION);
	vtkSmartPointer<vtkTetra> tetra =
			vtkSmartPointer<vtkTetra>::New();
	for (int i=0; i<vtkToIMesh.size(); i++) {
		vector<entHandle> verts =
				this->getAdjacentEntities(vtkToIMesh[i],
						iBase_VERTEX);
		tetra->GetPointIds()->SetId(0, indexOfVertex[verts[0]]);
		tetra->GetPointIds()->SetId(1, indexOfVertex[verts[1]]);
		tetra->GetPointIds()->SetId(2, indexOfVertex[verts[2]]);
		tetra->GetPointIds()->SetId(3, indexOfVertex[verts[3]]);
		cellArray->InsertNextCell(tetra);
	}
	mesh->SetCells(tetra->GetCellType(), cellArray);
	mesh->Update();

	return mesh;
}
