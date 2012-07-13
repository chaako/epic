/*
 * Mesh.h
 *
 *  Created on: Feb 10, 2012
 *      Author: chaako
 */

#ifndef MESH_H_
#define MESH_H_

#include "typesAndDefinitions.h"

//#include <map>
//#include <vector>
//
//#include <vtkVersion.h>
//#include <vtkSmartPointer.h>
//#include <vtkCellArray.h>
//#include <vtkTetra.h>
//#include <vtkUnstructuredGrid.h>
//#include <vtkPoints.h>
//#include <vtkCellType.h>
//#include <vtkDataSetMapper.h>


template <class T> class Field;

class Mesh {
public:
	Mesh(string inputMeshFile);
	virtual ~Mesh();

	void printElementNumbers();
	void save(string outputMeshFile);

	vect3d getCoordinates(entHandle node, bool useMap=false);
	entHandle findTet(vect3d oldPosition, vect3d position,
			entHandle adjacentTet, bool *tetFound, bool isTet=true);
	int findTet(vect3d oldPosition, vect3d position,
			int adjacentTetIndex, bool *tetFound, bool isTet=true);
//	entHandle findStartingTet(vect3d const &position,
//			vect3d const &velocity, entHandle vertex);
	vector<entHandle> getVertices(entHandle element);
	vector<entHandle> getAdjacentEntities(entHandle element,
			int dimension);
	vector<int> getAdjacentEntitiesIndices(int entityIndex, int dimension,
			int adjacentsDimension);
	vector<entHandle> getEntities(int dimension);
	entHandle getRandomVertex();
	vector<entHandle> getFaces(entHandle element);
	map<entHandle,vector<entHandle> >
	getSurroundingVerticesMap();
	map<entHandle,vector<entHandle> >
	getAdjacentsMap(int keyEntityType, int valueEntityType,
			int bridgeEntityType=-1);
	int getEntityDimension(entHandle entity);
	vector<vect3d> getVertexVectors(entHandle entity,
			bool useMap=true);
	vector<vect3d> getVertexVectors(int index, int dimension);
	void getVertexVectors(int index, int dimension,
			vector<vect3d> *vertexVectors);
	bool checkIfInTet(vect3d currentPosition,
			entHandle element);
	bool checkIfInTet(vect3d currentPosition,
			int elementIndex, int *nodeWithNegativeWeight=NULL);
	bool checkIfInTet(vect3d currentPosition,
			vector<vect3d> vertexVectors, int *nodeWithNegativeWeight=NULL);
	bool checkIfIntersectsTriangle(vect3d previousPosition,
			vect3d currentPosition,
			vector<vect3d> vertexVectors);
	entHandle findFaceCrossed(entHandle previousElement,
			vect3d previousPosition, vect3d currentPosition);
	int findFaceCrossed(int previousElementIndex,
			vect3d previousPosition, vect3d currentPosition);
	iBase_TagHandle getTagHandle(string tagName);
	iBase_TagHandle createTagHandle(string tagName, int size, int type);
	vect3d getSurfaceVector(entHandle face,
			vect3d point=vect3d(0.,0.,0.));
	vect3d getNormalVector(entHandle face,
			vect3d point=vect3d(0.,0.,0.));
	vect3d getVertexNormalVector(entHandle vertex,
			Field<int> faceTypeField);
	vector<entHandle> getSuperCellFaces(entHandle vertex);
	double getTetVolume(vect3d point, entHandle face);
	double getTetVolume(vector<vect3d> vertexVectors);
	vector<double> getTetSubVolumes(vect3d point,
			vector<vect3d> vertexVectors) ;
	vector<double> getVertexWeights(vect3d point,
			vector<vect3d> vertexVectors);
	Eigen::Vector4d evaluateLinearBasisFunctions(vect3d position,
			entHandle element);
	Eigen::Vector4d evaluateLinearBasisFunctions(vect3d position,
			vector<vect3d> vVs);
	Eigen::Vector4d evaluateLinearBasisFunctions(vect3d position,
			int regionIndex);
	void evaluateLinearBasisFunctions(const vect3d &position,
			int regionIndex, Eigen::Vector4d *basisFunctions);
	Eigen::Matrix4d calculatePositionToBasesMatrix(vector<vect3d> vVs);
	Eigen::VectorXd evaluateQuadraticErrorBases(
			Eigen::Vector4d linearBasisFunctions);
	void evaluateQuadraticErrorBases(
			Eigen::Vector4d linearBasisFunctions,
			Eigen::VectorXd *quadraticBasisFunctions);
	Eigen::VectorXd evaluateCubicErrorBases(
			Eigen::Vector4d linearBasisFunctions);
	void evaluateCubicErrorBases(
			Eigen::Vector4d linearBasisFunctions,
			Eigen::VectorXd *cubicBasisFunctions);

	vtkSmartPointer<vtkUnstructuredGrid> createVtkMesh();

	vtkSmartPointer<vtkUnstructuredGrid> vtkMesh_ptr;
	vtkSmartPointer<vtkCellTreeLocator> vtkCellTree_ptr;
	vector<entHandle> vtkToIMesh;

	iMesh_Instance meshInstance;
	// TODO: rename vtkInputMesh something like inputMeshIsVtk
	bool vtkInputMesh;
	iBase_EntitySetHandle rootEntitySet;
//	iBase_EntitySetHandle vextexEntitySet;
//	iBase_EntitySetHandle surfaceEntitySet;
//	iBase_EntitySetHandle volumeEntitySet;

	vector<vector<entHandle> > entitiesVectors;

	// adjacentEntities[dimension][index][adjacentDimension][adjacentIndices]
	vector<vector<vector<vector<int> > > > adjacentEntitiesVectors;

	vector<vector<int> > verticesSurroundingRegions;
	vector<vector<int> > regionsSurroundingRegions;
	// regionsOppositeVertices[regionIndex][vertexIndex]
	vector<vector<int> > regionsOppositeVertices;

	vector<vect3d> verticesPositions;
	vector<Eigen::Matrix4d> positionsToBases;

//	vector<entHandle> allVertices;
//	vector<entHandle> allFaces;
//	vector<entHandle> allElements;
//
//	map<entHandle,int> indexOfVertices;
//	map<entHandle,int> indexOfFaces;
//	map<entHandle,int> indexOfElements;

	map<entHandle,int> indicesOfEntities;
	map<entHandle,int> dimensionsOfEntities;

	entHandle previousCoordsToBasisElement;
	Eigen::Matrix4d previousCoordsToBasis;

	map<entHandle,vector<entHandle> > adjacentTetsMap;
	map<entHandle,vector<entHandle> > adjacentFacesMap;
	map<entHandle,vector<entHandle> > adjacentVertsMap;
	map<entHandle,vector<entHandle> > surroundingVertsMap;
	map<entHandle,vector<vect3d> > vertexVectorsMap;
	map<entHandle,vector<entHandle> > adjacentTetsToFaceMap;
};

#endif /* MESH_H_ */
