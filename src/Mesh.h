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
	Mesh(std::string inputMeshFile);
	virtual ~Mesh();

	void printElementNumbers();
	void save(std::string outputMeshFile);

	vect3d getCoordinates(entHandle node, bool useMap=false);
	entHandle findTet(vect3d oldPosition, vect3d position,
			entHandle adjacentTet, bool *tetFound, bool isTet=true);
//	entHandle findStartingTet(vect3d const &position,
//			vect3d const &velocity, entHandle vertex);
	std::vector<entHandle> getVertices(entHandle element);
	std::vector<entHandle> getAdjacentEntities(entHandle element,
			int dimension);
	std::vector<entHandle> getEntities(int dimension);
	entHandle getRandomVertex();
	std::vector<entHandle> getFaces(entHandle element);
	std::map<entHandle,std::vector<entHandle> >
	getSurroundingVerticesMap();
	std::map<entHandle,std::vector<entHandle> >
	getAdjacentsMap(int keyEntityType, int valueEntityType,
			int bridgeEntityType=-1);
	int getEntityDimension(entHandle entity);
	std::vector<vect3d> getVertexVectors(entHandle entity,
			bool useMap=true);
	bool checkIfInTet(vect3d currentPosition,
			entHandle element);
	bool checkIfInTet(vect3d currentPosition,
			std::vector<vect3d> vertexVectors);
	bool checkIfIntersectsTriangle(vect3d previousPosition,
			vect3d currentPosition,
			std::vector<vect3d> vertexVectors);
	entHandle findFaceCrossed(entHandle previousElement,
			vect3d previousPosition, vect3d currentPosition);
	iBase_TagHandle getTagHandle(std::string tagName);
	iBase_TagHandle createTagHandle(std::string tagName, int size, int type);
	vect3d getSurfaceVector(entHandle face,
			vect3d point=vect3d(0.,0.,0.));
	vect3d getNormalVector(entHandle face,
			vect3d point=vect3d(0.,0.,0.));
	vect3d getVertexNormalVector(entHandle vertex,
			Field<int> faceTypeField);
	std::vector<entHandle> getSuperCellFaces(entHandle vertex);
	double getTetVolume(vect3d point, entHandle face);
	double getTetVolume(std::vector<vect3d> vertexVectors);
	std::vector<double> getTetSubVolumes(vect3d point,
			std::vector<vect3d> vertexVectors) ;
	std::vector<double> getVertexWeights(vect3d point,
			std::vector<vect3d> vertexVectors);
	Eigen::VectorXd getErrorCoefficients(vect3d position,
			entHandle element, int interpolationOrder=1);
	Eigen::VectorXd getErrorCoefficients(vect3d position,
			std::vector<vect3d> vVs, int interpolationOrder=1);
	Eigen::Vector4d evaluateLinearBasisFunctions(vect3d position,
			entHandle element);
	Eigen::Vector4d evaluateLinearBasisFunctions(vect3d position,
			std::vector<vect3d> vVs);
	Eigen::VectorXd evaluateQuadraticErrorBases(
			Eigen::Vector4d linearBasisFunctions);
	Eigen::VectorXd evaluateCubicErrorBases(
			Eigen::Vector4d linearBasisFunctions);

	vtkSmartPointer<vtkUnstructuredGrid> createVtkMesh();

	vtkSmartPointer<vtkUnstructuredGrid> vtkMesh_ptr;
	vtkSmartPointer<vtkCellTreeLocator> vtkCellTree_ptr;
	std::vector<entHandle> vtkToIMesh;

	iMesh_Instance meshInstance;
	// TODO: rename vtkInputMesh something like inputMeshIsVtk
	bool vtkInputMesh;
	iBase_EntitySetHandle rootEntitySet;
//	iBase_EntitySetHandle vextexEntitySet;
//	iBase_EntitySetHandle surfaceEntitySet;
//	iBase_EntitySetHandle volumeEntitySet;

	std::vector<entHandle> allVertices;
	std::vector<entHandle> allFaces;
	std::vector<entHandle> allElements;

	std::map<entHandle,int> indexOfVertices;
	std::map<entHandle,int> indexOfFaces;
	std::map<entHandle,int> indexOfElements;

	entHandle previousCoordsToBasisElement;
	Eigen::Matrix4d previousCoordsToBasis;

	std::map<entHandle,std::vector<entHandle> > adjacentTetsMap;
	std::map<entHandle,std::vector<entHandle> > adjacentFacesMap;
	std::map<entHandle,std::vector<entHandle> > adjacentVertsMap;
	std::map<entHandle,std::vector<entHandle> > surroundingVertsMap;
	std::map<entHandle,std::vector<vect3d> > vertexVectorsMap;
	std::map<entHandle,std::vector<entHandle> > adjacentTetsToFaceMap;
};

#endif /* MESH_H_ */
