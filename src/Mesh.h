/*
 * Mesh.h
 *
 *  Created on: Feb 10, 2012
 *      Author: chaako
 */

#ifndef MESH_H_
#define MESH_H_

#include <map>

template <class T> class Field;

class Mesh {
public:
	Mesh(std::string inputMeshFile);
	virtual ~Mesh();

	void printElementNumbers();
	void save(std::string outputMeshFile);

	Eigen::Vector3d getCoordinates(iBase_EntityHandle node, bool useMap=false);
	iBase_EntityHandle findTet(Eigen::Vector3d oldPosition, Eigen::Vector3d position,
			iBase_EntityHandle adjacentTet, bool *tetFound, bool isTet=true);
	std::vector<iBase_EntityHandle> getVertices(iBase_EntityHandle element);
	std::vector<iBase_EntityHandle> getAdjacentEntities(iBase_EntityHandle element,
			int dimension);
	std::vector<iBase_EntityHandle> getEntities(int dimension);
	iBase_EntityHandle getRandomVertex();
	std::vector<iBase_EntityHandle> getFaces(iBase_EntityHandle element);
	bool vtkMesh;

	std::map<iBase_EntityHandle,std::vector<iBase_EntityHandle> >
	getAdjacentsMap(int keyEntityType, int valueEntityType,
			int bridgeEntityType=-1);
	int getEntityDimension(iBase_EntityHandle entity);
	std::vector<Eigen::Vector3d> getVertexVectors(iBase_EntityHandle entity,
			bool useMap=true);
	bool checkIfInTet(Eigen::Vector3d currentPosition,
			iBase_EntityHandle element);
	bool checkIfInTet(Eigen::Vector3d currentPosition,
			std::vector<Eigen::Vector3d> vertexVectors);
	bool checkIfIntersectsTriangle(Eigen::Vector3d previousPosition,
			Eigen::Vector3d currentPosition,
			std::vector<Eigen::Vector3d> vertexVectors);
	iBase_EntityHandle findFaceCrossed(iBase_EntityHandle previousElement,
			Eigen::Vector3d previousPosition, Eigen::Vector3d currentPosition);
	iBase_TagHandle getTagHandle(std::string tagName);
	iBase_TagHandle createTagHandle(std::string tagName, int size, int type);
	Eigen::Vector3d getSurfaceVector(iBase_EntityHandle face,
			Eigen::Vector3d point=Eigen::Vector3d(0.,0.,0.));
	Eigen::Vector3d getNormalVector(iBase_EntityHandle face,
			Eigen::Vector3d point=Eigen::Vector3d(0.,0.,0.));
	Eigen::Vector3d getVertexNormalVector(iBase_EntityHandle vertex,
			Field<int> faceTypeField);
	std::vector<iBase_EntityHandle> getSuperCellFaces(iBase_EntityHandle vertex);
	double getTetVolume(Eigen::Vector3d point, iBase_EntityHandle face);
	double getTetVolume(std::vector<Eigen::Vector3d> vertexVectors);
	std::vector<double> getTetSubVolumes(Eigen::Vector3d point,
			std::vector<Eigen::Vector3d> vertexVectors) ;
	std::vector<double> getVertexWeights(Eigen::Vector3d point,
			std::vector<Eigen::Vector3d> vertexVectors);
	Eigen::Vector4d getInterpolationCoeffs(Eigen::Vector3d position,
			iBase_EntityHandle element);
	Eigen::Vector4d getInterpolationCoeffs(Eigen::Vector3d position,
			std::vector<Eigen::Vector3d> vVs);

	iMesh_Instance meshInstance;
	iBase_EntitySetHandle rootEntitySet;
//	iBase_EntitySetHandle vextexEntitySet;
//	iBase_EntitySetHandle surfaceEntitySet;
//	iBase_EntitySetHandle volumeEntitySet;

	std::map<iBase_EntityHandle,std::vector<iBase_EntityHandle> > adjacentTetsMap;
	std::map<iBase_EntityHandle,std::vector<iBase_EntityHandle> > adjacentFacesMap;
	std::map<iBase_EntityHandle,std::vector<iBase_EntityHandle> > adjacentVertsMap;
	std::map<iBase_EntityHandle,std::vector<Eigen::Vector3d> > vertexVectorsMap;
	std::map<iBase_EntityHandle,std::vector<iBase_EntityHandle> > adjacentTetsToFaceMap;
};

#endif /* MESH_H_ */
