/*
 * Mesh.h
 *
 *  Created on: Feb 10, 2012
 *      Author: chaako
 */

#ifndef MESH_H_
#define MESH_H_

class Mesh {
public:
	Mesh(std::string inputMeshFile);
	virtual ~Mesh();

	void printElementNumbers();
	void save(std::string outputMeshFile);

	Eigen::Vector3d getCoordinates(iBase_EntityHandle node);
	iBase_EntityHandle findTet(Eigen::Vector3d position,
			iBase_EntityHandle adjacentTet, bool *tetFound, bool isTet=true);
	std::vector<iBase_EntityHandle> getVertices(iBase_EntityHandle element);
	iBase_EntityHandle getRandomVertex();

	iMesh_Instance meshInstance;
	iBase_EntitySetHandle rootEntitySet;
//	iBase_EntitySetHandle vextexEntitySet;
//	iBase_EntitySetHandle surfaceEntitySet;
//	iBase_EntitySetHandle volumeEntitySet;


};

#endif /* MESH_H_ */
