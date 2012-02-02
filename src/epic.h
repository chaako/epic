/*
 * epic.h
 *
 *  Created on: Jan 20, 2012
 *      Author: chaako
 */

#ifndef EPIC_H_
#define EPIC_H_

#endif /* EPIC_H_ */

#include <vector>


class mMesh;
int custom_importVTK(mMesh *, const char *);

Eigen::Vector3d getSurfaceVector(iMesh_Instance mesh, Eigen::Vector3d point,
		iBase_EntityHandle face);
std::vector<iBase_EntityHandle> getSuperCellFaces(iMesh_Instance mesh,
		iBase_EntityHandle vertex);
double getAverageDblData(iMesh_Instance mesh, iBase_EntityHandle entity,
		iBase_TagHandle dblData_tag);
double getVolumeBetweenPointAndFace(iMesh_Instance mesh, Eigen::Vector3d point,
		iBase_EntityHandle face);
