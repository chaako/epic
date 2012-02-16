/*
 * epic.h
 *
 *  Created on: Jan 20, 2012
 *      Author: chaako
 */

#ifndef EPIC_H_
#define EPIC_H_

#define VOLUME_TOLERANCE 1.e-8

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "iMesh.h"
#include <math.h>

#include "Eigen/Dense"

#include <numeric>
#include <vector>
//#include <type_traits>
#include <boost/type_traits.hpp>

#include "Mesh.h"
#include "Field.h"
#include "Orbit.h"

#define CHECK(a) if (iBase_SUCCESS != ierr) printf("%s\n", a), exit(ierr)

class mMesh;
int custom_importVTK(mMesh *, const char *);

Eigen::Vector3d getSurfaceVector(iMesh_Instance mesh, Eigen::Vector3d point,
		iBase_EntityHandle face);
std::vector<iBase_EntityHandle> getSuperCellFaces(iMesh_Instance mesh,
		iBase_EntityHandle vertex);
double getAverageDblData(iMesh_Instance mesh, iBase_EntityHandle entity,
		iBase_TagHandle dblData_tag);
double getTetVolume(iMesh_Instance mesh, Eigen::Vector3d point,
		iBase_EntityHandle face);
double getTetVolume(std::vector<Eigen::Vector3d> vertexVectors);
std::vector<Eigen::Vector3d> getVertexVectors(iMesh_Instance mesh,
		iBase_EntityHandle entity);
bool checkIfInTet(Eigen::Vector3d currentPosition, iMesh_Instance mesh,
		iBase_EntityHandle element);
bool checkIfInTet(Eigen::Vector3d currentPosition,
		std::vector<Eigen::Vector3d> vertexVectors);
std::vector<double> getTetSubVolumes(Eigen::Vector3d point,
		std::vector<Eigen::Vector3d> vertexVectors);
std::vector<double> getVertexWeights(Eigen::Vector3d point,
		std::vector<Eigen::Vector3d> vertexVectors);



#endif /* EPIC_H_ */
