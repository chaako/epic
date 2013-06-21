/*
 * SurfaceMesh.h
 *
 *  Created on: Dec 3, 2012
 *      Author: chaako
 */

#ifndef SURFACEMESH_H_
#define SURFACEMESH_H_

#include "typesAndDefinitions.h"
#include "variables.h"
#include "functions.h"

//#include "FMDB.h"
//#include "MeshAdapt.h"
//#include "AdaptUtil.h"
//#include "PWLinearSField.h"

namespace nglib {
#include "nglib.h"
}
//#include <ngexception.hpp>

class SurfaceMesh {
public:
	SurfaceMesh();
	SurfaceMesh(string inputFile);
	virtual ~SurfaceMesh();

	void load(string inputFile);
	void save(string outputFile);
	void saveVolumeMesh(string outputFile);

	void transformSurface(int surfaceId, vect3d origin,
			vect3d rotationAxis, double rotationAngle,
			vect3d scaleFactors, vect3d translation);
	void scaleVolumeMesh(vect3d origin, vect3d scaleFactors);
	void createVolumeMesh();
	vector<vect3d> getPoints(
			vect3dMap *vtkIdOfPoint,
			int surfaceCode=-1);

	vtkSmartPointer<vtkUnstructuredGrid> vtkMesh;
	vtkSmartPointer<vtkUnstructuredGrid> volumeMesh;

	vector<vect3d> vertexPositions;
	vect3dMap idOfVertex;
};

#endif /* SURFACEMESH_H_ */
