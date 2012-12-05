/*
 * SurfaceMesh.h
 *
 *  Created on: Dec 3, 2012
 *      Author: chaako
 */

#ifndef SURFACEMESH_H_
#define SURFACEMESH_H_

#include "typesAndDefinitions.h"

#include "FMDB.h"
#include "MeshAdapt.h"
#include "AdaptUtil.h"
#include "PWLinearSField.h"

namespace nglib {
#include "nglib.h"
}
//#include <ngexception.hpp>

class SurfaceMesh {
public:
	SurfaceMesh(string inputFile);
	virtual ~SurfaceMesh();

	void load(string inputFile);
	void save(string outputFile);
	void saveVolumeMesh(string outputFile);

	void rotateSurface(int surfaceId, vect3d origin,
			vect3d rotationAxis, double rotationAngle);
	void createVolumeMesh();

	vtkSmartPointer<vtkUnstructuredGrid> vtkMesh;
	vtkSmartPointer<vtkUnstructuredGrid> volumeMesh;
};

#endif /* SURFACEMESH_H_ */
