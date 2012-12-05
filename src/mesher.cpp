/*
 * mesher.cpp
 *
 *  Created on: Dec 4, 2012
 *      Author: chaako
 */

#include "mesher.h"

int setBSizeField(pMesh mesh, pSField field, void *)
{
  pVertex vt;
  double h[3], dirs[3][3], xyz[3], norm;
  VIter vit=M_vertexIter(mesh);
  while( vt=VIter_next(vit) ) {
    h[0] = 1.0;
    h[1] = .2;
    h[2] = 2.5;

    dirs[0][0]=1.;
    dirs[0][1]=0;
    dirs[0][2]=0;
    dirs[1][0]=0;
    dirs[1][1]=1.;
    dirs[1][2]=0;
    dirs[2][0]=0;
    dirs[2][1]=0;
    dirs[2][2]=1.;

    ((PWLsfield *)field)->setSize((pEntity)vt,dirs,h);
  }
  VIter_delete (vit);
  return 1;
}

int main(int argc, char *argv[]) {
	if (argc<2) {
		printf("usage: %s meshin\n", argv[0]);
		exit(1);
	}

	string inputMeshFile(argv[1]);
	stringstream volumeMeshFile;
	// TODO: distinguish between surface and volume mesh in less obscure way than format
	if (inputMeshFile.find(".vtu")!=string::npos) {
		SurfaceMesh surfaceMesh(inputMeshFile);
		double rotationAngle = M_PI/5;
		// TODO: don't hard code collector cell_code
		surfaceMesh.rotateSurface(4, vect3d(0.,0.,0.), vect3d::UnitY(), rotationAngle);
		int periodLocation = inputMeshFile.rfind(".vtu");
		stringstream rotatedSurfaceMeshFile;
		rotatedSurfaceMeshFile << inputMeshFile.substr(0,periodLocation)
						<< "_rotated_" << rotationAngle << ".vtu";
		surfaceMesh.save(rotatedSurfaceMeshFile.str());
		surfaceMesh.createVolumeMesh();
		volumeMeshFile << inputMeshFile.substr(0,periodLocation)
						<< "_meshed_" << rotationAngle << ".vtk";
		surfaceMesh.saveVolumeMesh(volumeMeshFile.str());
	}

	stringstream refinedMeshFile;
	// TODO: distinguish between surface and volume mesh in less obscure way than format
	if (inputMeshFile.find(".vtu")!=string::npos) {
		Mesh coarseMesh(volumeMeshFile.str());
		Field<int> faceType(&coarseMesh,string("cell_code"),iBase_FACE);
		Mesh refinedMesh(volumeMeshFile.str());
		Field<int> faceTypeRefined(&refinedMesh,string("cell_code"),iBase_FACE);

		refinedMesh.classifyBoundariesForMeshRefinement(faceTypeRefined);
		pMesh part;
		FMDB_Mesh_GetPart((mMesh*)refinedMesh.meshInstance, 0, part);
		pSField field=new PWLsfield(part);
		meshAdapt rdr(part,field,0,0);
		rdr.run(2,1, setBSizeField);

		int periodLocation = volumeMeshFile.str().rfind(".vtk");
		refinedMeshFile << volumeMeshFile.str().substr(0,periodLocation)
						<< "_refined" << ".sms";
		FMDB_Mesh_WriteToFile(part->getMesh(), refinedMeshFile.str().c_str(), 0);

		// Update cell_code field of refined mesh
		Mesh reopenedMesh(refinedMeshFile.str());
		Field<int> faceTypeReopened(&reopenedMesh,string("cell_code"),iBase_FACE);
		for (int i=0; i<faceTypeReopened.entities.size(); i++) {
			vect3d centroid(0.,0.,0.);
			vector<vect3d> vVs = reopenedMesh.getVertexVectors(
					faceTypeReopened.entities[i]);
			centroid = (vVs[0]+vVs[1]+vVs[2])/3.;
			bool foundTet=false;
			entHandle coarseEntity = coarseMesh.findTet(centroid, centroid,
					coarseMesh.entitiesVectors[iBase_REGION][0], &foundTet);
			vector<entHandle> coarseFaces =
					coarseMesh.getAdjacentEntities(coarseEntity, iBase_FACE);
			for (int j=0; j<coarseFaces.size(); j++) {
				// TODO: this doesn't work if two faces with different codes lie
				//       in the same plane
				double volume =
						coarseMesh.getTetVolume(centroid, coarseFaces[j]);
				if (volume<VOLUME_TOLERANCE) {
					int cellCode=faceType[coarseFaces[j]];
					// TODO: Don't hard-code cell codes
					if (cellCode==4 || cellCode==5) {
						faceTypeReopened.setField(faceTypeReopened.entities[i],
								cellCode);
					}
				}
			}
		}
		FMDB_Mesh_WriteToFile((mMesh*)reopenedMesh.meshInstance,
				refinedMeshFile.str().c_str(), 0);
	} else {
		refinedMeshFile << volumeMeshFile;
	}

	return 0;
}




clock_t extern_findTet=0, extern_checkIfInNewTet=0;

// TODO: find better way to distinguish orbits in output
int extern_orbitNumber = 0;
