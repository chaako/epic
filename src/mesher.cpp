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
	double scaleFactor=0.05;
	while( vt=VIter_next(vit) ) {
		V_coord(vt,xyz);
		h[0] = 1.*scaleFactor;
		h[1] = 1.*scaleFactor;
		h[2] = 1.*scaleFactor;
//		// Refine in near object
//		vect3d range(1.3,0.7,3.);
//		vect3d origin(-1,0.,0.);
//		vect3d coords(xyz[0],xyz[1],xyz[2]);
//		vect3d relPos=coords-origin;
//		if (fabs(relPos[0])<range[0] &&
//				fabs(relPos[1])<range[1] &&
//				fabs(relPos[2])<range[2]) {
//			h[0] /= (1+2.*fabs(fabs(relPos[0])-range[0])/range[0]);
//			h[1] /= (1+1.*fabs(fabs(relPos[1])-range[1])/range[1]);
//			h[2] /= (1+5.*fabs(fabs(relPos[2])-range[2])/range[2]);
//		}

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
	if (argc<3) {
		printf("usage: %s meshin Niter\n", argv[0]);
		exit(1);
	}

	//	// TODO: remove this testing ground
	//	{
	////		int playground_meshAdapt(int argc, char* argv[]);
	////		playground_meshAdapt(argc, argv);
	//		int playground_netgen(int argc, char* argv[]);
	//		playground_netgen(argc, argv);
	//		exit(0);
	//	}

	cout.precision(3);
//	for (double rotationAngle=0.; rotationAngle<M_PI; rotationAngle+=M_PI/12.) {
//	for (double rotationAngle=0.125*2.*M_PI; rotationAngle<0.126*2.*M_PI; rotationAngle+=M_PI/12.) {
	string inputMeshFile(argv[1]);
	stringstream volumeMeshFile;
	// TODO: distinguish between surface and volume mesh in less obscure way than format
	if (inputMeshFile.find(".vtu")!=string::npos) {
		SurfaceMesh surfaceMesh(inputMeshFile);
//		double rotationAngle = M_PI/5;
		vect3d scaleFactors(1./4.,1./1.,1./10.);
		vect3d inverseScaleFactors(1.,1.,1.);
		for (int i=0; i<NDIM; i++) {
			inverseScaleFactors[i] = 1./scaleFactors[i];
		}
		vect3d translation(-1.,0.,0.);
		// TODO: don't hard code collector cell_code
//		surfaceMesh.transformSurface(5, vect3d(0.,0.,0.), vect3d::UnitY(),
//				0., vect3d(0.7,0.5,1.3), vect3d(0.,0.,0.));
//		surfaceMesh.transformSurface(5, vect3d(0.,0.,0.), vect3d::UnitY(),
//				0., scaleFactors, vect3d(0.,0.,0.));
//		surfaceMesh.transformSurface(4, vect3d(0.,0.,0.), vect3d::UnitY(),
//				rotationAngle, scaleFactors, translation);
		surfaceMesh.createVolumeMesh();
//		surfaceMesh.scaleVolumeMesh(vect3d(0.,0.,0.), inverseScaleFactors);
		int periodLocation = inputMeshFile.rfind(".vtu");
		volumeMeshFile << inputMeshFile.substr(0,periodLocation)
				<< "_meshed" << ".vtk";
//		<< "_meshed_" << fixed << rotationAngle/(2.*M_PI) << ".vtk";
		surfaceMesh.saveVolumeMesh(volumeMeshFile.str());

		// TODO: don't hard code collector cell_code
//		surfaceMesh.transformSurface(4, vect3d(0.,0.,0.), vect3d::UnitY(),
//				0., inverseScaleFactors, vect3d(0.,0.,0.));
//		surfaceMesh.transformSurface(5, vect3d(0.,0.,0.), vect3d::UnitY(),
//				0., inverseScaleFactors, vect3d(0.,0.,0.));
		stringstream rotatedSurfaceMeshFile;
		rotatedSurfaceMeshFile << inputMeshFile.substr(0,periodLocation)
				<< "_rotated" << ".vtu";
//		<< "_rotated_" << fixed << rotationAngle/(2.*M_PI) << ".vtu";
		surfaceMesh.save(rotatedSurfaceMeshFile.str());
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
		rdr.run(atoi(argv[2]),1, setBSizeField);

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
//	}
	return 0;
}
