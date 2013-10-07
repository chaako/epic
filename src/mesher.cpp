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
	double scaleFactor=0.5;
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

//		dirs[0][0]=1.;
//		dirs[0][1]=0;
//		dirs[0][2]=0;
//		dirs[1][0]=0;
//		dirs[1][1]=1.;
//		dirs[1][2]=0;
//		dirs[2][0]=0;
//		dirs[2][1]=0;
//		dirs[2][2]=1.;

		// Refine in near object
		vect3d coords(xyz[0],xyz[1],xyz[2]);
		vect3d rHat=coords/coords.norm();
		vect3d perpX=rHat.cross(vect3d::UnitX());
		vect3d perpY=rHat.cross(vect3d::UnitY());
		vect3d perpZ=rHat.cross(vect3d::UnitZ());
//		h[0] *= pow(coords.norm(),2.)/10.;
		double transitionRadius=1.5;
		double wakeScaleFactor=0.5;
		double radialScaleLength =
				min(exp((coords.norm()-transitionRadius)/0.2)/4./scaleFactor,
						coords.norm()/transitionRadius*wakeScaleFactor);
//		double radialScaleLength=0.1;

		vect3d perp1, perp2;
		if (perpX.norm()<sqrt(LENGTH_TOLERANCE)) {
			perp1 = perpY;
			perp2 = perpZ;
		} else if (perpY.norm()<sqrt(LENGTH_TOLERANCE)) {
			perp1 = perpX;
			perp2 = perpZ;
		} else {
			perp1 = perpX;
			perp2 = perpY;
		}

//		double t=1./pow(coords.norm(),2.);
		double t;
		if (coords.norm()>transitionRadius) {
			t=0.;
		} else {
			t=1.;
		}
		vect3d scaleLengths(1.,1.,20.);
		double drift = 0.1;
		vect3d wakeDirection(0.,drift,1.);
		if (coords[2]<0)
			wakeDirection[2] *= -1.;
		wakeDirection /= wakeDirection.norm();
		vect3d wakeCenter = coords.dot(wakeDirection)*wakeDirection;
		double distanceFromCenterOfWake = (coords-wakeCenter).norm();
		if (distanceFromCenterOfWake<transitionRadius)
			scaleLengths *= wakeScaleFactor;
		vect3d tiltedY = wakeDirection.cross(vect3d::UnitX());
//		if (coords[2]<0)
//			tiltedY *= -1.;
		scaleLengths[0] = t*radialScaleLength + (1.-t)*scaleLengths[0];
		scaleLengths[1] = t + (1.-t)*scaleLengths[1];
		scaleLengths[2] = t + (1.-t)*scaleLengths[2];
		scaleLengths *= scaleFactor;

		// TODO: be more clever about which vectors to combine
//		vect3d d0 = t*rHat + (1.-t)*vect3d::UnitX();
//		vect3d d1 = t*perp1 + (1.-t)*vect3d::UnitY();
//		vect3d d2 = t*perp2 + (1.-t)*vect3d::UnitZ();
		vect3d d0 = t*rHat + (1.-t)*vect3d::UnitX();
		vect3d d1 = t*perp1 + (1.-t)*tiltedY;
		vect3d d2 = t*perp2 + (1.-t)*wakeDirection;
		d0 /= d0.norm();
		d1 /= d1.norm();
		d2 /= d2.norm();

		h[0] = scaleLengths[0];
		h[1] = scaleLengths[1];
		h[2] = scaleLengths[2];

		dirs[0][0]=d0[0];
		dirs[0][1]=d0[1];
		dirs[0][2]=d0[2];
		dirs[1][0]=d1[0];
		dirs[1][1]=d1[1];
		dirs[1][2]=d1[2];
		dirs[2][0]=d2[0];
		dirs[2][1]=d2[1];
		dirs[2][2]=d2[2];

		((PWLsfield *)field)->setSize((pEntity)vt,dirs,h);
	}
	VIter_delete (vit);
	return 1;
}

int main(int argc, char *argv[]) {
	string inputMeshFile;
	int numberOfIterations;
	bool isSurfaceMesh, refineMesh;
	double netgenFineness, netgenGrading, netgenMinh, netgenMaxh;
	// TODO: make names more transparent?
	{
		namespace po = boost::program_options;
		try {
			vector<double> evalPosX, evalPosY, evalPosZ;
			po::options_description desc("Allowed options");
			desc.add_options()
					("help", "produce help message")
					("inputFile", po::value<string>(&inputMeshFile), "input file")
					("isSurfaceMesh", po::value<bool>(&isSurfaceMesh)->default_value(true),
							"whether input is only surface mesh (true/false)")
					("refineMesh", po::value<bool>(&refineMesh)->default_value(true),
							"whether to do mesh refinement (true/false)")
					("numberOfIterations", po::value<int>(&numberOfIterations)->default_value(2),
							"number of refinement iterations")
					("netgenFineness", po::value<double>(&netgenFineness)->default_value(1.),
							"Netgen fineness parameter")
					("netgenGrading", po::value<double>(&netgenGrading)->default_value(0.25),
							"Netgen grading parameter")
					("netgenMinh", po::value<double>(&netgenMinh)->default_value(0.),
							"Netgen minh parameter")
					("netgenMaxh", po::value<double>(&netgenMaxh)->default_value(1000.),
							"Netgen maxh parameter")
			;

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);
			po::notify(vm);

			if (vm.count("help")) {
				cout << desc << "\n";
				exit(1);
			}

			if (vm.count("inputFile")) {
				cout << "Input file: " << inputMeshFile << endl;
			} else {
				cout << "Error: --inputFile was not set" << endl;
				cout << desc << "\n";
				exit(1);
			}
		} catch(exception& e) {
			cerr << "error: " << e.what() << "\n";
			exit(1);
		} catch(...) {
			cerr << "Exception of unknown type!\n";
		}
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
	stringstream volumeMeshFile;
	if (isSurfaceMesh) {
		if (inputMeshFile.find(".vtu")==string::npos) {
			cout << "ERROR: currently expect surface mesh to be .vtu file." << endl;
			exit(1);
		}
		SurfaceMesh surfaceMesh(inputMeshFile);
//		double rotationAngle = M_PI/5;
//		vect3d scaleFactors(1./4.,1./1.,1./10.);
		vect3d scaleFactors(1.,1.,1.);
		vect3d inverseScaleFactors(1.,1.,1.);
		for (int i=0; i<NDIM; i++) {
			inverseScaleFactors[i] = 1./scaleFactors[i];
		}
		vect3d origin(0.,0.,0.);
		vect3d noTranslation(0.,0.,0.);
		vect3d noScaling(1.,1.,1.);
		double noRotation=0;
//		vect3d translation(-2.,0.,0.);
		// TODO: don't hard code cell_codes
		surfaceMesh.transformSurface(5, origin, vect3d::UnitY(),
				noRotation, vect3d(0.4,0.6,5.), vect3d(0.,2.,0.));
		surfaceMesh.transformSurface(4, origin, vect3d::UnitY(),
				noRotation, scaleFactors, noTranslation);
		surfaceMesh.transformSurface(5, origin, vect3d::UnitY(),
				noRotation, scaleFactors, noTranslation);
//		surfaceMesh.transformSurface(5, vect3d(0.,0.,0.), vect3d::UnitY(),
//				0., vect3d(0.7,0.5,1.3), vect3d(0.,0.,0.));
//		surfaceMesh.transformSurface(5, vect3d(0.,0.,0.), vect3d::UnitY(),
//				0., scaleFactors, vect3d(0.,0.,0.));
//		surfaceMesh.transformSurface(4, vect3d(0.,0.,0.), vect3d::UnitY(),
//				rotationAngle, scaleFactors, translation);
		surfaceMesh.createVolumeMesh(netgenFineness,netgenGrading,netgenMinh,netgenMaxh);
		surfaceMesh.scaleVolumeMesh(origin, inverseScaleFactors);
		int periodLocation = inputMeshFile.rfind(".vtu");
		volumeMeshFile << inputMeshFile.substr(0,periodLocation)
				<< "_meshed" << ".vtk";
//		<< "_meshed_" << fixed << rotationAngle/(2.*M_PI) << ".vtk";
		surfaceMesh.saveVolumeMesh(volumeMeshFile.str());

		// TODO: don't hard code cell_codes
		surfaceMesh.transformSurface(4, origin, vect3d::UnitY(),
				0., inverseScaleFactors, noTranslation);
		surfaceMesh.transformSurface(5, origin, vect3d::UnitY(),
				0., inverseScaleFactors, noTranslation);
		stringstream transformedSurfaceMeshFile;
		transformedSurfaceMeshFile << inputMeshFile.substr(0,periodLocation)
				<< "_transformed" << ".vtu";
//		<< "_rotated_" << fixed << rotationAngle/(2.*M_PI) << ".vtu";
		surfaceMesh.save(transformedSurfaceMeshFile.str());
	} else {
		volumeMeshFile << inputMeshFile;
	}

	stringstream refinedMeshFile;
	if (refineMesh) {
		Mesh coarseMesh(volumeMeshFile.str());
		Field<int> faceType(&coarseMesh,string("cell_code"),iBase_FACE);
		Mesh refinedMesh(volumeMeshFile.str());
		Field<int> faceTypeRefined(&refinedMesh,string("cell_code"),iBase_FACE);

		refinedMesh.classifyBoundariesForMeshRefinement(faceTypeRefined);
		pMesh part;
		FMDB_Mesh_GetPart((mMesh*)refinedMesh.meshInstance, 0, part);
		pSField field=new PWLsfield(part);
		meshAdapt rdr(part,field,0,0);
		rdr.run(numberOfIterations,1, setBSizeField);

		int periodLocation = volumeMeshFile.str().rfind(".");
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
			// TODO: problem here if set VOLUME_TOLERANCE too small
			entHandle coarseEntity = coarseMesh.findTet(centroid, centroid,
					coarseMesh.entitiesVectors[iBase_REGION][0], &foundTet);
			vector<entHandle> coarseFaces =
					coarseMesh.getAdjacentEntities(coarseEntity, iBase_FACE);
			for (int j=0; j<coarseFaces.size(); j++) {
				// TODO: this doesn't work if two faces with different codes lie
				//       in the same plane
				double volume =
						coarseMesh.getTetVolume(centroid, coarseFaces[j]);
				// TODO: problem here if set VOLUME_TOLERANCE too small
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
