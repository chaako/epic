#include "typesAndDefinitions.h"

int playground_netgen(int argc, char* argv[]) {
	// Read the input .vtu file (only surface mesh)
	string inputFile = argv[1];
	if (inputFile.find(".vtu")==string::npos)
		throw; // Not a .vtu file
	vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
			vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
	reader->SetFileName(inputFile.c_str());
	reader->Update();
//	vtkSmartPointer<vtkUnstructuredGrid> vtkMesh = reader->GetOutput();
//	vtkMesh->Print(cout);
	vtkSmartPointer<vtkUnstructuredGrid> vtkMesh =
			vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkMesh->DeepCopy(reader->GetOutput());

	// Rotate and/or translate collecting surface
	vector<bool> transformNode(vtkMesh->GetNumberOfPoints(), false);
	string cell_code_name = "cell_code";
	vtkIntArray* orig_cell_codes = vtkIntArray::SafeDownCast(
			vtkMesh->GetCellData()->GetArray(cell_code_name.c_str()));
	for (vtkIdType id_cell = 0; id_cell < vtkMesh->GetNumberOfCells(); ++id_cell) {
		int cell_code = orig_cell_codes->GetValue(id_cell);
		if (cell_code==4) {
			vtkIdType *pts, N_pts;
			vtkMesh->GetCellPoints(id_cell, N_pts, pts);
			for (int i = 0; i < 3; ++i) {
				transformNode[pts[i]] = true;
			}
		}
	}
	for (vtkIdType id_node = 0; id_node < vtkMesh->GetNumberOfPoints(); ++id_node) {
		if (transformNode[id_node]) {
			double x[3];
			vtkMesh->GetPoints()->GetPoint(id_node, x);
			vect3d pos(x[0],x[1],x[2]);
			double angle=M_PI/3;
			pos = Eigen::AngleAxisd(-angle, vect3d::UnitY())*pos;
			pos += vect3d(1.,0.,0.);
			x[0] = pos[0];
			x[1] = pos[1];
			x[2] = pos[2];
			vtkMesh->GetPoints()->SetPoint(id_node, x);
		}
	}

	// Write the transformed .vtu file (only surface mesh)
	{
		vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
				vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
		stringstream outputFilename;
		int periodLocation = inputFile.rfind(".vtu");
		outputFilename << inputFile.substr(0,periodLocation)
						<< "_tranformed.vtu";
		writer->SetFileName(outputFilename.str().c_str());
		writer->SetInput(vtkMesh);
		writer->Write();
	}

	// Generate volume mesh (based on Engrid's createvolumemesh.cpp)
	using namespace nglib;
	Ng_Init();
	Ng_Meshing_Parameters mp;
	mp.fineness = 0.0;
	Ng_Mesh *ngMesh = Ng_NewMesh();
	mp.maxh = 1000.;
	mp.minh = 0.;
	mp.grading = 1;
	vector<vtkIdType> ng2eg(vtkMesh->GetNumberOfPoints()+1);
	vector<vtkIdType> eg2ng(vtkMesh->GetNumberOfPoints());
	{
		int N = 1;
		for (vtkIdType id_node = 0; id_node < vtkMesh->GetNumberOfPoints(); ++id_node) {
			double x[3];
			vtkMesh->GetPoints()->GetPoint(id_node, x);
			Ng_AddPoint(ngMesh, x);
			ng2eg[N] = id_node;
			eg2ng[id_node] = N;
			++N;
		}
	}
	vector<vector<vtkIdType> > tri;
	for (vtkIdType id_cell = 0; id_cell < vtkMesh->GetNumberOfCells(); ++id_cell) {
		vtkIdType type_cell = vtkMesh->GetCellType(id_cell);
		if (type_cell != VTK_TRIANGLE)
			throw; // Assuming all cells are triangular surface faces
		int trig[3];
		vtkIdType *pts, N_pts;
		vtkMesh->GetCellPoints(id_cell, N_pts, pts);
		for (int i = 0; i < 3; ++i) {
			trig[i] = eg2ng[pts[i]];
		}
		Ng_AddSurfaceElement(ngMesh, NG_TRIG, trig);
	}
	Ng_Result res;
	res = Ng_GenerateVolumeMesh (ngMesh, &mp);
	vtkSmartPointer<vtkUnstructuredGrid> volumeMesh =
			vtkSmartPointer<vtkUnstructuredGrid>::New();
//	volumeMesh->DeepCopy(vtkMesh);
	if (res == NG_OK) {
		int Npoints_ng = Ng_GetNP(ngMesh);
		int Ncells_ng  = Ng_GetNE(ngMesh);
		int Nscells_ng = Ng_GetNSE(ngMesh);

		vtkIdType new_point = 0;
//		vtkIdType new_point = vtkMesh->GetNumberOfPoints();
//		volumeMesh->GetPoints()->SetNumberOfPoints(Npoints_ng);

		// copy existing points
		vector<vtkIdType> old2new(vtkMesh->GetNumberOfPoints(), -1);
		vtkSmartPointer<vtkPoints> points =
				vtkSmartPointer<vtkPoints>::New();
		points->SetNumberOfPoints(Npoints_ng);
		for (vtkIdType id_point = 0; id_point < vtkMesh->GetNumberOfPoints(); ++id_point) {
			double x[3];
			vtkMesh->GetPoints()->GetPoint(id_point, x);
			//	      volumeMesh->GetPoints()->SetPoint(new_point, x);
			points->SetPoint(new_point, x);
			old2new[id_point] = new_point;
			++new_point;
		}
		volumeMesh->SetPoints(points);

		// mark all surface nodes coming from NETGEN
		vector<bool> ng_surf_node(Npoints_ng + 1, false);
		for (int i = 1; i <= Nscells_ng; ++i) {
			int pts[8];
			Ng_Surface_Element_Type ng_type;
			ng_type = Ng_GetSurfaceElement(ngMesh, i, pts);
			int N = 0;
			if (ng_type == NG_TRIG) {
				N = 3;
			} else {
				throw; // Only support triangular surface elements
			}
			for (int j = 0; j < N; ++j) {
				ng_surf_node[pts[j]] = true;
			}
		}

		// add new points from NETGEN
		vector<vtkIdType> ng2new(Npoints_ng+1, -1);
		for (int i = 1; i <= Npoints_ng; ++i) {
			if (!ng_surf_node[i]) {
				double x[3];
				Ng_GetPoint(ngMesh, i, x);
				volumeMesh->GetPoints()->SetPoint(new_point, x);
				ng2new[i] = new_point;
				++new_point;
			}
		}

		vtkSmartPointer<vtkIntArray> cell_codes =
				vtkSmartPointer<vtkIntArray>::New();
		cell_codes->SetNumberOfComponents(1);
		string cell_code_name = "cell_code";
		cell_codes->SetName(cell_code_name.c_str());
		vtkIntArray* orig_cell_codes = vtkIntArray::SafeDownCast(
				vtkMesh->GetCellData()->GetArray(cell_code_name.c_str()));

		// copy existing cells
		vector<vtkIdType> old2new_cell(vtkMesh->GetNumberOfCells(), -1);
		for (vtkIdType id_cell = 0; id_cell < vtkMesh->GetNumberOfCells(); ++id_cell) {
			bool ins_cell = false;
			vtkIdType type_cell = vtkMesh->GetCellType(id_cell);
			if (type_cell == VTK_TRIANGLE) ins_cell = true;
			if (ins_cell) {
				vtkIdType N_pts, *pts;
				vtkMesh->GetCellPoints(id_cell, N_pts, pts);
				for (int i = 0; i < N_pts; ++i) {
					pts[i] = old2new[pts[i]];
					if (pts[i] == -1) {
						throw;
					}
				}
				vtkIdType id_new = volumeMesh->InsertNextCell(type_cell, N_pts, pts);
				old2new_cell[id_cell] = id_new;
				if (vtkMesh->GetCellData()->GetArray(cell_code_name.c_str())) {
//				    if (vtkMesh->GetCellData()->GetScalars(cell_code_name.c_str())) {
//				    	int cell_code = vtkMesh->GetCellData()->
//					    		GetArray(cell_code_name.c_str())->GetValue(id_cell);
					int cell_code = orig_cell_codes->GetValue(id_cell);
//				    	cell_codes->SetValue(id_new, cell_code);
					cell_codes->InsertNextValue(cell_code);
//				    }
				}
			}
		}

		// add new cells
		vtkIdType id_new_cell;
		for (vtkIdType cellId = 0; cellId < Ncells_ng; ++cellId) {
			int       pts[8];
			vtkIdType new_pts[4];
			for (int i = 0; i < 8; ++i) {
				pts[i] = 0;
			}
			Ng_Volume_Element_Type ng_type;
			ng_type = Ng_GetVolumeElement(ngMesh, cellId + 1, pts);
			if (ng_type != NG_TET) {
				throw; // Only support tetrahedral meshes
			}
			for (int i = 0; i < 4; ++i) {
				if (!ng_surf_node[pts[i]]) {
					new_pts[i] = ng2new[pts[i]];
				} else {
					new_pts[i] = ng2eg[pts[i]];
				}
			}
			if (ng_type == NG_TET) {
				vtkIdType tet[4];
				tet[0] = new_pts[0];
				tet[1] = new_pts[1];
				tet[2] = new_pts[3];
				tet[3] = new_pts[2];
				id_new_cell = volumeMesh->InsertNextCell(VTK_TETRA, 4, tet);
				cell_codes->InsertNextValue(1);
			} else {
				throw; // Only support tetrahedral meshes
			}
		}

		volumeMesh->GetCellData()->AddArray(cell_codes);

	} else {
		cout << "NETGEN did not succeed.\nPlease check if the surface mesh is oriented correctly" <<
				" (red edges on the outside of the domain)";
		throw; // Netgen failed
	}

	// Write the output .vtk file (full volume mesh)
	vtkSmartPointer<vtkUnstructuredGridWriter> writer =
			vtkSmartPointer<vtkUnstructuredGridWriter>::New();
	stringstream outputFilename;
	int periodLocation = inputFile.rfind(".vtu");
	outputFilename << inputFile.substr(0,periodLocation)
					<< "_meshed.vtk";
	writer->SetFileName(outputFilename.str().c_str());
	writer->SetInput(volumeMesh);
	writer->Write();

	Ng_DeleteMesh(ngMesh);
	Ng_Exit();
	cout << "\n\nNETGEN call finished" << endl;
	cout << endl;
	return 0;
}













#include "MeshAdapt.h"
#include "AdaptUtil.h"
#include "PWLinearSField.h"
#include <stdlib.h>
#include <stdio.h>

#include "DiscUtil.h"

#include <math.h>
#include <fstream>
#include <iostream>


// anisotropic field one - no rotation, but with jumps
//  void setSizeField(pMesh mesh, pMField field, double lower, double up)
//  {
//    pVertex vt;
//    double xyz[3],h[3], dirs[3][3];
//    VIter vit=M_vertexIter(mesh);
//    while( vt=VIter_next(vit) ) {
//      V_coord(vt,xyz);

//      dirs[0][0]=1.0;
//      dirs[0][1]=0.0;
//      dirs[0][2]=0.0;
//      dirs[1][0]=0.0;
//      dirs[1][1]=1.0;
//      dirs[1][2]=0.0;
//      dirs[2][0]=0;
//      dirs[2][1]=0;
//      dirs[2][2]=1.;
    
//      if( xyz[0]<lower  ||  xyz[0]>up )
//        h[0]=0.1;
//      else
//        h[0]=0.005;
    
//      if( xyz[1]<lower  ||  xyz[1]>up )
//        h[1]=0.1;
//      else
//        h[1]=0.005;
    
//      h[2] = 1.0;

//      field->setMetric((pEntity)vt,dirs,h);

//    }
//    VIter_delete (vit);
//  }


// anisotropic size field two: with rotation
double R; 
double L;
int setSizeField(pMesh mesh, pSField field, void *)
{
  pVertex vt;
  double h[3], dirs[3][3], xyz[3], norm;
  VIter vit=M_vertexIter(mesh);
  while( vt=VIter_next(vit) ) {
//    V_coord(vt,xyz);
//    double circle= fabs(xyz[0] * xyz[0] +
//      xyz[1] * xyz[1] +
//      xyz[2] * xyz[2] - R*R);
//    h[0] = .5 * fabs(1. - exp (-circle*L)) + 1.e-3;
//    h[1] = .5;
//    h[2] = 1.;
    h[0] = .1;
    h[1] = .02;
    h[2] = .5;

//    norm=sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
//    dirs[0][0]=xyz[0]/norm;
//    dirs[0][1]=xyz[1]/norm;
//    dirs[0][2]=0;
//    dirs[1][0]=-1.0*xyz[1]/norm;
//    dirs[1][1]=xyz[0]/norm;
//    dirs[1][2]=0;
//    dirs[2][0]=0;
//    dirs[2][1]=0;
//    dirs[2][2]=1.;
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


int playground_meshAdapt(int argc, char* argv[])
{  

 if ( argc!=3 ) { 
   printf("Usage: %s  meshfile <# of iterations>\n",argv[0]);
   return 0;
 }
 
  char model_file[256];
  char mesh_file[256];
  char outmesh[256];
  char without_extension[256];

  snprintf(without_extension,strlen(argv[1])-3,"%s",argv[1]);
  sprintf(mesh_file,"%s",argv[1]);

  MS_init();

//  sprintf(model_file,"%s.xmt_txt",without_extension);
//  pGModel model=GM_createFromParasolidFile(model_file);
//  pMesh mesh=MS_newMesh(model);
  mMesh* meshInstance;
//  pGeomMdl geomModel=new meshModel::DiscreteModel();
//  pGeomMdl geomModel=new DiscreteModel(string(mesh_file),3,45,45,0,
//		  "newMeshName.sms","dmgFileName.dmg");
//  pMesh mesh=MS_newMesh(geomModel);
//  pGeomMdl geomModel=new SGModel("asdf");
//  FMDB_Mesh_Create(geomModel, meshInstance);
  FMDB_Mesh_Create(NULL, meshInstance);
//  mPart part(1,NULL,meshInstance);
//  pMesh mesh=&part;
//  pMesh mesh;
//  FMDB_Mesh_GetPart(meshInstance, 0, mesh);
  pMesh part;
  FMDB_Mesh_GetPart(meshInstance, 0, part);
//  DM_readSMS(mesh, string(mesh_file), 0);
//  pMesh mesh=MS_newMesh(NULL);
//  pSField field=new PWLsfield(mesh);
  pSField field=new PWLsfield(part);
//  FMDB_Mesh_LoadFromFile (meshInstance, mesh_file, 0);
//  int custom_importVTK(mMesh *, const char *);
//  custom_importVTK(mesh->getMesh(), mesh_file);
//  M_load(mesh,mesh_file);

  // Build tetrahedral domain from scratch
  typedef GEntity* (*entityBy_FP)(SGModel*,int,int);
  entityBy_FP fp = GM_entityByTag;

  vector<mEntity*> vertices;
  mVertex* vv;
  vv = part->createVertex(0,0.,0.,0.,part->getGEntity(0,3));
  vv->classify(part->getGEntity(0,0,fp));
  part->incTopo(FMDB_VERTEX);
  vertices.push_back((mEntity*)vv);
  vv = part->createVertex(1,1.,0.,0.,part->getGEntity(0,3));
  vv->classify(part->getGEntity(0,0,fp));
  part->incTopo(FMDB_VERTEX);
  vertices.push_back((mEntity*)vv);
  vv = part->createVertex(2,0.,1.,0.,part->getGEntity(0,3));
  vv->classify(part->getGEntity(0,0,fp));
  part->incTopo(FMDB_VERTEX);
  vertices.push_back((mEntity*)vv);
  vv = part->createVertex(3,0.,0.,1.,part->getGEntity(0,3));
  vv->classify(part->getGEntity(0,0,fp));
  part->incTopo(FMDB_VERTEX);
  vertices.push_back((mEntity*)vv);

  mEntity* ent;
  FMDB_Edge_Create(part, part->getGEntity(0,3), vertices[0], vertices[1],ent);
  ent->classify(part->getGEntity(0,1,fp));
  FMDB_Edge_Create(part, part->getGEntity(0,3), vertices[1], vertices[2],ent);
  ent->classify(part->getGEntity(0,1,fp));
  FMDB_Edge_Create(part, part->getGEntity(0,3), vertices[0], vertices[2],ent);
  ent->classify(part->getGEntity(0,1,fp));
  FMDB_Edge_Create(part, part->getGEntity(0,3), vertices[0], vertices[3],ent);
  ent->classify(part->getGEntity(0,1,fp));
  FMDB_Edge_Create(part, part->getGEntity(0,3), vertices[1], vertices[3],ent);
  ent->classify(part->getGEntity(0,1,fp));
  FMDB_Edge_Create(part, part->getGEntity(0,3), vertices[2], vertices[3],ent);
  ent->classify(part->getGEntity(0,1,fp));

  pMeshEnt pVtx[8];
  pVtx[0]=vertices[0];
  pVtx[1]=vertices[1];
  pVtx[2]=vertices[2];
  FMDB_Face_Create(part, part->getGEntity(0,3), FMDB_TRI, pVtx, 0, ent);
  ent->classify(part->getGEntity(0,2,fp));
  pVtx[0]=vertices[0];
  pVtx[1]=vertices[1];
  pVtx[2]=vertices[3];
  FMDB_Face_Create(part, part->getGEntity(0,3), FMDB_TRI, pVtx, 0, ent);
  ent->classify(part->getGEntity(0,2,fp));
  pVtx[0]=vertices[1];
  pVtx[1]=vertices[2];
  pVtx[2]=vertices[3];
  FMDB_Face_Create(part, part->getGEntity(0,3), FMDB_TRI, pVtx, 0, ent);
  ent->classify(part->getGEntity(0,2,fp));
  pVtx[0]=vertices[0];
  pVtx[1]=vertices[3];
  pVtx[2]=vertices[2];
  FMDB_Face_Create(part, part->getGEntity(0,3), FMDB_TRI, pVtx, 0, ent);
  ent->classify(part->getGEntity(0,2,fp));

  pVtx[0]=vertices[0];
  pVtx[1]=vertices[1];
  pVtx[2]=vertices[2];
  pVtx[3]=vertices[3];
  FMDB_Rgn_Create(part, part->getGEntity(0,3), FMDB_TET, 4, pVtx, ent);
  ent->classify(part->getGEntity(0,3,fp));

//  vector<mEntity*> vertices;
//  mVertex* vv;
//  vv = part->createVertex(0,0.,0.,0.,part->getGEntity(0,0,fp));
//  part->incTopo(FMDB_VERTEX);
//  vertices.push_back((mEntity*)vv);
//  vv = part->createVertex(1,1.,0.,0.,part->getGEntity(1,0,fp));
//  part->incTopo(FMDB_VERTEX);
//  vertices.push_back((mEntity*)vv);
//  vv = part->createVertex(2,0.,1.,0.,part->getGEntity(2,0,fp));
//  part->incTopo(FMDB_VERTEX);
//  vertices.push_back((mEntity*)vv);
//  vv = part->createVertex(3,0.,0.,1.,part->getGEntity(3,0,fp));
//  part->incTopo(FMDB_VERTEX);
//  vertices.push_back((mEntity*)vv);
//
//  mEntity* ent;
//  FMDB_Edge_Create(part, part->getGEntity(0,1,fp), vertices[0], vertices[1],ent);
//  FMDB_Edge_Create(part, part->getGEntity(1,1,fp), vertices[1], vertices[2],ent);
//  FMDB_Edge_Create(part, part->getGEntity(2,1,fp), vertices[0], vertices[2],ent);
//  FMDB_Edge_Create(part, part->getGEntity(3,1,fp), vertices[0], vertices[3],ent);
//  FMDB_Edge_Create(part, part->getGEntity(4,1,fp), vertices[1], vertices[3],ent);
//  FMDB_Edge_Create(part, part->getGEntity(5,1,fp), vertices[2], vertices[3],ent);
//
//  pMeshEnt pVtx[8];
//  pVtx[0]=vertices[0];
//  pVtx[1]=vertices[1];
//  pVtx[2]=vertices[2];
//  FMDB_Face_Create(part, part->getGEntity(0,2,fp), FMDB_TRI, pVtx, 0, ent);
//  pVtx[0]=vertices[0];
//  pVtx[1]=vertices[1];
//  pVtx[2]=vertices[3];
//  FMDB_Face_Create(part, part->getGEntity(1,2,fp), FMDB_TRI, pVtx, 0, ent);
//  pVtx[0]=vertices[1];
//  pVtx[1]=vertices[2];
//  pVtx[2]=vertices[3];
//  FMDB_Face_Create(part, part->getGEntity(2,2,fp), FMDB_TRI, pVtx, 0, ent);
//  pVtx[0]=vertices[0];
//  pVtx[1]=vertices[3];
//  pVtx[2]=vertices[2];
//  FMDB_Face_Create(part, part->getGEntity(3,2,fp), FMDB_TRI, pVtx, 0, ent);
//
//  pVtx[0]=vertices[0];
//  pVtx[1]=vertices[1];
//  pVtx[2]=vertices[2];
//  pVtx[3]=vertices[3];
//  FMDB_Rgn_Create(part, part->getGEntity(0,3,fp), FMDB_TET, 4, pVtx, ent);
//
  meshAdapt rdr(part,field,0,0);  // snap off; do refinement only
  //  rdr.setCallback(myCallback,0);
  L = 2;
  R = 0.45;
  rdr.run(atoi(argv[2]),1, setSizeField);
//  rdr.run(3,1, setSizeField);

  sprintf(outmesh,"%s-refined.sms",without_extension);
  FMDB_Mesh_WriteToFile(part->getMesh(), outmesh, 0);
//  M_writeSMS(mesh,outmesh,2);
  
//  delete field;
//  FMDB_Mesh_Del(mesh->getMesh());
//  M_delete(mesh);
//  GM_delete(model);
  MS_exit();

  return 1;
} 


