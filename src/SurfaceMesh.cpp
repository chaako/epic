/*
 * SurfaceMesh.cpp
 *
 *  Created on: Dec 3, 2012
 *      Author: chaako
 */

#include "SurfaceMesh.h"

SurfaceMesh::SurfaceMesh() {
	vtkMesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
	volumeMesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
}

SurfaceMesh::SurfaceMesh(string inputFile) {
	vtkMesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
	volumeMesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
	this->load(inputFile);
}

SurfaceMesh::~SurfaceMesh() {
	// TODO Auto-generated destructor stub
}

void SurfaceMesh::load(string inputFile) {
	if (inputFile.find(".vtu")==string::npos)
		throw; // Not a .vtu file
	vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
			vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
	reader->SetFileName(inputFile.c_str());
	reader->Update();
//	vtkMesh->Print(cout);
	vtkMesh->DeepCopy(reader->GetOutput());
}

void SurfaceMesh::save(string outputFile) {
	vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
			vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	writer->SetFileName(outputFile.c_str());
	writer->SetInput(vtkMesh);
	writer->Write();
}

void SurfaceMesh::saveVolumeMesh(string outputFile) {
	vtkSmartPointer<vtkUnstructuredGridWriter> writer =
			vtkSmartPointer<vtkUnstructuredGridWriter>::New();
	stringstream outputFilename;
	writer->SetFileName(outputFile.c_str());
	writer->SetInput(volumeMesh);
	writer->Write();
}

void SurfaceMesh::transformSurface(int surfaceId, vect3d origin,
		vect3d rotationAxis, double rotationAngle, vect3d scaleFactors,
		vect3d translation) {
	vector<bool> transformNode(vtkMesh->GetNumberOfPoints(), false);
	string cell_code_name = "cell_code";
	vtkIntArray* orig_cell_codes = vtkIntArray::SafeDownCast(
			vtkMesh->GetCellData()->GetArray(cell_code_name.c_str()));
	for (vtkIdType id_cell = 0; id_cell < vtkMesh->GetNumberOfCells(); ++id_cell) {
		int cell_code = orig_cell_codes->GetValue(id_cell);
		if (cell_code==surfaceId) {
			vtkIdType *pts, N_pts;
			vtkMesh->GetCellPoints(id_cell, N_pts, pts);
			for (int i = 0; i < 3; ++i) {
				transformNode[pts[i]] = true;
			}
		}
	}
	for (vtkIdType id_node = 0; id_node < vtkMesh->GetNumberOfPoints(); ++id_node) {
		if (transformNode[id_node]) {
			double x[NDIM];
			vtkMesh->GetPoints()->GetPoint(id_node, x);
			vect3d pos(x[0],x[1],x[2]);
			pos -= origin;
			pos = Eigen::AngleAxisd(-rotationAngle, rotationAxis)*pos;
			for (int i=0; i<NDIM; i++) {
				pos[i] *= scaleFactors[i];
				pos[i] += origin[i];
				pos[i] += translation[i]*scaleFactors[i];
				x[i] = pos[i];
			}
			vtkMesh->GetPoints()->SetPoint(id_node, x);
		}
	}
}

void SurfaceMesh::scaleVolumeMesh(vect3d origin, vect3d scaleFactors) {
	for (vtkIdType id_node = 0; id_node < volumeMesh->GetNumberOfPoints(); ++id_node) {
		double x[NDIM];
		volumeMesh->GetPoints()->GetPoint(id_node, x);
		vect3d pos(x[0],x[1],x[2]);
		pos -= origin;
		for (int i=0; i<NDIM; i++) {
			pos[i] *= scaleFactors[i];
			pos[i] += origin[i];
			x[i] = pos[i];
		}
		volumeMesh->GetPoints()->SetPoint(id_node, x);
	}
}

void SurfaceMesh::createVolumeMesh() {
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
	if (res == NG_OK) {
		int Npoints_ng = Ng_GetNP(ngMesh);
		int Ncells_ng  = Ng_GetNE(ngMesh);
		int Nscells_ng = Ng_GetNSE(ngMesh);
		vtkIdType new_point = 0;

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
					int cell_code = orig_cell_codes->GetValue(id_cell);
					cell_codes->InsertNextValue(cell_code);
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
}

vector<vect3d> SurfaceMesh::getPoints(int surfaceCode) {
	vector<vect3d> points;
	string cell_code_name = "cell_code";
	vtkIntArray* cell_codes = vtkIntArray::SafeDownCast(
			vtkMesh->GetCellData()->GetArray(cell_code_name.c_str()));
	int numberOfVtkPoints = vtkMesh->GetNumberOfPoints();
	vector<bool> processedPoint(numberOfVtkPoints);
	for (int i=0; i<processedPoint.size(); i++)
		processedPoint[i] = false;
	for (vtkIdType id_cell = 0; id_cell < vtkMesh->GetNumberOfCells(); ++id_cell) {
		if (vtkMesh->GetCellType(id_cell) == VTK_TRIANGLE) {
			// TODO: not very clean to use surfaceCode<0 to mean any code
			if (cell_codes->GetValue(id_cell)==surfaceCode || surfaceCode<0) {
				// TODO: figure out if supposed to free pts
				vtkIdType N_pts, *pts;
				vtkMesh->GetCellPoints(id_cell, N_pts, pts);
				for (int i = 0; i < N_pts; ++i) {
					if (processedPoint[pts[i]] == false) {
						double coords[3];
						vtkMesh->GetPoint(pts[i],coords);
						vect3d coordinates(coords[0],coords[1],coords[2]);
						points.push_back(coordinates);
						processedPoint[pts[i]] = true;
					}
				}
			}
		}
	}
	return points;
}
