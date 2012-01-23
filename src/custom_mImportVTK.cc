/****************************************************************************** 

  (c) 2004-2011 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE-SCOREC file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <string.h>
#include "FMDB_Internals.h"
#include "oldFMDB.h"
#include "FMDB.h"
#include "SCUtil.h"
#include "mPart.h"
#include "mEntity.h"
#include "mFMDB.h"
#include "mFlexDB.h"

using std::cout;
using std::ofstream;
using std::ostream;
using std::istream;
using std::endl;
using std::vector;
 
void importTags(pMeshMdl mesh, vector<pMeshEnt> data, FILE* in);
void importLookUpTable(pMeshMdl mesh, FILE* in, int numPnts);
void custom_importLookUpTable(pMeshMdl, FILE*, vector<pMeshEnt>, int);
// **********************************************************
int custom_importVTK(mMesh *mesh, const char *fName)
// **********************************************************
{
  FILE *in = fopen (fName,"r");
   if(!in)
    return SCUtil_FILE_NOT_FOUND; 
  int numStrPerLine[3] = {5, 1, 2};
  char dummy_str[64];
  char entity_type[64];
  bool point_flag=false, cell_flag=false;
  double x,y,z;
  int i, j, num_cells, num_cells2, num_pnts, num_cell_data;  
  int edges_count=0, faces_count=0, regs_count=0, cells_count=0;
  int cell_type, mDim;

  pPart part;
  FMDB_Mesh_GetPart(mesh, 0, part); 
  FMDB_Part_GetDim(part, &mDim);
#ifdef FLEXDB 
  set_MRM_oneLevel(part,3);
#endif
  
  while(1)
  {
     fscanf(in,"%s",dummy_str);
     if (strcmp(dummy_str, "DATASET")==0 )
     {
       break; 
     }
  }

  fscanf(in,"%s",dummy_str);
  if (strcmp(dummy_str, "UNSTRUCTURED_GRID")==0 )
     point_flag=true; 

  vector<mEntity*> vertices;    
  vector<pMeshEnt> faceVec;
  vector<pMeshEnt> edgeVec;
  vector<pMeshEnt> regVec;
  vector<pMeshEnt> cellVec;

  int vt1, vt2;   
  mVertex* vv;        
  mEntity* ent;
  list<int*> list_cells; 
  mEntity* tmpEdges[4]; 
  int* int_vids; 

  pMeshEnt pVtx[8];

  if(point_flag)
  {
   fscanf(in,"%s",entity_type);
   if (strcmp(entity_type, "POINTS")==0)
    {
      fscanf(in,"%d %s", &num_pnts, dummy_str);

      for (int i=0; i<num_pnts; ++i)
      {
        fscanf(in,"%lf %lf %lf", &x, &y, &z);
        vv = part->createVertex(i+1,x,y,z, part->getGEntity(0,3)); 
	part->incTopo(FMDB_VERTEX);
        vertices.push_back((mEntity*)vv);
      }
   }
   
   fscanf(in,"%s",entity_type);
   if (strcmp(entity_type, "CELLS")==0)
    {
      fscanf(in,"%d %d", &num_cells, &num_cell_data); 
      for (j=0; j<num_cells; ++j) {
	fscanf(in,"%d ", &num_pnts);
	
	int_vids = new int[num_pnts]; 
        for(int icount=0; icount<num_pnts; ++icount){
	  fscanf(in, "%d", &int_vids[icount]); 
	}
	list_cells.push_back(int_vids); 
      }    	
     }
  
    fscanf(in,"%s %d",entity_type, &num_cells2);
    assert(num_cells==num_cells2);
    for(j=0; j<num_cells; ++j)
    {
      if(!list_cells.empty())
      {
	 int_vids = list_cells.front(); 
	 list_cells.pop_front(); 
      }
      else 
      {
	 break; 
	 cout<<"Error in loading vtk files!"<<endl; 
      }
	 
      fscanf(in, "%d", &cell_type); 
      switch(cell_type)
      {
	case 3: 	 
	 vt1 = int_vids[0]; 
	 vt2 = int_vids[1]; 
         FMDB_Edge_Create(part, part->getGEntity(0,3), vertices[int_vids[0]], vertices[int_vids[1]],ent);
	 edgeVec.push_back(ent);
	 edges_count++; 
	 break; 
	
	case 5: // triangle
          for(int iV=0; iV<3; iV++)
           pVtx[iV]=vertices[int_vids[iV]];
          
	  FMDB_Face_Create(part, part->getGEntity(0,3), FMDB_TRI, pVtx, 0, ent); 
	  faceVec.push_back(ent);
	  faces_count++; 
	  break;
	  
        case 9: // quad

	 for(int iV=0; iV<4; iV++)
           pVtx[iV]=vertices[int_vids[iV]];

          FMDB_Face_Create(part, part->getGEntity(0,3), FMDB_QUAD, pVtx, 0, ent);
	  faceVec.push_back(ent);
          faces_count++; 
	  break;
        
	case 10: // tet 

  	 for(int iV=0; iV<4; iV++)
           pVtx[iV]=vertices[int_vids[iV]];
         FMDB_Rgn_Create(part, part->getGEntity(0,3), FMDB_TET, 4, pVtx, ent);
	 regVec.push_back(ent);
	 regs_count++;  
	 break; 	  
 
	case 12: // hex

	 for(int iV=0; iV<8; iV++)
           pVtx[iV]=vertices[int_vids[iV]];
         FMDB_Rgn_Create(part, part->getGEntity(0,3), FMDB_HEX, 8, pVtx, ent);
	 regVec.push_back(ent);
	 regs_count++; 
         break; 
	  
	case 13: // wedge
	  for(int iV=0; iV<6; iV++)
           pVtx[iV]=vertices[int_vids[iV]];
          FMDB_Rgn_Create(part, part->getGEntity(0,3), FMDB_PRISM, 6, pVtx, ent);
	  regVec.push_back(ent); 
          regs_count++; 
	  break; 
	 
	case 14: // pyramid
          for(int iV=0; iV<5; iV++)
           pVtx[iV]=vertices[int_vids[iV]];
          FMDB_Rgn_Create(part, part->getGEntity(0,3), FMDB_PYRAMID, 5, pVtx, ent);	
	  regVec.push_back(ent);
	  regs_count++; 
	  break; 

	default: 
         {
	  cout<<"FMDB Warning: the element type is not supported now.\n";
	  return SCUtil_INVALID_ENTITY_TYPE;
         }
     } // switch
   cells_count++;
   cellVec.push_back(ent);
   delete[] int_vids; 

   } // for j
   char pt_str[64], att_str[64];
   while(fscanf(in,"%s",pt_str) != EOF)
   {
   if (strcmp(pt_str, "POINT_DATA")!=0 && strcmp(pt_str, "CELL_DATA")!=0)
    return SCUtil_SUCCESS;

   int data_count;
   fscanf(in,"%d", &data_count);

   fscanf(in, "%s", dummy_str);
   strcpy(att_str, dummy_str);
   // right now, the support is only for field attributes
   if(strcmp(att_str, "FIELD")!=0 && 
      strcmp(att_str, "SCALARS")!=0 && 
      strcmp(att_str, "LOOKUP_TABLE") !=0)
    return SCUtil_SUCCESS;

   if(strcmp(pt_str, "POINT_DATA")==0)
    {
    if(strcmp(att_str, "SCALARS")==0)
     importLookUpTable(mesh, in, vertices.size());
    if(strcmp(att_str, "FIELD")==0)
     importTags(mesh, vertices, in);
    }
   else if(strcmp(pt_str, "CELL_DATA")==0)
    {
    if(strcmp(att_str, "SCALARS")==0)
//     importLookUpTable(mesh, in, regs_count);
     custom_importLookUpTable(mesh, in, cellVec, cells_count);
    if(strcmp(att_str, "FIELD")==0)
    {
    if(mDim==2)
     importTags(mesh, faceVec, in);
    else
     importTags(mesh, regVec, in);
    }
    }
  } 
 } // point_flag  
 fclose (in);
 return SCUtil_SUCCESS;
}
void custom_importLookUpTable(pMeshMdl mesh, FILE* in, vector<pMeshEnt> entVec, int numPnts)
{
 char scalar_type[64], scalar_name[64], dummy_str[64];
 fscanf(in, "%s %s", scalar_name, scalar_type);
 int tag_size=1, tag_type, nIntTags=0;
 pTag new_tag;
 // identify tag type
 if((strcmp(scalar_type, "int")==0) ||
    (strcmp(scalar_type, "unsigned_int")==0) ||
    (strcmp(scalar_type, "unsigned_short")==0) ||
    (strcmp(scalar_type, "short")==0))
    tag_type = SCUtil_INT;
 else if(strcmp(scalar_type, "float")==0 || strcmp(scalar_type, "double")==0)
    tag_type=SCUtil_DBL;
 else if(strcmp(scalar_type, "byte")==0)
    tag_type=SCUtil_BYTE;
 else if(strcmp(scalar_type, "entity")==0)
    tag_type=SCUtil_ENT;
 else if(strcmp(scalar_type, "entity_set")==0)
    tag_type=SCUtil_SET;
 else
    // data type not supported
   return;
 if (FMDB_Mesh_CreateTag (mesh, scalar_name, tag_type, tag_size, new_tag))
       FMDB_Mesh_FindTag (mesh, scalar_name, new_tag);
 if(fscanf(in,"%s",dummy_str) == EOF)
  return;

 if(strcmp(dummy_str, "LOOKUP_TABLE")==0)
 {
  fscanf(in, "%s", dummy_str);
  if(strcmp(scalar_type, "int")==0 || strcmp(dummy_str, "default")==0)
  {
   int int_data;
   for(int i=0; i<numPnts; i++)
   {
    fscanf(in, "%d", &int_data);
    FMDB_Ent_SetIntTag (mesh, entVec[i], new_tag, int_data);
    nIntTags++;
   }
  }
  else if(strcmp(scalar_type, "float")==0)
   {
    float float_data;
    for(int i=0; i<numPnts; i++)
     fscanf(in, "%f", &float_data);
   }
 }
 printf("nIntTags = %d\n", nIntTags);
}
void importLookUpTable(pMeshMdl mesh, FILE* in, int numPnts)
{
 char scalar_type[64], scalar_name[64],dummy_str[64];
 fscanf(in, "%s %s", scalar_name, scalar_type);
 if(fscanf(in,"%s",dummy_str) == EOF)
  return;

 if(strcmp(dummy_str, "LOOKUP_TABLE")==0)
 {
  fscanf(in, "%s", dummy_str);
  if(strcmp(scalar_type, "int")==0 || strcmp(dummy_str, "default")==0)
  {
   int int_data;
   for(int i=0; i<numPnts; i++)
    fscanf(in, "%d", &int_data);
  }
  else if(strcmp(scalar_type, "float")==0)
   {
    float float_data;
    for(int i=0; i<numPnts; i++)
     fscanf(in, "%f", &float_data);
   }
 }
}
void importTags(pMeshMdl mesh, vector<pMeshEnt> data, FILE* in)
{
  int numTag, opq_data_size;
  char dummy_str[64];
  float tmp_data;
  std::vector<pTag> tag_vector;
  fscanf(in,"%s %d",dummy_str, &numTag);
  
  for (int i=0; i<numTag; i++) //restore tags
  {
   int tag_num_pnts, tag_size, tag_type;
   pTag new_tag;
   char tagName[256], tagTypeChar[24];
   fscanf(in, "%s %d %d %s", tagName, &tag_size, &tag_num_pnts, tagTypeChar);
   /*if(tag_num_pnts!=data.size())
   {
    cout << "\n FMDB ERROR: Incorrect tagging information. Check VTK file format." << endl;
    return;
   }*/

   // identify tag type
   if((strcmp(tagTypeChar, "int")==0) ||
      (strcmp(tagTypeChar, "unsigned_int")==0) ||
      (strcmp(tagTypeChar, "unsigned_short")==0) ||
      (strcmp(tagTypeChar, "short")==0))
      tag_type = SCUtil_INT;
   else if(strcmp(tagTypeChar, "float")==0 || strcmp(tagTypeChar, "double")==0)
      tag_type=SCUtil_DBL;
   else if(strcmp(tagTypeChar, "byte")==0)
      tag_type=SCUtil_BYTE;
   else if(strcmp(tagTypeChar, "entity")==0)
      tag_type=SCUtil_ENT;
   else if(strcmp(tagTypeChar, "entity_set")==0)
      tag_type=SCUtil_SET;
   else
      // data type not supported
     return;
    
  if (FMDB_Mesh_CreateTag (mesh, tagName, tag_type, tag_size, new_tag))
        FMDB_Mesh_FindTag (mesh, tagName, new_tag);
  tag_vector.push_back(new_tag);
     
  switch(tag_type)
  {
   case SCUtil_BYTE:
   {
    char* opq_data=(char*)calloc(256, sizeof(char));
    for(int i=0; i<tag_num_pnts; i++)
    {
     fscanf(in, "%d %s", &opq_data_size, opq_data);
     FMDB_Ent_SetByteTag (mesh, data[i], new_tag, opq_data, opq_data_size);
    }
    free(opq_data);
   }
   case SCUtil_INT://integer
   {
    if(tag_size==1)
    {
     for(int i=0; i<tag_num_pnts; i++)
     {
      int int_data;
      if(strcmp(tagTypeChar, "unsigned_int") ==0)
      {
       unsigned int ui_data;
       fscanf(in, "%u", &ui_data);
       int_data=ui_data; // implicit conversion
      } 
      else if(strcmp(tagTypeChar, "short") ==0)
      {
       short int si_data;
       fscanf(in, "%hd", &si_data);
       int_data=si_data; //implicit conversion
      }
      else if(strcmp(tagTypeChar, "unsigned_short") ==0)
      {
       unsigned short us_data;
       fscanf(in, "%hu", &us_data);
       int_data=us_data;//implicit conversion
      }
     else if(strcmp(tagTypeChar, "int") ==0)
      fscanf(in, "%d", &int_data);
     
     FMDB_Ent_SetIntTag (mesh, data[i], new_tag, int_data);
     }
   }
  else
   {
   for(int i=0; i<tag_num_pnts; i++)
   {
    int* arr_data=new int[tag_size];
    for (int j=0; j<tag_size; j++)
    {
    if(strcmp(tagTypeChar, "unsigned_int") ==0)
     {
      unsigned int ui_data;
      fscanf(in, "%u", &ui_data);
      arr_data[j]=ui_data;
     }
   else if(strcmp(tagTypeChar, "short") ==0)
    {
    short int si_data;
    fscanf(in, "%hd", &si_data);
    arr_data[j]=si_data;
    }
   else if(strcmp(tagTypeChar, "unsigned_short") ==0)
    {
     unsigned short us_data;
     fscanf(in, "%hu", &us_data);
     arr_data[j]=us_data;
    }
   else if(strcmp(tagTypeChar, "int") ==0)
     fscanf(in, "%d", &arr_data[j]);
   }
   FMDB_Ent_SetIntArrTag (mesh, data[i], new_tag, arr_data, tag_size);	
   delete[] arr_data;
  }//end for
 }//end else
 }//end case
 break;
 case SCUtil_DBL://double
 {
 if(tag_size==1)
 {
  double dbl_data;
  if(strcmp(tagTypeChar, "float") ==0)
  {
   for(int i=0; i<tag_num_pnts; i++)
    {
     float float_data;
     fscanf(in, "%f", &float_data);
     dbl_data=float_data; //implicit conversion
    }
  }
 else if(strcmp(tagTypeChar, "double") ==0)
  {
   for(int i=0; i<tag_num_pnts; i++)
    fscanf(in, "%lf", &dbl_data);
  }
 FMDB_Ent_SetDblTag (mesh, data[i], new_tag, dbl_data);
 }//end if tag_size
 else
 {
  for(int i=0; i<tag_num_pnts; i++)
  {
   double* arr_data=new double[tag_size];
   if(strcmp(tagTypeChar, "float") ==0)
   {
    float float_data;
    for (int j=0; j<tag_size; j++)
    {
      fscanf(in, "%f", &float_data);
      arr_data[j]=float_data;
    }
   }
   else if(strcmp(tagTypeChar, "double") ==0)
   {
   for (int j=0; j<tag_size; j++)
    fscanf(in, "%lf", &arr_data[j]);
   }
   FMDB_Ent_SetDblArrTag (mesh, data[i], new_tag, arr_data, tag_size);
   delete[] arr_data;
  }// end for
 }//end else
 }
  break;
  case SCUtil_ENT:
  {
   // For now, not required for vtk file format
  }
  break;
  case SCUtil_SET:
  {
   // For now, not required for vtk file format
  }
  break;
  default:
   cout << "\n FMDB Warning: Field data type not supported" << endl;
  }//end switch
 }// end for
}
