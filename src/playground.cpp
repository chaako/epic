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


int playground(int argc, char* argv[])
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


