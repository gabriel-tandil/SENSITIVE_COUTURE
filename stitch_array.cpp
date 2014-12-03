/*
 *  stitch_array.cpp
 *  sensitive couture
 *
 *  Created by Nobuyuki Umetani on 12/20/10.
 *  Copyright 2010 The University of Tokyo and Columbia University. All rights reserved.
 *
 */


#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem/field.h"
#include "stitch_array.h"
#include "eqn_glue.h"

using namespace Fem::Field;



////////////////

bool CStitchAry::IsSeamLine(unsigned int id_e_in, unsigned int& id_e_out, bool& is_same_dir) const
{
  for(unsigned int ist=0;ist<aStitch.size();ist++){
    unsigned int id_e1 = aStitch[ist].id_e1;
    unsigned int id_e2 = aStitch[ist].id_e2;
    if( id_e1 == id_e_in ){
      id_e_out = id_e2;
      is_same_dir = aStitch[ist].is_same_dir;
      return true;
    }
    if( id_e2 == id_e_in ){
      id_e_out = id_e1;
      is_same_dir = aStitch[ist].is_same_dir;
      return true;
    }    
  }
  return false;
}

void CStitchAry::CopyCADSituation(const CStitchAry& sa)
{
//  std::cout << "CopyCADSituation " << sa.aStitch.size() << std::endl;
  this->aStitch.clear();
  for(unsigned int ist=0;ist<sa.aStitch.size();ist++){
    unsigned int id_e1 = sa.aStitch[ist].id_e1;
    unsigned int id_e2 = sa.aStitch[ist].id_e2;
//    std::cout << "CAD edge stitch: " << id_e1 << " " << id_e2 << std::endl;
    bool dir = sa.aStitch[ist].is_same_dir;
    aStitch.push_back( CStitch(id_e1,id_e2,dir) );    
  }
  this->stitch_stiff = sa.stitch_stiff;
  this->stitch_damp_coeff = sa.stitch_damp_coeff;
}

void CStitchAry::AddStitch(const Cad::CCadObj2D& cad_2d, const Msh::CMesher2D& mesh_2d, unsigned int id_e1, unsigned int id_e2)
{  
  bool is_left1 = true;
  {
    unsigned int id_l1_l,id_l1_r;
    cad_2d.GetIdLoop_Edge(id_l1_l,id_l1_r,id_e1);
    if( !mesh_2d.IsIdLCad_CutMesh(id_l1_l) ){ id_l1_l = 0; }
    if( !mesh_2d.IsIdLCad_CutMesh(id_l1_r) ){ id_l1_r = 0; }    
    assert( id_l1_l*id_l1_r == 0 );    
    if(      id_l1_l == 0 ){ is_left1 = false; }
    else if( id_l1_r == 0 ){ is_left1 = true;  }
  }
  bool is_left2 = true;
  {
    unsigned int id_l2_l,id_l2_r;
    cad_2d.GetIdLoop_Edge(id_l2_l,id_l2_r,id_e2);
    if( !mesh_2d.IsIdLCad_CutMesh(id_l2_l) ){ id_l2_l = 0; }
    if( !mesh_2d.IsIdLCad_CutMesh(id_l2_r) ){ id_l2_r = 0; }        
    assert( id_l2_l*id_l2_r == 0 );
    if(      id_l2_l == 0 ){ is_left2 = false; }
    else if( id_l2_r == 0 ){ is_left2 = true;  }
  }
  bool dir = ( is_left1 != is_left2 );    
  if( dir ){
    std::cout << "same dir " << id_e1 << " " << id_e2 << std::endl;
  }
  else{
    std::cout << "opp dir " << id_e1 << " " << id_e2 << std::endl;    
  }
  aStitch.push_back( CStitch(id_e1,id_e2,dir) );
}

void CStitchAry::AddLinSys_BackwardEular(double dt, double& se,                                
                                         Fem::Ls::CLinearSystem_Field& ls,
                                         unsigned int id_field_disp,
                                         const Fem::Field::CFieldWorld& world) const
{  
  for(unsigned int ist=0;ist<aStitch.size();ist++)
  {
    AddLinSys_Glut_Penalty_BackwardEular
    (dt, stitch_damp_coeff,
     ls,
     stitch_stiff,
     id_field_disp, aStitch[ist].id_field,
     world,
     se);
  }		
}

void CStitchAry::AddLinSys_Sensitivity(Fem::Ls::CLinearSystem_Field& ls,
                                       unsigned int id_field_disp,
                                       const Fem::Field::CFieldWorld& world) const
{
  for(unsigned int idart=0;idart<aStitch.size();idart++)
  {
    AddLinSys_Glut_Penalty_Sensitivity
    (ls,
     stitch_stiff,
     id_field_disp, aStitch[idart].id_field,
     world );
  }	  
}

void CStitchAry::MakeField(unsigned int id_field_disp, unsigned int id_base, Fem::Field::CFieldWorld& world){ 
  const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);    
  for(unsigned int ist=0;ist<aStitch.size();ist++){
    unsigned int id_e1 = aStitch[ist].id_e1;
    unsigned int id_e2 = aStitch[ist].id_e2;
    bool dir = aStitch[ist].is_same_dir;
    ////      
    bool is_reordering;
    unsigned int id_ea1 = conv.GetIdEA_fromCad(id_e1,Cad::EDGE);
    unsigned int id_ea2 = conv.GetIdEA_fromCad(id_e2,Cad::EDGE);
    unsigned int id_field_dart = MakePartialField_GlueEdge_Penalty
    (world,id_field_disp,
     id_ea1,id_ea2,dir,
     is_reordering);
    ////
    aStitch[ist].id_field = id_field_dart;
    aStitch[ist].need_reordering = is_reordering;
  }
}

void CStitchAry::AddStitchField(const Cad::CCadObj2D& cad_2d, const Msh::CMesher2D& mesh_2d, unsigned int id_e1, unsigned int id_e2,
                                unsigned int id_field_disp, unsigned int id_base, Fem::Field::CFieldWorld& world)
{
  assert( cad_2d.IsElemID(Cad::EDGE,id_e1) );
  assert( cad_2d.IsElemID(Cad::EDGE,id_e2) );
  assert( mesh_2d.GetElemID_FromCadID(id_e1,Cad::EDGE) != 0 );
  assert( mesh_2d.GetElemID_FromCadID(id_e2,Cad::EDGE) != 0 );    
  unsigned int ist0 = aStitch.size();
  this->AddStitch(cad_2d,mesh_2d,id_e1,id_e2);
  const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);    
  assert( id_e1 == aStitch[ist0].id_e1 );
  assert( id_e2 == aStitch[ist0].id_e2 );
  bool dir = aStitch[ist0].is_same_dir;
  ////      
  bool is_reordering;
  unsigned int id_ea1 = conv.GetIdEA_fromCad(id_e1,Cad::EDGE);
  unsigned int id_ea2 = conv.GetIdEA_fromCad(id_e2,Cad::EDGE);
  assert( world.IsIdEA(id_ea1) );
  assert( world.IsIdEA(id_ea2) );  
  unsigned int id_field_dart = MakePartialField_GlueEdge_Penalty
  (world,id_field_disp,
   id_ea1,id_ea2,dir,
   is_reordering);
  ////
  aStitch[ist0].id_field = id_field_dart;
  aStitch[ist0].need_reordering = is_reordering;                  
}

void CStitchAry::RemoveStitchField(unsigned int id_e1_in, unsigned int id_e2_in,
                                   unsigned int id_field_disp, unsigned int id_base, Fem::Field::CFieldWorld& world)
{
  int istitch = -1;
  for(unsigned int ist=0;ist<aStitch.size();ist++){
    unsigned int id_e1 = aStitch[ist].id_e1;
    unsigned int id_e2 = aStitch[ist].id_e2;
    if( (id_e1==id_e1_in && id_e2==id_e2_in ) || (id_e1==id_e2_in && id_e2==id_e1_in ) ){
      istitch = ist;
      std::vector<unsigned int> aIdFieldDel;
      aIdFieldDel.push_back(aStitch[ist].id_field);
      world.DeleteField(aIdFieldDel);      
    }
  }      
  if( istitch != -1 ){
    aStitch.erase(aStitch.begin()+istitch);
  }                   
}


void CStitchAry::AddPatternLinearSystem
(Fem::Ls::CLinearSystem_Field& ls, const Fem::Field::CFieldWorld& world)
{
  for(unsigned int idart=0;idart<aStitch.size();idart++){
    unsigned int id_field_dart = aStitch[idart].id_field;
    ls.AddPattern_Field(id_field_dart,world);
  }    
}

void DrawSeamLine(unsigned int id_field_disp, unsigned int id_field_dart, const CFieldWorld& world)
{
  const CField& field = world.GetField(id_field_disp);
	const Fem::Field::CNodeAry::CNodeSeg& ns_c = field.GetNodeSeg(CORNER,false,world,VALUE);
	const Fem::Field::CNodeAry::CNodeSeg& ns_u = field.GetNodeSeg(CORNER,true, world,VALUE);
	const Fem::Field::CNodeAry::CNodeSeg& ns_v = field.GetNodeSeg(CORNER,true, world,VELOCITY);
  const CField& field_dart = world.GetField(id_field_dart);
  const Fem::Field::CNodeAry::CNodeSeg& ns_r = field_dart.GetNodeSeg(BUBBLE,false,world,VALUE);
	const unsigned int id_ea = field_dart.GetAryIdEA()[0];
  const CElemAry& ea = world.GetEA(id_ea);
	const unsigned int nelem = ea.Size();  
	const Fem::Field::CElemAry::CElemSeg& es_c = field_dart.GetElemSeg(id_ea,CORNER,true,world);
	const Fem::Field::CElemAry::CElemSeg& es_b = field_dart.GetElemSeg(id_ea,BUBBLE,false,world);  
  bool is_lighting = ::glIsEnabled(GL_LIGHTING);
  ::glDisable(GL_LIGHTING);
  ::glLineWidth(1);
  ::glColor3d(0,0,1);
  ::glBegin(GL_LINES);
	for(unsigned int ielem=0;ielem<nelem;ielem++){
		unsigned int no[3];	es_c.GetNodes(ielem,no);
		double C[3][3];	ns_c.GetValue(no[0],C[0]);	ns_c.GetValue(no[1],C[1]);  ns_c.GetValue(no[2],C[2]);
		double u[3][3];	ns_u.GetValue(no[0],u[0]);	ns_u.GetValue(no[1],u[1]);  ns_u.GetValue(no[2],u[2]);
		double v[3][3];	ns_v.GetValue(no[0],v[0]);	ns_v.GetValue(no[1],v[1]);  ns_v.GetValue(no[2],v[2]);
    unsigned int ino_b; es_b.GetNodes(ielem,&ino_b);
    double r; ns_r.GetValue(ino_b,&r);
    const double c[3][3] = {
      { C[0][0]+u[0][0], C[0][1]+u[0][1], C[0][2]+u[0][2] },
      { C[1][0]+u[1][0], C[1][1]+u[1][1], C[1][2]+u[1][2] },    
      { C[2][0]+u[2][0], C[2][1]+u[2][1], C[2][2]+u[2][2] } };
    double m[3] = {
      c[1][0]*(1-r)+c[2][0]*r,
      c[1][1]*(1-r)+c[2][1]*r,
      c[1][2]*(1-r)+c[2][2]*r };
    ::glVertex3dv(c[0]);
    ::glVertex3dv(m);
  }
  ::glEnd();
  if( is_lighting ){ ::glEnable(GL_LIGHTING); }
}



void CStitchAry::Draw(const unsigned int id_field_disp,const Fem::Field::CFieldWorld& world) const
{    
  bool is_lighting =::glIsEnabled(GL_LIGHTING);
  ::glDisable(GL_LIGHTING);
  bool is_texture = ::glIsEnabled(GL_TEXTURE_2D);
  ::glDisable(GL_TEXTURE_2D);    
  for(unsigned int istitch=0;istitch<aStitch.size();istitch++){
    unsigned int id_field_stitch = aStitch[istitch].id_field;
    DrawSeamLine(id_field_disp, id_field_stitch, world);
  }    
  if( is_lighting ){ ::glEnable(GL_LIGHTING); }
  if( is_texture  ){ ::glEnable(GL_TEXTURE_2D); }    
}

void CStitchAry::SetReorderingPreconditionerIfNeeded
(LsSol::CPreconditioner_ILU& prec, 
 unsigned int id_field_disp,          
 const Fem::Field::CFieldWorld& world) const
{
//  return;
  bool is_reordering = false;
  for(unsigned int istitch=0;istitch<aStitch.size();istitch++){
    if( aStitch[istitch].need_reordering ){ is_reordering = true; break; }
  }
  if( !is_reordering ){
    std::vector<int> ord0;
    prec.SetOrdering(ord0);
    return;
  }
  ////
  const CField& field = world.GetField(id_field_disp);
  const unsigned int nno = field.GetNodeSeg(CORNER,true,world,VALUE).Size();
  std::set<unsigned int> set1;
  std::set<unsigned int> set2;
  unsigned int min = nno+1;
  unsigned int max = 0;
  for(unsigned int istitch=0;istitch<aStitch.size();istitch++){
    if( !aStitch[istitch].need_reordering ) continue;
    const CField& field_st = world.GetField(aStitch[istitch].id_field);
    unsigned int id_ea = field_st.GetAryIdEA()[0];
    const CElemAry::CElemSeg& es1 = field.GetElemSeg(id_ea,CORNER,true,world);    
    for(unsigned int ie=0;ie<es1.Size();ie++){
      unsigned int no[3]; es1.GetNodes(ie,no);
      set1.insert(no[0]);
      set2.insert(no[1]);
      set2.insert(no[2]);
      min = ( no[0] < min ) ? no[0] : min;
      min = ( no[1] < min ) ? no[1] : min;
      min = ( no[2] < min ) ? no[2] : min;
      max = ( no[0] > max ) ? no[0] : max;
      max = ( no[1] > max ) ? no[1] : max;
      max = ( no[2] > max ) ? no[2] : max;
    }
  }
  std::vector<int> ord0;
  ord0.resize(nno);
  for(unsigned int ino=0;ino<nno;ino++){ ord0[ino] = ino; }
  for(unsigned int ino=min;ino<=max;ino++){ ord0[ino] = -3; }
  for(std::set<unsigned int>::iterator itr=set2.begin();itr!=set2.end();itr++){ ord0[*itr] = -2; }
  for(std::set<unsigned int>::iterator itr=set1.begin();itr!=set1.end();itr++){ 
    unsigned int ipo1 = *itr;
    std::set<unsigned int>::iterator itr2 = set2.find(ipo1);
    if( itr2 != set2.end() ){
      set2.erase(itr2);
    }
    ord0[ipo1] = -1; 
  }
  std::vector<unsigned int> vec3; vec3.reserve(nno);
  for(unsigned int ino=min;ino<=max;ino++){ 
    if( ord0[ino] == -3 ){ vec3.push_back(ino); }  
  }  
  //  std::cout << vec3.size() << " " << max-min+1-vec1.size()-vec2.size() << std::endl;
  assert( max-min+1 == set1.size()+set2.size()+vec3.size() );
  unsigned int icnt = min;
  for(std::set<unsigned int>::iterator itr=set1.begin();itr!=set1.end();itr++){ 
    ord0[icnt] = *itr;
    icnt++;
  }  
  for(std::set<unsigned int>::iterator itr=set2.begin();itr!=set2.end();itr++){ 
    ord0[icnt] = *itr;
    icnt++;
  }    
  for(unsigned int i=0;i<vec3.size();i++){ ord0[min+set1.size()+set2.size()+i] = vec3[i]; }
  for(unsigned int i=min;i<=max;i++){
    //    std::cout << i << " " << ord0[i] << std::endl;
  }
  prec.SetOrdering(ord0);
}

void CStitchAry::DrawBoundaryCondition2D(const Cad::CCadObj2D& cad_2d) const
{
  bool is_lighting =::glIsEnabled(GL_LIGHTING);
  ::glDisable(GL_LIGHTING);
  bool is_texture =::glIsEnabled(GL_TEXTURE_2D);
  ::glDisable(GL_TEXTURE_2D);  
  ::glEnable(GL_LINE_SMOOTH);
  ::glEnable(GL_BLEND);
  ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  ::glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
  ::glColor3d(1,0,0);
  ::glLineWidth(2);
  for(unsigned int idart=0;idart<aStitch.size();idart++){
    const unsigned int id_e0 = aStitch[idart].id_e1;
    const unsigned int id_e1 = aStitch[idart].id_e2;
    bool is_same_dir = aStitch[idart].is_same_dir;
    unsigned int id_v0s,id_v0e;  cad_2d.GetIdVertex_Edge(id_v0s,id_v0e,id_e0);
    unsigned int id_v1s,id_v1e;  cad_2d.GetIdVertex_Edge(id_v1s,id_v1e,id_e1);    
    const Com::CVector2D& v0s = cad_2d.GetVertexCoord(id_v0s);
    const Com::CVector2D& v1e = cad_2d.GetVertexCoord(id_v1e);      
    const Com::CVector2D& v1s = cad_2d.GetVertexCoord(id_v1s);
    const Com::CVector2D& v0e = cad_2d.GetVertexCoord(id_v0e);      
    if( is_same_dir ){
      ::glBegin(GL_LINES);
      ::glVertex3d(v0s.x,v0s.y,0.2);      
      ::glVertex3d(v1s.x,v1s.y,0.2);
      ::glVertex3d(v1e.x,v1e.y,0.2);   
      ::glVertex3d(v0e.x,v0e.y,0.2);       
      ::glEnd();      
    }
    else{
      ::glBegin(GL_LINES);
      ::glVertex3d(v0s.x,v0s.y,0.2);
      ::glVertex3d(v1e.x,v1e.y,0.2);       
      ::glVertex3d(v1s.x,v1s.y,0.2);
      ::glVertex3d(v0e.x,v0e.y,0.2);       
      ::glEnd();
    }    
  }  
  if( is_lighting ){ ::glEnable(GL_LIGHTING); }
  if( is_texture  ){ ::glEnable(GL_TEXTURE_2D); }  
  ::glDisable(GL_LINE_SMOOTH);
  ::glDisable(GL_BLEND);  
}
