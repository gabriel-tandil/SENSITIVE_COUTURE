/*
 *  cloth_handler.cpp
 *  sensitive couture
 *
 *  Created by Nobuyuki Umetani on 9/1/10.
 *  Copyright 2010 The University of Tokyo and Columbia University. All rights reserved.
 *
 */

#include <math.h>
#include <vector>
#include <queue>
#include <stack>

#include "delfem/field_world.h"
#include "delfem/field.h"
#include "delfem/vector3d.h"
#include "delfem/vector2d.h"

#include "contact_target.h"
#include "cloth_handler.h"

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

using namespace Fem::Field;

void CClothHandler::Clear()
{
  for(unsigned int ip=0;ip<apPiece.size();ip++){ delete apPiece[ip]; }
  apPiece.clear();
  this->id_ea_hilight = 0; 
  this->itype = 0;
  if( pXYZs_ != 0 ){ delete pXYZs_; pXYZs_ = 0; }
  if( aNorm_ != 0 ){ delete aNorm_; aNorm_ = 0; }
  if( aTri_  != 0 ){ delete aTri_;  aTri_  = 0; }
  nnode_ = 0;
  ntri_ = 0;
}

void CClothHandler::AddClothPiece(unsigned int id_l_new, unsigned int id_l_old)
{
  
  int icp1 = -1;
  for(unsigned int i=0;i<apPiece.size();i++){
    if( apPiece[i]->id_l == id_l_old ){
      icp1 = i;
      break;
    }
  }  
  if( icp1 == -1 ){ return; }
  {
    CClothPiece* cp0 = new CClothPiece;
    CClothPiece* cp1 = apPiece[icp1];
    cp0->id_l = id_l_new;
    cp0->id_ea = 0;
    cp0->cent[0] = cp1->cent[0];
    cp0->cent[1] = cp1->cent[1];
    cp0->p[0] = cp1->p[0];  cp0->p[1] = cp1->p[1];  cp0->p[2] = cp1->p[2];
    cp0->n[0] = cp1->n[0];  cp0->n[1] = cp1->n[1];  cp0->n[2] = cp1->n[2];
    cp0->h[0] = cp1->h[0];  cp0->h[1] = cp1->h[1];  cp0->h[2] = cp1->h[2]; 
    cp0->radius = cp1->radius;
    apPiece.push_back(cp0);    
  }  
}

void CClothHandler::AddClothPiece(unsigned int id_l, double cent_x, double cent_y)
{
  int icp0 = -1;
  for(unsigned int i=0;i<apPiece.size();i++){
    if( apPiece[i]->id_l == id_l ){
      icp0 = i;
      break;
    }
  }  
  if( icp0 == -1 ){    
    CClothPiece* cp = new CClothPiece;
    cp->id_l = id_l;
    cp->id_ea = 0;
    cp->cent[0] = cent_x;
    cp->cent[1] = cent_y;
    cp->p[0] = cp->cent[0];
    cp->p[1] = cp->cent[1];      
    cp->p[2] = 0;
    cp->n[0] = 0;   cp->n[1] = 0;   cp->n[2] = 1;
    cp->h[0] = 1;   cp->h[1] = 0;   cp->h[2] = 0; 
    cp->radius = -1;
    apPiece.push_back(cp);    
  }
  else{
    CClothPiece* cp = apPiece[icp0];
    cp->cent[0] = cent_x;
    cp->cent[1] = cent_y;
  }
}

void CClothHandler::Transform_Cloth_Pan(unsigned int id_l, double anc_x, double anc_y, double anc_z)
{
  int icp0 = -1;
  for(unsigned int i=0;i<apPiece.size();i++){
    if( apPiece[i]->id_l == id_l ){
      icp0 = i;
      break;
    }
  }    
  if( icp0 == -1 ){ return; }
  CClothPiece* cp = apPiece[icp0];
  cp->p[0] = anc_x;
  cp->p[1] = anc_y;
  cp->p[2] = anc_z;    
}

void CClothHandler::SetRadius(unsigned int id_l, double r)
{
  int icp0 = -1;
  for(unsigned int i=0;i<apPiece.size();i++){
    if( apPiece[i]->id_l == id_l ){
      icp0 = i;
      break;
    }
  }    
  if( icp0 == -1 ){ return; }
  CClothPiece* cp = apPiece[icp0];
  cp->radius = r;
}

static void MatVec3(const double m[9], const double x[3], double y[3]){
  y[0] = m[0]*x[0] + m[1]*x[1] + m[2]*x[2];
  y[1] = m[3]*x[0] + m[4]*x[1] + m[5]*x[2];
  y[2] = m[6]*x[0] + m[7]*x[1] + m[8]*x[2];
}

static void VecMat3(const double x[3], const double m[9],  double y[3]){
  y[0] = m[0]*x[0] + m[3]*x[1] + m[6]*x[2];
  y[1] = m[1]*x[0] + m[4]*x[1] + m[7]*x[2];
  y[2] = m[2]*x[0] + m[5]*x[1] + m[8]*x[2];
}

void CClothHandler::Transform_Cloth_RotBryantAngle(unsigned int id_l, double phi, double theta, double psi)
{  
  int icp0 = -1;
  for(unsigned int i=0;i<apPiece.size();i++){
    if( apPiece[i]->id_l == id_l ){
      icp0 = i;
      break;
    }
  }    
  if( icp0 == -1 ){ return; }
  ////
  phi   *= 3.1416/180.0;
  theta *= 3.1416/180.0;
  psi   *= 3.1416/180.0;  
  ////
  const double mat[9] = {
    cos(psi)*cos(theta),	cos(psi)*sin(theta)*sin(phi)-sin(psi)*cos(phi), cos(psi)*sin(theta)*cos(phi)+sin(psi)*sin(phi),
    sin(psi)*cos(theta),	sin(psi)*sin(theta)*sin(phi)+cos(psi)*cos(phi), sin(psi)*sin(theta)*cos(phi)-cos(psi)*sin(phi),
    -sin(theta),		    	cos(theta)*sin(phi),							              cos(theta)*cos(phi)};  
  ////
  CClothPiece* cp = apPiece[icp0];  
  double res[3];
  ////
  MatVec3(mat,cp->n,res);
  cp->n[0] = res[0];
  cp->n[1] = res[1];
  cp->n[2] = res[2];      
  ////
  MatVec3(mat,cp->h,res);
  cp->h[0] = res[0];
  cp->h[1] = res[1];
  cp->h[2] = res[2];        
}


void CClothHandler::BuildClothMeshTopology(unsigned int id_base, unsigned int id_field_disp, const CFieldWorld& world)
{
  const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);
  const CField& field = world.GetField(id_field_disp);
  const CNodeAry::CNodeSeg& ns_c = field.GetNodeSeg(CORNER,false,world);
  const std::vector<unsigned int>& aIdEA = field.GetAryIdEA();
  ////
  for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
    const unsigned int id_ea = aIdEA[iiea];
//    std::cout << "try ea" << id_ea << std::endl;
    unsigned int id_l;  
    Cad::CAD_ELEM_TYPE cad_elem_type;  
    conv.GetIdCad_fromIdEA(id_ea,id_l,cad_elem_type);
    if( cad_elem_type != Cad::LOOP ){ continue; }
    if( id_l == 0 ) continue;
//    std::cout << "buid topology : " << id_ea << " " << id_l << std::endl;
    const Fem::Field::CElemAry::CElemSeg& es = field.GetElemSeg(id_ea,CORNER,false,world);
    unsigned int ntri = es.Size();
    unsigned int* lnods = new unsigned int [ntri*3];
    for(unsigned int itri=0;itri<ntri;itri++){ es.GetNodes(itri,lnods+itri*3); }
    ////
    int piece_ary_ind = -1;
    for(unsigned int i=0;i<apPiece.size();i++){
      if( apPiece[i]->id_l == id_l ){
        piece_ary_ind = i;
        break;
      }
    }
    if( piece_ary_ind == -1 ){    
      CClothPiece* cp = new CClothPiece;
      cp->topo.SetTriAry(ntri,lnods,ns_c.Size());    
      double x_min,x_max, y_min,y_max;
      for(unsigned int ielem=0;ielem<es.Size();ielem++){
        unsigned int no[3]; es.GetNodes(ielem, no);
        for(unsigned int inoel=0;inoel<3;inoel++){
          unsigned int ino1 = no[inoel];
          double c[3];  ns_c.GetValue(ino1,c);
          if( ielem==0 && inoel==0 ){ 
            x_min = x_max = c[0];
            y_min = y_max = c[1];
          }
          else{
            x_min = ( x_min < c[0] ) ? x_min : c[0];
            x_max = ( x_max > c[0] ) ? x_max : c[0];            
            y_min = ( y_min < c[1] ) ? y_min : c[1];
            y_max = ( y_max > c[1] ) ? y_max : c[1];
          }
        }
      }
      cp->id_l = id_l;
      cp->id_ea = id_ea;
      cp->cent[0] = (x_min+x_max)*0.5;
      cp->cent[1] = (y_min+y_max)*0.5;    
      cp->p[0] = cp->cent[0];
      cp->p[1] = cp->cent[1];      
      cp->p[2] = 0;
      cp->n[0] = 0;   cp->n[1] = 0;   cp->n[2] = 1;
      cp->h[0] = 1;   cp->h[1] = 0;   cp->h[2] = 0;
      cp->radius = -1;
      apPiece.push_back(cp);
    }
    else{
      CClothPiece* cp = apPiece[piece_ary_ind];
      cp->id_ea = id_ea;
      cp->topo.SetTriAry(ntri,lnods,ns_c.Size());    
    }
    delete[] lnods;
  }
}


bool CClothHandler::SetClothLocation(unsigned int id_field_disp, CFieldWorld& world)
{
  CField& field = world.GetField(id_field_disp);
  const std::vector<unsigned int>& aIdEA = field.GetAryIdEA();
  for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
    unsigned int id_ea = aIdEA[iiea];
    int piece_ary_ind = -1;
    for(unsigned int i=0;i<apPiece.size();i++){
      if( apPiece[i]->id_ea == id_ea ){
        piece_ary_ind = i;
        break;
      }
    }
    if( piece_ary_ind == -1 ) continue;    
    {
      const Fem::Field::CNodeAry::CNodeSeg& ns_c = field.GetNodeSeg(CORNER,false,world,VALUE);
      Fem::Field::CNodeAry::CNodeSeg& ns_u = field.GetNodeSeg(CORNER,true, world,VALUE);
      const Fem::Field::CElemAry::CElemSeg& es = field.GetElemSeg(id_ea,CORNER,true,world);
      const double* cent = apPiece[piece_ary_ind]->cent;
      const double* h    = apPiece[piece_ary_ind]->h;
      const double* n    = apPiece[piece_ary_ind]->n;
      const double* p    = apPiece[piece_ary_ind]->p;      
      const double  r    = apPiece[piece_ary_ind]->radius;
      double v[3]; Com::Cross3D(v,n,h);
      assert( es.Length() == 3 );
      for(unsigned int ielem=0;ielem<es.Size();ielem++){
        unsigned int no[3];	es.GetNodes(ielem,no);
        for(unsigned int inoel=0;inoel<3;inoel++){
          double C[3];	ns_c.GetValue(no[inoel],C);
          const double X = C[0] - cent[0];
          const double Y = C[1] - cent[1];
          double u[3];
          if( r <= 0.0 ){
            u[0] = p[0] + X*h[0] + Y*v[0] - C[0];
            u[1] = p[1] + X*h[1] + Y*v[1] - C[1];
            u[2] = p[2] + X*h[2] + Y*v[2] - C[2];
          }
          else {            
            const double theta = Y/r;
            double Y2 = sin(theta)*r;
            double Z  = (1-cos(theta))*r;
            u[0] = p[0] + X*h[0] + Y2*v[0] - Z*n[0] - C[0];
            u[1] = p[1] + X*h[1] + Y2*v[1] - Z*n[1] - C[1];
            u[2] = p[2] + X*h[2] + Y2*v[2] - Z*n[2] - C[2];            
          }          
          ns_u.SetValue(no[inoel],0,u[0]);
          ns_u.SetValue(no[inoel],1,u[1]);
          ns_u.SetValue(no[inoel],2,u[2]);          
        }
      }          
    }
  }  
  return true;
}

bool CClothHandler::GetAnchor_2D(double r[2], unsigned int id_ea) const
{
  for(unsigned int i=0;i<apPiece.size();i++){
    if( apPiece[i]->id_ea == id_ea ){
      r[0] = apPiece[i]->cent[0];
      r[1] = apPiece[i]->cent[1];
      return true;
    }
  }    
  return false;
}

bool CClothHandler::GetAnchor_2D_Loop(double r[2], unsigned int id_l) const
{
  for(unsigned int i=0;i<apPiece.size();i++){
    if( apPiece[i]->id_l == id_l ){
      r[0] = apPiece[i]->cent[0];
      r[1] = apPiece[i]->cent[1];
      return true;
    }
  }    
  return false;
  
}

bool CClothHandler::MoveAnchor_2D(double dx, double dy, unsigned int id_ea)
{
  for(unsigned int i=0;i<apPiece.size();i++){
    if( apPiece[i]->id_ea == id_ea ){
      apPiece[i]->cent[0] += dx;
      apPiece[i]->cent[1] += dy;
      return true;
    }
  }    
  return false;
}

bool CClothHandler::GetAnchor_3D(double p[3], double n[3], double h[3], unsigned int id_ea) const
{
  for(unsigned int i=0;i<apPiece.size();i++){
    if( apPiece[i]->id_ea == id_ea ){
      p[0] = apPiece[i]->p[0];
      p[1] = apPiece[i]->p[1];
      p[2] = apPiece[i]->p[2];
      ////
      n[0] = apPiece[i]->n[0];
      n[1] = apPiece[i]->n[1];
      n[2] = apPiece[i]->n[2];
      ////
      h[0] = apPiece[i]->h[0];
      h[1] = apPiece[i]->h[1];
      h[2] = apPiece[i]->h[2];      
      return true;
    }
  }    
  return false;  
}

void CClothHandler::Draw(unsigned int imode) const
{
  if(imode == 1 ){
    bool is_texture  = ::glIsEnabled(GL_TEXTURE_2D);
    bool is_lighting = ::glIsEnabled(GL_LIGHTING);
    ::glDisable(GL_LIGHTING);
    ::glDisable(GL_TEXTURE_2D);
    ::glPointSize(10);
    ::glLineWidth(3);
    ////
    //  const Fem::Field::CField& field_disp = world.GetField(id_field_disp);
    //  const std::vector<unsigned int>& aIdEA = field_disp.GetAryIdEA();
    for(unsigned int ip=0;ip<apPiece.size();ip++){
      unsigned int id_ea = apPiece[ip]->id_ea;
      double* p = apPiece[ip]->p;
      double* n = apPiece[ip]->n;    
      double* h = apPiece[ip]->h;
      double g1 = 0.03;
      const double v1[3] = { p[0]+g1*n[0], p[1]+g1*n[1], p[2]+g1*n[2] };
      double g2 = 0.3;
      const double v2[3] = { p[0]+g2*n[0], p[1]+g2*n[1], p[2]+g2*n[2] };
      const double v3[3] = { p[0]+g2*h[0], p[1]+g2*h[1], p[2]+g2*h[2] };
      if( id_ea == id_ea_hilight && itype == 1 ){ ::glColor3d(1,1,0);  }
      else{ ::glColor3d(0,1,0); }
      ::glBegin(GL_POINTS);          
      ::glVertex3dv(v1);       
      ::glEnd();
      ::glBegin(GL_LINES);
      ::glColor3d(1,0,0);
      ::glVertex3dv(v1);
      ::glVertex3dv(v3);
      ::glEnd();
      if( id_ea != id_ea_hilight ) continue;
      ////
      ::glBegin(GL_LINES);
      ::glColor3d(0,0,1);          
      ::glVertex3dv(v1);
      ::glVertex3dv(v2);
      ::glEnd();
      {
        if( itype == 2 ){ ::glColor3d(1,1,0); }
        else{ ::glColor3d(0,1,0); }          
        const double l3 = 0.2;
        double v[3]; Com::Cross3D(v, n, h);
        const unsigned int ndiv = 20;
        ::glBegin(GL_LINE_LOOP);
        for(unsigned int idiv=0;idiv<ndiv;idiv++){              
          double theta = 3.14*2.0*idiv/ndiv;
          double p[3];
          p[0] = v1[0]+l3*cos(theta)*h[0]+l3*sin(theta)*v[0];
          p[1] = v1[1]+l3*cos(theta)*h[1]+l3*sin(theta)*v[1];
          p[2] = v1[2]+l3*cos(theta)*h[2]+l3*sin(theta)*v[2]; 
          ::glVertex3dv(p);
        }
        ::glEnd();
      }          
    }
    if( is_lighting ){ glEnable(GL_LIGHTING); }
    if( is_texture  ){ glEnable(GL_TEXTURE_2D); }
  }
  
  bool is_lighting = ::glIsEnabled(GL_LIGHTING);
  ::glEnable(GL_LIGHTING);
  bool is_texture  = ::glIsEnabled(GL_TEXTURE_2D);
  ::glDisable(GL_TEXTURE_2D);    
  {
    float gray[4] = {1,1,0.5,1};
//    float gray[4] = {0.3,0.3,0.3,1};    
    ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, gray);
    float shine[4] = {0,0,0,0};
    ::glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, shine);
    ::glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 127.0);
    //    ::glColor3d(1,1,1);
  }
  ::glBegin(GL_TRIANGLES);
  for(unsigned int itri=0;itri<ntri_;itri++){
    const unsigned int i1 = aTri_[itri*3+0];
    const unsigned int i2 = aTri_[itri*3+1];
    const unsigned int i3 = aTri_[itri*3+2];
    ::glNormal3dv(aNorm_+i1*3);   ::glVertex3dv(pXYZs_+i1*3);
    ::glNormal3dv(aNorm_+i2*3);   ::glVertex3dv(pXYZs_+i2*3);
    ::glNormal3dv(aNorm_+i3*3);   ::glVertex3dv(pXYZs_+i3*3);
  }
  ::glEnd();
  if( !is_lighting ){ ::glDisable(GL_LIGHTING); }
  if(  is_texture  ){ ::glEnable(GL_TEXTURE_2D); } 
}

void Translate3D(double out[3], 
                 const double in[3],
                 const double trans0[3], const double rot[3], const double trans1[3])
{  
  const double c2[3] = { in[0]+trans0[0], in[1]+trans0[1], in[2]+trans0[2] };
  const double c3[3] = {
    rot[0]*c2[0] + rot[1]*c2[1] + rot[2]*c2[2],
    rot[3]*c2[0] + rot[4]*c2[1] + rot[5]*c2[2],
    rot[6]*c2[0] + rot[7]*c2[1] + rot[8]*c2[2] };
  out[0] = c3[0] + trans1[0];
  out[1] = c3[1] + trans1[1];
  out[2] = c3[2] + trans1[2];  
}

bool CClothHandler::Pick
(double scrx, double scry,
 const double trans0[3], const double rot[3], const double trans1[3], 
 const double dir[3], const double org[3],
 unsigned int id_field_disp, const CFieldWorld& world)
{  
  this->id_ea_hilight = 0;
  const Fem::Field::CField& field_disp = world.GetField(id_field_disp);
  const Fem::Field::CNodeAry::CNodeSeg& ns_c = field_disp.GetNodeSeg(CORNER,false,world,VALUE);
  const Fem::Field::CNodeAry::CNodeSeg& ns_u = field_disp.GetNodeSeg(CORNER,true, world,VALUE);
  std::vector<unsigned int> aIdEA = field_disp.GetAryIdEA();
  for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
    const unsigned int id_ea = aIdEA[iiea];
    const Fem::Field::CElemAry::CElemSeg& es = field_disp.GetElemSeg(id_ea,CORNER,true,world);
    assert( es.Length() == 3 );
    for(unsigned int ielem=0;ielem<es.Size();ielem++){
      unsigned int no[3];	es.GetNodes(ielem,no);
      double C[3][3];	ns_c.GetValue(no[0],C[0]);	ns_c.GetValue(no[1],C[1]);	ns_c.GetValue(no[2],C[2]);
      double u[3][3];	ns_u.GetValue(no[0],u[0]);	ns_u.GetValue(no[1],u[1]);	ns_u.GetValue(no[2],u[2]);
      double c0[3][3] = {
        { C[0][0]+u[0][0], C[0][1]+u[0][1], C[0][2]+u[0][2] },
        { C[1][0]+u[1][0], C[1][1]+u[1][1], C[1][2]+u[1][2] },
        { C[2][0]+u[2][0], C[2][1]+u[2][1], C[2][2]+u[2][2] } };
      double area0 = Com::TriArea3D(c0[0],c0[1],c0[2]);
      double c1[3][3];
      for(unsigned int i=0;i<3;i++){
        Translate3D(c1[i],c0[i],trans0,rot,trans1);
/*
        const double c2[3] = { c0[i][0]+trans0[0], c0[i][1]+trans0[1], c0[i][2]+trans0[2] };
        const double c3[3] = {
          rot[0]*c2[0] + rot[1]*c2[1] + rot[2]*c2[2],
          rot[3]*c2[0] + rot[4]*c2[1] + rot[5]*c2[2],
          rot[6]*c2[0] + rot[7]*c2[1] + rot[8]*c2[2] };
        c1[i][0] = c3[0] + trans1[0];
        c1[i][1] = c3[1] + trans1[1];
        c1[i][2] = c3[2] + trans1[2];
*/ 
      }
      double area1 = Com::TriArea2D(c1[0],c1[1],c1[2]);
      if( area1 < area0*0.1 ) continue;
      const double scr[2] = {scrx,scry};
      double area2 = Com::TriArea2D(scr,  c1[1],c1[2]);
      double area3 = Com::TriArea2D(c1[0],scr,  c1[2]);
      double area4 = Com::TriArea2D(c1[0],c1[1],scr  );
      if( area2 < -area1*0.01 ) continue;
      if( area3 < -area1*0.01 ) continue;
      if( area4 < -area1*0.01 ) continue;
//      std::cout << "hit" << std::endl;
      id_ea_hilight = id_ea;
    }
    if( id_ea_hilight != 0 ) break;
  }  
  if( id_ea_hilight == 0 ) return false;
  /////
  this->itype = 0;
//  if( id_ea_hilight_pre == id_ea_hilight ){    
  {
    int ip = -1;
    for(unsigned int i=0;i<apPiece.size();i++){
      if( apPiece[i]->id_ea == id_ea_hilight ){
        ip = i;
        break;
      }
    }
    if( ip == -1 ) return false;
    double* p = apPiece[ip]->p;
    double* n = apPiece[ip]->n;    
    double* h = apPiece[ip]->h;
    double g1 = 0.03;
    const double v1[3] = { p[0]+g1*n[0], p[1]+g1*n[1], p[2]+g1*n[2] };
    double l3 = 0.1;
    double v[3]; Com::Cross3D(v, n, h);
    double c2[3][3];
    Translate3D(c2[0],v1, trans0, rot, trans1);
    const unsigned int ndiv = 20;    
    double p1[3], q1[3];
    for(unsigned int idiv=0;idiv<ndiv;idiv++){              
      double theta1 = 3.14*2.0*idiv/ndiv;
      p1[0] = v1[0]+l3*cos(theta1)*h[0]+l3*sin(theta1)*v[0];
      p1[1] = v1[1]+l3*cos(theta1)*h[1]+l3*sin(theta1)*v[1];
      p1[2] = v1[2]+l3*cos(theta1)*h[2]+l3*sin(theta1)*v[2]; 
      double theta2 = 3.14*2.0*(idiv+1)/ndiv;
      q1[0] = v1[0]+l3*cos(theta2)*h[0]+l3*sin(theta2)*v[0];
      q1[1] = v1[1]+l3*cos(theta2)*h[1]+l3*sin(theta2)*v[1];
      q1[2] = v1[2]+l3*cos(theta2)*h[2]+l3*sin(theta2)*v[2];       
      ////
      Translate3D(c2[1],p1, trans0, rot, trans1);      
      Translate3D(c2[2],q1, trans0, rot, trans1); 
      double area1 = Com::TriArea2D(c2[0],c2[1],c2[2]);
      const double scr[2] = {scrx,scry};
      double area2 = Com::TriArea2D(scr,  c2[1],c2[2]);
      double area3 = Com::TriArea2D(c2[0],scr,  c2[2]);
      double area4 = Com::TriArea2D(c2[0],c2[1],scr  );
      if( area2 < -area1*0.01 ) continue;
      if( area3 < -area1*0.01 ) continue;
      if( area4 < -area1*0.01 ) continue;      
      this->itype = 1; // pick origin
    }            
    if( itype != 1 ){ 
      double tmp1 = Com::Dot3D(dir, n);
      if( fabs(tmp1) > 1.0e-10 ){
        const double tmp2 = (p[0]-org[0])*n[0] + (p[1]-org[1])*n[1] + (p[2]-org[2])*n[2];
        const double tmp3 = tmp2/tmp1;
        hit_pos[0] = org[0]+tmp3*dir[0];
        hit_pos[1] = org[1]+tmp3*dir[1];
        hit_pos[2] = org[2]+tmp3*dir[2];
        itype = 2; // rotate        
      }
      else{
        itype = 0;
      }
    }
  }
  return true;
}

void CClothHandler::SetObjectMesh(const std::vector<unsigned int>& aTri,
                                  const std::vector<double>& aXYZ)
{  
  nnode_ = aXYZ.size()/3;
  if( pXYZs_ != 0 ){ delete[] pXYZs_; } 
	pXYZs_ = new double [nnode_*3];
	for(unsigned int ino=0;ino<nnode_;ino++){
		pXYZs_[ino*3+0] = aXYZ[ino*3+0];
		pXYZs_[ino*3+1] = aXYZ[ino*3+1];
		pXYZs_[ino*3+2] = aXYZ[ino*3+2];
	}		

  ntri_ = aTri.size()/3;
  if( aTri_ != 0 ){ delete[] aTri_; }  
	aTri_ = new unsigned int [ntri_*3];
	for(unsigned int itri=0;itri<ntri_;itri++){			
		aTri_[itri*3+0] = aTri[itri*3+0];
		aTri_[itri*3+1] = aTri[itri*3+1];
		aTri_[itri*3+2] = aTri[itri*3+2];
	}
  
  if( aNorm_ != 0 ){ delete[] aNorm_; }
  aNorm_ = new double [nnode_*3];
  for(unsigned int i=0;i<nnode_*3;i++){ aNorm_[i] = 0; }
  for(unsigned int itri=0;itri<ntri_;itri++){
    unsigned int i1 = aTri_[itri*3+0];
    unsigned int i2 = aTri_[itri*3+1];
    unsigned int i3 = aTri_[itri*3+2];
    double un[3], area;    
    Com::UnitNormalAreaTri3D(un,area, pXYZs_+i1*3, pXYZs_+i2*3, pXYZs_+i3*3);    
    aNorm_[i1*3+0] += un[0];  aNorm_[i1*3+1] += un[1];  aNorm_[i1*3+2] += un[2];
    aNorm_[i2*3+0] += un[0];  aNorm_[i2*3+1] += un[1];  aNorm_[i2*3+2] += un[2];    
    aNorm_[i3*3+0] += un[0];  aNorm_[i3*3+1] += un[1];  aNorm_[i3*3+2] += un[2];    
  }
  for(unsigned int ino=0;ino<nnode_;ino++){
    double invlen = 1.0/Com::Length3D(aNorm_+ino*3);
    aNorm_[ino*3+0] *= invlen;
    aNorm_[ino*3+1] *= invlen;
    aNorm_[ino*3+2] *= invlen;    
  }
}
