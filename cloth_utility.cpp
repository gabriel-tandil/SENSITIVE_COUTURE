/*
 *  cloth_utility.cpp
 *  sensitive couture
 *
 *  Created by Nobuyuki Umetani on 9/17/10.
 *  Copyright 2010 The University of Tokyo and Columbia University. All rights reserved.
 *
 */


#if defined(__VISUALC__)
#pragma warning ( disable : 4996 )
#pragma warning ( disable : 4786 )
#endif

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <memory>

#include "delfem/ls/solver_ls_iter.h"

#include "delfem/field_world.h"
#include "delfem/field.h"
#include "delfem/drawer_field.h"
#include "delfem/drawer_field_face.h"
#include "delfem/drawer_field_edge.h"
#include "delfem/ls/solver_ls_iter.h"
#include "delfem/spatial_hash_grid2d.h"
#include "delfem/vector2d.h"
#include "delfem/vector3d.h"

#include "eqn_quad_bend.h"
#include "eqn_glue.h"
#include "eqn_contact3d.h"
#include "cloth_utility.h"

using namespace Fem::Field;

bool KineticDamping(double ke0, double ke1, double ke2,
                    unsigned int id_field_disp, Fem::Field::CFieldWorld& world)
{
  double c = (ke0-2*ke1+ke2)*0.5;
  //  std::cout << c << std::endl;
  if( fabs(c) < 1.0e-20 ) return false;
  double a = (ke0-ke2)*0.25/c;
  //  std::cout << c << " " << a << std::endl;
  if( a < -1.0 || a > 1.0 ) return false;
  std::cout << " ######## EXTREMUM ####### " << a << " " << c << std::endl;  
  if( c > 0 ) return false;
  //  return false;
  Fem::Field::CField& field = world.GetField(id_field_disp);  
	Fem::Field::CNodeAry::CNodeSeg& ns_v = field.GetNodeSeg(CORNER,true, world,VELOCITY); 
  ns_v.SetZero();
  return true;
}

void CopyValueVelo(unsigned int id0, unsigned int id1, Fem::Field::CFieldWorld& world)
{
  Fem::Field::CField& field0 = world.GetField(id0);  
	Fem::Field::CNodeAry::CNodeSeg& ns_u0 = field0.GetNodeSeg(CORNER,true, world,VALUE); 
	Fem::Field::CNodeAry::CNodeSeg& ns_v0 = field0.GetNodeSeg(CORNER,true, world,VELOCITY); 
  
  const Fem::Field::CField& field1 = world.GetField(id1);  
	const Fem::Field::CNodeAry::CNodeSeg& ns_u1 = field1.GetNodeSeg(CORNER,true, world,VALUE); 
	const Fem::Field::CNodeAry::CNodeSeg& ns_v1 = field1.GetNodeSeg(CORNER,true, world,VELOCITY); 
  
  assert( ns_u1.Size() == ns_u0.Size() );
  
  for(unsigned int ino=0;ino<ns_u0.Size();ino++){
    double u[3]; ns_u1.GetValue(ino,u);
    double v[3]; ns_v1.GetValue(ino,v);    
    ns_u0.SetValue(ino, 0, u[0]);
    ns_u0.SetValue(ino, 1, u[1]);
    ns_u0.SetValue(ino, 2, u[2]);
    ns_v0.SetValue(ino, 0, v[0]);
    ns_v0.SetValue(ino, 1, v[1]);
    ns_v0.SetValue(ino, 2, v[2]);        
  }  
}


void SetDeformedValue(unsigned int id0, unsigned int id1, Fem::Field::CFieldWorld& world)
{
  assert( world.IsIdField(id0) );
  Fem::Field::CField& field0 = world.GetField(id0);  
	Fem::Field::CNodeAry::CNodeSeg& ns_u0 = field0.GetNodeSeg(CORNER,true, world,VALUE); 
  
  const Fem::Field::CField& field1 = world.GetField(id1);  
	const Fem::Field::CNodeAry::CNodeSeg& ns_u1 = field1.GetNodeSeg(CORNER,true, world,VALUE); 
	const Fem::Field::CNodeAry::CNodeSeg& ns_X1 = field1.GetNodeSeg(CORNER,false,world,VALUE);   
  
  assert( ns_u1.Size() == ns_u0.Size() );
  
  for(unsigned int ino=0;ino<ns_u0.Size();ino++){
    double u[3]; ns_u1.GetValue(ino,u);
    double C[3]; ns_X1.GetValue(ino,C);
    ns_u0.SetValue(ino, 0, C[0]+u[0]);
    ns_u0.SetValue(ino, 1, C[1]+u[1]);
    ns_u0.SetValue(ino, 2, C[2]+u[2]);  
  }  
}

// Interpolate Value from field1 to field0
void InterpField(unsigned int id_field_disp0, unsigned int id_base0,
                 unsigned int id_field_disp1, unsigned int id_base1,
                 Fem::Field::CFieldWorld& world)
{
  // fine
  Fem::Field::CField& field0 = world.GetField(id_field_disp0);  
	Fem::Field::CNodeAry::CNodeSeg& ns_c0 = field0.GetNodeSeg(CORNER,false,world,VALUE);   
  Fem::Field::CNodeAry::CNodeSeg& ns_u0 = field0.GetNodeSeg(CORNER,true, world,VALUE); 
  
  // coarse
  const Fem::Field::CField& field1 = world.GetField(id_field_disp1);  
	const Fem::Field::CNodeAry::CNodeSeg& ns_c1 = field1.GetNodeSeg(CORNER,false,world,VALUE);   
  const Fem::Field::CNodeAry::CNodeSeg& ns_u1 = field1.GetNodeSeg(CORNER,true, world,VALUE); 
  
  std::vector<unsigned int> aFlg;
  aFlg.resize(ns_c0.Size(),0);
  const std::vector<unsigned int> aIdEA0 = field0.GetAryIdEA();
  const std::vector<unsigned int> aIdEA1 = field1.GetAryIdEA();  
  const Fem::Field::CIDConvEAMshCad& conv0 = world.GetIDConverter(id_base0);
  const Fem::Field::CIDConvEAMshCad& conv1 = world.GetIDConverter(id_base1);
  for(unsigned int iiea0=0;iiea0<aIdEA0.size();iiea0++){
    const unsigned int id_ea0 = aIdEA0[iiea0];
    unsigned int id_ea1 = 0;
    {
      unsigned int id_cad0;  Cad::CAD_ELEM_TYPE type_cad0; conv0.GetIdCad_fromIdEA(id_ea0, id_cad0,type_cad0);
      if( type_cad0 != Cad::LOOP ) continue;
      for(unsigned int iiea1=0;iiea1<aIdEA1.size();iiea1++){
        unsigned int iea1 = aIdEA1[iiea1];
        assert( world.IsIdEA(iea1) );
        unsigned int id_cad1;  Cad::CAD_ELEM_TYPE type_cad1; conv1.GetIdCad_fromIdEA(iea1, id_cad1,type_cad1); 
        if( type_cad1 != Cad::LOOP ) continue;
        if( id_cad1 == id_cad0 ){
          id_ea1 = iea1;
          break;
        }
      }
      assert( world.IsIdEA(id_ea1) );
    }
    if( !field1.GetIdElemSeg(id_ea1,CORNER,false,world) ) return;
    const Fem::Field::CElemAry::CElemSeg& es1 = field1.GetElemSeg(id_ea1,CORNER,false,world);    
    assert( es1.Length() == 3 );        
    CSpatialHash_Grid2D shg;
    {
      double min_x,max_x, min_y,max_y;
      for(unsigned int ielem1=0;ielem1<es1.Size();ielem1++){
        unsigned int no1[3];	es1.GetNodes(ielem1,no1);
        for(unsigned int inoes=0;inoes<3;inoes++){
          double C1[2];	ns_c1.GetValue(no1[inoes],C1);
          if( ielem1==0 && inoes==0 ){
            min_x = max_x = C1[0];
            min_y = max_y = C1[1];            
          }
          else{
            min_x = ( C1[0] < min_x ) ? C1[0] : min_x;
            max_x = ( C1[0] > max_x ) ? C1[0] : max_x;            
            min_y = ( C1[1] < min_y ) ? C1[1] : min_y;
            max_y = ( C1[1] > max_y ) ? C1[1] : max_y;
          }
        }
      }
      const double cent[2] = {(min_x+max_x)*0.5,(min_y+max_y)*0.5};
      double hw =  ( (max_x-min_x) > (max_y-min_y) ) ? (max_x-min_x)*0.5 : (max_y-min_y)*0.5;
      shg = CSpatialHash_Grid2D(64,cent,hw*1.154432165);
      for(unsigned int ielem1=0;ielem1<es1.Size();ielem1++){
        unsigned int no1[3];	es1.GetNodes(ielem1,no1);
        double C1[3][3];	ns_c1.GetValue(no1[0],C1[0]);	ns_c1.GetValue(no1[1],C1[1]); ns_c1.GetValue(no1[2],C1[2]);
        shg.AddTri(ielem1, C1[0], C1[1], C1[2]);
      }        
    }    
    const Fem::Field::CElemAry::CElemSeg& es0 = field0.GetElemSeg(id_ea0,CORNER,false,world);
    assert( es0.Length() == 3 );        
    std::vector<unsigned int> aIndTriCand;
    for(unsigned int ielem0=0;ielem0<es0.Size();ielem0++){
      unsigned int no0[3];	es0.GetNodes(ielem0,no0);      
      for(unsigned int inoes0=0;inoes0<3;inoes0++){
        unsigned int ino0 = no0[inoes0];        
        if( aFlg[ino0] == 1 ) continue;
        double C0[3]; ns_c0.GetValue(ino0,C0); 
        /////
//        double minsqlen = -1;
//        unsigned int minsqlen_ino = 0;
        shg.Find_IncludeTriCand(C0,aIndTriCand);
        /*
        for(unsigned int icand=0;icand<aIndTriCand.size();icand++){
          const unsigned int ielem1 = aIndTriCand[icand];
          unsigned int no1[3];	es1.GetNodes(ielem1,no1);
          double C1[3][3];	ns_c1.GetValue(no1[0],C1[0]);	ns_c1.GetValue(no1[1],C1[1]); ns_c1.GetValue(no1[2],C1[2]);
          if( minsqlen < 0 || minsqlen > Com::SqDistance2D(C0, C1[0]) ){ minsqlen = Com::SqDistance2D(C0,C1[0]); minsqlen_ino = no1[0]; }
          if( minsqlen < 0 || minsqlen > Com::SqDistance2D(C0, C1[1]) ){ minsqlen = Com::SqDistance2D(C0,C1[1]); minsqlen_ino = no1[1]; }
          if( minsqlen < 0 || minsqlen > Com::SqDistance2D(C0, C1[2]) ){ minsqlen = Com::SqDistance2D(C0,C1[2]); minsqlen_ino = no1[2]; }
          double area4 = Com::TriArea2D(C1[0],C1[1],C1[2]);
          double area0 = Com::TriArea2D(C0,   C1[1],C1[2]);
          double area1 = Com::TriArea2D(C1[0],C0,   C1[2]);
          double area2 = Com::TriArea2D(C1[0],C1[1],C0   );
          if( area0 < -area4*0.01 ) continue;
          if( area1 < -area4*0.01 ) continue;
          if( area2 < -area4*0.01 ) continue;
          const double r0 = area0/area4, r1 = area1/area4, r2 = area2/area4;
          aFlg[ino0] = 1;
          double u1[3][3];	ns_u1.GetValue(no1[0],u1[0]);	ns_u1.GetValue(no1[1],u1[1]); ns_u1.GetValue(no1[2],u1[2]);  
          const double c1[3][3] = {
            { u1[0][0]+C1[0][0], u1[0][1]+C1[0][1], u1[0][2]+C1[0][2] },
            { u1[1][0]+C1[1][0], u1[1][1]+C1[1][1], u1[1][2]+C1[1][2] },
            { u1[2][0]+C1[2][0], u1[2][1]+C1[2][1], u1[2][2]+C1[2][2] } };  
          double c0[3] = {
            c1[0][0]*r0 + c1[1][0]*r1 + c1[2][0]*r2,
            c1[0][1]*r0 + c1[1][1]*r1 + c1[2][1]*r2,
            c1[0][2]*r0 + c1[1][2]*r1 + c1[2][2]*r2 };          
          double u0[3] = { c0[0]-C0[0], c0[1]-C0[1], c0[2]-C0[2] };    
          ns_u0.SetValue(ino0,0,u0[0]);
          ns_u0.SetValue(ino0,1,u0[1]);
          ns_u0.SetValue(ino0,2,u0[2]);                
          break;
        }                     
        if( aFlg[ino0] == 0 ){
          aFlg[ino0] = 1;
          double C1[3];	ns_c1.GetValue(minsqlen_ino,C1);
          double u1[3]; ns_u1.GetValue(minsqlen_ino,u1);
          double c1[3] = { u1[0]+C1[0], u1[1]+C1[1], u1[2]+C1[2] };
          double u0[3] = { c1[0]-C0[0], c1[1]-C0[1], c1[2]-C0[2] };              
          ns_u0.SetValue(ino0,0,u0[0]);
          ns_u0.SetValue(ino0,1,u0[1]);
          ns_u0.SetValue(ino0,2,u0[2]);                          
        }*/
        unsigned int ielem_in = -1;
        for(unsigned int icand=0;icand<aIndTriCand.size();icand++){
          const unsigned int ielem1 = aIndTriCand[icand];
          unsigned int no1[3];	es1.GetNodes(ielem1,no1);
          double C1[3][3];	ns_c1.GetValue(no1[0],C1[0]);	ns_c1.GetValue(no1[1],C1[1]); ns_c1.GetValue(no1[2],C1[2]);
          double area4 = Com::TriArea2D(C1[0],C1[1],C1[2]);
          double area0 = Com::TriArea2D(C0,   C1[1],C1[2]);
          double area1 = Com::TriArea2D(C1[0],C0,   C1[2]);
          double area2 = Com::TriArea2D(C1[0],C1[1],C0   );
          if( area0 < -area4*0.01 ) continue;
          if( area1 < -area4*0.01 ) continue;
          if( area2 < -area4*0.01 ) continue;
          ielem_in = ielem1;
        }
        if( ielem_in == -1 ){
          unsigned int ielem_min = -1;
          double dist_min = -1;
          for(unsigned int ielem1=0;ielem1<es1.Size();ielem1++){
            unsigned int no1[3];	es1.GetNodes(ielem1,no1);
            double C1[3][3];	ns_c1.GetValue(no1[0],C1[0]);	ns_c1.GetValue(no1[1],C1[1]); ns_c1.GetValue(no1[2],C1[2]);
            const double GC1[3] = {
              (C1[0][0]+C1[1][0]+C1[2][0])/3.0,
              (C1[0][1]+C1[1][1]+C1[2][1])/3.0,
              (C1[0][2]+C1[1][2]+C1[2][2])/3.0 };
            double dist = Com::Distance3D(C0,GC1);
            if( ielem_min == -1 || dist < dist_min ){ dist_min = dist; ielem_min = ielem1; }
          }  
          ielem_in = ielem_min;
        }
        {
          unsigned int no1[3];	es1.GetNodes(ielem_in,no1);
          double C1[3][3];	ns_c1.GetValue(no1[0],C1[0]);	ns_c1.GetValue(no1[1],C1[1]); ns_c1.GetValue(no1[2],C1[2]);
          double area4 = Com::TriArea2D(C1[0],C1[1],C1[2]);
          double area0 = Com::TriArea2D(C0,   C1[1],C1[2]);
          double area1 = Com::TriArea2D(C1[0],C0,   C1[2]);
          double area2 = Com::TriArea2D(C1[0],C1[1],C0   );          
          const double r0 = area0/area4, r1 = area1/area4, r2 = area2/area4;
          aFlg[ino0] = 1;
          double u1[3][3];	ns_u1.GetValue(no1[0],u1[0]);	ns_u1.GetValue(no1[1],u1[1]); ns_u1.GetValue(no1[2],u1[2]);  
          const double c1[3][3] = {
            { u1[0][0]+C1[0][0], u1[0][1]+C1[0][1], u1[0][2]+C1[0][2] },
            { u1[1][0]+C1[1][0], u1[1][1]+C1[1][1], u1[1][2]+C1[1][2] },
            { u1[2][0]+C1[2][0], u1[2][1]+C1[2][1], u1[2][2]+C1[2][2] } };  
          double c0[3] = {
            c1[0][0]*r0 + c1[1][0]*r1 + c1[2][0]*r2,
            c1[0][1]*r0 + c1[1][1]*r1 + c1[2][1]*r2,
            c1[0][2]*r0 + c1[1][2]*r1 + c1[2][2]*r2 };          
          double u0[3] = { c0[0]-C0[0], c0[1]-C0[1], c0[2]-C0[2] };    
          ns_u0.SetValue(ino0,0,u0[0]);
          ns_u0.SetValue(ino0,1,u0[1]);
          ns_u0.SetValue(ino0,2,u0[2]);                
//          break;
        }                             
      }
    }
  }  
}


// Build Interp Relation
void FindBaseInterp(unsigned int id_base0,  // fine
                    unsigned int id_base1,  // coarse
                    Fem::Field::CFieldWorld& world,
                    std::vector<CInterpBarycentric>& aInterp)
{
  Fem::Field::CField& field0 = world.GetField(id_base0);  
	Fem::Field::CNodeAry::CNodeSeg& ns_c0 = field0.GetNodeSeg(CORNER,false,world,VALUE);   
  
  const Fem::Field::CField& field1 = world.GetField(id_base1);  
	const Fem::Field::CNodeAry::CNodeSeg& ns_c1 = field1.GetNodeSeg(CORNER,false,world,VALUE);   
  
  aInterp.resize(ns_c0.Size());
  
  std::vector<unsigned int> aFlg;
  aFlg.resize(ns_c0.Size(),0);
  const std::vector<unsigned int> aIdEA0 = field0.GetAryIdEA();
  const std::vector<unsigned int> aIdEA1 = field1.GetAryIdEA();  
  const Fem::Field::CIDConvEAMshCad& conv0 = world.GetIDConverter(id_base0);
  const Fem::Field::CIDConvEAMshCad& conv1 = world.GetIDConverter(id_base1);
  for(unsigned int iiea0=0;iiea0<aIdEA0.size();iiea0++){
    const unsigned int id_ea0 = aIdEA0[iiea0];
    unsigned int id_ea1 = 0;
    {
      unsigned int id_cad0;  Cad::CAD_ELEM_TYPE type_cad0; conv0.GetIdCad_fromIdEA(id_ea0, id_cad0,type_cad0);
      if( type_cad0 != Cad::LOOP ) continue;
      for(unsigned int iiea1=0;iiea1<aIdEA1.size();iiea1++){
        unsigned int iea1 = aIdEA1[iiea1];
        assert( world.IsIdEA(iea1) );
        unsigned int id_cad1;  Cad::CAD_ELEM_TYPE type_cad1; conv1.GetIdCad_fromIdEA(iea1, id_cad1,type_cad1); 
        if( type_cad1 != Cad::LOOP ) continue;
        if( id_cad1 == id_cad0 ){
          id_ea1 = iea1;
          break;
        }
      }
      assert( world.IsIdEA(id_ea1) );
    }
    const Fem::Field::CElemAry::CElemSeg& es0 = field0.GetElemSeg(id_ea0,CORNER,false,world);
    const Fem::Field::CElemAry::CElemSeg& es1 = field1.GetElemSeg(id_ea1,CORNER,false,world);    
    assert( es0.Length() == 3 );
    assert( es1.Length() == 3 );
    CSpatialHash_Grid2D shg;
    {
      double min_x,max_x, min_y,max_y;
      for(unsigned int ielem1=0;ielem1<es1.Size();ielem1++){
        unsigned int no1[3];	es1.GetNodes(ielem1,no1);
        for(unsigned int inoes=0;inoes<3;inoes++){
          double C1[2];	ns_c1.GetValue(no1[inoes],C1);
          if( ielem1==0 && inoes==0 ){
            min_x = max_x = C1[0];
            min_y = max_y = C1[1];            
          }
          else{
            min_x = ( C1[0] < min_x ) ? C1[0] : min_x;
            max_x = ( C1[0] > max_x ) ? C1[0] : max_x;            
            min_y = ( C1[1] < min_y ) ? C1[1] : min_y;
            max_y = ( C1[1] > max_y ) ? C1[1] : max_y;
          }
        }
      }
      const double cent[2] = {(min_x+max_x)*0.5,(min_y+max_y)*0.5};
      double hw =  ( (max_x-min_x) > (max_y-min_y) ) ? (max_x-min_x)*0.5 : (max_y-min_y)*0.5;
      shg = CSpatialHash_Grid2D(64,cent,hw*1.154432165);
      for(unsigned int ielem1=0;ielem1<es1.Size();ielem1++){
        unsigned int no1[3];	es1.GetNodes(ielem1,no1);
        double C1[3][3];	ns_c1.GetValue(no1[0],C1[0]);	ns_c1.GetValue(no1[1],C1[1]); ns_c1.GetValue(no1[2],C1[2]);
        shg.AddTri(ielem1, C1[0], C1[1], C1[2]);
      }        
    }
    std::vector<unsigned int> aIndTriCand;
    for(unsigned int ielem0=0;ielem0<es0.Size();ielem0++){
      unsigned int no0[3];	es0.GetNodes(ielem0,no0);      
      for(unsigned int inoes0=0;inoes0<3;inoes0++){
        unsigned int ino0 = no0[inoes0];        
        if( aFlg[ino0] == 1 ) continue;
        double C0[3]; ns_c0.GetValue(ino0,C0); 
        shg.Find_IncludeTriCand(C0,aIndTriCand);
        for(unsigned int icand=0;icand<aIndTriCand.size();icand++){
          const unsigned int ielem1 = aIndTriCand[icand];
          unsigned int no1[3];	es1.GetNodes(ielem1,no1);
          double C1[3][3];	ns_c1.GetValue(no1[0],C1[0]);	ns_c1.GetValue(no1[1],C1[1]); ns_c1.GetValue(no1[2],C1[2]);
          double area0 = Com::TriArea2D(C1[0],C1[1],C1[2]);
          double area1 = Com::TriArea2D(C0,   C1[1],C1[2]);
          double area2 = Com::TriArea2D(C1[0],C0,   C1[2]);
          double area3 = Com::TriArea2D(C1[0],C1[1],C0   );
          if( area1 < -area0*0.01 ) continue;
          if( area2 < -area0*0.01 ) continue;
          if( area3 < -area0*0.01 ) continue;
          const double r1 = area1/area0, r2 = area2/area0;
          aFlg[ino0] = 1;
          aInterp[ino0].ino1 = no1[0];
          aInterp[ino0].ino2 = no1[1];
          aInterp[ino0].ino3 = no1[2];
          aInterp[ino0].r1_ini = r1;
          aInterp[ino0].r2_ini = r2;          
          aInterp[ino0].r1 = r1;          
          aInterp[ino0].r2 = r2;
          aInterp[ino0].height = 0;
          break;
        }                     
        if( aFlg[ino0] != 0 ) continue;
        // can't find interpolation above
//        std::cout << "Cannot find interpolation" << " " << aIndTriCand.size() << std::endl;
        double height_min = -1;
        int ielem_min = -1;
        for(unsigned int icand=0;icand<aIndTriCand.size();icand++){
          const unsigned int ielem1 = aIndTriCand[icand];
          unsigned int no1[3];	es1.GetNodes(ielem1,no1);
          double C1[3][3];	ns_c1.GetValue(no1[0],C1[0]);	ns_c1.GetValue(no1[1],C1[1]); ns_c1.GetValue(no1[2],C1[2]);
          double h0 = -Com::TriArea2D(C0,   C1[1],C1[2]) / sqrt( Com::SqDistance2D(C1[1],C1[2]) );
          double h1 = -Com::TriArea2D(C1[0],C0,   C1[2]) / sqrt( Com::SqDistance2D(C1[0],C1[2]) );
          double h2 = -Com::TriArea2D(C1[0],C1[1],C0   ) / sqrt( Com::SqDistance2D(C1[0],C1[1]) );
          if( h0>0 && (ielem_min < 0 || h0<height_min) ){ ielem_min = ielem1; height_min = h0; }
          if( h1>0 && (ielem_min < 0 || h1<height_min) ){ ielem_min = ielem1; height_min = h1; }
          if( h2>0 && (ielem_min < 0 || h2<height_min) ){ ielem_min = ielem1; height_min = h2; }
        }
        if( ielem_min == -1 ){
          for(unsigned int ielem1=0;ielem1<es1.Size();ielem1++){
            unsigned int no1[3];	es1.GetNodes(ielem1,no1);
            double C1[3][3];	ns_c1.GetValue(no1[0],C1[0]);	ns_c1.GetValue(no1[1],C1[1]); ns_c1.GetValue(no1[2],C1[2]);
            double h0 = -Com::TriArea2D(C0,   C1[1],C1[2]) / sqrt( Com::SqDistance2D(C1[1],C1[2]) );
            double h1 = -Com::TriArea2D(C1[0],C0,   C1[2]) / sqrt( Com::SqDistance2D(C1[0],C1[2]) );
            double h2 = -Com::TriArea2D(C1[0],C1[1],C0   ) / sqrt( Com::SqDistance2D(C1[0],C1[1]) );
            if( h0>0 && (ielem_min < 0 || h0<height_min) ){ ielem_min = ielem1; height_min = h0; }
            if( h1>0 && (ielem_min < 0 || h1<height_min) ){ ielem_min = ielem1; height_min = h1; }
            if( h2>0 && (ielem_min < 0 || h2<height_min) ){ ielem_min = ielem1; height_min = h2; }
          }                  
        }
        assert( ielem_min >= 0 && ielem_min < es1.Size() );
//        std::cout << ielem_min << " " << height_min << std::endl;
        {
          unsigned int no1[3];	es1.GetNodes(ielem_min,no1);  // this sometime stack
          double C1[3][3];	ns_c1.GetValue(no1[0],C1[0]);	ns_c1.GetValue(no1[1],C1[1]); ns_c1.GetValue(no1[2],C1[2]);
          double area0 = Com::TriArea2D(C1[0],C1[1],C1[2]);
          double area1 = Com::TriArea2D(C0,   C1[1],C1[2]);
          double area2 = Com::TriArea2D(C1[0],C0,   C1[2]);
          const double r1 = area1/area0, r2 = area2/area0;
//          std::cout << r1 << " " << r2 << " " << r3 << std::endl;
          aFlg[ino0] = 1;
          aInterp[ino0].ino1 = no1[0];
          aInterp[ino0].ino2 = no1[1];
          aInterp[ino0].ino3 = no1[2];
          aInterp[ino0].r1_ini = r1;
          aInterp[ino0].r2_ini = r2;          
          aInterp[ino0].r1 = r1;          
          aInterp[ino0].r2 = r2;
          aInterp[ino0].height = 0;
          break;
        }                             
      }    
    }
  }
}


void MoveFineBaseCoord(unsigned int id_base0,
                       unsigned int id_base1,
                       Fem::Field::CFieldWorld& world,
                       const std::vector<CInterpBarycentric>& aInterp)
{
  // fine
  Fem::Field::CField& field0 = world.GetField(id_base0);  
	Fem::Field::CNodeAry::CNodeSeg& ns_c0 = field0.GetNodeSeg(CORNER,false,world,VALUE);   
  
  // coarse
  const Fem::Field::CField& field1 = world.GetField(id_base1);  
	const Fem::Field::CNodeAry::CNodeSeg& ns_c1 = field1.GetNodeSeg(CORNER,false,world,VALUE);   
  
  for(unsigned int ino0=0;ino0<ns_c0.Size();ino0++)
  {
    unsigned int ino1_0 = aInterp[ino0].ino1;
    unsigned int ino1_1 = aInterp[ino0].ino2;
    unsigned int ino1_2 = aInterp[ino0].ino3; 
    const double r0 = aInterp[ino0].r1_ini;
    const double r1 = aInterp[ino0].r2_ini;
    const double r2 = 1-r0-r1;
    double C1[3][3];	ns_c1.GetValue(ino1_0,C1[0]);	ns_c1.GetValue(ino1_1,C1[1]); ns_c1.GetValue(ino1_2,C1[2]);    
    const double C0[3] = {
      C1[0][0]*r0 + C1[1][0]*r1 + C1[2][0]*r2,
      C1[0][1]*r0 + C1[1][1]*r1 + C1[2][1]*r2,
      C1[0][2]*r0 + C1[1][2]*r1 + C1[2][2]*r2 };
    ns_c0.SetValue(ino0,0,C0[0]);
    ns_c0.SetValue(ino0,1,C0[1]);
    ns_c0.SetValue(ino0,2,C0[2]);        
  }  
}

void MoveFineDeformedCoord(unsigned int id_field_disp0, unsigned int id_base0,
                           unsigned int id_field_disp1, unsigned int id_base1,
                           Fem::Field::CFieldWorld& world,
                           const std::vector<CInterpBarycentric>& aInterp)
{
  // fine
  Fem::Field::CField& field0 = world.GetField(id_field_disp0);  
	Fem::Field::CNodeAry::CNodeSeg& ns_c0 = field0.GetNodeSeg(CORNER,false,world,VALUE);   
  Fem::Field::CNodeAry::CNodeSeg& ns_u0 = field0.GetNodeSeg(CORNER,true, world,VALUE); 
  
  // coarse
  const Fem::Field::CField& field1 = world.GetField(id_field_disp1);  
	const Fem::Field::CNodeAry::CNodeSeg& ns_c1 = field1.GetNodeSeg(CORNER,false,world,VALUE);   
  const Fem::Field::CNodeAry::CNodeSeg& ns_u1 = field1.GetNodeSeg(CORNER,true, world,VALUE); 
  
  for(unsigned int ino0=0;ino0<ns_c0.Size();ino0++)
  {
    const unsigned int ino1_0 = aInterp[ino0].ino1;
    const unsigned int ino1_1 = aInterp[ino0].ino2;
    const unsigned int ino1_2 = aInterp[ino0].ino3;
    const double r0 = aInterp[ino0].r1;
    const double r1 = aInterp[ino0].r2;    
    const double r2 = 1-r0-r1;
    double u1[3][3];	ns_u1.GetValue(ino1_0,u1[0]);	ns_u1.GetValue(ino1_1,u1[1]); ns_u1.GetValue(ino1_2,u1[2]);  
    double C1[3][3];	ns_c1.GetValue(ino1_0,C1[0]);	ns_c1.GetValue(ino1_1,C1[1]); ns_c1.GetValue(ino1_2,C1[2]);    
    const double c1[3][3] = {
      { u1[0][0]+C1[0][0], u1[0][1]+C1[0][1], u1[0][2]+C1[0][2] },
      { u1[1][0]+C1[1][0], u1[1][1]+C1[1][1], u1[1][2]+C1[1][2] },
      { u1[2][0]+C1[2][0], u1[2][1]+C1[2][1], u1[2][2]+C1[2][2] } };
    double c0[3] = {
      c1[0][0]*r0 + c1[1][0]*r1 + c1[2][0]*r2,
      c1[0][1]*r0 + c1[1][1]*r1 + c1[2][1]*r2,
      c1[0][2]*r0 + c1[1][2]*r1 + c1[2][2]*r2 };
    double height = aInterp[ino0].height;
    double n1[3],area1;  Com::UnitNormalAreaTri3D(n1,area1,c1[0],c1[1],c1[2]);
    c0[0] += height*n1[0];
    c0[1] += height*n1[1];
    c0[2] += height*n1[2];    
    double C0[3]; ns_c0.GetValue(ino0,C0);    
    double u0[3] = { c0[0]-C0[0], c0[1]-C0[1], c0[2]-C0[2] };    
    ns_u0.SetValue(ino0,0,u0[0]);
    ns_u0.SetValue(ino0,1,u0[1]);
    ns_u0.SetValue(ino0,2,u0[2]);      
  }
}

void UpdateFineDeformedInterp(unsigned int id_field_disp0, unsigned int id_base0,
                              unsigned int id_field_disp1, unsigned int id_base1,
                              const Fem::Field::CFieldWorld& world,
                              std::vector<CInterpBarycentric>& aInterp)
{
  // fine
  const Fem::Field::CField& field0 = world.GetField(id_field_disp0);  
	const Fem::Field::CNodeAry::CNodeSeg& ns_c0 = field0.GetNodeSeg(CORNER,false,world,VALUE);   
  const Fem::Field::CNodeAry::CNodeSeg& ns_u0 = field0.GetNodeSeg(CORNER,true, world,VALUE);
  
  // coarse
  const Fem::Field::CField& field1 = world.GetField(id_field_disp1);  
	const Fem::Field::CNodeAry::CNodeSeg& ns_c1 = field1.GetNodeSeg(CORNER,false,world,VALUE);   
  const Fem::Field::CNodeAry::CNodeSeg& ns_u1 = field1.GetNodeSeg(CORNER,true, world,VALUE); 
  
  for(unsigned int ino0=0;ino0<ns_c0.Size();ino0++)
  {
    double C0[3]; ns_c0.GetValue(ino0,C0);
    double u0[3]; ns_u0.GetValue(ino0,u0);
    double c0[3] = { C0[0]+u0[0], C0[1]+u0[1], C0[2]+u0[2] };
    const unsigned int ino1_0 = aInterp[ino0].ino1;
    const unsigned int ino1_1 = aInterp[ino0].ino2;
    const unsigned int ino1_2 = aInterp[ino0].ino3;
    double u1[3][3];	ns_u1.GetValue(ino1_0,u1[0]);	ns_u1.GetValue(ino1_1,u1[1]); ns_u1.GetValue(ino1_2,u1[2]);  
    double C1[3][3];	ns_c1.GetValue(ino1_0,C1[0]);	ns_c1.GetValue(ino1_1,C1[1]); ns_c1.GetValue(ino1_2,C1[2]);    
    const double c1[3][3] = {
      { u1[0][0]+C1[0][0], u1[0][1]+C1[0][1], u1[0][2]+C1[0][2] },
      { u1[1][0]+C1[1][0], u1[1][1]+C1[1][1], u1[1][2]+C1[1][2] },
      { u1[2][0]+C1[2][0], u1[2][1]+C1[2][1], u1[2][2]+C1[2][2] } };    
    double n1[3],a1; Com::UnitNormalAreaTri3D(n1, a1, c1[0], c1[1], c1[2]);
    double nh0[3]; Com::NormalTri3D(nh0, c0,    c1[1], c1[2]);
    double nh1[3]; Com::NormalTri3D(nh1, c1[0], c0,    c1[2]);
    double nh2[3]; Com::NormalTri3D(nh2, c1[0], c1[1], c0   );
    double r0 = Com::Dot3D(nh0, n1);
    double r1 = Com::Dot3D(nh1, n1);
    double r2 = Com::Dot3D(nh2, n1);
    double inv_sum_r = 1.0/(r0+r1+r2);
    r0 *= inv_sum_r;
    r1 *= inv_sum_r;
    r2 *= inv_sum_r;    
    aInterp[ino0].r1 = r0;
    aInterp[ino0].r2 = r1;
    double a[3] = {
      c1[0][0]*r0 + c1[1][0]*r1 + c1[2][0]*r2,
      c1[0][1]*r0 + c1[1][1]*r1 + c1[2][1]*r2,
      c1[0][2]*r0 + c1[1][2]*r1 + c1[2][2]*r2 };
    double hv[3] = { c0[0]-a[0], c0[1]-a[1], c0[2]-a[2] };
    aInterp[ino0].height = Com::Dot3D(hv,n1);
  }      
}

void InitFineDeformInterp(std::vector<CInterpBarycentric>& aInterp)
{
  for(unsigned int itp=0;itp<aInterp.size();itp++)
  {
    aInterp[itp].height = 0;
    aInterp[itp].r1 = aInterp[itp].r1_ini;
    aInterp[itp].r2 = aInterp[itp].r2_ini;
  }
}

bool StepTime_Static
(double& dt,
 double& total_energy,
 double torelance_static,
 bool& is_ilufrac_success,
 Fem::Ls::CLinearSystem_Field& ls, LsSol::CPreconditioner_ILU& prec,
 ////
 const CClothParam& cloth_param,
 double g_x, double g_y, double g_z,
 ////
 const CContactParam& contact_param, 
 const CContactTarget3D& CT,
 std::vector<CFrictionPoint>& aFrictionPoint, 
 /////
 const CStitchAry& stitch_ary,  
 /////
 unsigned int id_field_disp, 
 unsigned int id_field_hinge, 
 unsigned int id_field_disp_buffer,
 Fem::Field::CFieldWorld& world)
{
  //	std::cout << "Step Time BackwardEular" << std::endl;  
	ls.InitializeMarge();
  
  static int icnt = 0;
  static double ke0=0,ke1=0,ke2=0;
  static double te0=0,te1=0,te2=0;
  double ke=0,se=0,pe=0;
//	AddLinSys_DiscreteShell_CST_BackWardEular
	AddLinSys_QuadBend_CST_BackWardEular  
	(dt, 0, ls,
	 cloth_param.stiff_bend, cloth_param.stiff_lambda, cloth_param.stiff_myu, 
	 cloth_param.rho, 
   g_x, g_y, g_z,
	 id_field_disp, id_field_hinge,
	 world,
   ke,se,pe);
  stitch_ary.AddLinSys_BackwardEular(dt,se,ls,id_field_disp,world);  
  //  std::cout << "kine energy" << ke << " " << icnt << "   strain:" << se << "    tot:" << ke+se+pe << std::endl;
  double te = ke+se+pe;
  ke0=ke1;  ke1=ke2; ke2=ke;  
  te0=te1;  te1=te2; te2=te;
  total_energy = te;  
  icnt++;
  if( icnt > 10 && KineticDamping(ke0, ke1, ke2, id_field_disp, world) ){
    icnt=0;
    return false;
  }  
  if( te2 < te1 && te1 < te0 ){
    CopyValueVelo(id_field_disp_buffer, id_field_disp, world);        
  }
  
	{		
		AddLinSys_FrictionalContact_Penalty_NonStatic_BackwardEular
		(dt,
		 ls,
		 CT,
     contact_param.stiff_n, contact_param.stiff_f, contact_param.myu_s, contact_param.myu_k,
		 contact_param.offset,
		 id_field_disp,
		 world,
		 aFrictionPoint);
	}
  
//	const double res_norm = ls.FinalizeMarge();
	is_ilufrac_success = prec.SetValue(ls.m_ls);
  if( !is_ilufrac_success ){ 
    std::cout << "ilu frac false in update solution" << std::endl;
//    getchar();
    return false;
  }
  
  //	std::cout << "Residual : " << res_norm << std::endl;
  unsigned int iter;
  const unsigned int max_iter = 50;  
	{
		double tol = 1.0e-6;
		iter = max_iter;
		LsSol::CLinearSystemPreconditioner lsprec(ls.m_ls,prec);
		LsSol::Solve_PCG(tol,iter,lsprec);
		//		LsSol::Solve_CG(tol,iter,ls.m_ls);
//    std::cout << "  dt " << dt << "   "  << iter << " " << tol << " total eng  " << te2 << std::endl;
    if( iter >= max_iter ){
      //      SaveValue(id_field_disp,id_field_disp_buffer,world);
      dt = dt*0.5;
      std::cout << "&&&&&&&&&&&&&&&&&" << std::endl;
      std::cout << "cut time step half" << dt << std::endl;      
      return false;
    }
	}
	ls.UpdateValueOfField_BackwardEular(dt,   id_field_disp,world,true);
  //	std::cout << "Step Time BackwardEular 10" << std::endl;
	{
		Update_FrictionalContact
		(CT, 
     contact_param.stiff_n, contact_param.stiff_f, contact_param.myu_s, contact_param.myu_k,
		 contact_param.offset,
		 id_field_disp,
		 world,
		 aFrictionPoint);
	}		
  if( iter < 10
     && te0 > te1  
     && te1 > te2 
     && dt < 0.01 )
  {    
    std::cout << "#################" << std::endl;
    std::cout << "add time step double" << dt << std::endl;          
    dt=dt*1.2;
  }  
  const unsigned int nnode = aFrictionPoint.size();
	const double upd_per_node = sqrt(ls.DOT(-2,-2))/nnode;
//  std::cout << "hoge : " << upd_per_node << " " << torelance_static << std::endl;
  if( upd_per_node < torelance_static ){
    std::cout << "converge" << std::endl;
    return true;
  }
  return false;
}

bool GetSensitivity_fictbend
(Fem::Ls::CLinearSystem_Field& ls, LsSol::CPreconditioner_ILU& prec,
 ////
 const CClothParam& cloth_param,
 double g_x, double g_y, double g_z,
 ////
 const CContactParam& contact_param, 
 const CContactTarget3D& CT,
 std::vector<CFrictionPoint>& aFrictionPoint, 
 /////
 const CStitchAry& stitch_ary,    
// double stitch_stiff, double stitch_damp_coeff,
// const std::vector<unsigned int>& aIdFieldDart,
 /////
 unsigned int id_field_disp, 
 unsigned int id_field_hinge, 
 /////
 bool is_xy,
 unsigned int id_field_senseX, unsigned int id_field_senseY,
 unsigned int id_field_lamX,  unsigned int id_field_lamY,
 Fem::Field::CFieldWorld& world)
{  
  double fict_bend = 0.0001*(cloth_param.stiff_lambda+cloth_param.stiff_myu);
  fict_bend = ( fict_bend > cloth_param.stiff_bend ) ? fict_bend : cloth_param.stiff_bend;
//  fict_bend = 0;
  ////
	ls.InitializeMarge();
//	AddLinSys_DiscreteShell_CST_Sensitivity_FictitousBending
	AddLinSys_QuadBend_CST_Sensitivity_FictitousBending  
	(ls,
	 fict_bend, cloth_param.stiff_lambda, cloth_param.stiff_myu,
	 cloth_param.rho, g_x, g_y, g_z,
	 id_field_disp, id_field_hinge,
   id_field_senseX, id_field_senseY,
   id_field_lamX, id_field_lamY,
	 world );
  
  stitch_ary.AddLinSys_Sensitivity(ls,id_field_disp,world);  
  AddLinSys_FrictionalContact_Penalty_NonStatic_Sensitivity
  (ls,
   CT,
   contact_param.stiff_n, contact_param.stiff_f*0, contact_param.myu_s, contact_param.myu_k,
   contact_param.offset,
   id_field_disp,
   world,
   aFrictionPoint);

//  const double res_norm = ls.FinalizeMarge();
  MatVec::CVector_Blk rvec = ls.m_ls.GetVector(-1,0); 
  const bool is_ilufrac_success = prec.SetValue(ls.m_ls);
  if( !is_ilufrac_success ){
    std::cout << "ilu frac false in sensitivity analysis" << std::endl;
    return false;
  }
  ////    
  { // XSensitive 
    { // set x residual
      Fem::Field::CField& field_sx = world.GetField(id_field_senseX);
      Fem::Field::CNodeAry::CNodeSeg& ns_sx = field_sx.GetNodeSeg(CORNER,true, world,VALUE);      
      MatVec::CVector_Blk& res = ls.GetResidual(id_field_disp,CORNER,world);
      const unsigned int nnode = ns_sx.Size();
      assert( nnode == res.NBlk() );      
      for(unsigned int ino=0;ino<nnode;ino++){
        double v[3]; ns_sx.GetValue(ino,v);
        res.SetValue(ino,0,v[0]);
        res.SetValue(ino,1,v[1]);
        res.SetValue(ino,2,v[2]);        
      }
    }
//    const double res_norm = ls.FinalizeMarge();    
    //	std::cout << "Residual : " << res_norm << std::endl;
    {
      double tol = 1.0e-4;
      unsigned int iter = 1000;
      LsSol::CLinearSystemPreconditioner lsprec(ls.m_ls,prec);
      LsSol::Solve_PCG(tol,iter,lsprec);
//      		LsSol::Solve_CG(tol,iter,ls.m_ls);
			std::cout << iter << " " << tol << std::endl;      
    }	
    {
      const MatVec::CVector_Blk& upd = ls.GetUpdate(id_field_disp,CORNER,world);      
      Fem::Field::CField& field_sx = world.GetField(id_field_senseX);
      Fem::Field::CNodeAry::CNodeSeg& ns_sx = field_sx.GetNodeSeg(CORNER,true, world,VALUE);
      const unsigned int nnode = ns_sx.Size();
      assert( nnode == upd.NBlk() );
//      std::cout << upd.GetValue(0,2) << " " << upd.GetValue(1,2) << std::endl;
      for(unsigned int ino=0;ino<nnode;ino++){
        ns_sx.SetValue(ino,0,upd.GetValue(ino,0));
        ns_sx.SetValue(ino,1,upd.GetValue(ino,1));
        ns_sx.SetValue(ino,2,upd.GetValue(ino,2));        
      }                  
    }
  }
  if( is_xy ){ // YSensitive 
    { // set y residual      
      Fem::Field::CField& field_sy = world.GetField(id_field_senseY);
      Fem::Field::CNodeAry::CNodeSeg& ns_sy = field_sy.GetNodeSeg(CORNER,true, world,VALUE);
      MatVec::CVector_Blk& res = ls.GetResidual(id_field_disp,CORNER,world);      
      const unsigned int nnode = ns_sy.Size();
      assert( nnode == res.NBlk() );
      for(unsigned int ino=0;ino<nnode;ino++){
        double v[3]; ns_sy.GetValue(ino,v);
//        res.SetValue(ino,0,v[0]+rvec.GetValue(ino,0));
//        res.SetValue(ino,1,v[1]+rvec.GetValue(ino,1));
//        res.SetValue(ino,2,v[2]+rvec.GetValue(ino,2));
        res.SetValue(ino,0,v[0]);
        res.SetValue(ino,1,v[1]);
        res.SetValue(ino,2,v[2]);        
      }      
    }
    //	std::cout << "Residual : " << res_norm << std::endl;
//    const double res_norm = ls.FinalizeMarge();        
    {
      double tol = 1.0e-4;
      unsigned int iter = 1000;
      LsSol::CLinearSystemPreconditioner lsprec(ls.m_ls,prec);
      LsSol::Solve_PCG(tol,iter,lsprec);
      //		LsSol::Solve_CG(tol,iter,ls.m_ls);
			std::cout << iter << " " << tol << std::endl;            
    }	
    {
      const MatVec::CVector_Blk& upd = ls.GetUpdate(id_field_disp,CORNER,world);            
      Fem::Field::CField& field_sy = world.GetField(id_field_senseY);
      Fem::Field::CNodeAry::CNodeSeg& ns_sy = field_sy.GetNodeSeg(CORNER,true, world,VALUE);      
      const unsigned int nnode = ns_sy.Size();
      assert( nnode == upd.NBlk() );
      for(unsigned int ino=0;ino<nnode;ino++){
        ns_sy.SetValue(ino,0,upd.GetValue(ino,0));
        ns_sy.SetValue(ino,1,upd.GetValue(ino,1));
        ns_sy.SetValue(ino,2,upd.GetValue(ino,2));        
      }            
    }
  }  
  return true;
}




void NoResponse
(double mov_x, double mov_y,
 unsigned int id_field_disp,
 unsigned int id_field_reference,
 unsigned int id_field_lamX,  unsigned int id_field_lamY,
 Fem::Field::CFieldWorld& world )
{
  Fem::Field::CField& field_ref  = world.GetField(id_field_reference);  
  Fem::Field::CNodeAry::CNodeSeg& ns_ref = field_ref.GetNodeSeg( CORNER,true, world,VALUE); 
  Fem::Field::CField& field_lx  = world.GetField(id_field_lamX); 
  Fem::Field::CField& field_ly  = world.GetField(id_field_lamY);   
  const Fem::Field::CNodeAry::CNodeSeg& ns_lx = field_lx.GetNodeSeg( CORNER,true, world,VALUE);  
  const Fem::Field::CNodeAry::CNodeSeg& ns_ly = field_ly.GetNodeSeg( CORNER,true, world,VALUE);    
  unsigned int nno = ns_lx.Size();
  assert( ns_ly.Size() == nno );
  for(unsigned int ino=0;ino<nno;ino++){
    double lx[2];  ns_lx.GetValue(ino,lx);
    double ly[2];  ns_lx.GetValue(ino,ly);    
    double cr[3];  ns_ref.GetValue(ino,cr);
    ns_ref.SetValue(ino, 0, cr[0]-lx[0]*mov_x-ly[0]*mov_y);
    ns_ref.SetValue(ino, 1, cr[1]-lx[1]*mov_x-ly[1]*mov_y); 
  }    
  Fem::Field::CField& field_disp = world.GetField(id_field_disp);
  Fem::Field::CNodeAry::CNodeSeg& ns_u   = field_disp.GetNodeSeg(CORNER,true, world,VALUE);     
  for(unsigned int ino=0;ino<ns_u.Size();ino++){
    double cr[3];  ns_u.GetValue(ino,cr);
    double lx[2];  ns_lx.GetValue(ino,lx);    
    double ly[2];  ns_ly.GetValue(ino,ly);        
    ns_u.SetValue(ino,0,cr[0]-lx[0]*mov_x-ly[0]*mov_y);
    ns_u.SetValue(ino,1,cr[1]-lx[1]*mov_x-ly[1]*mov_y);
    ns_u.SetValue(ino,2,cr[2]);    
  }
}


void NoResponse_Slider
(double mov_v,
 unsigned int id_field_disp,
 unsigned int id_field_reference,
 unsigned int id_field_lamX,
 Fem::Field::CFieldWorld& world )
{
  Fem::Field::CField& field_ref  = world.GetField(id_field_reference);  
  Fem::Field::CNodeAry::CNodeSeg& ns_ref = field_ref.GetNodeSeg( CORNER,true, world,VALUE); 
  Fem::Field::CField& field_lx  = world.GetField(id_field_lamX); 
  const Fem::Field::CNodeAry::CNodeSeg& ns_lx = field_lx.GetNodeSeg( CORNER,true, world,VALUE);  
  unsigned int nno = ns_lx.Size();
  for(unsigned int ino=0;ino<nno;ino++){
    double lx[2];  ns_lx.GetValue(ino,lx);
    double cr[3];  ns_ref.GetValue(ino,cr);
    ns_ref.SetValue(ino, 0, cr[0]-lx[0]*mov_v);
    ns_ref.SetValue(ino, 1, cr[1]-lx[1]*mov_v); 
  }    
  Fem::Field::CField& field_disp = world.GetField(id_field_disp);
  Fem::Field::CNodeAry::CNodeSeg& ns_u   = field_disp.GetNodeSeg(CORNER,true, world,VALUE);     
  for(unsigned int ino=0;ino<ns_u.Size();ino++){
    double cr[3];  ns_u.GetValue(ino,cr);
    double lx[2];  ns_lx.GetValue(ino,lx);    
    ns_u.SetValue(ino,0,cr[0]-lx[0]*mov_v);
    ns_u.SetValue(ino,1,cr[1]-lx[1]*mov_v);
    ns_u.SetValue(ino,2,cr[2]);    
  }
}




void SensitiveResponse
(double mov_x, double mov_y,
 unsigned int id_field_disp,
 unsigned int id_field_reference,
 unsigned int id_field_lamX,  unsigned int id_field_lamY,
 unsigned int id_field_senseX, unsigned int id_field_senseY,
 Fem::Field::CFieldWorld& world )
{
  Fem::Field::CField& field_ref  = world.GetField(id_field_reference);  
  Fem::Field::CNodeAry::CNodeSeg& ns_ref = field_ref.GetNodeSeg( CORNER,true, world,VALUE); 
  Fem::Field::CField& field_lx  = world.GetField(id_field_lamX); 
  Fem::Field::CField& field_ly  = world.GetField(id_field_lamY);   
  const Fem::Field::CNodeAry::CNodeSeg& ns_lx = field_lx.GetNodeSeg( CORNER,true, world,VALUE);  
  const Fem::Field::CNodeAry::CNodeSeg& ns_ly = field_ly.GetNodeSeg( CORNER,true, world,VALUE);    
  unsigned int nno = ns_lx.Size();
  assert( ns_ly.Size() == nno );
  for(unsigned int ino=0;ino<nno;ino++){
    double lx[2];  ns_lx.GetValue(ino,lx);
    double ly[2];  ns_lx.GetValue(ino,ly);    
    double cr[3];  ns_ref.GetValue(ino,cr);
    ns_ref.SetValue(ino, 0, cr[0]-lx[0]*mov_x-ly[0]*mov_y);
    ns_ref.SetValue(ino, 1, cr[1]-lx[1]*mov_x-ly[1]*mov_y); 
  }  
  Fem::Field::CField& field_disp = world.GetField(id_field_disp);
  Fem::Field::CNodeAry::CNodeSeg& ns_u   = field_disp.GetNodeSeg(CORNER,true, world,VALUE);     
  Fem::Field::CField& field_sx = world.GetField(id_field_senseX);
  const Fem::Field::CNodeAry::CNodeSeg& ns_sx  = field_sx.GetNodeSeg(  CORNER,true, world,VALUE);
  Fem::Field::CField& field_sy = world.GetField(id_field_senseY);
  const Fem::Field::CNodeAry::CNodeSeg& ns_sy  = field_sy.GetNodeSeg(  CORNER,true, world,VALUE);        
  for(unsigned int ino=0;ino<ns_u.Size();ino++){
    double sx[3];  ns_sx.GetValue(ino,sx);  
    double sy[3];  ns_sy.GetValue(ino,sy);      
    double cr[3];  ns_u.GetValue(ino,cr);
    double lx[2];  ns_lx.GetValue(ino,lx);    
    double ly[2];  ns_ly.GetValue(ino,ly);        
    ns_u.SetValue(ino,0,cr[0]+sx[0]*mov_x+sy[0]*mov_y-lx[0]*mov_x-ly[0]*mov_y);
    ns_u.SetValue(ino,1,cr[1]+sx[1]*mov_x+sy[1]*mov_y-lx[1]*mov_x-ly[1]*mov_y);
    ns_u.SetValue(ino,2,cr[2]+sx[2]*mov_x+sy[2]*mov_y);
  }
}



void SensitiveResponse_Slider
(double mov_v,
 unsigned int id_field_disp,
 unsigned int id_field_reference,
 unsigned int id_field_lamX,
 unsigned int id_field_senseX,
 Fem::Field::CFieldWorld& world )
{
  Fem::Field::CField& field_ref  = world.GetField(id_field_reference);  
  Fem::Field::CNodeAry::CNodeSeg& ns_ref = field_ref.GetNodeSeg( CORNER,true, world,VALUE); 
  Fem::Field::CField& field_lx  = world.GetField(id_field_lamX); 
  const Fem::Field::CNodeAry::CNodeSeg& ns_lx = field_lx.GetNodeSeg( CORNER,true, world,VALUE);  
  unsigned int nno = ns_lx.Size();
  for(unsigned int ino=0;ino<nno;ino++){
    double lx[2];  ns_lx.GetValue(ino,lx);
    double cr[3];  ns_ref.GetValue(ino,cr);
    ns_ref.SetValue(ino, 0, cr[0]-lx[0]*mov_v);
    ns_ref.SetValue(ino, 1, cr[1]-lx[1]*mov_v); 
  }  
  Fem::Field::CField& field_disp = world.GetField(id_field_disp);
  Fem::Field::CNodeAry::CNodeSeg& ns_u   = field_disp.GetNodeSeg(CORNER,true, world,VALUE);     
  Fem::Field::CField& field_sx = world.GetField(id_field_senseX);
  const Fem::Field::CNodeAry::CNodeSeg& ns_sx  = field_sx.GetNodeSeg(  CORNER,true, world,VALUE);
//  Fem::Field::CField& field_sy = world.GetField(id_field_senseY);
//  const Fem::Field::CNodeAry::CNodeSeg& ns_sy  = field_sy.GetNodeSeg(  CORNER,true, world,VALUE);        
  for(unsigned int ino=0;ino<ns_u.Size();ino++){
    double sx[3];  ns_sx.GetValue(ino,sx);  
//    double sy[3];  ns_sy.GetValue(ino,sy);      
    double cr[3];  ns_u.GetValue(ino,cr);
    double lx[2];  ns_lx.GetValue(ino,lx);    
//    double ly[2];  ns_ly.GetValue(ino,ly);        
    ns_u.SetValue(ino,0,cr[0]+sx[0]*mov_v-lx[0]*mov_v);
    ns_u.SetValue(ino,1,cr[1]+sx[1]*mov_v-lx[1]*mov_v);
    ns_u.SetValue(ino,2,cr[2]+sx[2]*mov_v);
  }
}

static inline void CalcInvMat2(double a[])
{
	const double det = a[0]*a[3] - a[1]*a[2];
	const double inv_det = 1.0/det;
  
  const double t[4] = { a[0], a[1], a[2], a[3] };
  
	a[0] = inv_det*(+t[3]);
	a[1] = inv_det*(-t[1]);  
	a[2] = inv_det*(-t[2]);  
	a[3] = inv_det*(+t[0]);
}
  


static inline void CalcInvMat3(double a[])
{
	const double det = a[0]*a[4]*a[8] + a[3]*a[7]*a[2] + a[6]*a[1]*a[5]
  - a[0]*a[7]*a[5] - a[6]*a[4]*a[2] - a[3]*a[1]*a[8];
	const double inv_det = 1.0/det;
  
  double t[9];  for(int i=0;i<9;i++){ t[i] = a[i]; }  // copy a to t
  
	a[0] = inv_det*(t[4]*t[8]-t[5]*t[7]);
	a[1] = inv_det*(t[2]*t[7]-t[1]*t[8]);
	a[2] = inv_det*(t[1]*t[5]-t[2]*t[4]);
  
	a[3] = inv_det*(t[5]*t[6]-t[3]*t[8]);
	a[4] = inv_det*(t[0]*t[8]-t[2]*t[6]);
	a[5] = inv_det*(t[2]*t[3]-t[0]*t[5]);
  
	a[6] = inv_det*(t[3]*t[7]-t[4]*t[6]);
	a[7] = inv_det*(t[1]*t[6]-t[0]*t[7]);
	a[8] = inv_det*(t[0]*t[4]-t[1]*t[3]);
}

void GuessSolution_GMLS_Slider
(double pre_v, double pos_v,
 unsigned int id_field_disp,
 unsigned int id_field_reference,
 unsigned int id_field_lamX,
 const std::vector<CSolutionSensitivity>& aSolSens,
 Fem::Field::CFieldWorld& world )
{
  double mov_v = pos_v-pre_v;
  Fem::Field::CField& field_ref  = world.GetField(id_field_reference);  
  Fem::Field::CNodeAry::CNodeSeg& ns_ref = field_ref.GetNodeSeg( CORNER,true, world,VALUE); 
  Fem::Field::CField& field_lx  = world.GetField(id_field_lamX); 
  const Fem::Field::CNodeAry::CNodeSeg& ns_lx = field_lx.GetNodeSeg( CORNER,true, world,VALUE);  
  unsigned int nno = ns_lx.Size();
  for(unsigned int ino=0;ino<nno;ino++){
    double lx[2];  ns_lx.GetValue(ino,lx);
    double cr[3];  ns_ref.GetValue(ino,cr);
    ns_ref.SetValue(ino, 0, cr[0]-lx[0]*mov_v);
    ns_ref.SetValue(ino, 1, cr[1]-lx[1]*mov_v); 
  }  
  double nr;
  std::vector<double> aNxdxdy;
  {
    aNxdxdy.resize(aSolSens.size()*2,0);
    std::vector<double> aTmp;
    aTmp.resize(aSolSens.size()*3,0);
    double pc[2] = { 1, pos_v };
    double pr[2] = { 1, pre_v };
    const double sqdist_r = (pos_v-pre_v)*(pos_v-pre_v);
    double wr = 1.0/(sqdist_r+0.001);
    for(unsigned int iss=0;iss<aSolSens.size();iss++){
      if( !aSolSens[iss].is_active ) continue;
      double obj_v = aSolSens[iss].val_slider;
      aTmp[iss*3+0] = 1;
      aTmp[iss*3+1] = obj_v;
      const double sqdist = (pos_v-obj_v)*(pos_v-obj_v);
      aTmp[iss*3+2] = 1.0/(sqdist+0.001);
    }
    double deriv_ratio = 1.0;
    double invG[4];
    {
      double G[4] = { 0,0,0,0 };
      for(unsigned int i=0;i<2;i++){
        for(unsigned int j=0;j<2;j++){
          G[i*2+j] += pr[i]*pr[j]*wr;
        }
      }      
      for(unsigned int iss=0;iss<aSolSens.size();iss++){
        if( !aSolSens[iss].is_active ) continue;
        double w0 = aTmp[iss*3+2];
        for(unsigned int i=0;i<2;i++){
          for(unsigned int j=0;j<2;j++){
            G[i*2+j] += aTmp[iss*3+i]*aTmp[iss*3+j]*w0;
          }
        }
        G[3] += w0*deriv_ratio;
      }
      for(unsigned int i=0;i<4;i++){ invG[i] = G[i]; }
      CalcInvMat2(invG);
    }    
    nr = 0;
    for(unsigned int i=0;i<2;i++){
      for(unsigned int j=0;j<2;j++){      
        nr += pc[i]*invG[i*2+j]*pr[j]*wr;
      }
    }        
    for(unsigned int iss=0;iss<aSolSens.size();iss++){
      if( !aSolSens[iss].is_active ) continue;
      const double w0 = aTmp[iss*3+2];   
      for(unsigned int i=0;i<2;i++){
        for(unsigned int j=0;j<2;j++){      
          aNxdxdy[iss*2+0] += pc[i]*invG[i*2+j]*aTmp[iss*3+j]*w0;
        }
      }
      for(unsigned int i=0;i<2;i++){    
        aNxdxdy[iss*2+1] += pc[i]*invG[i*2+1]*w0*deriv_ratio;
      }
    }
  }
  //  std::cout << nr << " " << n0 << " " << n0x << " " << n0y << std::endl;
  //  getchar();
  
  {
    Fem::Field::CField& field_disp = world.GetField(id_field_disp);
    Fem::Field::CNodeAry::CNodeSeg& ns_u   = field_disp.GetNodeSeg(CORNER,true, world,VALUE);     
    Fem::Field::CNodeAry::CNodeSeg& ns_X   = field_disp.GetNodeSeg(CORNER,false,world,VALUE);        
    for(unsigned int ino=0;ino<ns_u.Size();ino++){
      double cu[3];  ns_u.GetValue( ino,cu);
      double cX[3];  ns_X.GetValue( ino,cX);    
      double lx[2];  ns_lx.GetValue(ino,lx);    
      double cx[3] = {      
        cX[0]-lx[0]*mov_v+cu[0],
        cX[1]-lx[1]*mov_v+cu[1],
        cX[2]            +cu[2] };
      double x[3] = {
        nr*cx[0],
        nr*cx[1],
        nr*cx[2] };
      double u[3] = {
        x[0]-cX[0],
        x[1]-cX[1],
        x[2]-cX[2] };
      ns_u.SetValue(ino,0,u[0]);
      ns_u.SetValue(ino,1,u[1]);
      ns_u.SetValue(ino,2,u[2]);
    }          
  }
  
  for(unsigned int iss=0;iss<aSolSens.size();iss++){
    if( !aSolSens[iss].is_active ) continue;
    Fem::Field::CField& field_disp = world.GetField(id_field_disp);
    Fem::Field::CNodeAry::CNodeSeg& ns_u   = field_disp.GetNodeSeg(CORNER,true, world,VALUE);     
    
    unsigned int id_field_senseX = aSolSens[iss].id_field_dudpx;
    unsigned int id_field_x      = aSolSens[iss].id_field_x;
    ////
    Fem::Field::CField& field_sx = world.GetField(id_field_senseX);
    const Fem::Field::CNodeAry::CNodeSeg& ns_sx  = field_sx.GetNodeSeg( CORNER,true, world,VALUE);
    Fem::Field::CField& field_x = world.GetField(id_field_x);
    const Fem::Field::CNodeAry::CNodeSeg& ns_x0  = field_x.GetNodeSeg(  CORNER,true, world,VALUE);        
    ////    
    for(unsigned int ino=0;ino<ns_x0.Size();ino++){
      double x0[3];  ns_x0.GetValue(ino,x0);    
      double sx[3];  ns_sx.GetValue(ino,sx);  
      double x[3] = {
        aNxdxdy[iss*2+0]*x0[0] + aNxdxdy[iss*2+1]*sx[0],
        aNxdxdy[iss*2+0]*x0[1] + aNxdxdy[iss*2+1]*sx[1],
        aNxdxdy[iss*2+0]*x0[2] + aNxdxdy[iss*2+1]*sx[2] };
      ns_u.AddValue(ino,0,x[0]);
      ns_u.AddValue(ino,1,x[1]);
      ns_u.AddValue(ino,2,x[2]);      
    }
  }    
}


void GuessSolution_GMLS
(double pre_x, double pre_y, 
 double pos_x, double pos_y,
 unsigned int id_field_disp,
 unsigned int id_field_reference,
 unsigned int id_field_lamX,  unsigned int id_field_lamY,
 const std::vector<CSolutionSensitivity>& aSolSens,
 Fem::Field::CFieldWorld& world )
{
//  std::cout << "GuessSolution_GMLS" << std::endl;
  double mov_x = pos_x-pre_x;
  double mov_y = pos_y-pre_y;
  Fem::Field::CField& field_ref  = world.GetField(id_field_reference);  
  Fem::Field::CNodeAry::CNodeSeg& ns_ref = field_ref.GetNodeSeg( CORNER,true, world,VALUE); 
  Fem::Field::CField& field_lx  = world.GetField(id_field_lamX); 
  Fem::Field::CField& field_ly  = world.GetField(id_field_lamY);   
  const Fem::Field::CNodeAry::CNodeSeg& ns_lx = field_lx.GetNodeSeg( CORNER,true, world,VALUE);  
  const Fem::Field::CNodeAry::CNodeSeg& ns_ly = field_ly.GetNodeSeg( CORNER,true, world,VALUE);    
  unsigned int nno = ns_lx.Size();
  assert( ns_ly.Size() == nno );
  for(unsigned int ino=0;ino<nno;ino++){
    double lx[2];  ns_lx.GetValue(ino,lx);
    double ly[2];  ns_lx.GetValue(ino,ly);    
    double cr[3];  ns_ref.GetValue(ino,cr);
    ns_ref.SetValue(ino, 0, cr[0]-lx[0]*mov_x-ly[0]*mov_y);
    ns_ref.SetValue(ino, 1, cr[1]-lx[1]*mov_x-ly[1]*mov_y); 
  }
  double nr;
  std::vector<double> aNxdxdy;
  {
    aNxdxdy.resize(aSolSens.size()*3,0);
    std::vector<double> aTmp;
    aTmp.resize(aSolSens.size()*4,0);
    double pc[3] = { 1, pos_x, pos_y };
    double pr[3] = { 1, pre_x, pre_y };
    const double sqdist_r = (pos_x-pre_x)*(pos_x-pre_x) + (pos_y-pre_y)*(pos_y-pre_y);
    double wr = 1.0/(sqdist_r+0.001);
    for(unsigned int iss=0;iss<aSolSens.size();iss++){
      if( !aSolSens[iss].is_active ) continue;
      double obj_x = aSolSens[iss].obj_x;
      double obj_y = aSolSens[iss].obj_y;
      aTmp[iss*4+0] = 1;
      aTmp[iss*4+1] = obj_x;
      aTmp[iss*4+2] = obj_y;
      const double sqdist = (pos_x-obj_x)*(pos_x-obj_x) + (pos_y-obj_y)*(pos_y-obj_y);
      aTmp[iss*4+3] = 1.0/(sqdist+0.001);
    }
    double deriv_ratio = 1.0;
    double invG[9];
    {
      double G[9] = { 0,0,0,0,0,0,0,0,0 };
      for(unsigned int i=0;i<3;i++){
      for(unsigned int j=0;j<3;j++){
        G[i*3+j] += pr[i]*pr[j]*wr;
      }
      }      
      for(unsigned int iss=0;iss<aSolSens.size();iss++){
        if( !aSolSens[iss].is_active ) continue;
        double w0 = aTmp[iss*4+3];
        for(unsigned int i=0;i<3;i++){
        for(unsigned int j=0;j<3;j++){
          G[i*3+j] += aTmp[iss*4+i]*aTmp[iss*4+j]*w0;
        }
        }
        G[4] += w0*deriv_ratio;
        G[8] += w0*deriv_ratio;      
      }
      for(unsigned int i=0;i<9;i++){ invG[i] = G[i]; }
      CalcInvMat3(invG);
    }    
    nr = 0;
    for(unsigned int i=0;i<3;i++){
    for(unsigned int j=0;j<3;j++){      
      nr += pc[i]*invG[i*3+j]*pr[j]*wr;
    }
    }        
    for(unsigned int iss=0;iss<aSolSens.size();iss++){
      if( !aSolSens[iss].is_active ) continue;
      const double w0 = aTmp[iss*4+3];                          
      for(unsigned int i=0;i<3;i++){
      for(unsigned int j=0;j<3;j++){      
        aNxdxdy[iss*3+0] += pc[i]*invG[i*3+j]*aTmp[iss*4+j]*w0;
      }
      }
      for(unsigned int i=0;i<3;i++){    
        aNxdxdy[iss*3+1] += pc[i]*invG[i*3+1]*w0*deriv_ratio;
      }
      for(unsigned int i=0;i<3;i++){    
        aNxdxdy[iss*3+2] += pc[i]*invG[i*3+2]*w0*deriv_ratio;
      }    
    }
  }
//  std::cout << nr << " " << n0 << " " << n0x << " " << n0y << std::endl;
//  getchar();
  
  {
    Fem::Field::CField& field_disp = world.GetField(id_field_disp);
    Fem::Field::CNodeAry::CNodeSeg& ns_u   = field_disp.GetNodeSeg(CORNER,true, world,VALUE);     
    Fem::Field::CNodeAry::CNodeSeg& ns_X   = field_disp.GetNodeSeg(CORNER,false,world,VALUE);        
    for(unsigned int ino=0;ino<ns_u.Size();ino++){
      double cu[3];  ns_u.GetValue( ino,cu);
      double cX[3];  ns_X.GetValue( ino,cX);    
      double lx[2];  ns_lx.GetValue(ino,lx);    
      double ly[2];  ns_ly.GetValue(ino,ly);        
      double cx[3] = {      
        cX[0]-lx[0]*mov_x-ly[0]*mov_y+cu[0],
        cX[1]-lx[1]*mov_x-ly[1]*mov_y+cu[1],
        cX[2]                        +cu[2] };
      double x[3] = {
        nr*cx[0],
        nr*cx[1],
        nr*cx[2] };
      double u[3] = {
        x[0]-cX[0],
        x[1]-cX[1],
        x[2]-cX[2] };
      ns_u.SetValue(ino,0,u[0]);
      ns_u.SetValue(ino,1,u[1]);
      ns_u.SetValue(ino,2,u[2]);
    }          
  }
  
  for(unsigned int iss=0;iss<aSolSens.size();iss++){
    if( !aSolSens[iss].is_active ) continue;
    Fem::Field::CField& field_disp = world.GetField(id_field_disp);
    Fem::Field::CNodeAry::CNodeSeg& ns_u   = field_disp.GetNodeSeg(CORNER,true, world,VALUE);     
    
    unsigned int id_field_senseX = aSolSens[iss].id_field_dudpx;
    unsigned int id_field_senseY = aSolSens[iss].id_field_dudpy;  
    unsigned int id_field_x      = aSolSens[iss].id_field_x;
    ////
    Fem::Field::CField& field_sx = world.GetField(id_field_senseX);
    const Fem::Field::CNodeAry::CNodeSeg& ns_sx  = field_sx.GetNodeSeg( CORNER,true, world,VALUE);
    Fem::Field::CField& field_sy = world.GetField(id_field_senseY);
    const Fem::Field::CNodeAry::CNodeSeg& ns_sy  = field_sy.GetNodeSeg( CORNER,true, world,VALUE);        
    Fem::Field::CField& field_x = world.GetField(id_field_x);
    const Fem::Field::CNodeAry::CNodeSeg& ns_x0  = field_x.GetNodeSeg(  CORNER,true, world,VALUE);        
    ////    
    for(unsigned int ino=0;ino<ns_x0.Size();ino++){
      double x0[3];  ns_x0.GetValue(ino,x0);    
      double sx[3];  ns_sx.GetValue(ino,sx);  
      double sy[3];  ns_sy.GetValue(ino,sy);      
      double x[3] = {
        aNxdxdy[iss*3+0]*x0[0] + aNxdxdy[iss*3+1]*sx[0] + aNxdxdy[iss*3+2]*sy[0],
        aNxdxdy[iss*3+0]*x0[1] + aNxdxdy[iss*3+1]*sx[1] + aNxdxdy[iss*3+2]*sy[1],
        aNxdxdy[iss*3+0]*x0[2] + aNxdxdy[iss*3+1]*sx[2] + aNxdxdy[iss*3+2]*sy[2] };
      ns_u.AddValue(ino,0,x[0]);
      ns_u.AddValue(ino,1,x[1]);
      ns_u.AddValue(ino,2,x[2]);      
    }
  }  
}

void GuessSolution_MLS
(double pre_x, double pre_y, 
 double pos_x, double pos_y,
 unsigned int id_field_disp,
 unsigned int id_field_reference,
 unsigned int id_field_lamX,  unsigned int id_field_lamY,
 const std::vector<CSolutionSensitivity>& aSolSens,
 Fem::Field::CFieldWorld& world )
{
  double mov_x = pos_x-pre_x;
  double mov_y = pos_y-pre_y;
  Fem::Field::CField& field_ref  = world.GetField(id_field_reference);  
  Fem::Field::CNodeAry::CNodeSeg& ns_ref = field_ref.GetNodeSeg( CORNER,true, world,VALUE); 
  Fem::Field::CField& field_lx  = world.GetField(id_field_lamX); 
  Fem::Field::CField& field_ly  = world.GetField(id_field_lamY);   
  const Fem::Field::CNodeAry::CNodeSeg& ns_lx = field_lx.GetNodeSeg( CORNER,true, world,VALUE);  
  const Fem::Field::CNodeAry::CNodeSeg& ns_ly = field_ly.GetNodeSeg( CORNER,true, world,VALUE);    
  unsigned int nno = ns_lx.Size();
  assert( ns_ly.Size() == nno );
  for(unsigned int ino=0;ino<nno;ino++){
    double lx[2];  ns_lx.GetValue(ino,lx);
    double ly[2];  ns_lx.GetValue(ino,ly);    
    double cr[3];  ns_ref.GetValue(ino,cr);
    ns_ref.SetValue(ino, 0, cr[0]-lx[0]*mov_x-ly[0]*mov_y);
    ns_ref.SetValue(ino, 1, cr[1]-lx[1]*mov_x-ly[1]*mov_y); 
  }  
  double nr;
  std::vector<double> aNxdxdy;
  {
    aNxdxdy.resize(aSolSens.size()*3,0);
    std::vector<double> aTmp;
    aTmp.resize(aSolSens.size()*4,0);
    double pc[3] = { 1, pos_x, pos_y };
    double pr[3] = { 1, pre_x, pre_y };
    const double sqdist_r = (pos_x-pre_x)*(pos_x-pre_x) + (pos_y-pre_y)*(pos_y-pre_y);
    double wr = 1.0/(sqdist_r+0.001);
    for(unsigned int iss=0;iss<aSolSens.size();iss++){
      if( !aSolSens[iss].is_active ) continue;
      double obj_x = aSolSens[iss].obj_x;
      double obj_y = aSolSens[iss].obj_y;
      aTmp[iss*4+0] = 1;
      aTmp[iss*4+1] = obj_x;
      aTmp[iss*4+2] = obj_y;
      const double sqdist = (pos_x-obj_x)*(pos_x-obj_x) + (pos_y-obj_y)*(pos_y-obj_y);
      aTmp[iss*4+3] = 1.0/(sqdist+0.001);
    }
//    double deriv_ratio = 1.0;
    double invG[9];
    {
      double G[9] = { 0,0,0,0,0,0,0,0,0 };
      for(unsigned int i=0;i<3;i++){
        for(unsigned int j=0;j<3;j++){
          G[i*3+j] += pr[i]*pr[j]*wr;
        }
      }      
      for(unsigned int iss=0;iss<aSolSens.size();iss++){
        if( !aSolSens[iss].is_active ) continue;
        double w0 = aTmp[iss*4+3];
        for(unsigned int i=0;i<3;i++){
          for(unsigned int j=0;j<3;j++){
            G[i*3+j] += aTmp[iss*4+i]*aTmp[iss*4+j]*w0;
          }
        }
//        G[4] += w0*deriv_ratio;
//        G[8] += w0*deriv_ratio;      
      }
      G[4] += 2.0e-53;
      G[8] += 2.0e-53; 
      for(unsigned int i=0;i<9;i++){ invG[i] = G[i]; }
      CalcInvMat3(invG);
    }    
    nr = 0;
    for(unsigned int i=0;i<3;i++){
      for(unsigned int j=0;j<3;j++){      
        nr += pc[i]*invG[i*3+j]*pr[j]*wr;
      }
    }        
    for(unsigned int iss=0;iss<aSolSens.size();iss++){
      if( !aSolSens[iss].is_active ) continue;
      const double w0 = aTmp[iss*4+3];                          
      for(unsigned int i=0;i<3;i++){
        for(unsigned int j=0;j<3;j++){      
          aNxdxdy[iss*3+0] += pc[i]*invG[i*3+j]*aTmp[iss*4+j]*w0;
        }
      }
      for(unsigned int i=0;i<3;i++){    
//        aNxdxdy[iss*3+1] += pc[i]*invG[i*3+1]*w0*deriv_ratio;
        aNxdxdy[iss*3+1] = 0;
      }
      for(unsigned int i=0;i<3;i++){    
//        aNxdxdy[iss*3+2] += pc[i]*invG[i*3+2]*w0*deriv_ratio;
        aNxdxdy[iss*3+2] = 0;
      }    
    }
  }
  //  std::cout << nr << " " << n0 << " " << n0x << " " << n0y << std::endl;
  //  getchar();
  
  {
    Fem::Field::CField& field_disp = world.GetField(id_field_disp);
    Fem::Field::CNodeAry::CNodeSeg& ns_u   = field_disp.GetNodeSeg(CORNER,true, world,VALUE);     
    Fem::Field::CNodeAry::CNodeSeg& ns_X   = field_disp.GetNodeSeg(CORNER,false,world,VALUE);        
    for(unsigned int ino=0;ino<ns_u.Size();ino++){
      double cu[3];  ns_u.GetValue( ino,cu);
      double cX[3];  ns_X.GetValue( ino,cX);    
      double lx[2];  ns_lx.GetValue(ino,lx);    
      double ly[2];  ns_ly.GetValue(ino,ly);        
      double cx[3] = {      
        cX[0]-lx[0]*mov_x-ly[0]*mov_y+cu[0],
        cX[1]-lx[1]*mov_x-ly[1]*mov_y+cu[1],
        cX[2]                        +cu[2] };
      double x[3] = {
        nr*cx[0],
        nr*cx[1],
        nr*cx[2] };
      double u[3] = {
        x[0]-cX[0],
        x[1]-cX[1],
        x[2]-cX[2] };
      ns_u.SetValue(ino,0,u[0]);
      ns_u.SetValue(ino,1,u[1]);
      ns_u.SetValue(ino,2,u[2]);
    }          
  }
  
  for(unsigned int iss=0;iss<aSolSens.size();iss++){
    if( !aSolSens[iss].is_active ) continue;
    Fem::Field::CField& field_disp = world.GetField(id_field_disp);
    Fem::Field::CNodeAry::CNodeSeg& ns_u   = field_disp.GetNodeSeg(CORNER,true, world,VALUE);     
    
    unsigned int id_field_senseX = aSolSens[iss].id_field_dudpx;
    unsigned int id_field_senseY = aSolSens[iss].id_field_dudpy;  
    unsigned int id_field_x      = aSolSens[iss].id_field_x;
    ////
    Fem::Field::CField& field_sx = world.GetField(id_field_senseX);
    const Fem::Field::CNodeAry::CNodeSeg& ns_sx  = field_sx.GetNodeSeg( CORNER,true, world,VALUE);
    Fem::Field::CField& field_sy = world.GetField(id_field_senseY);
    const Fem::Field::CNodeAry::CNodeSeg& ns_sy  = field_sy.GetNodeSeg( CORNER,true, world,VALUE);        
    Fem::Field::CField& field_x = world.GetField(id_field_x);
    const Fem::Field::CNodeAry::CNodeSeg& ns_x0  = field_x.GetNodeSeg(  CORNER,true, world,VALUE);        
    ////    
    for(unsigned int ino=0;ino<ns_x0.Size();ino++){
      double x0[3];  ns_x0.GetValue(ino,x0);    
      double sx[3];  ns_sx.GetValue(ino,sx);  
      double sy[3];  ns_sy.GetValue(ino,sy);      
      double x[3] = {
        aNxdxdy[iss*3+0]*x0[0] + aNxdxdy[iss*3+1]*sx[0] + aNxdxdy[iss*3+2]*sy[0],
        aNxdxdy[iss*3+0]*x0[1] + aNxdxdy[iss*3+1]*sx[1] + aNxdxdy[iss*3+2]*sy[1],
        aNxdxdy[iss*3+0]*x0[2] + aNxdxdy[iss*3+1]*sx[2] + aNxdxdy[iss*3+2]*sy[2] };
      ns_u.AddValue(ino,0,x[0]);
      ns_u.AddValue(ino,1,x[1]);
      ns_u.AddValue(ino,2,x[2]);      
    }
  }  
}


