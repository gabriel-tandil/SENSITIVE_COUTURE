/*
 *  contact_target_cad3d.cpp
 *  cad_view
 *
 *  Created by Nobuyuki Umetani on 3/5/11.
 *  Copyright 2011 The University of Tokyo. All rights reserved.
 *
 */

#include "../include/contact_target_cad3d.h"



CContactTarget3D_Cad3D::CContactTarget3D_Cad3D(const Cad::CCadObj3D& cad3d)
{  
  SetCad3D(cad3d);
}

void CContactTarget3D_Cad3D::SetCad3D(const Cad::CCadObj3D& cad3d)
{
  aLoop.clear();
  const std::vector<unsigned int>& aIdL = cad3d.GetAryElemID(Cad::LOOP);
  for(unsigned int iil=0;iil<aIdL.size();iil++){
    const unsigned int id_l = aIdL[iil];    
    const Cad::CLoop3D& l = cad3d.GetLoop(id_l);
    aLoop.push_back(l);
  }        
}

double CContactTarget3D_Cad3D::Projection
(double px, double py, double pz,
 double n[3]) const
{
  const Com::CVector3D p(px,py,pz);
  bool is_inside = true;
  {
    for(unsigned int i=0;i<3;i++){
      Com::CVector3D dir;
      if(      i == 0 ){ dir = Com::CVector3D(1,0,0); }
      else if( i == 1 ){ dir = Com::CVector3D(0,1,0); }
      else if( i == 2 ){ dir = Com::CVector3D(0,1,0); }      
      unsigned int cross_counter = 0;      
      bool is_amb = false;
      for(unsigned int iil=0;iil<aLoop.size();iil++){
        const Cad::CLoop3D& l = aLoop[iil];      
        int num_sect = l.NumIntersecRay(p,dir);
        if( num_sect == -1 ){ is_amb = true; break; }
//        if( num_sect !=0 ){ std::cout << num_sect << std::endl; }
        cross_counter += num_sect;
      } 
      if( !is_amb ){
        is_inside = (cross_counter%2==1);
        break;
      }
    }
  }
  double min_dist = -1;
  Com::CVector3D p_near;
  for(unsigned int iil=0;iil<aLoop.size();iil++){
    const Cad::CLoop3D& l = aLoop[iil];
    const Com::CVector3D& pn = l.GetNearestPoint(p);
    const double d0 = Distance(pn,p);
    if( min_dist > 0 && d0 > min_dist ){ continue; }
    min_dist = d0; 
    p_near = pn;
  }
  Com::CVector3D dir(p-p_near); dir.Normalize(); 
  n[0] = dir.x;
  n[1] = dir.y;
  n[2] = dir.z;  
  double depth = Distance(p, p_near);
  if( !is_inside ){ depth *= -1; }
  return depth;
}


Com::CBoundingBox3D CContactTarget3D_Cad3D::GetBoundingBox() const
{
  Com::CBoundingBox3D bb;
  for(unsigned int iil=0;iil<aLoop.size();iil++){
    const Cad::CLoop3D& l = aLoop[iil];
    bb += l.GetBoundingBox();
  }
  return bb;
}
