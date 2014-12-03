/*
 *  contact_target_cad3d.h
 *  cad_view
 *
 *  Created by Nobuyuki Umetani on 3/5/11.
 *  Copyright 2011 The University of Tokyo. All rights reserved.
 *
 */

#include "cad_obj3d.h"
#include "contact_target.h"

class CContactTarget3D_Cad3D : public CContactTarget3D
{
public:
  CContactTarget3D_Cad3D(){}
  CContactTarget3D_Cad3D(const Cad::CCadObj3D& cad3d);
  void SetCad3D(const Cad::CCadObj3D& cad3d);
public:
	virtual void Draw() const{}
	virtual double Projection
	(double px, double py, double pz,
	 double n[3]) const;
  
  virtual bool IntersectionPoint
  (double p[3], 
   const double org[3], const double dir[3]) const{}
  
  virtual void GetMesh(std::vector<unsigned int>& aTri,
                       std::vector<double>& aXYZ,
                       double elen) const {}
  /////
  Com::CBoundingBox3D GetBoundingBox() const;
private:
  std::vector<Cad::CLoop3D> aLoop;
};

