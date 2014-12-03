/*
 *  cloth_handler.h
 *  sensitive couture
 *
 *  Created by Nobuyuki Umetani on 9/1/10.
 *  Copyright 2010 The University of Tokyo and Columbia University. All rights reserved.
 *
 */

#if !defined(CLOTH_HANDLER_H)
#define CLOTH_HANDLER_H

#include "delfem/tri_ary_topology.h"


class CClothHandler
{    
  class CNode;
  class CClothPiece;
public:
  CClothHandler(){
    this->id_ea_hilight = 0;
    pXYZs_ = 0;
    aTri_ = 0;
    aNorm_ = 0;
    nnode_ = 0;
    ntri_ = 0;    
  }
  ~CClothHandler(){
    Clear();
  }

  void Clear();
  bool MoveAnchor_2D(double x, double y, unsigned int id_ea);
  bool GetAnchor_2D(double r[2], unsigned int id_ea) const;
  bool GetAnchor_2D_Loop(double r[2], unsigned int id_l) const;  
  bool GetAnchor_3D(double p[3], double n[3], double h[3], unsigned int id_ea) const; 
  void BuildClothMeshTopology(unsigned int id_base, unsigned int id_field_disp, const Fem::Field::CFieldWorld& world);
  ////
  void Draw(unsigned int imode) const;
  ////
  void SetObjectMesh(const std::vector<unsigned int>& aTri,
                     const std::vector<double>& aXYZ);                  
  /////
  void AddClothPiece(unsigned int id_l, double cent_x, double cent_y);
  void AddClothPiece(unsigned int id_l_new, unsigned int id_l_old);  
  void Transform_Cloth_Pan(unsigned int id_l, double anc_x, double anc_y, double anc_z);  
  void Transform_Cloth_RotBryantAngle(unsigned int id_l, double phi, double theta, double psi); 
  void SetRadius(unsigned int id_l, double r);
    /////
  bool SetClothLocation(unsigned int id_field_disp, Fem::Field::CFieldWorld& world);
  bool Pick(double scrx, double scry,
            const double trans0[3], const double rot[3], const double trans1[3],
            const double dir[3], const double org[3],            
            unsigned int id_field_disp, const Fem::Field::CFieldWorld& world );
private:
  unsigned int id_ea_hilight; // hilighted elemary
  unsigned int itype; // itype operation : 0:pan, 1:rot
  double hit_pos[3];
  std::vector<CClothPiece*> apPiece;
  
  /////
  unsigned int nnode_;
  double* pXYZs_;   // nno*3
  double* aNorm_;   // nno*3
  unsigned int ntri_;
  unsigned int* aTri_;  
  
private:
  class CClothPiece
  {
  public:
    // consistient anytime
    unsigned int id_l;
    double cent[2]; // center of 2D cloth piece
    double p[3];    
    double n[3];
    double h[3];
    /// update when field change
    CTriAryTopology topo;
    unsigned int id_ea;
    double radius;
  };
};

#endif

