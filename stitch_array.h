/*
 *  stitch_array.h
 *  sensitive couture
 *
 *  Created by Nobuyuki Umetani on 12/20/10.
 *  Copyright 2010 The University of Toky and Columbia University. All rights reserved.
 *
 */

#if !defined(STITCH_ARRAY_H)
#define STITCH_ARRAY_H

#include <vector>

#include "delfem/femls/linearsystem_field.h"
#include "delfem/ls/preconditioner.h"
#include "delfem/field_world.h"
#include "delfem/cad_obj2d_move.h"
#include "delfem/mesher2d_edit.h"


class CStitchAry{
public:
  CStitchAry(){
    stitch_stiff = 1000;
    stitch_damp_coeff = 1.0;
  }
  void Clear(){ aStitch.clear(); }  
  void SetStiff(double val){ this->stitch_stiff = val; }
  void SetDampingCoeff(double val){ this->stitch_damp_coeff = val; }
  void AddStitch(const Cad::CCadObj2D& cad_2d, const Msh::CMesher2D& mesh_2d, unsigned int id_e1, unsigned int id_e2);
  void CopyCADSituation(const CStitchAry& sa);
  void AddLinSys_BackwardEular(double dt, double& se,                                
                               Fem::Ls::CLinearSystem_Field& ls,
                               unsigned int id_field_disp,
                               const Fem::Field::CFieldWorld& world) const;
  void AddLinSys_Sensitivity(Fem::Ls::CLinearSystem_Field& ls,
                             unsigned int id_field_disp,
                             const Fem::Field::CFieldWorld& world) const;  
  void MakeField(unsigned int id_field_disp, unsigned int id_base, Fem::Field::CFieldWorld& world);
  void AddStitchField(const Cad::CCadObj2D& cad_2d, const Msh::CMesher2D& mesh_2d, unsigned int id_e1, unsigned int id_e2,
                      unsigned int id_field_disp, unsigned int id_base, Fem::Field::CFieldWorld& world);
  void RemoveStitchField(unsigned int id_e1, unsigned int id_e2,
                         unsigned int id_field_disp, unsigned int id_base, Fem::Field::CFieldWorld& world);
  void AddPatternLinearSystem(Fem::Ls::CLinearSystem_Field& ls, const Fem::Field::CFieldWorld& world);
  void Draw(const unsigned int id_field_disp,const Fem::Field::CFieldWorld& world) const;
  void DrawBoundaryCondition2D(const Cad::CCadObj2D& cad_2d) const;
  void SetReorderingPreconditionerIfNeeded(LsSol::CPreconditioner_ILU& prec,
                                           unsigned int id_field_disp,          
                                           const Fem::Field::CFieldWorld& world) const;
  std::vector<unsigned int> GetIdFieldAry(){
    std::vector<unsigned int> aIdField;
    for(unsigned int ist=0;ist<aStitch.size();ist++){ aIdField.push_back( aStitch[ist].id_field ); }
    return aIdField;
  }
  bool IsSeamLine(unsigned int id_e_in, unsigned int& id_e_out, bool& is_same_dir) const;
  //private:
public:
  class CStitch{
  public:
    CStitch(unsigned int id_e1, unsigned int id_e2, bool dir ){
      this->id_e1 = id_e1;
      this->id_e2 = id_e2;
      this->is_same_dir = dir;
    }
  public:
    unsigned int id_e1;
    unsigned int id_e2;
    bool is_same_dir;  
    /////
    unsigned int id_ea1;
    unsigned int id_ea2;
    /////
    unsigned int id_field;
    bool need_reordering;    
  };
  //private:
public:
  std::vector<CStitch> aStitch;
  double stitch_stiff;
  double stitch_damp_coeff;
};

#endif
