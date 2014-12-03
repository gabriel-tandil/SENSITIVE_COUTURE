/*
 *  cloth_utility.h
 *  sensitive couture
 *
 *  Created by Nobuyuki Umetani on 9/17/10.
 *  Copyright 2010 The University of Tokyo and Columbia University. All rights reserved.
 *
 */

#if !defined(CLOTH_UTILITY_H)
#define CLOTH_UTILITY_H

#include "delfem/field_world.h"

#include "contact_target.h"
#include "eqn_contact3d.h"

#include "stitch_array.h"

class CClothParam{
public:
  CClothParam(){
    stiff_bend   = 0;
    stiff_myu    = 1;
    stiff_lambda = 0;
    rho          = 1;
  }
  CClothParam(const CClothParam& cp){
    this->stiff_bend   = cp.stiff_bend;
    this->stiff_myu    = cp.stiff_myu;
    this->stiff_lambda = cp.stiff_lambda;
    this->rho          = cp.rho;
  }
public:
  double stiff_bend;
  double stiff_myu;
  double stiff_lambda;
  double rho;      
};

class CContactParam{    
public:    
  double myu_k;
  double myu_s;
  double stiff_n;
  double stiff_f;
  double offset;
};      

bool KineticDamping(double ke0, double ke1, double ke2,
                    unsigned int id_field_disp, Fem::Field::CFieldWorld& world);

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
 ////
 unsigned int id_field_disp, 
 unsigned int id_field_hinge, 
 unsigned int id_field_disp_buffer,
 Fem::Field::CFieldWorld& world);

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
 /////
 unsigned int id_field_disp, 
 unsigned int id_field_hinge,  
 ////
 bool is_xy,
 unsigned int id_field_senseX, unsigned int id_field_senseY,
 unsigned int id_field_lamX, unsigned int id_field_lamY,
 Fem::Field::CFieldWorld& world);


void DrawSeamLine(unsigned int id_field_disp, unsigned int id_field_dart, const Fem::Field::CFieldWorld& world);

class CInterpBarycentric
{
public:
  CInterpBarycentric(){
    ino1 = 0;
    ino2 = 0;
    ino3 = 0;
    r1_ini = 0;
    r2_ini = 0;    
    height = 0;
    r1 = 0;
    r2 = 0;
  }
public:
  unsigned int ino1, ino2, ino3;
  double r1_ini;
  double r2_ini;
  ////
  double height;
  double r1;
  double r2;
};

// copy id1 into id0
void CopyValueVelo(unsigned int id0, unsigned int id1, Fem::Field::CFieldWorld& world);

// copy id1 into id0
void SetDeformedValue(unsigned int id0, unsigned int id1, Fem::Field::CFieldWorld& world);

void InterpField(unsigned int id_base_to,   unsigned int id_field_to,
                 unsigned int id_base_from, unsigned int id_field_from,
                 Fem::Field::CFieldWorld& world);

void FindBaseInterp(unsigned int id_base_to,
                    unsigned int id_base_from,
                    Fem::Field::CFieldWorld& world,
                    std::vector<CInterpBarycentric>& aInterp);

void MoveFineBaseCoord(unsigned int id_base_to,
                       unsigned int id_base_from,
                       Fem::Field::CFieldWorld& world,
                       const std::vector<CInterpBarycentric>& aInterp);

void MoveFineDeformedCoord(unsigned int id_field_disp_to,   unsigned int id_base_to,
                           unsigned int id_field_disp_from, unsigned int id_base_from,
                           Fem::Field::CFieldWorld& world,
                           const std::vector<CInterpBarycentric>& aInterp);

void UpdateFineDeformedInterp(unsigned int id_field_disp_to,   unsigned int id_base_to,
                            unsigned int id_field_disp_from, unsigned int id_base_from,
                            const Fem::Field::CFieldWorld& world,
                            std::vector<CInterpBarycentric>& aInterp);

void InitFineDeformInterp(std::vector<CInterpBarycentric>& aInterp);

//////

void NoResponse
(double mov_x, double mov_y,
 unsigned int id_field_disp,
 unsigned int id_field_reference,
 unsigned int id_field_lamX,  unsigned int id_field_lamY,
 Fem::Field::CFieldWorld& world );

void NoResponse_Slider
(double mov_v,
 unsigned int id_field_disp,
 unsigned int id_field_reference,
 unsigned int id_field_lamX,
 Fem::Field::CFieldWorld& world );


void SensitiveResponse
(double mov_x, double mov_y,
 unsigned int id_field_disp,
 unsigned int id_field_reference,
 unsigned int id_field_lamX,  unsigned int id_field_lamY,
 unsigned int id_field_senseX, unsigned int id_field_senseY,
 Fem::Field::CFieldWorld& world );

void SensitiveResponse_Slider
(double mov_v,
 unsigned int id_field_disp,
 unsigned int id_field_reference,
 unsigned int id_field_lamX,
 unsigned int id_field_senseX,
 Fem::Field::CFieldWorld& world );


class CSolutionSensitivity{
public:
  CSolutionSensitivity(){
    is_active = false;
    is_xy = true;
    obj_x = 0;
    obj_y = 0;
    val_slider = 0;
    id_field_x = 0;
    id_field_dudpx = 0;
    id_field_dudpy = 0;
  }
public:
  bool is_active;
  bool is_xy;
  double obj_x;
  double obj_y;
  double val_slider;
  unsigned int id_field_x;
  unsigned int id_field_dudpx;
  unsigned int id_field_dudpy;  
};

void GuessSolution_GMLS
(double pre_x, double pre_y, 
 double pos_x, double pos_y,
 unsigned int id_field_disp,
 unsigned int id_field_reference,
 unsigned int id_field_lamX,  unsigned int id_field_lamY,
 const std::vector<CSolutionSensitivity>& aSolSens,
 Fem::Field::CFieldWorld& world );

void GuessSolution_GMLS_Slider
(double pre_v, double pos_v,
 unsigned int id_field_disp,
 unsigned int id_field_reference,
 unsigned int id_field_lamX,
 const std::vector<CSolutionSensitivity>& aSolSens,
 Fem::Field::CFieldWorld& world );


void GuessSolution_MLS
(double pre_x, double pre_y, 
 double pos_x, double pos_y,
 unsigned int id_field_disp,
 unsigned int id_field_reference,
 unsigned int id_field_lamX,  unsigned int id_field_lamY,
 const std::vector<CSolutionSensitivity>& aSolSens,
 Fem::Field::CFieldWorld& world );


#endif



