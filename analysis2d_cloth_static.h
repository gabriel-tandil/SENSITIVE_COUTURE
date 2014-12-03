/*
 *  analysis2d_cloth_static.h
 *  sensitive couture
 *
 *  Created by Nobuyuki Umetani on 9/14/10.
 *  Copyright 2010 The University of Tokyo and Columbia University. All rights reserved.
 *
 */


#if !defined(ANALYSIS_2D_CLOTH_STATIC_H)
#define ANALYSIS_2D_CLOTH_STATIC_H

#include "delfem/drawer_field.h"
#include "delfem/femls/linearsystem_field.h"
#include "delfem/ls/preconditioner.h"
#include "delfem/field_world.h"
#include "delfem/msh/surface_mesh_reader.h"
#include "delfem/drawer_cad.h"

#include "eqn_contact3d.h"
#include "cloth_handler.h"
#include "cloth_utility.h"
#include "slider_deform.h"
#include "contact_target.h"


enum SOLVER_FLAG{ 
	SUCCESS = 0,
	NOTCONV_LSSOL = 1,
	NOTCONV_NR = 2,
	ILU_FAIL = 3,
};

enum CLOTH_DESIGN_SYSTEM_MODE
{
  CLOTH_INITIAL_LOCATION, // 0
  SIMULATION_STATIC,      // 1 
  WAITING_PICK,           // 2
  SENSITIVITY_DONE,       // 3
  SIMULATION_DETAIL,      // 4
};

namespace Fem{
  namespace Field{
    namespace View{
      class CDrawerFace;
    }    
  }
}

class CAnalysis2D_Cloth_Static_TimeStamp{
public:
  Cad::CCadObj2D cad;
  Msh::CMesher2D msh;
  ////
  Fem::Field::CFieldWorld world;  
  
  // basic
  unsigned int id_field_base_;
  unsigned int id_field_disp_;
  unsigned int id_field_disp_buffer_;  
  unsigned int id_field_hinge_;  
  CStitchAry stitch_ary_;    // stitch parameters  
  std::vector<CFrictionPoint> aFrictionPoint;  
  
  // detail
  std::vector<CInterpBarycentric> aInterp_detail;  
  unsigned int id_field_base_detail_;
  unsigned int id_field_disp_detail_;   
  unsigned int id_field_disp_buffer_detail;  
  unsigned int id_field_hinge_detail;  
  CStitchAry stitch_ary_detail;  
  std::vector<CFrictionPoint> aFrictionPoint_detail;    

  // sensitivity analysis
  unsigned int id_field_lamX, id_field_lamY;          
  std::vector<CSolutionSensitivity> aSolSens;
  ////
  double cur_time;
  double dt;  
  CClothParam cloth_param;  // cloth parameters
  CContactParam contact_param;  // contact parameters
  bool is_detail;
};

//! interactive integration of shell and modeling
class CAnalysis2D_Cloth_Static
{
public:
public:
  CAnalysis2D_Cloth_Static();
  ~CAnalysis2D_Cloth_Static(){ delete pCT; }
  
  //// virtual 
  virtual void SetModelProblem(Cad::CCadObj2D& cad_2d, Msh::CMesher2D& mesh_2d, unsigned int iprob){}
  virtual void Draw() const;
  virtual void DrawBoundaryCondition(const Cad::CCadObj2D& cad_2d) const;
  virtual void UpdateAnimationDrawer(){}
  
  // Called from CDesigner
  virtual SOLVER_FLAG UpdateMeshAndSolve(const Msh::CMesher2D& mesh_2d, 
                                         bool is_upd_coord, bool is_upd_edge, 
                                         bool is_solve_mesh_change, bool is_update_disp); //!< solve fem
  void ConnectEdge(const Cad::CCadObj2D& cad, const Msh::CMesher2D& mesh_2d, unsigned int id_e1, unsigned int id_e2);
  void DisconnectEdge(unsigned int id_e1,unsigned int id_e2);  
  
  // rebuild every FEM rerated things (e.g. Field LS Drawer), and solve from initial value
  virtual void BuildFEM_ClearValueField( const Cad::CCadObj2D& cad_2d, const Msh::CMesher2D& mesh_2d);
  virtual void BuildFEM_InterpValueField(const Cad::CCadObj2D& cad_2d, const Msh::CMesher2D& mesh_2d, 
                                         const std::vector< std::pair<unsigned int,unsigned int> >& aNewL);
  
  void Update_Boundary_Condition_Cad_Move(Cad::CAD_ELEM_TYPE itype, unsigned int cad_elem_id,
                                          double del_x, double del_y, 
                                          CSliderDeform& slider_deform);
  
  virtual void Update_Boundary_Condition_Cad_Move(Cad::CAD_ELEM_TYPE itype, unsigned int cad_elem_id,
                                                  double del_x, double del_y){}  
  
  virtual void Serialize( Com::CSerializer& arch );  
  void MeshQualityInfo(double& max_aspect, bool& is_inverted) const;

  //// non-virtual
  void SetModelProblem_Cloth(Cad::CCadObj2D_Move& cad_2d, Msh::CMesher2D& mesh_2d, unsigned int iprob,                               
                             CSliderDeform& slider_deform,
                             std::vector< std::pair<unsigned int,unsigned int> >& aSymIdVPair);
  Com::CBoundingBox3D GetBoundingBox(double rot[]) const;
  void PerformStaticSolver();
  void SetClothPiecePlacingMode();  
  CLOTH_DESIGN_SYSTEM_MODE GetMode(){ return imode_; }
  void SetColor_FaceFEM(double r, double g, double b){
    face_color[0] = r;
    face_color[1] = g;
    face_color[2] = b;    
    InitDrawer();
  }
  void MoveClothLoopInitialPosition(unsigned int id_l, double der[3]);
  void SetHilight(Cad::CAD_ELEM_TYPE itype, unsigned int cad_elem_id);
  void SetTextureCenter(double cent_x, double cent_y);    
  void SetTextureScale_FaceFEM(double scale);
  bool IsSeamLine(unsigned int id_e_in, unsigned int& id_e_out, bool& is_same_dir) const;
  
  void SetSensitivity(Cad::CAD_ELEM_TYPE itype_cad_part, unsigned int id_cad_part,
                      const std::vector<double>& har,
                      double objx, double objy);  

  void GuessSolutionMouseMove(double pre_x, double pre_y,
                              double pos_x, double pos_y);


  // Time Stamp
  void SaveTimeStamp(const Cad::CCadObj2D& cad_2d, const Msh::CMesher2D& mesh_2d);
  void LoadTimeStamp(Cad::CCadObj2D& cad_2d, Msh::CMesher2D& mesh_2d);
  bool IsBlowUp() const {
    return false;
    if(total_energy_ref < 0 ){ return false; }
    if(total_energy_ >= total_energy_ref * 1.0e+10 ) return true; 
    return false;
  }

  //Slider
  void SetSensitivity_Slider(unsigned int islider,
                             const std::vector<double>& har,
                             double val);    
  void GetSensitivityElem_Slider(const unsigned int no[3], const double r[3], 
                                 double sns[3], double sns2[3]) const;
  void GuessSolutionSliderMove(double pre_v, double pos_v);  
    
  ///////////
  // Set Parameter
  
  void SetParam_Cloth(const CClothParam& p){ this->cloth_param = p; }
  void SetParam_Contact(const CContactParam& p){ this->contact_param = p; }  
  void SetTimeStep(double dt){ this->dt_ = dt; }
  void SetIsDetail(bool is_detail,const Cad::CCadObj2D& cad_2d, const Msh::CMesher2D& mesh_2d);
  void SetIsShowEdge(bool is_show);
  void SetIsLighting(bool is_lighting);
  void SetIsDrawPatternBoundary(bool is_draw){ this->is_draw_pattern_boundary = is_draw; }

  ///////////
  // Get Parameter
        
  CClothParam GetParam_Cloth() const { return this->cloth_param; }
  CContactParam GetParam_Contact() const { return this->contact_param; }  
  double GetCurTime() const { return this->cur_time_; } 
  double GetStampedTime() const { return this->time_stamp.cur_time; }
  double GetTimeStep() const { return this->dt_; }  
  bool IsDetail() const { return is_detail_; }  
  bool IsShowEdge() const { return is_show_edge_!= 0; }  
  bool IsLighting() const { return this->is_lighting_; }  
  
  bool Pick(double scrx, double scry,
            const double trans0[], const double rot[], const double trans1[],
            unsigned int picked_elem_nodes[3], double picked_elem_ratio[3],
            unsigned int& id_l);
  
  bool WriteObjMeshSTL(const std::string& fname,double scale){ return obj_mesh.WriteSTL(fname,scale); }
  virtual SOLVER_FLAG Solve();  
private:
  void ClearDetailField();
  void MakeDetailField(const Cad::CCadObj2D& cad_2d, const Msh::CMesher2D& mesh_2d);
  virtual void ClearLinearSystemPreconditioner();
  virtual void InitFieldEqn_fromMsh(const Cad::CCadObj2D& cad_2d, const Msh::CMesher2D& mesh_2d);
  virtual void InitDrawer();
  virtual void CreatedLoop(unsigned int id_l_created, unsigned int id_l_from){}
private:    
  double face_color[3];
  double tex_scale;
  ///
  double total_energy_;
  double total_energy_ref;
  double gx_,gy_,gz_;
  unsigned int is_show_edge_;  
  bool is_lighting_;
  bool is_draw_pattern_boundary;
  unsigned int imode_sensitivity_guess_;  // (0)noguess, (1)sens, (2)sens+upd, (3)gmls, (4)prg+gmls, (5)prg+mls
  Cad::CAD_ELEM_TYPE itype_hilight_;
  unsigned int id_hilight_;    
  double mouse_obj_x, mouse_obj_y;
  double mouse_slider_v;
  
  CContactTarget3D* pCT;    
  CClothHandler clothHandler_;      
  CSurfaceMeshReader obj_mesh;
  
  // parameters
  double cur_time_;    // current time
  double dt_;
  CClothParam cloth_param;  // cloth parameters
  CContactParam contact_param;  // contact parameters
  CLOTH_DESIGN_SYSTEM_MODE imode_;  // type of mode      
  
  Fem::Field::CFieldWorld world;  
  
  // sensitivity analysis
  unsigned int id_field_lamX, id_field_lamY;      
  std::vector<CSolutionSensitivity> aSolSens;
  
  // coarse simulation
  unsigned int id_field_base; // The handle of mesh imported to FEM
  unsigned int id_field_disp;
  unsigned int id_field_hinge;
  unsigned int id_field_disp_buffer;
  CStitchAry stitch_ary_;    // stitch parameters
  std::vector<CFrictionPoint> aFrictionPoint;  
  
  // detail analysis
  bool is_detail_;    
  std::vector<CInterpBarycentric> aInterp_detail;
  unsigned int id_field_base_detail;
  unsigned int id_field_disp_detail;
  unsigned int id_field_disp_buffer_detail;  
  unsigned int id_field_hinge_detail;
  CStitchAry stitch_ary_detail;    
  std::vector<CFrictionPoint> aFrictionPoint_detail;  
  
  CAnalysis2D_Cloth_Static_TimeStamp time_stamp;  
  
  ////
  // linear system setting
  Fem::Ls::CLinearSystem_Field ls;
  LsSol::CPreconditioner_ILU prec;
  bool is_cleared_ls_prec;
  
  std::vector<unsigned int> aIdECad_Fix;
  std::vector<unsigned int> aIdField_Fix;
  std::vector<unsigned int> aIdField_Fix_detail;  
   
  Fem::Field::View::CDrawerArrayField m_aDrawerField;
  Fem::Field::View::CDrawerFace *pDF;  

  Fem::Field::View::CDrawerArrayField m_aDrawerField_detail;  
  Fem::Field::View::CDrawerFace *pDF_detail; 
};

#endif
