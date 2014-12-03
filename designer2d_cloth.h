/*
 *  design2d_cloth.h
 *  sensitive couture
 *
 *  Created by Nobuyuki Umetani on 9/14/10.
 *  Copyright 2010 The University of Tokyo and Columbia University. All rights reserved.
 *
 */

#if !defined(DESIGNER_2D_CLOTH_H)
#define DESIGNER_2D_CLOTH_H

#if defined(__APPLE__) && defined(__MACH__)
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include <set>
#include <map>
#include <vector>

#include "delfem/cad_obj2d_move.h"
#include "delfem/mesher2d_edit.h"
#include "delfem/drawer_msh.h"

#include "analysis2d_cloth_static.h"

//! class for combining simulation and design
class CDesigner2D_Cloth
{
public:
	//! constructor
	CDesigner2D_Cloth(){
    tex_scale = 1;
    face_color[0] = 1;
    face_color[1] = 1;
    face_color[2] = 1;    
		is_updated_coord = false;
		is_updated_edge = false;
		m_is_solve_cad_change = true;
		pDrawerCAD = 0;
    pDrawerMsh = 0;
    pAnalysis = 0;
	}
	//! destructor
  virtual ~CDesigner2D_Cloth(){
		if( pDrawerCAD != 0 ) delete pDrawerCAD;
		if( pDrawerMsh != 0 ) delete pDrawerMsh;
	}
	////////////////
	// public virtual function
	void SetAnalysisInitialize(CAnalysis2D_Cloth_Static* pAnalysis, unsigned int iprob);
	void SetAnalysis(CAnalysis2D_Cloth_Static* pAnalysis);
	void Draw(unsigned int idaw_obj = 0);
	void DrawSelection();
	void Solve_fromCad_InitValue();	//!< function to solve from shape
  void Solve_fromCad_InterpValue(){
    std::vector< std::pair<unsigned int,unsigned int> > aNewL;
    return this->Solve_fromCad_InterpValue(aNewL);
  }
  void Solve_fromCad_InterpValue(const std::vector< std::pair<unsigned int,unsigned int> >& aNewL);
	void Solve_ifNeeded();	//!< function to solve if needed
	void Serialize( Com::CSerializer& arch );	//!< save/load to file
	//! set if solve when shape change
	void Enable_SolveCadChange(bool flg){ m_is_solve_cad_change = flg; }
  
  // Slider
  unsigned int GetNumberOfSlider() const { return slider_deform.count(); }
  std::string GetNameSlider(unsigned int islider){
    std::string name;
    slider_deform.GetSliderProperty(islider,name);
    return name;
  }
  double GetValueSlider(unsigned int islider, double& min, double& max) const { 
    return slider_deform.GetValueSlider(islider, min,max); 
  }
  void SetValueSlider(unsigned int islider, double val);
  
	////////////////////////////////
	// gui related functions

	// bounding box of shape
	// rot : 3x3 rotation matrix
	Com::CBoundingBox3D GetBoundingBox(double rot[]) const;

	// make mesh follow CAD (please call in idle loop)
	void FollowMshToCad_ifNeeded();

	////////////////////////////////
	// pick related functions

	void SetSelection(const std::vector<Com::View::SSelectedObject>& aSelecObj);
  void Cad_SetPicked(Cad::CAD_ELEM_TYPE itype, unsigned int id, double x, double y){
    m_picked_x = x; m_picked_y = y;
    m_id_cad_part = id; m_itype_cad_part = itype;
  }
	void Cad_GetPicked(Cad::CAD_ELEM_TYPE& itype, unsigned int& id, double& x, double& y) const ;
  
	////////////////////////////////
	// cad related functions

  void HilightCadTypeID(Cad::CAD_ELEM_TYPE itype, unsigned int id){
    pDrawerCAD->ClearSelected();    
    if( !this->cad_2d.IsElemID(itype,id) ) return;
    pDrawerCAD->AddSelected(itype,id);
	}
  void SetTextureScale_CadFace(double tex_scale);
  void SetTextureCenter_FaceCAD(double cent_x,double cent_y);
  void SetColor_CadFace(double r, double g, double b);

	////////////////////////////////
	// mesh related function
  
  void Msh_SetMeshingMode_ElemSize(unsigned int esize){
    if( esize == 0 ) return;
    mesh_2d.SetMeshingMode_ElemSize(esize);
    this->Solve_fromCad_InitValue();
  }
  void Msh_SetMeshingMode_Length(double elen){
    mesh_2d.SetMeshingMode_ElemLength(elen);
    this->Solve_fromCad_InitValue();
  }
	void Msh_PrecompDrag(Cad::CAD_ELEM_TYPE itype_cad_elem, unsigned id_cad_elem);
  bool IsntMeshDeformSensitivity() const { return isnt_mesh_deform_sensitivity; } 
  void Msh_PrecompSlider(unsigned int islider);
  void Msh_GetHarmonicPrecomp(std::vector<double>& har) const{
    mesh_2d.GetXYHarmonicFunction(har);
  }
	//! store CAD model to DXF
	bool File_WriteDXF(const std::string& fname, double scale) const {
		return cad_2d.WriteToFile_dxf(fname, scale);
  }
  const Cad::CCadObj2D& GetCad() const{ return cad_2d; }
  const Msh::CMesher2D& GetMesh() const{ return mesh_2d; }
public:
  bool Cad_AddHole(const Com::CVector2D& pos, const Com::CVector2D& poe);
  bool Cad_AddCutLine(const Com::CVector2D& pos, const Com::CVector2D& poe);  
  bool Cad_Remove(Cad::CAD_ELEM_TYPE itype, unsigned int id);  
  bool Cad_Move(Cad::CAD_ELEM_TYPE itype_cad, unsigned int id_part_cad,
                const Com::CVector2D& pos_obj_pre, const Com::CVector2D& pos_obj_cur, double tor = -1.0);   
	bool Cad_SetCurveType(unsigned int id_e, unsigned int itype);	//! 0:line, 1:arc, 2:polyline
  bool Cad_DragArc(unsigned int id_part_cad,const Com::CVector2D& pos_obj0);  
  //! if edge (ID:id_e) is an polyline make it go through point(dist)
  bool Cad_DragPolyline(unsigned int id_e, const Com::CVector2D& dist);  
  bool Cad_PreCompDragPolyline(unsigned int id_e, const Com::CVector2D& pick_pos); 
  // smoothing edge (id_e) if radius is negative smooth whole edge 
  bool Cad_SmoothingPolylineEdge(unsigned int id_e, unsigned int niter,
                                 const Com::CVector2D& pos, double radius);  
  void LoadTimeStamp();
  void SaveTimeStamp();
  
protected:
	void InitDrawer();
  bool MoveLoop(unsigned int id_l, const Com::CVector2D& dir);
  bool MoveEdge(unsigned int id_e, const Com::CVector2D& dir);
  bool MoveVertex(unsigned int id_v, const Com::CVector2D& dir);  
  unsigned int AddDartDiamond(unsigned int id_l,const Com::CVector2D& pos, const Com::CVector2D& poe);  
  bool AddDartDiamond_Sym(unsigned int id_l,const Com::CVector2D& pos, const Com::CVector2D& poe);    
  unsigned int AddDartEdge(unsigned int id_l, unsigned int id_e0, bool is_same_dir0, const Com::CVector2D& pos, const Com::CVector2D& p_nearest);  
  bool AddDartEdge_Sym(unsigned id_l, const Com::CVector2D& pos0, const Com::CVector2D& poe0);
protected:
  class CDartOnEdge{
  public:
    unsigned int id_l;    
    unsigned int id_vc; 
    unsigned int id_v1, id_v2;
    unsigned int id_vo;        
    unsigned int id_e1, id_e2;
    unsigned int id_l_in;
  };
  class CDartDiamond{
  public:
    unsigned int id_l;   
    unsigned int id_l_in;
    unsigned int id_vu, id_vd; 
    unsigned int id_vl, id_vr;
    unsigned int id_e_ul, id_e_ld, id_e_dr, id_e_ru;
  };  
protected:
  double face_color[3];
  double tex_scale;
  
  CSliderDeform slider_deform;
  std::vector<CDartDiamond> aDartDiamond;
  std::vector<CDartOnEdge> aDartOnEdge;
  std::vector< std::pair<unsigned int,unsigned int> > aSymIdVPair;
  /////
  // cad related functions
  Cad::CCadObj2D_Move cad_2d;
	std::set<unsigned int> setIdVCad_NeedFollow;
	Cad::View::CDrawer_Cad2D* pDrawerCAD;

	// mesh related functions
	Msh::CMesher2D_Edit mesh_2d;
	Msh::View::CDrawerMsh2D* pDrawerMsh;
	
  // picked object
	double m_picked_x, m_picked_y;
	unsigned int m_id_cad_part;
	Cad::CAD_ELEM_TYPE m_itype_cad_part;
  //
  std::vector<double> aLamXY_Vtx;
  Cad::CAD_ELEM_TYPE itype_cad_elem_prec;
  unsigned id_cad_elem_prec;  
  bool isnt_mesh_deform_sensitivity;
  
	bool m_is_solve_cad_change;
	bool is_updated_cad;
	bool is_updated_coord;	// is updated the cad element
	bool is_updated_edge;	// is updated mesh edge

  CAnalysis2D_Cloth_Static* pAnalysis;
};

#endif
