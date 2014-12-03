/*
DelFEM (Finite Element Analysis)
Copyright (C) 2009  Nobuyuki Umetani    n.umetani@gmail.com

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/*! @file
@brief Interactive design and analysis
@author Nobuyuki Umetani
*/

#if !defined(GUI_LISTENER_INTERACTIVE_ANALYSIS_2D_H)
#define GUI_LISTENER_INTERACTIVE_ANALYSIS_2D_H

#if defined(_WIN32)
#  include <windows.h>
#if defined(__VISUALC__)
#  pragma comment (lib, "winmm.lib")      /* link with Windows MultiMedia lib */
#  pragma comment (lib, "opengl32.lib")  /* link with Microsoft OpenGL lib */
#  pragma comment (lib, "glu32.lib")     /* link with Microsoft OpenGL Utility lib */
#endif
#endif  /* _WIN32 */

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

#include "cad_obj2d_move.h"
#include "mesher2d_edit.h"
#include "analysis2d_interface.h"

#include "delfem/drawer_msh.h"

//! class for combining simulation and design
class CDesigner2D_Analysis
{
public:
	//! constructor
	CDesigner2D_Analysis(){
		is_updated_coord = false;
		is_updated_edge = false;
		m_is_solve_cad_change = true;
		pDrawerCAD = 0;
    pDrawerMsh = 0;
    pAnalysis = 0;
	}
	//! destructor
  virtual ~CDesigner2D_Analysis(){
		if( pDrawerCAD != 0 ) delete pDrawerCAD;
		if( pDrawerMsh != 0 ) delete pDrawerMsh;
	}
	////////////////
	// public virtual function

	void SetAnalysisInitialize(IAnalysis2D* pAnalysis,unsigned int iprob);
	void SetAnalysis(IAnalysis2D* pAnalysis);
	void Draw(unsigned int idaw_obj = 0);
	void DrawSelection();
	void Solve_fromCad_InitValue();	//!< solve all from shape (including mesh,field,ls generation)
	void Solve_fromCad_InterpValue();	//!< solve all from shape (including mesh,field,ls generation) 
	void Solve_ifNeeded();	//!< solve if needed
	void Serialize( Com::CSerializer& arch );	//!< io function (save/load)

	//! 形状が変わったときに，方程式を解くかどうか
	void Enable_SolveCadChange(bool flg){ m_is_solve_cad_change = flg; }
  
	////////////////////////////////
	// GUIに関するクラス

	// get bounding box of the shape
	// rot : 3x3 matrix
	Com::CBoundingBox3D GetBoundingBox(double rot[]) const;

	// CADに追従していないメッシュを追従させる(Idle時に実行する)
	void FollowMshToCad_ifNeeded();

	////////////////////////////////
	// pick related function

	void SetSelection(const std::vector<Com::View::SSelectedObject>& aSelecObj);
	void Cad_GetPicked(Cad::CAD_ELEM_TYPE& itype, unsigned int& id, double& x, double& y) const ;

	////////////////////////////////
	// CADに関する関数

  void HilightCadTypeID(Cad::CAD_ELEM_TYPE itype, unsigned int id){
    if( !this->cad_2d.IsElemID(itype,id) ) return;
    pDrawerCAD->ClearSelected();
    pDrawerCAD->AddSelected(itype,id);
	}

	////////////////////////////////
	// Mshに関する関数

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
  void Msh_GetHarmonicPrecomp(std::vector<double>& har) const{
    mesh_2d.GetXYHarmonicFunction(har);
  }

	// ループを穴に設定する
  void SetLoopType(LOOP_TYPE loop_type, unsigned int id_l)
  {
		if( !cad_2d.IsElemID(Cad::LOOP,id_l) ) return;
    if(      loop_type == LOOP_MESH  ){
      if( mesh_2d.IsIdLCad_CutMesh(id_l)  && !pAnalysis->IsRigid(id_l) ) return;
      mesh_2d.AddIdLCad_CutMesh(id_l);
      pAnalysis->SetRigid(false,id_l);
		}
    else if( loop_type == LOOP_RIGID ){	
			mesh_2d.RemoveIdLCad_CutMesh(id_l);
			pAnalysis->SetRigid(true,id_l);
		}
    else if( loop_type == LOOP_HOLE  ){
			mesh_2d.RemoveIdLCad_CutMesh(id_l);
			pAnalysis->SetRigid(false,id_l);
    }
		this->Solve_fromCad_InitValue();
	}
  LOOP_TYPE GetLoopType(unsigned int id_l){
    if( !cad_2d.IsElemID(Cad::LOOP,id_l) ) return NOT_LOOP;
		if( mesh_2d.IsIdLCad_CutMesh(id_l) ) return LOOP_MESH;
    if( pAnalysis == 0 ){ return LOOP_HOLE; }
    if( pAnalysis->IsRigid(id_l) ){ return LOOP_RIGID; }
    return LOOP_HOLE;
  }
	//! save shape to dxf file
	void File_WriteDXF(const std::string& fname) const {
		cad_2d.WriteToFile_dxf(fname,1);
  }
  const Cad::CCadObj2D& GetCad() const{ return cad_2d; }
  const Msh::CMesher2D& GetMesh() const{ return mesh_2d; }
public:
	unsigned int Cad_MakeRectRegion(unsigned int id_l, const Com::CVector2D& pos_obj0, const Com::CVector2D& pos_obj1 );
  unsigned int Cad_ConnectVertex_Line(unsigned int id_v0, unsigned int id_v1);
  unsigned int Cad_ConnectVertex(const Cad::CEdge2D& e);
  bool Cad_Remove(Cad::CAD_ELEM_TYPE itype, unsigned int id);
  unsigned int Cad_AddVertex(Cad::CAD_ELEM_TYPE itype, unsigned int id, const Com::CVector2D& pos_obj0);
  void Cad_SetColor_Loop(unsigned int id_l, double color[3]);	
	bool Cad_SetCurveType(unsigned int id_e, unsigned int itype);	//! 既存のカーブに出来るだけ沿うように新しいカーブを設定する
  bool Cad_SetCurve_Polyline(unsigned int id_e, const std::vector<Com::CVector2D>& aVec);
  bool Cad_DragArc(unsigned int id_part_cad,const Com::CVector2D& pos_obj0);
  bool Cad_Move(Cad::CAD_ELEM_TYPE itype_cad, unsigned int id_part_cad,
                const Com::CVector2D& pos_obj_pre, const Com::CVector2D& pos_obj_cur);
  //! id_eがメッシュ辺なら，スムージングをかけて滑らかにする(radiusが負なら全体をスムージングする)
  bool Cad_SmoothingPolylineEdge(unsigned int id_e, unsigned int niter, const Com::CVector2D& pos, double radius);
protected:
	void CreatedLoop(unsigned int id_l_created, unsigned int id_l_from){}
	void InitDrawer();
protected:
  // Mesh deformation precomputation
  Cad::CAD_ELEM_TYPE itype_cad_elem_prec;
  unsigned int id_cad_elem_prec;
  std::vector<double> aLamXY_Vtx;  
  
  // cad related function
  Cad::CCadObj2D_Move cad_2d;
	std::set<unsigned int> setIdVCad_NeedFollow;
	Cad::View::CDrawer_Cad2D* pDrawerCAD;

	// mesh related fuction
	Msh::CMesher2D_Edit mesh_2d;
	Msh::View::CDrawerMsh2D* pDrawerMsh;
	
  // picked element
	double m_picked_x, m_picked_y;
	unsigned int m_id_cad_part;
	Cad::CAD_ELEM_TYPE m_itype_cad_part;

	bool m_is_solve_cad_change;
	bool is_updated_cad;
	bool is_updated_coord;	// 要素の座標がUpdateされたかどうか
	bool is_updated_edge;	// 要素の辺がUpdateされたかどうか

  IAnalysis2D* pAnalysis;
};

#endif
