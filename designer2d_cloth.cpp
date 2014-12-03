/*
 *  design2d_cloth.cpp
 *  sensitive couture
 *
 *  Created by Nobuyuki Umetani on 9/14/10.
 *  Copyright 2010 The University of Tokyo and Columbia University. All rights reserved.
 *
 */


#if defined(__VISUALC__)
#pragma warning ( disable : 4996 )
#pragma warning ( disable : 4786 )
#endif

#include "delfem/serialize.h"
#include "designer2d_cloth.h"

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void CDesigner2D_Cloth::DrawSelection()
{
	pDrawerCAD->DrawSelection(0);
}

Com::CBoundingBox3D CDesigner2D_Cloth::GetBoundingBox(double rot[]) const
{
  //return pDrawerCAD->GetBoundingBox(rot);
//  cad_2d.Get
  const std::vector<unsigned int>& aIdE =cad_2d.GetAryElemID(Cad::EDGE);
  Com::CBoundingBox2D bb2d;
  for(unsigned int iide=0;iide<aIdE.size();iide++){
    unsigned int id_e = aIdE[iide];
    const Cad::CEdge2D& e = cad_2d.GetEdge(id_e);
    bb2d += e.GetBoundingBox();
  }
  return Com::CBoundingBox3D(bb2d.x_min,bb2d.x_max, bb2d.y_min,bb2d.y_max,-1,1);  
}

void CDesigner2D_Cloth::SetTextureScale_CadFace(double tex_scale)
{
  this->tex_scale = tex_scale;
  this->InitDrawer();  
}

void CDesigner2D_Cloth::SetColor_CadFace(double r, double g, double b)
{
  face_color[0] = r;
  face_color[1] = g;
  face_color[2] = b;  
  const std::vector<unsigned int>& aIdL = cad_2d.GetAryElemID(Cad::LOOP);
  for(unsigned int iidl=0;iidl<aIdL.size();iidl++){
    unsigned int id_l = aIdL[iidl];
    cad_2d.SetColor_Loop(id_l,face_color);
  }  
  this->InitDrawer();
}

void CDesigner2D_Cloth::SetTextureCenter_FaceCAD(double cx, double cy)
{
  //  std::cout << cx << " " << cy << std::endl;

  if( pDrawerCAD == 0 ) return;
  pDrawerCAD->SetTexCenter(cx, cy);
}



void CDesigner2D_Cloth::FollowMshToCad_ifNeeded()
{
  if( setIdVCad_NeedFollow.empty() ){ return; }
  //    std::cout << "FollowMshToCad_ifNeeded()" << std::endl;
  std::vector<unsigned int> aIdV_del;
  std::set<unsigned int>::iterator itr = setIdVCad_NeedFollow.begin();
  for(;itr!=setIdVCad_NeedFollow.end();itr++){
    const unsigned int id_v = *itr;
//    std::cout << "follow " << id_v << std::endl;
    if( !cad_2d.IsElemID(Cad::VERTEX,id_v) ) continue;
    unsigned int itype_ope = 0;
    if( mesh_2d.FitMeshToCad_Vertex(cad_2d, id_v, itype_ope) ){ aIdV_del.push_back(id_v); }
    if( itype_ope == 0 ){
      std::cout << "fail folow mesh to cad rebuild mesh" << std::endl;
      this->Solve_fromCad_InterpValue();
      this->Msh_PrecompDrag(itype_cad_elem_prec, id_cad_elem_prec);      
      is_updated_coord = false;
      is_updated_edge = false;
      setIdVCad_NeedFollow.clear();
      return;
    }
    is_updated_coord = is_updated_coord || (itype_ope&1);	// êﬂì_ÇÃà⁄ìÆ
    is_updated_edge  = is_updated_edge  || (itype_ope&2);	// óvëfÇÃêÿÇËë÷Ç¶
  }
  for(unsigned int iid_v=0;iid_v<aIdV_del.size();iid_v++){
    setIdVCad_NeedFollow.erase( aIdV_del[iid_v] );
  } 
}

void CDesigner2D_Cloth::SetSelection(const std::vector<Com::View::SSelectedObject>& aSelecObj)
{
	pDrawerCAD->ClearSelected();
	m_id_cad_part = 0;  m_itype_cad_part = Cad::NOT_SET;
	if( aSelecObj.size() == 0 ) return;
	pDrawerCAD->GetCadPartID(aSelecObj[0].name,m_itype_cad_part,m_id_cad_part);
	pDrawerCAD->AddSelected( aSelecObj[0].name );
	m_picked_x = aSelecObj[0].picked_pos.x;
	m_picked_y = aSelecObj[0].picked_pos.y;
}

void CDesigner2D_Cloth::Cad_GetPicked
(Cad::CAD_ELEM_TYPE& itype, unsigned int& id,
 double& x, double& y) const
{
  id = m_id_cad_part;
  itype = m_itype_cad_part;
  x = m_picked_x;
  y = m_picked_y;  
}

bool CDesigner2D_Cloth::Cad_Move
(Cad::CAD_ELEM_TYPE itype_cad_part, unsigned int id_cad_part, 
 const Com::CVector2D& pos_pre, const Com::CVector2D& pos_cur, double tor)
{
  if( !cad_2d.IsElemID(itype_cad_part,id_cad_part) ){ return false; }
  if( pAnalysis != 0 ){
    pAnalysis->Update_Boundary_Condition_Cad_Move(itype_cad_part,id_cad_part,
                                                  pos_cur.x-pos_pre.x, pos_cur.y-pos_pre.y,
                                                  slider_deform);
  }  
  unsigned int itype_ope = 0;  
	if( itype_cad_part == Cad::VERTEX ){
    if( !this->MoveVertex(id_cad_part,pos_cur) ){ return false; }
		is_updated_cad = true;
//		if( mesh_2d.GetElemID_FromCadID(id_cad_part,Cad::VERTEX ) == 0 ) return true;
    // if this vertex doesn't touch any loop that have mesh, we don't have to move mesh and sensitivity
		if( itype_cad_elem_prec!=Cad::VERTEX || id_cad_elem_prec!=id_cad_part){ this->Msh_PrecompDrag(Cad::VERTEX, id_cad_part); }        
		if( !mesh_2d.FitMeshToCad_UsingPrecomp(cad_2d, Cad::VERTEX,id_cad_part, itype_ope) ){ goto FAIL_FIT_MESH; }
	}
	if( itype_cad_part == Cad::EDGE ){
		if( !this->MoveEdge(id_cad_part, pos_cur-pos_pre) ){ return false; }
		is_updated_cad = true;
//		if( mesh_2d.GetElemID_FromCadID(id_cad_part,Cad::EDGE ) == 0 ) return true;
    // if this edge doesn't touch any loop that have mesh, we don't have to move mesh and sensitivity
		if( itype_cad_elem_prec!=Cad::EDGE || id_cad_elem_prec!=id_cad_part){ this->Msh_PrecompDrag(Cad::EDGE, id_cad_part); }    
		if( !mesh_2d.FitMeshToCad_UsingPrecomp(cad_2d, Cad::EDGE,id_cad_part, itype_ope) ){ goto FAIL_FIT_MESH; }
	}
	if( itype_cad_part == Cad::LOOP ){
		if( !this->MoveLoop(id_cad_part,pos_cur-pos_pre) ){ return false; }
		is_updated_cad = true;
		if( itype_cad_elem_prec!=Cad::LOOP || id_cad_elem_prec!=id_cad_part){ this->Msh_PrecompDrag(Cad::LOOP, id_cad_part); }
		if( !mesh_2d.FitMeshToCad_UsingPrecomp(cad_2d, Cad::LOOP,id_cad_part, itype_ope) ){ goto FAIL_FIT_MESH; }
	}
  is_updated_coord = is_updated_coord || (itype_ope&1);
  is_updated_edge  = is_updated_edge  || (itype_ope&2);  
  return true;
FAIL_FIT_MESH:
  this->Solve_fromCad_InterpValue();
  this->Msh_PrecompDrag(itype_cad_part, id_cad_part);  
  is_updated_coord = false;
  is_updated_edge = false;
  is_updated_cad = false;  
//  return true;
//  for(unsigned int i=0;i<1;i++){ this->Solve_ifNeeded(); }
  if(pAnalysis != 0 ){
    /*
    if( pAnalysis->GetMode() != CLOTH_INITIAL_LOCATION ){    
      std::vector<double> har;
      this->Msh_GetHarmonicPrecomp(har);  // Get mesh harmonic function
      if( har.size() == 0 ) return true;
      pAnalysis->SetSensitivity(itype_cad_part,id_cad_part,har,pos_cur.x,pos_cur.y);
    } 
     */
  }
  return true;
//  */
}

////////////////////////////////////////////////////////////////
// à ëäïœâª


bool CDesigner2D_Cloth::Cad_Remove(Cad::CAD_ELEM_TYPE itype, unsigned int id)
{
	if( itype == Cad::LOOP ){
    if( mesh_2d.IsIdLCad_CutMesh(id) ) return false;    
    std::vector<unsigned int> aIdE, aIdV;
    for(Cad::CBRepSurface::CItrLoop pItr = cad_2d.GetItrLoop(id);!pItr.IsEnd();pItr++){
      unsigned int id_e;  bool is_same_dir;
      pItr.GetIdEdge(id_e,is_same_dir);
      aIdE.push_back(id_e);
      aIdV.push_back(pItr.GetIdVertex());
    }        
    for(unsigned int iie=0;iie<aIdE.size();iie++){
      cad_2d.RemoveElement(Cad::EDGE,aIdE[iie]);
    }
    for(unsigned int iiv=0;iiv<aIdV.size();iiv++){
      cad_2d.RemoveElement(Cad::VERTEX,aIdV[iiv]);
    }    
  }
	if( itype == Cad::VERTEX ){  
		if( !cad_2d.RemoveElement(itype,id) ){ return false; }
  }  
  if( itype == Cad::EDGE ){
    unsigned int id_l_l, id_l_r;
    cad_2d.GetIdLoop_Edge(id_l_l,id_l_r, id);
		const bool is_l_l_cut = mesh_2d.IsIdLCad_CutMesh(id_l_l);
		const bool is_l_r_cut = mesh_2d.IsIdLCad_CutMesh(id_l_r);
		if( !cad_2d.RemoveElement(itype,id) ) return false;
		if( is_l_l_cut != is_l_r_cut ){	// àÍï˚Ç™ÉÅÉbÉVÉÖÇ™Ç ÇËÅCàÍï˚Ç…ÉÅÉbÉVÉÖÇ™Ç»Ç¢èÍçáÇÕÅCÇ≈Ç´ÇÈÇæÇØÉÅÉbÉVÉÖÇçÏÇÈ
			if( cad_2d.IsElemID(Cad::LOOP,id_l_l) ){ mesh_2d.AddIdLCad_CutMesh(id_l_l); }
			if( cad_2d.IsElemID(Cad::LOOP,id_l_r) ){ mesh_2d.AddIdLCad_CutMesh(id_l_r); }
		}
		// ÉÅÉbÉVÉÖÇêÿÇÈÉãÅ[ÉvÇ™ë∂ç›ÇµÇ»ÇØÇÍÇŒè¡Ç∑
		const std::vector<unsigned int>& aIdL = mesh_2d.GetIdLCad_CutMesh();
		for(unsigned int iil=0;iil<aIdL.size();iil++){
			unsigned int id_l = aIdL[iil];
			if( cad_2d.IsElemID(Cad::LOOP,id_l) ){ continue; }
			mesh_2d.RemoveIdLCad_CutMesh(id_l);
			break;	// Ç±Ç±Ç≈breakÇµÇ»Ç¢Ç∆ÉGÉâÅ[Ç…Ç»ÇÈ
		}  
  }
  this->Solve_fromCad_InterpValue();
  return true;
}


bool CDesigner2D_Cloth::Cad_SetCurveType(unsigned int id_e, unsigned int itype )
{
  //    std::cout << "CGuiListner_Analysis2D_Interactive::Cad_SetCurveType" << std::endl;
  if( !cad_2d.IsElemID(Cad::EDGE,id_e) ) return false;
  if(      itype == 0 ){ is_updated_cad = cad_2d.SetCurve_Line(id_e); }
  else if( itype == 1 ){ is_updated_cad = cad_2d.SetCurve_Arc( id_e); }
  else if( itype == 2 ){ is_updated_cad = cad_2d.SetCurve_Polyline(id_e); }
	if( !is_updated_cad ) return false;
  ////////////////
	if(  m_is_solve_cad_change ){	// ÉÅÉbÉVÉÖÇìÆÇ©Ç∑
		unsigned int itype_ope;
		if( !mesh_2d.FitMeshToCad_Edge(cad_2d,id_e,   itype_ope) ){
			unsigned int id_vs, id_ve;
			cad_2d.GetIdVertex_Edge(id_vs,id_ve,id_e);
			setIdVCad_NeedFollow.insert( id_vs );
			setIdVCad_NeedFollow.insert( id_ve );
		}
		is_updated_coord = is_updated_coord || (itype_ope&1);	// êﬂì_ÇÃà⁄ìÆ
		is_updated_edge  = is_updated_edge  || (itype_ope&2);	// óvëfÇÃêÿÇËë÷Ç¶
	}
	return true;  
}



bool CDesigner2D_Cloth::Cad_DragArc(unsigned int id_e, const Com::CVector2D& pos_obj)
{
	if( !cad_2d.IsElemID(Cad::EDGE,id_e) ) return false;
	if( cad_2d.GetEdgeCurveType(id_e) != 1 ) return true;
  if( !cad_2d.DragArc(id_e,pos_obj) ) return false;
  is_updated_cad = true;
  if( m_is_solve_cad_change ){
    unsigned int itype_ope;
    if( !mesh_2d.FitMeshToCad_Edge(cad_2d,id_e,   itype_ope) ){
      unsigned int id_vs, id_ve;
      cad_2d.GetIdVertex_Edge(id_vs,id_ve,id_e);
      setIdVCad_NeedFollow.insert( id_vs );
      setIdVCad_NeedFollow.insert( id_ve );
    }
    is_updated_coord = is_updated_coord || (itype_ope&1);	// êﬂì_ÇÃà⁄ìÆ
    is_updated_edge  = is_updated_edge  || (itype_ope&2);	// óvëfÇÃêÿÇËë÷Ç¶
    is_updated_cad = true;
  }
  return true;  
}

bool CDesigner2D_Cloth::Cad_DragPolyline(unsigned int id_e, const Com::CVector2D& dist)
{
	if( !cad_2d.IsElemID(Cad::EDGE,id_e) ) return false;
	if( cad_2d.GetEdgeCurveType(id_e) != 2 ) return true;
  if( !cad_2d.DragPolyline(id_e,dist) ) return false;
  std::cout << "Cad_DragPolyline" << std::endl;
  is_updated_cad = true;
  if( m_is_solve_cad_change ){
    unsigned int itype_ope;
    if( !mesh_2d.FitMeshToCad_Edge(cad_2d,id_e,   itype_ope) ){
      this->Solve_fromCad_InterpValue();
      is_updated_coord = false;
      is_updated_edge = false;
      is_updated_cad = false;  
      return true;    
    }
    is_updated_coord = is_updated_coord || (itype_ope&1);	// êﬂì_ÇÃà⁄ìÆ
    is_updated_edge  = is_updated_edge  || (itype_ope&2);	// óvëfÇÃêÿÇËë÷Ç¶
    is_updated_cad = true;
  }
  return true;    
}

bool CDesigner2D_Cloth::Cad_PreCompDragPolyline(unsigned int id_e, const Com::CVector2D& pick_pos)
{  
	if( !cad_2d.IsElemID(Cad::EDGE,id_e) ) return false;
	if( cad_2d.GetEdgeCurveType(id_e) != 2 ) return true;
  if( !cad_2d.PreCompDragPolyline(id_e,pick_pos) ) return false;
  is_updated_cad = true;
  itype_cad_elem_prec = Cad::NOT_SET;
  id_cad_elem_prec = 0;
  isnt_mesh_deform_sensitivity = true;
  if( m_is_solve_cad_change ){
    unsigned int itype_ope;
    if( !mesh_2d.FitMeshToCad_Edge(cad_2d,id_e,   itype_ope) ){
      unsigned int id_vs, id_ve;
      cad_2d.GetIdVertex_Edge(id_vs,id_ve,id_e);
      setIdVCad_NeedFollow.insert( id_vs );
      setIdVCad_NeedFollow.insert( id_ve );
    }
    is_updated_coord = is_updated_coord || (itype_ope&1);	// êﬂì_ÇÃà⁄ìÆ
    is_updated_edge  = is_updated_edge  || (itype_ope&2);	// óvëfÇÃêÿÇËë÷Ç¶
    is_updated_cad = true;
  }
  return true;  
}

// smoothing edge (id_e) if radius is negative smooth whole edge 
bool CDesigner2D_Cloth::Cad_SmoothingPolylineEdge(unsigned int id_e, unsigned int niter,
                                                  const Com::CVector2D& pos, double radius)
{
	if( !cad_2d.IsElemID(Cad::EDGE,id_e) ) return false;
	if( cad_2d.GetEdgeCurveType(id_e) != 2 ) return true;
  if( !cad_2d.SmoothingPolylineEdge(id_e,niter,pos,radius) ) return false;
  is_updated_cad = true;
  if( m_is_solve_cad_change ){
    unsigned int itype_ope;
    if( !mesh_2d.FitMeshToCad_Edge(cad_2d,id_e,   itype_ope) ){
      unsigned int id_vs, id_ve;
      cad_2d.GetIdVertex_Edge(id_vs,id_ve,id_e);
      setIdVCad_NeedFollow.insert( id_vs );
      setIdVCad_NeedFollow.insert( id_ve );
    }
    is_updated_coord = is_updated_coord || (itype_ope&1);	// êﬂì_ÇÃà⁄ìÆ
    is_updated_edge  = is_updated_edge  || (itype_ope&2);	// óvëfÇÃêÿÇËë÷Ç¶
    is_updated_cad = true;
  }
  return true;  
}



void CDesigner2D_Cloth::Draw(unsigned int idraw_object)
{
	if( this->is_updated_cad ){
		this->is_updated_cad = false;
		pDrawerCAD->UpdateCAD_TopologyGeometry(cad_2d);
	}
	if(      idraw_object == 0 ){ // analysis and cad_rigid_loop, cad_edge and cad_vertex
    if( pAnalysis != 0 ){
      pAnalysis->UpdateAnimationDrawer();
      pAnalysis->Draw(); 
    }
		else{ 
      pDrawerMsh->Draw(); 
    }
		////////////////
		{
			const std::vector<unsigned int>& aIdL = this->cad_2d.GetAryElemID(Cad::LOOP);
			for(unsigned int i=0;i<aIdL.size();i++){
				const unsigned int id_l = aIdL[i];
				if( pAnalysis != 0 ){ pDrawerCAD->SetIsShow(false,Cad::LOOP,id_l); }
				else{                 pDrawerCAD->SetIsShow(true,Cad::LOOP,id_l);  }
			}
		}
		{
			const std::vector<unsigned int>& aIdL = mesh_2d.GetIdLCad_CutMesh();
			for(unsigned int iil=0;iil<aIdL.size();iil++){
				const unsigned int id_l = aIdL[iil];
				pDrawerCAD->SetIsShow(false,Cad::LOOP,id_l);
			}
		}
    pDrawerCAD->Draw();
	}
	else if( idraw_object == 1 ){ // cad
    pDrawerCAD->SetIsShow(false,Cad::LOOP, this->cad_2d.GetAryElemID(Cad::LOOP) );
    const std::vector<unsigned int>& aIdL = mesh_2d.GetIdLCad_CutMesh();
    for(unsigned int iil=0;iil<aIdL.size();iil++){
      const unsigned int id_l = aIdL[iil];
      pDrawerCAD->SetIsShow(true,Cad::LOOP,id_l);
    }
		pDrawerCAD->Draw();
	}
	else if( idraw_object == 2 ){ // mesh
    pDrawerMsh->Draw();
	}
	else if( idraw_object == 3 ){ // mesh and cad edge & vertex
    pDrawerMsh->Draw();
    pDrawerCAD->SetIsShow(false,Cad::LOOP, this->cad_2d.GetAryElemID(Cad::LOOP) );    
    /*
    pDrawerCAD->SetIsShow(true,Cad::LOOP, this->cad_2d.GetAryElemID(Cad::LOOP) );
		const std::vector<unsigned int>& aIdL = mesh_2d.GetIdLCad_CutMesh();
		for(unsigned int iil=0;iil<aIdL.size();iil++){
      const unsigned int id_l = aIdL[iil];
      pDrawerCAD->SetIsShow(false,Cad::LOOP,id_l);
    }
     */
    pDrawerCAD->Draw();
	}
	else if( idraw_object == 4 ){
    if( pAnalysis != 0 ){     
      pAnalysis->UpdateAnimationDrawer();
      pAnalysis->Draw(); 
    }
	}
  else if( idraw_object == 5 ){
    if( pAnalysis != 0 ){ pAnalysis->DrawBoundaryCondition(cad_2d); }
  }
}

void CDesigner2D_Cloth::InitDrawer()
{
//  std::cout << "Init Drawer Cloth Designer" << std::endl;
  double tex_cent_x=0, tex_cent_y=0;
	if( pDrawerCAD != 0 ){ 
    pDrawerCAD->GetTexCenter(tex_cent_x, tex_cent_y);
    delete pDrawerCAD; 
  }
  pDrawerCAD = new Cad::View::CDrawer_Cad2D(cad_2d);
  pDrawerCAD->EnableUVMap(true);  
  pDrawerCAD->SetAntiAliasing(true);
  pDrawerCAD->SetTextureScale(tex_scale);  
  pDrawerCAD->SetTexCenter(tex_cent_x,tex_cent_y);
	if( pDrawerMsh != 0 ){ delete pDrawerMsh; }
	pDrawerMsh = new Msh::View::CDrawerMsh2D(mesh_2d);
//	if( pAnalysis != 0 ){ pAnalysis->InitDrawer(); }
}


void CDesigner2D_Cloth::SetAnalysisInitialize(CAnalysis2D_Cloth_Static* pAnalysis,unsigned int iprob)
{
  this->setIdVCad_NeedFollow.clear();
	this->pAnalysis = pAnalysis;
	if( pAnalysis != 0 ){        
		pAnalysis->SetModelProblem_Cloth(cad_2d,mesh_2d,iprob,slider_deform,aSymIdVPair); 
	}
	else{
    aSymIdVPair.clear();
    if( iprob == 0 ){
      unsigned int id_l0;
      {	// Make model
        cad_2d.Clear();
        std::vector<Com::CVector2D> vec_ary;
        vec_ary.push_back( Com::CVector2D(0.0,0.0) );
        vec_ary.push_back( Com::CVector2D(1.0,0.0) );
        vec_ary.push_back( Com::CVector2D(1.0,1.0) );
        vec_ary.push_back( Com::CVector2D(0.0,1.0) );
        id_l0 = cad_2d.AddPolygon( vec_ary ).id_l_add;
      }
      {   // ÉÅÉbÉVÉÖÇÃê›íË
        mesh_2d.Clear();			
        mesh_2d.AddIdLCad_CutMesh(id_l0);
        mesh_2d.SetMeshingMode_ElemSize(2000);
        mesh_2d.Meshing(cad_2d);
      }
    }
    else if( iprob == 1 ){
      unsigned int id_l0;
      {	// Make model
        cad_2d.Clear();
        std::vector<Com::CVector2D> vec_ary;
        vec_ary.push_back( Com::CVector2D(+0.00,0.00) );  // 0
        vec_ary.push_back( Com::CVector2D(+1.00,0.00) );  // 1
        vec_ary.push_back( Com::CVector2D(+1.00,1.00) );  // 2
        vec_ary.push_back( Com::CVector2D(+1.35,0.83) );  // 3
        vec_ary.push_back( Com::CVector2D(+1.50,1.2) );   // 4
        vec_ary.push_back( Com::CVector2D(+0.80,1.5) );   // 5
        ////
        vec_ary.push_back( Com::CVector2D(+0.20,1.50) );  // 6
        vec_ary.push_back( Com::CVector2D(-0.50,1.20) );  // 7
        vec_ary.push_back( Com::CVector2D(-0.35,0.83) );  // 8
        vec_ary.push_back( Com::CVector2D(+0.00,1.00) );  // 9
        Cad::CCadObj2D::CResAddPolygon res = cad_2d.AddPolygon( vec_ary );
        id_l0 = res.id_l_add;
        cad_2d.SetCurve_Polyline(res.aIdE[5]);
        cad_2d.PreCompDragPolyline(res.aIdE[5],vec_ary[5]*0.5+vec_ary[6]*0.5);
        cad_2d.DragPolyline(res.aIdE[5],vec_ary[5]*0.5+vec_ary[6]*0.5+Com::CVector2D(-0.0,-0.08));
        ////
        aSymIdVPair.push_back( std::make_pair(res.aIdV[1],res.aIdV[0]) );
        aSymIdVPair.push_back( std::make_pair(res.aIdV[2],res.aIdV[9]) );        
        aSymIdVPair.push_back( std::make_pair(res.aIdV[3],res.aIdV[8]) );                
        aSymIdVPair.push_back( std::make_pair(res.aIdV[4],res.aIdV[7]) );        
        aSymIdVPair.push_back( std::make_pair(res.aIdV[5],res.aIdV[6]) );                
      }
      {   // ÉÅÉbÉVÉÖÇÃê›íË
        mesh_2d.Clear();			
        mesh_2d.AddIdLCad_CutMesh(id_l0);
        mesh_2d.SetMeshingMode_ElemSize(4000);
        mesh_2d.Meshing(cad_2d);
      }
    }
    else if( iprob == 2 ){
      unsigned int id_l0;
      {	// Make model
        cad_2d.Clear();
        std::vector<Com::CVector2D> vec_ary;
        vec_ary.push_back( Com::CVector2D(+0.00,0.00) );  // 0
        vec_ary.push_back( Com::CVector2D(+1.00,0.00) );  // 1
        vec_ary.push_back( Com::CVector2D(+1.00,1.00) );  // 2
        //      vec_ary.push_back( Com::CVector2D(+0.90,1.10) );      
        vec_ary.push_back( Com::CVector2D(+0.90,1.50) );  // 3
        vec_ary.push_back( Com::CVector2D(+0.80,1.50) );  // 4
        ////
        vec_ary.push_back( Com::CVector2D(+0.20,1.50) );  // 5
        vec_ary.push_back( Com::CVector2D(+0.10,1.50) );  // 6
        //      vec_ary.push_back( Com::CVector2D(+0.10,1.10) );
        vec_ary.push_back( Com::CVector2D(+0.00,1.00) );  // 7
        Cad::CCadObj2D::CResAddPolygon res = cad_2d.AddPolygon( vec_ary );
        id_l0 = res.id_l_add;
        cad_2d.SetCurve_Polyline(res.aIdE[4]);
        cad_2d.PreCompDragPolyline(res.aIdE[4],vec_ary[4]*0.5+vec_ary[5]*0.5);
        cad_2d.DragPolyline(res.aIdE[4],vec_ary[4]*0.5+vec_ary[5]*0.5+Com::CVector2D(-0.0,-0.08));
        ////
        aSymIdVPair.push_back( std::make_pair(res.aIdV[1],res.aIdV[0]) );
        aSymIdVPair.push_back( std::make_pair(res.aIdV[2],res.aIdV[7]) );        
        aSymIdVPair.push_back( std::make_pair(res.aIdV[3],res.aIdV[6]) );                
        aSymIdVPair.push_back( std::make_pair(res.aIdV[4],res.aIdV[5]) );
      }
      {   // ÉÅÉbÉVÉÖÇÃê›íË
        mesh_2d.Clear();			
        mesh_2d.AddIdLCad_CutMesh(id_l0);
        mesh_2d.SetMeshingMode_ElemSize(1000);
        mesh_2d.Meshing(cad_2d);
      }
      Cad_AddCutLine(Com::CVector2D(0.5,1), Com::CVector2D(0.5,0.2));
    }      
	}  
  this->InitDrawer();
  /*
	if( pDrawerCAD != 0 ){ delete pDrawerCAD; }
  pDrawerCAD = new Cad::View::CDrawer_Cad2D(cad_2d);
	if( pDrawerMsh != 0 ){ delete pDrawerMsh; }
	pDrawerMsh = new Msh::View::CDrawerMsh2D(mesh_2d);
   */
}

void CDesigner2D_Cloth::SetAnalysis(CAnalysis2D_Cloth_Static* pAnalysis){
	this->pAnalysis = pAnalysis;
  this->InitDrawer();
  /*
	if( pDrawerCAD != 0 ){ delete pDrawerCAD; }
    pDrawerCAD = new Cad::View::CDrawer_Cad2D(cad_2d);
	if( pDrawerMsh != 0 ){ delete pDrawerMsh; }
	pDrawerMsh = new Msh::View::CDrawerMsh2D(mesh_2d);
   */
}

void CDesigner2D_Cloth::Solve_fromCad_InitValue()
{
	mesh_2d.Meshing(cad_2d);
	if( pAnalysis != 0 ){ pAnalysis->BuildFEM_ClearValueField(cad_2d,mesh_2d); }
	this->InitDrawer();
}

void CDesigner2D_Cloth::Solve_fromCad_InterpValue(const std::vector< std::pair<unsigned int,unsigned int> >& aNewL)
{
	mesh_2d.Meshing(cad_2d);
	if( pAnalysis != 0 ){ pAnalysis->BuildFEM_InterpValueField(cad_2d,mesh_2d,aNewL); }
	this->InitDrawer();
}


void CDesigner2D_Cloth::Solve_ifNeeded(){
	/*
	{	// bulding mesh for each computation
		mesh_2d.Meshing(cad_2d);
		if( pAnalysis != 0 ){ pAnalysis->SolveInitial(cad_2d,mesh_2d); }
		this->is_updated_cad = false;
		this->is_updated_coord = false;
		this->is_updated_edge = false;
		if( pDrawerMsh != 0 ){ delete pDrawerMsh; }
		pDrawerMsh = new Msh::View::CDrawerMsh2D(mesh_2d);
		return;
	}	
	*/
	if( is_updated_cad && !(is_updated_coord||is_updated_edge) ){
    std::cout << "fail mesh update" << std::endl;
		mesh_2d.Meshing(cad_2d); 
    const std::vector< std::pair<unsigned int,unsigned int> > aNewL;
		if( pAnalysis != 0 ){ pAnalysis->BuildFEM_InterpValueField(cad_2d,mesh_2d,aNewL); }    
//		if( pAnalysis != 0 ){ pAnalysis->BuildFEM_ClearValueField(cad_2d,mesh_2d); }        
    this->InitDrawer();    
		is_updated_cad = false;    
    is_updated_coord = false;
    is_updated_edge = false;    
    this->Msh_PrecompDrag(itype_cad_elem_prec, id_cad_elem_prec);    
		return;
	}
	if( is_updated_coord || is_updated_edge ){
		if( pDrawerMsh != 0 ){ delete pDrawerMsh; }
		pDrawerMsh = new Msh::View::CDrawerMsh2D(mesh_2d);
	}
//	if( !m_is_solve_cad_change ){ return; }
	if( pAnalysis == 0 ) return;
	SOLVER_FLAG solver_flg = this->pAnalysis->UpdateMeshAndSolve(mesh_2d,
                                                               is_updated_coord,is_updated_edge,
                                                               m_is_solve_cad_change,isnt_mesh_deform_sensitivity);
	if( solver_flg != SUCCESS ){
		mesh_2d.Meshing(cad_2d);
    const std::vector< std::pair<unsigned int,unsigned int> > aNewL;
		this->pAnalysis->BuildFEM_InterpValueField(cad_2d,mesh_2d,aNewL);
	}
	is_updated_coord = false;
	is_updated_edge  = false;
}

void CDesigner2D_Cloth::Serialize( Com::CSerializer& arch )
{
//	if( !arch.IsOpen() ) return;
	if( arch.IsLoading() ){
		const unsigned int buff_size = 256;
		char class_name[buff_size];
		arch.ReadDepthClassName(class_name,buff_size);
		std::cout << class_name << std::endl;
		assert( strncmp(class_name,"CGuiListner_Analysis2D_Interactive",34) == 0 );
    arch.ShiftDepth(true);
		cad_2d.Serialize(arch);
		mesh_2d.Serialize(arch,true);
		arch.ShiftDepth(false);
		mesh_2d.Meshing(cad_2d);
	}
	else{
		arch.WriteDepthClassName("CGuiListner_Analysis2D_Interactive");
    arch.ShiftDepth(true);
		cad_2d.Serialize(arch);
		mesh_2d.Serialize(arch,true);
    arch.ShiftDepth(false);
	}
}

unsigned int CDesigner2D_Cloth::AddDartDiamond
(unsigned int id_l,const Com::CVector2D& pos, const Com::CVector2D& poe)
{
  if( !mesh_2d.IsIdLCad_CutMesh(id_l) ) return 0;
  Com::CVector2D ve(-(pos-poe).y,(pos-poe).x);
  ve.Normalize();
  ve *= 0.025;
  Com::CVector2D c = (pos+poe)*0.5;
  Com::CVector2D po1 = c+ve;
  Com::CVector2D po2 = c-ve;
  if( !cad_2d.CheckIsPointInsideLoop(id_l,po1)  ) return 0;
  if( !cad_2d.CheckIsPointInsideLoop(id_l,po2)  ) return 0;        
  std::vector< Com::CVector2D > aVec;
  aVec.push_back(pos);
  aVec.push_back(po1);
  aVec.push_back(poe);
  aVec.push_back(po2);
  //    const unsigned int id_l_add = cad_2d.AddPolygon(aVec,id_l).id_l_add;
  const Cad::CCadObj2D::CResAddPolygon& res = cad_2d.AddPolygon(aVec,id_l);
  const unsigned int id_l_add = res.id_l_add;
  if( id_l_add == 0 ){ return 0; }
  CDartDiamond dart;
  {
    dart.id_l = res.id_l_add;
    dart.id_vu = res.aIdV[0];
    dart.id_vl = res.aIdV[1];
    dart.id_vd = res.aIdV[2];
    dart.id_vr = res.aIdV[3];
    dart.id_e_ul = res.aIdE[0];
    dart.id_e_ld = res.aIdE[1];
    dart.id_e_dr = res.aIdE[2];
    dart.id_e_ru = res.aIdE[3];          
  }
  aDartDiamond.push_back(dart);
  Solve_fromCad_InterpValue(); 	
  this->InitDrawer();     
  if( pAnalysis != 0 ){
    pAnalysis->PerformStaticSolver();
    pAnalysis->ConnectEdge(this->cad_2d,this->mesh_2d,res.aIdE[0],res.aIdE[3]);
    pAnalysis->ConnectEdge(this->cad_2d,this->mesh_2d,res.aIdE[1],res.aIdE[2]);    
  }
  return id_l_add;  
}

unsigned int CDesigner2D_Cloth::AddDartEdge
(unsigned int id_l, unsigned int id_e0, bool is_same_dir0,
 const Com::CVector2D& pos, const Com::CVector2D& p_nearest)
{
  if( !cad_2d.IsElemID(Cad::EDGE,id_e0) ) return 0;
  if( cad_2d.GetEdgeCurveType(id_e0) != 0 ) return 0;
  Com::CVector2D out0,out1;    
  double ratio0, ratio1;
  {
    const Cad::CEdge2D& edge = cad_2d.GetEdge(id_e0);
    bool is_exceed0, is_exceed1;
    edge.GetPointOnCurve_OnCircle(p_nearest, 0.025, true, is_exceed0, out0);    
    edge.GetPointOnCurve_OnCircle(p_nearest, 0.025,false, is_exceed1, out1);
    if( is_exceed0 ) return 0;
    if( is_exceed1 ) return 0;          
    ratio0 = Com::Distance(p_nearest,edge.po_s);
    ratio1 = Com::Distance(p_nearest,edge.po_e);
    const double invsum01 = 1.0/(ratio0+ratio1);
    ratio0 *= invsum01;
    ratio1 *= invsum01;      
  }
  CDartOnEdge dart;        
  unsigned int id_e_add0 = 0;
  {
    dart.id_vc = cad_2d.AddVertex(Cad::LOOP,id_l,pos).id_v_add;
    Cad::CCadObj2D::CResAddVertex res_add_out0 = cad_2d.AddVertex(Cad::EDGE,id_e0,out0);
    dart.id_v1 = res_add_out0.id_v_add;
    id_e_add0 = res_add_out0.id_e_add;
    dart.id_v2 = cad_2d.AddVertex(Cad::EDGE,id_e0,out1).id_v_add;
    dart.id_e1 = cad_2d.ConnectVertex_Line(dart.id_v1,dart.id_vc).id_e_add;
    Cad::CBRepSurface::CResConnectVertex res;
    if( is_same_dir0 ){ res = cad_2d.ConnectVertex_Line(dart.id_v2,dart.id_vc); }
    else{ res = cad_2d.ConnectVertex_Line(dart.id_vc,dart.id_v2); }
    dart.id_e2 = res.id_e_add;
    dart.id_l = res.id_l_add;
    dart.id_l_in = id_l;
    dart.id_vo = 0;
  }    
  if( pAnalysis != 0 ){
    unsigned int id_e_add1 = 0;
    unsigned int id_e_opp;  bool is_same_dir_opp;
    bool is_seam_line_cutted = false;
    if( pAnalysis->IsSeamLine(id_e0,id_e_opp,is_same_dir_opp) ){
      is_seam_line_cutted = true;
      assert( cad_2d.IsElemID(Cad::EDGE,id_e_opp) );
      unsigned int id_v_opp0, id_v_opp1;
      cad_2d.GetIdVertex_Edge(id_v_opp0,id_v_opp1,id_e_opp);
      Com::CVector2D opp0 = cad_2d.GetVertexCoord(id_v_opp0);
      Com::CVector2D opp1 = cad_2d.GetVertexCoord(id_v_opp1);        
      Cad::CCadObj2D::CResAddVertex res1 = cad_2d.AddVertex(Cad::EDGE,id_e_opp,opp0*ratio0+opp1*ratio1);
      id_e_add1 = res1.id_e_add;
      dart.id_vo = res1.id_v_add;
    }
    aDartOnEdge.push_back(dart);       
    Solve_fromCad_InterpValue(); 	          
    if( is_seam_line_cutted ){
      pAnalysis->DisconnectEdge(id_e0,id_e_opp);
      pAnalysis->ConnectEdge(this->cad_2d,this->mesh_2d,id_e0,id_e_add1); 
      pAnalysis->ConnectEdge(this->cad_2d,this->mesh_2d,id_e_add0,id_e_opp);         
    }
    pAnalysis->ConnectEdge(this->cad_2d,this->mesh_2d,dart.id_e2,dart.id_e1);
    this->InitDrawer(); 
    if( pAnalysis->GetMode() != CLOTH_INITIAL_LOCATION ){
      pAnalysis->PerformStaticSolver();    
    }
  }
  else{
    aDartOnEdge.push_back(dart);       
    Solve_fromCad_InterpValue(); 	          
  }
  return dart.id_l;  
}

bool FindIntersectionEdge
(const Cad::CCadObj2D& cad_2d,
 unsigned int id_l,const Com::CVector2D& pos, const Com::CVector2D& poe,
 unsigned int& id_e_nearest, bool& is_same_dir_e_near,
 Com::CVector2D& p_nearest)
{  
  id_e_nearest = 0;
  Com::CVector2D dir = poe-pos; dir.Normalize();
  double sqdist = -1;
  for(Cad::CBRepSurface::CItrLoop pItr = cad_2d.GetItrLoop(id_l);!pItr.IsEndChild();pItr.ShiftChildLoop()){
    for(;!pItr.IsEnd();pItr++){
      unsigned int id_e;  bool is_same_dir;
      pItr.GetIdEdge(id_e,is_same_dir);
      const Cad::CEdge2D& edge = cad_2d.GetEdge(id_e);
      Com::CVector2D sec;
      if( !edge.GetNearestIntersectionPoint_AgainstHalfLine(sec,pos,dir) ) continue;
      if( sqdist < 0 || (sec-pos).SqLength() < sqdist ){
        id_e_nearest = id_e;
        is_same_dir_e_near = is_same_dir;
        p_nearest = sec;
        sqdist = (sec-pos).SqLength();
      }
    }
  }        
  if( !cad_2d.IsElemID(Cad::EDGE,id_e_nearest) ){ return false; }
  return true;
}


unsigned int GetSymPairIdV(unsigned int id_v, const std::vector< std::pair<unsigned int,unsigned int> >& aPair){
  for(unsigned int isym=0;isym<aPair.size();isym++){
    if(      aPair[isym].first  == id_v ){ return aPair[isym].second; }
    else if( aPair[isym].second == id_v ){ return aPair[isym].first;  } 
  }
  return 0;
}

bool CDesigner2D_Cloth::AddDartDiamond_Sym
(unsigned int id_l,const Com::CVector2D& pos0, const Com::CVector2D& poe0)
{
  unsigned int id_l_dart0 = AddDartDiamond(id_l,pos0,poe0);
  if( !cad_2d.IsElemID(Cad::LOOP,id_l_dart0) ) return false;
  /////
  double cnt_x;
  {
    double sum;
    unsigned int n=0;
    for(Cad::CBRepSurface::CItrLoop itrl = cad_2d.GetItrLoop(id_l);!itrl.IsEnd();itrl++){
      unsigned int id_v  = itrl.GetIdVertex();
      unsigned int id_vo = GetSymPairIdV(id_v,aSymIdVPair);
      if( id_vo == 0 ) continue;
      sum += cad_2d.GetVertexCoord(id_v).x;
      sum += cad_2d.GetVertexCoord(id_vo).x;      
      n+=2;
    }      
    cnt_x = sum/n;
  }
  Com::CVector2D pos1(2*cnt_x-pos0.x,pos0.y);
  Com::CVector2D poe1(2*cnt_x-poe0.x,poe0.y);
  if( (pos1.x-cnt_x)*(poe1.x-cnt_x) < 0 ) return true;
  unsigned int id_l_dart1 = AddDartDiamond(id_l,pos1,poe1);
  if( !cad_2d.IsElemID(Cad::LOOP,id_l_dart1) ) return false;
  
  int idart0=-1, idart1=-1;
  for(unsigned int idart=0;idart<aDartDiamond.size();idart++){
    if( aDartDiamond[idart].id_l == id_l_dart0 ){ idart0 = idart; }
    if( aDartDiamond[idart].id_l == id_l_dart1 ){ idart1 = idart; }
  }
  if( idart0 == -1 && idart1 == -1 ) return true;
  const CDartDiamond& dart0 = aDartDiamond[idart0];
  const CDartDiamond& dart1 = aDartDiamond[idart1];  
  aSymIdVPair.push_back( std::make_pair(dart0.id_vu,dart1.id_vu) );
  aSymIdVPair.push_back( std::make_pair(dart0.id_vd,dart1.id_vd) );
  aSymIdVPair.push_back( std::make_pair(dart0.id_vl,dart1.id_vr) );
  aSymIdVPair.push_back( std::make_pair(dart0.id_vr,dart1.id_vl) );  
  return true;
}

///////////////////////
bool CDesigner2D_Cloth::AddDartEdge_Sym
(unsigned id_l0,
 const Com::CVector2D& pos0, 
 const Com::CVector2D& poe0)
{  
  unsigned int id_e0=0; bool is_same_dir0;
  Com::CVector2D p0;
  if( !FindIntersectionEdge(cad_2d,id_l0,pos0,poe0, id_e0,is_same_dir0, p0) ) return 0;
  unsigned int id_v10,id_v20;
  cad_2d.GetIdVertex_Edge(id_v10,id_v20,id_e0);
  const unsigned int id_v1o = GetSymPairIdV(id_v10,aSymIdVPair);
  const unsigned int id_v2o = GetSymPairIdV(id_v20,aSymIdVPair);
  
  unsigned int id_l_dart0 = this->AddDartEdge(id_l0, id_e0,is_same_dir0, pos0,p0); // Add Dart0
  if( !cad_2d.IsElemID(Cad::EDGE,id_l_dart0) ){ return false; }
  
  if( id_v1o == 0 || id_v2o == 0 ) return true;
  double cnt_x;  
  unsigned int id_l_dart1 = 0;
  if( id_v1o == id_v20 ){ // dart in centered edge
    assert( id_v2o == id_v10 );
    {
      const Com::CVector2D& vec10 = cad_2d.GetVertexCoord(id_v10);
      const Com::CVector2D& vec20 = cad_2d.GetVertexCoord(id_v20);
      cnt_x = (vec10.x + vec20.x)*0.5;
    }        
    {
      Com::CVector2D pos1(2*cnt_x-pos0.x,pos0.y);
      Com::CVector2D poe1(2*cnt_x-poe0.x,poe0.y);
//      if( (pos1.x-cnt_x)*(poe1.x-cnt_x) < -0.001 ) return true;      
      {
        if( fabs(pos1.x-cnt_x) < 0.03 ) return true;
        if( fabs(poe1.x-cnt_x) < 0.03 ) return true;        
        if( (pos1.x-cnt_x)*(poe1.x-cnt_x) < 0 ) return true;        
      }
      unsigned int id_e1=0; bool is_same_dir1;
      Com::CVector2D p1;
      if( !FindIntersectionEdge(cad_2d,id_l0,pos1,poe1, id_e1,is_same_dir1, p1) ) return true;
      id_l_dart1 = this->AddDartEdge(id_l0, id_e1,is_same_dir1,  pos1,p1);  // Add Dart1
    }
  }
  else {  // dart in side edge
    {
      const Com::CVector2D& vec10 = cad_2d.GetVertexCoord(id_v10);
      const Com::CVector2D& vec20 = cad_2d.GetVertexCoord(id_v20);
      const Com::CVector2D& vec1o = cad_2d.GetVertexCoord(id_v1o);
      const Com::CVector2D& vec2o = cad_2d.GetVertexCoord(id_v2o);
      cnt_x = (vec10.x + vec20.x + vec1o.x + vec2o.x)*0.25;
    }
    {
      Com::CVector2D pos1(2*cnt_x-pos0.x,pos0.y);
      Com::CVector2D poe1(2*cnt_x-poe0.x,poe0.y);
      if( fabs(pos1.x-cnt_x) < 0.03 ) return true;
      if( fabs(poe1.x-cnt_x) < 0.03 ) return true;        
      if( (pos1.x-cnt_x)*(poe1.x-cnt_x) < 0 ) return true;              
      unsigned int id_l1 = 0;
      {
        const std::vector<unsigned int>& aIdL = cad_2d.GetAryElemID(Cad::LOOP);  
        for(unsigned int iidl=0;iidl<aIdL.size();iidl++){
          unsigned int id_l = aIdL[iidl];
          if( cad_2d.CheckIsPointInsideLoop(id_l,pos1) ){ id_l1 = id_l; break; }
        }    
      }
      unsigned int id_e1=0; bool is_same_dir1;
      Com::CVector2D p1;
      if( !FindIntersectionEdge(cad_2d,id_l1,pos1,poe1, id_e1,is_same_dir1, p1) ) return true;
      id_l_dart1 = this->AddDartEdge(id_l1, id_e1,is_same_dir1, pos1,p1);  // Add Dart1
    }
  }
  
  if( !cad_2d.IsElemID(Cad::EDGE,id_l_dart1) ){ return true; }
  
  int idart0=-1, idart1=-1;
  for(unsigned int idart=0;idart<aDartOnEdge.size();idart++){
    if( aDartOnEdge[idart].id_l == id_l_dart0 ){ idart0 = idart; }
    if( aDartOnEdge[idart].id_l == id_l_dart1 ){ idart1 = idart; }
  }
  if( idart0 == -1 && idart1 == -1 ) return true;
  const CDartOnEdge& dart0 = aDartOnEdge[idart0];
  const CDartOnEdge& dart1 = aDartOnEdge[idart1];  
  aSymIdVPair.push_back( std::make_pair(dart0.id_vc,dart1.id_vc) );
  aSymIdVPair.push_back( std::make_pair(dart0.id_vo,dart1.id_vo) );
  {
    const Com::CVector2D& vec10 = cad_2d.GetVertexCoord(dart0.id_v1);
    const Com::CVector2D& vec20 = cad_2d.GetVertexCoord(dart0.id_v2);
    const Com::CVector2D& vec11 = cad_2d.GetVertexCoord(dart1.id_v1);
    const Com::CVector2D& vec21 = cad_2d.GetVertexCoord(dart1.id_v2); 
    const Com::CVector2D vec10o(2*cnt_x-vec10.x,vec10.y);
    const Com::CVector2D vec20o(2*cnt_x-vec20.x,vec20.y);
    double sum0 = Com::Distance(vec10o,vec11) + Com::Distance(vec20o,vec21);    
    double sum1 = Com::Distance(vec10o,vec21) + Com::Distance(vec20o,vec11);
    if( sum0 < sum1 ){
      aSymIdVPair.push_back( std::make_pair(dart0.id_v1,dart1.id_v1) );
      aSymIdVPair.push_back( std::make_pair(dart0.id_v2,dart1.id_v2) );      
    }
    else{
      aSymIdVPair.push_back( std::make_pair(dart0.id_v1,dart1.id_v2) );      
      aSymIdVPair.push_back( std::make_pair(dart0.id_v2,dart1.id_v1) );            
    }
  }    
  return true;
}




bool CDesigner2D_Cloth::Cad_AddHole(const Com::CVector2D& pos, const Com::CVector2D& poe)
{
  //  pAnalysis->SetClothPiecePlacingMode();
  unsigned int id_ls = 0, id_le = 0;
  {
    const std::vector<unsigned int>& aIdL = cad_2d.GetAryElemID(Cad::LOOP);  
    for(unsigned int iidl=0;iidl<aIdL.size();iidl++){
      unsigned int id_l = aIdL[iidl];
      if( cad_2d.CheckIsPointInsideLoop(id_l,pos) ){ id_ls = id_l; break; }
    }    
    for(unsigned int iidl=0;iidl<aIdL.size();iidl++){
      unsigned int id_l = aIdL[iidl];
      if( cad_2d.CheckIsPointInsideLoop(id_l,poe) ){ id_le = id_l; break; }
    }      
  }    
  ////
  if( id_ls == 0 || id_le == 0 || id_ls != id_le ){ return true; }  
  const unsigned int id_l = id_ls;
  if( !mesh_2d.IsIdLCad_CutMesh(id_l) ) return 0;
  Com::CVector2D po1,po2,po3,po4;
  {
    double min_x = ( pos.x < poe.x ) ? pos.x : poe.x;
    double max_x = ( pos.x > poe.x ) ? pos.x : poe.x;
    double min_y = ( pos.y < poe.y ) ? pos.y : poe.y;
    double max_y = ( pos.y > poe.y ) ? pos.y : poe.y;    
    po1 = Com::CVector2D(min_x,min_y);
    po2 = Com::CVector2D(max_x,min_y);
    po3 = Com::CVector2D(max_x,max_y);
    po4 = Com::CVector2D(min_x,max_y);    
  }
  if( !cad_2d.CheckIsPointInsideLoop(id_l,po1)  ) return 0;
  if( !cad_2d.CheckIsPointInsideLoop(id_l,po2)  ) return 0;       
  if( !cad_2d.CheckIsPointInsideLoop(id_l,po3)  ) return 0;
  if( !cad_2d.CheckIsPointInsideLoop(id_l,po4)  ) return 0;        
  std::vector< Com::CVector2D > aVec;
  aVec.push_back(po1);
  aVec.push_back(po2);
  aVec.push_back(po3);
  aVec.push_back(po4);
  //    const unsigned int id_l_add = cad_2d.AddPolygon(aVec,id_l).id_l_add;
  const Cad::CCadObj2D::CResAddPolygon& res = cad_2d.AddPolygon(aVec,id_l);
  const unsigned int id_l_add = res.id_l_add;
  if( id_l_add == 0 ){ return 0; }
  
  Solve_fromCad_InterpValue(); 	
  this->InitDrawer();     
  return id_l_add;    
}

bool CDesigner2D_Cloth::Cad_AddCutLine(const Com::CVector2D& pos, const Com::CVector2D& poe)
{
//  pAnalysis->SetClothPiecePlacingMode();
  unsigned int id_ls = 0, id_le = 0;
  {
    const std::vector<unsigned int>& aIdL = cad_2d.GetAryElemID(Cad::LOOP);  
    for(unsigned int iidl=0;iidl<aIdL.size();iidl++){
      unsigned int id_l = aIdL[iidl];
      if( cad_2d.CheckIsPointInsideLoop(id_l,pos) ){ id_ls = id_l; break; }
    }    
    for(unsigned int iidl=0;iidl<aIdL.size();iidl++){
      unsigned int id_l = aIdL[iidl];
      if( cad_2d.CheckIsPointInsideLoop(id_l,poe) ){ id_le = id_l; break; }
    }      
  }    
  ////
  if( id_ls > 0 && id_le > 0 && id_ls == id_le ){ // diamond dart
    return this->AddDartDiamond_Sym(id_ls,pos,poe);
  }
  if( id_ls > 0 && id_le == 0 )
  {
    if( !mesh_2d.IsIdLCad_CutMesh(id_ls) ) return 0;
    return this->AddDartEdge_Sym(id_ls,pos,poe);
  }
  if( id_le > 0 && id_ls == 0 )
  {
    if( !mesh_2d.IsIdLCad_CutMesh(id_le) ) return 0;    
    return this->AddDartEdge_Sym(id_le,poe,pos);      
  }  
  if( id_ls == 0 && id_le == 0 )
  {
    return false;
    unsigned int id_l_cut = 0;
    unsigned int id_e_near_s; bool is_same_dir_e_near_s;
    unsigned int id_e_near_e; bool is_same_dir_e_near_e;
    Com::CVector2D p_near_s, p_near_e;
    double ratio_s, ratio_e;
    const std::vector<unsigned int>& aIdL = cad_2d.GetAryElemID(Cad::LOOP);  
    for(unsigned int iidl=0;iidl<aIdL.size();iidl++){
      const unsigned int id_l = aIdL[iidl];            
      if( !mesh_2d.IsIdLCad_CutMesh(id_l) ) continue;
      bool bres_s = FindIntersectionEdge(cad_2d,id_l,pos,poe, id_e_near_s,is_same_dir_e_near_s, p_near_s);
      bool bres_e = FindIntersectionEdge(cad_2d,id_l,poe,pos, id_e_near_e,is_same_dir_e_near_e, p_near_e);
      if( !bres_s || !bres_e ){ continue; }
      if( Dot(poe-pos,p_near_e-pos) < 0 ) continue;
      if( Dot(pos-poe,p_near_s-poe) < 0 ) continue;      
      id_l_cut = id_l;
      {
        const Cad::CEdge2D& edge = cad_2d.GetEdge(id_e_near_s);
        const double len_s = Com::Distance(p_near_s,edge.po_s);
        const double len_e = Com::Distance(p_near_s,edge.po_e);
        const double invsum01 = 1.0/(len_s+len_e);
        ratio_s = len_s*invsum01;
      }      
      {
        const Cad::CEdge2D& edge = cad_2d.GetEdge(id_e_near_e);
        const double len_s = Com::Distance(p_near_e,edge.po_s);
        const double len_e = Com::Distance(p_near_e,edge.po_e);
        const double invsum01 = 1.0/(len_s+len_e);
        ratio_e = len_s*invsum01;
      }
      break;
    }
    if( id_l_cut == 0 ) return 0;
    const double eps_cad = cad_2d.GetMinClearance();
    assert( cad_2d.IsElemID(Cad::EDGE,id_e_near_s) );
    assert( cad_2d.IsElemID(Cad::EDGE,id_e_near_e) );
    Cad::CCadObj2D::CResAddVertex res0_s, res1_s;      
    Cad::CCadObj2D::CResAddVertex res0_e, res1_e; 
    {
      Com::CVector2D out0_s,out1_s;
      const Cad::CEdge2D& edge = cad_2d.GetEdge(id_e_near_s);
      bool is_exceed0, is_exceed1;
      edge.GetPointOnCurve_OnCircle(p_near_s, eps_cad, true, is_exceed0, out0_s);    
      edge.GetPointOnCurve_OnCircle(p_near_s, eps_cad,false, is_exceed1, out1_s);
      if( is_exceed0 || is_exceed1 ) return 0;
      res0_s = cad_2d.AddVertex(Cad::EDGE,id_e_near_s,out0_s);              
      res1_s = cad_2d.AddVertex(Cad::EDGE,id_e_near_s,out1_s);
    }      
    {
      Com::CVector2D out0_e,out1_e;
      const Cad::CEdge2D& edge = cad_2d.GetEdge(id_e_near_e);
      bool is_exceed0, is_exceed1;
      edge.GetPointOnCurve_OnCircle(p_near_e, eps_cad, true, is_exceed0, out0_e);    
      edge.GetPointOnCurve_OnCircle(p_near_e, eps_cad,false, is_exceed1, out1_e);
      if( is_exceed0 || is_exceed1 ) return 0;
      res0_e = cad_2d.AddVertex(Cad::EDGE,id_e_near_e,out0_e);
      res1_e = cad_2d.AddVertex(Cad::EDGE,id_e_near_e,out1_e);
    }    
    unsigned int id_l_add = 0;
    unsigned int id_e_add0, id_e_add1;
    if( is_same_dir_e_near_e != is_same_dir_e_near_s ){
      Cad::CBRepSurface::CResConnectVertex res_cv0 = cad_2d.ConnectVertex_Line(res0_s.id_v_add,res0_e.id_v_add);
      Cad::CBRepSurface::CResConnectVertex res_cv1 = cad_2d.ConnectVertex_Line(res1_s.id_v_add,res1_e.id_v_add);  
      cad_2d.RemoveElement(Cad::EDGE,res1_s.id_e_add);
      cad_2d.RemoveElement(Cad::EDGE,res1_e.id_e_add);
      if( is_same_dir_e_near_s ){ id_l_add = res_cv0.id_l_add; }
      else{                       id_l_add = res_cv1.id_l_add; }
      id_e_add0 = res_cv0.id_e_add;
      id_e_add1 = res_cv1.id_e_add;
    }
    else{
      Cad::CBRepSurface::CResConnectVertex res_cv0 = cad_2d.ConnectVertex_Line(res0_s.id_v_add,res1_e.id_v_add);                
      Cad::CBRepSurface::CResConnectVertex res_cv1 = cad_2d.ConnectVertex_Line(res1_s.id_v_add,res0_e.id_v_add);  
      cad_2d.RemoveElement(Cad::EDGE,res1_s.id_e_add);
      cad_2d.RemoveElement(Cad::EDGE,res1_e.id_e_add);      
      if( is_same_dir_e_near_s ){ id_l_add = res_cv0.id_l_add; }
      else{                       id_l_add = res_cv1.id_l_add; }
      id_e_add0 = res_cv0.id_e_add;
      id_e_add1 = res_cv1.id_e_add;      
    }    
    {
      Com::CVector2D dir;
      {
        dir = pos-poe; dir.Normalize();
        dir = Com::CVector2D(dir.y,-dir.x); dir*=-0.1;      
      }
    }
    {
      mesh_2d.AddIdLCad_CutMesh(id_l_add);       
      std::vector< std::pair<unsigned int,unsigned int> > aNewL;
      aNewL.push_back( std::make_pair(id_l_add,id_l_cut) );
      Solve_fromCad_InterpValue(aNewL);       
    }
    if( pAnalysis != 0 ){
      unsigned int id_e_add_opp_s = 0;
      unsigned int id_e_opp_s;  bool is_same_dir_opp_s;
      bool is_seam_line_cutted_s = false;
      if( pAnalysis->IsSeamLine(id_e_near_s,id_e_opp_s,is_same_dir_opp_s) ){
        is_seam_line_cutted_s = true;
        assert( cad_2d.IsElemID(Cad::EDGE,id_e_opp_s) );
        unsigned int id_v_opp0, id_v_opp1;
        cad_2d.GetIdVertex_Edge(id_v_opp0,id_v_opp1,id_e_opp_s);
        Com::CVector2D opp0 = cad_2d.GetVertexCoord(id_v_opp0);
        Com::CVector2D opp1 = cad_2d.GetVertexCoord(id_v_opp1);        
        Cad::CCadObj2D::CResAddVertex res1 = cad_2d.AddVertex(Cad::EDGE,id_e_opp_s,opp0*ratio_s+opp1*(1-ratio_s));
        id_e_add_opp_s = res1.id_e_add;
      }      
      ////
      unsigned int id_e_add_opp_e = 0;
      unsigned int id_e_opp_e;  bool is_same_dir_opp_e;
      bool is_seam_line_cutted_e = false;
      if( pAnalysis->IsSeamLine(id_e_near_e,id_e_opp_e,is_same_dir_opp_e) ){
        is_seam_line_cutted_e = true;
        assert( cad_2d.IsElemID(Cad::EDGE,id_e_opp_e) );
        unsigned int id_v_opp0, id_v_opp1;
        cad_2d.GetIdVertex_Edge(id_v_opp0,id_v_opp1,id_e_opp_e);
        Com::CVector2D opp0 = cad_2d.GetVertexCoord(id_v_opp0);
        Com::CVector2D opp1 = cad_2d.GetVertexCoord(id_v_opp1);        
        Cad::CCadObj2D::CResAddVertex res1 = cad_2d.AddVertex(Cad::EDGE,id_e_opp_e,opp0*ratio_e+opp1*(1-ratio_e));
        id_e_add_opp_e = res1.id_e_add;
      }            
      if( is_seam_line_cutted_s || is_seam_line_cutted_e ){ Solve_fromCad_InterpValue(); }
      if( is_seam_line_cutted_s ){
        pAnalysis->DisconnectEdge(id_e_near_s,id_e_opp_s);
        assert( mesh_2d.GetElemID_FromCadID(id_e_near_s,    Cad::EDGE) != 0 );
        assert( mesh_2d.GetElemID_FromCadID(id_e_add_opp_s, Cad::EDGE) != 0 );           
        assert( mesh_2d.GetElemID_FromCadID(res0_s.id_e_add,Cad::EDGE) != 0 );
        assert( mesh_2d.GetElemID_FromCadID(id_e_opp_s,     Cad::EDGE) != 0 );
        pAnalysis->ConnectEdge(this->cad_2d,this->mesh_2d,id_e_near_s,    id_e_add_opp_s); 
        pAnalysis->ConnectEdge(this->cad_2d,this->mesh_2d,res0_s.id_e_add,id_e_opp_s    );
      }
      if( is_seam_line_cutted_e ){
        pAnalysis->DisconnectEdge(id_e_near_e,id_e_opp_e);
        assert( mesh_2d.GetElemID_FromCadID(id_e_near_e,    Cad::EDGE) != 0 );
        assert( mesh_2d.GetElemID_FromCadID(id_e_add_opp_e, Cad::EDGE) != 0 );           
        assert( mesh_2d.GetElemID_FromCadID(res0_e.id_e_add,Cad::EDGE) != 0 );
        assert( mesh_2d.GetElemID_FromCadID(id_e_opp_e,     Cad::EDGE) != 0 );
        pAnalysis->ConnectEdge(this->cad_2d,this->mesh_2d,id_e_near_e,    id_e_add_opp_e); 
        pAnalysis->ConnectEdge(this->cad_2d,this->mesh_2d,res0_e.id_e_add,id_e_opp_e    );
      }      
      pAnalysis->ConnectEdge(this->cad_2d,this->mesh_2d,id_e_add0,id_e_add1);
    }
    
      
    this->InitDrawer();     
    return 0;    
  }
  
  return 0;
}

void CDesigner2D_Cloth::Msh_PrecompSlider(unsigned int islider)
{
  this->itype_cad_elem_prec = Cad::NOT_SET;
  this->id_cad_elem_prec = 0;
  slider_deform.GetLamVtx(islider,cad_2d,aLamXY_Vtx);
  mesh_2d.Precomp_FitMeshToCad(this->cad_2d,aLamXY_Vtx);
}

void SetSymTnsrIfSym
(unsigned int id_v, const std::vector< std::pair<unsigned int,unsigned int> >& aPair,
 std::vector<double>& aLamXY_Vtx)
{
  unsigned int id_vo = GetSymPairIdV(id_v,aPair);
  if( id_vo == 0 ) return;
  aLamXY_Vtx[id_vo*4+0] = -aLamXY_Vtx[id_v*4+0];
  aLamXY_Vtx[id_vo*4+1] = -aLamXY_Vtx[id_v*4+1];
  aLamXY_Vtx[id_vo*4+2] = +aLamXY_Vtx[id_v*4+2];
  aLamXY_Vtx[id_vo*4+3] = +aLamXY_Vtx[id_v*4+3];
}

void CDesigner2D_Cloth::Msh_PrecompDrag
(Cad::CAD_ELEM_TYPE itype_cad_elem, unsigned id_cad_elem)
{
  this->itype_cad_elem_prec = itype_cad_elem;
  this->id_cad_elem_prec = id_cad_elem;
  this->isnt_mesh_deform_sensitivity = false;
  {
    unsigned int max_id_v = 0;    
    const std::vector<unsigned int>& aIdV = cad_2d.GetAryElemID(Cad::VERTEX);
    for(unsigned int iiv=0;iiv<aIdV.size();iiv++){
      max_id_v = ( aIdV[iiv] > max_id_v ) ? aIdV[iiv] : max_id_v;
    }
    aLamXY_Vtx.clear();
    aLamXY_Vtx.resize((max_id_v+1)*4,0);    
  }
  if( itype_cad_elem_prec == Cad::VERTEX ){  
    const unsigned int id_v = id_cad_elem;
    bool flg_dart = false;
    for(unsigned int idart=0;idart<aDartOnEdge.size();idart++){
      const unsigned int id_v1 = aDartOnEdge[idart].id_v1;
      const unsigned int id_v2 = aDartOnEdge[idart].id_v2;      
      const unsigned int id_vc = aDartOnEdge[idart].id_vc;            
      if( id_v == id_v1 || id_v == id_v2 ){
        unsigned int id_l_in = aDartOnEdge[idart].id_l_in;
        unsigned int id_v1a = 0, id_v2a = 0;
        {
          for(Cad::CBRepSurface::CItrLoop itrl=cad_2d.GetItrLoop(id_l_in);!itrl.IsEnd();itrl++){
            if( itrl.GetIdVertex() == id_v1 ){
              if( itrl.GetIdVertex_Ahead() == id_vc ){ id_v1a = itrl.GetIdVertex_Behind(); }
              else{
                assert( id_vc == itrl.GetIdVertex_Behind() );
                id_v1a = itrl.GetIdVertex_Ahead();           
              }
            }
            if( itrl.GetIdVertex() == id_v2 ){
              if( itrl.GetIdVertex_Ahead() == id_vc ){ id_v2a = itrl.GetIdVertex_Behind(); }
              else{
                assert( id_vc == itrl.GetIdVertex_Behind() );
                id_v2a = itrl.GetIdVertex_Ahead();           
              }
            }            
          }          
        }
        const Com::CVector2D& vec1a = cad_2d.GetVertexCoord(id_v1a);      
        const Com::CVector2D& vec2a = cad_2d.GetVertexCoord(id_v2a);              
        const Com::CVector2D& vec1  = cad_2d.GetVertexCoord(id_v1);      
        const Com::CVector2D& vec2  = cad_2d.GetVertexCoord(id_v2);       
        Com::CVector2D evm1 = vec1-vec1a;  evm1.Normalize();        
        Com::CVector2D evm2 = vec2-vec2a;  evm2.Normalize();        
        flg_dart = true;        
        if( id_v == id_v1 ){
          aLamXY_Vtx[id_v1*4+0] = evm1.x*evm1.x;
          aLamXY_Vtx[id_v1*4+1] = evm1.y*evm1.x;
          aLamXY_Vtx[id_v1*4+2] = evm1.x*evm1.y;        
          aLamXY_Vtx[id_v1*4+3] = evm1.y*evm1.y;
          aLamXY_Vtx[id_v2*4+0] = evm1.x*evm2.x;
          aLamXY_Vtx[id_v2*4+1] = evm1.y*evm2.x;
          aLamXY_Vtx[id_v2*4+2] = evm1.x*evm2.y;
          aLamXY_Vtx[id_v2*4+3] = evm1.y*evm2.y;          
        }
        if( id_v == id_v2 ){        
          aLamXY_Vtx[id_v1*4+0] = evm1.x*evm2.x;
          aLamXY_Vtx[id_v1*4+1] = evm1.x*evm2.y;
          aLamXY_Vtx[id_v1*4+2] = evm1.y*evm2.x;
          aLamXY_Vtx[id_v1*4+3] = evm1.y*evm2.y;
          aLamXY_Vtx[id_v2*4+0] = evm2.x*evm2.x;
          aLamXY_Vtx[id_v2*4+1] = evm2.y*evm2.x;
          aLamXY_Vtx[id_v2*4+2] = evm2.x*evm2.y;        
          aLamXY_Vtx[id_v2*4+3] = evm2.y*evm2.y;
        }  
        SetSymTnsrIfSym(id_v1,aSymIdVPair,aLamXY_Vtx);
        SetSymTnsrIfSym(id_v2,aSymIdVPair,aLamXY_Vtx);        
      }    
    }    
    for(unsigned int idart=0;idart<aDartDiamond.size();idart++){
      const unsigned int id_vl = aDartDiamond[idart].id_vr;
      const unsigned int id_vr = aDartDiamond[idart].id_vl;
      const unsigned int id_vu = aDartDiamond[idart].id_vu;
      const unsigned int id_vd = aDartDiamond[idart].id_vd;            
      if( id_v == id_vl ){
        const Com::CVector2D& vecu = cad_2d.GetVertexCoord(id_vu);
        const Com::CVector2D& vecd = cad_2d.GetVertexCoord(id_vd);
        Com::CVector2D evm1 = vecu-vecd;  evm1.Normalize();        
        Com::CVector2D evm2(evm1.y,-evm1.x);
        aLamXY_Vtx[id_vl*4+0] = +1;
        aLamXY_Vtx[id_vl*4+3] = +1;
        aLamXY_Vtx[id_vr*4+0] = evm1.x*evm1.x - evm2.x*evm2.x;
        aLamXY_Vtx[id_vr*4+1] = evm1.y*evm1.x - evm2.y*evm2.x;
        aLamXY_Vtx[id_vr*4+2] = evm1.x*evm1.y - evm2.x*evm2.y;
        aLamXY_Vtx[id_vr*4+3] = evm1.y*evm1.y - evm2.y*evm2.y;
        flg_dart = true;
      }
      if( id_v == id_vr ){
        const Com::CVector2D& vecu = cad_2d.GetVertexCoord(id_vu);
        const Com::CVector2D& vecd = cad_2d.GetVertexCoord(id_vd);     
        Com::CVector2D evm1 = vecu-vecd;  evm1.Normalize();        
        Com::CVector2D evm2(evm1.y,-evm1.x);
        aLamXY_Vtx[id_vl*4+0] = evm1.x*evm1.x - evm2.x*evm2.x;
        aLamXY_Vtx[id_vl*4+1] = evm1.x*evm1.y - evm2.x*evm2.y;
        aLamXY_Vtx[id_vl*4+2] = evm1.y*evm1.x - evm2.y*evm2.x;
        aLamXY_Vtx[id_vl*4+3] = evm1.y*evm1.y - evm2.y*evm2.y;
        aLamXY_Vtx[id_vr*4+0] = +1;
        aLamXY_Vtx[id_vr*4+3] = +1;            
        flg_dart = true;        
      }    
      SetSymTnsrIfSym(id_vl,aSymIdVPair,aLamXY_Vtx);
      SetSymTnsrIfSym(id_vr,aSymIdVPair,aLamXY_Vtx);        
    }        
    if( !flg_dart ){
      for(unsigned int isym=0;isym<aSymIdVPair.size();isym++){
        unsigned int id_v1 = GetSymPairIdV(id_v,aSymIdVPair);        
        if( id_v1 == 0 ) continue;
        aLamXY_Vtx[id_v *4+0] = +1;
        aLamXY_Vtx[id_v *4+3] = +1;
        aLamXY_Vtx[id_v1*4+0] = -1;
        aLamXY_Vtx[id_v1*4+3] = +1;    
        flg_dart = true;
      }      
    }
    if( !flg_dart ){
      aLamXY_Vtx[id_v*4+0] = 1;
      aLamXY_Vtx[id_v*4+3] = 1;    
    }
  }
  if( itype_cad_elem_prec == Cad::EDGE ){
    const unsigned int id_e = id_cad_elem;
    bool flg_dart = false;
    for(unsigned int idart=0;idart<aDartOnEdge.size();idart++){
      if( id_e == aDartOnEdge[idart].id_e1 || id_e == aDartOnEdge[idart].id_e2 ){ 
        flg_dart = true;
        break;
      }
    }        
    unsigned int id_v1,id_v2;
    cad_2d.GetIdVertex_Edge(id_v1,id_v2,id_e);
    unsigned int id_v1o = GetSymPairIdV(id_v1,aSymIdVPair);
    unsigned int id_v2o = GetSymPairIdV(id_v2,aSymIdVPair);    
    if( !flg_dart ){
      if( id_v1o == id_v2 && id_v2o == id_v1 ){
        aLamXY_Vtx[id_v1*4+0] = 0;
        aLamXY_Vtx[id_v1*4+3] = 1;
        aLamXY_Vtx[id_v2*4+0] = 0;
        aLamXY_Vtx[id_v2*4+3] = 1;                
      }
      else if( id_v1o !=0 && id_v2o != 0 ){
        aLamXY_Vtx[id_v1 *4+0] = 1;
        aLamXY_Vtx[id_v1 *4+3] = 1;
        aLamXY_Vtx[id_v2 *4+0] = 1;
        aLamXY_Vtx[id_v2 *4+3] = 1;          
        aLamXY_Vtx[id_v1o*4+0] = -1;
        aLamXY_Vtx[id_v1o*4+3] = +1;
        aLamXY_Vtx[id_v2o*4+0] = -1;
        aLamXY_Vtx[id_v2o*4+3] = +1;
      }    
      else{
        aLamXY_Vtx[id_v1*4+0] = 1;
        aLamXY_Vtx[id_v1*4+3] = 1;
        aLamXY_Vtx[id_v2*4+0] = 1;
        aLamXY_Vtx[id_v2*4+3] = 1;          
      }      
    }
  }
  if( itype_cad_elem_prec == Cad::LOOP ){
    unsigned int id_l = id_cad_elem;
    bool isnt_dart = true;
    for(unsigned int idart=0;idart<aDartOnEdge.size();idart++){
      if( id_l == aDartOnEdge[idart].id_l ){
        unsigned int id_v1 = aDartOnEdge[idart].id_v1;
        unsigned int id_v2 = aDartOnEdge[idart].id_v2;                
        unsigned int id_vc = aDartOnEdge[idart].id_vc;           
        const Com::CVector2D& v1 = cad_2d.GetVertexCoord(id_v1);
        const Com::CVector2D& v2 = cad_2d.GetVertexCoord(id_v2);
        Com::CVector2D dir12 = v2-v1; dir12.Normalize();        
        aLamXY_Vtx[id_v1*4+0] = dir12.x*dir12.x;
        aLamXY_Vtx[id_v1*4+1] = dir12.x*dir12.y;
        aLamXY_Vtx[id_v1*4+2] = dir12.y*dir12.x;
        aLamXY_Vtx[id_v1*4+3] = dir12.y*dir12.y;
        aLamXY_Vtx[id_v2*4+0] = dir12.x*dir12.x;
        aLamXY_Vtx[id_v2*4+1] = dir12.x*dir12.y;
        aLamXY_Vtx[id_v2*4+2] = dir12.y*dir12.x;
        aLamXY_Vtx[id_v2*4+3] = dir12.y*dir12.y;
        aLamXY_Vtx[id_vc*4+0] = dir12.x*dir12.x;
        aLamXY_Vtx[id_vc*4+1] = dir12.x*dir12.y;
        aLamXY_Vtx[id_vc*4+2] = dir12.y*dir12.x;
        aLamXY_Vtx[id_vc*4+3] = dir12.y*dir12.y;
        SetSymTnsrIfSym(id_v1,aSymIdVPair,aLamXY_Vtx);
        SetSymTnsrIfSym(id_v2,aSymIdVPair,aLamXY_Vtx);        
        SetSymTnsrIfSym(id_vc,aSymIdVPair,aLamXY_Vtx);                        
        isnt_dart = false;        
      }
    }
    if( isnt_dart ){            
      for(unsigned int idart=0;idart<aDartDiamond.size();idart++){
        if( id_l == aDartDiamond[idart].id_l ){
          const unsigned int aIdV[4] = {
            aDartDiamond[idart].id_vu,
            aDartDiamond[idart].id_vr,
            aDartDiamond[idart].id_vd,
            aDartDiamond[idart].id_vl };
          for(unsigned int i=0;i<4;i++){
            unsigned int id_v0 = aIdV[i];
            aLamXY_Vtx[id_v0*4+0] = 1;
            aLamXY_Vtx[id_v0*4+1] = 0;
            aLamXY_Vtx[id_v0*4+2] = 0;
            aLamXY_Vtx[id_v0*4+3] = 1;
            SetSymTnsrIfSym(id_v0,aSymIdVPair,aLamXY_Vtx);
          }
          isnt_dart = false;        
        }
      }
    }
    if( isnt_dart ){
      for(Cad::CBRepSurface::CItrLoop itrl=cad_2d.GetItrLoop(id_l);!itrl.IsEnd();itrl.ShiftChildLoop()){
        for(itrl.Begin();!itrl.IsEnd();itrl++){        
          const unsigned int id_v = itrl.GetIdVertex();
          aLamXY_Vtx[id_v*4+0] = 1;
          aLamXY_Vtx[id_v*4+3] = 1;        
        }
      }
    }
    {
      /*
      bool iflg = true;
      for(unsigned int idart=0;idart<aDartDiamond.size();idart++){
        if( id_l == aDartDiamond[idart].id_l ){ iflg = false; }
      }
      for(unsigned int idart=0;idart<aDartOnEdge.size();idart++){
        if( id_l == aDartOnEdge[idart].id_l ){ iflg = false; }
      }        
      if( iflg ){ isnt_mesh_deform_sensitivity = true; }
       */
      isnt_mesh_deform_sensitivity = mesh_2d.IsIdLCad_CutMesh(id_l);
    }
    

  }
  mesh_2d.Precomp_FitMeshToCad(this->cad_2d,aLamXY_Vtx);
}

bool CDesigner2D_Cloth::MoveLoop(unsigned int id_l, const Com::CVector2D& dir)
{
  std::vector<unsigned int> aIdV_Move;
  
  {
    bool iflg = false;
    int idart0=-1;
    for(unsigned int idart=0;idart<aDartOnEdge.size();idart++){
      if( id_l == aDartOnEdge[idart].id_l ){
        idart0 = idart;
        break;
      }
    }
    if( idart0 != -1 ){ 
      aIdV_Move.push_back( aDartOnEdge[idart0].id_v1 );
      aIdV_Move.push_back( aDartOnEdge[idart0].id_v2 );
      aIdV_Move.push_back( aDartOnEdge[idart0].id_vc );      
      iflg = true;    
    }    
    if( !iflg ){
      for(unsigned int idart=0;idart<aDartDiamond.size();idart++){
        if( id_l == aDartDiamond[idart].id_l ){
          idart0 = idart;
          break;
        }
      }
      if( idart0 != -1 ){ 
        aIdV_Move.push_back( aDartDiamond[idart0].id_vu );
        aIdV_Move.push_back( aDartDiamond[idart0].id_vd );
        aIdV_Move.push_back( aDartDiamond[idart0].id_vr );      
        aIdV_Move.push_back( aDartDiamond[idart0].id_vl );            
        iflg = true;    
      }          
    }
    if( !iflg ){
      return cad_2d.MoveLoop(id_l,dir);
    }
  }      
  
  std::vector<unsigned int> aIdV_Add;
  for(unsigned int iiv=0;iiv<aIdV_Move.size();iiv++){
    unsigned int id_v = aIdV_Move[iiv];
    unsigned int id_vo = GetSymPairIdV(id_v,aSymIdVPair);
    if( id_vo == 0 ) continue;
    for(unsigned int jiv=0;jiv<aIdV_Move.size();jiv++){
      if( aIdV_Move[jiv] == id_vo ){ id_vo = 0; break; }
    }
    if( id_vo == 0 ) continue;
    aIdV_Add.push_back(id_vo);
  }
  for(unsigned int iiv=0;iiv<aIdV_Add.size();iiv++){ aIdV_Move.push_back(aIdV_Add[iiv]); }
  
  std::vector< std::pair<unsigned int,Com::CVector2D> > aIdDist;
  for(unsigned int iiv=0;iiv<aIdV_Move.size();iiv++){
    unsigned int id_v1 = aIdV_Move[iiv];
    const Com::CVector2D der1(aLamXY_Vtx[id_v1*4+0]*dir.x+aLamXY_Vtx[id_v1*4+1]*dir.y,
                              aLamXY_Vtx[id_v1*4+2]*dir.x+aLamXY_Vtx[id_v1*4+3]*dir.y);
    const Com::CVector2D vec1 = cad_2d.GetVertexCoord(id_v1);      
    aIdDist.push_back( std::make_pair(id_v1, vec1+der1) );    
  }
  return cad_2d.MoveVertex(aIdDist);
}

bool CDesigner2D_Cloth::MoveEdge(unsigned int id_e, const Com::CVector2D& dir)
{
  for(unsigned int idart=0;idart<aDartOnEdge.size();idart++){
    if( id_e == aDartOnEdge[idart].id_e1 || id_e == aDartOnEdge[idart].id_e2 ){ return false; }
  }
  unsigned int id_v1,id_v2;
  cad_2d.GetIdVertex_Edge(id_v1,id_v2,id_e);
  const unsigned int id_v1o = GetSymPairIdV(id_v1,aSymIdVPair);
  const unsigned int id_v2o = GetSymPairIdV(id_v2,aSymIdVPair);
  if( id_v1o == id_v2 && id_v2o == id_v1 ){
    return cad_2d.MoveEdge(id_e,Com::CVector2D(0,dir.y)); 
  }
  if( id_v1o == 0 && id_v2o == 0 ){
    return cad_2d.MoveEdge(id_e,dir); 
  }
  if( id_v1o !=0 && id_v2o != 0 ){
    const Com::CVector2D& vec1  = cad_2d.GetVertexCoord(id_v1 );
    const Com::CVector2D& vec2  = cad_2d.GetVertexCoord(id_v2 );
    const Com::CVector2D& vec1o = cad_2d.GetVertexCoord(id_v1o);
    const Com::CVector2D& vec2o = cad_2d.GetVertexCoord(id_v2o);
    std::vector< std::pair<unsigned int,Com::CVector2D> > aIdDist;                  
    aIdDist.push_back( std::make_pair(id_v1, vec1+dir) );
    aIdDist.push_back( std::make_pair(id_v2, vec2+dir) );
    Com::CVector2D diro = Com::CVector2D(-dir.x,dir.y);
    aIdDist.push_back( std::make_pair(id_v1o,vec1o+diro) );
    aIdDist.push_back( std::make_pair(id_v2o,vec2o+diro) );
    return cad_2d.MoveVertex(aIdDist);    
  }    
  return cad_2d.MoveEdge(id_e,dir); 
}

bool CDesigner2D_Cloth::MoveVertex(unsigned int id_v, const Com::CVector2D& mouse)
{
  std::vector<unsigned int> aIdV_Move;
  aIdV_Move.push_back(id_v);
  for(unsigned int idart=0;idart<aDartOnEdge.size();idart++){
    unsigned int id_v1 = aDartOnEdge[idart].id_v1;
    unsigned int id_v2 = aDartOnEdge[idart].id_v2;    
    if( id_v == id_v1 ){ aIdV_Move.push_back(id_v2); break; }
    if( id_v == id_v2 ){ aIdV_Move.push_back(id_v1); break; }
  }
  
  for(unsigned int idart=0;idart<aDartDiamond.size();idart++){
    const unsigned int id_vl = aDartDiamond[idart].id_vr;
    const unsigned int id_vr = aDartDiamond[idart].id_vl;
    if( id_v == id_vl ){ aIdV_Move.push_back(id_vr); break; }
    if( id_v == id_vr ){ aIdV_Move.push_back(id_vl); break; }
  }          
  
  std::vector<unsigned int> aIdV_Add;
  for(unsigned int iiv=0;iiv<aIdV_Move.size();iiv++){
    unsigned int id_v = aIdV_Move[iiv];
    unsigned int id_vo = GetSymPairIdV(id_v,aSymIdVPair);
    if( id_vo == 0 ) continue;
    for(unsigned int jiv=0;jiv<aIdV_Move.size();jiv++){
      if( aIdV_Move[jiv] == id_vo ){ id_vo = 0; break; }
    }
    if( id_vo == 0 ) continue;
    aIdV_Add.push_back(id_vo);
  }
  for(unsigned int iiv=0;iiv<aIdV_Add.size();iiv++){ aIdV_Move.push_back(aIdV_Add[iiv]); }
    if( aIdV_Move.size() > 1 ){
    std::vector< std::pair<unsigned int,Com::CVector2D> > aIdDist;
    const Com::CVector2D der0 = mouse-cad_2d.GetVertexCoord(id_v);
    for(unsigned int iiv=0;iiv<aIdV_Move.size();iiv++){
      unsigned int id_v1 = aIdV_Move[iiv];
      const Com::CVector2D der1(aLamXY_Vtx[id_v1*4+0]*der0.x+aLamXY_Vtx[id_v1*4+1]*der0.y,
                                aLamXY_Vtx[id_v1*4+2]*der0.x+aLamXY_Vtx[id_v1*4+3]*der0.y);
      const Com::CVector2D vec1 = cad_2d.GetVertexCoord(id_v1);      
      aIdDist.push_back( std::make_pair(id_v1, vec1+der1) );    
    }
    return cad_2d.MoveVertex(aIdDist);
  }  
  return cad_2d.MoveVertex(id_v,mouse);
}


void CDesigner2D_Cloth::SetValueSlider(unsigned int islider, double val)
{
  double min,max;
  const double val0 = slider_deform.GetValueSlider(islider, min,max);
  slider_deform.SetValueSlider(islider,val);
  const double val1 = slider_deform.GetValueSlider(islider, min,max);
  if( fabs(val0-val1) < 1.0e-20 ) return;
  if( !slider_deform.MoveCad(islider,val1, val0,cad_2d) ) return;
  is_updated_cad = true;
  if(  !m_is_solve_cad_change  ) return;
  unsigned int itype_ope = 0;    
  if( !mesh_2d.FitMeshToCad_Slider_UsingPrecomp(cad_2d,itype_ope,val1-val0) ){ goto FAIL_FIT_MESH; }  
  is_updated_coord = is_updated_coord || (itype_ope&1);
  is_updated_edge  = is_updated_edge  || (itype_ope&2);  
  return;
FAIL_FIT_MESH:
  this->Solve_fromCad_InterpValue();
//  this->Msh_PrecompDrag(itype_cad_part, id_cad_part);  
//  slider_deform.GetLamVtx(islider,cad_2d,aLamXY_Vtx);
  mesh_2d.Precomp_FitMeshToCad(this->cad_2d,aLamXY_Vtx);  
  is_updated_coord = false;
  is_updated_edge = false;
  is_updated_cad = false;  
}


void CDesigner2D_Cloth::LoadTimeStamp()
{
	this->m_id_cad_part = 0;
	this->m_itype_cad_part = Cad::NOT_SET;
  if( pAnalysis != 0 ){
    this->pAnalysis->LoadTimeStamp(this->cad_2d,this->mesh_2d);
  }
  this->InitDrawer();
}

void CDesigner2D_Cloth::SaveTimeStamp(){
  if( pAnalysis != 0 ){  
    this->pAnalysis->SaveTimeStamp(this->cad_2d,this->mesh_2d);  
  }
}




