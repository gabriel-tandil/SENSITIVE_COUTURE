
//#pragma comment(linker, "/subsystem:\"windows\" /entry:\"mainCRTStartup\"")

#if defined(__VISUALC__)
#pragma warning( disable : 4786 ) 
#endif
#define for if(0);else for

#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include <math.h>
#include <cstdlib> //(exit)

#if defined(__APPLE__) && defined(__MACH__)
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif


#include "delfem/cad_obj2d.h"
#include "delfem/mesher2d.h"
#include "delfem/mesh3d.h"
#include "delfem/drawer_msh.h"

#include "delfem/camera.h"
#include "delfem/drawer_gl_utility.h"

Com::View::CCamera camera;
double mov_begin_x, mov_begin_y;
int imodifier;

void myGlutMotion( int x, int y ){
	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
	const int win_w = viewport[2];
	const int win_h = viewport[3];
	const double mov_end_x = (2.0*x-win_w)/win_w;
	const double mov_end_y = (win_h-2.0*y)/win_h;
  ////
  if(      imodifier == GLUT_ACTIVE_SHIFT ){ // pan
    camera.MousePan(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 		
	}
	else if( imodifier == GLUT_ACTIVE_CTRL ){ //  rotation
    camera.MouseRotation(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 		
	}
  ////
	mov_begin_x = mov_end_x;
	mov_begin_y = mov_end_y;
	::glutPostRedisplay();
}

void myGlutIdle(){
	::glutPostRedisplay();
}

void myGlutResize(int w, int h)
{
	camera.SetWindowAspect((double)w/h);

	::glViewport(0, 0, w, h);
	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
	Com::View::SetProjectionTransform(camera);
}

Com::View::CDrawerArray drawer_ary;

void myGlutSpecial(int Key, int x, int y)
{
	switch(Key)
	{
	case GLUT_KEY_PAGE_UP:
		if( ::glutGetModifiers() && GLUT_ACTIVE_SHIFT ){
			if( camera.IsPers() ){
				const double tmp_fov_y = camera.GetFovY() + 10.0;
				camera.SetFovY( tmp_fov_y );
			}
		}
		else{
			const double tmp_scale = camera.GetScale() * 0.9;
			camera.SetScale( tmp_scale );
		}
		break;
	case GLUT_KEY_PAGE_DOWN:
		if( ::glutGetModifiers() && GLUT_ACTIVE_SHIFT ){
			if( camera.IsPers() ){
				const double tmp_fov_y = camera.GetFovY() - 10.0;
				camera.SetFovY( tmp_fov_y );
			}
		}
		else{
			const double tmp_scale = camera.GetScale() * 1.111;
			camera.SetScale( tmp_scale );
		}
		break;
	case GLUT_KEY_HOME :
		drawer_ary.InitTrans(camera);
		camera.Fit();
		break;
	case GLUT_KEY_END :
		if( camera.IsPers() ) camera.SetIsPers(false);
		else{ camera.SetIsPers(true); }
		break;
	default:
		break;
	}
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	Com::View::SetProjectionTransform(camera);
	::glutPostRedisplay();
}

void myGlutDisplay(void)
{
	::glClearColor(1.0, 1.0, 1.0, 1.0);
	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	::glEnable(GL_DEPTH_TEST);

	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 1.1f, 4.0f );
   
	::glMatrixMode(GL_MODELVIEW);
	::glLoadIdentity();
	Com::View::SetModelViewTransform(camera);

	drawer_ary.Draw();
	::glutSwapBuffers();
}

void myGlutMouse(int button, int state, int x, int y){
	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
	const int win_w = viewport[2];
	const int win_h = viewport[3];
	mov_begin_x = (2.0*x-win_w)/(double)win_w;
	mov_begin_y = (win_h-2.0*y)/(double)win_h;
	imodifier = glutGetModifiers();
	if( button == GLUT_LEFT_BUTTON && state == GLUT_DOWN ){
		const unsigned int size_buffer = 128;
		unsigned int select_buffer[size_buffer];
		Com::View::PickPre(size_buffer,select_buffer,   x,y,  5,5,   camera);
		drawer_ary.DrawSelection();
		std::vector<Com::View::SSelectedObject> aSelecObj
			= Com::View::PickPost(select_buffer,   x,y,   camera );
		if( aSelecObj.size() > 0 ){ drawer_ary.AddSelected( aSelecObj[0].name ); }
		else{                       drawer_ary.ClearSelected(); }
	}
}

bool SetNewProblem()
{
	const unsigned int nprob = 14;
	static unsigned int iprob = 0;

	if( iprob == 0 )
	{
		Cad::CCadObj2D cad_2d;
		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.2) );
			vec_ary.push_back( Com::CVector2D(0.5,0.2) );
			vec_ary.push_back( Com::CVector2D(0.5,0.8) );
			vec_ary.push_back( Com::CVector2D(1.0,0.8) );
			vec_ary.push_back( Com::CVector2D(1.0,1.0) );
			vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			cad_2d.AddPolygon( vec_ary );
			cad_2d.ConnectVertex_Line(6,3);
		}
		Msh::CMesher2D mesh_2d(cad_2d,0.05);
    ////
    {
      Msh::CMesher2D msh_tmp(mesh_2d);
      mesh_2d.Clear();
      mesh_2d = msh_tmp;
    }
		{	// write file
			Com::CSerializer fout("hoge.txt",false);
			mesh_2d.Serialize(fout);
		}
		{	// load file
			Com::CSerializer fin( "hoge.txt",true);
			mesh_2d.Serialize(fin);
		}
    ////
		drawer_ary.Clear();
		drawer_ary.PushBack( new Msh::View::CDrawerMsh2D(mesh_2d) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 1 )
	{
		Msh::CMesher2D mesh_2d;
		mesh_2d.ReadFromFile_GiDMsh("../input_file/hexa_tri.msh");
    {
      Msh::CMesher2D msh_tmp(mesh_2d);
      mesh_2d.Clear();
      mesh_2d = msh_tmp;
    }    
		{	// write file
			Com::CSerializer fout("hoge.txt",false);
			mesh_2d.Serialize(fout);
		}
		{	// read file
			Com::CSerializer fin( "hoge.txt",true);
			mesh_2d.Serialize(fin);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new Msh::View::CDrawerMsh2D(mesh_2d) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 2 )
	{
		Msh::CMesher2D mesh_2d;
		mesh_2d.ReadFromFile_GiDMsh("../input_file/rect_quad.msh");
    {
      Msh::CMesher2D msh_tmp(mesh_2d);
      mesh_2d.Clear();
      mesh_2d = msh_tmp;
    }    
		{	// write file
			Com::CSerializer fout("hoge.txt",false);
			mesh_2d.Serialize(fout);
		}
		{	// read file
			Com::CSerializer fin( "hoge.txt",true);
			mesh_2d.Serialize(fin);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new Msh::View::CDrawerMsh2D(mesh_2d) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 3 )
	{
		Cad::CCadObj2D cad_2d;
		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.5) );
			vec_ary.push_back( Com::CVector2D(0.5,0.5) );
			vec_ary.push_back( Com::CVector2D(0.5,1.0) );
			vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			cad_2d.AddPolygon( vec_ary );
		}
		Msh::CMesher2D mesh_2d(cad_2d, 0.1);
		Msh::CMesh3D_Extrude mesh_3d;
		mesh_3d.Extrude(mesh_2d, 0.5, 0.1 );
    {
      Msh::CMesher2D msh_tmp(mesh_2d);
      mesh_2d.Clear();
      mesh_2d = msh_tmp;
    }    
		{	// write file
			Com::CSerializer fout("hoge.txt",false);
			mesh_3d.Serialize(fout);
		}
		{	// read file
			Com::CSerializer fin( "hoge.txt",true);
			mesh_3d.Serialize(fin);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new Msh::View::CDrawerMsh3D(mesh_3d) );
		drawer_ary.InitTrans(camera);
	}
	else if( iprob == 4 )	// load mesh of GiD
	{
		Msh::CMesher2D mesh_2d;
		mesh_2d.ReadFromFile_GiDMsh("../input_file/hexa_tri.msh");
		Msh::CMesh3D_Extrude mesh_3d;
		mesh_3d.Extrude(mesh_2d, 5.0, 0.5 );
    {
      Msh::CMesher2D msh_tmp(mesh_2d);
      mesh_2d.Clear();
      mesh_2d = msh_tmp;
    }    
		{	// write file
			Com::CSerializer fout("hoge.txt",false);
			mesh_2d.Serialize(fout);
		}
		{	// read file
			Com::CSerializer fin( "hoge.txt",true);
			mesh_2d.Serialize(fin);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new Msh::View::CDrawerMsh3D(mesh_3d) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 5 )
	{
		Msh::CMesher3D mesh_3d;
		mesh_3d.ReadFromFile_GiDMsh("../input_file/cylinder_hex.msh");
		{	// write file
			Com::CSerializer fout("hoge.txt",false);
			mesh_3d.Serialize(fout);
		}
		{	// read files
			Com::CSerializer fin( "hoge.txt",true);
			mesh_3d.Serialize(fin);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new Msh::View::CDrawerMsh3D(mesh_3d) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 6 )
	{
		Msh::CMesher3D mesh_3d;
		mesh_3d.ReadFromFile_GiDMsh("../input_file/cylinder_tet.msh");
		{	// write file
			Com::CSerializer fout("hoge.txt",false);
			mesh_3d.Serialize(fout);
		}
		{	// read file
			Com::CSerializer fin( "hoge.txt",true);
			mesh_3d.Serialize(fin);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new Msh::View::CDrawerMsh3D(mesh_3d) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 7 )
	{
		Msh::CMesher2D mesh_2d;
		mesh_2d.ReadFromFile_GiDMsh("../input_file/rect_quad.msh");
		Msh::CMesh3D_Extrude mesh_3d;
		mesh_3d.Extrude(mesh_2d, 5.0, 0.5 );
    {
      Msh::CMesher2D msh_tmp(mesh_2d);
      mesh_2d.Clear();
      mesh_2d = msh_tmp;
    }    
		{	// write file
			Com::CSerializer fout("hoge.txt",false);
			mesh_2d.Serialize(fout);
		}
		{	// read file
			Com::CSerializer fin( "hoge.txt",true);
			mesh_2d.Serialize(fin);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new Msh::View::CDrawerMsh3D(mesh_3d) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 8 )
	{
		Cad::CCadObj2D cad_2d;
		unsigned int id_l = 0;
		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,1.0) );
			vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			id_l = cad_2d.AddPolygon( vec_ary ).id_l_add;
		}
		cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.8,0.6) );
		cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.6,0.6) );
		cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.4,0.6) );
		cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.2,0.6) );
		cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.8,0.4) );
		cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.6,0.4) );
		cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.4,0.4) );
		cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.2,0.4) );
		Msh::CMesher2D mesh_2d(cad_2d,0.02);
    {
      Msh::CMesher2D msh_tmp(mesh_2d);
      mesh_2d.Clear();
      mesh_2d = msh_tmp;
    }    
		{	// write file
			Com::CSerializer fout("hoge.txt",false);
			mesh_2d.Serialize(fout);
		}
		{	// read file
			Com::CSerializer fin( "hoge.txt",true);
			mesh_2d.Serialize(fin);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new Msh::View::CDrawerMsh2D(mesh_2d) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 9 )
	{
		Cad::CCadObj2D cad_2d;
		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
      vec_ary.push_back( Com::CVector2D(1.0,0.0) );
      vec_ary.push_back( Com::CVector2D(1.0,1.0) );
      vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			cad_2d.AddPolygon( vec_ary );
		}
		cad_2d.SetCurve_Arc(1,true, -0.5);
		cad_2d.SetCurve_Arc(2,false,-0.5);
		cad_2d.SetCurve_Arc(3,true, -0.5);
		cad_2d.SetCurve_Arc(4,false,-0.5);
		Msh::CMesher2D mesh_2d(cad_2d,0.05);
    {
      Msh::CMesher2D msh_tmp(mesh_2d);
      mesh_2d.Clear(); 
      mesh_2d = msh_tmp;
    }    
		{	// write file
			Com::CSerializer fout("hoge.txt",false);
			mesh_2d.Serialize(fout);
		}
		{	// read file
			Com::CSerializer fin( "hoge.txt",true);
			mesh_2d.Serialize(fin);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new Msh::View::CDrawerMsh2D(mesh_2d) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 10 )
	{
		Cad::CCadObj2D cad_2d;
		{	// define shape
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,1.0) );
			vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			const unsigned int id_l = cad_2d.AddPolygon( vec_ary ).id_l_add;
			const unsigned int id_v1 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.3,0.2)).id_v_add;
			const unsigned int id_v2 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.7,0.2)).id_v_add;
			const unsigned int id_v3 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.7,0.8)).id_v_add;
			const unsigned int id_v4 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.3,0.8)).id_v_add;
			cad_2d.ConnectVertex_Line(id_v1,id_v2);
			cad_2d.ConnectVertex_Line(id_v2,id_v3);
			cad_2d.ConnectVertex_Line(id_v3,id_v4);
			cad_2d.ConnectVertex_Line(id_v4,id_v1);
		}
		Msh::CMesher2D mesh_2d;
		{
		  mesh_2d.AddIdLCad_CutMesh(1);	// cut mesh to loop whitch have id 1
		  mesh_2d.SetMeshingMode_ElemLength(0.05);
		  mesh_2d.Meshing(cad_2d);
		}
    {
      Msh::CMesher2D msh_tmp(mesh_2d);
      mesh_2d.Clear();
      mesh_2d = msh_tmp;
    }    
		{	// write file
			Com::CSerializer fout("hoge.txt",false);
			mesh_2d.Serialize(fout);
		}
		{	// read file
			Com::CSerializer fin( "hoge.txt",true);
			mesh_2d.Serialize(fin);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new Msh::View::CDrawerMsh2D(mesh_2d) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 11 )	// mesh with cut
	{
		Cad::CCadObj2D cad_2d;
		unsigned int id_l;
		unsigned int id_e1, id_e2, id_e3, id_e4, id_e5;
		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(0.3,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,1.0) );
			vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			id_l = cad_2d.AddPolygon( vec_ary ).id_l_add;
			unsigned int id_v1 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(0.3,0.5) ).id_v_add;
			id_e1 = cad_2d.ConnectVertex_Line(2,id_v1).id_e_add;
			unsigned int id_v2 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(0.7,0.5) ).id_v_add;
			unsigned int id_v3 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(0.7,0.2) ).id_v_add;
			unsigned int id_v4 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(0.7,0.8) ).id_v_add;
			unsigned int id_v5 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(0.5,0.5) ).id_v_add;
			unsigned int id_v6 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(0.9,0.5) ).id_v_add;
			id_e2 = cad_2d.ConnectVertex_Line(id_v2,id_v3).id_e_add;
			id_e3 = cad_2d.ConnectVertex_Line(id_v2,id_v4).id_e_add;
			id_e4 = cad_2d.ConnectVertex_Line(id_v2,id_v5).id_e_add;
			id_e5 = cad_2d.ConnectVertex_Line(id_v2,id_v6).id_e_add;
		}
		Msh::CMesher2D mesh_2d(cad_2d,0.2);
//		mesh_2d.Tesselation(cad_2d);

		{
			std::vector<unsigned int> aIdMsh_Inc;
			aIdMsh_Inc.push_back( mesh_2d.GetElemID_FromCadID(id_l,Cad::LOOP) );
			std::vector<unsigned int> aIdMshBar_Cut;
			aIdMshBar_Cut.push_back( mesh_2d.GetElemID_FromCadID(id_e1,Cad::EDGE) );
			aIdMshBar_Cut.push_back( mesh_2d.GetElemID_FromCadID(id_e2,Cad::EDGE) );
			aIdMshBar_Cut.push_back( mesh_2d.GetElemID_FromCadID(id_e3,Cad::EDGE) );
			aIdMshBar_Cut.push_back( mesh_2d.GetElemID_FromCadID(id_e4,Cad::EDGE) );
			aIdMshBar_Cut.push_back( mesh_2d.GetElemID_FromCadID(id_e5,Cad::EDGE) );
			////////////////
			std::vector< std::vector<int> > aLnods;
			std::vector<unsigned int> mapVal2Co;
			mesh_2d.GetClipedMesh(aLnods,mapVal2Co,   aIdMsh_Inc,aIdMshBar_Cut);
		}
    {
      Msh::CMesher2D msh_tmp(mesh_2d);
      mesh_2d.Clear();
      mesh_2d = msh_tmp;
    }    
		drawer_ary.Clear();
		drawer_ary.PushBack( new Msh::View::CDrawerMsh2D(mesh_2d) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 12 )	// mesh with cut
	{
		Cad::CCadObj2D cad_2d;
		unsigned int id_l1, id_l2;
		unsigned int id_e3, id_e4;
		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0, 0.0) );
			vec_ary.push_back( Com::CVector2D(1.0, 0.0) );
			vec_ary.push_back( Com::CVector2D(1.5, 0.0) );
			vec_ary.push_back( Com::CVector2D(2.0, 0.0) );
			vec_ary.push_back( Com::CVector2D(2.0, 1.0) );
			vec_ary.push_back( Com::CVector2D(1.5, 1.0) );
			vec_ary.push_back( Com::CVector2D(1.0, 1.0) );
			vec_ary.push_back( Com::CVector2D(0.0, 1.0) );
			unsigned int id_l0 = cad_2d.AddPolygon( vec_ary ).id_l_add;
			unsigned int id_v1 = cad_2d.AddVertex(Cad::LOOP,id_l0, Com::CVector2D(0.5,0.5) ).id_v_add;
			unsigned int id_e1 = cad_2d.ConnectVertex_Line(2,7).id_e_add;
			unsigned int id_e2 = cad_2d.ConnectVertex_Line(3,6).id_e_add;
			unsigned int id_v2 = cad_2d.AddVertex(Cad::EDGE,id_e1, Com::CVector2D(1.0,0.5) ).id_v_add;
			unsigned int id_v3 = cad_2d.AddVertex(Cad::EDGE,1,     Com::CVector2D(0.5,0.0) ).id_v_add;
			id_e3 = cad_2d.ConnectVertex_Line(id_v1,id_v2).id_e_add;
			id_e4 = cad_2d.ConnectVertex_Line(id_v1,id_v3).id_e_add;
			id_l1 = 1;
			id_l2 = 2;
		}
		Msh::CMesher2D mesh_2d(cad_2d,0.2);
//		mesh_2d.Tesselation(cad_2d);
		{
			std::vector<unsigned int> aIdMsh_Inc;
			aIdMsh_Inc.push_back( mesh_2d.GetElemID_FromCadID(id_l1,Cad::LOOP) );
			aIdMsh_Inc.push_back( mesh_2d.GetElemID_FromCadID(id_l2,Cad::LOOP) );
			std::vector<unsigned int> aIdMshBar_Cut;
//			aIdMshBar_Cut.push_back( mesh_2d.GetElemID_FromCadID(id_e3,Cad::EDGE) );
			aIdMshBar_Cut.push_back( mesh_2d.GetElemID_FromCadID(id_e4,Cad::EDGE) );
			////////////////
			std::vector< std::vector<int> > aLnods;
			std::vector<unsigned int> mapVal2Co;
			mesh_2d.GetClipedMesh(aLnods,mapVal2Co,   aIdMsh_Inc,aIdMshBar_Cut);
		}
    {
      Msh::CMesher2D msh_tmp(mesh_2d);
      mesh_2d.Clear();
      mesh_2d = msh_tmp;
    }    
		drawer_ary.Clear();
		drawer_ary.PushBack( new Msh::View::CDrawerMsh2D(mesh_2d) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 13 ){
		Cad::CCadObj2D cad_2d;
		unsigned int id_e3, id_e4;
		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0, 0.0) );	// 1
			vec_ary.push_back( Com::CVector2D(1.5, 0.0) );	// 2
			vec_ary.push_back( Com::CVector2D(1.5, 0.4) );	// 3
			vec_ary.push_back( Com::CVector2D(1.0, 0.4) );	// 4
			vec_ary.push_back( Com::CVector2D(1.0, 0.5) );	// 5
			vec_ary.push_back( Com::CVector2D(2.0, 0.5) );	// 6
			vec_ary.push_back( Com::CVector2D(2.0, 1.0) );	// 7
			vec_ary.push_back( Com::CVector2D(0.0, 1.0) );	// 8
			vec_ary.push_back( Com::CVector2D(0.0, 0.5) );	// 9
			unsigned int id_l0 = cad_2d.AddPolygon( vec_ary ).id_l_add;
			unsigned int id_e1 = cad_2d.ConnectVertex_Line(5,9).id_e_add;
			cad_2d.ShiftLayer_Loop(id_l0,true);
			const double col[3] = { 0.9, 0.4, 0.4 };
			cad_2d.SetColor_Loop(id_l0, col);
			cad_2d.AddVertex(Cad::EDGE,3, Com::CVector2D(1.3,0.5) );
		}
		Msh::CMesher2D mesh_2d(cad_2d,0.05);
    {
      Msh::CMesher2D msh_tmp(mesh_2d);
      mesh_2d.Clear();
      mesh_2d = msh_tmp;
    }    
		drawer_ary.Clear();
		drawer_ary.PushBack( new Msh::View::CDrawerMsh2D(mesh_2d) );
		drawer_ary.InitTrans( camera );
	}

	iprob++;
	if( iprob == nprob ) iprob = 0;

	return true;
}

void myGlutKeyboard(unsigned char key, int x, int y)
{
  switch (key) {
  case 'q':
  case 'Q':
  case '\033': // ´033 is esc
	  exit(0);
	  break;
  case ' ':
	  SetNewProblem();
	  ::glMatrixMode(GL_PROJECTION);
	  ::glLoadIdentity();
	  Com::View::SetProjectionTransform(camera);
	  break;
  default:
    break;
  }
}

int main(int argc,char* argv[])
{
	::glutInitWindowPosition(200,200);
	::glutInitWindowSize(400, 400);
	::glutInit(&argc, argv);
	::glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
	::glutCreateWindow("DelFEM demo");

    // Set call back function
	::glutMotionFunc(myGlutMotion);
	::glutMouseFunc(myGlutMouse);
	::glutKeyboardFunc(myGlutKeyboard);
	::glutSpecialFunc(myGlutSpecial);
	::glutDisplayFunc(myGlutDisplay);
	::glutReshapeFunc(myGlutResize);
	::glutIdleFunc(myGlutIdle);
  
  camera.SetRotationMode( Com::View::ROT_3D );
	
	SetNewProblem();

    // Enter Main Loop
	::glutMainLoop();
	return 0;
}
