
//#pragma comment(linker, "/subsystem:\"windows\" /entry:\"mainCRTStartup\"")

#if defined(__VISUALC__)
#pragma warning( disable : 4786 ) 
#endif
#define for if(0);else for


#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <cstdlib> //(exit)

#if defined(__APPLE__) && defined(__MACH__)
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif

#include "delfem/camera.h"
#include "delfem/drawer_gl_utility.h"
#include "delfem/glut_utility.h"

#include "delfem/serialize.h"
#include "delfem/cad_obj2d.h"
#include "delfem/drawer_cad.h"

Com::View::CCamera mvp_trans;
double mov_begin_x, mov_begin_y;
int imodifier;
Com::View::CDrawerArray drawer_ary;


void myGlutResize(int w, int h)
{
	mvp_trans.SetWindowAspect((double)w/h);
	glViewport(0, 0, w, h);
	::glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	Com::View::SetProjectionTransform(mvp_trans);
	glutPostRedisplay();
}

void myGlutDisplay(void)
{
//	::glClearColor(0.2, .7, 0.7, 1.0);
	::glClearColor(1.0, 1.0, 1.0, 1.0);
	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	::glEnable(GL_DEPTH_TEST);

	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 1.1f, 4.0f );

	::glMatrixMode(GL_MODELVIEW);
	::glLoadIdentity();
	Com::View::SetModelViewTransform(mvp_trans);

	drawer_ary.Draw();
	ShowFPS();

	glutSwapBuffers();
}

void myGlutIdle(){
	::glutPostRedisplay();
}

void myGlutMotion( int x, int y ){
	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
	const int win_w = viewport[2];
	const int win_h = viewport[3];
	const double mov_end_x = (2.0*x-win_w)/win_w;
	const double mov_end_y = (win_h-2.0*y)/win_h;
	if(      imodifier == GLUT_ACTIVE_CTRL ){
		mvp_trans.MouseRotation(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 
	}
	else if( imodifier == GLUT_ACTIVE_SHIFT ){
		mvp_trans.MousePan(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 
	}
	mov_begin_x = mov_end_x;
	mov_begin_y = mov_end_y;
	::glutPostRedisplay();
}

void myGlutMouse(int button, int state, int x, int y){
	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
	const int win_w = viewport[2];
	const int win_h = viewport[3];
	mov_begin_x = (2.0*x-win_w)/win_w;
	mov_begin_y = (win_h-2.0*y)/win_h;
  imodifier = ::glutGetModifiers();
	if( state == GLUT_DOWN && imodifier == 0 ){
		const unsigned int size_buffer = 2048;
		unsigned int selec_buffer[size_buffer];
		Com::View::PickPre(size_buffer,selec_buffer, x,y,5,5, mvp_trans);
		drawer_ary.DrawSelection();
		std::vector<Com::View::SSelectedObject> aSelecObj = Com::View::PickPost(selec_buffer, x,y, mvp_trans );
		drawer_ary.ClearSelected();
		if( aSelecObj.size() > 0 ){
			drawer_ary.AddSelected( aSelecObj[0].name );
		}
	}
}

bool SetNewProblem()
{
	const unsigned int nprob = 12;
	static unsigned int iprob = 0;
  
  std::cout << "SetNewProblem() " << iprob << std::endl;

	Cad::CCadObj2D cad_2d;
	if( iprob == 0 )
	{
		unsigned int id_l0;
 		{	// define shape
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
      vec_ary.push_back( Com::CVector2D(1.0,0.0) );
      vec_ary.push_back( Com::CVector2D(1.0,1.0) );
			vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			id_l0 = cad_2d.AddPolygon( vec_ary,0 ).id_l_add;
		}
		const unsigned int id_v5 = cad_2d.AddVertex(Cad::LOOP,id_l0,Com::CVector2D(0.5,0.5)).id_v_add;
		const unsigned int id_v6 = cad_2d.AddVertex(Cad::LOOP,id_l0,Com::CVector2D(0.5,0.8)).id_v_add;
		const unsigned int id_v7 = cad_2d.AddVertex(Cad::LOOP,id_l0,Com::CVector2D(0.8,0.5)).id_v_add;
		const unsigned int id_v8 = cad_2d.AddVertex(Cad::LOOP,id_l0,Com::CVector2D(0.5,0.2)).id_v_add;
		const unsigned int id_v9 = cad_2d.AddVertex(Cad::LOOP,id_l0,Com::CVector2D(0.2,0.5)).id_v_add;
		cad_2d.ConnectVertex_Line(id_v5,id_v6);	
		cad_2d.ConnectVertex_Line(id_v5,id_v7);	
		cad_2d.ConnectVertex_Line(id_v5,id_v8);	
		cad_2d.ConnectVertex_Line(id_v5,id_v9);
    {
      Cad::CCadObj2D cad_tmp(cad_2d);
      cad_2d.Clear();
      cad_2d = cad_tmp;
    }
		{ // export to the file
			Com::CSerializer fout("hoge.txt",false);
			cad_2d.Serialize(fout);
		}
		{ // import form the file
			Com::CSerializer fin( "hoge.txt",true);
			cad_2d.Serialize(fin);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new Cad::View::CDrawer_Cad2D(cad_2d) );
		drawer_ary.InitTrans(mvp_trans);
	}
	else if( iprob == 1 )
	{
    Cad::CCadObj2D::CResAddPolygon res;
 		{	// define shape
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.0) );
			vec_ary.push_back( Com::CVector2D(2.0,1.0) );
			vec_ary.push_back( Com::CVector2D(1.0,1.0) );
			vec_ary.push_back( Com::CVector2D(1.0,2.0) );
			vec_ary.push_back( Com::CVector2D(0.0,2.0) );
			res = cad_2d.AddPolygon( vec_ary );
		}
    cad_2d.RemoveElement(Cad::EDGE,res.aIdE[1]);    
    cad_2d.RemoveElement(Cad::EDGE,res.aIdE[2]);        
    cad_2d.RemoveElement(Cad::VERTEX,res.aIdV[2]);        
    unsigned int id_e_new = cad_2d.ConnectVertex_Line(res.aIdE[3],res.aIdE[1]).id_e_add;
    Cad::CCadObj2D::CResAddVertex res0 = cad_2d.AddVertex(Cad::EDGE,id_e_new,Com::CVector2D(1,0.4));
    Cad::CCadObj2D::CResAddVertex res1 = cad_2d.AddVertex(Cad::EDGE,id_e_new,Com::CVector2D(1,0.6));
    Cad::CCadObj2D::CResAddVertex res2 = cad_2d.AddVertex(Cad::EDGE,res.aIdE[5],Com::CVector2D(0,0.4));
    Cad::CCadObj2D::CResAddVertex res3 = cad_2d.AddVertex(Cad::EDGE,res.aIdE[5],Com::CVector2D(0,0.6));
    cad_2d.ConnectVertex_Line(res0.id_v_add,res2.id_v_add);    
    cad_2d.ConnectVertex_Line(res1.id_v_add,res3.id_v_add);
    cad_2d.RemoveElement(Cad::EDGE,res1.id_e_add);
    cad_2d.RemoveElement(Cad::EDGE,res3.id_e_add);
    { // copy test
      Cad::CCadObj2D cad_tmp(cad_2d);
      cad_2d.Clear();
      cad_2d = cad_tmp;
    }
		{	// export to the file
			Com::CSerializer fout("hoge.txt",false);
			cad_2d.Serialize(fout);
		}
		{	// import from the file
			Com::CSerializer fin( "hoge.txt",true);
			cad_2d.Serialize(fin);
		}
		drawer_ary.Clear();
		drawer_ary.PushBack( new Cad::View::CDrawer_Cad2D(cad_2d) );
		drawer_ary.InitTrans(mvp_trans);
	}
	else if( iprob == 2 )
	{
		unsigned int id_l0;
 		{	// Make model
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.0) );
      vec_ary.push_back( Com::CVector2D(1.0,0.4) );
      vec_ary.push_back( Com::CVector2D(0.0,0.4) );
			id_l0 = cad_2d.AddPolygon( vec_ary ).id_l_add;
		}
//		cad_2d.SetCurve_Arc(1,true,-0.2);
//		cad_2d.SetCurve_Arc(2,true, -0.5);
		cad_2d.SetCurve_Arc(3,false, 0);
//		cad_2d.SetCurve_Arc(4,true, -0.5);
    {
      Cad::CCadObj2D cad_tmp(cad_2d);
      cad_2d.Clear();
      cad_2d = cad_tmp;
    }
		drawer_ary.Clear();
		drawer_ary.PushBack( new Cad::View::CDrawer_Cad2D(cad_2d) );
		drawer_ary.InitTrans(mvp_trans);
	}
	else  if( iprob == 3 )
  {
		const unsigned int id_v1 = cad_2d.AddVertex(Cad::NOT_SET,0,Com::CVector2D(0,0)).id_v_add;
		const unsigned int id_v2 = cad_2d.AddVertex(Cad::NOT_SET,0,Com::CVector2D(1,0)).id_v_add;
		cad_2d.ConnectVertex_Line(id_v1,id_v2);
		const unsigned int id_v3 = cad_2d.AddVertex(Cad::NOT_SET,0,Com::CVector2D(1,1)).id_v_add;
		cad_2d.ConnectVertex_Line(id_v2,id_v3);
		const unsigned int id_v4 = cad_2d.AddVertex(Cad::NOT_SET,0,Com::CVector2D(0,1)).id_v_add;
		cad_2d.ConnectVertex_Line(id_v4,id_v3);
		const unsigned int id_v5 = cad_2d.AddVertex(Cad::EDGE,2,Com::CVector2D(1,0.5)).id_v_add;
		const unsigned int id_v7 = cad_2d.AddVertex(Cad::NOT_SET,0,Com::CVector2D(0.5,0.5)).id_v_add;
		cad_2d.ConnectVertex_Line(id_v5,id_v7);
		const unsigned int id_v6 = cad_2d.AddVertex(Cad::NOT_SET,0,Com::CVector2D(1.5,0.5)).id_v_add;
		cad_2d.ConnectVertex_Line(id_v5,id_v6);
		cad_2d.ConnectVertex_Line(id_v4,id_v1);
    {
      Cad::CCadObj2D cad_tmp(cad_2d);
      cad_2d.Clear();
      cad_2d = cad_tmp;
    }    
		drawer_ary.Clear();
		drawer_ary.PushBack( new Cad::View::CDrawer_Cad2D(cad_2d) );
		drawer_ary.InitTrans(mvp_trans);
	}
	else if( iprob == 4 )
	{
		unsigned int id_l0;
 		{	// define shape
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,1.0) );
			vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			id_l0 = cad_2d.AddPolygon( vec_ary ).id_l_add;
		}
//		const unsigned int id_v1 = cad_2d.AddVertex(Cad::LOOP,id_l0,Com::CVector2D(0.5,0.5));
		const unsigned int id_v2 = cad_2d.AddVertex(Cad::LOOP,id_l0,Com::CVector2D(0.5,0.2)).id_v_add;
		const unsigned int id_v3 = cad_2d.AddVertex(Cad::LOOP,id_l0,Com::CVector2D(0.8,0.5)).id_v_add;
		const unsigned int id_v4 = cad_2d.AddVertex(Cad::LOOP,id_l0,Com::CVector2D(0.5,0.8)).id_v_add;
		const unsigned int id_v5 = cad_2d.AddVertex(Cad::LOOP,id_l0,Com::CVector2D(0.2,0.5)).id_v_add;
		cad_2d.ConnectVertex_Line(id_v2,id_v3);
		cad_2d.ConnectVertex_Line(id_v3,id_v4);
		cad_2d.ConnectVertex_Line(id_v4,id_v5);
		unsigned int id_e1 = cad_2d.ConnectVertex_Line(id_v5,id_v2).id_e_add;
    cad_2d.RemoveElement(Cad::EDGE,id_e1);
    cad_2d.ConnectVertex_Line(id_v5,id_v2);    
    {
      Cad::CCadObj2D cad_tmp(cad_2d);
      cad_2d.Clear();
      cad_2d = cad_tmp;
    }    
		drawer_ary.Clear();
		drawer_ary.PushBack( new Cad::View::CDrawer_Cad2D(cad_2d) );
		drawer_ary.InitTrans(mvp_trans);
	}
	else if( iprob == 5 )
	{
		unsigned int id_l0;
 		{	// define shape
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			vec_ary.push_back( Com::CVector2D(1.0,1.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.0) );
			id_l0 = cad_2d.AddPolygon( vec_ary ).id_l_add;
		}
    { // copy test
      Cad::CCadObj2D cad_tmp(cad_2d);
      cad_2d.Clear();
      cad_2d = cad_tmp;
    }    
		drawer_ary.Clear();
		drawer_ary.PushBack( new Cad::View::CDrawer_Cad2D(cad_2d) );
		drawer_ary.InitTrans(mvp_trans);
	}
	else if( iprob == 6 )
	{
		unsigned int id_l0;
 		{	// define initial loop
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,1.0) );
			vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			id_l0 = cad_2d.AddPolygon( vec_ary ).id_l_add;
		}
		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.3,0.1) );
			vec_ary.push_back( Com::CVector2D(0.9,0.1) );
			vec_ary.push_back( Com::CVector2D(0.9,0.7) );
			cad_2d.AddPolygon( vec_ary, id_l0 );
		}
		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.1,0.9) );
			vec_ary.push_back( Com::CVector2D(0.7,0.9) );
			vec_ary.push_back( Com::CVector2D(0.1,0.3) );
			cad_2d.AddPolygon( vec_ary, id_l0 );
		}
		unsigned int id_e0 = cad_2d.ConnectVertex_Line(1,3).id_e_add;
		cad_2d.RemoveElement(Cad::EDGE,id_e0);
    { // copy test
      Cad::CCadObj2D cad_tmp(cad_2d);
      cad_2d.Clear();
      cad_2d = cad_tmp;
    }    
		drawer_ary.Clear();
		drawer_ary.PushBack( new Cad::View::CDrawer_Cad2D(cad_2d) );
		drawer_ary.InitTrans(mvp_trans);
	}
	else if( iprob == 7 )
	{
		unsigned int id_l0;
 		{	// define shape
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,1.0) );
			vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			id_l0 = cad_2d.AddPolygon( vec_ary ).id_l_add;
		}
		unsigned int id_v1 = cad_2d.AddVertex(Cad::LOOP,id_l0,Com::CVector2D(0.9,0.7) ).id_v_add;
		unsigned int id_v2 = cad_2d.AddVertex(Cad::LOOP,id_l0,Com::CVector2D(0.9,0.1) ).id_v_add;
		unsigned int id_v3 = cad_2d.AddVertex(Cad::LOOP,id_l0,Com::CVector2D(0.3,0.1) ).id_v_add;
		cad_2d.ConnectVertex_Line(id_v1,id_v2);
		cad_2d.ConnectVertex_Line(id_v2,id_v3);
		cad_2d.ConnectVertex_Line(id_v3,id_v1);
		{
			unsigned int id_e0 = cad_2d.ConnectVertex_Line(id_v2,2).id_e_add;
			cad_2d.RemoveElement(Cad::EDGE,id_e0);
		}
		{
			unsigned int id_e0 = cad_2d.ConnectVertex_Line(2,id_v2).id_e_add;
			cad_2d.RemoveElement(Cad::EDGE,id_e0);
		}
    { // copy test
      Cad::CCadObj2D cad_tmp(cad_2d);
      cad_2d.Clear();
      cad_2d = cad_tmp;
    }    
		drawer_ary.Clear();
		drawer_ary.PushBack( new Cad::View::CDrawer_Cad2D(cad_2d) );
		drawer_ary.InitTrans(mvp_trans);
	}
	else if( iprob == 8 ){
		Cad::CCadObj2D cad_2d;
		unsigned int id_e3, id_e4;
		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0, 0.0) );	// 1
			vec_ary.push_back( Com::CVector2D(1.0, 0.0) );	// 2
			vec_ary.push_back( Com::CVector2D(1.5, 0.0) );	// 3
			vec_ary.push_back( Com::CVector2D(2.0, 0.0) );	// 4
			vec_ary.push_back( Com::CVector2D(2.0, 1.0) );	// 5
			vec_ary.push_back( Com::CVector2D(1.5, 1.0) );	// 6
			vec_ary.push_back( Com::CVector2D(1.0, 1.0) );	// 7
			vec_ary.push_back( Com::CVector2D(0.0, 1.0) );	// 8
			unsigned int id_l0 = cad_2d.AddPolygon( vec_ary ).id_l_add;
			unsigned int id_v1 = cad_2d.AddVertex(Cad::LOOP,id_l0, Com::CVector2D(0.5,0.5) ).id_v_add;
			unsigned int id_e1 = cad_2d.ConnectVertex_Line(2,7).id_e_add;
			unsigned int id_e2 = cad_2d.ConnectVertex_Line(3,6).id_e_add;
			unsigned int id_v2 = cad_2d.AddVertex(Cad::EDGE,id_e1, Com::CVector2D(1.0,0.5) ).id_v_add;
			unsigned int id_v3 = cad_2d.AddVertex(Cad::EDGE,1,     Com::CVector2D(0.5,0.0) ).id_v_add;
			id_e3 = cad_2d.ConnectVertex_Line(id_v1,id_v2).id_e_add;
			id_e4 = cad_2d.ConnectVertex_Line(id_v1,id_v3).id_e_add;
		}
    { // copy test
      Cad::CCadObj2D cad_tmp(cad_2d);
      cad_2d.Clear();
      cad_2d = cad_tmp;
    }    
		drawer_ary.Clear();
		drawer_ary.PushBack( new Cad::View::CDrawer_Cad2D(cad_2d) );
		drawer_ary.InitTrans(mvp_trans);
	}
	else if( iprob == 9 )
	{
		unsigned int id_l0;
 		{	// define shape
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,1.0) );
			vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			id_l0 = cad_2d.AddPolygon( vec_ary ).id_l_add;
		}
		unsigned int id_v1 = cad_2d.AddVertex(Cad::LOOP,id_l0,Com::CVector2D(0.2,0.7) ).id_v_add;
		unsigned int id_v2 = cad_2d.AddVertex(Cad::LOOP,id_l0,Com::CVector2D(0.2,0.3) ).id_v_add;
		unsigned int id_e0 = cad_2d.ConnectVertex_Line(id_v1,id_v2).id_e_add;
		unsigned int id_v3 = cad_2d.AddVertex(Cad::EDGE,id_e0, Com::CVector2D(0.2,0.5) ).id_v_add;
		unsigned int id_v4 = cad_2d.AddVertex(Cad::LOOP,id_l0,Com::CVector2D(0.5,0.5) ).id_v_add;
		unsigned int id_v5 = cad_2d.AddVertex(Cad::LOOP,id_l0,Com::CVector2D(0.7,0.5) ).id_v_add;
		unsigned int id_e1 = cad_2d.ConnectVertex_Line(id_v3,id_v4).id_e_add;
		unsigned int id_e2 = cad_2d.ConnectVertex_Line(id_v4,id_v5).id_e_add;
		cad_2d.RemoveElement(Cad::EDGE,id_e1);
    { // copy test
      Cad::CCadObj2D cad_tmp(cad_2d);
      cad_2d.Clear();
      cad_2d = cad_tmp;
    }    
		drawer_ary.Clear();
		drawer_ary.PushBack( new Cad::View::CDrawer_Cad2D(cad_2d) );
		drawer_ary.InitTrans(mvp_trans);
	}
	else if( iprob == 10 )
	{
		unsigned int id_l0;
 		{	// define shape
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,1.0) );
			vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			id_l0 = cad_2d.AddPolygon( vec_ary ).id_l_add;
		}
		unsigned int id_v1 = cad_2d.AddVertex(Cad::LOOP,0,Com::CVector2D(1.1,0)).id_v_add;
		cad_2d.RemoveElement(Cad::VERTEX,id_v1);
    { // copy test
      Cad::CCadObj2D cad_tmp(cad_2d);
      cad_2d.Clear();
      cad_2d = cad_tmp;
    }    
		drawer_ary.Clear();
		drawer_ary.PushBack( new Cad::View::CDrawer_Cad2D(cad_2d) );
		drawer_ary.InitTrans(mvp_trans);
	}
	else if( iprob == 11 ){
		{	// define initial loop
      std::vector<Com::CVector2D> vec_ary;
      vec_ary.push_back( Com::CVector2D(0.0,0.0) );
      vec_ary.push_back( Com::CVector2D(1.0,0.0) );
      vec_ary.push_back( Com::CVector2D(1.0,1.0) );
      vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			cad_2d.AddPolygon( vec_ary );
		}
		{
			std::vector<Com::CVector2D> aRelCo;
			aRelCo.push_back( Com::CVector2D(0.25, -0.1) );
			aRelCo.push_back( Com::CVector2D( 0.5, -0.0) );
			aRelCo.push_back( Com::CVector2D(0.75, -0.1) );
			cad_2d.SetCurve_Polyline(1,aRelCo);
		}
		cad_2d.SetCurve_Arc(2,true,-1);
		{
			std::vector<Com::CVector2D> aRelCo;
			aRelCo.push_back( Com::CVector2D(+0.02, 0.75) );
			aRelCo.push_back( Com::CVector2D(-0.20, 0.50) );
			aRelCo.push_back( Com::CVector2D(+0.02, 0.25) );
			cad_2d.SetCurve_Polyline(4,aRelCo);
		}
		cad_2d.RemoveElement(Cad::VERTEX,1);
    { // copy shape
      Cad::CCadObj2D cad_tmp(cad_2d);
      cad_2d.Clear();
      cad_2d = cad_tmp;
    }    
		drawer_ary.Clear();
		drawer_ary.PushBack( new Cad::View::CDrawer_Cad2D(cad_2d) );
		drawer_ary.InitTrans(mvp_trans);
	}
	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
	Com::View::SetProjectionTransform(mvp_trans);

	iprob++;
	if( iprob == nprob ) iprob=0;

	return true;
}

void myGlutKeyboard(unsigned char key, int x, int y)
{
  switch (key) {
  case 'q':
  case 'Q':
  case '\033':  /* '\033' ÇÕ ESC ÇÃ ASCII ÉRÅ[Éh */
	  exit(0);
	  break;
  case ' ':
	  SetNewProblem();
	  break;
  default:
    break;
  }
}

int main(int argc,char* argv[])
{
	// Initailze GLUT
	::glutInitWindowPosition(200,200);
	::glutInitWindowSize(400, 300);
	::glutInit(&argc, argv);	
    ::glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
	::glutCreateWindow("Cad View");

	// Set callback function
	::glutMotionFunc(myGlutMotion);
	::glutMouseFunc(myGlutMouse);
	::glutDisplayFunc(myGlutDisplay);
	::glutReshapeFunc(myGlutResize);
	::glutKeyboardFunc(myGlutKeyboard);
	::glutIdleFunc(myGlutIdle);

	SetNewProblem();

	// Enter main loop
	::glutMainLoop();
	return 0;
}
