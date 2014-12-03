////////////////////////////////////////////////////////////////
//                                                            //
//		DelFEM Test_glut of Fluid 2D                            //
//                                                            //
//          Copy Rights (c) Nobuyuki Umetani 2009             //
//          e-mail : numetani@gmail.com                       //
////////////////////////////////////////////////////////////////

//#pragma comment(linker, "/subsystem:\"windows\" /entry:\"mainCRTStartup\"")

#if defined(__VISUALC__)
#pragma warning( disable : 4786 ) 
#endif
#define for if(0);else for

#include <iostream>
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

#include "delfem/camera.h"	// camera parameter class
#include "delfem/glut_utility.h"  // FPS()

#include "delfem/cad_obj2d.h"	// Cad::CCadObj2D
#include "delfem/mesher2d.h"	// Mesh::CMesher2D

#include "delfem/field.h" // Fem::Field::CField
#include "delfem/field_value_setter.h"
#include "delfem/field_world.h"	// Fem::Field::CFieldWorld
#include "delfem/drawer_field_face.h"
#include "delfem/drawer_field_edge.h"
#include "delfem/drawer_field_vector.h"
#include "delfem/drawer_field_image_based_flow_vis.h"
#include "delfem/drawer_field_streamline.h"

#include "delfem/eqnsys_fluid.h"		// Fem::Eqn::CEqnSystem_Fluid2D

void myGlutIdle(){
	glutPostRedisplay();
}

Com::View::CCamera camera;
double mov_begin_x, mov_begin_y;
bool is_animation = true;
int imodifier;

void myGlutResize(int w, int h)
{
	camera.SetWindowAspect((double)w/h);
	glViewport(0, 0, w, h);
	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
	Com::View::SetProjectionTransform(camera);
	glutPostRedisplay();
}

Fem::Field::View::CDrawerArrayField drawer_ary;
double cur_time = 0.0;
double dt = 0.8;
unsigned int id_base;
Fem::Field::CFieldWorld world;
Fem::Eqn::CEqnSystem_Fluid2D fluid;
Fem::Field::CFieldValueSetter field_value_setter;

void myGlutMotion( int x, int y ){
	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
	const int win_w = viewport[2];
	const int win_h = viewport[3];
	const double mov_end_x = (2.0*x-win_w)/win_w;
	const double mov_end_y = (win_h-2.0*y)/win_h;
  if(      imodifier == GLUT_ACTIVE_SHIFT ){
    camera.MousePan(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 
  }
  else if( imodifier == GLUT_ACTIVE_CTRL ){
    camera.MouseRotation(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 
  }
	mov_begin_x = mov_end_x;
	mov_begin_y = mov_end_y;
	::glutPostRedisplay();
}

void myGlutMouse(int button, int state, int x, int y){
	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
  imodifier = ::glutGetModifiers();
	const int win_w = viewport[2];
	const int win_h = viewport[3];
	mov_begin_x = (2.0*x-win_w)/win_w;
	mov_begin_y = (win_h-2.0*y)/win_h;
}

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
		camera.Fit();
		break;
	case GLUT_KEY_END :
		if( camera.IsPers() ) camera.SetIsPers(false);
		else{ camera.SetIsPers(true); }
		break;
	default:
		break;
	}
	
	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
	Com::View::SetProjectionTransform(camera);
	::glutPostRedisplay();
}

void myGlutDisplay(void)
{
//	::glClearColor(0.2, .7, .7 ,1.0);
	::glClearColor(1, 1, 1 ,1.0);
	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	::glEnable(GL_DEPTH_TEST);

	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 1.1f, 4.0f );

	::glMatrixMode(GL_MODELVIEW);
	::glLoadIdentity();
	Com::View::SetModelViewTransform(camera);
	
	if( is_animation )
	{
		cur_time += dt;
    field_value_setter.ExecuteValue(cur_time,world);
		fluid.Solve(world);
		if( fluid.GetAry_ItrNormRes().size() > 0 ){
//			std::cout << "Iter : " << fluid.GetAry_ItrNormRes()[0].first << " ";
//			std::cout << "Res : " << fluid.GetAry_ItrNormRes()[0].second << std::endl;
		}
		drawer_ary.Update(world);
	}
	drawer_ary.Draw();
	ShowFPS();
	glutSwapBuffers();
}


void SetNewProblem();
void myGlutKeyboard(unsigned char Key, int x, int y)
{
	switch(Key)
	{
	case 'q':
	case 'Q':
	case '\033':
		exit(0);  /* '\033' ? ESC ? ASCII ??? */
	case 'a':
		if( is_animation ){ is_animation = false; }
		else{ is_animation = true; }
		break;
	case ' ':
		SetNewProblem();
		break;
	}
	::glutPostRedisplay();
}

    
void SetNewProblem()
{
	const unsigned int nprob = 15;	// number of problems
	static unsigned int iprob = 0;

	if( iprob == 0 )	// cavity flow (stationaly stokes)
	{
		Cad::CCadObj2D cad_2d;
 		{	// define shape
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(-0.5,-0.5) );
			vec_ary.push_back( Com::CVector2D( 0.5,-0.5) );
			vec_ary.push_back( Com::CVector2D( 0.5, 0.5) );
			vec_ary.push_back( Com::CVector2D(-0.5, 0.5) );
			unsigned int id_l = cad_2d.AddPolygon( vec_ary ).id_l_add;
			cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.0,0.0));
		}
		cur_time = 0;
		world.Clear();
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.04) );
		const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);
		fluid.Clear();
		fluid.UnSetInterpolationBubble();
		fluid.UpdateDomain_Field(id_base,world);

		unsigned int id_field_press = fluid.GetIdField_Press();
//		unsigned int id_field_press_bc0 = world.GetPartialField(id_field_press,10);
//		fluid.AddFixField(id_field_press_bc0,world);

		unsigned int id_field_velo  = fluid.GetIdField_Velo();
		unsigned int id_field_bc0 = fluid.AddFixElemAry(conv.GetIdEA_fromCad(3,Cad::EDGE),world);
    field_value_setter = Fem::Field::CFieldValueSetter(id_field_bc0,world);
    field_value_setter.SetMathExp("0.5*sin(0.05*t)", 0,Fem::Field::VELOCITY, world);
		unsigned int id_field_bc1;
		{
			std::vector<unsigned int> id_ea_bc1;
			id_ea_bc1.push_back(conv.GetIdEA_fromCad(1,Cad::EDGE));
			id_ea_bc1.push_back(conv.GetIdEA_fromCad(2,Cad::EDGE));
			id_ea_bc1.push_back(conv.GetIdEA_fromCad(4,Cad::EDGE));
			id_field_bc1 = fluid.AddFixElemAry(id_ea_bc1,world);
		}
		fluid.SetRho(0.1);
		fluid.SetMyu(0.0002);
		fluid.SetStokes();
		fluid.SetIsStationary(true);
		dt = 0.5;
		fluid.SetTimeIntegrationParameter(dt);

		// registration of visualization objects
		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerVector(id_field_velo,world) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_press,true,world, id_field_press) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerEdge(id_field_velo,true,world) );
//		drawer_ary.PushBack( new Fem::Field::View::CDrawerImageBasedFlowVis(id_field_velo,world,1) );
//		drawer_ary.PushBack( new Fem::Field::View::CDrawerStreamline(id_field_velo,world) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 1 )	// cavity flow (non-stationary storks)
	{
		fluid.SetIsStationary(false);
	}
	else if( iprob == 2 )	// cavity flow (non-stationary storks with larger rho )
	{
		fluid.SetRho(0.5);
	}
	else if( iprob == 3 )	// cavity flowÅCnon-static Naiver-Stokes flow
	{
		fluid.SetRho(0.02);
		fluid.SetMyu(0.00001);
		fluid.SetNavierStokes();
	}
	else if( iprob == 4 )	// cavity flowÅCbubble interpolationÅCconstant storks
	{
		fluid.Clear();
		fluid.SetInterpolationBubble();
		fluid.UpdateDomain_Field(id_base,world);
		fluid.SetStokes();
		fluid.SetIsStationary(true);
		const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);
		unsigned int id_field_press = fluid.GetIdField_Press();
        std::cout << "press : " << id_field_press << std::endl;
		unsigned int id_field_press_bc0 = world.GetPartialField(id_field_press,conv.GetIdEA_fromCad(1,Cad::VERTEX));
		fluid.AddFixField(id_field_press_bc0,world);

		unsigned int id_field_velo  = fluid.GetIdField_Velo();
        std::cout << "velo : " << id_field_velo << std::endl;
		unsigned int id_field_bc0 = fluid.AddFixElemAry(conv.GetIdEA_fromCad(3,Cad::EDGE),world);
    field_value_setter = Fem::Field::CFieldValueSetter(id_field_bc0,world);
    field_value_setter.SetMathExp("0.5*sin(0.05*t)", 0,Fem::Field::VELOCITY, world);
		unsigned int id_field_bc1;
		{
			std::vector<unsigned int> id_ea_bc1;
			id_ea_bc1.push_back(conv.GetIdEA_fromCad(1,Cad::EDGE));
			id_ea_bc1.push_back(conv.GetIdEA_fromCad(2,Cad::EDGE));
			id_ea_bc1.push_back(conv.GetIdEA_fromCad(4,Cad::EDGE));
			id_field_bc1 = fluid.AddFixElemAry(id_ea_bc1,world);
		}
		
		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerVector(id_field_velo,world) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_press,true,world, id_field_press) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 5 ) // l shaped flow
	{
		Cad::CCadObj2D cad_2d;
 		{	// define shape
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0, 0.0) );
			vec_ary.push_back( Com::CVector2D(1.0, 0.0) );
			vec_ary.push_back( Com::CVector2D(1.0, 0.5) );
			vec_ary.push_back( Com::CVector2D(0.5, 0.5) );
			vec_ary.push_back( Com::CVector2D(0.5, 1.0) );
			vec_ary.push_back( Com::CVector2D(0.0, 1.0) );
			cad_2d.AddPolygon( vec_ary );
		}
		cur_time = 0;
		world.Clear();
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.04) );
		const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);
		fluid.UnSetInterpolationBubble();
		fluid.Clear();
		fluid.UpdateDomain_Field(id_base,world);
		fluid.SetTimeIntegrationParameter(dt);

		unsigned int id_field_press = fluid.GetIdField_Press();
//		unsigned int id_field_press_bc0 = world.GetPartialField(id_field_press,10);
//		fluid.AddFixField(id_field_press_bc0,world);

		unsigned int id_field_velo  = fluid.GetIdField_Velo();
		unsigned int id_field_bc0 = fluid.AddFixElemAry(conv.GetIdEA_fromCad(2,Cad::EDGE),world);
    field_value_setter = Fem::Field::CFieldValueSetter(id_field_bc0,world);
    field_value_setter.SetMathExp("0.1*sin(0.1*t)", 0,Fem::Field::VELOCITY, world);    
		unsigned int id_field_bc1;
		{
			std::vector<unsigned int> id_ea_bc1;
			id_ea_bc1.push_back( conv.GetIdEA_fromCad(1,Cad::EDGE) );
			id_ea_bc1.push_back( conv.GetIdEA_fromCad(3,Cad::EDGE) );
			id_ea_bc1.push_back( conv.GetIdEA_fromCad(4,Cad::EDGE) );
			id_ea_bc1.push_back( conv.GetIdEA_fromCad(6,Cad::EDGE) );
			id_field_bc1 = fluid.AddFixElemAry(id_ea_bc1,world);
		}
		fluid.SetRho(0.1);
		fluid.SetMyu(0.0002);
//		fluid.UnSetStationary(world);
		fluid.SetIsStationary(true);
		fluid.SetStokes();
		fluid.SetTimeIntegrationParameter(dt);

		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerVector(id_field_velo,world) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_press,true,world, id_field_press) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerEdge(id_field_velo,true,world) );
//		drawer_ary.PushBack( new Fem::Field::View::CDrawerStreamline(id_field_velo,world) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 6 )
	{
		fluid.SetIsStationary(false);
		fluid.SetNavierStokes();
		fluid.SetRho(1.3);
		fluid.SetMyu(0.0002);
	}
	else if( iprob == 7 )
	{
		Cad::CCadObj2D cad_2d;
		{	// define shape ( square devided vertically in half )
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,1.0) );
			vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			cad_2d.AddPolygon( vec_ary );
			const unsigned int id_v1 = cad_2d.AddVertex(Cad::EDGE,1,Com::CVector2D(0.5,0.0)).id_v_add;
			const unsigned int id_v2 = cad_2d.AddVertex(Cad::EDGE,3,Com::CVector2D(0.5,1.0)).id_v_add;
			cad_2d.ConnectVertex_Line(id_v1,id_v2);
		}

		cur_time = 0;
		world.Clear();
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.03) );
		const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);
		fluid.Clear();
		fluid.UpdateDomain_FieldElemAry(id_base, conv.GetIdEA_fromCad(2,Cad::LOOP) ,world);
		unsigned int id_field_press = fluid.GetIdField_Press();
		unsigned int id_field_velo  = fluid.GetIdField_Velo();

		unsigned int id_field_bc1 = fluid.AddFixElemAry( conv.GetIdEA_fromCad(3,Cad::EDGE) ,world);
    field_value_setter = Fem::Field::CFieldValueSetter(id_field_bc1,world);
    field_value_setter.SetMathExp("0.3*sin(0.5*t)", 1,Fem::Field::VELOCITY, world);        
		unsigned int id_field_bc2 = fluid.AddFixElemAry( conv.GetIdEA_fromCad(5,Cad::EDGE) ,world);

		fluid.SetRho(0.1);
		fluid.SetMyu(0.0002);
		fluid.SetStokes();
		fluid.SetIsStationary(true);
		fluid.SetTimeIntegrationParameter(dt);

		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerVector(id_field_velo,world) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_press,true,world,id_field_press) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerEdge(id_base,true,world) );
//		drawer_ary.PushBack( new Fem::Field::View::CDrawerVector(id_field_velo,world) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 8 ){
		fluid.SetIsStationary(false);
	}
	else if( iprob == 9 ){
		fluid.SetNavierStokes();
	}
	else if( iprob == 10 ){	// Karman vortex sheet problem
		Cad::CCadObj2D cad_2d;
		{	// define shape (hole in rectangle)
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(2.0,0.0) );
			vec_ary.push_back( Com::CVector2D(2.0,0.6) );
			vec_ary.push_back( Com::CVector2D(0.0,0.6) );
			const unsigned int id_l = cad_2d.AddPolygon( vec_ary ).id_l_add;
			const unsigned int id_v1 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.2,0.2)).id_v_add;
			const unsigned int id_v2 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.3,0.2)).id_v_add;
			const unsigned int id_v3 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.3,0.4)).id_v_add;
			const unsigned int id_v4 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.2,0.4)).id_v_add;
			cad_2d.ConnectVertex_Line(id_v1,id_v2);
			cad_2d.ConnectVertex_Line(id_v2,id_v3);
			cad_2d.ConnectVertex_Line(id_v3,id_v4);
			cad_2d.ConnectVertex_Line(id_v4,id_v1);
		}

		cur_time = 0;
		world.Clear();
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.05) );
		const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);

		fluid.Clear();
		fluid.UpdateDomain_FieldElemAry(id_base,conv.GetIdEA_fromCad(1,Cad::LOOP),world);
		unsigned int id_field_press = fluid.GetIdField_Press();
		unsigned int id_field_velo  = fluid.GetIdField_Velo();

		unsigned int id_field_bc1 = fluid.AddFixElemAry( conv.GetIdEA_fromCad(4,Cad::EDGE) ,world);
    field_value_setter = Fem::Field::CFieldValueSetter();
    Fem::Field::SetFieldValue_Constant(id_field_bc1,0,Fem::Field::VELOCITY,world,0.1);
		{
			std::vector<unsigned int> aIdEAFixVelo;
			aIdEAFixVelo.push_back( conv.GetIdEA_fromCad(1,Cad::EDGE) );
			aIdEAFixVelo.push_back( conv.GetIdEA_fromCad(3,Cad::EDGE) );
			unsigned int id_field_bc2 = fluid.AddFixElemAry(aIdEAFixVelo,world);
		}
		{
			std::vector<unsigned int> aIdEAFixVelo;
			aIdEAFixVelo.push_back( conv.GetIdEA_fromCad(5,Cad::EDGE) );
			aIdEAFixVelo.push_back( conv.GetIdEA_fromCad(6,Cad::EDGE) );
			aIdEAFixVelo.push_back( conv.GetIdEA_fromCad(7,Cad::EDGE) );
			aIdEAFixVelo.push_back( conv.GetIdEA_fromCad(8,Cad::EDGE) );
			unsigned int id_field_bc3 = fluid.AddFixElemAry(aIdEAFixVelo,world);
		}
		dt = 0.13;
		fluid.SetRho(200.0);
		fluid.SetMyu(0.0001);
		fluid.SetNavierStokes();
		fluid.SetTimeIntegrationParameter(dt);

		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerVector(id_field_velo,world) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_press,true,world,id_field_press) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerEdge(id_field_velo,true,world) );
//		drawer_ary.PushBack( new Fem::Field::View::CDrawerImageBasedFlowVis(id_field_velo,world,1) );
//		drawer_ary.PushBack( new Fem::Field::View::CDrawerStreamline(id_field_velo,world) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 11 )	// two element array problem
	{
		Cad::CCadObj2D cad_2d;
 		{	// define shape
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0, 0.0) );
			vec_ary.push_back( Com::CVector2D(1.0, 0.0) );
			vec_ary.push_back( Com::CVector2D(1.0, 0.5) );
			vec_ary.push_back( Com::CVector2D(1.0, 1.0) );
			vec_ary.push_back( Com::CVector2D(0.0, 1.0) );
			vec_ary.push_back( Com::CVector2D(0.0, 0.5) );
			cad_2d.AddPolygon( vec_ary );
			cad_2d.ConnectVertex_Line(3,6);
		}
		cur_time = 0;
		world.Clear();
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.04) );
		const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);

		fluid.Clear();
		fluid.UnSetInterpolationBubble();
		fluid.UpdateDomain_Field(id_base,world);
		unsigned int id_field_press = fluid.GetIdField_Press();
		unsigned int id_field_velo  = fluid.GetIdField_Velo();
		unsigned int id_field_bc0 = fluid.AddFixElemAry(conv.GetIdEA_fromCad(2,Cad::EDGE),world);
    field_value_setter = Fem::Field::CFieldValueSetter(id_field_bc0,world);
    field_value_setter.SetMathExp("0.1*sin(0.1*t)", 0,Fem::Field::VELOCITY, world);            
		unsigned int id_field_bc1;
		{
			std::vector<unsigned int> id_ea_bc1;
			id_ea_bc1.push_back( conv.GetIdEA_fromCad(1,Cad::EDGE) );
			id_ea_bc1.push_back( conv.GetIdEA_fromCad(3,Cad::EDGE) );
			id_ea_bc1.push_back( conv.GetIdEA_fromCad(4,Cad::EDGE) );
			id_ea_bc1.push_back( conv.GetIdEA_fromCad(6,Cad::EDGE) );
			id_field_bc1 = fluid.AddFixElemAry(id_ea_bc1,world);
		}
		dt = 0.8;
		fluid.SetRho(1);
		fluid.SetMyu(0.0002);
		fluid.SetIsStationary(false);
//		fluid.SetStationary(world);
		fluid.SetNavierStokes();
		fluid.SetTimeIntegrationParameter(dt);

		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerVector(id_field_velo,world) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_press,true,world,id_field_press) );
//		drawer_ary.PushBack( new Fem::Field::View::CDrawerStreamline(id_field_velo,world) );
//		drawer_ary.PushBack( new Fem::Field::View::CDrawerFaceContour(id_field_press,world) );    
		drawer_ary.InitTrans( camera );
	}
  else if( iprob == 12 ){ // back-step flow
		Cad::CCadObj2D cad_2d;
 		{	// define shape
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D( 0.0, 0.0) );
			vec_ary.push_back( Com::CVector2D( 1.4, 0.0) );
			vec_ary.push_back( Com::CVector2D( 1.5, 0.0) );
			vec_ary.push_back( Com::CVector2D( 1.5, 1.0) );
			vec_ary.push_back( Com::CVector2D( 1.4, 1.0) );
			vec_ary.push_back( Com::CVector2D(-0.5, 1.0) );
			vec_ary.push_back( Com::CVector2D(-0.5, 0.7) );
			vec_ary.push_back( Com::CVector2D( 0.0, 0.7) );
			cad_2d.AddPolygon( vec_ary );
			cad_2d.ConnectVertex_Line(2,5);	
		}
		cur_time = 0;
		world.Clear();
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.04) );
		const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);

		fluid.Clear();
		fluid.UnSetInterpolationBubble();
		fluid.UpdateDomain_Field(id_base,world);

		unsigned int id_field_press = fluid.GetIdField_Press();
		unsigned int id_field_velo  = fluid.GetIdField_Velo();
		unsigned int id_field_bc0 = fluid.AddFixElemAry(conv.GetIdEA_fromCad(6,Cad::EDGE),world);
    field_value_setter = Fem::Field::CFieldValueSetter();
    Fem::Field::SetFieldValue_Constant(id_field_bc0,0,Fem::Field::VELOCITY, world, 0.2);
		unsigned int id_field_bc1;
		{
			std::vector<unsigned int> id_ea_bc1;
			id_ea_bc1.push_back( conv.GetIdEA_fromCad(1,Cad::EDGE) );
			id_ea_bc1.push_back( conv.GetIdEA_fromCad(2,Cad::EDGE) );
			id_ea_bc1.push_back( conv.GetIdEA_fromCad(4,Cad::EDGE) );
			id_ea_bc1.push_back( conv.GetIdEA_fromCad(5,Cad::EDGE) );
			id_ea_bc1.push_back( conv.GetIdEA_fromCad(7,Cad::EDGE) );
			id_ea_bc1.push_back( conv.GetIdEA_fromCad(8,Cad::EDGE) );
			id_field_bc1 = fluid.AddFixElemAry(id_ea_bc1,world);
		}
    field_value_setter = Fem::Field::CFieldValueSetter();    
		fluid.SetRho(5);
		fluid.SetMyu(0.0002);
		fluid.SetIsStationary(false);
//		fluid.SetStationary(world);
		fluid.SetNavierStokes();
		fluid.SetTimeIntegrationParameter(dt);
		{
			Fem::Eqn::CEqn_Fluid2D eqn_fluid = fluid.GetEquation( conv.GetIdEA_fromCad(2,Cad::LOOP) );
			eqn_fluid.SetMyu(0.01);
			fluid.SetEquation(eqn_fluid);
		}

		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerVector(id_field_velo,world) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_press,true,world,id_field_press) );
    //		drawer_ary.PushBack( new Fem::Field::View::CDrawerFaceContour(id_field_press,world) );
    //		drawer_ary.PushBack( new Fem::Field::View::CDrawerImageBasedFlowVis(id_field_velo,world,1) );
    //		drawer_ary.PushBack( new Fem::Field::View::CDrawerStreamline(id_field_velo,world) );        
		drawer_ary.InitTrans( camera );
	}
	if( iprob == 13 ){ // buoyant problem
		Cad::CCadObj2D cad_2d;
 		{	// define shape
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D( 0.0, 0.0) );
			vec_ary.push_back( Com::CVector2D( 3.0, 0.0) );
			vec_ary.push_back( Com::CVector2D( 3.0, 3.0) );
			vec_ary.push_back( Com::CVector2D( 0.0, 3.0) );
			const unsigned int id_l = cad_2d.AddPolygon( vec_ary ).id_l_add;
			vec_ary.clear();
			vec_ary.push_back( Com::CVector2D( 1.0, 1.0) );
			vec_ary.push_back( Com::CVector2D( 1.5, 1.0) );
			vec_ary.push_back( Com::CVector2D( 1.5, 1.5) );
			vec_ary.push_back( Com::CVector2D( 1.0, 1.5) );
			cad_2d.AddPolygon( vec_ary, id_l );
		}
		cur_time = 0;
		world.Clear();
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.1) );
		const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);

		fluid.Clear();
		fluid.UnSetInterpolationBubble();
		fluid.UpdateDomain_Field(id_base,world);

		unsigned int id_field_press = fluid.GetIdField_Press();
		unsigned int id_field_velo  = fluid.GetIdField_Velo();
		unsigned int id_field_bc0;
		{
			std::vector<unsigned int> id_ea_bc;
			id_ea_bc.push_back( conv.GetIdEA_fromCad(1,Cad::EDGE) );
			id_ea_bc.push_back( conv.GetIdEA_fromCad(2,Cad::EDGE) );
			id_ea_bc.push_back( conv.GetIdEA_fromCad(3,Cad::EDGE) );
			id_ea_bc.push_back( conv.GetIdEA_fromCad(4,Cad::EDGE) );
			id_field_bc0 = fluid.AddFixElemAry(id_ea_bc,world);
		}
    field_value_setter = Fem::Field::CFieldValueSetter();    
		fluid.SetRho(5);
		fluid.SetMyu(0.005);
		fluid.SetIsStationary(false);
		fluid.SetNavierStokes();
		fluid.SetTimeIntegrationParameter(dt);
		{
			Fem::Eqn::CEqn_Fluid2D eqn_fluid = fluid.GetEquation( conv.GetIdEA_fromCad(2,Cad::LOOP) );
			eqn_fluid.SetBodyForce(0.0,0.5);
			fluid.SetEquation(eqn_fluid);
		}

		// registration of visualization objects
		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerVector(id_field_velo,world) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_press,true,world,id_field_press) );    
    //		drawer_ary.PushBack( new Fem::Field::View::CDrawerFaceContour(id_field_press,world);
    //		drawer_ary.PushBack( new Fem::Field::View::CDrawerStreamline(id_field_velo,world) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 14 )
	{
		Cad::CCadObj2D cad_2d;
		unsigned int id_l;
		unsigned int id_e1, id_e2, id_e3,id_e4,id_e5,id_e6;
		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D( 0.0,0.0) );
			vec_ary.push_back( Com::CVector2D( 0.5,0.0) );
			vec_ary.push_back( Com::CVector2D( 2.0,0.0) );
			vec_ary.push_back( Com::CVector2D( 2.0,1.0) );
			vec_ary.push_back( Com::CVector2D( 0.0,1.0) );
			id_l = cad_2d.AddPolygon( vec_ary ).id_l_add;
			unsigned int id_v1 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(0.5,0.5) ).id_v_add;
			id_e1 = cad_2d.ConnectVertex_Line(2,id_v1).id_e_add;
			unsigned int id_v2 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(1.0,0.3) ).id_v_add;
			unsigned int id_v3 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(1.0,0.9) ).id_v_add;
			id_e2 = cad_2d.ConnectVertex_Line(id_v2,id_v3).id_e_add;
			unsigned int id_v4 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(1.5,0.4) ).id_v_add;
			unsigned int id_v5 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(1.5,0.1) ).id_v_add;
			unsigned int id_v6 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(1.5,0.7) ).id_v_add;
			unsigned int id_v7 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(1.2,0.4) ).id_v_add;
			unsigned int id_v8 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(1.8,0.4) ).id_v_add;
			id_e3 = cad_2d.ConnectVertex_Line(id_v4,id_v5).id_e_add;
			id_e4 = cad_2d.ConnectVertex_Line(id_v4,id_v6).id_e_add;
			id_e5 = cad_2d.ConnectVertex_Line(id_v4,id_v7).id_e_add;
			id_e6 = cad_2d.ConnectVertex_Line(id_v4,id_v8).id_e_add;
		}
		Msh::CMesher2D mesh_2d(cad_2d,0.05);

		cur_time = 0;
		world.Clear();
		id_base = world.AddMesh(mesh_2d);
		const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);
		unsigned int id_base2 = 0;
		{
			std::vector<unsigned int> mapVal2Co;
			std::vector< std::vector<int> > aLnods;
			{
				std::vector<unsigned int> aIdMsh_Inc;
				aIdMsh_Inc.push_back( mesh_2d.GetElemID_FromCadID(id_l,Cad::LOOP) );
				std::vector<unsigned int> aIdMshBar_Cut;
				aIdMshBar_Cut.push_back( mesh_2d.GetElemID_FromCadID(id_e1,Cad::EDGE) );
				aIdMshBar_Cut.push_back( mesh_2d.GetElemID_FromCadID(id_e2,Cad::EDGE) );
				aIdMshBar_Cut.push_back( mesh_2d.GetElemID_FromCadID(id_e3,Cad::EDGE) );
				aIdMshBar_Cut.push_back( mesh_2d.GetElemID_FromCadID(id_e4,Cad::EDGE) );
				aIdMshBar_Cut.push_back( mesh_2d.GetElemID_FromCadID(id_e5,Cad::EDGE) );
				aIdMshBar_Cut.push_back( mesh_2d.GetElemID_FromCadID(id_e6,Cad::EDGE) );
				mesh_2d.GetClipedMesh(aLnods,mapVal2Co, aIdMsh_Inc,aIdMshBar_Cut);
			}
			std::vector<unsigned int> aIdEA_Inc;
			aIdEA_Inc.push_back( conv.GetIdEA_fromCad(1,Cad::LOOP) );
			id_base2 = world.SetCustomBaseField(id_base,aIdEA_Inc,aLnods,mapVal2Co);
		}

		fluid.Clear();
		fluid.UnSetInterpolationBubble();
		fluid.UpdateDomain_FieldVeloPress(id_base,id_base2,world);
		unsigned int id_field_press = fluid.GetIdField_Press();
		unsigned int id_field_velo  = fluid.GetIdField_Velo();
		unsigned int id_field_bc0;
		{
			std::vector<unsigned int> id_ea_bc;
			id_ea_bc.push_back( conv.GetIdEA_fromCad(1,Cad::EDGE) );
			id_ea_bc.push_back( conv.GetIdEA_fromCad(2,Cad::EDGE) );
			id_ea_bc.push_back( conv.GetIdEA_fromCad(4,Cad::EDGE) );
			id_ea_bc.push_back( conv.GetIdEA_fromCad(id_e1,Cad::EDGE) );
			id_ea_bc.push_back( conv.GetIdEA_fromCad(id_e2,Cad::EDGE) );
			id_ea_bc.push_back( conv.GetIdEA_fromCad(id_e3,Cad::EDGE) );
			id_ea_bc.push_back( conv.GetIdEA_fromCad(id_e4,Cad::EDGE) );
			id_ea_bc.push_back( conv.GetIdEA_fromCad(id_e5,Cad::EDGE) );
			id_ea_bc.push_back( conv.GetIdEA_fromCad(id_e6,Cad::EDGE) );
			id_field_bc0 = fluid.AddFixElemAry(id_ea_bc,world);
		}
		unsigned int id_field_bc1 = fluid.AddFixElemAry(conv.GetIdEA_fromCad(5,Cad::EDGE),world);
    field_value_setter = Fem::Field::CFieldValueSetter(id_field_bc1,world);
    field_value_setter.SetMathExp("0.1*sin(t*PI*0.1+0.01)", 0,Fem::Field::VELOCITY, world);            

		dt = 0.8;
		fluid.SetRho(10);
		fluid.SetMyu(0.005);
		fluid.SetIsStationary(false);
//		fluid.SetStokes();
		fluid.SetNavierStokes();
		fluid.SetTimeIntegrationParameter(dt);

		// registration of visualization objects
		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerVector(     id_field_velo, world) );    
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_press,true,world,id_field_press) );    
//		drawer_ary.PushBack( new Fem::Field::View::CDrawerImageBasedFlowVis(id_field_velo,world,1) );
//		drawer_ary.PushBack( new Fem::Field::View::CDrawerStreamline(id_field_velo,world) );
//		drawer_ary.PushBack( new Fem::Field::View::CDrawerFaceContour(id_field_press,world) );
		drawer_ary.InitTrans( camera );
	}
	
	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
	Com::View::SetProjectionTransform(camera);

	iprob++;
	if( iprob == nprob ){ iprob = 0; }
}

int main(int argc,char* argv[])
{
	// initialize GLUT
	glutInitWindowPosition(200,200);
	glutInitWindowSize(400, 300);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
	glutCreateWindow("DelFEM Demo");

	// set call-back functions
	glutIdleFunc(myGlutIdle);
	glutKeyboardFunc(myGlutKeyboard);
	glutDisplayFunc(myGlutDisplay);
	glutReshapeFunc(myGlutResize);
	glutSpecialFunc(myGlutSpecial);;
	glutMotionFunc(myGlutMotion);
	glutMouseFunc(myGlutMouse);
	
	SetNewProblem();
	glutMainLoop();
	return 0;
}
