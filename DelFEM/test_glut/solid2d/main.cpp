////////////////////////////////////////////////////////////////
//                                                            //
//		DelFEM demo : Solid2D                                 //
//                                                            //
//          Copy Rights (c) Nobuyuki Umetani 2008             //
//          e-mail : numetani@gmail.com                       //
////////////////////////////////////////////////////////////////

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif
#define for if(0); else for

#include <cassert>
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <math.h>
#include <fstream>
#include <time.h>
#include <cstdlib> //(exit)

#if defined(__APPLE__) && defined(__MACH__)
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif

#include "delfem/camera.h"
#include "delfem/glut_utility.h"

#include "delfem/cad_obj2d.h"
#include "delfem/mesh3d.h"

#include "delfem/field.h"
#include "delfem/field_world.h"
#include "delfem/field_value_setter.h"

#include "delfem/drawer_field.h"
#include "delfem/drawer_field_face.h"
#include "delfem/drawer_field_edge.h"
#include "delfem/drawer_field_vector.h"

#include "delfem/eqnsys_solid.h"

using namespace Fem::Ls;
using namespace Fem::Field;

Com::View::CCamera camera;
double mov_begin_x, mov_begin_y;
bool is_animation = true;


void myGlutResize(int w, int h)
{
	camera.SetWindowAspect((double)w/h);
	::glViewport(0, 0, w, h);
	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
	Com::View::SetProjectionTransform(camera);
	::glutPostRedisplay();
}

void myGlutMotion( int x, int y ){
	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
	const int win_w = viewport[2];
	const int win_h = viewport[3];
	const double mov_end_x = (2.0*x-win_w)/win_w;
	const double mov_end_y = (win_h-2.0*y)/win_h;
//	camera.MouseRotation(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 
	camera.MousePan(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 
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
}

void SetNewProblem();
void myGlutKeyboard(unsigned char Key, int x, int y)
{
	switch(Key)
	{
	case 'q':
	case 'Q':
	case '\033':
		exit(0);  // '\033' is ESC key in ASCII
	case 'a':
		if( is_animation ){ is_animation = false; }
		else{ is_animation = true; }
		break;
	case ' ':
		SetNewProblem();
		break;
	}
	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
	Com::View::SetProjectionTransform(camera);
	::glutPostRedisplay();
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

void myGlutIdle(){
	::glutPostRedisplay();
}

////////////////////////////////

Fem::Field::CFieldWorld world;
Fem::Field::View::CDrawerArrayField drawer_ary;
Fem::Field::CFieldValueSetter field_value_setter;
double cur_time = 0.0;
double dt = 0.05;
Fem::Eqn::CEqnSystem_Solid2D solid;
unsigned int id_field_disp;
unsigned int id_field_equiv_stress;	// handle of equiv stress (scalar value)
unsigned int id_field_stress;	// handle of stress field (sym-tensor field)

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

	if( is_animation ){
		cur_time += dt;
//		world.FieldValueExec(cur_time);
    field_value_setter.ExecuteValue(cur_time,world);
		solid.Solve(world);    
		if( id_field_equiv_stress != 0 ){ solid.SetEquivStressValue(id_field_equiv_stress,world); }
		if( id_field_stress       != 0 ){ solid.SetStressValue(     id_field_stress,      world); }
		if( solid.GetAry_ItrNormRes().size() > 0 ){
			std::cout << "Iter : " << solid.GetAry_ItrNormRes()[0].first << " ";
			std::cout << "Res : " << solid.GetAry_ItrNormRes()[0].second << std::endl;
		}
//		world.FieldValueDependExec();
		drawer_ary.Update(world);
	}

	drawer_ary.Draw();
	ShowFPS();
	::glutSwapBuffers();
}

void SetNewProblem()
{
	const unsigned int nprob = 30;	// number of problem settings
	static unsigned int iprob = 0;

	static unsigned int id_field_disp_fix0 = 0;
	static unsigned int id_field_temp = 0;

	if( iprob == 0 )	// linear solid stationary analysis
	{
		id_field_disp_fix0 = 0;
		id_field_temp = 0;
		id_field_stress = 0;
		id_field_equiv_stress = 0;
		////////////////
		Cad::CCadObj2D cad_2d;
		{	// define shape
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(5.0,0.0) );
			vec_ary.push_back( Com::CVector2D(5.0,1.0) );
			vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			cad_2d.AddPolygon( vec_ary );
		}
		world.Clear();
		const unsigned int id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.1) );
		const CIDConvEAMshCad conv = world.GetIDConverter(id_base);

		solid.Clear();
    solid.UpdateDomain_Field(id_base, world);
		solid.SetSaveStiffMat(false);	
		solid.SetStationary(true);
		// Setting Material Parameter
		solid.SetYoungPoisson(10.0,0.3,true);	// planter stress
		solid.SetGeometricalNonlinear(false);
		solid.SetGravitation(0.0,0.0);
		solid.SetTimeIntegrationParameter(dt,0.7);

		unsigned int id_field_bc0 = solid.AddFixElemAry(conv.GetIdEA_fromCad(2,Cad::EDGE),world);
		unsigned int id_field_bc1 = solid.AddFixElemAry(conv.GetIdEA_fromCad(4,Cad::EDGE),world);
    
    field_value_setter = CFieldValueSetter(id_field_bc0,world);
    field_value_setter.SetMathExp("sin(t*PI*2*0.1)", 1,Fem::Field::VALUE, world);	// oscilate bc1_field y axis

		// Setting Visualiziation
		drawer_ary.Clear();
		id_field_disp = solid.GetIdField_Disp();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_disp,false,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,false,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,true ,world) );
		drawer_ary.InitTrans(camera);	// set view transformation
	}
	else if( iprob == 1 )	// save stiffness matrix for the efficency of computation
	{
		solid.SetSaveStiffMat(true);
	}
	else if( iprob == 2 )	// non-stationary analysis
	{
		solid.SetSaveStiffMat(true);
		solid.SetStationary(false);
	}
	else if( iprob == 3 )	// set stiffer material
	{
		solid.SetYoungPoisson(50,0.3,true);
	}
	else if( iprob == 4 )	// set more stiffer material
	{
		solid.SetYoungPoisson(100,0.3,true);
	}
	else if( iprob == 5 )	// geometrical non-linear stationaly
	{
		solid.SetStationary(true);
		solid.SetGeometricalNonlinear(true);
	}
	else if( iprob == 6 )	// geometrical non-linear non-stationary
	{
		solid.SetYoungPoisson(10,0.0,true);
		solid.SetStationary(false);
		solid.SetGeometricalNonlinear(true);
	}
	else if( iprob == 7 )	// display equivalent stress field in deformedconfigulation
	{
		id_field_equiv_stress = world.MakeField_FieldElemDim(id_field_disp,2,SCALAR,VALUE,BUBBLE);
		solid.SetGeometricalNonlinear(false);
		solid.SetStationary(true);
		// set up visualization
		drawer_ary.Clear();
		id_field_disp = solid.GetIdField_Disp();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_equiv_stress,false,world,id_field_equiv_stress, 0,0.5) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,false,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,true ,world) );
		drawer_ary.InitTrans(camera);	// init view transformation
	}
	else if( iprob == 8 )	// display equivalent stress field in initial configulation
	{
		// set up visualization
		drawer_ary.Clear();
		id_field_disp = solid.GetIdField_Disp();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_disp,false,world,id_field_equiv_stress, 0,0.5) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,false,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,true ,world) );
		drawer_ary.InitTrans(camera);	// init view transformation
	}
	else if( iprob == 9 )	// thermal-solid analysis
	{
		id_field_equiv_stress = 0;
		id_field_stress = 0;
		////////////////
		Cad::CCadObj2D cad_2d;
 		{	// define shape
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(3.0,0.0) );
			vec_ary.push_back( Com::CVector2D(3.0,1.0) );
			vec_ary.push_back( Com::CVector2D(2.0,1.0) );
			vec_ary.push_back( Com::CVector2D(1.0,1.0) );
			vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			cad_2d.AddPolygon( vec_ary );
		}
		world.Clear();
		const unsigned int id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.1) );
		CIDConvEAMshCad conv = world.GetIDConverter(id_base);

		solid.UpdateDomain_Field(id_base,world);
		solid.SetSaveStiffMat(false);
		solid.SetStationary(true);
		// set material property
		solid.SetYoungPoisson(10.0,0.3,true);	// set planer stress
		solid.SetGeometricalNonlinear(false);	// geometricaly linear model
		solid.SetGravitation(0.0,-0.1);
		solid.SetTimeIntegrationParameter(dt);	// set time setp
		
		unsigned int id_field_bc0 = solid.AddFixElemAry(conv.GetIdEA_fromCad(2,Cad::EDGE),world);
		unsigned int id_field_bc1 = solid.AddFixElemAry(conv.GetIdEA_fromCad(6,Cad::EDGE),world);

		// set temparature field
		id_field_temp = world.MakeField_FieldElemDim(id_field_disp,2,SCALAR,VALUE,CORNER);
    field_value_setter = CFieldValueSetter(id_field_temp,world);
    field_value_setter.SetMathExp("sin(6.28*y)*sin(x)*sin(t)", 0,Fem::Field::VALUE, world);	// oscilate bc1_field y axis
    
		solid.SetThermalStress(id_field_temp);
		solid.ClearFixElemAry(3,world);

		drawer_ary.Clear();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_disp,false, world, id_field_temp, -1,1) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,false,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,true ,world) );
		drawer_ary.InitTrans(camera);	// set view transformation
	}
	else if( iprob == 10 )	// show contour in undeformed configuration
	{
		drawer_ary.Clear();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_temp,true, world, id_field_temp, -1,1) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,false,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,true ,world) );
		drawer_ary.InitTrans(camera);	// set view transformation
	}
	else if( iprob == 11 )	// stop considering thermal-effect
	{
		solid.SetThermalStress(0);
	}
	else if( iprob == 12 )
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
		world.Clear();
		const unsigned int id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.05) );
		CIDConvEAMshCad conv = world.GetIDConverter(id_base);		// get ID converter

		solid.SetDomain_FieldEA(id_base,conv.GetIdEA_fromCad(1,Cad::LOOP),world);
		solid.SetSaveStiffMat(true);
		solid.SetStationary(true);
		solid.SetTimeIntegrationParameter(dt);	// set time step
		solid.SetYoungPoisson(2.5,0.3,true);	// planer stress
		solid.SetGeometricalNonlinear(false);	// set geometrical liner
		solid.SetGravitation(0.0,0.0);	// set gravitation

		unsigned int id_field_bc1 = solid.AddFixElemAry(conv.GetIdEA_fromCad(3,Cad::EDGE),world);
		unsigned int id_field_bc2 = solid.AddFixElemAry(conv.GetIdEA_fromCad(1,Cad::EDGE),world);    
    field_value_setter = CFieldValueSetter(id_field_bc1,world);
    field_value_setter.SetMathExp("0.3*sin(1.5*t)", 0,Fem::Field::VALUE, world);	// oscilate bc1_field x axis
    field_value_setter.SetMathExp("0.1*(cos(t)+1)", 1,Fem::Field::VALUE, world);	// oscilate bc1_field y axis

		// set visualization
		drawer_ary.Clear();
		id_field_disp = solid.GetIdField_Disp();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_disp,false,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,false,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,true ,world) );
		drawer_ary.InitTrans(camera);	// set view transformation
	}
	else if( iprob == 13 )
	{
		solid.SetSaveStiffMat(true);
	}
	else if( iprob == 14 )
	{
		solid.SetSaveStiffMat(false);
		solid.SetStationary(false);
	}
	else if( iprob == 15 ){
		solid.SetStationary(true);
		solid.SetGeometricalNonlinear(true);
	}
	else if( iprob == 16 ){
		solid.SetStationary(false);
		solid.SetGeometricalNonlinear(true);
	}
	else if( iprob == 17 )	// hard and soft solid are connected
	{
		Cad::CCadObj2D cad_2d;
		{	// define shape
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

		world.Clear();
		const unsigned int id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.05) );
		CIDConvEAMshCad conv = world.GetIDConverter(id_base);  // get ID converter

		solid.SetDomain_FieldEA(id_base,conv.GetIdEA_fromCad(2,Cad::LOOP),world);
		solid.SetTimeIntegrationParameter(dt);
		solid.SetSaveStiffMat(false);
		solid.SetStationary(true);

		solid.SetYoungPoisson(3.0,0.3,true);
		solid.SetGeometricalNonlinear(false);	// set geometrically linear model
		solid.SetGravitation(0.0,0.0);

		unsigned int id_field_bc1 = solid.AddFixElemAry(conv.GetIdEA_fromCad(3,Cad::EDGE),world);
		unsigned int id_field_bc2 = solid.AddFixElemAry(conv.GetIdEA_fromCad(5,Cad::EDGE),world);    
    field_value_setter = CFieldValueSetter(id_field_bc1,world);
    field_value_setter.SetMathExp("0.3*sin(1.5*t)",     0,Fem::Field::VALUE, world);	// oscilate bc1_field x axis
    field_value_setter.SetMathExp("0.1*(cos(t)+1)+0.1", 1,Fem::Field::VALUE, world);	// oscilate bc1_field y axis
    
		// set up visualization
		drawer_ary.Clear();
		id_field_disp = solid.GetIdField_Disp();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_disp,false,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,false,world) );
//		drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,true ,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_base,true,world) );
		drawer_ary.InitTrans(camera);	// set view trnsformation
	}
	else if( iprob == 18 )
	{
		solid.SetSaveStiffMat(true);
	}
	else if( iprob == 19 )
	{
		solid.SetSaveStiffMat(false);
		solid.SetStationary(false);
	}
	else if( iprob == 20 )
	{
		solid.SetStationary(true);
		solid.SetGeometricalNonlinear(true);
	}
	else if( iprob == 21 )
	{
		solid.SetStationary(false);
		solid.SetGeometricalNonlinear(true);
	}
	else if( iprob == 22 )	// 4 different type of solid combined
	{
		Cad::CCadObj2D cad_2d;
		{	// define shape
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(2.0,0.0) );
			vec_ary.push_back( Com::CVector2D(2.0,0.5) );
			vec_ary.push_back( Com::CVector2D(0.0,0.5) );
			cad_2d.AddPolygon( vec_ary );
			const unsigned int id_v5 = cad_2d.AddVertex(Cad::EDGE,1,Com::CVector2D(1.5,0.0)).id_v_add;
			const unsigned int id_v3 = cad_2d.AddVertex(Cad::EDGE,1,Com::CVector2D(1.0,0.0)).id_v_add;
			const unsigned int id_v1 = cad_2d.AddVertex(Cad::EDGE,1,Com::CVector2D(0.5,0.0)).id_v_add;
			const unsigned int id_v2 = cad_2d.AddVertex(Cad::EDGE,3,Com::CVector2D(0.5,0.5)).id_v_add;
			const unsigned int id_v4 = cad_2d.AddVertex(Cad::EDGE,3,Com::CVector2D(1.0,0.5)).id_v_add;
			const unsigned int id_v6 = cad_2d.AddVertex(Cad::EDGE,3,Com::CVector2D(1.5,0.5)).id_v_add;
			cad_2d.ConnectVertex_Line(id_v1,id_v2);
			cad_2d.ConnectVertex_Line(id_v3,id_v4);
			cad_2d.ConnectVertex_Line(id_v5,id_v6);
		}

		world.Clear();
		const unsigned int id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.05) );	
		const CIDConvEAMshCad& conv = world.GetIDConverter(id_base);  // get ID converter

		solid.UpdateDomain_Field(id_base,world);	// set domain of solid analysis
		solid.SetTimeIntegrationParameter(dt);	// set time step
		solid.SetSaveStiffMat(false);
		solid.SetStationary(true);
		// set material property
		solid.SetYoungPoisson(1.0,0.3,true);
		solid.SetGeometricalNonlinear(false);
		solid.SetGravitation(0.0,-0.0);

		{	// St.Venant-Kirchhoff material
			Fem::Eqn::CEqn_Solid2D eqn = solid.GetEquation(conv.GetIdEA_fromCad(1,Cad::LOOP));
			eqn.SetGeometricalNonlinear(true);
			solid.SetEquation(eqn);
		}
		{	// soft elastic material
			Fem::Eqn::CEqn_Solid2D eqn = solid.GetEquation(conv.GetIdEA_fromCad(2,Cad::LOOP));
			eqn.SetYoungPoisson(0.1,0.3,true);
			solid.SetEquation(eqn);
		}
		unsigned int id_field_temp = world.MakeField_FieldElemAry(id_base, conv.GetIdEA_fromCad(3,Cad::LOOP), SCALAR,VALUE,CORNER);
    field_value_setter = CFieldValueSetter(id_field_temp,world);
    field_value_setter.SetMathExp("0.1*sin(3.14*4*y)*sin(2*t)", 0,Fem::Field::VALUE, world);
		{	// linear elastic material concidering thrmal-solid
			Fem::Eqn::CEqn_Solid2D eqn = solid.GetEquation(conv.GetIdEA_fromCad(3,Cad::LOOP));
			eqn.SetThermalStress(id_field_temp);
			solid.SetEquation(eqn);
		}
		{	// hard elastic material
			Fem::Eqn::CEqn_Solid2D eqn = solid.GetEquation(conv.GetIdEA_fromCad(4,Cad::LOOP));
			eqn.SetYoungPoisson(10,0.3,true);
			solid.SetEquation(eqn);
		}

		id_field_disp_fix0 = solid.AddFixElemAry(conv.GetIdEA_fromCad(2,Cad::EDGE),world);

		// set up visualization
		drawer_ary.Clear();
		id_field_disp = solid.GetIdField_Disp();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_disp,false,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,false,world) );
//		drawer_ary.PushBack( new View::CDrawerEdge(id_base,false,world) );
		drawer_ary.InitTrans(camera);	// set view transformation
	}
	else if( iprob == 23 )
	{
		solid.SetRho(0.0001);
		solid.SetStationary(false);
    field_value_setter = CFieldValueSetter(id_field_disp_fix0,world);
    field_value_setter.SetMathExp("0.5*cos(2*t)", 1,Fem::Field::VALUE, world);    
	}
	else if( iprob == 24 ){
		Cad::CCadObj2D cad_2d;
		{	// define shape
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(2.0,0.0) );
			vec_ary.push_back( Com::CVector2D(2.0,1.0) );
			vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			cad_2d.AddPolygon( vec_ary );
			const unsigned int id_v1 = cad_2d.AddVertex(Cad::EDGE,1,Com::CVector2D(1.0,0.0)).id_v_add;
			const unsigned int id_v2 = cad_2d.AddVertex(Cad::EDGE,3,Com::CVector2D(1.0,1.0)).id_v_add;
			cad_2d.ConnectVertex_Line(id_v1,id_v2);	
		}

		world.Clear();
		const unsigned int id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.05) );
		const CIDConvEAMshCad& conv = world.GetIDConverter(id_base);  // ID converter

		solid.UpdateDomain_Field(id_base,world);	// set displacement field
		solid.SetTimeIntegrationParameter(dt);	// set time step
		solid.SetSaveStiffMat(false);
		solid.SetStationary(false);
		// set material property
		solid.SetYoungPoisson(1.0,0.3,true);	// planter stress
		solid.SetGeometricalNonlinear(false);
		solid.SetGravitation(0.0,-0.0);
        solid.SetRho(0.001);

		{	// soft elastic material
			Fem::Eqn::CEqn_Solid2D eqn = solid.GetEquation(conv.GetIdEA_fromCad(1,Cad::LOOP));
			eqn.SetYoungPoisson(0.1,0.3,true);
			solid.SetEquation(eqn);
		}
		{	// hard elastic material
			Fem::Eqn::CEqn_Solid2D eqn = solid.GetEquation(conv.GetIdEA_fromCad(2,Cad::LOOP));
			eqn.SetYoungPoisson(100000000,0.3,true);
			solid.SetEquation(eqn);
		}

//		id_field_disp_fix0 = solid.AddFixElemAry(conv.GetIdEA_fromCad(2,1),world);
		const unsigned int id_field_bc1 = solid.AddFixElemAry(conv.GetIdEA_fromCad(4,Cad::EDGE),world);
    field_value_setter = CFieldValueSetter(id_field_bc1,world);
    field_value_setter.SetMathExp("0.3*sin(1.5*t)",     0,Fem::Field::VALUE, world);	// oscilate bc1_field x axis
    field_value_setter.SetMathExp("0.1*(cos(t)+1)+0.1", 1,Fem::Field::VALUE, world);	// oscilate bc1_field y axis    

		// set up visualization
		drawer_ary.Clear();
		id_field_disp = solid.GetIdField_Disp();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_disp,false,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,false,world) );
		drawer_ary.InitTrans(camera);	// initialize view transmation
	}
	else if( iprob == 25 )
	{
		Cad::CCadObj2D cad_2d;
		unsigned int id_l;
		unsigned int id_e1,id_e2,id_e3,id_e4,id_e5;
		{	// define shape
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(0.2,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,1.0) );
			vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			id_l = cad_2d.AddPolygon( vec_ary ).id_l_add;
			unsigned int id_v1 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(0.2,0.5) ).id_v_add;
			id_e1 = cad_2d.ConnectVertex_Line(2,id_v1).id_e_add;
			unsigned int id_v2 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(0.5,0.2) ).id_v_add;
			unsigned int id_v3 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(0.5,0.5) ).id_v_add;
			unsigned int id_v4 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(0.5,0.8) ).id_v_add;
			unsigned int id_v5 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(0.8,0.5) ).id_v_add;
			unsigned int id_v6 = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(0.3,0.5) ).id_v_add;
			id_e2 = cad_2d.ConnectVertex_Line(id_v2,id_v3).id_e_add;
			id_e3 = cad_2d.ConnectVertex_Line(id_v3,id_v4).id_e_add;
			id_e4 = cad_2d.ConnectVertex_Line(id_v3,id_v5).id_e_add;
			id_e5 = cad_2d.ConnectVertex_Line(id_v3,id_v6).id_e_add;
		}
		Msh::CMesher2D mesh_2d(cad_2d,0.1);

		world.Clear();
		const unsigned int id_base = world.AddMesh(mesh_2d);
		const CIDConvEAMshCad& conv = world.GetIDConverter(id_base);  // get ID converter
		unsigned int id_base2 = 0;
		{	// cut mesh
			std::vector<unsigned int> mapVal2Co;
			std::vector< std::vector<int> > aLnods;
			{
				std::vector<unsigned int> aIdMsh_Inc;
				aIdMsh_Inc.push_back( mesh_2d.GetElemID_FromCadID(id_l,Cad::LOOP) );
				std::vector<unsigned int> aIdMshCut;
				aIdMshCut.push_back( mesh_2d.GetElemID_FromCadID(id_e1,Cad::EDGE) );
				aIdMshCut.push_back( mesh_2d.GetElemID_FromCadID(id_e2,Cad::EDGE) );
				aIdMshCut.push_back( mesh_2d.GetElemID_FromCadID(id_e3,Cad::EDGE) );
				aIdMshCut.push_back( mesh_2d.GetElemID_FromCadID(id_e4,Cad::EDGE) );
				aIdMshCut.push_back( mesh_2d.GetElemID_FromCadID(id_e5,Cad::EDGE) );
				mesh_2d.GetClipedMesh(aLnods,mapVal2Co, aIdMsh_Inc,aIdMshCut);
			}
			std::vector<unsigned int> aIdEA_Inc;
			aIdEA_Inc.push_back( conv.GetIdEA_fromCad(1,Cad::LOOP) );
			id_base2 = world.SetCustomBaseField(id_base,aIdEA_Inc, aLnods,mapVal2Co);
		}
   
	    solid.UpdateDomain_Field(id_base2, world);
		solid.SetSaveStiffMat(false);	
		solid.SetStationary(true);
		// set material parameter
		solid.SetYoungPoisson(10.0,0.3,true);
		solid.SetGeometricalNonlinear(false);
		solid.SetGravitation(0.0,0.0);
		solid.SetTimeIntegrationParameter(dt,0.7);

		unsigned int id_field_bc0 = solid.AddFixElemAry(conv.GetIdEA_fromCad(3,Cad::EDGE),world);
		unsigned int id_field_bc1 = solid.AddFixElemAry(conv.GetIdEA_fromCad(5,Cad::EDGE),world);
    field_value_setter = CFieldValueSetter(id_field_bc0,world);
    field_value_setter.SetMathExp("0.1*sin(t*PI*2*0.1)",    0,Fem::Field::VALUE, world);	// oscilate bc1_field x axis
    field_value_setter.SetMathExp("0.1*(1-cos(t*PI*2*0.1))",1,Fem::Field::VALUE, world);	// oscilate bc1_field y axis
    
		// set visualization
		drawer_ary.Clear();
		id_field_disp = solid.GetIdField_Disp();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_disp,false,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,false,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,true ,world) );
		drawer_ary.InitTrans(camera);	// set view transformation
	}
	else if( iprob == 26 )
	{
		Cad::CCadObj2D cad_2d;
		unsigned int id_e;
		unsigned int id_l;
		{	// define shape
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(5.0,0.0) );
			vec_ary.push_back( Com::CVector2D(5.0,2.0) );
			vec_ary.push_back( Com::CVector2D(0.0,2.0) );
			id_l = cad_2d.AddPolygon( vec_ary ).id_l_add;
			unsigned int id_v1 = cad_2d.AddVertex(Cad::EDGE,3,   Com::CVector2D(2.5,2.0)).id_v_add;
			unsigned int id_v2 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(2.5,1.0)).id_v_add;
			id_e = cad_2d.ConnectVertex_Line(id_v1,id_v2).id_e_add;
		}
		Msh::CMesher2D mesh_2d(cad_2d,0.2);
		world.Clear();
		cur_time = 0;
		const unsigned int id_base = world.AddMesh(mesh_2d);
		const CIDConvEAMshCad conv = world.GetIDConverter(id_base);
		unsigned int id_base2 = 0;
		{	// cut mesh
			std::vector<unsigned int> mapVal2Co;
			std::vector< std::vector<int> > aLnods;
			{
				std::vector<unsigned int> aIdMsh_Inc;
				aIdMsh_Inc.push_back( mesh_2d.GetElemID_FromCadID(id_l,Cad::LOOP) );
				std::vector<unsigned int> aIdMshCut;
				aIdMshCut.push_back( mesh_2d.GetElemID_FromCadID(id_e,Cad::EDGE) );
				mesh_2d.GetClipedMesh(aLnods,mapVal2Co, aIdMsh_Inc,aIdMshCut);
			}
			std::vector<unsigned int> aIdEA_Inc;
			aIdEA_Inc.push_back( conv.GetIdEA_fromCad(id_l,Cad::LOOP) );
			id_base2 = world.SetCustomBaseField(id_base,aIdEA_Inc, aLnods,mapVal2Co);
		}

	    solid.UpdateDomain_Field(id_base2, world);
		solid.SetSaveStiffMat(false);	
		solid.SetStationary(true);
		// set material parameter
		solid.SetYoungPoisson(10.0,0.3,true);	// set planer stress
		solid.SetGeometricalNonlinear(false);
		solid.SetGravitation(0.0,0.0);
		solid.SetTimeIntegrationParameter(dt,0.7);

		unsigned int id_field_bc0 = solid.AddFixElemAry(conv.GetIdEA_fromCad(2,Cad::EDGE),world);
		unsigned int id_field_bc1 = solid.AddFixElemAry(conv.GetIdEA_fromCad(4,Cad::EDGE),world);
    field_value_setter = CFieldValueSetter(id_field_bc0,world);
    field_value_setter.SetMathExp("0.5*(1-cos(t*PI*2*0.1))", 0,Fem::Field::VALUE, world);	// oscilate bc1_field x axis
    field_value_setter.SetMathExp("0.2*sin(t*PI*2*0.1)",     1,Fem::Field::VALUE, world);	// oscilate bc1_field y axis
		
		id_field_disp = solid.GetIdField_Disp();
		id_field_equiv_stress = world.MakeField_FieldElemDim(id_field_disp,2,SCALAR,VALUE,BUBBLE);

		// set up visualization
		drawer_ary.Clear();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_disp,false,world,id_field_equiv_stress) );
//		drawer_ary.PushBack( new View::CDrawerFace(id_field_disp,false,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,false,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,true ,world) );
		drawer_ary.InitTrans(camera);	// set transformation
	}	
	else if( iprob == 27 )
	{	
		Cad::CCadObj2D cad_2d;
		{	// define shape
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(3.0,0.0) );
			vec_ary.push_back( Com::CVector2D(3.0,1.0) );
			vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			cad_2d.AddPolygon( vec_ary );
		}
		world.Clear();
		const unsigned int id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.3) );
		const CIDConvEAMshCad conv = world.GetIDConverter(id_base);	// Get ID converter
		
		solid.Clear();
	    solid.UpdateDomain_Field(id_base, world);
		solid.SetSaveStiffMat(false);	
		solid.SetStationary(true);
		// Setting Material Parameter
		solid.SetYoungPoisson(2,0.3,true);	// planter stress
		solid.SetGeometricalNonlinear(false);
		solid.SetGravitation(0.0,0.0);
		solid.SetTimeIntegrationParameter(dt,0.7);
		
		unsigned int id_field_bc0 = solid.AddFixElemAry(conv.GetIdEA_fromCad(2,Cad::EDGE),world);
		unsigned int id_field_bc1 = solid.AddFixElemAry(conv.GetIdEA_fromCad(4,Cad::EDGE),world);
    field_value_setter = CFieldValueSetter(id_field_bc0,world);
    field_value_setter.SetMathExp("0.5*sin(t*PI*2*0.1)",    0,Fem::Field::VALUE, world);	// oscilate bc1_field x axis
    field_value_setter.SetMathExp("0.3*(1-cos(t*PI*2*0.1))",1,Fem::Field::VALUE, world);	// oscilate bc1_field y axis
    
		id_field_disp = solid.GetIdField_Disp();
		id_field_stress = world.MakeField_FieldElemDim(id_field_disp,2,STSR2,VALUE,BUBBLE);
//		id_field_equiv_stress = world.MakeField_FieldElemDim(id_field_disp,2,SCALAR,VALUE,BUBBLE);
		
		// set up visualization
		drawer_ary.Clear();
//		drawer_ary.PushBack( new View::CDrawerFace(id_field_equiv_stress,false,world,id_field_equiv_stress, 0,0.5) );
		drawer_ary.PushBack( new View::CDrawerFace(id_field_disp,false,world) );
		drawer_ary.PushBack( new View::CDrawerVector(id_field_stress,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,false,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,true ,world) );
		drawer_ary.InitTrans(camera);	// init view transformation
	}
	else if( iprob == 28 )
	{	
		Cad::CCadObj2D cad_2d;
		{	// define shape
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(0.2,0.0) );
			vec_ary.push_back( Com::CVector2D(0.2,0.5) );
			vec_ary.push_back( Com::CVector2D(0.8,0.5) );
			vec_ary.push_back( Com::CVector2D(0.8,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.7) );
			vec_ary.push_back( Com::CVector2D(0.6,0.7) );
			vec_ary.push_back( Com::CVector2D(0.4,0.7) );
			vec_ary.push_back( Com::CVector2D(0.0,0.7) );
			cad_2d.AddPolygon( vec_ary );
		}
		world.Clear();
		const unsigned int id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.1) );
		const CIDConvEAMshCad conv = world.GetIDConverter(id_base);	// Get ID converter
		
		solid.Clear();
	    solid.UpdateDomain_Field(id_base, world);
		solid.SetSaveStiffMat(false);	
		solid.SetStationary(true);
		// Setting Material Parameter
		solid.SetYoungPoisson(2,0.1,true);	// planter stress
		solid.SetGeometricalNonlinear(false);
		solid.SetGravitation(0.0,0.0);
		solid.SetTimeIntegrationParameter(dt,0.7);
		
		unsigned int id_field_bc0;
		{
			std::vector<unsigned int> aIdEA;
			aIdEA.push_back(conv.GetIdEA_fromCad(1,Cad::EDGE));
			aIdEA.push_back(conv.GetIdEA_fromCad(5,Cad::EDGE));
			id_field_bc0 = solid.AddFixElemAry(aIdEA,world);
		}
		unsigned int id_field_bc1 = solid.AddFixElemAry(conv.GetIdEA_fromCad(8,Cad::EDGE),world);
    field_value_setter = CFieldValueSetter(id_field_bc1,world);
    field_value_setter.SetMathExp("-0.03*(1-cos(t*PI*2*0.1))", 1,Fem::Field::VALUE, world);	// oscilate bc1_field x axis    
		
		id_field_disp = solid.GetIdField_Disp();
		id_field_stress = world.MakeField_FieldElemDim(id_field_disp,2,STSR2,VALUE,BUBBLE);
//		id_field_equiv_stress = world.MakeField_FieldElemDim(id_field_disp,2,SCALAR,VALUE,BUBBLE);

		// set up visualization
		drawer_ary.Clear();		
//		drawer_ary.PushBack( new View::CDrawerFace(id_field_equiv_stress,false,world,id_field_equiv_stress, 0,0.5) );
		drawer_ary.PushBack( new View::CDrawerFace(id_field_disp,false,world) );
		drawer_ary.PushBack( new View::CDrawerVector(id_field_stress,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,false,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,true ,world) );
		drawer_ary.InitTrans(camera);	// init view transformation
	}
	else if( iprob == 29 )
	{	
		Cad::CCadObj2D cad_2d;
		{	// define shape
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(0.2,0.0) );
			vec_ary.push_back( Com::CVector2D(0.2,0.5) );
			vec_ary.push_back( Com::CVector2D(0.8,0.5) );
			vec_ary.push_back( Com::CVector2D(0.8,0.3) );
			vec_ary.push_back( Com::CVector2D(1.0,0.3) );
			vec_ary.push_back( Com::CVector2D(1.0,0.7) );
			vec_ary.push_back( Com::CVector2D(0.6,0.7) );
			vec_ary.push_back( Com::CVector2D(0.4,0.7) );
			vec_ary.push_back( Com::CVector2D(0.0,0.7) );
			cad_2d.AddPolygon( vec_ary );
		}
		world.Clear();
		const unsigned int id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.1) );
		const CIDConvEAMshCad conv = world.GetIDConverter(id_base);	// Get ID converter
		
		solid.Clear();
	    solid.UpdateDomain_Field(id_base, world);
		solid.SetSaveStiffMat(false);	
		solid.SetStationary(true);
		// Setting Material Parameter
		solid.SetYoungPoisson(2,0.1,true);	// planter stress
		solid.SetGeometricalNonlinear(false);
		solid.SetGravitation(0.0,0.0);
		solid.SetTimeIntegrationParameter(dt,0.7);
		
		unsigned int id_field_bc0;
		{
			std::vector<unsigned int> aIdEA;
			aIdEA.push_back(conv.GetIdEA_fromCad(1,Cad::EDGE));
//			aIdEA.push_back(conv.GetIdEA_fromCad(5,Cad::EDGE));
			id_field_bc0 = solid.AddFixElemAry(aIdEA,world);
		}
		unsigned int id_field_bc1 = solid.AddFixElemAry(conv.GetIdEA_fromCad(8,Cad::EDGE),world);
    field_value_setter = CFieldValueSetter(id_field_bc1,world);
    field_value_setter.SetMathExp("-0.03*(1-cos(t*PI*2*0.1))", 1,Fem::Field::VALUE, world);	// oscilate bc1_field x axis        
		
		id_field_disp = solid.GetIdField_Disp();
		id_field_stress = world.MakeField_FieldElemDim(id_field_disp,2,STSR2,VALUE,BUBBLE);		
//		id_field_equiv_stress = world.MakeField_FieldElemDim(id_field_disp,2,SCALAR,VALUE,BUBBLE);
		
		// set up visualization
		drawer_ary.Clear();		
//		drawer_ary.PushBack( new View::CDrawerFace(id_field_equiv_stress,false,world,id_field_equiv_stress, 0,0.5) );
		drawer_ary.PushBack( new View::CDrawerFace(id_field_disp,false,world) );
		drawer_ary.PushBack( new View::CDrawerVector(id_field_stress,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,false,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,true ,world) );
		drawer_ary.InitTrans(camera);	// init view transformation
	}
	
	iprob++;
	if( iprob == nprob ){
		iprob = 0;
	}
}


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

int main(int argc,char* argv[])
{
	// initialize GLUT
	glutInitWindowPosition(200,200);
	glutInitWindowSize(400, 300);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
	glutCreateWindow("FEM View");

	// set call back functions
	glutDisplayFunc(myGlutDisplay);
	glutReshapeFunc(myGlutResize);
	glutMotionFunc(myGlutMotion);
	glutMouseFunc(myGlutMouse);
	glutKeyboardFunc(myGlutKeyboard);
	glutSpecialFunc(myGlutSpecial);
	glutIdleFunc(myGlutIdle);
	
	SetNewProblem();
	glutMainLoop();
	return 0;
}

