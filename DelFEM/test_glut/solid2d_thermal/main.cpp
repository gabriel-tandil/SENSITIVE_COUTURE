////////////////////////////////////////////////////////////////
//                                                            //
//		Solid2d thermal                                       //
//                                                            //
//          Copy Rights (c) Nobuyuki Umetani 2007             //
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

#if defined(__APPLE__) && defined(__MACH__)
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif

#include "delfem/camera.h"
#include "delfem/glut_utility.h"

#include "delfem/cad_obj2d.h"
#include "delfem/mesh3d.h"
#include "delfem/mesher2d.h"
#include "delfem/field_world.h"
#include "delfem/field.h"
#include "delfem/field_value_setter.h"
#include "delfem/drawer_field_face.h"
#include "delfem/drawer_field_edge.h"
#include "delfem/eqnsys_solid.h"
#include "delfem/eqnsys_scalar.h"

using namespace Fem::Ls;
using namespace Fem::Field;

Com::View::CCamera camera;
double mov_begin_x, mov_begin_y;
bool is_animation = true;



// ÉäÉTÉCÉYéûÇÃÉRÅ[ÉãÉoÉbÉNä÷êî
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
		exit(0);  /* '\033' ? ESC ? ASCII ??? */
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

void myGlutIdle(){	// idle call back function
	::glutPostRedisplay();
}

////////////////////////////////


Fem::Field::CFieldWorld world;
View::CDrawerArrayField drawer_ary;
Fem::Field::CFieldValueSetter field_value_setter;
double cur_time = 0.0;
double dt = 0.05;
Fem::Eqn::CEqnSystem_Solid2D solid;
Fem::Eqn::CEqnSystem_Scalar2D scalar;

void myGlutDisplay(void)	// display call back function
{
//	::glClearColor(0.2f, 0.7f, 0.7f,1.0f);
	::glClearColor(1.0f, 1.0f, 1.0f,1.0f);
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
		scalar.Solve(world);
		solid.Solve(world);
//		world.FieldValueDependExec();
		drawer_ary.Update(world);
	}

	drawer_ary.Draw();
	ShowFPS();
	::glutSwapBuffers();
}

void SetNewProblem()
{
	Cad::CCadObj2D cad_2d;
	{	// å`ÇçÏÇÈ
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
	unsigned int id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.1) );
	const CIDConvEAMshCad& conv = world.GetIDConverter(id_base);
	// îMägéUï˚íˆéÆÇ…ê›íË
	scalar.SetDomain_Field(id_base,world);
	scalar.SetCapacity(5);
	scalar.SetAlpha(1.0);
	scalar.SetStationary(false);
	scalar.SetTimeIntegrationParameter(dt);
	const unsigned int id_field_tmp_fix0 = scalar.AddFixElemAry( conv.GetIdEA_fromCad(4,Cad::EDGE),world);
	{	// îMÇÃéûä‘óöóÇó^Ç¶ÇÈ
//		Fem::Field::CField& field_fix = world.GetField(id_field_tmp_fix0);
//		field_fix.SetValue("floor(0.5*sin(t)+1.0)",0,Fem::Field::VALUE,world,true);
    field_value_setter = Fem::Field::CFieldValueSetter(id_field_tmp_fix0,world);
    field_value_setter.SetMathExp("floor(0.5*sin(t)+1.0)",0,Fem::Field::VALUE,world);
	}

	solid.UpdateDomain_Field(id_base,world);
	solid.SetRho(0.01);
	solid.SetYoungPoisson(10.0,0.0,true);	// ÉÑÉìÉOó¶Ç∆É|ÉAÉ\Éìî‰ÇÃê›íË(ïΩñ âûóÕ)
	solid.SetSaveStiffMat(false);
	solid.SetStationary(true);
	solid.SetTimeIntegrationParameter(dt);
	// îMâûóÕñ‚ëËÇ…ê›íË
	const unsigned int id_field_temp = scalar.GetIdField_Value();
	solid.SetThermalStress(id_field_temp);

	// ó¿ÇÃï–ï˚Çå≈íË
	solid.AddFixElemAry(7,world);

	drawer_ary.Clear();
	const unsigned int id_field_disp = solid.GetIdField_Disp();
  Fem::Field::SetFieldValue_Constant(id_field_disp,0,Fem::Field::VALUE,world,0);
  Fem::Field::SetFieldValue_Constant(id_field_disp,1,Fem::Field::VALUE,world,0);  
	drawer_ary.PushBack( new View::CDrawerFace(id_field_temp,true,world, id_field_temp,-1,1) );
	drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,false,world) );
	drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,true ,world) );

	// Viewç¿ïWïœä∑ÇÃê›íË
	drawer_ary.InitTrans(camera);
}


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

int main(int argc,char* argv[])
{
	// initialization of glut
	glutInitWindowPosition(200,200);
	glutInitWindowSize(400, 300);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
	glutCreateWindow("FEM View");

	// setting call back function
	glutDisplayFunc(myGlutDisplay);
	glutReshapeFunc(myGlutResize);
	glutMotionFunc(myGlutMotion);
	glutMouseFunc(myGlutMouse);
	glutKeyboardFunc(myGlutKeyboard);
	glutSpecialFunc(myGlutSpecial);
	glutIdleFunc(myGlutIdle);

	SetNewProblem();	// setting problem

	glutMainLoop();	// main loop
	return 0;
}
