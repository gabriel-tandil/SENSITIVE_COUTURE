
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

#include "delfem/camera.h"
#include "delfem/glut_utility.h"  // FPS()

#include "delfem/cad_obj2d.h"
#include "delfem/mesher2d.h"

#include "delfem/femls/zlinearsystem.h"
#include "delfem/femls/zsolver_ls_iter.h"
#include "delfem/field_world.h"
#include "delfem/field_value_setter.h"
#include "delfem/drawer_field.h"
#include "delfem/drawer_field_face.h"
#include "delfem/drawer_field_edge.h"
#include "delfem/drawer_field_vector.h"
#include "delfem/femeqn/eqn_helmholtz.h"

using namespace Fem::Field;
using namespace Fem::Ls;

Fem::Field::CFieldWorld world;
double dt = 0.02;
View::CDrawerArrayField drawer_ary;
Com::View::CCamera camera;
double mov_begin_x, mov_begin_y;

bool SetNewProblem()
{
	const unsigned int nprob = 1;
	static unsigned int iprob = 0;
	
	if( iprob == 0 )
  {
		////////////////
		Cad::CCadObj2D cad_2d;
		unsigned int id_v;
 		{	// define shape
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(2.0,0.0) );
			vec_ary.push_back( Com::CVector2D(2.0,2.0) );
			vec_ary.push_back( Com::CVector2D(0.0,2.0) );
			const unsigned int id_l = cad_2d.AddPolygon(vec_ary).id_l_add;
			id_v = cad_2d.AddVertex(Cad::LOOP, id_l, Com::CVector2D(0.5,0.05) ).id_v_add;
		}
    
		world.Clear();
		const unsigned int id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.04) );
		Fem::Field::CIDConvEAMshCad conv = world.GetIDConverter(id_base);
		const unsigned int id_field_val = world.MakeField_FieldElemDim(id_base,2,ZSCALAR,VALUE,CORNER);
//		unsigned int id_field_bc0 = world.GetPartialField(id_field_val,conv.GetIdEA_fromCad(2,1));
		unsigned int id_field_bc1;
		{
			std::vector<unsigned int> aEA;
			aEA.push_back( conv.GetIdEA_fromCad(1,Cad::EDGE) );
			aEA.push_back( conv.GetIdEA_fromCad(2,Cad::EDGE) );
			aEA.push_back( conv.GetIdEA_fromCad(3,Cad::EDGE) );
			aEA.push_back( conv.GetIdEA_fromCad(4,Cad::EDGE) );
			id_field_bc1 = world.GetPartialField(id_field_val,aEA);
		}
    Fem::Field::SetFieldValue_Constant(id_field_bc1,0,Fem::Field::VALUE,world,0);

		CZLinearSystem ls;
		CZPreconditioner_ILU prec;
		ls.AddPattern_Field(id_field_val,world);
//		ls.SetFixedBoundaryCondition_Field(id_field_bc0,world);
//		ls.SetFixedBoundaryCondition_Field(id_field_bc1,world);
		prec.SetFillInLevel(1);
		prec.SetLinearSystem(ls);

		double wave_length = 0.4;
		ls.InitializeMarge();
		Fem::Eqn::AddLinSys_Helmholtz(ls,wave_length,world,id_field_val);
		Fem::Eqn::AddLinSys_SommerfeltRadiationBC(ls,wave_length,world,id_field_bc1);
		double res = ls.FinalizeMarge();
		prec.SetValue(ls);

		{
			const unsigned int id_ea_v = conv.GetIdEA_fromCad(id_v,Cad::VERTEX);
			std::cout << id_ea_v << std::endl;
			const CElemAry& ea = world.GetEA(id_ea_v);
			const CElemAry::CElemSeg& es = ea.GetSeg(1);
			assert( ea.ElemType() == Fem::Field::POINT );
			unsigned int noes[1];
			es.GetNodes(0,noes);
			std::cout << noes[0] << std::endl;
			ls.GetResidualPtr(id_field_val,CORNER,world)->AddValue(noes[0],0,Com::Complex(1,0));
		}

		std::cout << "Residual : " << res << std::endl;
		{
			double tol = 1.0e-6;
			unsigned int iter = 2000;
//			Fem::Ls::Solve_CG(tol,iter,ls);
//			Fem::Ls::Solve_PCG(tol,iter,ls,prec);
			Fem::Ls::Solve_PCOCG(tol,iter,ls,prec);
//			Fem::Ls::Solve_CGNR(tol,iter,ls);
//			Fem::Ls::Solve_BiCGSTAB(tol,iter,ls);
//			Fem::Ls::Solve_BiCGStabP(tol,iter,ls,prec);
			std::cout << iter << " " << tol << std::endl;
		}
		ls.UpdateValueOfField(id_field_val,world,VALUE);

		drawer_ary.Clear();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_val,true,world, id_field_val,-0.05,0.05) );
//		drawer_ary.PushBack( new View::CDrawerFaceContour(id_field_val,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_val,true,world) );
		drawer_ary.InitTrans(camera);	
	}
	iprob++;
	if( iprob == nprob ) iprob=0;

	return true;
}


////////////////////////////////////////////////////////////////

double cur_time = 0.0;
bool is_animation = true;

void myGlutResize(int w, int h)
{
	camera.SetWindowAspect((double)w/h);
	::glViewport(0, 0, w, h);
	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
  SetProjectionTransform(camera);
	::glutPostRedisplay();
}

void myGlutDisplay(void)
{
	::glClearColor(1.0, 1.0, 1.0 ,1.0);
	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	::glEnable(GL_DEPTH_TEST);

	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 1.1f, 4.0f );

  ::glMatrixMode(GL_MODELVIEW);
	::glLoadIdentity();
	SetModelViewTransform(camera);

//	ShowFPS();
	drawer_ary.Draw();
	glutSwapBuffers();
}

void myGlutMotion( int x, int y ){
	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
	const int win_w = viewport[2];
	const int win_h = viewport[3];
	const double mov_end_x = (2.0*x-win_w)/win_w;
	const double mov_end_y = (win_h-2.0*y)/win_h;
	camera.MouseRotation(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 
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
      ::glMatrixMode(GL_PROJECTION);
      ::glLoadIdentity();
      SetProjectionTransform(camera);
		break;
	}
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
  SetProjectionTransform(camera);
	::glutPostRedisplay();
}

void myGlutIdle(){
	::glutPostRedisplay();
}


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

int main(int argc,char* argv[])
{

	SetNewProblem();

	glutInitWindowPosition(200,200);
	glutInitWindowSize(250, 250);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
	glutCreateWindow("FEM View");

	glutDisplayFunc(myGlutDisplay);
	glutReshapeFunc(myGlutResize);
	glutMotionFunc(myGlutMotion);
	glutMouseFunc(myGlutMouse);
	glutKeyboardFunc(myGlutKeyboard);
	glutSpecialFunc(myGlutSpecial);
	glutIdleFunc(myGlutIdle);

	glutMainLoop();
	return 0;
}
