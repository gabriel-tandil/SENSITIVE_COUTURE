

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
#include "delfem/glut_utility.h"

#include "delfem/cad_obj2d.h"
#include "delfem/mesher2d.h"

#include "delfem/field_world.h"
#include "delfem/field.h"
#include "delfem/field_value_setter.h"
#include "delfem/drawer_field_face.h"
#include "delfem/drawer_field_edge.h"
#include "delfem/drawer_field_vector.h"

#include "delfem/ls/preconditioner.h"
#include "delfem/ls/solver_ls_iter.h"

#include "delfem/femls/linearsystem_field.h"

#include "delfem/femeqn/eqn_dkt.h"

using namespace Fem::Ls;
using namespace Fem::Field;


void myGlutIdle(){	// call back function for idle
	glutPostRedisplay();
}

Com::View::CCamera camera;

void myGlutResize(int w, int h){ // call vack function for resize
	camera.SetWindowAspect((double)w/h);
	glViewport(0, 0, w, h);
	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
	Com::View::SetProjectionTransform(camera);
	glutPostRedisplay();
}

View::CDrawerArrayField drawer_ary;
Fem::Field::CFieldWorld world;
double mov_begin_x, mov_begin_y;
bool is_animation = true;

void myGlutMotion( int x, int y ){
	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
	const int win_w = viewport[2];
	const int win_h = viewport[3];
	const double mov_end_x = (2.0*x-win_w)/win_w;
	const double mov_end_y = (win_h-2.0*y)/win_h;
	camera.MouseRotation(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 
//	camera.MousePan(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 
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

void myGlutDisplay(void){	// call back function for display
	::glClearColor(1.0, 1.0, 1.0, 1.0);
	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	::glEnable(GL_DEPTH_TEST);

	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 1.1f, 4.0f );

	::glMatrixMode(GL_MODELVIEW);
	::glLoadIdentity();
	Com::View::SetModelViewTransform(camera);

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
		::glMatrixMode(GL_PROJECTION);
		::glLoadIdentity();
		Com::View::SetProjectionTransform(camera);
		break;
	}
	::glutPostRedisplay();
}


void SetNewProblem()
{
	const unsigned int nprob = 2;
	static unsigned int iprob = 0;

	static int id_val_bc0=0, id_val_bc1=0, id_val_bc2=0;
	
	if( iprob == 0 ){	// ２次元問題の設定
		Cad::CCadObj2D cad_2d;
		{	// 形を作る
            std::vector<Com::CVector2D> vec_ary;
            vec_ary.push_back( Com::CVector2D(-0.5,-0.5) );
            vec_ary.push_back( Com::CVector2D( 0.5,-0.5) );
            vec_ary.push_back( Com::CVector2D( 0.5, 0.5) );
            vec_ary.push_back( Com::CVector2D(-0.5, 0.5) );
			cad_2d.AddPolygon( vec_ary );
            cad_2d.AddVertex(Cad::LOOP,1,Com::CVector2D(0.0,0.0));
		}
		world.Clear();
		const unsigned int id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.05) );
		const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);

		unsigned int id_field_rot = world.MakeField_FieldElemDim(id_base,2,VECTOR2,VALUE,CORNER);
		unsigned int id_field_deflect = world.MakeField_FieldElemDim(id_base,2,SCALAR,VALUE,CORNER);

		Fem::Ls::CLinearSystem_Field ls;
		LsSol::CPreconditioner_ILU prec;
		ls.AddPattern_Field(id_field_deflect,world);
		ls.AddPattern_Field(id_field_rot,id_field_deflect,world);
		prec.SetLinearSystem(ls.m_ls);

		unsigned int id_field_rot_fix0 = world.GetPartialField(id_field_rot,conv.GetIdEA_fromCad(2,Cad::EDGE));
		unsigned int id_field_rot_fix1 = world.GetPartialField(id_field_rot,conv.GetIdEA_fromCad(4,Cad::EDGE));
		ls.SetFixedBoundaryCondition_Field(id_field_rot_fix0,world);
		ls.SetFixedBoundaryCondition_Field(id_field_rot_fix1,world);

		unsigned int id_field_def_fix0 = world.GetPartialField(id_field_deflect,conv.GetIdEA_fromCad(2,Cad::EDGE));
		unsigned int id_field_def_fix1 = world.GetPartialField(id_field_deflect,conv.GetIdEA_fromCad(4,Cad::EDGE));
    Fem::Field::SetFieldValue_Constant(id_field_rot_fix1,1,Fem::Field::VALUE,world, -1);
    
		ls.SetFixedBoundaryCondition_Field(id_field_def_fix0,world);
		ls.SetFixedBoundaryCondition_Field(id_field_def_fix1,world);    
    
		ls.InitializeMarge();
		Fem::Eqn::AddLinearSystem_DKT2D_Static(ls,world,id_field_deflect,id_field_rot);
		double res = ls.FinalizeMarge();
		prec.SetValue(ls.m_ls);

		std::cout << "Residual : " << res << std::endl;
		{
			double tol = 1.0e-6;
			unsigned int iter = 10000;
			LsSol::Solve_CG(tol,iter,ls);
//			Fem::Sol::Solve_PCG(tol,iter,ls,prec);
			std::cout << iter << " " << tol << std::endl;
		}
		ls.UpdateValueOfField(id_field_deflect,world,VALUE);
		ls.UpdateValueOfField(id_field_rot,world,VALUE);

		drawer_ary.Clear();
	//	drawer_ary.PushBack( new View::CDrawerVector(id_field_rot,world) );
		drawer_ary.PushBack( new View::CDrawerFace(id_field_deflect,false,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_deflect,false,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_deflect,true,world) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 1 ){
		Msh::CMesher2D msh;
		msh.ReadFromFile_GiDMsh("../input_file/rect_tri.msh");
		world.Clear();
		const unsigned int id_base = world.AddMesh( msh );

		unsigned int id_field_rot = world.MakeField_FieldElemDim(id_base,2,  VECTOR2,VALUE,CORNER);
		unsigned int id_field_deflect = world.MakeField_FieldElemDim(id_base,2,  SCALAR,VALUE,CORNER);

		unsigned int id_field_rot_fix0 = world.GetPartialField(id_field_rot,3);
		unsigned int id_field_rot_fix1 = world.GetPartialField(id_field_rot,4);
		unsigned int id_field_def_fix0 = world.GetPartialField(id_field_deflect,3);
		unsigned int id_field_def_fix1 = world.GetPartialField(id_field_deflect,4);
    
    Fem::Field::SetFieldValue_Constant(id_field_def_fix1,0,Fem::Field::VALUE,world, -3);    

		Fem::Ls::CLinearSystem_Field ls;
		LsSol::CPreconditioner_ILU prec;
		ls.AddPattern_Field(id_field_deflect,world);
		ls.AddPattern_Field(id_field_rot,id_field_deflect,world);
		prec.SetLinearSystem(ls.m_ls);

		ls.SetFixedBoundaryCondition_Field(id_field_rot_fix0,world);
		ls.SetFixedBoundaryCondition_Field(id_field_rot_fix1,world);
		ls.SetFixedBoundaryCondition_Field(id_field_def_fix0,world);
		ls.SetFixedBoundaryCondition_Field(id_field_def_fix1,world);
		ls.InitializeMarge();
		Fem::Eqn::AddLinearSystem_DKT2D_Static(ls,world,id_field_deflect,id_field_rot);
		double res = ls.FinalizeMarge();
		prec.SetValue(ls.m_ls);

		std::cout << "Residual : " << res << std::endl;
		{
			double tol = 1.0e-6;
			unsigned int iter = 10000;
			LsSol::Solve_CG(tol,iter,ls);
//			Fem::Sol::Solve_PCG(tol,iter,ls,prec);
			std::cout << iter << " " << tol << std::endl;
		}
		ls.UpdateValueOfField(id_field_deflect,world,VALUE);
		ls.UpdateValueOfField(id_field_rot,world,VALUE);

		drawer_ary.Clear();
	//	drawer_ary.PushBack( new View::CDrawerVector(id_field_rot,world) );
		drawer_ary.PushBack( new View::CDrawerFace(id_field_deflect,false,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_deflect,false,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_deflect,true,world) );
		drawer_ary.InitTrans( camera );
	}

	iprob++;
	if( iprob == nprob ) iprob=0;
}

int main(int argc,char* argv[])
{
	SetNewProblem();

	// GLUTの初期設定
	glutInitWindowPosition(200,200);
	glutInitWindowSize(400, 300);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
	glutCreateWindow("DelFEM Demo");

	// コールバック関数の設定
	glutIdleFunc(myGlutIdle);
	glutKeyboardFunc(myGlutKeyboard);
	glutDisplayFunc(myGlutDisplay);
	glutReshapeFunc(myGlutResize);
	glutSpecialFunc(myGlutSpecial);;
	glutMotionFunc(myGlutMotion);
	glutMouseFunc(myGlutMouse);
	
	// メインループ
	glutMainLoop();
	return 0;
}
