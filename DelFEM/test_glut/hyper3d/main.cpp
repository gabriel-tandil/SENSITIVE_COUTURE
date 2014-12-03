////////////////////////////////////////////////////////////////
//                                                            //
//		Test of Hyper Elastic Solid 3D                        //
//                                                            //
//          Copy Rights (c) Nobuyuki Umetani 2008             //
//          e-mail : numetani@gmail.com                       //
////////////////////////////////////////////////////////////////

#if defined(__VISUALC__)
#pragma warning ( disable : 4786 )
#endif
#define for if(0);else for

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
#include "delfem/glut_utility.h"  // FPS()

#include "delfem/cad_obj2d.h"
#include "delfem/mesh_primitive.h"

#include "delfem/ls/preconditioner.h"
#include "delfem/ls/solver_ls_iter.h"
#include "delfem/femls/linearsystem_field.h"
#include "delfem/femeqn/eqn_linear_solid2d.h"
#include "delfem/femeqn/eqn_linear_solid3d.h"
#include "delfem/femeqn/eqn_hyper.h"
#include "delfem/field.h"
#include "delfem/field_world.h"
#include "delfem/field_value_setter.h"
#include "delfem/drawer_field_face.h"
#include "delfem/drawer_field_edge.h"

using namespace Fem::Ls;
using namespace Fem::Field;

Com::View::CCamera camera;
double mov_begin_x, mov_begin_y;
int imodifier;
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
  imodifier = glutGetModifiers();
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
		::glMatrixMode(GL_PROJECTION);
		::glLoadIdentity();
		Com::View::SetProjectionTransform(camera);
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
	Com::View::SetProjectionTransform(camera);
	::glutPostRedisplay();
}

void myGlutIdle(){
	::glutPostRedisplay();
}

////////////////////////////////


Fem::Field::CFieldWorld world;
Fem::Field::CFieldValueSetter field_value_setter;
View::CDrawerArrayField drawer_ary;
double cur_time = 0.0;
const double dt = 0.06;
const double newmarkb_gamma = 0.59;
const double newmarkb_beta = 0.25*(0.5+newmarkb_gamma)*(0.5+newmarkb_gamma);
unsigned int id_field_disp;
unsigned int id_field_lambda;
Fem::Ls::CLinearSystem_Field ls;
LsSol::CPreconditioner_ILU prec;

void StepTime(){
  const double c1 = 200;
  const double c2 = 200;
  const double rho = 1.8;
  const double g[3] = { 0,0,0 };
  for(unsigned int iitr=0;iitr<2;iitr++){
    ls.InitializeMarge();
    Fem::Eqn::AddLinSys_Hyper3D_NonStatic_NewmarkBeta
    (dt, newmarkb_gamma, newmarkb_beta, ls,
     c1, c2,
     rho, g[0], g[1], g[2],    
     id_field_disp, id_field_lambda, world, 
     iitr==0);
    /*        Fem::Eqn::AddLinSys_LinearSolid3D_NonStatic_NewmarkBeta(
     dt, newmarkb_gamma, newmarkb_beta, ls,
     lambda, myu,
     rho, g[0], g[1], g[2],
     world,
     id_field_disp );*/
    const double res = ls.FinalizeMarge();
    prec.SetValue(ls.m_ls);
    {
      double conv_ratio = 1.0e-6;
      unsigned int iteration = 400;
      LsSol::CLinearSystemPreconditioner lsp(ls.m_ls,prec);
      //            Sol::Solve_PCG(conv_ratio,iteration, lsp);
      LsSol::Solve_PBiCGSTAB(conv_ratio,iteration, lsp);
      std::cout << "Iter NR : " << iitr << "   Res : " << res << "   Iter : " << iteration << "   Conv : " << conv_ratio << std::endl;
    }
    ls.UpdateValueOfField_NewmarkBeta(newmarkb_gamma,newmarkb_beta,dt, id_field_disp,  world,iitr==0);
    ls.UpdateValueOfField_NewmarkBeta(newmarkb_gamma,newmarkb_beta,dt, id_field_lambda,world,iitr==0);
  }  
}

void myGlutDisplay(void)
{
	::glClearColor(0.2f, 0.7f, 0.7f,1.0f);
	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	::glEnable(GL_DEPTH_TEST);

	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 1.1f, 4.0f );

	::glMatrixMode(GL_MODELVIEW);
	::glLoadIdentity();
	Com::View::SetModelViewTransform(camera);
    
	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
    Com::View::SetProjectionTransform(camera);

	if( is_animation ){
		cur_time += dt;
//		world.FieldValueExec(cur_time);
    field_value_setter.ExecuteValue(cur_time,world);    
    StepTime();
		drawer_ary.Update(world);
	}

  ShowBackGround();
	drawer_ary.Draw();
	ShowFPS();



	::glutSwapBuffers();
}

void SetNewProblem()
{
	const unsigned int nprob = 1;	// –â‘è”
	static unsigned int iprob = 0;

	if( iprob == 0 ){
    cur_time = 0;
    Msh::CMesh_Primitive_Hexahedra mesh_3d(0.5, 4, 6,  1,8,8);
    //        Msh::CMesh_Primitive_ThickCylinder mesh_3d(2,2.5,20, 1, 12, 6);
    world.Clear();
		const unsigned int id_base = world.AddMesh( mesh_3d );
    const CIDConvEAMshCad& conv = world.GetIDConverter(id_base);
    id_field_disp   = world.MakeField_FieldElemDim(id_base,3,VECTOR3,VALUE|VELOCITY|ACCELERATION,CORNER);
    id_field_lambda = world.MakeField_FieldElemDim(id_base,3,SCALAR, VALUE|VELOCITY|ACCELERATION,BUBBLE);
		unsigned int id_field_bc1 = world.GetPartialField(id_field_disp,conv.GetIdEA_fromMsh(2));
    field_value_setter = Fem::Field::CFieldValueSetter(id_field_bc1,world);
    field_value_setter.SetMathExp("1*sin(2*t)", 0,Fem::Field::VALUE, world);
    field_value_setter.SetMathExp("1*sin(t  )", 0,Fem::Field::VALUE, world);    
    ls.Clear();
    ls.AddPattern_Field(id_field_disp,world);
    ls.AddPattern_Field(id_field_lambda,id_field_disp,world);
    ls.SetFixedBoundaryCondition_Field(id_field_bc1,world);
    //        ls.SetFixedBoundaryCondition_Field(id_field_bc2,world);
    
    prec.Clear();
    prec.SetFillInLevel(0);
    prec.SetLinearSystem(ls.m_ls);
    
		// •`‰æƒIƒuƒWƒFƒNƒg‚Ì“o˜^
		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_disp,false,world) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerEdge(id_field_disp,false,world) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerEdge(id_field_disp,true ,world) );
		drawer_ary.InitTrans(camera);    
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
	// Initialization of GLUT
	glutInitWindowPosition(200,200);
	glutInitWindowSize(400, 300);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
	glutCreateWindow("FEM View");

	// Setting of call back function
	glutDisplayFunc(myGlutDisplay);
	glutReshapeFunc(myGlutResize);
	glutMotionFunc(myGlutMotion);
	glutMouseFunc(myGlutMouse);
	glutKeyboardFunc(myGlutKeyboard);
	glutSpecialFunc(myGlutSpecial);
	glutIdleFunc(myGlutIdle);
	
	SetNewProblem();	// set new problem
	glutMainLoop();		// enter main loop
	return 0;
}
