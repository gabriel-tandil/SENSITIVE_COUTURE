
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
#include "delfem/mesh3d.h"

#include "delfem/field_world.h"
#include "delfem/field.h"
#include "delfem/field_value_setter.h"
#include "delfem/drawer_field_face.h"
#include "delfem/drawer_field_edge.h"
#include "delfem/drawer_field_vector.h"

#include "delfem/eqnsys_scalar.h"

using namespace Fem::Field;
using namespace Fem::Ls;

Fem::Field::CFieldWorld world;
Fem::Eqn::CEqn_Scalar3D eqn_scalar;
unsigned int id_base;
double dt = 0.02;
View::CDrawerArrayField drawer_ary;
std::vector<Fem::Field::CFieldValueSetter> field_value_setter_ary;
Com::View::CCamera camera;
double mov_begin_x, mov_begin_y;

bool SetNewProblem()
{
	const unsigned int nprob = 4;
	static unsigned int iprob = 0;
	
	if( iprob == 0 )	// ２次元問題の設定
	{
		////////////////
		Cad::CCadObj2D cad_2d;
 		{	// 形を作る
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,1.0) );
			vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			cad_2d.AddPolygon( vec_ary );
		}
		// メッシュを作る
		Msh::CMesh3D_Extrude msh_3d;
		msh_3d.Extrude( Msh::CMesher2D(cad_2d,0.1),1.0,0.1 );
		world.Clear();
		id_base = world.AddMesh( msh_3d );
        const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);
		// 方程式の設定
		eqn_scalar.SetDomain(id_base,world);
		eqn_scalar.SetAlpha(1.0);
		eqn_scalar.SetTimeIntegrationParameter(dt);
		unsigned int id_val_bc0 = eqn_scalar.AddFixElemAry(conv.GetIdEA_fromCad(1,Cad::EDGE,1),world);
		unsigned int id_val_bc1 = eqn_scalar.AddFixElemAry(conv.GetIdEA_fromCad(2,Cad::EDGE,2),world);
		unsigned int id_val_bc2 = eqn_scalar.AddFixElemAry(conv.GetIdEA_fromCad(1,Cad::EDGE,2),world);    
    field_value_setter_ary.clear();    
    {
      Fem::Field::CFieldValueSetter fvs(id_val_bc0,world);
      fvs.SetMathExp("cos(2*PI*t+0.1)", 0, Fem::Field::VALUE,world);
      field_value_setter_ary.push_back(fvs);
    }
    Fem::Field::SetFieldValue_Constant(id_val_bc1,0,Fem::Field::VALUE,world, 1.0);    
    {
      Fem::Field::CFieldValueSetter fvs(id_val_bc2,world);
      fvs.SetMathExp("sin(2*PI*t+0.1)", 0, Fem::Field::VALUE,world);
      field_value_setter_ary.push_back(fvs);
    }    
		// 描画オブジェクトの登録
		drawer_ary.Clear();
		const unsigned int id_field_val = eqn_scalar.GetIdField_Value();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_val,true,world, id_field_val,-1.0,1.0) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_val,true,world) );
		drawer_ary.InitTrans(camera);
	}
	else if( iprob == 1 )
	{
		Msh::CMesher2D mesh_2d;
		mesh_2d.ReadFromFile_GiDMsh("../input_file/rect_quad.msh");
		Msh::CMesh3D_Extrude mesh_3d;
		mesh_3d.Extrude(mesh_2d, 5.0, 0.5 );
		world.Clear();
		unsigned int id_base = world.AddMesh( mesh_3d );
    const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);
		// 方程式の設定
		eqn_scalar.SetDomain(id_base,world);
		eqn_scalar.SetAlpha(1.0);
		eqn_scalar.SetTimeIntegrationParameter(dt);
		unsigned int id_val_bc0 = eqn_scalar.AddFixElemAry(conv.GetIdEA_fromMshExtrude(1,1),world);
		unsigned int id_val_bc1 = eqn_scalar.AddFixElemAry(conv.GetIdEA_fromMshExtrude(3,2),world);
		unsigned int id_val_bc2 = eqn_scalar.AddFixElemAry(conv.GetIdEA_fromMshExtrude(1,3),world);
    
    field_value_setter_ary.clear();    
    {
      Fem::Field::CFieldValueSetter fvs(id_val_bc0,world);
      fvs.SetMathExp("cos(2*PI*t+0.1)", 0, Fem::Field::VALUE,world);
      field_value_setter_ary.push_back(fvs);
    }
    Fem::Field::SetFieldValue_Constant(id_val_bc1,0,Fem::Field::VALUE,world, 1.0);
    {
      Fem::Field::CFieldValueSetter fvs(id_val_bc2,world);
      fvs.SetMathExp("sin(2*PI*t+0.1)", 0, Fem::Field::VALUE,world);
      field_value_setter_ary.push_back(fvs);
    } 
		// 描画オブジェクトの登録
		drawer_ary.Clear();
		const unsigned int id_field_val = eqn_scalar.GetIdField_Value();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_val,true,world, id_field_val,-1.0,1.0) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_val,true,world) );
		drawer_ary.InitTrans(camera);
	}
	else if( iprob == 2 )
	{
		Msh::CMesher2D mesh_2d;
		mesh_2d.ReadFromFile_GiDMsh("../input_file/hexa_tri.msh");
		Msh::CMesh3D_Extrude mesh_3d;
		mesh_3d.Extrude(mesh_2d, 5.0, 0.5 );
		world.Clear();
		id_base = world.AddMesh( mesh_3d );
        const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);
		// 方程式の設定
		eqn_scalar.SetDomain(id_base,world);
		eqn_scalar.SetAlpha(1.0);
		eqn_scalar.SetTimeIntegrationParameter(dt);
		unsigned int id_val_bc0 = eqn_scalar.AddFixElemAry(conv.GetIdEA_fromMshExtrude(10,2),world);
		unsigned int id_val_bc1 = eqn_scalar.AddFixElemAry(conv.GetIdEA_fromMshExtrude(1,1),world);
		unsigned int id_val_bc2 = eqn_scalar.AddFixElemAry(conv.GetIdEA_fromMshExtrude(11,2),world);        
    field_value_setter_ary.clear();    
    {
      Fem::Field::CFieldValueSetter fvs(id_val_bc0,world);
      fvs.SetMathExp("cos(2*PI*t+0.1)", 0, Fem::Field::VALUE,world);
      field_value_setter_ary.push_back(fvs);
    }
    Fem::Field::SetFieldValue_Constant(id_val_bc1,0,Fem::Field::VALUE,world, 1.0);    
    {
      Fem::Field::CFieldValueSetter fvs(id_val_bc2,world);
      fvs.SetMathExp("sin(2*PI*t+0.1)", 0, Fem::Field::VALUE,world);
      field_value_setter_ary.push_back(fvs);
    } 
		// 描画オブジェクトの登録
		drawer_ary.Clear();
		const unsigned int id_field_val = eqn_scalar.GetIdField_Value();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_val,true,world, id_field_val,-1.0,1.0) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_val,true,world) );
		drawer_ary.InitTrans(camera);
	}
	else if( iprob == 3 )
	{
		Msh::CMesher3D mesh_3d;
		mesh_3d.ReadFromFile_GiDMsh("../input_file/cylinder_tet.msh");
		world.Clear();
		id_base = world.AddMesh( mesh_3d );
		// 方程式の設定
		eqn_scalar.SetDomain(id_base, world);
		eqn_scalar.SetAlpha(1.0);
		eqn_scalar.SetTimeIntegrationParameter(dt);
		unsigned int id_val_bc0 = eqn_scalar.AddFixElemAry(3,world);
		unsigned int id_val_bc1 = eqn_scalar.AddFixElemAry(4,world);    
    field_value_setter_ary.clear();    
    {
      Fem::Field::CFieldValueSetter fvs(id_val_bc0,world);
      fvs.SetMathExp("cos(2*PI*t+0.1)", 0, Fem::Field::VALUE,world);
      field_value_setter_ary.push_back(fvs);
    }
    Fem::Field::SetFieldValue_Constant(id_val_bc1,0,Fem::Field::VALUE,world, 0.0);    
		// 描画オブジェクトの登録
		drawer_ary.Clear();
		const unsigned int id_field_val = eqn_scalar.GetIdField_Value();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_val, true,world, id_field_val,-1.0,1.0) );
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

	if( is_animation ){
		cur_time += dt;
//		world.FieldValueExec(cur_time);
    for(unsigned int iset=0;iset<field_value_setter_ary.size();iset++){
      field_value_setter_ary[iset].ExecuteValue(cur_time,world);
    }    
		eqn_scalar.Solve(world);
		if( eqn_scalar.GetAry_ItrNormRes().size() > 0 ){
			std::cout << "Iter : " << eqn_scalar.GetAry_ItrNormRes()[0].first << " ";
			std::cout << "Res : " << eqn_scalar.GetAry_ItrNormRes()[0].second << std::endl;
		}
		drawer_ary.Update(world);
	}

	drawer_ary.Draw();
	ShowFPS();
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
	
	::glMatrixMode(GL_MODELVIEW);
	::glLoadIdentity();
	Com::View::SetProjectionTransform(camera);
	::glutPostRedisplay();
}

void myGlutIdle(){
	::glutPostRedisplay();
}


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

int main(int argc,char* argv[])
{
	// initialization of glut
	glutInitWindowPosition(200,200);
	glutInitWindowSize(400,300);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
	glutCreateWindow("FEM View");

	// setting of call back function
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
