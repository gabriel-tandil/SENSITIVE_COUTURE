
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
#include "delfem/drawer_gl_utility.h"
#include "delfem/glut_utility.h"

#include "delfem/cad_obj2d.h" // CCadObj2D
#include "delfem/mesher2d.h"  // CMesher2D

#include "delfem/field.h"	// 
#include "delfem/field_world.h"		// 
#include "delfem/field_value_setter.h"
#include "delfem/drawer_field.h"
#include "delfem/drawer_field_face.h"
#include "delfem/drawer_field_edge.h"
#include "delfem/drawer_field_vector.h"

#include "delfem/eqnsys_scalar.h"

using namespace Fem::Field;
using namespace Fem::Ls;

Fem::Field::CFieldWorld world;
Fem::Eqn::CEqnSystem_Scalar2D eqn_scalar;
double dt = 0.001;
View::CDrawerArrayField drawer_ary;
std::vector<Fem::Field::CFieldValueSetter> field_value_setter_ary;
Com::View::CCamera camera;
double mov_begin_x, mov_begin_y;
unsigned int id_base;


bool SetNewProblem()
{
	const unsigned int nprob = 20;
	static unsigned int iprob = 0;

	static int id_val_bc0=0, id_val_bc1=0, id_val_bc2=0;
	
	if( iprob == 0 )	// ２次元問題の設定
	{
		world.Clear();
		drawer_ary.Clear();
		////////////////
		Cad::CCadObj2D cad_2d;
 		{	// 形を作る
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,1.0) );
			vec_ary.push_back( Com::CVector2D(0.0,1.0) );
			const unsigned int id_l = cad_2d.AddPolygon( vec_ary ).id_l_add;
			const unsigned int id_v1 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.7,0.5)).id_v_add;
			const unsigned int id_v2 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.7,0.9)).id_v_add;
			const unsigned int id_v3 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.8,0.9)).id_v_add;
			const unsigned int id_v4 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(0.8,0.5)).id_v_add;
			cad_2d.ConnectVertex_Line(id_v1,id_v2);
			cad_2d.ConnectVertex_Line(id_v2,id_v3);
			cad_2d.ConnectVertex_Line(id_v3,id_v4);
			cad_2d.ConnectVertex_Line(id_v4,id_v1);
		}
		// メッシュを作る
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.02) );
		Fem::Field::CIDConvEAMshCad conv = world.GetIDConverter(id_base);
		eqn_scalar.SetDomain_Field(id_base,world);    // 方程式の設定
		dt = 0.02;
		eqn_scalar.SetTimeIntegrationParameter(dt);
		eqn_scalar.SetSaveStiffMat(false);
		eqn_scalar.SetStationary(false);
		eqn_scalar.SetAxialSymmetry(false);
		// 全体の方程式の係数設定
		eqn_scalar.SetAlpha(1.0);
		eqn_scalar.SetCapacity(30.0);
		eqn_scalar.SetAdvection(0);

    id_val_bc0 = eqn_scalar.AddFixElemAry(conv.GetIdEA_fromCad(2,Cad::LOOP),world);
		id_val_bc1 = eqn_scalar.AddFixElemAry(conv.GetIdEA_fromCad(1,Cad::EDGE),world);
		id_val_bc2 = eqn_scalar.AddFixElemAry(conv.GetIdEA_fromCad(3,Cad::EDGE),world);
    
    field_value_setter_ary.clear();
    {
      Fem::Field::CFieldValueSetter fvs(id_val_bc0,world);
      fvs.SetMathExp("floor(1+0.8*cos(2*PI*t+0.1))",0,Fem::Field::VALUE,world);
      field_value_setter_ary.push_back(fvs);
    }
    Fem::Field::SetFieldValue_Constant(id_val_bc1,0,Fem::Field::VALUE,world,+1.0);
    Fem::Field::SetFieldValue_Constant(id_val_bc2,0,Fem::Field::VALUE,world,-1.0);    

		// 描画オブジェクトの登録
		const unsigned int id_field_val = eqn_scalar.GetIdField_Value();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_val,true,world,id_field_val,-1,1) );
//		drawer_ary.PushBack( new View::CDrawerFaceContour(id_field_val,world) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_val,true,world) );
		drawer_ary.InitTrans(camera);	// 視線座標変換行列の初期化
	}
	else if( iprob == 1 ){
		eqn_scalar.SetCapacity(10);
	}
	else if( iprob == 2 ){
		eqn_scalar.SetCapacity(5);
	}
	else if( iprob == 3 ){
		eqn_scalar.SetCapacity(1);
	}
	else if( iprob == 4 ){
		eqn_scalar.SetStationary(true);
	}
	else if( iprob == 5 ){
		eqn_scalar.SetSaveStiffMat(true);
	}
	else if( iprob == 6 ){
		eqn_scalar.SetStationary(false);
		eqn_scalar.SetCapacity(10);
		dt = 0.02;
		eqn_scalar.SetTimeIntegrationParameter(dt);
	}
	else if( iprob == 7 ){
		const unsigned int id_field_velo = world.MakeField_FieldElemDim(id_base,2,Fem::Field::VECTOR2,VELOCITY);
		std::cout << "Velo : " << id_field_velo << std::endl;
    Fem::Field::SetFieldValue_MathExp(id_field_velo,0,Fem::Field::VELOCITY,world," (y-0.5)",0);    
    Fem::Field::SetFieldValue_MathExp(id_field_velo,1,Fem::Field::VELOCITY,world,"-(x-0.5)",0);    
    field_value_setter_ary.clear();
		{	// 固定境界条件の設定
      Fem::Field::CFieldValueSetter fvs(id_val_bc0,world);
      fvs.SetMathExp("floor(1+0.8*cos(3*t))",0,Fem::Field::VALUE,world);
      field_value_setter_ary.push_back(fvs);      
		}
		{	// 周囲の固定境界条件の設定
			const CIDConvEAMshCad conv = world.GetIDConverter(id_base);
			std::vector<unsigned int> m_aIDEA;
			m_aIDEA.push_back(conv.GetIdEA_fromCad(1,Cad::EDGE));
			m_aIDEA.push_back(conv.GetIdEA_fromCad(2,Cad::EDGE));
			m_aIDEA.push_back(conv.GetIdEA_fromCad(3,Cad::EDGE));
			m_aIDEA.push_back(conv.GetIdEA_fromCad(4,Cad::EDGE));
			id_val_bc1 = eqn_scalar.AddFixElemAry(m_aIDEA,world);
    }
    Fem::Field::SetFieldValue_Constant(id_val_bc1,0,Fem::Field::VALUE,world,-1.0);
    
		eqn_scalar.SetSaveStiffMat(false);
		eqn_scalar.SetStationary(false);
		eqn_scalar.SetAxialSymmetry(false);
		eqn_scalar.SetTimeIntegrationParameter(dt);
		// 方程式の係数の設定
		eqn_scalar.SetAlpha(0.00001);
		eqn_scalar.SetCapacity(1.0);
		eqn_scalar.SetAdvection(id_field_velo);
	}
	else if( iprob == 8 ){
		eqn_scalar.SetSaveStiffMat(true);
	}
	else if( iprob == 9 ){
		eqn_scalar.SetStationary(true);
	}
	else if( iprob == 10 ){
		eqn_scalar.SetSaveStiffMat(false);
	}
	else if( iprob == 11 )
	{
		////////////////
		Cad::CCadObj2D cad_2d;
 		{	// 形を作る
			{
				std::vector<Com::CVector2D> vec_ary;
				vec_ary.push_back( Com::CVector2D(0.0,0.0) );
				vec_ary.push_back( Com::CVector2D(1.0,0.0) );
				vec_ary.push_back( Com::CVector2D(1.0,1.0) );
				vec_ary.push_back( Com::CVector2D(0.0,1.0) );
				cad_2d.AddPolygon( vec_ary );
			}
			const unsigned int id_v1 = cad_2d.AddVertex(Cad::EDGE,1,Com::CVector2D(0.5,0.0)).id_v_add;
			const unsigned int id_v2 = cad_2d.AddVertex(Cad::EDGE,3,Com::CVector2D(0.5,1.0)).id_v_add;
			cad_2d.ConnectVertex_Line(id_v1,id_v2);
		}
		// メッシュを作る
		world.Clear();
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.05) );	// メッシュで表される場のハンドルを得る
		const CIDConvEAMshCad& conv = world.GetIDConverter(id_base);
		eqn_scalar.SetDomain_FieldElemAry(id_base,conv.GetIdEA_fromCad(2,Cad::LOOP),world);
		eqn_scalar.SetSaveStiffMat(false);
		eqn_scalar.SetStationary(true);
		eqn_scalar.SetTimeIntegrationParameter(dt);
		// 方程式の設定
		dt = 0.02;
		eqn_scalar.SetTimeIntegrationParameter(dt);
		eqn_scalar.SetAlpha(1.0);
		eqn_scalar.SetCapacity(30.0);
		eqn_scalar.SetAdvection(false);
		// 境界条件の設定
		id_val_bc0 = eqn_scalar.AddFixElemAry(conv.GetIdEA_fromCad(2,Cad::EDGE),world);
    field_value_setter_ary.clear();
		{
      Fem::Field::CFieldValueSetter fvs(id_val_bc0,world);
      fvs.SetMathExp("cos(2*PI*t+0.1)",0,Fem::Field::VALUE,world);
      field_value_setter_ary.push_back(fvs);            
		}
		id_val_bc1 = eqn_scalar.AddFixElemAry(conv.GetIdEA_fromCad(3,Cad::EDGE),world);
    Fem::Field::SetFieldValue_Constant(id_val_bc1,0,Fem::Field::VALUE,world,+1.0);    

		// 描画オブジェクトの登録
		drawer_ary.Clear();
		const unsigned int id_field_val = eqn_scalar.GetIdField_Value();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_val,true,world,id_field_val,-1.0,1.0) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_base,true,world) );
		drawer_ary.InitTrans(camera);
	}
	else if( iprob == 12 ){
		eqn_scalar.SetSaveStiffMat(true);
	}
	else if( iprob == 13 ){
		eqn_scalar.SetSaveStiffMat(false);
		eqn_scalar.SetStationary(false);
	}
	else if( iprob == 14 ){
		eqn_scalar.SetSaveStiffMat(true);
	}
	else if( iprob == 15 ){		
		Cad::CCadObj2D cad_2d;
		{	// 正方形に矩形の穴
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
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.02) );
		const CIDConvEAMshCad conv = world.GetIDConverter(id_base);
		eqn_scalar.SetDomain_Field(id_base,world);
		eqn_scalar.SetStationary(true);
		eqn_scalar.SetTimeIntegrationParameter(dt);
		// 方程式の設定
		eqn_scalar.SetAlpha(1.0);
		eqn_scalar.SetCapacity(30.0);
		eqn_scalar.SetAdvection(false);
		eqn_scalar.SetSaveStiffMat(false);
		{
			Fem::Eqn::CEqn_Scalar2D eqn1 = eqn_scalar.GetEquation(conv.GetIdEA_fromCad(2,Cad::LOOP));
			eqn1.SetAlpha(10.0);
			eqn_scalar.SetEquation(eqn1);
		}

		id_val_bc0 = eqn_scalar.AddFixElemAry(conv.GetIdEA_fromCad(3,Cad::EDGE),world);
    field_value_setter_ary.clear();
		{
      Fem::Field::CFieldValueSetter fvs(id_val_bc0,world);
      fvs.SetMathExp("cos(2*PI*t+0.1)",0,Fem::Field::VALUE,world);
      field_value_setter_ary.push_back(fvs);            
		}
		id_val_bc1 = eqn_scalar.AddFixElemAry(conv.GetIdEA_fromCad(6,Cad::EDGE),world);
    Fem::Field::SetFieldValue_Constant(id_val_bc1,0,Fem::Field::VALUE,world,+1.0);    

		// 描画オブジェクトの登録
		drawer_ary.Clear();
		const unsigned int id_field_val = eqn_scalar.GetIdField_Value();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_val,true,world,id_field_val,-1.0,1.0) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_val,true,world) );
		drawer_ary.InitTrans(camera);
	}
	else if( iprob == 16 ){
		eqn_scalar.SetSaveStiffMat(true);
	}
	else if( iprob == 17 ){
		eqn_scalar.SetSaveStiffMat(false);
		eqn_scalar.SetStationary(false);
	}
	else if( iprob == 18 ){
		eqn_scalar.SetSaveStiffMat(true);
	}
	else if( iprob == 19 ){
    Cad::CCadObj2D cad_2d;
		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(0.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.0) );
			vec_ary.push_back( Com::CVector2D(1.0,0.1) );
			vec_ary.push_back( Com::CVector2D(0.0,0.1) );
			cad_2d.AddPolygon(vec_ary);
		}
		world.Clear();
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.05) );
		const CIDConvEAMshCad& conv = world.GetIDConverter(id_base);

		eqn_scalar.SetDomain_Field(id_base,world);
		eqn_scalar.SetStationary(false);
		eqn_scalar.SetSaveStiffMat(false);
		dt = 1;
		eqn_scalar.SetTimeIntegrationParameter(dt,0.5);
		eqn_scalar.SetAxialSymmetry(true);
		// 全体の方程式の係数設定
		eqn_scalar.SetAlpha(48.0);
		eqn_scalar.SetCapacity(480*7.86*1000);
		eqn_scalar.SetAdvection(0);
		eqn_scalar.SetSource(0);

		const unsigned int id_field_val = eqn_scalar.GetIdField_Value();
    Fem::Field::SetFieldValue_Constant(id_field_val,0,Fem::Field::VALUE,world,500.0);    

		id_val_bc1 = eqn_scalar.AddFixElemAry(conv.GetIdEA_fromCad(2,Cad::EDGE),world);
    Fem::Field::SetFieldValue_Constant(id_val_bc1,0,Fem::Field::VALUE,world,+0.0);    

		// 描画オブジェクトの登録
		drawer_ary.Clear();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_val,true,world,id_field_val,0,500) );
		drawer_ary.PushBack( new View::CDrawerEdge(id_field_val,true,world) );
		drawer_ary.InitTrans(camera);	// 視線座標変換行列の初期化
	}

	iprob++;
	if( iprob == nprob ) iprob=0;

	return true;
}


////////////////////////////////////////////////////////////////

double cur_time = 0.0;
bool is_animation = true;

// リサイズ時のコールバック関数
void myGlutResize(int w, int h)
{
	camera.SetWindowAspect((double)w/h);
	::glViewport(0, 0, w, h);

	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
	Com::View::SetProjectionTransform(camera);
	::glutPostRedisplay();
}

// 描画時のコールバック関数
void myGlutDisplay(void)
{
//	::glClearColor(0.2, 0.7, 0.7 ,1.0);
	::glClearColor(1.0, 1.0, 1.0 ,1.0);
	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	::glEnable(GL_DEPTH_TEST);

	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 1.1f, 4.0f );

	::glMatrixMode(GL_MODELVIEW);
	::glLoadIdentity();
	Com::View::SetModelViewTransform(camera);

	if( is_animation ){
		cur_time += dt;
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

	ShowFPS();
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
	// glutの初期設定
	glutInitWindowPosition(200,200);
	glutInitWindowSize(400, 300);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
	glutCreateWindow("DelFEM demo");

	// コールバック関数の設定
	glutDisplayFunc(myGlutDisplay);
	glutReshapeFunc(myGlutResize);
	glutMotionFunc(myGlutMotion);
	glutMouseFunc(myGlutMouse);
	glutKeyboardFunc(myGlutKeyboard);
	glutSpecialFunc(myGlutSpecial);
	glutIdleFunc(myGlutIdle);

    // 問題設定
	SetNewProblem();

	// メインループ
	glutMainLoop();
	return 0;
}
