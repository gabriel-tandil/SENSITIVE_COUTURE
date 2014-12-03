
#include <QtGui>
#include <QtOpenGL>

#include <math.h>

#include "glwidget.h"

#include "delfem/camera.h"
#include "delfem/cad_obj2d.h"
#include "delfem/mesh3d.h"
#include "delfem/field.h"
#include "delfem/field_world.h"
#include "delfem/drawer_field.h"
#include "delfem/drawer_field_face.h"
#include "delfem/drawer_field_edge.h"
#include "delfem/drawer_field_vector.h"
#include "delfem/eqnsys_solid.h"

using namespace Fem::Field;

GLWidget::GLWidget(QWidget *parent)
    : QGLWidget(parent)
{
  cur_time = 0;
  dt = 0.05;
  SetNewProblem();

  QTimer *timer = new QTimer(this);
  connect(timer, SIGNAL(timeout()), this, SLOT(StepTime()));
  timer->start(20);
}

GLWidget::~GLWidget()
{
  makeCurrent();
}

void GLWidget::initializeGL()
{

}

void GLWidget::paintGL()
{
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  ::glMatrixMode(GL_MODELVIEW);
  ::glLoadIdentity();
  Com::View::SetModelViewTransform(camera);

  ::glMatrixMode(GL_PROJECTION);
  ::glLoadIdentity();
  Com::View::SetProjectionTransform(camera);

  drawer_ary.Draw();
}

void GLWidget::resizeGL(int w, int h)
{
  camera.SetWindowAspect((double)w/h);
  ::glViewport(0, 0, w, h);
  ::glMatrixMode(GL_PROJECTION);
  ::glLoadIdentity();
  Com::View::SetProjectionTransform(camera);
  updateGL();
}

void GLWidget::mousePressEvent(QMouseEvent *event)
{
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
}

void GLWidget::StepTime()
{
  cur_time += dt;
  field_value_setter.ExecuteValue(cur_time,world);
  eqn_scalar.Solve(world);
  drawer_ary.Update(world);
  updateGL();
}

void GLWidget::SetNewProblem()
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

    field_value_setter = Fem::Field::CFieldValueSetter(id_val_bc0,world);
    field_value_setter.SetMathExp("floor(1+0.8*cos(2*PI*t+0.1))",0,Fem::Field::VALUE,world);
    Fem::Field::SetFieldValue_Constant(id_val_bc1,0,Fem::Field::VALUE,world,+1.0);
    Fem::Field::SetFieldValue_Constant(id_val_bc2,0,Fem::Field::VALUE,world,-1.0);
				
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
    Fem::Field::SetFieldValue_MathExp(id_field_velo,0,Fem::Field::VELOCITY,world," (y-0.5)");
    Fem::Field::SetFieldValue_MathExp(id_field_velo,1,Fem::Field::VELOCITY,world,"-(x-0.5)");
    field_value_setter = Fem::Field::CFieldValueSetter(id_val_bc0,world);
    field_value_setter.SetMathExp("floor(1+0.8*cos(3*t))",0,Fem::Field::VALUE,world);
		{	// 周囲の固定境界条件の設定
			const CIDConvEAMshCad conv = world.GetIDConverter(id_base);
			std::vector<unsigned int> m_aIDEA;
			m_aIDEA.push_back(conv.GetIdEA_fromCad(1,Cad::EDGE));
			m_aIDEA.push_back(conv.GetIdEA_fromCad(2,Cad::EDGE));
			m_aIDEA.push_back(conv.GetIdEA_fromCad(3,Cad::EDGE));
			m_aIDEA.push_back(conv.GetIdEA_fromCad(4,Cad::EDGE));
			id_val_bc1 = eqn_scalar.AddFixElemAry(m_aIDEA,world);
      Fem::Field::SetFieldValue_Constant(id_val_bc1,0,Fem::Field::VALUE,world,-1.0);
		}
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
    field_value_setter = Fem::Field::CFieldValueSetter(id_val_bc0,world);
    field_value_setter.SetMathExp("cos(2*PI*t+0.1)",0,Fem::Field::VALUE,world);

//    Fem::Field::SetFieldValue_MathExp(id_val_bc0,0,Fem::Field::VALUE,world,);
/*		{
			CField& field = world.GetField(id_val_bc0);
			field.SetValue("cos(2*PI*t+0.1)",0,Fem::Field::VALUE,world,true);
			//			field.SetValue("floor(1+0.8*cos(2*PI*t+0.1))",0,world,true);
    }*/
		id_val_bc1 = eqn_scalar.AddFixElemAry(conv.GetIdEA_fromCad(3,Cad::EDGE),world);
    Fem::Field::SetFieldValue_Constant(id_val_bc1,0,Fem::Field::VALUE,world,1.0);
/*		{
			CField& field = world.GetField(id_val_bc1);
			field.SetValue(1.0,0,Fem::Field::VALUE,world,false);
    }*/
				
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
    {	// quadrateral hole in the square
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
    id_val_bc1 = eqn_scalar.AddFixElemAry(conv.GetIdEA_fromCad(6,Cad::EDGE),world);

    field_value_setter = Fem::Field::CFieldValueSetter(id_val_bc0,world);
    field_value_setter.SetMathExp("cos(2*PI*t+0.1)",0,Fem::Field::VALUE,world);
    Fem::Field::SetFieldValue_Constant(id_val_bc1,0,Fem::Field::VALUE,world,1.0);
			
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


    field_value_setter = Fem::Field::CFieldValueSetter();
    Fem::Field::SetFieldValue_Constant(id_val_bc1,0,Fem::Field::VALUE,world,0.0);
		
		drawer_ary.Clear();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_val,true,world,id_field_val,0,500) );
    drawer_ary.PushBack( new View::CDrawerEdge(id_field_val,true,world) );
    drawer_ary.InitTrans(camera);	// initialize camera
	}
	
	iprob++;
	if( iprob == nprob ) iprob=0;

  return;
}
