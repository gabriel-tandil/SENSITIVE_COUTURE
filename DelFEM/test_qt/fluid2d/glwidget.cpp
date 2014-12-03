
#include <QtGui>
#include <QtOpenGL>

#include <math.h>

#include "glwidget.h"

#include "delfem/camera.h"
#include "delfem/cad_obj2d.h"
#include "delfem/mesh3d.h"
#include "delfem/field.h"
#include "delfem/field_world.h"
#include "delfem/drawer_field.h"	// fem visualization class
#include "delfem/drawer_field_face.h"
#include "delfem/drawer_field_edge.h"
#include "delfem/drawer_field_vector.h"
#include "delfem/drawer_field_image_based_flow_vis.h"
#include "delfem/drawer_field_streamline.h"

#include "delfem/eqnsys_fluid.h"		// class CEqnSystem_Fluid2D


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
//    world.FieldValueExec(cur_time);
  fluid.Solve(world);
//    world.FieldValueDependExec();
  drawer_ary.Update(world);
  updateGL();
}

void GLWidget::SetNewProblem()
{
  const unsigned int nprob = 11;	// number of problems
	static unsigned int iprob = 0;
	
  if( iprob == 0 )	// cavity flow， static storks fluid
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
		
		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerVector(id_field_velo,world) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_press,true,world, id_field_press) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerEdge(id_field_velo,true,world) );
		//		drawer_ary.PushBack( new Fem::Field::View::CDrawerImageBasedFlowVis(id_field_velo,world,1) );
		//		drawer_ary.PushBack( new Fem::Field::View::CDrawerStreamline(id_field_velo,world) );
    drawer_ary.InitTrans( camera );
	}
  else if( iprob == 1 )	// cavity flow，nonstatic storks
	{
		fluid.SetIsStationary(false);
	}
  else if( iprob == 2 )	// cavity flow，nonstatic storks(big Rho)
	{
		fluid.SetRho(0.5);
	}
  else if( iprob == 3 )	// cavity flow, Naiver-Stokes flow
	{
		fluid.SetRho(0.02);
		fluid.SetMyu(0.00001);
		fluid.SetNavierStokes();
	}
  else if( iprob == 4 )	// cavity flow，bubble interpolation，nonstatic navier-storks flow
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
    field_value_setter.SetMathExp("0.5*sin(0.1*t)", 0,Fem::Field::VELOCITY, world);
/*		{
			Fem::Field::CField& bc0_field_velo = world.GetField(id_field_bc0);
			bc0_field_velo.SetValue("0.5*sin(0.1*t)", 0,Fem::Field::VELOCITY, world,true);
			//			bc0_field_velo.SetVelocity("0.1", 0, world,true);
    }*/
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
  else if( iprob == 5 )	// L shaped channel
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
/*		{
			Fem::Field::CField& bc0_field_velo = world.GetField(id_field_bc0);
			bc0_field_velo.SetValue("0.1*sin(0.1*t)", 0,Fem::Field::VELOCITY, world,true);
			//			bc0_field_velo.SetVelocity("0.1", 0, world,true);
    }*/
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
		{	// 正方形にが２つに分割
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
/*		{
			Fem::Field::CField& field = world.GetField(id_field_bc1);
			field.SetValue("0.3*sin(0.5*t)", 1,Fem::Field::VELOCITY, world,true);
    }*/
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
  else if( iprob == 10 ){	// karman vortex sheets
		Cad::CCadObj2D cad_2d;
    {	// quadrateral hole in the channel
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
/*		{
			Fem::Field::CField& field = world.GetField(id_field_bc1);
			field.SetValue(0.1,0,Fem::Field::VELOCITY,world,true);
    }      */
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
	
	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
  Com::View::SetProjectionTransform(camera);
	//	::glutPostRedisplay();
	
	iprob++;
	if( iprob == nprob ){ iprob = 0; }
}
