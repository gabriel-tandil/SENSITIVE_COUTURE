
#include <QtGui>
#include <QtOpenGL>

#include <math.h>

#include "glwidget_solid2d.h"

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

GLWidget_Solid2d::GLWidget_Solid2d(QWidget *parent)
    : QGLWidget(parent)
{
    // View Default Setting
    this->is_show_coord = false;
    this->is_draw_edge_deformed = true;
    this->is_draw_edge_undeformed = true;
    this->itype_face_draw = 1;

    cur_time = 0;
    dt = 0.05;
    SetNewProblem();   

    QTimer *timer = new QTimer(this);
    connect(timer, SIGNAL(timeout()), this, SLOT(StepTime()));
    timer->start(20);      
}

GLWidget_Solid2d::~GLWidget_Solid2d()
{
    makeCurrent();
}

void GLWidget_Solid2d::initializeGL()
{

}

void GLWidget_Solid2d::paintGL()
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

void GLWidget_Solid2d::resizeGL(int w, int h)
{
    camera.SetWindowAspect((double)w/h);
    ::glViewport(0, 0, w, h);
    ::glMatrixMode(GL_PROJECTION);
    ::glLoadIdentity();
    Com::View::SetProjectionTransform(camera);
    updateGL();
}

void GLWidget_Solid2d::mousePressEvent(QMouseEvent *event)
{
}

void GLWidget_Solid2d::mouseMoveEvent(QMouseEvent *event)
{
}

void GLWidget_Solid2d::StepTime()
{
    cur_time += dt;
    for(unsigned int ifvs=0;ifvs<field_value_setter_ary.size();ifvs++){
      field_value_setter_ary[ifvs].ExecuteValue(cur_time,world);
    }
    //world.FieldValueExec(cur_time);
    solid.Solve(world);
    if( id_field_equiv_stress != 0 ){ solid.SetEquivStressValue(id_field_equiv_stress,world); } // Åe?Åg?ÅÒ??I?I?e?d?X?V
    if( id_field_stress       != 0 ){ solid.SetStressValue(     id_field_stress,      world); } // ÅÒ??I?I?e?d?X?V
//    world.FieldValueDependExec();
    drawer_ary.Update(world);
    updateGL();
}

void GLWidget_Solid2d::SetDrawerArray(){
    this->drawer_ary.Clear();
    if( this->is_show_coord ){
    }
    const unsigned int id_field_disp = solid.GetIdField_Disp();
//    drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_disp,false,world) );
    if(      this->itype_face_draw==1 ){    // plane deform
        drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_disp,false,world) );
    }
    else if( this->itype_face_draw==2 ){    // equiv stress deform
       drawer_ary.PushBack( new View::CDrawerFace(id_field_disp,false,world,id_field_equiv_stress, 0,0.5) );
    }
    else if( this->itype_face_draw==3 ){    // equiv stress undeform
        drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_equiv_stress,false,world,id_field_equiv_stress, 0,0.5) );
    }
    else if( this->itype_face_draw==4 ){    // temp deform
        assert( world.IsIdField(id_field_temp) );
        drawer_ary.PushBack( new View::CDrawerFace(id_field_disp,false, world, id_field_temp, -1,1) );
    }
    else if( this->itype_face_draw==5 ){    // temp undeform
        assert( world.IsIdField(id_field_temp) );
        drawer_ary.PushBack( new View::CDrawerFace(id_field_temp,true, world, id_field_temp, -1,1) );
    }
    ////////
    if( this->is_draw_edge_deformed ){
        drawer_ary.PushBack( new Fem::Field::View::CDrawerEdge(id_field_disp,false,world) );
    }
    if( this->is_draw_edge_undeformed ){
        drawer_ary.PushBack( new Fem::Field::View::CDrawerEdge(id_field_disp,true ,world) );
    }
}



void GLWidget_Solid2d::SetNewProblem()
{
    const unsigned int nprob = 13;	// number of problem settings
    static unsigned int iprob = 0;

    static unsigned int id_field_disp_fix0 = 0;

    if( iprob == 0 )
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
        const Fem::Field::CIDConvEAMshCad conv = world.GetIDConverter(id_base);

        cur_time = 0;
        dt = 0.05;
        solid.Clear();
        solid.UpdateDomain_Field(id_base, world);
        solid.SetSaveStiffMat(false);
        solid.SetStationary(true);
        // Set Material Property
        solid.SetYoungPoisson(10.0,0.3,true);	// ???Åg?O?|?A?|?A?\?ÅgÅha?I?YÅfe(???EÅÒ??I)
        solid.SetGeometricalNonlinear(false);
        solid.SetGravitation(0.0,0.0);
        solid.SetTimeIntegrationParameter(dt,0.7);

        unsigned int id_field_bc0 = solid.AddFixElemAry(conv.GetIdEA_fromCad(2,Cad::EDGE),world);
        unsigned int id_field_bc1 = solid.AddFixElemAry(conv.GetIdEA_fromCad(4,Cad::EDGE),world);
        Fem::Field::CFieldValueSetter fvs(id_field_bc1,world);
        fvs.SetMathExp("sin(t*PI*2*0.1)", 1,Fem::Field::VALUE,world);
        field_value_setter_ary.clear();
        field_value_setter_ary.push_back(fvs);

        // View Setting
        this->itype_face_draw = 1;  // plain deform
        this->SetDrawerArray();
        drawer_ary.InitTrans(camera);	// View Fitting
    }
    else if( iprob == 1 )	// Disp Stress Value
    {
        const unsigned int id_field_disp = this->solid.GetIdField_Disp();
        id_field_equiv_stress = world.MakeField_FieldElemDim(id_field_disp,2,SCALAR,VALUE,BUBBLE);
        solid.SetGeometricalNonlinear(false);
        solid.SetStationary(true);

        // View Setting
        this->itype_face_draw = 3;  // stress undeform
        this->SetDrawerArray();
        drawer_ary.InitTrans(camera);
    }
    else if( iprob == 2 )	// Disp Stress Value
    {
        // View Setting
        this->itype_face_draw = 2;  // stress deform
        this->SetDrawerArray();
        drawer_ary.InitTrans(camera);
    }
    else if( iprob == 3 )	// ÅhMÅÒ??I?aÅee
    {
        id_field_equiv_stress = 0;
        id_field_stress = 0;
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

        cur_time = 0;
        dt = 0.05;
        solid.UpdateDomain_Field(id_base,world);
        solid.SetSaveStiffMat(false);	// ???Å·?s?n?UÅeÅ˜?É ?E?Åë
        solid.SetStationary(true);	// ?AÅgI?aÅee?E?Z?b?g
        // Set Material Property
        solid.SetYoungPoisson(10.0,0.3,true);	// ???Åg?O?|?A?|?A?\?ÅgÅha?I?YÅfe(???EÅÒ??I)
        solid.SetGeometricalNonlinear(false);	// ?oÅÒ??wÅgIÅhn?u?`?Å·?3??
        solid.SetGravitation(0.0,-0.1);	// ?d?I?I?YÅfe
        solid.SetTimeIntegrationParameter(dt);	// ???O?X?e?b?v?YÅfe

        unsigned int id_field_bc0 = solid.AddFixElemAry(conv.GetIdEA_fromCad(2,Cad::EDGE),world);
        unsigned int id_field_bc1 = solid.AddFixElemAry(conv.GetIdEA_fromCad(6,Cad::EDGE),world);

        // Set Temparatuer Distribution
        const unsigned int id_field_disp = solid.GetIdField_Disp();
        id_field_temp = world.MakeField_FieldElemDim(id_field_disp,2,SCALAR,VALUE,CORNER);

        Fem::Field::CFieldValueSetter fvs(id_field_temp,world);
        fvs.SetMathExp("sin(6.28*y)*sin(x)*sin(t)", 0,Fem::Field::VALUE, world);
        field_value_setter_ary.clear();
        field_value_setter_ary.push_back(fvs);

        solid.SetThermalStress(id_field_temp);
        solid.ClearFixElemAry(3,world);

        this->itype_face_draw = 4;  // temp deform
        this->SetDrawerArray();
        drawer_ary.InitTrans(camera);	// View Fitting
    }
    else if( iprob == 4 )	// ÅhMÅÒ??I?d?l?Å˜?ÅE?e?Å}?A?d?a?s?e
    {
        this->itype_face_draw = 5;  // temp undeform
        this->SetDrawerArray();
        drawer_ary.InitTrans(camera);	// View Fitting
    }
    else if( iprob == 5 )	// ÅhMÅÒ??I?d?l?Å˜?ÅE?e?Å}?A?d?a?s?e
    {
        solid.SetThermalStress(0);
    }
    else if( iprob == 6 )
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
        CIDConvEAMshCad conv = world.GetIDConverter(id_base);		// ID?I?ÅE?N?ÅÒ?X

        cur_time = 0;
        dt = 0.05;
        solid.SetDomain_FieldEA(id_base,conv.GetIdEA_fromCad(1,Cad::LOOP),world);
        solid.SetSaveStiffMat(true);
        solid.SetStationary(true);
        solid.SetTimeIntegrationParameter(dt);	// ?^?C???X?e?b?v?d?YÅfe
        solid.SetYoungPoisson(2.5,0.3,true);	// ???Åg?O?|?A?|?A?\?ÅgÅha?I?YÅfe(???EÅÒ??I)
        solid.SetGeometricalNonlinear(false);	// ?oÅÒ??wÅgIÅhn?u?`?Å·?d?l?Å˜?É ?E?Åë
        solid.SetGravitation(0.0,0.0);	// ?d?I?O

        unsigned int id_field_bc1 = solid.AddFixElemAry(conv.GetIdEA_fromCad(3,Cad::EDGE),world);
        unsigned int id_field_bc2 = solid.AddFixElemAry(conv.GetIdEA_fromCad(1,Cad::EDGE),world);

        Fem::Field::CFieldValueSetter fvs(id_field_bc1,world);
        fvs.SetMathExp("0.3*sin(1.5*t)", 0,Fem::Field::VALUE, world);
        fvs.SetMathExp("0.1*(cos(t)+1)", 1,Fem::Field::VALUE, world);
        field_value_setter_ary.clear();
        field_value_setter_ary.push_back(fvs);

        // View Setting
        this->itype_face_draw = 1;  // plain deform
        this->SetDrawerArray();
        const unsigned int id_field_disp = solid.GetIdField_Disp();
        drawer_ary.InitTrans(camera);	// View Fitting
    }
    else if( iprob == 7 ){
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
        CIDConvEAMshCad conv = world.GetIDConverter(id_base);  // ID?I?ÅE?N?ÅÒ?X

        cur_time = 0;
        dt = 0.05;
        solid.SetDomain_FieldEA(id_base,conv.GetIdEA_fromCad(2,Cad::LOOP),world);
        solid.SetTimeIntegrationParameter(dt);
        solid.SetSaveStiffMat(false);
        solid.SetStationary(true);

        solid.SetYoungPoisson(3.0,0.3,true);	// ???Åg?O?|?A?|?A?\?ÅgÅha?I?YÅfe(???EÅÒ??I)
        solid.SetGeometricalNonlinear(false);	// ?oÅÒ??wÅgIÅhn?u?`?Å·?d?l?Å˜?É ?E?Åë
        solid.SetGravitation(0.0,0.0);

        unsigned int id_field_bc1 = solid.AddFixElemAry(conv.GetIdEA_fromCad(3,Cad::EDGE),world);
        unsigned int id_field_bc2 = solid.AddFixElemAry(conv.GetIdEA_fromCad(5,Cad::EDGE),world);

        Fem::Field::CFieldValueSetter fvs(id_field_bc1,world);
        fvs.SetMathExp("0.3*sin(1.5*t)",     0,Fem::Field::VALUE, world);
        fvs.SetMathExp("0.1*(cos(t)+1)+0.1", 1,Fem::Field::VALUE, world);
        field_value_setter_ary.clear();
        field_value_setter_ary.push_back(fvs);

        // View Setting
        this->itype_face_draw = 1;  // plain deform
        this->SetDrawerArray();
        const unsigned int id_field_disp = solid.GetIdField_Disp();
        drawer_ary.InitTrans(camera);	// View Fitting
    }
    else if( iprob == 8 ){
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
        const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);  // ID?I?ÅE?N?ÅÒ?X

        cur_time = 0;
        dt = 0.05;
        solid.UpdateDomain_Field(id_base,world);	// ÅÒd?I?I?a?d?Z?b?g
        solid.SetTimeIntegrationParameter(dt);	// ???O???Y?d?Z?b?g
        solid.SetSaveStiffMat(false);	// ???Å·?s?n?d?UÅeÅ˜?É ?E?Åë
        solid.SetStationary(true);	// ?AÅgI?aÅee
        // ÅeSÅeI?I?ÅN?Å·Åfn?d?Z?b?g
        solid.SetYoungPoisson(1.0,0.3,true);	// ???Åg?O?|?A?|?A?\?ÅgÅha?I?YÅfe(???EÅÒ??I)
        solid.SetGeometricalNonlinear(false);
        solid.SetGravitation(0.0,-0.0);

        {	// St.Venant-KirchhoffÅeI
            Fem::Eqn::CEqn_Solid2D eqn = solid.GetEquation(conv.GetIdEA_fromCad(1,Cad::LOOP));
            eqn.SetGeometricalNonlinear(true);
            solid.SetEquation(eqn);
        }
        {	// ?_?c?c?ÅëÅfe?Å·ÅeI
            Fem::Eqn::CEqn_Solid2D eqn = solid.GetEquation(conv.GetIdEA_fromCad(2,Cad::LOOP));
            eqn.SetYoungPoisson(0.1,0.3,true);
            solid.SetEquation(eqn);
        }
        unsigned int id_field_temp = world.MakeField_FieldElemAry(id_base, conv.GetIdEA_fromCad(3,Cad::LOOP),SCALAR,VALUE,CORNER);
        Fem::Field::CFieldValueSetter fvs(id_field_temp,world);
        fvs.SetMathExp("0.1*sin(3.14*4*y)*sin(2*t)", 0,Fem::Field::VALUE, world);
        field_value_setter_ary.clear();
        field_value_setter_ary.push_back(fvs);

        {	// ÅhMÅÒ??I?d?l?Å˜?É ???u?`Åfe?Å·ÅeI
            Fem::Eqn::CEqn_Solid2D eqn = solid.GetEquation(conv.GetIdEA_fromCad(3,Cad::LOOP));
            eqn.SetThermalStress(id_field_temp);
            solid.SetEquation(eqn);
        }
        {	// ?d?Åë?u?`Åfe?Å·ÅeI
            Fem::Eqn::CEqn_Solid2D eqn = solid.GetEquation(conv.GetIdEA_fromCad(4,Cad::LOOP));
            eqn.SetYoungPoisson(10,0.3,true);
            solid.SetEquation(eqn);
        }

        id_field_disp_fix0 = solid.AddFixElemAry(conv.GetIdEA_fromCad(2,Cad::EDGE),world);

        // View Setting
        this->itype_face_draw = 1;  // plain deform
        this->SetDrawerArray();
        drawer_ary.InitTrans(camera);	// View Fitting
    }
    else if( iprob == 9 ){
        solid.SetRho(0.0001);
        solid.SetStationary(false);
        Fem::Field::CFieldValueSetter fvs(id_field_disp_fix0,world);
        fvs.SetMathExp("0.5*cos(2*t)", 1,Fem::Field::VALUE, world);
//        field_value_setter_ary.clear();
        field_value_setter_ary.push_back(fvs);
    }
    else if( iprob == 10 ){
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
        const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);  // ID?I?ÅE?N?ÅÒ?X

        cur_time = 0;
        dt = 0.05;
        solid.UpdateDomain_Field(id_base,world);	// ÅÒd?I?I?a?d?Z?b?g
        solid.SetTimeIntegrationParameter(dt);	// ???O???Y?d?Z?b?g
        solid.SetSaveStiffMat(false);	// ???Å·?s?n?d?UÅeÅ˜?É ?E?Åë
        solid.SetStationary(false);	// ?AÅgI?aÅee
        // ÅeSÅeI?I?ÅN?Å·Åfn?d?Z?b?g
        solid.SetYoungPoisson(1.0,0.3,true);	// ???Åg?O?|?A?|?A?\?ÅgÅha?I?YÅfe(???EÅÒ??I)
        solid.SetGeometricalNonlinear(false);
        solid.SetGravitation(0.0,-0.0);
        solid.SetRho(0.001);

        {	// ?_?c?c?ÅëÅfe?Å·ÅeI
            Fem::Eqn::CEqn_Solid2D eqn = solid.GetEquation(conv.GetIdEA_fromCad(1,Cad::LOOP));
            eqn.SetYoungPoisson(0.1,0.3,true);
            solid.SetEquation(eqn);
        }
        {	// ?d?Åë?u?`Åfe?Å·ÅeI
            Fem::Eqn::CEqn_Solid2D eqn = solid.GetEquation(conv.GetIdEA_fromCad(2,Cad::LOOP));
            eqn.SetYoungPoisson(100000000,0.3,true);
            solid.SetEquation(eqn);
        }

//		id_field_disp_fix0 = solid.AddFixElemAry(conv.GetIdEA_fromCad(2,1),world);
        const unsigned int id_field_bc1 = solid.AddFixElemAry(conv.GetIdEA_fromCad(4,Cad::EDGE),world);

        Fem::Field::CFieldValueSetter fvs(id_field_bc1,world);
        fvs.SetMathExp("0.3*sin(1.5*t)",     0,Fem::Field::VALUE, world);
        fvs.SetMathExp("0.1*(cos(t)+1)+0.1", 1,Fem::Field::VALUE, world);
        field_value_setter_ary.clear();
        field_value_setter_ary.push_back(fvs);

        // Set Visualization
        this->itype_face_draw = 1;  // plain deform
        this->SetDrawerArray();
        const unsigned int id_field_disp = solid.GetIdField_Disp();
        drawer_ary.InitTrans(camera);	// View Fitting
    }
    else if( iprob == 11 )
    {
      Cad::CCadObj2D cad_2d;      
      unsigned int id_l;
        unsigned int id_e1,id_e2,id_e3,id_e4,id_e5;
        {   // define shape
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
        const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);  // ID?I?ÅE?N?ÅÒ?X
        unsigned int id_base2 = 0;
        {
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
        cur_time = 0;
        dt = 0.05;
        solid.UpdateDomain_Field(id_base2, world);
        solid.SetSaveStiffMat(false);
        solid.SetStationary(true);
        // Set Material Property
        solid.SetYoungPoisson(10.0,0.3,true);	// ???Åg?O?|?A?|?A?\?ÅgÅha?I?YÅfe(???EÅÒ??I)
//		  solid.SetRho(10);
        solid.SetGeometricalNonlinear(false);
        solid.SetGravitation(0.0,0.0);
        solid.SetTimeIntegrationParameter(dt,0.7);

        unsigned int id_field_bc0 = solid.AddFixElemAry(conv.GetIdEA_fromCad(3,Cad::EDGE),world);
        unsigned int id_field_bc1 = solid.AddFixElemAry(conv.GetIdEA_fromCad(5,Cad::EDGE),world);
        Fem::Field::CFieldValueSetter fvs(id_field_bc1,world);
        fvs.SetMathExp("0.1*sin(t*PI*2*0.1)",     0,Fem::Field::VALUE, world);
        fvs.SetMathExp("0.1*(1-cos(t*PI*2*0.1))", 1,Fem::Field::VALUE, world);
        field_value_setter_ary.clear();
        field_value_setter_ary.push_back(fvs);

        // Set Visualization
        this->itype_face_draw = 1;  // plane deform
        this->SetDrawerArray();
        drawer_ary.InitTrans(camera);	// View Fitting
    }
    else if( iprob == 12 )
    {
        cur_time = 0;
        dt = 0.05;
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
            unsigned int id_v1 = cad_2d.AddVertex(Cad::EDGE,3,Com::CVector2D(2.5,2.0)).id_v_add;
            unsigned int id_v2 = cad_2d.AddVertex(Cad::LOOP,id_l,Com::CVector2D(2.5,1.0)).id_v_add;
            id_e = cad_2d.ConnectVertex_Line(id_v1,id_v2).id_e_add;
        }
        Msh::CMesher2D mesh_2d(cad_2d,0.2);
        world.Clear();
        cur_time = 0;
        dt = 0.05;
        const unsigned int id_base = world.AddMesh(mesh_2d);
        const Fem::Field::CIDConvEAMshCad conv = world.GetIDConverter(id_base);
        unsigned int id_base2 = 0;
        {
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
        // Setting Material Parameter
        solid.SetYoungPoisson(10.0,0.3,true);	// ???Åg?O?|?A?|?A?\?ÅgÅha?I?YÅfe(???EÅÒ??I)
        solid.SetGeometricalNonlinear(false);
        solid.SetGravitation(0.0,0.0);
        solid.SetTimeIntegrationParameter(dt,0.7);

        unsigned int id_field_bc0 = solid.AddFixElemAry(conv.GetIdEA_fromCad(2,Cad::EDGE),world);
        unsigned int id_field_bc1 = solid.AddFixElemAry(conv.GetIdEA_fromCad(4,Cad::EDGE),world);
        Fem::Field::CFieldValueSetter fvs(id_field_bc1,world);
        fvs.SetMathExp("0.5*(1-cos(t*PI*2*0.1))",0,Fem::Field::VALUE, world);
        fvs.SetMathExp("0.2*sin(t*PI*2*0.1)",    1,Fem::Field::VALUE, world);
        field_value_setter_ary.clear();
        field_value_setter_ary.push_back(fvs);

        const unsigned int id_field_disp = solid.GetIdField_Disp();
        id_field_equiv_stress = world.MakeField_FieldElemDim(id_field_disp,2,SCALAR,VALUE,BUBBLE);

        // Setting of Visualizatoin
        this->itype_face_draw = 2;  // stress defrom
        this->SetDrawerArray();
        drawer_ary.InitTrans(camera);   // View Fitting
    }

    iprob++;
    if( iprob == nprob ){
        iprob = 0;
    }
}

