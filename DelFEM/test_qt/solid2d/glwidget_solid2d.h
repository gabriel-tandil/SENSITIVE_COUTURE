#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QGLWidget>

#include "delfem/camera.h"
#include "delfem/field.h"
#include "delfem/field_world.h"
#include "delfem/field_value_setter.h"
#include "delfem/drawer_field.h"
#include "delfem/eqnsys_solid.h"

class GLWidget_Solid2d : public QGLWidget
{
    Q_OBJECT

public:
    GLWidget_Solid2d(QWidget *parent = 0);
    ~GLWidget_Solid2d();

    void setShowCoordinate(bool is_show){
        if( is_show == this->is_show_coord ) return;
        this->is_show_coord = is_show;
        this->SetDrawerArray();
    }

protected:
    void initializeGL();
    void paintGL();
    void resizeGL(int width, int height);
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);

private slots:
    void StepTime();

public:
    void SetNewProblem();
public:
    Com::View::CCamera camera;

    Fem::Field::CFieldWorld world;
    std::vector<Fem::Field::CFieldValueSetter> field_value_setter_ary;
    Fem::Field::View::CDrawerArrayField drawer_ary;
    double cur_time;
    double dt;
    Fem::Eqn::CEqnSystem_Solid2D solid;
    unsigned int id_field_equiv_stress;
    unsigned int id_field_stress;
    unsigned int id_field_temp;
private:
    void SetDrawerArray();
private:
    bool is_show_coord;
    bool is_draw_edge_deformed;
    bool is_draw_edge_undeformed;
    // 0:none
    // 1:plane_deformed
    // 2:equiv_stress_deformed
    // 3:equiv_stress_undeformed
    // 4:temp_deformed
    // 5:temp_undeformed
    unsigned int itype_face_draw;

    Com::View::CDrawerCoord* pDrawerCoord;
};

#endif
