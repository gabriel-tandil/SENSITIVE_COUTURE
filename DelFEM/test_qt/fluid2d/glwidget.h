#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QGLWidget>

#include "delfem/camera.h"
#include "delfem/field.h"
#include "delfem/field_world.h"
#include "delfem/field_value_setter.h"
#include "delfem/drawer_field.h"
#include "delfem/eqnsys_fluid.h"


class GLWidget : public QGLWidget
{
    Q_OBJECT

public:
    GLWidget(QWidget *parent = 0);
    ~GLWidget();

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
    Fem::Eqn::CEqnSystem_Fluid2D fluid;
    Fem::Field::CFieldWorld world;
    Fem::Field::CFieldValueSetter field_value_setter;
private:
    Com::View::CCamera camera;

    Fem::Field::View::CDrawerArrayField drawer_ary;
    double cur_time;
    double dt;
    unsigned int id_base;
};

#endif
