#ifndef DIALOG_FLUID2D_H
#define DIALOG_FLUID2D_H

#include <QDialog>

#include "delfem/field.h"
#include "delfem/field_world.h"
#include "delfem/eqnsys_fluid.h"

namespace Ui {
    class Dialog_Fluid2D;
}

class QAbstractButton;
class Dialog_Fluid2D : public QDialog {
    Q_OBJECT
public:
    Dialog_Fluid2D(const Fem::Field::CFieldWorld& world_input,
                   Fem::Eqn::CEqnSystem_Fluid2D& fluid_input,
                   QWidget *parent = 0);
    ~Dialog_Fluid2D();
private slots:
    void comboBox_currentIndexChenged(int);
    void matPorpChanged();
    void buttonBox_clicked(QAbstractButton* pbtn);

protected:
    void changeEvent(QEvent *e);

private:
    Ui::Dialog_Fluid2D *ui;

    const Fem::Field::CFieldWorld& world;
    Fem::Eqn::CEqnSystem_Fluid2D& fluid;
    std::vector<Fem::Eqn::CEqn_Fluid2D> aEqn;
};

#endif // DIALOG_FLUID2D_H
