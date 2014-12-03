#ifndef DIALOG_SOLID2D_H
#define DIALOG_SOLID2D_H

#include <QDialog>

#include "delfem/field.h"
#include "delfem/field_world.h"
#include "delfem/eqnsys_solid.h"

namespace Ui {
    class Dialog_Solid2D;
}

class QAbstractButton;
class Dialog_Solid2D : public QDialog {
    Q_OBJECT
public:
    Dialog_Solid2D(const Fem::Field::CFieldWorld& world_input,
                   Fem::Eqn::CEqnSystem_Solid2D& solid_input,
                   QWidget *parent = 0);
    ~Dialog_Solid2D();
private slots:
    void comboBox_currentIndexChenged(int);
    void matPorpChanged();
    void buttonBox_clicked(QAbstractButton* pbtn);
protected:
    void changeEvent(QEvent *e);

private:
    Ui::Dialog_Solid2D *ui;

    Fem::Eqn::CEqnSystem_Solid2D& solid;
    const Fem::Field::CFieldWorld& world;
    std::vector<Fem::Eqn::CEqn_Solid2D> aEqn;
};

#endif // DIALOG_SOLID2D_H
