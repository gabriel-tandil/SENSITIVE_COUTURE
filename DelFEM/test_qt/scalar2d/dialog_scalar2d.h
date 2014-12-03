#ifndef DIALOG_SCALAR2D_H
#define DIALOG_SCALAR2D_H

#include <QDialog>

#include "delfem/field.h"
#include "delfem/field_world.h"
#include "delfem/eqnsys_scalar.h"

class QAbstractButton;
namespace Ui {
    class Dialog_Scalar2D;
}

class Dialog_Scalar2D : public QDialog {
    Q_OBJECT
public:
    Dialog_Scalar2D(const Fem::Field::CFieldWorld& world_input,
                   Fem::Eqn::CEqnSystem_Scalar2D& scalar_input,
                   QWidget *parent = 0);
    ~Dialog_Scalar2D();
private slots:
    void comboBox_currentIndexChenged(int);
    void matPorpChanged();
    void buttonBox_clicked(QAbstractButton* pbtn);
protected:
    void changeEvent(QEvent *e);

private:
    Ui::Dialog_Scalar2D *ui;

    Fem::Eqn::CEqnSystem_Scalar2D& scalar;
    const Fem::Field::CFieldWorld& world;
    std::vector<Fem::Eqn::CEqn_Scalar2D> aEqn;
};

#endif // DIALOG_SCALAR2D_H
