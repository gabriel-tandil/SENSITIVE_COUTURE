#include <QPushButton>
#include "dialog_solid2d.h"
#include "ui_dialog_solid2d.h"

Dialog_Solid2D::Dialog_Solid2D(const Fem::Field::CFieldWorld& world_input,
                               Fem::Eqn::CEqnSystem_Solid2D& solid_input,
                               QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Dialog_Solid2D),
    world(world_input),
    solid(solid_input)
{
    ui->setupUi(this);

    const unsigned int id_field_disp = solid.GetIdField_Disp();
    const Fem::Field::CField& disp = world.GetField(id_field_disp);
    const std::vector<unsigned int>& aIdEA = disp.GetAry_IdElemAry();
    for(unsigned int iiea=0;iiea<aIdEA.size();iiea++)
    {
        unsigned int id_ea = aIdEA[iiea];
        char str_id[16];
        sprintf(str_id,"%d",id_ea);
        ui->comboBox_IdEA->addItem(str_id);
        aEqn.push_back(solid.GetEquation(id_ea));
    }

    connect(ui->comboBox_IdEA,          SIGNAL(currentIndexChanged(int)),   this,SLOT(comboBox_currentIndexChenged(int)));
    connect(ui->doubleSpinBox_Young,    SIGNAL(valueChanged(double)),       this,SLOT(matPorpChanged()));
    connect(ui->doubleSpinBox_Poisson,  SIGNAL(valueChanged(double)),       this,SLOT(matPorpChanged()));
    connect(ui->doubleSpinBox_rho,      SIGNAL(valueChanged(double)),       this,SLOT(matPorpChanged()));
    connect(ui->checkBox_nonlin,        SIGNAL(clicked(bool)),              this,SLOT(matPorpChanged()));

    connect(ui->buttonBox,SIGNAL(clicked(QAbstractButton*)),    this,SLOT(buttonBox_clicked(QAbstractButton*)));

    this->ui->buttonBox->button(QDialogButtonBox::Apply)->setEnabled(false);

    this->comboBox_currentIndexChenged(0);
}

Dialog_Solid2D::~Dialog_Solid2D()
{
    delete ui;
}

void Dialog_Solid2D::changeEvent(QEvent *e)
{
    QDialog::changeEvent(e);
    switch (e->type()) {
    case QEvent::LanguageChange:
        ui->retranslateUi(this);
        break;
    default:
        break;
    }
}

void Dialog_Solid2D::comboBox_currentIndexChenged(int ind)
{
    if( ind < 0 || ind >= this->ui->comboBox_IdEA->count() ) return;

    disconnect(ui->doubleSpinBox_Young,     SIGNAL(valueChanged(double)),   this,SLOT(matPorpChanged()));
    disconnect(ui->doubleSpinBox_Poisson,   SIGNAL(valueChanged(double)),   this,SLOT(matPorpChanged()));
    disconnect(ui->doubleSpinBox_rho,       SIGNAL(valueChanged(double)),   this,SLOT(matPorpChanged()));
    disconnect(ui->checkBox_nonlin,         SIGNAL(clicked(bool)),          this,SLOT(matPorpChanged()));

    const Fem::Eqn::CEqn_Solid2D& eqn = aEqn[ind];
    double young, poisson;
    eqn.GetYoungPoisson(young,poisson);
    ui->doubleSpinBox_Young->setValue(young);
    ui->doubleSpinBox_Poisson->setValue(poisson);
    ui->doubleSpinBox_rho->setValue(eqn.GetRho());
    ui->checkBox_nonlin->setChecked(eqn.IsGeometricalNonlinear());

    connect(ui->doubleSpinBox_Young,    SIGNAL(valueChanged(double)),   this,SLOT(matPorpChanged()));
    connect(ui->doubleSpinBox_Poisson,  SIGNAL(valueChanged(double)),   this,SLOT(matPorpChanged()));
    connect(ui->doubleSpinBox_rho,      SIGNAL(valueChanged(double)),   this,SLOT(matPorpChanged()));
    connect(ui->checkBox_nonlin,        SIGNAL(clicked(bool)),          this,SLOT(matPorpChanged()));
}

void Dialog_Solid2D::matPorpChanged()
{
    unsigned int id_ea0;
    {
        unsigned int ind = ui->comboBox_IdEA->currentIndex();
        QString itemData = ui->comboBox_IdEA->itemText(ind);
        id_ea0 = itemData.toInt();
    }
    unsigned int ieqn=0;
    for(ieqn=0;ieqn<aEqn.size();ieqn++){
        if( aEqn[ieqn].GetIdEA() == id_ea0 ) break;
    }
    if( ieqn == aEqn.size() ) return;
    Fem::Eqn::CEqn_Solid2D& eqn = aEqn[ieqn];
    double young = ui->doubleSpinBox_Young->value();
    double poisson = ui->doubleSpinBox_Poisson->value();
    eqn.SetYoungPoisson(young,poisson,true);
    eqn.SetRho(ui->doubleSpinBox_rho->value());
    eqn.SetGeometricalNonlinear(ui->checkBox_nonlin->isChecked());

    this->ui->buttonBox->button(QDialogButtonBox::Apply)->setEnabled(true);
}

void Dialog_Solid2D::buttonBox_clicked(QAbstractButton* pbtn)
{
    if( this->ui->buttonBox->buttonRole(pbtn) != this->ui->buttonBox->ApplyRole ||
        this->ui->buttonBox->buttonRole(pbtn) != this->ui->buttonBox->AcceptRole ){
        for(unsigned int ieqn=0;ieqn<aEqn.size();ieqn++){
            solid.SetEquation(aEqn[ieqn]);
        }
    }    
    this->ui->buttonBox->button(QDialogButtonBox::Apply)->setEnabled(false);
}
