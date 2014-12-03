#include "dialog_scalar2d.h"
#include "ui_dialog_scalar2d.h"

Dialog_Scalar2D::Dialog_Scalar2D(const Fem::Field::CFieldWorld& world_input,
                               Fem::Eqn::CEqnSystem_Scalar2D& scalar_input,
                               QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Dialog_Scalar2D),
    world(world_input),
    scalar(scalar_input)
{
    ui->setupUi(this);

    const unsigned int id_field_val = scalar.GetIdField_Value();
    const Fem::Field::CField& field = world.GetField(id_field_val);
    const std::vector<unsigned int>& aIdEA = field.GetAry_IdElemAry();
    for(unsigned int iiea=0;iiea<aIdEA.size();iiea++)
    {
        unsigned int id_ea = aIdEA[iiea];
        char str_id[16];
        sprintf(str_id,"%d",id_ea);
        ui->comboBox_IdEA->addItem(str_id);
        aEqn.push_back(scalar.GetEquation(id_ea));
    }

   connect(ui->comboBox_IdEA,           SIGNAL(currentIndexChanged(int)),   this,SLOT(comboBox_currentIndexChenged(int)));
   connect(ui->doubleSpinBox_alpha,     SIGNAL(valueChanged(double)),       this,SLOT(matPorpChanged()));
   connect(ui->doubleSpinBox_capacity,  SIGNAL(valueChanged(double)),       this,SLOT(matPorpChanged()));
   connect(ui->doubleSpinBox_source,    SIGNAL(valueChanged(double)),       this,SLOT(matPorpChanged()));

   connect(ui->buttonBox,SIGNAL(clicked(QAbstractButton*)), this,SLOT(buttonBox_clicked(QAbstractButton*)));

   this->comboBox_currentIndexChenged(0);
}

Dialog_Scalar2D::~Dialog_Scalar2D()
{
    delete ui;
}

void Dialog_Scalar2D::changeEvent(QEvent *e)
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


void Dialog_Scalar2D::comboBox_currentIndexChenged(int ind)
{    
    if( ind < 0 || ind >= this->ui->comboBox_IdEA->count() ) return;

    disconnect(ui->comboBox_IdEA,           SIGNAL(currentIndexChanged(int)),   this,SLOT(comboBox_currentIndexChenged(int)));
    disconnect(ui->doubleSpinBox_alpha,     SIGNAL(valueChanged(double)),       this,SLOT(matPorpChanged()));
    disconnect(ui->doubleSpinBox_capacity,  SIGNAL(valueChanged(double)),       this,SLOT(matPorpChanged()));
    disconnect(ui->doubleSpinBox_source,    SIGNAL(valueChanged(double)),       this,SLOT(matPorpChanged()));

    const Fem::Eqn::CEqn_Scalar2D& eqn = aEqn[ind];
    ui->doubleSpinBox_alpha->setValue(eqn.GetAlpha());
    ui->doubleSpinBox_capacity->setValue(eqn.GetCapacity());
    ui->doubleSpinBox_source->setValue(eqn.GetSource());

    connect(ui->comboBox_IdEA,           SIGNAL(currentIndexChanged(int)),   this,SLOT(comboBox_currentIndexChenged(int)));
    connect(ui->doubleSpinBox_alpha,     SIGNAL(valueChanged(double)),       this,SLOT(matPorpChanged()));
    connect(ui->doubleSpinBox_capacity,  SIGNAL(valueChanged(double)),       this,SLOT(matPorpChanged()));
    connect(ui->doubleSpinBox_source,    SIGNAL(valueChanged(double)),       this,SLOT(matPorpChanged()));
}

void Dialog_Scalar2D::matPorpChanged()
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
    Fem::Eqn::CEqn_Scalar2D& eqn = aEqn[ieqn];
    eqn.SetAlpha(ui->doubleSpinBox_alpha->value());
    eqn.SetCapacity(ui->doubleSpinBox_capacity->value());
    eqn.SetSource(ui->doubleSpinBox_source->value());
}

void Dialog_Scalar2D::buttonBox_clicked(QAbstractButton* pbtn)
{
    if( this->ui->buttonBox->buttonRole(pbtn) != this->ui->buttonBox->ApplyRole ||
        this->ui->buttonBox->buttonRole(pbtn) != this->ui->buttonBox->AcceptRole ){
        for(unsigned int ieqn=0;ieqn<aEqn.size();ieqn++){
            scalar.SetEquation(aEqn[ieqn]);
        }
    }
}
