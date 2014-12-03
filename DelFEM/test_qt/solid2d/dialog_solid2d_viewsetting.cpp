#include "dialog_solid2d_viewsetting.h"
#include "ui_dialog_solid2d_viewsetting.h"

Dialog_Solid2D_ViewSetting::Dialog_Solid2D_ViewSetting(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Dialog_Solid2D_ViewSetting)
{
    ui->setupUi(this);
}

Dialog_Solid2D_ViewSetting::~Dialog_Solid2D_ViewSetting()
{
    delete ui;
}

void Dialog_Solid2D_ViewSetting::changeEvent(QEvent *e)
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
