/********************************************************************************
** Form generated from reading UI file 'dialog_solid2d_viewsetting.ui'
**
** Created: Thu Jan 27 18:47:54 2011
**      by: Qt User Interface Compiler version 4.6.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_DIALOG_SOLID2D_VIEWSETTING_H
#define UI_DIALOG_SOLID2D_VIEWSETTING_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QDialog>
#include <QtGui/QDialogButtonBox>
#include <QtGui/QHeaderView>

QT_BEGIN_NAMESPACE

class Ui_Dialog_Solid2D_ViewSetting
{
public:
    QDialogButtonBox *buttonBox;

    void setupUi(QDialog *Dialog_Solid2D_ViewSetting)
    {
        if (Dialog_Solid2D_ViewSetting->objectName().isEmpty())
            Dialog_Solid2D_ViewSetting->setObjectName(QString::fromUtf8("Dialog_Solid2D_ViewSetting"));
        Dialog_Solid2D_ViewSetting->resize(400, 300);
        buttonBox = new QDialogButtonBox(Dialog_Solid2D_ViewSetting);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setGeometry(QRect(30, 240, 341, 32));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        retranslateUi(Dialog_Solid2D_ViewSetting);
        QObject::connect(buttonBox, SIGNAL(accepted()), Dialog_Solid2D_ViewSetting, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), Dialog_Solid2D_ViewSetting, SLOT(reject()));

        QMetaObject::connectSlotsByName(Dialog_Solid2D_ViewSetting);
    } // setupUi

    void retranslateUi(QDialog *Dialog_Solid2D_ViewSetting)
    {
        Dialog_Solid2D_ViewSetting->setWindowTitle(QApplication::translate("Dialog_Solid2D_ViewSetting", "Dialog", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class Dialog_Solid2D_ViewSetting: public Ui_Dialog_Solid2D_ViewSetting {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_DIALOG_SOLID2D_VIEWSETTING_H
