/********************************************************************************
** Form generated from reading UI file 'dialog_fluid2d.ui'
**
** Created: Thu Jan 27 18:47:54 2011
**      by: Qt User Interface Compiler version 4.6.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_DIALOG_FLUID2D_H
#define UI_DIALOG_FLUID2D_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QComboBox>
#include <QtGui/QDialog>
#include <QtGui/QDialogButtonBox>
#include <QtGui/QDoubleSpinBox>
#include <QtGui/QGridLayout>
#include <QtGui/QGroupBox>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_Dialog_Fluid2D
{
public:
    QDialogButtonBox *buttonBox;
    QComboBox *comboBox_IdEA;
    QLabel *label_4;
    QGroupBox *MatProp;
    QWidget *gridLayoutWidget;
    QGridLayout *gridLayout;
    QDoubleSpinBox *doubleSpinBox_myu;
    QLabel *label;
    QLabel *label_2;
    QDoubleSpinBox *doubleSpinBox_rho;

    void setupUi(QDialog *Dialog_Fluid2D)
    {
        if (Dialog_Fluid2D->objectName().isEmpty())
            Dialog_Fluid2D->setObjectName(QString::fromUtf8("Dialog_Fluid2D"));
        Dialog_Fluid2D->resize(387, 281);
        Dialog_Fluid2D->setSizeGripEnabled(false);
        Dialog_Fluid2D->setModal(true);
        buttonBox = new QDialogButtonBox(Dialog_Fluid2D);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setGeometry(QRect(30, 240, 341, 32));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Apply|QDialogButtonBox::Cancel|QDialogButtonBox::Ok);
        buttonBox->setCenterButtons(false);
        comboBox_IdEA = new QComboBox(Dialog_Fluid2D);
        comboBox_IdEA->setObjectName(QString::fromUtf8("comboBox_IdEA"));
        comboBox_IdEA->setGeometry(QRect(120, 20, 111, 26));
        label_4 = new QLabel(Dialog_Fluid2D);
        label_4->setObjectName(QString::fromUtf8("label_4"));
        label_4->setGeometry(QRect(40, 20, 81, 31));
        MatProp = new QGroupBox(Dialog_Fluid2D);
        MatProp->setObjectName(QString::fromUtf8("MatProp"));
        MatProp->setGeometry(QRect(40, 60, 291, 171));
        gridLayoutWidget = new QWidget(MatProp);
        gridLayoutWidget->setObjectName(QString::fromUtf8("gridLayoutWidget"));
        gridLayoutWidget->setGeometry(QRect(20, 30, 251, 121));
        gridLayout = new QGridLayout(gridLayoutWidget);
        gridLayout->setContentsMargins(10, 10, 10, 10);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        gridLayout->setContentsMargins(0, 0, 0, 0);
        doubleSpinBox_myu = new QDoubleSpinBox(gridLayoutWidget);
        doubleSpinBox_myu->setObjectName(QString::fromUtf8("doubleSpinBox_myu"));
        doubleSpinBox_myu->setDecimals(5);
        doubleSpinBox_myu->setMaximum(1e+06);
        doubleSpinBox_myu->setValue(0);

        gridLayout->addWidget(doubleSpinBox_myu, 0, 2, 1, 1);

        label = new QLabel(gridLayoutWidget);
        label->setObjectName(QString::fromUtf8("label"));
        QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(label->sizePolicy().hasHeightForWidth());
        label->setSizePolicy(sizePolicy);
        label->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label, 0, 0, 1, 1);

        label_2 = new QLabel(gridLayoutWidget);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        sizePolicy.setHeightForWidth(label_2->sizePolicy().hasHeightForWidth());
        label_2->setSizePolicy(sizePolicy);
        label_2->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_2, 1, 0, 1, 1);

        doubleSpinBox_rho = new QDoubleSpinBox(gridLayoutWidget);
        doubleSpinBox_rho->setObjectName(QString::fromUtf8("doubleSpinBox_rho"));
        doubleSpinBox_rho->setDecimals(5);
        doubleSpinBox_rho->setMinimum(0);
        doubleSpinBox_rho->setMaximum(1e+06);
        doubleSpinBox_rho->setSingleStep(0.1);
        doubleSpinBox_rho->setValue(1);

        gridLayout->addWidget(doubleSpinBox_rho, 1, 2, 1, 1);


        retranslateUi(Dialog_Fluid2D);
        QObject::connect(buttonBox, SIGNAL(accepted()), Dialog_Fluid2D, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), Dialog_Fluid2D, SLOT(reject()));

        QMetaObject::connectSlotsByName(Dialog_Fluid2D);
    } // setupUi

    void retranslateUi(QDialog *Dialog_Fluid2D)
    {
        Dialog_Fluid2D->setWindowTitle(QApplication::translate("Dialog_Fluid2D", "Dialog", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("Dialog_Fluid2D", "Element ID", 0, QApplication::UnicodeUTF8));
        MatProp->setTitle(QApplication::translate("Dialog_Fluid2D", "Material Property", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("Dialog_Fluid2D", "Myu", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("Dialog_Fluid2D", "Mass Deinsity", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class Dialog_Fluid2D: public Ui_Dialog_Fluid2D {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_DIALOG_FLUID2D_H
