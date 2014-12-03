/********************************************************************************
** Form generated from reading UI file 'dialog_solid2d.ui'
**
** Created: Thu Jul 15 15:46:06 2010
**      by: Qt User Interface Compiler version 4.6.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_DIALOG_SOLID2D_H
#define UI_DIALOG_SOLID2D_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
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

class Ui_Dialog_Solid2D
{
public:
    QDialogButtonBox *buttonBox;
    QComboBox *comboBox_IdEA;
    QLabel *label_4;
    QGroupBox *MatProp;
    QWidget *gridLayoutWidget;
    QGridLayout *gridLayout;
    QLabel *label_2;
    QLabel *label_3;
    QDoubleSpinBox *doubleSpinBox_rho;
    QDoubleSpinBox *doubleSpinBox_Poisson;
    QLabel *label;
    QDoubleSpinBox *doubleSpinBox_Young;
    QCheckBox *checkBox_nonlin;
    QLabel *label_5;

    void setupUi(QDialog *Dialog_Solid2D)
    {
        if (Dialog_Solid2D->objectName().isEmpty())
            Dialog_Solid2D->setObjectName(QString::fromUtf8("Dialog_Solid2D"));
        Dialog_Solid2D->resize(387, 281);
        Dialog_Solid2D->setSizeGripEnabled(false);
        Dialog_Solid2D->setModal(true);
        buttonBox = new QDialogButtonBox(Dialog_Solid2D);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setGeometry(QRect(30, 240, 341, 32));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Apply|QDialogButtonBox::Cancel|QDialogButtonBox::Ok);
        comboBox_IdEA = new QComboBox(Dialog_Solid2D);
        comboBox_IdEA->setObjectName(QString::fromUtf8("comboBox_IdEA"));
        comboBox_IdEA->setGeometry(QRect(120, 20, 111, 26));
        label_4 = new QLabel(Dialog_Solid2D);
        label_4->setObjectName(QString::fromUtf8("label_4"));
        label_4->setGeometry(QRect(40, 20, 81, 31));
        MatProp = new QGroupBox(Dialog_Solid2D);
        MatProp->setObjectName(QString::fromUtf8("MatProp"));
        MatProp->setGeometry(QRect(40, 50, 291, 181));
        gridLayoutWidget = new QWidget(MatProp);
        gridLayoutWidget->setObjectName(QString::fromUtf8("gridLayoutWidget"));
        gridLayoutWidget->setGeometry(QRect(0, 22, 291, 161));
        gridLayout = new QGridLayout(gridLayoutWidget);
        gridLayout->setContentsMargins(10, 10, 10, 10);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        gridLayout->setContentsMargins(0, 0, 0, 0);
        label_2 = new QLabel(gridLayoutWidget);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(label_2->sizePolicy().hasHeightForWidth());
        label_2->setSizePolicy(sizePolicy);
        label_2->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_2, 2, 0, 1, 1);

        label_3 = new QLabel(gridLayoutWidget);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        sizePolicy.setHeightForWidth(label_3->sizePolicy().hasHeightForWidth());
        label_3->setSizePolicy(sizePolicy);
        label_3->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_3, 3, 0, 1, 1);

        doubleSpinBox_rho = new QDoubleSpinBox(gridLayoutWidget);
        doubleSpinBox_rho->setObjectName(QString::fromUtf8("doubleSpinBox_rho"));
        doubleSpinBox_rho->setDecimals(3);
        doubleSpinBox_rho->setMaximum(100000);
        doubleSpinBox_rho->setValue(0);

        gridLayout->addWidget(doubleSpinBox_rho, 3, 2, 1, 1);

        doubleSpinBox_Poisson = new QDoubleSpinBox(gridLayoutWidget);
        doubleSpinBox_Poisson->setObjectName(QString::fromUtf8("doubleSpinBox_Poisson"));
        doubleSpinBox_Poisson->setDecimals(3);
        doubleSpinBox_Poisson->setMinimum(-1);
        doubleSpinBox_Poisson->setMaximum(0.5);
        doubleSpinBox_Poisson->setSingleStep(0.1);

        gridLayout->addWidget(doubleSpinBox_Poisson, 2, 2, 1, 1);

        label = new QLabel(gridLayoutWidget);
        label->setObjectName(QString::fromUtf8("label"));
        sizePolicy.setHeightForWidth(label->sizePolicy().hasHeightForWidth());
        label->setSizePolicy(sizePolicy);
        label->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label, 1, 0, 1, 1);

        doubleSpinBox_Young = new QDoubleSpinBox(gridLayoutWidget);
        doubleSpinBox_Young->setObjectName(QString::fromUtf8("doubleSpinBox_Young"));
        doubleSpinBox_Young->setDecimals(3);
        doubleSpinBox_Young->setMaximum(1e+06);
        doubleSpinBox_Young->setValue(0);

        gridLayout->addWidget(doubleSpinBox_Young, 1, 2, 1, 1);

        checkBox_nonlin = new QCheckBox(gridLayoutWidget);
        checkBox_nonlin->setObjectName(QString::fromUtf8("checkBox_nonlin"));

        gridLayout->addWidget(checkBox_nonlin, 4, 2, 1, 1);

        label_5 = new QLabel(gridLayoutWidget);
        label_5->setObjectName(QString::fromUtf8("label_5"));
        label_5->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_5, 4, 0, 1, 1);


        retranslateUi(Dialog_Solid2D);
        QObject::connect(buttonBox, SIGNAL(accepted()), Dialog_Solid2D, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), Dialog_Solid2D, SLOT(reject()));

        QMetaObject::connectSlotsByName(Dialog_Solid2D);
    } // setupUi

    void retranslateUi(QDialog *Dialog_Solid2D)
    {
        Dialog_Solid2D->setWindowTitle(QApplication::translate("Dialog_Solid2D", "Dialog", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("Dialog_Solid2D", "Element ID", 0, QApplication::UnicodeUTF8));
        MatProp->setTitle(QApplication::translate("Dialog_Solid2D", "Material Property", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("Dialog_Solid2D", "Poisson's ratio", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("Dialog_Solid2D", "Mass Density", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("Dialog_Solid2D", "young's modulus", 0, QApplication::UnicodeUTF8));
        checkBox_nonlin->setText(QString());
        label_5->setText(QApplication::translate("Dialog_Solid2D", "Nonlinear", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class Dialog_Solid2D: public Ui_Dialog_Solid2D {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_DIALOG_SOLID2D_H
