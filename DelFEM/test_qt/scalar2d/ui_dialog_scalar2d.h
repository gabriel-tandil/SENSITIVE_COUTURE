/********************************************************************************
** Form generated from reading UI file 'dialog_scalar2d.ui'
**
** Created: Sat Jan 29 16:27:39 2011
**      by: Qt User Interface Compiler version 4.6.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_DIALOG_SCALAR2D_H
#define UI_DIALOG_SCALAR2D_H

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

class Ui_Dialog_Scalar2D
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
    QDoubleSpinBox *doubleSpinBox_source;
    QDoubleSpinBox *doubleSpinBox_capacity;
    QLabel *label;
    QDoubleSpinBox *doubleSpinBox_alpha;

    void setupUi(QDialog *Dialog_Scalar2D)
    {
        if (Dialog_Scalar2D->objectName().isEmpty())
            Dialog_Scalar2D->setObjectName(QString::fromUtf8("Dialog_Scalar2D"));
        Dialog_Scalar2D->resize(387, 281);
        Dialog_Scalar2D->setSizeGripEnabled(false);
        Dialog_Scalar2D->setModal(true);
        buttonBox = new QDialogButtonBox(Dialog_Scalar2D);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setGeometry(QRect(30, 240, 341, 32));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Apply|QDialogButtonBox::Cancel|QDialogButtonBox::Ok);
        buttonBox->setCenterButtons(false);
        comboBox_IdEA = new QComboBox(Dialog_Scalar2D);
        comboBox_IdEA->setObjectName(QString::fromUtf8("comboBox_IdEA"));
        comboBox_IdEA->setGeometry(QRect(120, 20, 111, 26));
        label_4 = new QLabel(Dialog_Scalar2D);
        label_4->setObjectName(QString::fromUtf8("label_4"));
        label_4->setGeometry(QRect(40, 20, 81, 31));
        MatProp = new QGroupBox(Dialog_Scalar2D);
        MatProp->setObjectName(QString::fromUtf8("MatProp"));
        MatProp->setGeometry(QRect(40, 60, 291, 171));
        gridLayoutWidget = new QWidget(MatProp);
        gridLayoutWidget->setObjectName(QString::fromUtf8("gridLayoutWidget"));
        gridLayoutWidget->setGeometry(QRect(30, 30, 251, 121));
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

        doubleSpinBox_source = new QDoubleSpinBox(gridLayoutWidget);
        doubleSpinBox_source->setObjectName(QString::fromUtf8("doubleSpinBox_source"));
        doubleSpinBox_source->setDecimals(3);
        doubleSpinBox_source->setMinimum(-1e+06);
        doubleSpinBox_source->setMaximum(1e+06);

        gridLayout->addWidget(doubleSpinBox_source, 3, 2, 1, 1);

        doubleSpinBox_capacity = new QDoubleSpinBox(gridLayoutWidget);
        doubleSpinBox_capacity->setObjectName(QString::fromUtf8("doubleSpinBox_capacity"));
        doubleSpinBox_capacity->setDecimals(3);
        doubleSpinBox_capacity->setMinimum(0);
        doubleSpinBox_capacity->setMaximum(1e+06);
        doubleSpinBox_capacity->setSingleStep(0.1);
        doubleSpinBox_capacity->setValue(1);

        gridLayout->addWidget(doubleSpinBox_capacity, 2, 2, 1, 1);

        label = new QLabel(gridLayoutWidget);
        label->setObjectName(QString::fromUtf8("label"));
        sizePolicy.setHeightForWidth(label->sizePolicy().hasHeightForWidth());
        label->setSizePolicy(sizePolicy);
        label->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label, 1, 0, 1, 1);

        doubleSpinBox_alpha = new QDoubleSpinBox(gridLayoutWidget);
        doubleSpinBox_alpha->setObjectName(QString::fromUtf8("doubleSpinBox_alpha"));
        doubleSpinBox_alpha->setDecimals(3);
        doubleSpinBox_alpha->setMaximum(1e+06);
        doubleSpinBox_alpha->setValue(0);

        gridLayout->addWidget(doubleSpinBox_alpha, 1, 2, 1, 1);


        retranslateUi(Dialog_Scalar2D);
        QObject::connect(buttonBox, SIGNAL(accepted()), Dialog_Scalar2D, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), Dialog_Scalar2D, SLOT(reject()));

        QMetaObject::connectSlotsByName(Dialog_Scalar2D);
    } // setupUi

    void retranslateUi(QDialog *Dialog_Scalar2D)
    {
        Dialog_Scalar2D->setWindowTitle(QApplication::translate("Dialog_Scalar2D", "Dialog", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("Dialog_Scalar2D", "Element ID", 0, QApplication::UnicodeUTF8));
        MatProp->setTitle(QApplication::translate("Dialog_Scalar2D", "Material Property", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("Dialog_Scalar2D", "Calpacity", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("Dialog_Scalar2D", "Source", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("Dialog_Scalar2D", "Alpha", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class Dialog_Scalar2D: public Ui_Dialog_Scalar2D {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_DIALOG_SCALAR2D_H
