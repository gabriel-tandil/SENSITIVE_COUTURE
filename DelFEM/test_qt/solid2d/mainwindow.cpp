#include <QMenuBar>
#include <QMessageBox>

#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "glwidget_solid2d.h"

#include "dialog_solid2d.h"
#include "dialog_solid2d_viewsetting.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    connect(ui->pushButton_prob,    SIGNAL(clicked()),  this,SLOT(setNewProblem()));
    connect(ui->pushButton_matprop, SIGNAL(clicked()),  this,SLOT(showDialogMatProp()));
    connect(ui->pushButton_viewset, SIGNAL(clicked()),  this,SLOT(showDialogViewSetting()));

    glWidget = new GLWidget_Solid2d;
    ui->scrollArea->setWidget(glWidget);

    createActions();
    createMenus();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::createActions()
{
    actMatPropDlg = new QAction(tr("&Material Property Dialog"),this);
    connect(actMatPropDlg,  SIGNAL(triggered()), this, SLOT(showDialogMatProp()));

    actAbout = new QAction(tr("&About"), this);
    actAbout->setStatusTip(tr("Show the application's About box"));
    connect(actAbout, SIGNAL(triggered()), this, SLOT(about()));

    actAboutQt = new QAction(tr("About &Qt"), this);
    actAboutQt->setStatusTip(tr("Show the Qt library's About box"));
    connect(actAboutQt, SIGNAL(triggered()), qApp, SLOT(aboutQt()));
}

void MainWindow::createMenus()
{
    menuEdit = menuBar()->addMenu(tr("&Edit"));
    menuEdit->addAction(actMatPropDlg);

    menuView = menuBar()->addMenu(tr("&View"));
//    fileMenu->addAction(actViewCoord);
//    fileMenu->addAction(openAct);
//    fileMenu->addAction(saveAct);
//    fileMenu->addAction(saveAsAct);
//    fileMenu->addSeparator();
//    fileMenu->addAction(exitAct);
    menuHelp = menuBar()->addMenu(tr("&Help"));
    menuHelp->addAction(actAbout);
    menuHelp->addAction(actAboutQt);
}

void MainWindow::changeEvent(QEvent *e)
{
    QMainWindow::changeEvent(e);
    switch (e->type()) {
    case QEvent::LanguageChange:
        ui->retranslateUi(this);
        break;
    default:
        break;
    }
}


void MainWindow::about()
{
   QMessageBox::about(this, tr("About SDI"),
            tr("The <b>SDI</b> example demonstrates how to write single "
               "document interface applications using Qt."));
}


void MainWindow::setNewProblem(){
    std::cout << "setNewProblem()" << std::endl;
    glWidget->SetNewProblem();
}

void MainWindow::showDialogMatProp(){
    std::cout << "show Dialog Materila Property()" << std::endl;
    Dialog_Solid2D* dialog = new Dialog_Solid2D(glWidget->world,glWidget->solid,this);
    ((QDialog*)dialog)->setModal(true);
    ((QDialog*)dialog)->exec();
    delete dialog;
}

void MainWindow::showDialogViewSetting(){
    std::cout << "show Dialog View Setting()" << std::endl;
    Dialog_Solid2D_ViewSetting* dialog = new Dialog_Solid2D_ViewSetting(this);
    ((QDialog*)dialog)->setModal(true);
    ((QDialog*)dialog)->exec();
    delete dialog;
}
