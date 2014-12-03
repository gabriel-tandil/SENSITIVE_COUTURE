#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

class GLWidget_Solid2d;
namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();

protected:
    void changeEvent(QEvent *e);
public slots:
    void setNewProblem();
    void showDialogMatProp();
    void showDialogViewSetting();

    void about();

private:
    void createActions();
    void createMenus();

private:
    Ui::MainWindow *ui;
    GLWidget_Solid2d* glWidget;

    QMenu *menuEdit;
    QMenu *menuView;
    QMenu *menuHelp;

    QAction *actMatPropDlg;
    QAction *actAbout;
    QAction *actAboutQt;
};

#endif // MAINWINDOW_H
