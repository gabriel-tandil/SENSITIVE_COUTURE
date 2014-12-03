#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

class GLWidget;
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

private:
    Ui::MainWindow *ui;
    GLWidget* glWidget;
};

#endif // MAINWINDOW_H
