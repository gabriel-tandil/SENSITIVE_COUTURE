#ifndef DIALOG_SOLID2D_VIEWSETTING_H
#define DIALOG_SOLID2D_VIEWSETTING_H

#include <QDialog>

namespace Ui {
    class Dialog_Solid2D_ViewSetting;
}

class Dialog_Solid2D_ViewSetting : public QDialog {
    Q_OBJECT
public:
    Dialog_Solid2D_ViewSetting(QWidget *parent = 0);
    ~Dialog_Solid2D_ViewSetting();

protected:
    void changeEvent(QEvent *e);

private:
    Ui::Dialog_Solid2D_ViewSetting *ui;
};

#endif // DIALOG_SOLID2D_VIEWSETTING_H
