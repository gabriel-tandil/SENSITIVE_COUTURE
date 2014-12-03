# -------------------------------------------------
# Project created by QtCreator 2010-05-21T13:41:55
# -------------------------------------------------
QT += opengl
TARGET = fluid2d
TEMPLATE = app
QMAKE_INCDIR += ../../include
QMAKE_LIBDIR += ../../lib
QMAKE_LIBS += -ldfm_core
SOURCES += main.cpp \
    mainwindow.cpp \
    glwidget.cpp \
    dialog_fluid2d.cpp \
    ../solid2d/dialog_solid2d_viewsetting.cpp
HEADERS += mainwindow.h \
    glwidget.h \
    dialog_fluid2d.h \
    ../solid2d/dialog_solid2d_viewsetting.h
FORMS += mainwindow.ui \
    dialog_fluid2d.ui \
    ../solid2d/dialog_solid2d_viewsetting.ui
