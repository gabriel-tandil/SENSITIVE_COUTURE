# -------------------------------------------------
# Project created by QtCreator 2010-05-21T13:41:55
# -------------------------------------------------
QT += opengl
TARGET = solid2d
TEMPLATE = app
QMAKE_INCDIR += ../../include
QMAKE_LIBDIR += ../../lib
QMAKE_LIBS += -ldfm_core
SOURCES += main.cpp \
    mainwindow.cpp \
    glwidget_solid2d.cpp \
    dialog_solid2d.cpp \
    dialog_solid2d_viewsetting.cpp
HEADERS += mainwindow.h \
    glwidget_solid2d.h \
    dialog_solid2d.h \
    dialog_solid2d_viewsetting.h
FORMS += mainwindow.ui \
    dialog_solid2d.ui \
    dialog_solid2d_viewsetting.ui
