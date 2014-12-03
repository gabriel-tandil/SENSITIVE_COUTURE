# -------------------------------------------------
# Project created by QtCreator 2010-05-21T13:41:55
# -------------------------------------------------
QT += opengl
TARGET = scalar2d
TEMPLATE = app
QMAKE_INCDIR += ../../include
QMAKE_LIBDIR += ../../lib
QMAKE_LIBS += -ldfm_core
SOURCES += main.cpp \
    mainwindow.cpp \
    glwidget.cpp \
    dialog_scalar2d.cpp
HEADERS += mainwindow.h \
    glwidget.h \
    dialog_scalar2d.h
FORMS += mainwindow.ui \
    dialog_scalar2d.ui
