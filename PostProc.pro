#-------------------------------------------------
#
# Project created by QtCreator 2015-04-01T12:03:13
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = PostProc
TEMPLATE = app
INCLUDEPATH += /usr/include/eigen3
QMAKE_CXXFLAGS = -O3 -DEIGEN_DONT_VECTORIZE -DEIGEN_DONT_ALIGN

SOURCES += main.cpp\
        mainwindow.cpp \
    mesh_element.cc \
    mesh_file.cc \
    mesh.cc \
    meshplot.cpp \
    assembly.cc \
    element.cc \
    partition.cc \
    physlist.cpp \
    region.cc

HEADERS  += mainwindow.h \
    mesh_element.h \
    mesh_file.h \
    mesh.h \
    meshplot.h \
    assembly.h \
    element.h \
    partition.h \
    physlist.h \
    region.h
