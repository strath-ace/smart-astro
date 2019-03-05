#-------------------------------------------------
#
# Project created by QtCreator 2019-01-14T12:30:27
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = smart-astro
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

CONFIG += c++11

SOURCES += \
        src/Controller/parameterinputcontroller.cpp \
	src/View/answerdialog.cpp \
	src/View/genericfunctiondialog.cpp \
	src/View/mainwindow.cpp \ 
        src/View/parameterinputdialog.cpp \
	src/main.cpp \
	src/Model/testModel.cpp \

HEADERS += \
	src/Controller/parameterinputcontroller.h \
        src/View/answerdialog.h \
	src/View/genericfunctiondialog.h \
	src/View/mainwindow.h \ 
	src/View/parameterinputdialog.h \
	src/Model/testModel.h \

FORMS += \
	src/View/answerdialog.ui \
	src/View/genericfunctiondialog.ui \
	src/View/mainwindow.ui \ 
        src/View/parameterinputdialog.ui \

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
