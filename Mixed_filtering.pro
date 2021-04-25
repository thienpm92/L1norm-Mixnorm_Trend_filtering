QT += core
QT -= gui

CONFIG += c++11

TARGET = Mixed_filtering
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app
INCLUDEPATH += "/usr/local/include"
INCLUDEPATH += "/usr/include/mlpack"
INCLUDEPATH += "/home/eric/Armadillo/armadillo-7.400.4"
INCLUDEPATH += "/home/eric/Armadillo/armadillo-7.400.4/include"
INCLUDEPATH += "/usr/local/include/eigen3"
INCLUDEPATH += "/home/eric/Ceres/ceres-solver-1.12.0new/config"
INCLUDEPATH += "/usr/include/suitesparse"

LIBS += `pkg-config --libs opencv`


LIBS += -L/home/eric/coinlib/lib
LIBS += -lmlpack

LIBS += -L/usr/local/lib
LIBS += -llapack -lblas -larmadillo
LIBS += -lmlpack
LIBS += -lrt
LIBS += -larmadillo
LIBS += -lboost_program_options
LIBS += -lboost_unit_test_framework
LIBS += -lboost_serialization

SOURCES += main.cpp \
    l1trend_filter.cpp \
    mixlp_filter.cpp \
    l1trend_admm.cpp

unix:!macx: LIBS += -L$$PWD/../../../../usr/lib/x86_64-linux-gnu/ -lsuperlu

INCLUDEPATH += $$PWD/../../../../usr/lib/x86_64-linux-gnu
DEPENDPATH += $$PWD/../../../../usr/lib/x86_64-linux-gnu

unix:!macx: LIBS += -L$$PWD/../../Armadillo/armadillo-7.400.4/ -larmadillo

INCLUDEPATH += $$PWD/../../Armadillo/armadillo-7.400.4
DEPENDPATH += $$PWD/../../Armadillo/armadillo-7.400.4

unix:!macx: LIBS += -L$$PWD/../../../../usr/lib/ -lopenblas

INCLUDEPATH += $$PWD/../../../../usr/include
DEPENDPATH += $$PWD/../../../../usr/include

unix:!macx: LIBS += -L$$PWD/../../../../usr/lib/ -llapack

INCLUDEPATH += $$PWD/../../../../usr/include
DEPENDPATH += $$PWD/../../../../usr/include

unix:!macx: LIBS += -L$$PWD/../../superLU/SuperLU_5.2.1/build/SRC/ -lsuperlu

INCLUDEPATH += $$PWD/../../superLU/SuperLU_5.2.1/build/SRC
DEPENDPATH += $$PWD/../../superLU/SuperLU_5.2.1/build/SRC

unix:!macx: PRE_TARGETDEPS += $$PWD/../../superLU/SuperLU_5.2.1/build/SRC/libsuperlu.a

unix:!macx: LIBS += -L$$PWD/../../superLU/SuperLU_5.2.1/build/CBLAS/ -lblas

INCLUDEPATH += $$PWD/../../superLU/SuperLU_5.2.1/build/CBLAS
DEPENDPATH += $$PWD/../../superLU/SuperLU_5.2.1/build/CBLAS

unix:!macx: PRE_TARGETDEPS += $$PWD/../../superLU/SuperLU_5.2.1/build/CBLAS/libblas.a


unix:!macx: LIBS += -L$$PWD/../../Ceres/ceres-solver-1.12.0new/build/lib/ -lceres

INCLUDEPATH += $$PWD/../../Ceres/ceres-solver-1.12.0new/build
DEPENDPATH += $$PWD/../../Ceres/ceres-solver-1.12.0new/build

unix:!macx: PRE_TARGETDEPS += $$PWD/../../Ceres/ceres-solver-1.12.0new/build/lib/libceres.a

unix:!macx: LIBS += -L$$PWD/../../../../usr/lib/x86_64-linux-gnu/ -lglog

INCLUDEPATH += $$PWD/../../../../usr/lib/x86_64-linux-gnu
DEPENDPATH += $$PWD/../../../../usr/lib/x86_64-linux-gnu

unix:!macx: LIBS += -L$$PWD/../../../../usr/lib/x86_64-linux-gnu/ -lgflags

INCLUDEPATH += $$PWD/../../../../usr/lib/x86_64-linux-gnu
DEPENDPATH += $$PWD/../../../../usr/lib/x86_64-linux-gnu

unix:!macx: LIBS += -L$$PWD/../../../../usr/lib/x86_64-linux-gnu/ -lsuitesparseconfig

INCLUDEPATH += $$PWD/../../../../usr/lib/x86_64-linux-gnu
DEPENDPATH += $$PWD/../../../../usr/lib/x86_64-linux-gnu

unix:!macx: LIBS += -L$$PWD/../../../../usr/lib/x86_64-linux-gnu/ -lspqr

INCLUDEPATH += $$PWD/../../../../usr/lib/x86_64-linux-gnu
DEPENDPATH += $$PWD/../../../../usr/lib/x86_64-linux-gnu

unix:!macx: LIBS += -L$$PWD/../../../../usr/lib/x86_64-linux-gnu/ -lcholmod

INCLUDEPATH += $$PWD/../../../../usr/lib/x86_64-linux-gnu
DEPENDPATH += $$PWD/../../../../usr/lib/x86_64-linux-gnu

HEADERS += \
    l1trend_filter.h \
    mixlp_filter.h \
    l1trend_admm.h
