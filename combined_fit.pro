TEMPLATE = app
CONFIG += c++14
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_LFLAGS += -O3 -std=c++1y
INCLUDEPATH = /usr/local/include/eigen3
OBJECTS_DIR = objects/
QMAKE_CXXFLAGS += -O3 -std=c++1y
TARGET = comb_fit
LIBS += -lecf

SOURCES += \
    src/RandomNumbers.cc \
    src/GaussMutOp.cc \
    src/MyFloatingPoint.cc \
    src/main.cpp \
    src/GaussianFitFunctor.cpp \
    src/EvalOp.cc \
    src/basis.cpp

HEADERS += \
    src/RandomNumbers.h \
    src/GaussMutOp.h \
    src/MyFloatingPoint.h \
    src/GaussianFitFunctor.h \
    src/EvalOp.h \
    src/basis.h
