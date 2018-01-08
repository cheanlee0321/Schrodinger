TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    schrodinger.cpp

HEADERS += \
    schrodinger.h
    Parameter.h

QMAKE_CXXFLAGS += -fopenmp
LIBS += -fopenmp
