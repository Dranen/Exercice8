TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

TARGET = ../Exercice8/Exercice8

QMAKE_CXXFLAGS += -std=gnu++11
QMAKE_CXXFLAGS += -O3
QMAKE_CXXFLAGS += -march=native

SOURCES += \
    SchroedingerMain.cpp

