noinst_LTLIBRARIES = libTIE.la
libTIE_la_FFLAGS = $(AM_FFLAGS) -fautomatic

libTIE_la_SOURCES = \
tensor_redux_Eij.F \
tensor_redux_Eij_new_Re_Com_1Mas.F\
tensor_redux_Eij_new_Re_Com.F 



if WITH_QUAD
libTIE_la_FFLAGS += -DWITH_QUAD
libTIE_la_SOURCES += \
tensor_redux_Eij_new_Re_Com_QUAD.F
endif
