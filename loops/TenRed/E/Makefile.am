noinst_LTLIBRARIES = libTIE.la
libTIE_la_FFLAGS = $(AM_FFLAGS) -fautomatic

libTIE_la_SOURCES = \
tensor_redux_Eij.F \
tensor_redux_Eij_new_Re_Com_1Mas.F\
tensor_redux_Eij_new_Re_Com.F 

if WITH_FERMIONLOOPS
libTIE_la_SOURCES += \
tensor_redux_Eij_new_Re_Com_R5_G.F\
tensor_redux_Eij_new_Re_Com_R5_G_Div.F
endif

if WITH_QUAD
libTIE_la_FFLAGS += -DWITH_QUAD
libTIE_la_SOURCES += \
tensor_redux_Eij_new_Re_Com_QUAD.F\
tensor_redux_Eij_new_Re_Com_R5_G_QUAD.F\
QUAD_tensor_redux_Eij_new_Re_Com_1Mas.F 
endif

$(top_builddir)/include/%.mod: libTIE_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
