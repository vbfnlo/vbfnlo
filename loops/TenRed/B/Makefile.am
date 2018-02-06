noinst_LTLIBRARIES = libTIB.la
libTIB_la_SOURCES = \
tensor_redux_BijLUde1M.F\
tensor_redux_BijLUde.F\
tensor_redux_BijLUde1M_Div.F\
tensor_redux_BijLUde_Div.F

if WITH_FERMIONLOOPS
libTIB_la_SOURCES += \
tensor_redux_BijLUde_G.F\
tensor_redux_BijLUde_G_Div.F\
tensor_redux_BijLUde_R3_G.F\
tensor_redux_BijLUde_R3_G_Div.F\
BijForGram.f\
BijForGram_M.f\
BijForGram_G_Num.f
endif


libTIB_la_FFLAGS = $(AM_FFLAGS)

if WITH_QUAD
libTIB_la_FFLAGS += -DWITH_QUAD
libTIB_la_SOURCES += \
tensor_redux_BijLUde_R3_G_QUAD.F
endif

$(top_builddir)/include/%.mod: libTIB_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
