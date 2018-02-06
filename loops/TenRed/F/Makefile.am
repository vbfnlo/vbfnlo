noinst_LTLIBRARIES = libTIF.la
libTIF_la_FFLAGS = $(AM_FFLAGS) -fautomatic
libTIF_la_SOURCES = \
tensor_redux_Fij_new_Re_Com.F 


if WITH_QUAD
libTIF_la_FFLAGS += -DWITH_QUAD
libTIF_la_SOURCES += \
tensor_redux_Fij_new_Re_Com_QUAD.F
endif

$(top_builddir)/include/%.mod: libTIF_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
