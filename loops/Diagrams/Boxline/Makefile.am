noinst_LTLIBRARIES = libDiagramsBoxL.la
libDiagramsBoxL_la_SOURCES = \
boxline.F \
boxlineABE.F\
boxlineABEDiv.F\
boxlineABETotal.F\
boxlineNoAbe.F\
boxlineNoAbe_p.F\
boxlineNoAbeDiv.F\
boxlineNoAbe_pDiv.F\
boxlineNoAbeTotal.F\
boxlineNoAbeTotal_p.F


libDiagramsBoxL_la_FFLAGS = $(AM_FFLAGS)

if WITH_QUAD
if WITH_QCD_JJ_PLUSX
libDiagramsBoxL_la_FFLAGS += -DWITH_QUAD
libDiagramsBoxL_la_SOURCES += \
XDGetM0_QUAD_Dble.F \
boxline_QUAD.F \
boxlineABE_QUAD.F \
boxlineNoAbe_QUAD.F
endif
endif

$(top_builddir)/include/%.mod: libDiagramsBoxL_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
