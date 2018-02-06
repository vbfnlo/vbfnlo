noinst_LTLIBRARIES = libUtilities.la
libUtilities_la_SOURCES = \
GaugeTest.inc \
Pent_to_boxes.F\
pent45_noAbe.F\
GaugeTest.F\
MatrixElements.F \
bra_ket_subroutines.F \
Wardtest.F \
GluonSelf.F \
ten2cc.F

libUtilities_la_FFLAGS = $(AM_FFLAGS)

if WITH_FERMIONLOOPS
libUtilities_la_SOURCES += \
determinants.f
endif

if WITH_QUAD
libUtilities_la_FFLAGS += -DWITH_QUAD
libUtilities_la_SOURCES += \
Utilities_QUAD.F\
dotprod_QUAD.F
if WITH_QCD_JJ_PLUSX
libUtilities_la_SOURCES += \
GaugeTest_QUAD.F
endif
endif

$(top_builddir)/include/%.mod: libUtilities_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
