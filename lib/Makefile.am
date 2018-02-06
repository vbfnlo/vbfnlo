libVBFNLO_la_SOURCES = vbfnlo_main.F90
libVBFNLO_la_FCFLAGS = $(AM_FFLAGS)

if WITH_NLO
libVBFNLO_la_FCFLAGS += -DWITH_NLO
endif
if WITH_MPI
libVBFNLO_la_FCFLAGS += -DWITH_MPI
endif
if WITH_LOOPTOOLS
libVBFNLO_la_FCFLAGS += -DWITH_LT 
endif
if WITH_QUAD
libVBFNLO_la_FCFLAGS += -DWITH_QUAD
endif
if WITH_KK
libVBFNLO_la_SOURCES += trigger_cxx_linking.cpp
endif
if WITH_ROOT
libVBFNLO_la_SOURCES += trigger_cxx_linking.cpp
endif

pkglib_LTLIBRARIES = libVBFNLO.la

libVBFNLO_la_LIBADD = \
$(top_builddir)/PDFsets/libVBFNLOPDFsets.la \
$(top_builddir)/phasespace/libVBFNLOPhasespace.la \
$(top_builddir)/utilities/libVBFNLOUtilities.la \
$(top_builddir)/loops/libVBFNLOLoops.la \
$(top_builddir)/amplitudes/libVBFNLOAmplitudes.la 

$(top_builddir)/include/%.mod: libVBFNLO_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
