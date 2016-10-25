pkglib_LTLIBRARIES = libVBFNLO.la
libVBFNLO_la_SOURCES = \
	vbfnlo_main.F
libVBFNLO_la_FFLAGS = $(AM_FFLAGS)

if WITH_NLO
libVBFNLO_la_FFLAGS += -DWITH_NLO
endif
if WITH_LOOPTOOLS
libVBFNLO_la_FFLAGS += -DWITH_LT 
endif
if WITH_QUAD
libVBFNLO_la_FFLAGS += -DWITH_QUAD
endif

if WITH_KK
libVBFNLO_la_SOURCES += trigger_cxx_linking.cpp
endif

if WITH_ROOT
libVBFNLO_la_SOURCES += trigger_cxx_linking.cpp
endif

libVBFNLO_la_LIBADD = \
$(top_builddir)/amplitudes/libVBFNLOAmplitudes.la \
$(top_builddir)/helas/libHELAS.la \
$(top_builddir)/PDFsets/libVBFNLOPDFsets.la \
$(top_builddir)/phasespace/libVBFNLOPhasespace.la \
$(top_builddir)/utilities/libVBFNLOUtilities.la \
$(top_builddir)/utilities/libVBFspecific.la \
$(top_builddir)/loops/libVBFNLOLoops.la

if WITH_GGFLO
pkglib_LTLIBRARIES += libGGFLO.la
libGGFLO_la_SOURCES = \
	ggflo_main.F
libGGFLO_la_FFLAGS = $(AM_FFLAGS)
if WITH_LOOPTOOLS
libGGFLO_la_FFLAGS += -DWITH_LT 
endif

if WITH_ROOT
libGGFLO_la_SOURCES += trigger_cxx_linking.cpp
endif

libGGFLO_la_LIBADD = \
$(top_builddir)/amplitudes/libGGFLOAmplitudes.la \
$(top_builddir)/helas/libHELAS.la \
$(top_builddir)/PDFsets/libVBFNLOPDFsets.la \
$(top_builddir)/phasespace/libVBFNLOPhasespace.la \
$(top_builddir)/utilities/libVBFNLOUtilities.la \
$(top_builddir)/utilities/libGGFspecific.la \
$(top_builddir)/loops/libVBFNLOLoops.la
endif
