SUBDIRS = 

pkglib_LTLIBRARIES = libVBFNLOAmplitudes.la 
libVBFNLOAmplitudes_la_SOURCES = \
amplitudes.F
libVBFNLOAmplitudes_la_FFLAGS = $(AM_FFLAGS)
libVBFNLOAmplitudes_la_LIBADD = 

if WITH_DIBOSON
SUBDIRS += diboson
libVBFNLOAmplitudes_la_FFLAGS += -DWITH_DIBOSON
libVBFNLOAmplitudes_la_LIBADD += diboson/libdiboson.la
endif

if WITH_DIBOSONJET
if !WITH_DIBOSON
SUBDIRS += diboson
libVBFNLOAmplitudes_la_LIBADD += diboson/libdiboson.la
endif
SUBDIRS += vvj
libVBFNLOAmplitudes_la_FFLAGS += -DWITH_DIBOSONJET
libVBFNLOAmplitudes_la_LIBADD += vvj/libvvj.la
endif

if WITH_QCDV
SUBDIRS += qcdwz
libVBFNLOAmplitudes_la_FFLAGS += -DWITH_QCDV
libVBFNLOAmplitudes_la_LIBADD += qcdwz/libqcdwz.la
if WITH_QCDVV
libVBFNLOAmplitudes_la_FFLAGS += -DWITH_QCDVV
endif
endif

if WITH_TRIBOSON
SUBDIRS += triboson
libVBFNLOAmplitudes_la_FFLAGS += -DWITH_TRIBOSON
libVBFNLOAmplitudes_la_LIBADD += triboson/libtriboson.la
endif

if WITH_TRIBOSONJET
if !WITH_TRIBOSON
SUBDIRS += triboson
libVBFNLOAmplitudes_la_LIBADD += triboson/libtriboson.la
endif
libVBFNLOAmplitudes_la_FFLAGS += -DWITH_TRIBOSONJET
endif

if WITH_VBF
SUBDIRS += vvjj
libVBFNLOAmplitudes_la_FFLAGS += -DWITH_VBF
libVBFNLOAmplitudes_la_LIBADD += vvjj/libvbf.la
endif

if WITH_HJJJ
SUBDIRS += hjjj
libVBFNLOAmplitudes_la_FFLAGS += -DWITH_HJJJ
libVBFNLOAmplitudes_la_LIBADD += hjjj/libhjjj.la
endif

if WITH_GGF 
SUBDIRS += ggf
libVBFNLOAmplitudes_la_FFLAGS += -DWITH_GGF
libVBFNLOAmplitudes_la_LIBADD += ggf/libggf.la
endif

if WITH_NLO
libVBFNLOAmplitudes_la_FFLAGS += -DWITH_NLO
endif

if WITH_KK
libVBFNLOAmplitudes_la_FFLAGS += -DWITH_KK
endif

SUBDIRS += Spin2
libVBFNLOAmplitudes_la_LIBADD += Spin2/libSpin2.la

libVBFNLOAmplitudes_la_LIBADD += $(top_builddir)/helas/libHELAS.la

$(top_builddir)/include/%.mod: libVBFNLOAmplitudes_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
