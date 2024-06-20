pkglib_LTLIBRARIES = libVBFNLOPhasespace.la
libVBFNLOPhasespace_la_SOURCES = \
multires.inc \
phasespace.F \
ps_tools.F \
ps_vnj.F \
ps_vvnj.F \
ps_qqHHj.F \
ps_qqVAj.F 

if WITH_DIBOSON
libVBFNLOPhasespace_la_SOURCES += \
ps_qqVVj.F \
ps_qqVHj.F \
ps_qqVj.F \
ps_qqWZ.F \
ps_qqAAj.F 
endif

if WITH_DIBOSONJET
if !WITH_DIBOSON
libVBFNLOPhasespace_la_SOURCES += \
ps_qqVVj.F \
ps_qqVHj.F \
ps_qqVj.F
endif
endif

if WITH_TRIBOSON
libVBFNLOPhasespace_la_SOURCES += \
ps_qqVVVj.F \
ps_ZZZj.F \
ps_WWW.F \
ps_qqVVAj.F \
ps_qqVAAj.F \
ps_qqAAAj.F
endif

if WITH_TRIBOSONJET
if !WITH_TRIBOSON
libVBFNLOPhasespace_la_SOURCES += \
ps_qqVAAj.F 
endif
endif

if WITH_GGF 
libVBFNLOPhasespace_la_SOURCES += \
ps_ggf.F
endif

libVBFNLOPhasespace_la_FFLAGS = $(AM_FFLAGS)
libVBFNLOPhasespace_la_CPPFLAGS = $(AM_CPPFLAGS)

if WITH_DIBOSON
libVBFNLOPhasespace_la_FFLAGS += -DWITH_DIBOSON
endif

if WITH_DIBOSONJET
libVBFNLOPhasespace_la_FFLAGS += -DWITH_DIBOSONJET
endif

if WITH_TRIBOSON
libVBFNLOPhasespace_la_FFLAGS += -DWITH_TRIBOSON
endif

if WITH_TRIBOSONJET
libVBFNLOPhasespace_la_FFLAGS += -DWITH_TRIBOSONJET
endif

if WITH_QCDVV
libVBFNLOPhasespace_la_FFLAGS += -DWITH_QCDVV
endif

if WITH_QCDV
libVBFNLOPhasespace_la_FFLAGS += -DWITH_QCDV
endif

if WITH_KK
libVBFNLOPhasespace_la_FFLAGS += -DWITH_KK
libVBFNLOPhasespace_la_CPPFLAGS += $(GSLINCLUDE)
libVBFNLOPhasespace_la_LIBADD = $(GSLLIBS)
libVBFNLOPhasespace_la_SOURCES += \
ps_tribos_kk.F
endif

$(top_builddir)/include/%.mod: libVBFNLOPhasespace_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
