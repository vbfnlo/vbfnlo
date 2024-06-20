noinst_LTLIBRARIES = libqcdwz.la
libqcdwz_la_FFLAGS = $(AM_FFLAGS)

libqcdwz_la_SOURCES = \
brakets.inc \
common.inc \
cacheZZ.inc \
m2s_WZ2j.F \
m2s_WZ3j.F \
nloStuff.F \
Wjj.F


if WITH_QCDVV
libqcdwz_la_SOURCES += \
m2s_WW2j.F \
m2s_WW3j.F \
uu2udwz.F \
ug2dgwz.F \
m2s_ZZ2j.F \
m2s_ZZ3j.F 
endif


if WITH_NLO
if WITH_QCDVV
libqcdwz_la_SOURCES += \
qqwwqqVirt.F\
qqwzqqVirt.F\
qqwzggVirt.F
libqcdwz_la_SOURCES += \
qqzzggVirt.F\
qqzzqqVirt.F
endif
endif

if WITH_QUAD
libqcdwz_la_FFLAGS += -DWITH_QUAD
endif

if WITH_NLO 
libqcdwz_la_FFLAGS += -DWITH_NLO
endif

if WITH_QCDV 
libqcdwz_la_FFLAGS += -DWITH_QCDV
endif

if WITH_QCDVV
libqcdwz_la_FFLAGS += -DWITH_QCDVV
endif

if WITH_VBF
libqcdwz_la_FFLAGS += -DWITH_VBF
endif

$(top_builddir)/include/%.mod: libqcdwz_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
