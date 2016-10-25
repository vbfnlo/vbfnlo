noinst_LTLIBRARIES = libqcdwz.la
libqcdwz_la_FFLAGS = $(AM_FFLAGS)

libqcdwz_la_SOURCES = \
m2s_WZ2j.F \
m2s_WZ3j.F \
nloStuff.F \
Wjj.F


if WITH_QCDVV
libqcdwz_la_SOURCES += \
m2s_WW2j.F \
m2s_WW3j.F \
ug2dgwz.F \
uu2udwz.F 
endif

if WITH_NLO
if WITH_QCDVV
libqcdwz_la_SOURCES += \
qqwwqqVirt.F\
qqwzqqVirt.F\
qqwzggVirt.F
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
