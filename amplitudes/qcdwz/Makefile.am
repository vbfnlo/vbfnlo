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

#EXTRA_DIST = temp.F.in
#BUILT_SOURCES = temp.F
#
#if WITH_LOOPTOOLS
#libqcdwz_la_FFLAGS  += -DWITH_LOOPTOOLS
## libqcdwz_la_FFLAGS  += -I$(LOOPTOOLS_DIR)/include -IFAFCme/
#endif
#
#if WITH_LOOPTOOLS
#temp.F: $(srcdir)/temp.F.in
#	@sed s%@LOOPTOOLSPATH@%"$(LOOPTOOLS_DIR)"% < $(srcdir)/temp.F.in > temp.F
#	@touch temp.F
#else
#temp.F: $(srcdir)/temp.F.in
#	@cp $(srcdir)/temp.F.in temp.F
#	@touch temp.F
#endif
#
#
#if WITH_LOOPTOOLS
#libqcdwz_la_FFLAGS+= -DWITH_LT
#libqcdwz_la_LIBADD = $(LOOPTOOLS_LDFLAGS) $(LOOPTOOLS_LIBS)
#endif
