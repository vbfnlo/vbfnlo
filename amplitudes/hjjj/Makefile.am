noinst_LTLIBRARIES = libhjjj.la
libhjjj_la_SOURCES = \
coupl_red.inc \
genps.inc \
color.inc \
m2s_qqh3j.F  \
qqhqq.F       \
m2s_qqh4j.F    \
qqh4q_sub.F  \
qqh4q_me.F \
qqh2q2g_me.F

if WITH_NLO
libhjjj_la_SOURCES += \
hjjj_amp_aux.F\
coeff.F        \
dipsub.F     \
dipolesub.F  
endif

if WITH_MADGRAPH
libhjjj_la_SOURCES += \
qqhqqgg_mgfast.F \
qqh4q_mg.F \
qqhqqgg_cc.F \
qqhqqgg_mg.F \
dc_usbbh.F \
us_usbbh.F \
us_usuuh.F
endif

libhjjj_la_FFLAGS = $(AM_FFLAGS)

if WITH_NLO
libhjjj_la_FFLAGS += -DWITH_NLO
endif

if WITH_MADGRAPH
libhjjj_la_FFLAGS += -DWITH_MADGRAPH
endif

$(top_builddir)/include/%.mod: libhjjj_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
