noinst_LTLIBRARIES = libvvj.la
libvvj_la_SOURCES = \
tensor.inc \
genps.inc \
nexternal.inc \
tensorlep.inc \
leptens.F \
m2s_qqVVj2.F \
m2s_qqZZj2.F \
m2s_qqVVjj.F \
m2s_qqZZjj.F \
m2s_qqwmag.F \
m2s_qqwmagj.F \
m2s_qqwpag.F \
m2s_qqwpagj.F \
m2s_qqwmzg.F \
m2s_qqwmzgj.F \
m2s_qqwpzg.F \
m2s_qqwpzgj.F \
qqVVgg.F \
qqVVqq.F \
qqZZgg.F \
qqZZqq.F \
ggWWg.F \
ggZZg.F \
qqWWj_curr.F \
qqZZj_curr.F 


if WITH_NLO
libvvj_la_SOURCES += \
m2s_qqwmagvirt.F \
m2s_qqwpagvirt.F \
m2s_qqwmzgvirt.F \
m2s_qqwpzgvirt.F 
endif

if WITH_MADGRAPH
libvvj_la_SOURCES += \
mg2_uuWWgg.F \
mg2_ddWWgg.F \
mg2_uuWWuu.F \
mg2_ddWWdd.F \
mg2_udWWud.F \
mg2_duWWdu.F \
mg2_uuZZgg.F \
mg2_ddZZgg.F \
mg2_uuZZuu.F \
mg2_ddZZdd.F \
mg2_udZZud.F \
mg2_duZZdu.F 
endif

libvvj_la_FFLAGS = $(AM_FFLAGS)

if WITH_NLO
libvvj_la_FFLAGS += -DWITH_NLO
endif

if WITH_MADGRAPH
libvvj_la_FFLAGS += -DWITH_MADGRAPH
endif

$(top_builddir)/include/%.mod: libvvj_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
