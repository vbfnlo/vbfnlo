noinst_LTLIBRARIES = libdiboson.la
libdiboson_la_SOURCES = \
tensor.inc \
vvv_anomal.inc \
genps.inc \
genpsr.inc \
DIBOStoww.F \
DIBOStoww_had.F \
DIBOStoww_anomal.F \
DIBOSwtowz.F\
DIBOSwtowz_had.F\
DIBOSwtowz_anomal.F\
ZZtozz.F \
ZZtozz_had.F \
DIBOSwtowa.F\
DIBOSwtowa_anomal.F\
ZAtoza.F \
AAtoaa.F \
m2s_qqVV.F \
m2s_qqWZ.F \
m2s_qqZZ.F \
m2s_qqWA.F \
m2s_qqZA.F \
m2s_qqAA.F \
m2s_qqVVj.F \
m2s_qqWZj.F \
m2s_qqZZj.F \
m2s_qqWAj.F \
m2s_qqZAj.F \
m2s_qqAAj.F \
qqVV.F \
qqWZ.F \
qqZZ.F \
qqWA.F \
qqZA.F \
qqAA.F \
ggWW.F \
ggZZ.F \
ggZA.F \
ggAA.F \
qqVVj.F \
qqWZj.F \
qqZZj.F \
qqWAj.F \
qqZAj.F \
qqAAj.F \
DIBOSwtowh.F \
DIBOSjwtowh.F \
DIBOSwtowh_anomal.F \
DIBOSjwtowh_anomal.F \
fl_vvjj.F \
m2s_qqWH.F \
m2s_qqWHj.F \
m2s_qqWHj2.F \
m2s_qqWHjj.F \
qqHWj_curr.F \
qqWH.F \
qqWHgg.F \
qqWHj.F \
qqWHj2.F \
qqWHqq.F \
DIBOStoqqll.F

if WITH_NLO
libdiboson_la_SOURCES += \
vv_amp_aux.F
endif

if WITH_MADGRAPH
libdiboson_la_SOURCES += \
mg2_ddww_had.F \
mg2_uuww_had.F \
mg2_udwz_had.F \
mg2_uuZZ_had.F \
mg2_ddZZ_had.F \
mg2_udwph.F \
mg2_udwmh.F \
mg2_udwhg.F \
mg2_duwhg.F \
mg2_udWHuu.F \
mg2_udWHdd.F \
mg2_udWHgg.F \
mg2_uuwwg.F \
mg2_ddwwg.F \
mg2_udwz.F \
mg2_duwz.F \
mg2_udwzg.F \
mg2_duwzg.F \
mg2_uuZZ.F \
mg2_ddZZ.F \
mg2_uuZZg.F \
mg2_ddZZg.F \
mg2_udwa.F \
mg2_duwa.F \
mg2_udwag.F \
mg2_duwag.F \
mg2_uuZA.F \
mg2_ddZA.F \
mg2_uuZAg.F \
mg2_ddZAg.F \
mg2_uuAA.F \
mg2_ddAA.F \
mg2_uuAAg.F \
mg2_ddAAg.F \
mg2_ddww.F \
mg2_uuww.F  
endif

libdiboson_la_FFLAGS = $(AM_FFLAGS)

if WITH_NLO
libdiboson_la_FFLAGS += -DWITH_NLO
endif

if WITH_MADGRAPH
libdiboson_la_FFLAGS += -DWITH_MADGRAPH
endif

if WITH_KK
libdiboson_la_FFLAGS += -DWITH_KK
endif

$(top_builddir)/include/%.mod: libdiboson_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
