noinst_LTLIBRARIES = libtriboson.la
libtriboson_la_SOURCES = \
ud_www.inc \
tensor.inc \
du_zzw.inc \
genps.inc \
ud_zzw.inc \
du_www.inc \
genpsr.inc \
dd_wwz.inc \
uu_wwz.inc \
m2s_qqZZWj.F  \
m2s_qqZZZ.F   \
m2s_qqWWZ.F   \
m2s_qqZZZj.F  \
m2s_qqWWZj.F  \
m2s_qqZZW.F   \
m2s_qqWWW.F   \
m2s_qqWWWj.F  \
m2s_qqWWA.F \
m2s_qqWWAj.F \
m2s_qqZZA.F \
m2s_qqZZAj.F \
m2s_qqWZA.F \
m2s_qqWZAj.F \
m2s_qqWAA.F \
m2s_qqWAAj.F \
m2s_qqZAA.F \
m2s_qqZAAj.F \
m2s_qqAAA.F \
m2s_qqAAAj.F \
qqWWZ.F   \
qqZZW.F   \
qqZZZ.F   \
qqZZZ1.F  \
qqZZWj.F  \
qqWWZj.F  \
qqZZZj.F  \
qqZZZj1.F  \
qqWWW.F   \
qqWWWj.F  \
qqVVA.F   \
qqVVAj.F \
qqZZA.F   \
qqZZAj.F \
qqWZA.F   \
qqWZAj.F \
qqAAW.F   \
qqAAWj.F   \
qqZAA.F   \
qqZAAj.F   \
qqAAA.F   \
qqAAAj.F   \
ZZWtozzw.F \
WWZtoww.F  \
WWZtowwz.F \
ZZZtozz.F  \
ZZWtowz.F  \
WWZtowz.F  \
ZZZtozzz.F \
WWWtoww.F  \
WWWtowww.F \
WWAtowa.F  \
WWAtowwa.F \
ZZAtozae.F \
ZZAtozam.F \
ZZAtozza.F \
ZZnAtozza.F \
ZZAtozz.F \
ZZnAtozz.F \
WZAtowa.F \
WZAtowza.F \
WZAtowz.F \
WZAtoza.F \
WAAtowaa.F \
WAAtowa.F \
ZAAtozaa.F \
ZAAtoza.F \
WWZ_ltens_had.F  \
WWZtoww_anomal.F  \
WWZtowz_anomal.F  \
WWZtowwz_anomal.F \
ZZZtozz_anomal.F  \
ZZW_ltens_had.F  \
ZZWtowz_anomal.F  \
ZZWtozzw_anomal.F \
WWW_ltens_had.F  \
WWWtoww_anomal.F  \
WWWtowww_anomal.F \
ZZZ_ltens_had.F  \
ZZZtozzz_anomal.F \
WWA_ltens_had.F  \
WWAtowa_anomal.F  \
WWAtowwa_anomal.F \
ZZAtozae_anomal.F \
ZZAtozam_anomal.F \
ZZAtozza_anomal.F \
ZZAtozz_anomal.F \
WZAtowa_anomal.F \
WZAtowza_anomal.F \
WZAtowz_anomal.F \
WZAtoza_anomal.F \
WAAtowaa_anomal.F \
WAAtowa_anomal.F \
ZAAtozaa_anomal.F \
ZAAtoza_anomal.F \
AAAtoaaa_anomal.F \
WZA_ltens_had.F \
ZZA_ltens_had.F


if WITH_MADGRAPH
libtriboson_la_SOURCES += \
mg2_ddZZZ.F   \
mg2_udZZWg.F  \
mg2_ddZZZg.F  \
mg2_uuWWZ.F   \
mg2_duZZW.F   \
mg2_uuWWZg.F  \
mg2_ddWWZ.F   \
mg2_duZZWg.F  \
mg2_uuZZZ.F   \
mg2_ddWWZg.F  \
mg2_udZZW.F   \
mg2_uuZZZg.F  \
mg2_udbarWWW.F  \
mg2_udbarWWWg.F \
mg2_dubarWWW.F  \
mg2_dubarWWWg.F \
mg2_uuWWA.F \
mg2_ddWWA.F \
mg2_uuWWAg.F \
mg2_ddWWAg.F \
mg2_uuZZA.F \
mg2_ddZZA.F \
mg2_uuZZnA.F \
mg2_ddZZnA.F \
mg2_uuZZAg.F \
mg2_ddZZAg.F \
mg2_udWZA.F \
mg2_udWZAg.F \
mg2_udWAA.F \
mg2_udWAAg.F \
mg2_uuZAA.F \
mg2_ddZAA.F \
mg2_uuZAAg.F \
mg2_ddZAAg.F \
mg2_uuZnAA.F \
mg2_ddZnAA.F \
mg2_uuZnAAg.F \
mg2_ddZnAAg.F \
mg2_uuAAA.F \
mg2_ddAAA.F \
mg2_uuAAAg.F \
mg2_ddAAAg.F \
mg2_WZA_ltens_had.F \
mg2_ZZA_ltens_had.F \
mg2_WWA_ltens_had.F \
mg2_WWW_ltens_had.F \
mg2_WWZ_ltens_had.F \
mg2_ZZW_ltens_had.F \
mg2_ZZZ_ltens_had.F
endif


libtriboson_la_FFLAGS = $(AM_FFLAGS)
libtriboson_la_CPPFLAGS = $(AM_CPPFLAGS)

if WITH_NLO
libtriboson_la_FFLAGS += -DWITH_NLO
endif

if WITH_MADGRAPH
libtriboson_la_FFLAGS += -DWITH_MADGRAPH
endif


if WITH_KK
libtriboson_la_SOURCES += \
WWZtowwz_kk.F \
ZZWtozzw_kk.F \
WWWtowww_kk.F \
jwwwxxtest.F \
jw3wxxtest.F \
w3w3xxtest.F \
wwwwxxtest.F \
mg2_udZZW_mhlm.F \
mg2_duZZW_mhlm.F \
mg2_uuWWZ_mhlm.F \
mg2_ddWWZ_mhlm.F \
mg2_udbarWWW_mhlm.F \
mg2_dubarWWW_mhlm.F

libtriboson_la_FFLAGS += -DWITH_KK
libtriboson_la_CPPFLAGS += $(GSLINCLUDE)
libtriboson_la_LIBADD = $(GSLLIBS)
endif


if WITH_TRIBOSONJET

libtriboson_la_FFLAGS += -DWITH_TRIBOSONJET

libtriboson_la_SOURCES += \
m2s_qqWAAj2.F \
m2s_qqWAAjj.F \
qqAAWj_curr.F \
qqWAAgg.F \
qqWAAqq.F 

if WITH_MADGRAPH
libtriboson_la_SOURCES += \
mg2_udWAAgg.F \
mg2_udWAAqq.F \
mg2_udgenuaa.F \
mg2_ugdenuaa.F \
mg2_dguenuaa.F 
endif

endif

$(top_builddir)/include/%.mod: libtriboson_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
