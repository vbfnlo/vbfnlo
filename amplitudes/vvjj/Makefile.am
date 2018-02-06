noinst_LTLIBRARIES = libvbf.la
libvbf_la_SOURCES = \
tensor.inc \
tensorz.inc \
tampr-lists.inc \
tensorsp2aa.inc \
tampz-lists.inc \
genps.inc \
tampr_wpp_list.inc \
gmunu.inc \
tensorw.inc \
tensorzn.inc \
toza.F \
toza_anomal.F \
qqzaqq.F\
m2s_qqwqq.F   \
qqwpzqq.F   \
qqzzqqj.F        \
m2s_qqZjqq.F  \
m2s_qqWPPjqq.F  \
qqbqq.F     \
qqwpzqqj.F  \
tozz.F \
tozz_had.F \
tozz_anomal.F \
towpm.F \
towpm_anomal.F \
towpm_had.F \
m2s_qqzqq.F   \
m2s_qqwppqq.F   \
m2s_qqHAqq.F   \
qqhqqa.F    \
qqwwqq.F    \
toww_anomal2.F \
toww_had.F \
qqhqq.F     \
qqhqq_pht.F     \
qqwwqqj.F   \
tbv.F            \
toww.F         \
qqwmzqq.F   \
qqZWjqqi.F  \
m2s_qqWjqq.F    \
qqwmzqqj.F  \
qqzzqq.F    \
HAjj_currents.F \
HAjj_ME.F \
HAjjj_ME.F \
towz.F \
towz_had.F \
towz_anomal.F \
qqwppqqj.F \
tohh.F \
qqhhqq.F \
qqhhqqj.F \
qqwppqq.F \
qqwpaqq.F \
qqwpaqqj.F \
qqwmaqq.F \
qqwmaqqj.F \
towa.F \
towa_anomal.F

libvbf_la_SOURCES += \
toww_spin2.F \
towz_spin2.F \
tozz_spin2.F \
spin2toaa.F \
vvsp2tovv.F \
qqsp2aaqq.F \
qqsp2aaqqj.F \
qqsp2wwqq.F  \
qqsp2wwqqj.F  \
qqsp2zzqq.F  \
qqsp2zzqqj.F \
toww_Hres.F \
tozz_Hres.F

if WITH_NLO
libvbf_la_SOURCES += \
vvjj_amp_aux.F \
jvirt5.F        \
gauge5-test-vv.F  \
jvirt5-isig.F \
VBF_PENT_ABE.F \
VBF_BOX_ABE.F \
VBF_BOX_ABE_CURRENT.F\
VBF_PENT_ABE_CURRENT.F
endif

libvbf_la_FFLAGS = $(AM_FFLAGS)
libvbf_la_CPPFLAGS = $(AM_CPPFLAGS)

if WITH_LOOPTOOLS
libvbf_la_FFLAGS += -DWITH_LT
endif

if WITH_MADGRAPH
libvbf_la_SOURCES += \
mg2_usus.F  \
mg2_ucuc.F    \
mg2_usdc.F    \
mg2_dsds.F \
mg2_toww.F \
mg2_towpwp.F \
mg2_towz.F \
mg2_tozz.F \
mg2_wpa_uucs.F \
mg2_wpa_ddcs.F \
mg2_wpa_udss.F \
mg2_wpa_udcc.F \
mg2_wma_uusc.F \
mg2_wma_ddsc.F \
mg2_wma_duss.F \
mg2_wma_ducc.F \
mg2_za_ucuc.F \
mg2_za_udsc.F
endif


if WITH_KK
libvbf_la_SOURCES += \
toww_kk.F \
towz_kk.F \
tozz_kk.F 
libvbf_la_FFLAGS += -DWITH_KK
libvbf_la_CPPFLAGS += $(GSLINCLUDE)
libvbf_la_LIBADD = $(GSLLIBS)
endif

if WITH_NLO
libvbf_la_FFLAGS += -DWITH_NLO
endif

if WITH_MADGRAPH
libvbf_la_FFLAGS += -DWITH_MADGRAPH
endif

$(top_builddir)/include/%.mod: libvbf_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
