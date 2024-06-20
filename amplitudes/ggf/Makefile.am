
noinst_LTLIBRARIES = libggf.la

libggf_la_FFLAGS = $(AM_FFLAGS) 

## fix memory leak with gfortran 4.8
libggf_la_FFLAGS +=$(GFORTRAN_48_FIX_GGF)

libggf_la_SOURCES = \
ggf_params.inc \
ggf_dotproducts.inc \
ggf_process.inc \
d_coefficients.inc \
d_functions.inc \
c_coefficients.inc \
c_functions.inc \
ggf_default.inc \
gf_higgs_ME.F    \
ggf_dotproducts.F    \
c_functions.F       \
d_functions.F       \
ggf_amp_aux.F    \
m2s_ggfh.F \
ggf_tri_box.F \
ggf_Apen.F \
ggf_Hpen.F \
ggf_Hpen_a.F \
ggf_Hpen_b.F \
ggf_Hpen_c.F \
m2s_ggf3jh.F \
gf_qqhqqg.F \
gf_qqhqqg-bottom.F \
gf_qqaqqg-eff.F \
gf_qqaqqg-top.F \
gf_qqaqqg-bottom.F \
gf_qghqgg.F \
gf_qghqgg-bottom.F \
gf_qgaqgg-eff.F \
gf_qgaqgg-top.F \
gf_qgaqgg-bottom.F \
gf_gghggg.F \
mgcheck.F \
ggfloMG.inc \
genps.inc \
dimension.inc	\
nexternal.inc \
gf_ucucgH.F \
gf_uuuugH.F \
gf_uuxccxgH.F \
gf_uuxuuxgH.F \
gf_guuccxH.F \
gf_guuuuxH.F \
matrixuguggH.F \
matrixgggggH.F \
matrixgggggA.F \
matrixuuuugH.F \
matrixuuuugH_MG5.F \
matrixucucgH.F \
matrixucucgH_MG5.F \
matrixguuccxH.F \
matrixguhccu_MG5.F \
matrixguuuuxH.F \
matrixguuuh_MG5.F \
matrixuuxuuxgH.F \
matrixubuhubug_MG5.F \
matrixuuxccxgH.F \
matrixubuhdbdg_MG5.F \
ggf_tri_box_pen_hex.F 

$(top_builddir)/include/%.mod: libggf_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
