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
ggf_Hpen_c.F

$(top_builddir)/include/%.mod: libggf_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
