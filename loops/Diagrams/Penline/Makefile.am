noinst_LTLIBRARIES = libDiagramsPenL.la

libDiagramsPenL_la_FFLAGS = $(AM_FFLAGS) 

## fix memory leak with gfortran 4.8
libDiagramsPenL_la_FFLAGS +=$(GFORTRAN_48_FIX_LOOPS)

libDiagramsPenL_la_SOURCES = \
penline.F\
penlineABE.F\
penlineABEDiv.F\
penlineABETotal.F\
penlineNoAbe.F\
penlineNoAbeDiv.F\
penlineNoAbeTotal.F\
penlineNoAbe_parts.F\
penlineNoAbeDiv_parts.F\
penlineNoAbeTotal_parts.F\
tri_box_penta_tens_red.F \
virtual1.F \
virtual2.F 

if WITH_QCD_JJ_PLUSX
libDiagramsPenL_la_SOURCES += \
NoAbe2gaPenLinetempFF10.F\
NoAbe2gaPenLinetempFa1.F\
NoAbe2gaPenLinetempFa2.F\
NoAbe2gaPenLinetempFa3.F\
NoAbe2gaPenLinetempFF1.F\
NoAbe2gaPenLinetempFF2.F\
NoAbe2gaPenLinetempFF3.F\
NoAbe2gaPenLinetemp.F	\
NoAbe2gaPenLinetempFF4.F\
NoAbe2gaPenLinetempFF5.F\
NoAbe2gaPenLinetempFF6.F\
NoAbe2gaPenLinetempDiv.F\
NoAbe2gaPenLinetempFF7.F\
NoAbe2gaPenLinetempFF8.F\
NoAbe2gaPenLinetempFF9.F\
penlineNoAbe2gaTotal.F\
NoAbe4vertBoxLinetemp.F\
NoAbe4vertBoxLinetempFF1.F\
NoAbe4vertBoxLinetempFa1.F\
NoAbe4vertBoxLinetempDiv.F\
Boxline4Vert_Total.F
endif

if WITH_QUAD
if WITH_QCD_JJ_PLUSX
libDiagramsPenL_la_FFLAGS += -DWITH_QUAD
libDiagramsPenL_la_SOURCES += \
NoAbe4vertBoxLinetemp_QUAD.F\
XEGetM0_QUAD_Dble.F \
XEGetM0_QUAD_Dble2g.F \
XEGetM0_QUAD_Dble_4g.F \
penlineNoAbe_parts_QUAD.F  \
penlineNoAbe_QUAD.F  \
NoAbe2gaPenLinetemp_QUAD.F	\
penline_QUAD.F  \
penlineABE_QUAD.F  
endif
endif

$(top_builddir)/include/%.mod: libDiagramsPenL_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
