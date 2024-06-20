noinst_LTLIBRARIES = libDiagramsAQDire.la

libDiagramsAQDire_la_FFLAGS = $(AM_FFLAGS)

## fix memory leak with gfortran 4.8
libDiagramsAQDire_la_FFLAGS +=$(GFORTRAN_48_FIX_LOOPS)

libDiagramsAQDire_la_SOURCES = \
CheckHexBox2temp66.F\
CheckHexBox2temp67.F\
CheckHexBox2temp76.F\
CheckHexBox2temp77.F\
HexBox2temp66.F\
HexBox2temp67.F\
HexBox2temp76.F\
HexBox2temp77.F\
HexBox2temp66Div.F\
HexBox2temp67Div.F\
HexBox2temp76Div.F\
HexBox2temp77Div.F\
HexBox2tempT.F\
HexBox2tempFa10.F\
HexBox2tempFa11.F\
HexBox2tempFa12.F\
HexBox2tempFa13.F\
HexBox2tempFa14.F\
HexBox2tempFa15.F\
HexBox2tempFa16.F\
HexBox2tempFa17.F\
HexBox2tempFa18.F\
HexBox2tempFa19.F\
HexBox2tempFa1.F\
HexBox2tempFa20.F\
HexBox2tempFa2.F\
HexBox2tempFa3.F\
HexBox2tempFa4.F\
HexBox2tempFa5.F\
HexBox2tempFa6.F\
HexBox2tempFa7.F\
HexBox2tempFa8.F\
HexBox2tempFa9.F\
HexBox2tempFF10.F\
HexBox2tempFF11.F\
HexBox2tempFF12.F\
HexBox2tempFF13.F\
HexBox2tempFF14.F\
HexBox2tempFF15.F\
HexBox2tempFF16.F\
HexBox2tempFF17.F\
HexBox2tempFF18.F\
HexBox2tempFF19.F\
HexBox2tempFF1.F\
HexBox2tempFF20.F\
HexBox2tempFF2.F\
HexBox2tempFF3.F\
HexBox2tempFF4.F\
HexBox2tempFF5.F\
HexBox2tempFF6.F\
HexBox2tempFF7.F\
HexBox2tempFF8.F\
HexBox2tempFF9.F

$(top_builddir)/include/%.mod: libDiagramsAQDire_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
