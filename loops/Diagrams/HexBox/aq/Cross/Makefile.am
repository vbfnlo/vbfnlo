noinst_LTLIBRARIES = libDiagramsAQCross.la

libDiagramsAQCross_la_FFLAGS = $(AM_FFLAGS)

## fix memory leak with gfortran 4.8
libDiagramsAQCross_la_FFLAGS +=$(GFORTRAN_48_FIX_LOOPS)

libDiagramsAQCross_la_SOURCES = \
CheckHexBoxCross2temp66.F\
CheckHexBoxCross2temp67.F\
CheckHexBoxCross2temp76.F\
CheckHexBoxCross2temp77.F\
HexBoxCross2temp66.F\
HexBoxCross2temp67.F\
HexBoxCross2temp76.F\
HexBoxCross2temp77.F\
HexBoxCross2temp66Div.F\
HexBoxCross2temp67Div.F\
HexBoxCross2temp76Div.F\
HexBoxCross2temp77Div.F\
HexBoxCross2tempT.F\
HexBoxCross2tempFa10.F\
HexBoxCross2tempFa11.F\
HexBoxCross2tempFa12.F\
HexBoxCross2tempFa13.F\
HexBoxCross2tempFa14.F\
HexBoxCross2tempFa15.F\
HexBoxCross2tempFa16.F\
HexBoxCross2tempFa17.F\
HexBoxCross2tempFa18.F\
HexBoxCross2tempFa19.F\
HexBoxCross2tempFa1.F\
HexBoxCross2tempFa20.F\
HexBoxCross2tempFa2.F\
HexBoxCross2tempFa3.F\
HexBoxCross2tempFa4.F\
HexBoxCross2tempFa5.F\
HexBoxCross2tempFa6.F\
HexBoxCross2tempFa7.F\
HexBoxCross2tempFa8.F\
HexBoxCross2tempFa9.F\
HexBoxCross2tempFF10.F\
HexBoxCross2tempFF11.F\
HexBoxCross2tempFF12.F\
HexBoxCross2tempFF13.F\
HexBoxCross2tempFF14.F\
HexBoxCross2tempFF15.F\
HexBoxCross2tempFF16.F\
HexBoxCross2tempFF17.F\
HexBoxCross2tempFF18.F\
HexBoxCross2tempFF19.F\
HexBoxCross2tempFF1.F\
HexBoxCross2tempFF20.F\
HexBoxCross2tempFF2.F\
HexBoxCross2tempFF3.F\
HexBoxCross2tempFF4.F\
HexBoxCross2tempFF5.F\
HexBoxCross2tempFF6.F\
HexBoxCross2tempFF7.F\
HexBoxCross2tempFF8.F\
HexBoxCross2tempFF9.F

$(top_builddir)/include/%.mod: libDiagramsAQCross_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
