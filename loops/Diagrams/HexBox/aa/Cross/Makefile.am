noinst_LTLIBRARIES = libDiagramsAACross.la

libDiagramsAACross_la_FFLAGS = $(AM_FFLAGS)

## fix memory leak with gfortran 4.8
libDiagramsAACross_la_FFLAGS +=$(GFORTRAN_48_FIX_LOOPS)

libDiagramsAACross_la_SOURCES = \
CheckHexBoxCross1temp66.F\
CheckHexBoxCross1temp67.F\
CheckHexBoxCross1temp76.F\
CheckHexBoxCross1temp77.F\
HexBoxCross1temp66.F\
HexBoxCross1temp67.F\
HexBoxCross1temp76.F\
HexBoxCross1temp77.F\
HexBoxCross1temp66Div.F\
HexBoxCross1temp67Div.F\
HexBoxCross1temp76Div.F\
HexBoxCross1temp77Div.F\
HexBoxCross1tempT.F\
HexBoxCross1tempFa10.F\
HexBoxCross1tempFa11.F\
HexBoxCross1tempFa12.F\
HexBoxCross1tempFa13.F\
HexBoxCross1tempFa14.F\
HexBoxCross1tempFa15.F\
HexBoxCross1tempFa16.F\
HexBoxCross1tempFa17.F\
HexBoxCross1tempFa18.F\
HexBoxCross1tempFa19.F\
HexBoxCross1tempFa1.F\
HexBoxCross1tempFa20.F\
HexBoxCross1tempFa2.F\
HexBoxCross1tempFa3.F\
HexBoxCross1tempFa4.F\
HexBoxCross1tempFa5.F\
HexBoxCross1tempFa6.F\
HexBoxCross1tempFa7.F\
HexBoxCross1tempFa8.F\
HexBoxCross1tempFa9.F\
HexBoxCross1tempFF10.F\
HexBoxCross1tempFF11.F\
HexBoxCross1tempFF12.F\
HexBoxCross1tempFF13.F\
HexBoxCross1tempFF14.F\
HexBoxCross1tempFF15.F\
HexBoxCross1tempFF16.F\
HexBoxCross1tempFF17.F\
HexBoxCross1tempFF18.F\
HexBoxCross1tempFF19.F\
HexBoxCross1tempFF1.F\
HexBoxCross1tempFF20.F\
HexBoxCross1tempFF2.F\
HexBoxCross1tempFF3.F\
HexBoxCross1tempFF4.F\
HexBoxCross1tempFF5.F\
HexBoxCross1tempFF6.F\
HexBoxCross1tempFF7.F\
HexBoxCross1tempFF8.F\
HexBoxCross1tempFF9.F

$(top_builddir)/include/%.mod: libDiagramsAACross_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
