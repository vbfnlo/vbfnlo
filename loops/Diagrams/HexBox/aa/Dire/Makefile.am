noinst_LTLIBRARIES = libDiagramsAADire.la

libDiagramsAADire_la_FFLAGS = $(AM_FFLAGS)

## fix memory leak with gfortran 4.8
libDiagramsAADire_la_FFLAGS +=$(GFORTRAN_48_FIX_LOOPS)

libDiagramsAADire_la_SOURCES = \
CheckHexBoxtemp66.F\
CheckHexBoxtemp67.F\
CheckHexBoxtemp76.F\
CheckHexBoxtemp77.F\
HexBoxtemp66.F\
HexBoxtemp67.F\
HexBoxtemp76.F\
HexBoxtemp77.F\
HexBoxtemp66Div.F\
HexBoxtemp67Div.F\
HexBoxtemp76Div.F\
HexBoxtemp77Div.F\
HexBoxtempT.F\
HexBoxtempFa10.F\
HexBoxtempFa11.F\
HexBoxtempFa12.F\
HexBoxtempFa13.F\
HexBoxtempFa14.F\
HexBoxtempFa15.F\
HexBoxtempFa16.F\
HexBoxtempFa17.F\
HexBoxtempFa18.F\
HexBoxtempFa19.F\
HexBoxtempFa1.F\
HexBoxtempFa20.F\
HexBoxtempFa2.F\
HexBoxtempFa3.F\
HexBoxtempFa4.F\
HexBoxtempFa5.F\
HexBoxtempFa6.F\
HexBoxtempFa7.F\
HexBoxtempFa8.F\
HexBoxtempFa9.F\
HexBoxtempFF10.F\
HexBoxtempFF11.F\
HexBoxtempFF12.F\
HexBoxtempFF13.F\
HexBoxtempFF14.F\
HexBoxtempFF15.F\
HexBoxtempFF16.F\
HexBoxtempFF17.F\
HexBoxtempFF18.F\
HexBoxtempFF19.F\
HexBoxtempFF1.F\
HexBoxtempFF20.F\
HexBoxtempFF2.F\
HexBoxtempFF3.F\
HexBoxtempFF4.F\
HexBoxtempFF5.F\
HexBoxtempFF6.F\
HexBoxtempFF7.F\
HexBoxtempFF8.F\
HexBoxtempFF9.F

$(top_builddir)/include/%.mod: libDiagramsAADire_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
