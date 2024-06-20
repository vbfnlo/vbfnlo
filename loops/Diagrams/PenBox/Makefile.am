noinst_LTLIBRARIES = libDiagramsPenBox.la

libDiagramsPenBox_la_FFLAGS = $(AM_FFLAGS)

## fix memory leak with gfortran 4.8
libDiagramsPenBox_la_FFLAGS +=$(GFORTRAN_48_FIX_LOOPS)

libDiagramsPenBox_la_SOURCES = \
Dir/PenBoxtempT.F\
Dir/PenBoxtemp66Div.F\
Dir/PenBoxtemp67Div.F\
Dir/PenBoxtemp76Div.F\
Dir/PenBoxtemp77Div.F\
Dir/PenBoxtempFF1.F\
Dir/PenBoxtemp66.F\
Dir/PenBoxtemp67.F\
Dir/PenBoxtemp76.F\
Dir/PenBoxtemp77.F\
Dir/PenBoxtempFF2.F\
Cross/PenBoxCrosstempT.F\
Cross/PenBoxCrosstemp66Div.F\
Cross/PenBoxCrosstemp67Div.F\
Cross/PenBoxCrosstemp76Div.F\
Cross/PenBoxCrosstemp77Div.F\
Cross/PenBoxCrosstempFF1.F\
Cross/PenBoxCrosstemp66.F\
Cross/PenBoxCrosstemp67.F\
Cross/PenBoxCrosstemp76.F\
Cross/PenBoxCrosstemp77.F\
Cross/PenBoxCrosstempFF2.F


$(top_builddir)/include/%.mod: libDiagramsPenBox_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
