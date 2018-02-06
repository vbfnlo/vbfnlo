SUBDIRS =
libIntegrals_la_LIBADD =

noinst_LTLIBRARIES = libIntegrals.la

SUBDIRS += FF
libIntegrals_la_LIBADD += FF/libFF.la
if WITH_QUAD
libIntegrals_la_LIBADD += FF/libFFquad.la
endif

libIntegrals_la_SOURCES = \
F0.F\
E0functions.F \
D0functions.F \
D0functionsDiv.F \
C0functions.F \
Aux_Functions.F\
A0.F\
B0functions.F \
cdilog.F\
dilog.F 

if WITH_FERMIONLOOPS
libIntegrals_la_SOURCES += \
B0C0D01MASG.F\
scalarints.F
endif

libIntegrals_la_FFLAGS = $(AM_FFLAGS)

if WITH_QUAD
libIntegrals_la_FFLAGS += -DWITH_QUAD
libIntegrals_la_SOURCES += \
B0QUAD.F\
C0QUAD.F\
D0QUAD.F\
D04QUAD.F\
E0QUAD.F\
F0_QUAD.F\
QUAD_poly_Re.F \
B0C0D01MASG_QUAD.F
endif

$(top_builddir)/include/%.mod: libIntegrals_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
