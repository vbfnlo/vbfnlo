noinst_LTLIBRARIES = libSpin2.la
libSpin2_la_SOURCES = \
aatxxx_sing.F \
aztxxx_sing.F \
wwtxxx_sing.F \
zztxxx_sing.F \
uaaxxx_sing.F \
uazxxx_sing.F \
uwwxxx_sing.F \
uzzxxx_sing.F \
aatxxx_tripn.F \
aztxxx_tripn.F \
watxxx_trippm.F \
wwtxxx_tripn.F \
wztxxx_trippm.F \
zztxxx_tripn.F \
uaaxxx_tripn.F \
uazxxx_tripn.F \
uwaxxx_trippm.F \
uwwxxx_tripn.F \
uwzxxx_trippm.F \
uzzxxx_tripn.F


libSpin2_la_FFLAGS = $(AM_FFLAGS)

$(top_builddir)/include/%.mod: libSpin2_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
