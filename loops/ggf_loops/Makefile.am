noinst_LTLIBRARIES = libggf_loops.la
libggf_loops_la_SOURCES = \
loops_aux.F 

$(top_builddir)/include/%.mod: libggf_loops_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
