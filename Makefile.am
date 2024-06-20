SUBDIRS = include \
          helas \
          PDFsets \
          utilities \
          loops \
          amplitudes \
          phasespace \
          lib src doc \
          regress

ACLOCAL_AMFLAGS = -I m4
EXTRA_DIST = CHANGELOG.md GUIDELINES README.md Makefile.dependencylist

all: all-recursive

DISTCHECK_CONFIGURE_FLAGS = --enable-processes=all --enable-kk --enable-madgraph --enable-quad --with-LHAPDF
