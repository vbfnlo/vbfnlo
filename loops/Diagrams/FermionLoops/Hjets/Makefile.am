noinst_LTLIBRARIES = libDiagramsFLoopsHjets.la

libDiagramsFLoopsHjets_la_FFLAGS = $(AM_FFLAGS)

## fix memory leak with gfortran 4.8
libDiagramsFLoopsHjets_la_FFLAGS +=$(GFORTRAN_48_FIX_LOOPS)

libDiagramsFLoopsHjets_la_SOURCES = \
HiggsEventempver1.F\
HiggsEventempver_Div.F\
HexHiggsEventemp.F\
HexHiggsEventempFa.F\
HexHiggsEventempFF1.F\
HexHiggsEventempFF2.F\
HexHiggsEventempFF3.F\
HexHiggsEventempFF4.F\
HexHiggsEventempFF5.F\
HexHiggsOddtemp.F\
HexHiggsOddtempFa.F\
HexHiggsOddtempFa_optimize.F\
HexHiggsOddtemp_optimize.F\
HexHiggsOddtempFF1.F\
HexHiggsOddtempFF2.F\
HexHiggsOddtempFF3.F\
HexHiggsOddtempFF4.F\
HexHiggsOddtempFF5.F\
HiggsEventempbox.F\
HiggsEventempboxT.F\
HiggsEventemp.F\
HiggsEventempver.F\
HiggsOddtempbox.F\
HiggsOddtempbox_old.F\
HiggsOddtempboxT.F\
HiggsOddtemppen.F\
HiggsOddtemppen_old.F\
HiggsOddtempvert.F \
HiggsEventemppenT.F \
HiggsEventemppenT_QUAD.F \
HiggsOddtemppenT.F \
HiggsOddtemppenT_QUAD.F \
QUAD_HiggsOddtempbox.F \
QUAD_HiggsOddtemppen.F \
QUAD_HiggsEventemp.F \
QUAD_HiggsEventempbox.F


if WITH_QUAD

libDiagramsFLoopsHjets_la_FFLAGS += -DWITH_QUAD
libDiagramsFLoopsHjets_la_SOURCES += \
QUAD_HexHiggsEventemp.F\
QUAD_HexHiggsOddtemp.F\
QUAD_HexHiggsOddtemp_optimize.F\
QUAD_HiggsBoxEvenT.F\
QUAD_HiggsBoxOddT.F\
QUAD_HiggsEventempbox.F\
QUAD_HiggsOddtempbox.F\
QUAD_HiggsOddtempbox_old.F\
QUAD_HiggsEventemp.F\
QUAD_HiggsEventempver.F\
QUAD_HiggsOddtemppen.F\
QUAD_HiggsOddtemppen_old.F\
QUAD_HiggsOddtempvert.F
endif
