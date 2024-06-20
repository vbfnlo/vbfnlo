noinst_LTLIBRARIES = libXGet.la
libXGet_la_SOURCES = \
XEgetM1M3.F\
XEgetM1M3Sim.F\
XEgetM1M3_Div.F\
XEgetOddM1M3.F\
XEgetOddSimM1M3.F\
XEgetOddM1M3_Div.F\
MyCget.F\
MyCget_M.F

libXGet_la_FFLAGS = $(AM_FFLAGS)

if WITH_FERMIONLOOPS
libXGet_la_SOURCES += \
MyCgetGram.F\
MyCgetGram_M.F
endif


if WITH_QUAD
libXGet_la_FFLAGS += -DWITH_QUAD
libXGet_la_SOURCES += \
XBoxGetM1M3_QUAD_check.F\
XBoxGetM1M3_QUAD.F \
XDGetM0_QUAD_Dble.F \
XEGetM0_QUAD_Dble2g.F \
XEGetM0_QUAD_Dble_4g.F \
XEGetM0_QUAD_Dble.F \
XEgetM1M3_QUAD.F\
XEgetM1M3_QUAD_noSim.F \
XEgetOddM1M3_QUAD.F \
XEgetOddSimM1M3_QUAD.F \
XFGetM0_QUAD_Dble2g.F \
XFGetM0_QUAD_Dble_4g.F \
XFGetM0_QUAD_Dble.F \
XGetC_QUAD_2Dble.F \
XGetCOdd_QUAD_2Dble.F \
XGetD_QUAD_2Dble.F\
XGetE_QUAD_2Dble.F\
XGetF_QUAD_2Dble.F
endif
