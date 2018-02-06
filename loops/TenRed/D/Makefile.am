noinst_LTLIBRARIES = libTID.la
libTID_la_SOURCES = \
DijGramM.inc \
DijGram.inc \
tensor_redux_Dij_new_Re_Com_1Mas.F\
tensor_redux_Dij.F \
tensor_redux_Dij_new_Re_Com.F \
tensor_redux_Dij_Div.F \
tensor_redux_DijLUde_1MRank4.F\
tensor_redux_DijLUdeRank4123_1M.F\
tensor_redux_DijLUdeRank4123.F\
tensor_redux_DijLUdeRank4132_1M.F\
tensor_redux_DijLUdeRank4132.F\
tensor_redux_DijLUdeRank4213_1M.F\
tensor_redux_DijLUdeRank4213.F\
tensor_redux_DijLUdeRank4231_1M.F\
tensor_redux_DijLUdeRank4231.F\
tensor_redux_DijLUdeRank4312_1M.F\
tensor_redux_DijLUdeRank4312.F\
tensor_redux_DijLUdeRank4321_1M.F\
tensor_redux_DijLUdeRank4321.F\
tensor_redux_DijLUdeRank4.F\
tensor_redux_DijLUde_1MRank4_Div.F\
tensor_redux_DijLUdeRank4123_1M_Div.F\
tensor_redux_DijLUdeRank4123_Div.F\
tensor_redux_DijLUdeRank4132_1M_Div.F\
tensor_redux_DijLUdeRank4132_Div.F\
tensor_redux_DijLUdeRank4213_1M_Div.F\
tensor_redux_DijLUdeRank4213_Div.F\
tensor_redux_DijLUdeRank4231_1M_Div.F\
tensor_redux_DijLUdeRank4231_Div.F\
tensor_redux_DijLUdeRank4312_1M_Div.F\
tensor_redux_DijLUdeRank4312_Div.F\
tensor_redux_DijLUdeRank4321_1M_Div.F\
tensor_redux_DijLUdeRank4321_Div.F\
tensor_redux_DijLUdeRank4_Div.F


if WITH_FERMIONLOOPS
libTID_la_SOURCES += \
tensor_redux_DijLUdeRank4_G.F\
tensor_redux_DijLUdeRank4_G_Div.F\
tensor_redux_DijLUdeRank5_G.F\
tensor_redux_DijLUdeRank5_G_Div.F\
tensor_redux_DijLUdeRank4123_G.F\
tensor_redux_DijLUdeRank4132_G.F\
tensor_redux_DijLUdeRank4213_G.F\
tensor_redux_DijLUdeRank4231_G.F\
tensor_redux_DijLUdeRank4312_G.F\
tensor_redux_DijLUdeRank4321_G.F\
tensor_redux_DijLUdeRank4123_G_Div.F\
tensor_redux_DijLUdeRank4132_G_Div.F\
tensor_redux_DijLUdeRank4213_G_Div.F\
tensor_redux_DijLUdeRank4231_G_Div.F\
tensor_redux_DijLUdeRank4312_G_Div.F\
tensor_redux_DijLUdeRank4321_G_Div.F\
tensor_redux_DijLUde_Rank5123_G.F\
tensor_redux_DijLUde_Rank5132_G.F\
tensor_redux_DijLUde_Rank5213_G.F\
tensor_redux_DijLUde_Rank5231_G.F\
tensor_redux_DijLUde_Rank5312_G.F\
tensor_redux_DijLUde_Rank5321_G.F\
DijGram1.f\
DijGramNewM1.f
endif


libTID_la_FFLAGS = $(AM_FFLAGS)

if WITH_QUAD
libTID_la_FFLAGS += -DWITH_QUAD
libTID_la_SOURCES += \
tensor_redux_DijLUde_QUAD.F \
tensor_redux_DijLUdeRank5_G_QUAD.F \
tensor_redux_DijLUde_Rank5123_G_QUAD.F \
tensor_redux_DijLUde_Rank5132_G_QUAD.F \
tensor_redux_DijLUde_Rank5213_G_QUAD.F \
tensor_redux_DijLUde_Rank5231_G_QUAD.F \
tensor_redux_DijLUde_Rank5312_G_QUAD.F \
tensor_redux_DijLUde_Rank5321_G_QUAD.F \
XBoxGetM1M3_QUAD.F
endif

$(top_builddir)/include/%.mod: libTID_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
