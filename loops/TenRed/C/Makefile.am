noinst_LTLIBRARIES = libTIC.la
libTIC_la_SOURCES = \
tensor_redux_Cij_new_Re_Com_1Mas.F\
tensor_redux_Cij.F \
tensor_redux_Cij_new_Re_Com.F \
tensor_redux_Cij_Div.F \
tensor_redux_CijLUde_1MasRank3.F\
tensor_redux_CijLUdeRank3.F\
tensor_redux_CijLUde_1MasRank3_Div.F\
tensor_redux_CijLUdeRank3_Div.F



if WITH_FERMIONLOOPS
libTIC_la_SOURCES += \
tensor_redux_CijLUde_1Mas_Div.F\
tensor_redux_CijLUde_Rank3_G.F\
tensor_redux_CijLUde_Rank3_G_Div.F\
tensor_redux_CijLUdeRank4_G.F\
tensor_redux_CijLUdeRank4_G_Div.F
endif

libTIC_la_FFLAGS = $(AM_FFLAGS)

if WITH_QUAD
libTIC_la_FFLAGS += -DWITH_QUAD
libTIC_la_SOURCES += \
tensor_redux_Cij_new_Re_Com_QUAD.F \
tensor_redux_CijLUdeRank4_G_QUAD.F

endif
