noinst_LTLIBRARIES = libDiagramsFLoops.la

libDiagramsFLoops_la_FFLAGS = $(AM_FFLAGS)

## fix memory leak with gfortran 4.8
libDiagramsFLoops_la_FFLAGS +=$(GFORTRAN_48_FIX_LOOPS)

libDiagramsFLoops_la_SOURCES = \
FermionBoxEventemp1.F\
FermionVertEventemp.F\
FermionBoxEventemp1_Div.F\
FermionVertEventemp_Div.F\
FermionVertCouplingT.F\
FermionVertCouplingT_Div.F\
FermionVertEventempCoupling.F\
FermionVertEventempCoupling_Div.F\
FermionVertOddtempCoupling.F\
FermionVertOddtempCoupling_Div.F\
FermionBoxCouplingT.F\
FermionBoxCouplingT_Div.F\
FermionBoxCouplingT_GM.F\
FermionBoxCouplingT_GM_Div.F\
FermionBoxEventempCoupling.F\
FermionBoxEventempCoupling_Div.F\
FermionBoxEvenG1tempCoupling1.F\
FermionBoxEvenG1tempCoupling1_Div.F\
FermionBoxOddtempCoupling.F\
FermionBoxOddtempCoupling_Div.F\
FermionBoxOddtempCoupling_G.F\
FermionBoxOddtempCoupling_G_Div.F\
HiggsEventempver1.F\
HiggsEventempver_Div.F


if WITH_QUAD

libDiagramsFLoops_la_FFLAGS += -DWITH_QUAD
libDiagramsFLoops_la_SOURCES += \
XBoxGetM1M3_QUAD.F\
FermionBoxEventempCoupling_QUAD.F

endif
