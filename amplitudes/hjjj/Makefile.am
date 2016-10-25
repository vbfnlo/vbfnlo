noinst_LTLIBRARIES = libhjjj.la
libhjjj_la_SOURCES = \
dc_usbbh.F \
m2s_qqh3j.F  \
qqhqq.F       \
qqhqqgg_mgfast.F \
qqh4q_sub.F  \
m2s_qqh4j.F    \
qqh4q_mg.F   \
qqhqqgg_cc.F  \
qqh4qcc_sub.F  \
qqhqqgg_mg.F  \
us_usbbh.F 

if WITH_NLO
libhjjj_la_SOURCES += \
hjjj_amp_aux.F\
coeff.F        \
dipsub.F     \
dipolesub.F  
endif

libhjjj_la_FFLAGS = $(AM_FFLAGS)

if WITH_NLO
libhjjj_la_FFLAGS += -DWITH_NLO
endif



