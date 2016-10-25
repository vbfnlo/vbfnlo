ALLDIRS = \
  100_Hjj   \
  100_Hjj_EW \
  100_Hjj_anom \
  101_Hjj_AA   \
  102_Hjj_mumu   \
  103_Hjj_tautau   \
  104_Hjj_bb   \
  105_Hjj_WW   \
  106_Hjj_ZZ_llll   \
  107_Hjj_ZZ_llnunu   \
  110_Hjjj   \
  111_Hjjj_AA   \
  112_Hjjj_mumu   \
  113_Hjjj_tautau   \
  114_Hjjj_bb   \
  115_Hjjj_WW   \
  116_Hjjj_ZZ_llll   \
  117_Hjjj_ZZ_llnunu   \
  120_Zjj_l \
  121_Zjj_nu \
  130_Wpjj \
  140_Wmjj \
  150_Ajj  \
  160_HHjj \
  200_WWjj \
  210_ZZjj_l \
  211_ZZjj_nu \
  220_WpZjj \
  230_WmZjj \
  250_WpWpjj \
  270_WpAjj \
  280_WmAjj \
  290_ZAjj_l \
  291_ZAjj_n \
  300_WW    \
  310_WpZ    \
  320_WmZ    \
  330_ZZ    \
  340_WpA    \
  350_WmA    \
  360_ZA    \
  370_AA    \
  400_WWZ   \
  400_WWZ_anom  \
  410_WpZZ  \
  420_WmZZ  \
  430_WpWW  \
  430_WpWW_anom  \
  430_WpWW_RS  \
  430_WpWW_MHLM  \
  440_WmWW  \
  450_ZZZ   \
  460_WWA   \
  470_ZZA   \
  480_WpZA  \
  490_WmZA  \
  500_WpAA  \
  510_WmAA  \
  520_ZAA   \
  521_ZnuAA \
  530_AAA   \
  610_WmAj  \
  620_WpAj  \
  620_WpAj_anom \
  630_WmZj  \
  640_WpZj  \
  2100_HAjj \
  2101_HAjj_AA \
  2102_HAjj_mu \
  2103_HAjj_tau \
  2104_HAjj_bbar \
  2105_HAjj_WW \
  2106_HAjj_ZZ_ll \
  2107_HAjj_ZZ_lnu \
  4100_ggf   \
  4101_ggf_AA   \
  4102_ggf_mu   \
  4103_ggf_tau   \
  4104_ggf_bb   \
  4105_ggf_WW   \
  4106_ggf_ZZ_l   \
  4107_ggf_ZZ_nu   \
  4300_GFWW   \
  4330_GFZZ    \
  4360_GFZA    \
  4370_GFAA    \
  191_AAjj_spin2 \
  195_Sp2jj_WW \
  196_Sp2jj_ZZ_llll \
  197_Sp2jj_ZZ_llnunu \
  200_WWjj_spin2 \
  210_ZZjj_l_spin2 \
  211_ZZjj_nu_spin2 \
  220_WpZjj_spin2 \
  230_WmZjj_spin2 \
  108_Hjj_WPhadWM \
  109_Hjj_WPWMhad \
  1010_Hjj_ZZhad \
  201_WPhadWMjj \
  202_WPWMhadjj \
  212_ZZhadjj \
  221_WPhadZjj \
  222_WPZhadjj \
  231_WMhadZjj \
  232_WMZhadjj \
  251_WphadWpjj \
  261_WmhadWmjj \
  301_WPhadWM \
  302_WPWMhad \
  312_WPhadZ \
  313_WPZhad \
  322_WMhadZ \
  323_WMZhad \
  331_ZZhad \
  401_WPhadWMZ \
  402_WPWMhadZ \
  403_WPWMZhad \
  411_ZZWPhad  \
  412_ZZhadWP  \
  421_ZZWMhad  \
  422_ZZhadWM  \
  431_WPhadWMWP  \
  432_WPWMhadWP  \
  441_WMWPhadWM  \
  442_WMhadWPWM  \
  451_ZZZhad \
  461_WPhadWMA \
  462_WPWMhadA \
  471_ZhadZA \
  481_WPhadZA  \
  482_WPZhadA  \
  491_WMhadZA  \
  492_WMZhadA \
  1300_WPH \
  1310_WMH \
  1330_WP \
  1340_WM \
  1600_WPHJ \
  1610_WMHJ \
  1630_WPJ \
  1640_WMJ \
  200_WWjj_TwoHiggs\
  800_WpAAj  \
  810_WmAAj  \
  3120_QCDZjj_l \
  3121_QCDZjj_nu \
  3130_QCDWPjj \
  3140_QCDWMjj \
  3210_QCDZZjj_l \
  3211_QCDZZjj_nu \
  3220_QCDWPZjj \
  3230_QCDWMZjj \
  3250_QCDWPWPjj \
  3260_QCDWMWMjj \
  3270_QCDWPAjj \
  3280_QCDWMAjj \
  3290_QCDZAjj_l \
  3291_QCDZAjj_nu \



INPUTFILES := `$(MAKE) -s -f ../src/Makefile printinputfiles`

ifeq ($(strip $(DIRS)),)
DIRS = $(ALLDIRS)
endif

all: 
	@for d in $(DIRS); do \
	   for i in $(INPUTFILES); do \
	     if ! test -r $$d/$$i; then \
	       cp ../src/$$i $$d/$$i; \
	       echo $$d/$$i >>deletedat; \
	     fi; \
	   done; \
	 done

clean:
	@rm -f `cat deletedat`
	@rm -f deletedat
	@for d in $(DIRS); do \
	   cd $$d; \
	   rm -rf check_* grid* GRID xsection* *.gp core kkcheck.dat; \
	   cd ..; \
	 done

run-all:
	@for d in $(DIRS); do \
           pwd=`pwd`; \
	   cd $$d; \
           num=`echo "$$d" | awk -F'_' '{print $$1}'`; \
           if ( [ $$num -ge 4000 ] && [ $$num -lt 5000 ] ); then \
	     qsub -cwd -b y -N check_$$d $$pwd/../src/ggflo; \
           else \
	     qsub -cwd -b y -N check_$$d $$pwd/../src/vbfnlo; \
           fi; \
	   cd ..; \
	 done

run-serial:
	@for d in $(DIRS); do \
	pwd=`pwd`; \
	   cd $$d; \
           num=`echo "$$d" | awk -F'_' '{print $$1}'`; \
           if ( [ $$num -ge 4000 ] && [ $$num -lt 5000 ] ); then \
	     $$pwd/../src/ggflo >check_$${d}.out; \
           else \
	     $$pwd/../src/vbfnlo >check_$${d}.out; \
           fi; \
	   cd ..; \
	done

printalldirs:
	@echo $(ALLDIRS)

.PHONY: all clean run-all run-serial printalldirs
