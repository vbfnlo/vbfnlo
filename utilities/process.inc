
!!! THIS .INC FILE IS USED IN FORTRAN 77 AND FORTRAN 90 FILES: 
!!! CONTINUATION LINES ARE NOT ALLOWED, COMMENTS HAVE TO START WITH AN !


      INTEGER lepgens
      LOGICAL alllep_flag
      common /alllep/ lepgens, alllep_flag

      INTEGER process
      common /proc/ process

      DOUBLE PRECISION symmetryfactor
      common /factor/ symmetryfactor

      INTEGER narrowwidth
      common /narrow/ narrowwidth

      LOGICAL LOplusjet
      common /proc2/ LOplusjet

      logical photoncut_vbf, zdecay_had
      common /techcuts/ photoncut_vbf, zdecay_had

      INTEGER finallep, finallep_default, N_lepton
      common /leptonnumbers/ finallep(6), finallep_default(6), N_lepton

      INTEGER finalquarks(2), finalquarks_psp(2), finalquarks_default(2), N_decayquarks
      common /quarknumbers/ finalquarks, finalquarks_psp, finalquarks_default, N_decayquarks

      INTEGER Hjj, Hjj_AA, Hjj_mu, Hjj_tau, Hjj_bbar
      INTEGER Hjj_WW, Hjj_ZZ_ll, Hjj_ZZ_lnu
      INTEGER HjjjLO, HjjjLO_AA, HjjjLO_mu, HjjjLO_tau, HjjjLO_bbar
      INTEGER HjjjLO_WW, HjjjLO_ZZ_ll, HjjjLO_ZZ_lnu
      INTEGER Hjj_WPhadWM, Hjj_WPWMhad, Hjj_ZZhad
      INTEGER HjjjLO_WPhadWM, HjjjLO_WPWMhad, HjjjLO_ZZhad
      PARAMETER( Hjj           = 100  )
      PARAMETER( Hjj_AA        = 101  )
      PARAMETER( Hjj_mu        = 102  )
      PARAMETER( Hjj_tau       = 103  )
      PARAMETER( Hjj_bbar      = 104  )
      PARAMETER( Hjj_WW        = 105  )
      PARAMETER( Hjj_ZZ_ll     = 106  )
      PARAMETER( Hjj_ZZ_lnu    = 107  )
      PARAMETER( Hjj_WPhadWM   = 108  )
      PARAMETER( Hjj_WPWMhad   = 109  )
      PARAMETER( Hjj_ZZhad     = 1010 )
      PARAMETER( HjjjLO        = 9100 )
      PARAMETER( HjjjLO_AA     = 9101 )
      PARAMETER( HjjjLO_mu     = 9102 )
      PARAMETER( HjjjLO_tau    = 9103 )
      PARAMETER( HjjjLO_bbar   = 9104 )
      PARAMETER( HjjjLO_WW     = 9105 )
      PARAMETER( HjjjLO_ZZ_ll  = 9106 )
      PARAMETER( HjjjLO_ZZ_lnu = 9107 )
      PARAMETER( HjjjLO_WPhadWM = 9108  )
      PARAMETER( HjjjLO_WPWMhad = 9109  )
      PARAMETER( HjjjLO_ZZhad   = 91010 )

      INTEGER HJJJ,HJJJJLO
      INTEGER Hjjj_AA, Hjjj_mu, Hjjj_tau, Hjjj_bbar
      INTEGER Hjjj_WW, Hjjj_ZZ_ll, Hjjj_ZZ_lnu
      INTEGER HjjjjLO_AA, HjjjjLO_mu, HjjjjLO_tau, HjjjjLO_bbar
      INTEGER HjjjjLO_WW, HjjjjLO_ZZ_ll, HjjjjLO_ZZ_lnu
      PARAMETER( Hjjj           = 110 )
      PARAMETER( Hjjj_AA        = 111  )
      PARAMETER( Hjjj_mu        = 112  )
      PARAMETER( Hjjj_tau       = 113  )
      PARAMETER( Hjjj_bbar      = 114  )
      PARAMETER( Hjjj_WW        = 115  )
      PARAMETER( Hjjj_ZZ_ll     = 116  )
      PARAMETER( Hjjj_ZZ_lnu    = 117  )
      PARAMETER( HjjjjLO        = 9110 )
      PARAMETER( HjjjjLO_AA     = 9111  )
      PARAMETER( HjjjjLO_mu     = 9112  )
      PARAMETER( HjjjjLO_tau    = 9113  )
      PARAMETER( HjjjjLO_bbar   = 9114  )
      PARAMETER( HjjjjLO_WW     = 9115  )
      PARAMETER( HjjjjLO_ZZ_ll  = 9116  )
      PARAMETER( HjjjjLO_ZZ_lnu = 9117  )

      INTEGER Zjj_l, Zjj_nu, WPjj, WMjj
      INTEGER ZjjjLO_l, ZjjjLO_nu, WPjjjLO, WMjjjLO
      INTEGER Ajj
      INTEGER AjjjLO
      PARAMETER( Zjj_l         = 120  )
      PARAMETER( Zjj_nu        = 121  )
      PARAMETER( WPjj          = 130  )
      PARAMETER( WMjj          = 140  )      
      PARAMETER( Ajj           = 150  )
      PARAMETER( ZjjjLO_l      = 9120 )
      PARAMETER( ZjjjLO_nu     = 9121 )
      PARAMETER( WPjjjLO       = 9130 )
      PARAMETER( WMjjjLO       = 9140 )
      PARAMETER( AjjjLO        = 9150 )

      INTEGER AAjj, Sp2jj_WW, AAjjjLO, Sp2jjjLO_WW
      INTEGER Sp2jj_ZZ_ll, Sp2jjjLO_ZZ_ll, Sp2jj_ZZ_lnu, Sp2jjjLO_ZZ_lnu
      PARAMETER( AAjj           = 191  )
      PARAMETER( Sp2jj_WW       = 195  )
      PARAMETER( Sp2jj_ZZ_ll    = 196  )
      PARAMETER( Sp2jj_ZZ_lnu   = 197  )
      PARAMETER( AAjjjLO        = 9191 )
      PARAMETER( Sp2jjjLO_WW    = 9195 )
      PARAMETER( Sp2jjjLO_ZZ_ll = 9196 )
      PARAMETER( Sp2jjjLO_ZZ_lnu= 9197 )

      INTEGER HHjj, HHjj_bbtau, HHjj_bbAA
      INTEGER HHjjjLO, HHjjjLO_bbtau, HHjjjLO_bbAA
      PARAMETER( HHjj          = 160  )
      PARAMETER( HHjj_bbtau    = 161  )
      PARAMETER( HHjj_bbAA     = 162  )
      PARAMETER( HHjjjLO       = 9160 )
      PARAMETER( HHjjjLO_bbtau = 9161 )
      PARAMETER( HHjjjLO_bbAA  = 9162 )

      INTEGER WPWMjj, WPWPjj, WMWMjj, WPZjj, WMZjj, WMAjj, WPAjj, ZAjj, ZAjj_n
      INTEGER WPhadWMjj, WPWMhadjj, WPhadZjj, WMhadZjj, WPZhadjj, WMZhadjj, WPhadWPjj, WMhadWMjj
      INTEGER ZZjj_ll, ZZjj_lnu, ZZhadjj
      INTEGER WPWMjjjLO,WPWPjjjLO,WMWMjjjLO,WPZjjjLO,WMZjjjLO,WPAjjjLO,WMAjjjLO,ZAjjjLO,ZAjjjLO_n
      INTEGER WPhadWMjjjLO, WPWMhadjjjLO, WPhadZjjjLO, WMhadZjjjLO, WPZhadjjjLO, WMZhadjjjLO,WPhadWPjjjLO,WMhadWMjjjLO
      INTEGER ZZjjjLO_ll,ZZjjjLO_lnu, ZZhadjjjLO
      PARAMETER( WPWMjj        = 200  )
      PARAMETER( WPhadWMjj     = 201  )
      PARAMETER( WPWMhadjj     = 202  )
      PARAMETER( ZZjj_ll       = 210  )
      PARAMETER( ZZjj_lnu      = 211  )
      PARAMETER( ZZhadjj       = 212  )
      PARAMETER( WPZjj         = 220  )
      PARAMETER( WPhadZjj      = 221  )
      PARAMETER( WPZhadjj      = 222  )
      PARAMETER( WMZjj         = 230  )
      PARAMETER( WMhadZjj      = 231  )
      PARAMETER( WMZhadjj      = 232  )
      PARAMETER( WPWPjj        = 250  )
      PARAMETER( WPhadWPjj     = 251  )
      PARAMETER( WMWMjj        = 260  )
      PARAMETER( WMhadWMjj     = 261  )
      PARAMETER( WPAjj         = 270  )
      PARAMETER( WMAjj         = 280  )
      PARAMETER( ZAjj          = 290  )
      PARAMETER( ZAjj_n        = 291  )
      PARAMETER( WPWMjjjLO     = 9200 )
      PARAMETER( WPhadWMjjjLO  = 9201 )
      PARAMETER( WPWMhadjjjLO  = 9202 )
      PARAMETER( ZZjjjLO_ll    = 9210 )
      PARAMETER( ZZjjjLO_lnu   = 9211 )
      PARAMETER( ZZhadjjjLO    = 9212 )
      PARAMETER( WPZjjjLO      = 9220 )
      PARAMETER( WPhadZjjjLO   = 9221 )
      PARAMETER( WPZhadjjjLO   = 9222 )
      PARAMETER( WMZjjjLO      = 9230 )
      PARAMETER( WMhadZjjjLO   = 9231 )
      PARAMETER( WMZhadjjjLO   = 9232 )
      PARAMETER( WPWPjjjLO     = 9250  )
      PARAMETER( WPhadWPjjjLO  = 9251  )
      PARAMETER( WMWMjjjLO     = 9260  )
      PARAMETER( WMhadWMjjjLO  = 9261  )
      PARAMETER( WPAjjjLO      = 9270 )
      PARAMETER( WMAjjjLO      = 9280 )
      PARAMETER( ZAjjjLO       = 9290 )
      PARAMETER( ZAjjjLO_n     = 9291 )


      INTEGER WW, WPZ, WMZ, ZZ, WPA, WMA, ZA, AA
      INTEGER WWjLO, WPZjLOb, WMZjLOb, ZZjLO, WPAjLOb, WMAjLOb, ZAjLO, AAjLO
      INTEGER WPhadWMlep, WPlepWMhad, WPhadWMlepjLO, WPlepWMhadjLO
      INTEGER ZZhad, ZZhadjLO
      INTEGER WPhadZ, WPZhad, WMhadZ, WMZhad, WPhadZjLO, WPZhadjLO, WMhadZjLO, WMZhadjLO
      PARAMETER( WW            = 300  )
      PARAMETER( WPhadWMlep    = 301  )
      PARAMETER( WPlepWMhad    = 302  )
      PARAMETER( WPZ           = 310  )
      PARAMETER( WPhadZ        = 312  )
      PARAMETER( WPZhad        = 313  )
      PARAMETER( WMZ           = 320  )
      PARAMETER( WMhadZ        = 322  )
      PARAMETER( WMZhad        = 323  )
      PARAMETER( ZZ            = 330  )
      PARAMETER( ZZhad         = 331  )
      PARAMETER( WPA           = 340  )
      PARAMETER( WMA           = 350  )
      PARAMETER( ZA            = 360  )
      PARAMETER( AA            = 370  )
      PARAMETER( WWjLO         = 9300 )
      PARAMETER( WPhadWMlepjLO = 9301 )
      PARAMETER( WPlepWMhadjLO = 9302 )
      PARAMETER( WPZjLOb       = 9310 )
      PARAMETER( WPhadZjLO     = 9312 )
      PARAMETER( WPZhadjLO     = 9313 )
      PARAMETER( WMZjLOb       = 9320 )
      PARAMETER( WMhadZjLO     = 9322 )
      PARAMETER( WMZhadjLO     = 9323 )
      PARAMETER( ZZjLO         = 9330 )
      PARAMETER( ZZhadjLO      = 9331 )
      PARAMETER( WPAjLOb       = 9340 )
      PARAMETER( WMAjLOb       = 9350 )       
      PARAMETER( ZAjLO         = 9360 )
      PARAMETER( AAjLO         = 9370 )

      INTEGER WPH, WPH_AA, WPH_mu, WPH_tau, WPH_bbar, WPH_WW, WPH_ZZ_ll, WPH_ZZ_lnu
      INTEGER WMH, WMH_AA, WMH_mu, WMH_tau, WMH_bbar, WMH_WW, WMH_ZZ_ll, WMH_ZZ_lnu
      INTEGER WPHjLOb, WPHjLOb_AA, WPHjLOb_mu, WPHjLOb_tau, WPHjLOb_bbar 
      INTEGER WPHjLOb_WW, WPHjLOb_ZZ_ll, WPHjLOb_ZZ_lnu
      INTEGER WMHjLOb, WMHjLOb_AA, WMHjLOb_mu, WMHjLOb_tau, WMHjLOb_bbar 
      INTEGER WMHjLOb_WW, WMHjLOb_ZZ_ll, WMHjLOb_ZZ_lnu

      PARAMETER( WPH           = 1300  )
      PARAMETER( WPH_AA        = 1301  )
      PARAMETER( WPH_mu        = 1302  )
      PARAMETER( WPH_tau       = 1303  )
      PARAMETER( WPH_bbar      = 1304  )
      PARAMETER( WPH_WW        = 1305  )
      PARAMETER( WPH_ZZ_ll     = 1306  )
      PARAMETER( WPH_ZZ_lnu    = 1307  )

      PARAMETER( WPHjLOb           = 91300  )
      PARAMETER( WPHjLOb_AA        = 91301  )
      PARAMETER( WPHjLOb_mu        = 91302  )
      PARAMETER( WPHjLOb_tau       = 91303  )
      PARAMETER( WPHjLOb_bbar      = 91304  )
      PARAMETER( WPHjLOb_WW        = 91305  )
      PARAMETER( WPHjLOb_ZZ_ll     = 91306  )
      PARAMETER( WPHjLOb_ZZ_lnu    = 91307  )

      PARAMETER( WMH           = 1310  )
      PARAMETER( WMH_AA        = 1311  )
      PARAMETER( WMH_mu        = 1312  )
      PARAMETER( WMH_tau       = 1313  )
      PARAMETER( WMH_bbar      = 1314  )
      PARAMETER( WMH_WW        = 1315  )
      PARAMETER( WMH_ZZ_ll     = 1316  )
      PARAMETER( WMH_ZZ_lnu    = 1317  )

      PARAMETER( WMHjLOb           = 91310  )
      PARAMETER( WMHjLOb_AA        = 91311  )
      PARAMETER( WMHjLOb_mu        = 91312  )
      PARAMETER( WMHjLOb_tau       = 91313  )
      PARAMETER( WMHjLOb_bbar      = 91314  )
      PARAMETER( WMHjLOb_WW        = 91315  )
      PARAMETER( WMHjLOb_ZZ_ll     = 91316  )
      PARAMETER( WMHjLOb_ZZ_lnu    = 91317  )

      INTEGER WP_only, WM_only, WPjLOb, WMjLOb
      PARAMETER( WP_only      = 1330  )
      PARAMETER( WM_only      = 1340  )
      PARAMETER( WPjLOb       = 91330  )
      PARAMETER( WMjLOb       = 91340  )


      INTEGER GFWW, GFZZ, GFZA, GFAA
      INTEGER GFWPhadWMlep, GFWPlepWMhad
      INTEGER GFZZhad
      PARAMETER( GFWW          = 4300 )
      PARAMETER( GFWPhadWMlep  = 4301 )
      PARAMETER( GFWPlepWMhad  = 4302 )
      PARAMETER( GFZZ          = 4330 )
      PARAMETER( GFZZhad       = 4331 )
      PARAMETER( GFZA          = 4360 )
      PARAMETER( GFAA          = 4370 )


      INTEGER WWZ, ZZWP, ZZWM, WWWP, WWWM, WWA, ZZA, ZZnA, WPZA, WMZA
      INTEGER WPAA, WMAA, ZAA,ZnAA, AAA
      INTEGER WWZjLO, ZZWPjLO, ZZWMjLO, WWWPjLO, WWWMjLO, WWAjLO, ZZAjLO, ZZnAjLO
      INTEGER WPZAjLO, WMZAjLO, WPAAjLO, WMAAjLO, ZAAjLO,ZnAAjLO, AAAjLO
      INTEGER WPAAj, WMAAj
      INTEGER WPAAjjLO, WMAAjjLO
      INTEGER ZZZ
      INTEGER ZZZjLO
      INTEGER WPhadWMZ, WPWMhadZ, WWZhad, ZZhadWP, ZZWPhad, ZZhadWM, ZZWMhad
      INTEGER WPhadWMWP, WPWMhadWP, WMhadWPWM, WMWPhadWM, ZZZhad
      INTEGER WPhadWMZjLO, WPWMhadZjLO, WWZhadjLO, ZZhadWPjLO, ZZWPhadjLO, ZZhadWMjLO, ZZWMhadjLO
      INTEGER WPhadWMWPjLO, WPWMhadWPjLO, WMhadWPWMjLO, WMWPhadWMjLO, ZZZhadjLO
      INTEGER WPhadWMA, WPWMhadA, ZZhadA, WPhadZA, WPZhadA, WMhadZA, WMZhadA
      INTEGER WPhadWMAjLO, WPWMhadAjLO, ZZhadAjLO, WPhadZAjLO, WPZhadAjLO, WMhadZAjLO, WMZhadAjLO

      PARAMETER( WWZ           = 400  )
      PARAMETER( WPhadWMZ      = 401  )
      PARAMETER( WPWMhadZ      = 402  )
      PARAMETER( WWZhad        = 403  )
      PARAMETER( ZZWP          = 410  )
      PARAMETER( ZZWPhad       = 411  )
      PARAMETER( ZZhadWP       = 412  )
      PARAMETER( ZZWM          = 420  )
      PARAMETER( ZZWMhad       = 421  )
      PARAMETER( ZZhadWM       = 422  )
      PARAMETER( WWWP          = 430  )
      PARAMETER( WPhadWMWP     = 431  )
      PARAMETER( WPWMhadWP     = 432  )
      PARAMETER( WWWM          = 440  )
      PARAMETER( WMWPhadWM     = 441  )
      PARAMETER( WMhadWPWM     = 442  )
      PARAMETER( ZZZ           = 450  )
      PARAMETER( ZZZhad        = 451  )

      PARAMETER( WWA           = 460  )
      PARAMETER( WPhadWMA      = 461  )
      PARAMETER( WPWMhadA      = 462  )
      PARAMETER( ZZA           = 470  )
      PARAMETER( ZZhadA        = 471  )
      PARAMETER( ZZnA          = 472  )
      PARAMETER( WPZA          = 480  )
      PARAMETER( WPhadZA       = 481  )
      PARAMETER( WPZhadA       = 482  )
      PARAMETER( WMZA          = 490  )
      PARAMETER( WMhadZA       = 491  )
      PARAMETER( WMZhadA       = 492  )
      PARAMETER( WPAA          = 500  )
      PARAMETER( WMAA          = 510  )
      PARAMETER( ZAA           = 520  )
      PARAMETER( ZnAA          = 521  )
      PARAMETER( AAA           = 530  )

      PARAMETER( WWZjLO        = 9400 )
      PARAMETER( WPhadWMZjLO   = 9401 )
      PARAMETER( WPWMhadZjLO   = 9402 )
      PARAMETER( WWZhadjLO     = 9403 )
      PARAMETER( ZZWPjLO       = 9410 )
      PARAMETER( ZZWPhadjLO    = 9411 )
      PARAMETER( ZZhadWPjLO    = 9412 )
      PARAMETER( ZZWMjLO       = 9420 )
      PARAMETER( ZZWMhadjLO    = 9421 )
      PARAMETER( ZZhadWMjLO    = 9422 )
      PARAMETER( WWWPjLO       = 9430 )
      PARAMETER( WPhadWMWPjLO  = 9431 )
      PARAMETER( WPWMhadWPjLO  = 9432 )
      PARAMETER( WWWMjLO       = 9440 )
      PARAMETER( WMWPhadWMjLO  = 9441 )
      PARAMETER( WMhadWPWMjLO  = 9442 )
      PARAMETER( ZZZjLO        = 9450 )
      PARAMETER( ZZZhadjLO     = 9451 )

      PARAMETER( WWAjLO        = 9460 )
      PARAMETER( WPhadWMAjLO   = 9461 )
      PARAMETER( WPWMhadAjLO   = 9462 )
      PARAMETER( ZZAjLO        = 9470 )
      PARAMETER( ZZhadAjLO     = 9471 )
      PARAMETER( ZZnAjLO       = 9472 )
      PARAMETER( WPZAjLO       = 9480 )
      PARAMETER( WPhadZAjLO    = 9481 )
      PARAMETER( WPZhadAjLO    = 9482 )
      PARAMETER( WMZAjLO       = 9490 )
      PARAMETER( WMhadZAjLO    = 9491 )
      PARAMETER( WMZhadAjLO    = 9492 )
      PARAMETER( WPAAjLO       = 9500 )
      PARAMETER( WMAAjLO       = 9510 )
      PARAMETER( ZAAjLO        = 9520 )
      PARAMETER( ZnAAjLO       = 9521 )
      PARAMETER( AAAjLO        = 9530 )

      PARAMETER( WPAAj         = 800  )
      PARAMETER( WMAAj         = 810  )
      PARAMETER( WPAAjjLO      = 9800 )
      PARAMETER( WMAAjjLO      = 9810 )


! Higgs plus photon plus jets
      INTEGER HAjj, HAjj_AA, HAjj_mu, HAjj_tau, HAjj_bbar
      INTEGER HAjj_WW, HAjj_ZZ_ll, HAjj_ZZ_lnu
      INTEGER HAjjjLO, HAjjjLO_AA, HAjjjLO_mu, HAjjjLO_tau, HAjjjLO_bbar
      INTEGER HAjjjLO_WW, HAjjjLO_ZZ_ll, HAjjjLO_ZZ_lnu

      PARAMETER( HAjj          = 2100 )
      PARAMETER( HAjj_AA       = 2101 )
      PARAMETER( HAjj_mu       = 2102 )
      PARAMETER( HAjj_tau      = 2103 )
      PARAMETER( HAjj_bbar     = 2104 )
      PARAMETER( HAjj_WW       = 2105 )
      PARAMETER( HAjj_ZZ_ll    = 2106 )
      PARAMETER( HAjj_ZZ_lnu   = 2107 )
      PARAMETER( HAjjjLO       = 92100)
      PARAMETER( HAjjjLO_AA    = 92101)
      PARAMETER( HAjjjLO_mu    = 92102)
      PARAMETER( HAjjjLO_tau   = 92103)
      PARAMETER( HAjjjLO_bbar  = 92104)
      PARAMETER( HAjjjLO_WW    = 92105)
      PARAMETER( HAjjjLO_ZZ_ll = 92106)
      PARAMETER( HAjjjLO_ZZ_lnu= 92107)
      
! diboson plus jet
      INTEGER WWJ, WMZJ, WPZJ, WMAJ, WPAJ, ZZJ
      INTEGER WWJJLO, WMZJLO, WPZJLO, WMAJLO, WPAJLO
      INTEGER ZZJJLO
      INTEGER WPhadWMJ, WPWMhadJ
      INTEGER WPhadZJ, WPZhadJ, WMhadZJ, WMZhadJ
      INTEGER WPhadWMJJLO, WPWMhadJJLO
      INTEGER WPhadZJjLO, WPZhadJjLO, WMhadZJjLO, WMZhadJjLO
      PARAMETER ( WWJ         = 600 )
      PARAMETER ( WPHADWMJ    = 601 )
      PARAMETER ( WPWMHADJ    = 602 )
      PARAMETER ( WMAJ        = 610 )
      PARAMETER ( WPAJ        = 620 )
      PARAMETER ( WMZJ        = 630 )
      PARAMETER ( WMhadZJ     = 631 )
      PARAMETER ( WMZhadJ     = 632 )
      PARAMETER ( WPZJ        = 640 )
      PARAMETER ( WPhadZJ     = 641 )
      PARAMETER ( WPZhadJ     = 642 )
      PARAMETER ( ZZJ         = 650 )
      PARAMETER ( WWJJLO      = 9600 )
      PARAMETER ( WPHADWMJJLO = 9601 )
      PARAMETER ( WPWMHADJJLO = 9602 )
      PARAMETER ( WMAJLO      = 9610 )
      PARAMETER ( WPAJLO      = 9620 )
      PARAMETER ( WMZJLO      = 9630 )
      PARAMETER ( WMhadZJjLO  = 9631 )
      PARAMETER ( WMZhadJjLO  = 9632 )
      PARAMETER ( WPZJLO      = 9640 )
      PARAMETER ( WPhadZJjLO  = 9641 )
      PARAMETER ( WPZhadJjLO  = 9642 )
      PARAMETER ( ZZJJLO      = 9650 )

! gluon fusion
      INTEGER GFWWj, GFZZj
      INTEGER GFWPHADWMJ, GFWPWMHADJ
      PARAMETER( GFWWj         = 4600 )
      PARAMETER( GFWPHADWMJ    = 4601 )
      PARAMETER( GFWPWMHADJ    = 4602 )
      PARAMETER( GFZZj         = 4650 )

! W+- plus jet
      INTEGER WPJ, WMJ, WPJJLO, WMJJLO

      PARAMETER ( WPJ   = 1630 )
      PARAMETER ( WMJ   = 1640 )
      PARAMETER ( WPJJLO   = 91630 )
      PARAMETER ( WMJJLO   = 91640 )

! W+- H plus jet
      INTEGER WPHJ, WPHJ_AA, WPHJ_mu, WPHJ_tau, WPHJ_bbar, WPHJ_WW, WPHJ_ZZ_ll, WPHJ_ZZ_lnu
      INTEGER WMHJ, WMHJ_AA, WMHJ_mu, WMHJ_tau, WMHJ_bbar, WMHJ_WW, WMHJ_ZZ_ll, WMHJ_ZZ_lnu
      INTEGER WPHJJLO, WPHJJLO_AA, WPHJJLO_mu, WPHJJLO_tau, WPHJJLO_bbar 
      INTEGER WPHJJLO_WW, WPHJJLO_ZZ_ll, WPHJJLO_ZZ_lnu
      INTEGER WMHJJLO, WMHJJLO_AA, WMHJJLO_mu, WMHJJLO_tau, WMHJJLO_bbar 
      INTEGER WMHJJLO_WW, WMHJJLO_ZZ_ll, WMHJJLO_ZZ_lnu

      PARAMETER( WPHJ           = 1600  )
      PARAMETER( WPHJ_AA        = 1601  )
      PARAMETER( WPHJ_mu        = 1602  )
      PARAMETER( WPHJ_tau       = 1603  )
      PARAMETER( WPHJ_bbar      = 1604  )
      PARAMETER( WPHJ_WW        = 1605  )
      PARAMETER( WPHJ_ZZ_ll     = 1606  )
      PARAMETER( WPHJ_ZZ_lnu    = 1607  )

      PARAMETER( WPHJJLO           = 91600  )
      PARAMETER( WPHJJLO_AA        = 91601  )
      PARAMETER( WPHJJLO_mu        = 91602  )
      PARAMETER( WPHJJLO_tau       = 91603  )
      PARAMETER( WPHJJLO_bbar      = 91604  )
      PARAMETER( WPHJJLO_WW        = 91605  )
      PARAMETER( WPHJJLO_ZZ_ll     = 91606  )
      PARAMETER( WPHJJLO_ZZ_lnu    = 91607  )

      PARAMETER( WMHJ           = 1610  )
      PARAMETER( WMHJ_AA        = 1611  )
      PARAMETER( WMHJ_mu        = 1612  )
      PARAMETER( WMHJ_tau       = 1613  )
      PARAMETER( WMHJ_bbar      = 1614  )
      PARAMETER( WMHJ_WW        = 1615  )
      PARAMETER( WMHJ_ZZ_ll     = 1616  )
      PARAMETER( WMHJ_ZZ_lnu    = 1617  )

      PARAMETER( WMHJJLO           = 91610  )
      PARAMETER( WMHJJLO_AA        = 91611  )
      PARAMETER( WMHJJLO_mu        = 91612  )
      PARAMETER( WMHJJLO_tau       = 91613  )
      PARAMETER( WMHJJLO_bbar      = 91614  )
      PARAMETER( WMHJJLO_WW        = 91615  )
      PARAMETER( WMHJJLO_ZZ_ll     = 91616  )
      PARAMETER( WMHJJLO_ZZ_lnu    = 91617  )


! gluon fusion
      INTEGER GFHjj, GFHjj_AA, GFHjj_mu, GFHjj_tau, GFHjj_bbar
      INTEGER GFHjj_WW, GFHjj_ZZ_ll, GFHjj_ZZ_lnu
      PARAMETER( GFHjj         = 4100 )
      PARAMETER( GFHjj_AA      = 4101 )
      PARAMETER( GFHjj_mu      = 4102 )
      PARAMETER( GFHjj_tau     = 4103 )
      PARAMETER( GFHjj_bbar    = 4104 )
      PARAMETER( GFHjj_WW      = 4105 )
      PARAMETER( GFHjj_ZZ_ll   = 4106 )
      PARAMETER( GFHjj_ZZ_lnu  = 4107 )

! QCD induced processes VVjj
      INTEGER QCDWPZjj, QCDWMZjj, QCDWPZjjjLO, QCDWMZjjjLO
      INTEGER QCDWPWPjj,QCDWPWPjjjLO,QCDWMWMjj,QCDWMWMjjjLO
      INTEGER QCDZZjj_ll,QCDZZjj_lnu,QCDZZjjjLO_ll,QCDZZjjjLO_lnu
      INTEGER QCDWPAjj, QCDWMAjj, QCDWPAjjjLO, QCDWMAjjjLO
      INTEGER QCDZAjj_l,  QCDZAjjjLO_l, QCDAAjj, QCDAAjjjLO
      INTEGER QCDZAjj_n,  QCDZAjjjLO_n
      !PARAMETER( QCDWPWMjj     =  3200 )
      PARAMETER( QCDZZjj_ll    =  3210 )
      PARAMETER( QCDZZjj_lnu   =  3211 )
      PARAMETER( QCDZZjjjLO_ll = 93210 )
      PARAMETER( QCDZZjjjLO_lnu= 93211 )
      PARAMETER( QCDWPZjj      =  3220 )
      PARAMETER( QCDWMZjj      =  3230 )
      PARAMETER( QCDWPZjjjLO   = 93220 )
      PARAMETER( QCDWMZjjjLO   = 93230 )
      PARAMETER( QCDWPAjj      =  3270  )
      PARAMETER( QCDWMAjj      =  3280  )
      PARAMETER( QCDWPAjjjLO   = 93270  )
      PARAMETER( QCDWMAjjjLO   = 93280  )
      PARAMETER( QCDZAjj_l    =   3290 )
      PARAMETER( QCDZAjjjLO_l =  93290 )
      PARAMETER( QCDZAjj_n    =   3291 )
      PARAMETER( QCDZAjjjLO_n =  93291 )
      PARAMETER( QCDAAjj       =  3299 )
      PARAMETER( QCDAAjjjLO    = 93299 )
      PARAMETER( QCDWPWPjj     =  3250 )
      PARAMETER( QCDWMWMjj     =  3260 )
      PARAMETER( QCDWPWPjjjLO  = 93233 )
      PARAMETER( QCDWMWMjjjLO  = 93244 )

! QCD induced processes Vjj
      INTEGER QCDZjj_l, QCDZjj_nu, QCDWPjj,QCDWPjjjLO,QCDWMjj,QCDWMjjjLO
      INTEGER QCDZjjjLO_l, QCDZjjjLO_nu
      PARAMETER( QCDZjj_l      = 3120 )
      PARAMETER( QCDZjjjLO_l   =93120 )
      PARAMETER( QCDZjj_nu     = 3121 )
      PARAMETER( QCDZjjjLO_nu  =93121 )
      PARAMETER( QCDWPjj       = 3130 )
      PARAMETER( QCDWPjjjLO    =93130 )
      PARAMETER( QCDWMjj       = 3140 )
      PARAMETER( QCDWMjjjLO    =93140 )


