c last modified by Vera January 2009
c
c  precalculated lepton tensors forr VVV --> 6 leptons

      integer it, ie, iu

c  lepton spinors, polarization vectors, momenta
      double complex wep, wve, wmu, wvm, ztp, ztm
      double complex zep, zem, zmm, zmp, ae, ze, amu, zmu
      double complex wp, wm, wpp, wmp, xp, xm, ztau, atau, wbos
      double complex wa, wz, zp, wap
      double complex ztaup, ataup, xz, xa, xaz
      double complex atau1, atau2, atau3, w0, atau1p, atau2p, w0p

      double precision qp, qm, qz, qww, qwwz, qze, qzmu, qztau, qwbos
      double precision qzzw, qzzz, qzz, qzza, qa, qwz, qwza
      double precision qw0, qa1, qa2, qa3, qaa, qaaw, qzaa
      double complex aww, zww, awwz, zwwz, wpwz, wmwz, aza, zaz
      double complex wpwa, wmwa, awwa, zwwa
      double complex wpwa1, wmwa1, wpwa2, wmwa2, wpwaa, wmwaa
      double complex azza,zzza,azz,zzztens,aaze,zzae,aazm,zzam
      double complex aaze1,aaze2,zzae1,zzae2
      double complex azzz, zzzz, azezm, azezt, azmzt, azaa,  zzaa
      double complex zzzz_higgs1, zzzz_higgs2, zzzz_higgs3
      double complex zzezm, zzezt, zzmzt
      double complex wpzzw, wmzzw, wzzw, wpzw1, wpzw2, wmzw1, wmzw2
      double complex wzw1, wzw2
      double complex zmup, aep, amup, zepp, wbosp, xze, xae, xzmu
      double complex xaon, xw, xamu,xatau1,xatau2
      double complex aaaa, zaaa

      common /wwzspinor/ wep(6,11), wve(6,11), wmu(6,11), wvm(6,11),
     1                   ztp(6,11), ztm(6,11), wp(6,11), wm(6,11),
     2                   ztau(6,11), atau(6,11),wpp(6,11), wmp(6,11), 
     3                   wap(6,11), xp, xm,ataup(6,11), ztaup(6,11), 
     4                   xz, xa, xaz, qp(0:4,11), qm(0:4,11),
     5                   qz(0:4,11), qww(0:4,11), qwwz(0:4,11)


      common /zzzspinor/ zep(6,11), zem(6,11), zmm(6,11), zmp(6,11),
     1                   ae(6,11), ze(6,11), amu(6,11), zmu(6,11),
     2                   qze(0:4,11), qzmu(0:4,11), qztau(0:4,11), 
     3                   qzzz(0:4,11)

      common /aawspinor/ atau1(6,11), atau2(6,11), w0(6,11), atau1p(6,11), 
     1	                 atau2p(6,11), w0p(6,11), xatau1, xatau2, atau3(6,11), 
     2                   qzaa(0:4,11), qa3(0:4,11), qw0(0:4,11), qa1(0:4,11), 
     3                   qa2(0:4,11), qaa(0:4,11), qaaw(0:4,11)


      common /zzwspinor/ wbos(6,11), zepp(6,11), zmup(6,11), aep(6,11), 
     1                   amup(6,11), wbosp(6,11), xze, xae, 
     2                   xzmu, xamu, xw, xaon, qzza(0:4,11),
     3                   qwbos(0:4,11), qzzw(0:4,11), qzz(0:4,11)

      common /wzaspinor/ wa(6,11), wz(6,11), zp(6,11), 
     1                   qa(0:4,11), qwz(0:4,11), qwza(0:4,11)

      common /tensww/ aww(0:5,11), zww(0:5,11),
     1                awwz(0:5,11), zwwz(0:5,11),
     2                wpwz(0:5,11), wmwz(0:5,11)


      common /tenszzz/ azzz(0:5,11), zzzz(0:5,11),
     1                 azezm(0:5,11), azezt(0:5,11), azmzt(0:5,11),
     2                 zzezm(0:5,11), zzezt(0:5,11), zzmzt(0:5,11),
     3                 zzzz_higgs1(0:5,11), zzzz_higgs2(0:5,11),
     4                 zzzz_higgs3(0:5,11) 

      common /tensaaw/ wpwa1(0:5,11), wpwa2(0:5,11), wmwa1(0:5,11),
     1	               wmwa2(0:5,11), wpwaa(0:5,11), wmwaa(0:5,11)

      common /tensvva/ awwa(0:5,11), zwwa(0:5,11),
     1                 wpwa(0:5,11), wmwa(0:5,11),
     2		       azza(0:5,11), zzza(0:5,11),
     3                 azz(0:5,11), zzztens(0:5,11),
     4                 aaze(0:5,11), zzae(0:5,11),
     5                 aaze1(0:5,11), zzae1(0:5,11),
     6                 aaze2(0:5,11), zzae2(0:5,11),
     7                 aazm(0:5,11), zzam(0:5,11),
     8                 aza(0:5,11), zaz(0:5,11),
     9                 azaa(0:5,11), zzaa(0:5,11), aaaa(0:5,11), zaaa(0:5,11)


      common /tenszzw/ wzzw(0:5,11), wzw1(0:5,11), wzw2(0:5,11)


c-------------------------www-momenta and wavefunctions--------------------
c  lepton spinors and W+,W-, Z/A polarization vectors

      double complex wep1, wve1, wmu1, wvm1, wtp, wvt
      double complex wp1, wp2, w1m
      double complex wpp1, wpp2, w1mp, w2mp,
     1               xp1, xp2, x1m, x2m
      double precision qp1, qp2, q1m, qwp1wm, qwp2wm,
     1                 qwpwmwp, qwwwp

      double complex wep2, wve2, wmu2, wvm2, wtm,wvtbar
      double complex w1p, wm1, wm2
      double complex wmm1, wmm2, w1pm, w2pm,
     1               xm1, xmm2, x1p, x2p
      double precision q1p, qm1, qm2, qwpwm1, qwpwm2,
     1                 qwpwmwm, qwwwm

      common /wwwpspinor/ wep1(6,11), wve1(6,11), wmu1(6,11), wvm1(6,11),
     1                    wtp(6,11), wvt(6,11), wp1(6,11), wp2(6,11),
     2                    w1m(6,11), wpp1(6,11), wpp2(6,11), w1mp(6,11),
     3                    w2mp(6,11), xp1, xp2, x1m, x2m, qp1(0:4,11), 
     4                    qp2(0:4,11), q1m(0:4,11), qwp1wm(0:4,11), 
     5                    qwp2wm(0:4,11), qwpwmwp(0:4,11), qwwwp(0:4,11)


      common /wwwmspinor/ wep2(6,11), wve2(6,11), wmu2(6,11), wvm2(6,11),
     1                    wtm(6,11), wvtbar(6,11), wm1(6,11), wm2(6,11),
     2                    w1p(6,11), wmm1(6,11), wmm2(6,11), w1pm(6,11),
     3                    w2pm(6,11), xm1, xmm2, x1p, x2p,q1p(0:4,11),
     4                    qm1(0:4,11),qm2(0:4,11),qwpwm1(0:4,11),qwpwm2(0:4,11),
     5                    qwpwmwm(0:4,11),qwwwm(0:4,11)


c---------------------leptonic tensors for w+w-w+/- --------------------------
      double complex WWW_aww, WWW_zww, wptowww, wmtowww
      common /tenswww/ WWW_aww(0:5,11,11), WWW_zww(0:5,11,11),
     1                wptowww(0:5,11),wmtowww(0:5,11)


c---------------------common block for testing purposes-----------------------
      double complex WWWtamp,WWWtampr
      common /WWWtest/ WWWtamp(6,-1:1,2)
      common /WWWrtest/ WWWtampr(6,-1:1,-1:1,2)

      double complex tamp,tampr,tamprr
      integer ip,ipr,iprr
      common /VVVtest/ tamp(10,-1:1,2), ip
      common /VVVrtest/ tampr(10,-1:1,-1:1,2), ipr
      common /VVVrrtest/ tamprr(10,-1:1,-1:1,-1:1,2), iprr

c---------------------common block for vvvjet-----------------------
c  lepton spinors and W+- polarization vectors
      double complex wai, waf, wem, wafpolhel, wavew,
     1               wpcurr,acurr,zcurr,wmcurr,currlepw,wafpolc,currlepwga,
     2               glupolc
      double precision wafpol,glupol
      common /wmaspinor/ wai(6), waf(6,-1:1), wem(6),
     1                   wafpol(0:3,-1:1),glupol(0:3,2,-1:1),
     2                   wafpolhel(0:3,-1:1),wavew(6),wpcurr(6,0:3,2),
     3                   acurr(6,0:3,2),zcurr(6,0:3,2),wmcurr(6),
     4                   currlepw(0:3),currlepwga(0:3,-1:1),
     5                   wafpolc(0:3,-1:1),glupolc(0:3,2,-1:1)

      double precision wemnumom,wemnuga1mom,wemnuga2mom,wemnugagamom,
     1                 phot1fourvec,phot2fourvec,
     2                 wemmom,wnumom
      common /vgagajmomenta/ wemnumom(0:4), wemnuga1mom(0:4),
     1                       wemnuga2mom(0:4), wemnugagamom(0:4),
     2                       phot1fourvec(0:4), phot2fourvec(0:4),
     3                       wemmom(0:4), wnumom(0:4)
