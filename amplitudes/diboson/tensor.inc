c Dieter's version from Apr. 12, 2005, extended for real emission contributions
c modified by Vera April 2006
c modified by Johannes July 2011
c modified by Michael R May 2012
c  lepton spinors and W+- polarization vectors
      double complex wep, wen, wve, wvebar, wmu, wmubar
      double complex wvm, wp, wm, wpp, wmp, xp, xm, zn, an, 
     1               we, wmup, wmum, wa, wz
      double complex atau1, atau2
      double precision qp, qm, qz, qa, qa1, qa2, qaa, qww
      common /wpmspinor/ wep(6,max_kin), wve(6,max_kin), wmu(6,max_kin),
     &      wvm(6,max_kin), wp(6,max_kin), wm(6,max_kin), 
     &      atau1(6,max_kin), atau2(6,max_kin),
     &      wpp(6,max_kin), wmp(6,max_kin), 
     &      zn(6,max_kin), wmubar(6,max_kin), 
     &      an(6,max_kin), 
     &      wen(6,max_kin),wvebar(6,max_kin),
     &      we(6,max_kin), wmup(6,max_kin), wmum(6,max_kin),
     &      wa(6,max_kin),wz(6,max_kin),
     &      qp(0:4,max_kin), qm(0:4,max_kin),
     &      qz(0:4,max_kin) ,qa(0:4,max_kin),
     &      qa1(0:4,max_kin), qa2(0:4,max_kin), qaa(0:4,max_kin),
     &      qww(0:4,max_kin),
     &      xp, xm
c  leptonic tensors
      double complex aww, zww, wtowz, wtowa, z1ww, w4l, z2ww
      integer ie,iu
      common /tensww/ aww(0:5,max_kin), zww(0:5,max_kin), wtowz(0:5,max_kin),
     1                wtowa(0:5,max_kin), z1ww(0:5,max_kin),z2ww(0:5,max_kin),  
     2                w4l(0:9,max_kin), ie, iu

c for testing purposes
      double complex tamp, tampr, tamprr
      integer ip, ipr, iprr
      common /wwtest/ tamp(10,-1:1,4), ip
      common /wwrtest/ tampr(20,-1:1,-1:1,4), ipr
      common /WWWrrtest/ tamprr(36,-1:1,-1:1,-1:1,4), iprr


c---------------------common block for vhj -----------------------
      double complex hig, wtowh, wtoonlyw
      double precision qh, qwh,qw0
      common /vhj/ hig(6,11), wtowh(0:5,max_kin), ! TODO: change max_kin to 1?!
     1       qh(0:4,11),qwh(0:4,11),qw0(0:4,11),
     2       wtoonlyw(0:5,max_kin)

      double precision mad_series_t, mad_series_p, mad_tau
      double complex mad_gh
      common /madtoploop/ mad_gh(2)

c---------------------common block for vhj@NLO-----------------------
c  lepton spinors and W+- polarization vectors
      double complex glupolc
      common /wmhspinor/ glupolc(0:3,2,-1:1)

      double precision wemnumom,wemnuhmom,
     2                 wemmom,wnumom,
     3                 higmom
      common /vhjmomenta/ wemnumom(0:4), wemnuhmom(0:4),
     3                       wemmom(0:4), wnumom(0:4),
     4                       higmom(0:4)

