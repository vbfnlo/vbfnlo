c Common block for spinor-information for the process wgammajet
c
c  lepton spinors and W+- polarization vectors
      double complex wai(6), waf(6,-1:1), wem(6), wve(6),wafpolhel(0:3,-1:1)
      double complex wavew(6), wpcurr(6,0:3,2),acurr(6,0:3,2)
      double complex zcurr(6,0:3,2),wmcurr(6),currlepw(0:3),wafpolc(0:3,-1:1)
      double complex wafpolcz(0:3,-1:1,2,-1:1),currlepwga(0:3,-1:1)
      double complex glupolc(0:3,2,-1:1),wmum(6),wmup(6),zcurrz(6)
      double complex currlepz(0:3,-1:1),wavez(6,-1:1),wavea(6,-1:1)
      double complex currlepa(0:3,-1:1)

      double precision wafpol(0:3,-1:1),glupol(0:3,2,-1:1)

      common /wmaspinor/ wai, waf, wem, wve
      common /wmaspinor/ wafpolhel,wavew,wpcurr
      common /wmaspinor/ acurr,zcurr,wmcurr,currlepw,currlepwga
      common /wmaspinor/ wafpolc,wafpolcz,glupolc,wmum
      common /wmaspinor/ wmup,zcurrz,currlepz
      common /wmaspinor/ wavez,wavea,currlepa,wafpol,glupol

      double precision wemnumom(0:4),wemnugamom(0:4),photfourvec(0:4)
      common /vgajmomenta/ wemnumom, wemnugamom,photfourvec

c b-quarks
      integer lbquark
      common /bquarks/ lbquark

