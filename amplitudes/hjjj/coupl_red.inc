c====================================================================
c	remove all stuff which doesn't seem to be needed 
c
c  Define common block containing all coupling constants and masses.
c
c====================================================================
c
      double complex       gg(2)
c      double precision     alpha,ee, sin2w, gfermi, alfas,g
c      common /COUPL_BASIC/ gg, alpha,ee, sin2w, gfermi, alfas,g   
c      
      double precision     alpha,ee, sin2w, gfermi, alfas,g
      common /COUPL_BASIC/ gg, alpha,ee, sin2w, gfermi, alfas,g   
c
      double precision     hmass, wmass, zmass, amass,
     &                     tmass, bmass, lmass, cmass
      common /COUPL_MASS/  hmass, wmass, zmass, amass,
     &                     tmass, bmass, lmass, cmass

      double precision     hwidth, wwidth, zwidth, 
     &                     twidth, lwidth, awidth
      common /COUPL_WIDTH/ hwidth, wwidth, zwidth, 
     &                     twidth, lwidth, awidth

      double complex       gal(2), gad(2), gau(2), gwf(2),
     &                     gzn(2), gzl(2), gzd(2), gzu(2)
      double precision     gw, gwwa, gwwz
      common /COUPL_GAUGE/ gal   , gad   , gau   , gwf   ,
     &                     gzn   , gzl   , gzd   , gzu   ,
     &                     gw, gwwa, gwwz

      double complex       gwwh, gzzh, gwwhh, gzzhh, ghhh, ghhhh
      common /COUPL_SCAL/  gwwh, gzzh, gwwhh, gzzhh, ghhh, ghhhh



