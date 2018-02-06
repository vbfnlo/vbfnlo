!====================================================================
!
!  Define common block containing all coupling constants and masses.
!
!====================================================================
!
      double complex       gg(2)
      double precision     alpha,ee, sin2w, gfermi, alfas,g
      common /COUPL_BASIC/ gg, alpha,ee, sin2w, gfermi, alfas,g  
!
      double precision     hmass, wmass, zmass, amass, tmass, bmass, lmass, cmass, h2mass
      common /COUPL_MASS/  hmass, wmass, zmass, amass, tmass, bmass, lmass, cmass, h2mass

      double precision     hwidth, wwidth, zwidth, twidth, lwidth, awidth, h2width
      common /COUPL_WIDTH/ hwidth, wwidth, zwidth, twidth, lwidth, awidth, h2width

      double complex       gal(2), gad(2), gau(2), gwf(2), gzn(2), gzl(2), gzd(2), gzu(2)
      double precision     gw, gwwa, gwwz
      common /COUPL_GAUGE/ gal   , gad   , gau   , gwf   , gzn   , gzl   , gzd   , gzu   , gw, gwwa, gwwz

      double complex       gwfc(2),  gwfs(2)
      common /coupl_ckm/   gwfc,     gwfs

      double complex       gwwh, gzzh, gwwhh, gzzhh, ghhh, ghhhh, gwwh2, gzzh2
      common /COUPL_SCAL/  gwwh, gzzh, gwwhh, gzzhh, ghhh, ghhhh, gwwh2, gzzh2

      double complex       ghtop(2), ghbot(2), ghtau(2), ghcha(2)
      common /COUPL_YUK/   ghtop   , ghbot   , ghtau   , ghcha

      double precision     xzmass, xwmass
      common /COUPL_XMASS/ xzmass, xwmass

      double complex       xzl(2) , xzb(2) , xzt(2) , xwpq(2), xwmq(2), xwpl(2), xwml(2)
      common /COUPL_XFFV/  xzl    , xzb    , xzt    , xwpq   , xwmq   , xwpl   , xwml

      double complex       xzhz, xwhwp, xwhwm
      common /COUPL_XVSS/  xzhz, xwhwp, xwhwm

      double complex       xwzwp, xwzwm, xwawp, xwawm
      common /COUPL_XVVS/  xwzwp, xwzwm, xwawp, xwawm

      double complex       xwzhwp, xwzhwm, xwahwp, xwahwm
      common /COUPL_XVVSS/ xwzhwp, xwzhwm, xwahwp, xwahwm

! settings for two-Higgs model
      integer gsba, gcba
      logical with_secondhiggs, use_secondhiggs
      common /COUPL_TWOHIGGS/ gsba, gcba, with_secondhiggs, use_secondhiggs

!  switch for MadGraph comparison in MHLM
      logical              mgtest
      parameter (mgtest = .false.)
