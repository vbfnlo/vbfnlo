      subroutine coupsm(idef)
c
c This subroutine sets up the coupling constants of the STANDARD MODEL.
c
c input:
c       real    sin2w             : square of sine of the weak angle
c       real    alpha             : EM fine structure constant
c       real    Gfermi            : G_Fermi
c       real    zmass,wmass       : weak boson masses
c                                 : ONLY zmass OR wmass is input;
c                                 : the other is calculated
c       real    hmass             : Higgs mass
c       real    zwidth            : Z boson width
c       real    tmass,bmass,lmass : fermion masses
c       real    xzmass, xwmass    : mass of Goldstone bosons
c
c output:
c       real    gw                : weak coupling constant
c       real    gwwa              : dimensionless WWA  coupling
c       real    gwwz              : dimensionless WWZ  coupling
c       complex gwwh              : dimensionful  WWH  coupling
c       complex gzzh              : dimensionful  ZZH  coupling
c       complex ghhh              : dimensionful  HHH  coupling
c       complex gwwhh             : dimensionful  WWHH coupling
c       complex gzzhh             : dimensionful  ZZHH coupling
c       complex ghhhh             : dimensionless HHHH coupling
c       complex ghtop(2)          : Yukawa coupling of (L,R) top   quark
c       complex ghbot(2)          : Yukawa coupling of (L,R) bottom quark
c       complex ghtau(2)          : Yukawa coupling of (L,R) tau  lepton
c       complex gal(2)            : coupling with A of charged leptons
c       complex gau(2)            : coupling with A of up-type quarks
c       complex gad(2)            : coupling with A of down-type quarks
c       complex gwf(2)            : coupling with W-,W+ of fermions
c       complex gzn(2)            : coupling with Z of neutrinos
c       complex gzl(2)            : coupling with Z of charged leptons
c       complex gzu(2)            : coupling with Z of up-type quarks
c       complex gzd(2)            : coupling with Z of down-type quarks
c       real    g                 : QCD 3-,4-gluon coupling
c       complex gg(2)             : QCD gqq coupling (L,R)
c       complex xzhz              : Goldstone z   coupling  to ZH
c       complex xwhwp,xwhwm       : Goldstone w+- couplings to WH
c       complex xwzwp,xwzwm       : Goldstone w+- couplings to WZ
c       complex xwawp,xwawm       : Goldstone w+- couplings to WA
c       complex xwzhwp,xwzhwm     : Goldstone w+- couplings to WZH
c       complex xwahwp,xwahwm     : Goldstone w+- couplings to WAH
c       complex xwpq,xwpl         : Goldstone w+  couplings to t,b,tau
c       complex xwmq,xwml         : Goldstone w-  couplings to t,b,tau
c       complex xzt,xzb,xzl       : Goldstone z   couplings to t,b,tau
c
c ----------------------------------------------------------------------
c
c input parameters
c
      implicit none
      integer  idef
c
c calculated couplings
c
#include "coupl.inc"
      double precision  ee2, ez, ey, sw, cw, sc2, v, rho, alpha_s
      double precision  gwne, gwud, lambda, lam4, xt, rew, rqcd
      double precision  decw, w_w_nl, w_w_tau, w_w_ud, w_w_cs
      double precision  decz, w_z_nn, w_z_ll, w_z_tau
      double precision  w_z_uu, w_z_dd, w_z_cc, w_z_bb
      double precision  heavy, a, b, c, dum
      double precision  alphas, alfa, alfaw, mfrun
      external          alphas, alfa, alfaw, mfrun
c
c predefined SM parameters (these values may change over time)
c
      double precision  mz, wz, mh
      double precision  mt, wt, mb, ml, mc
      double precision  hals, hale, hss, hgfermi
      parameter( mz = 91.188d0, wz = 2.495d0 )
      parameter( mh = 120.0d0 )
      parameter( mt = 174.3d0, mb = 4.7d0, ml = 1.777d0 )
      parameter( mc = 1.42d0)
      parameter( hals = 0.1185d0, hale = 1d0/128.9d0 )
      parameter( hss = 0.2312d0, hgfermi = 1.16639d-5 )
c
c HDECAY variables
c
      logical           lhdecay
      integer           ihiggs, nnlo, ipole
      common /flag/     ihiggs, nnlo, ipole
      integer           ionsh, ionwz, iofsusy
      common /onshell/  ionsh, ionwz, iofsusy
      double precision  gf, alph, almass, ammuon, amz, amw
      common /param/    gf, alph, almass, ammuon, amz, amw
      double precision  ams, amc, amb, amt
      common /massesC/   ams, amc, amb, amt
      double precision  gamc0, gamt0, gamt1, gamw, gamz
      common /wzwdth/   gamc0, gamt0, gamt1, gamw, gamz
      double precision  hdals, hdmhbeg, hdmhend, tgbet
      common /hdparms/  hdals, hdmhbeg, hdmhend, tgbet
      double precision  SMBRB,SMBRL,SMBRM,SMBRS,SMBRC,SMBRT,
     &                  SMBRG,SMBRGA,SMBRZGA,SMBRW,SMBRZ,SMWDTH
      common /widthsm/  SMBRB,SMBRL,SMBRM,SMBRS,SMBRC,SMBRT,
     &                  SMBRG,SMBRGA,SMBRZGA,SMBRW,SMBRZ,SMWDTH
c
c constants
c
      double complex  ci
      parameter( ci = ( 0.0d0, 1.0d0 ) )
      double precision  Zero, One, Two, Three, Four, Half, Rt2
      parameter( Zero = 0.0d0, One = 1.0d0, Two = 2.0d0 )
      parameter( Three = 3.0d0, Four = 4.0d0, Half = 0.5d0 )
      double precision  Pi, Fourpi
      parameter( Pi = 3.14159265358979323846d0 )
      parameter( Fourpi = Four * Pi )
c
c define phase space reduction factor for heavy fermionic decays
c
      heavy(a,b,c) = ( One - Half*(b**2+c**2) - Half*(b**2-c**2)**2
     &                 + Three*b*c*(a**2 - One)/(a**2 + One)        )
     &               * sqrt( (One-b**2-c**2)**2 - Four * b**2 * c**2 )
c
c convenient parameter
c
      Rt2 = sqrt(Two)
c
c     Photon Mass
c
      AMASS = zero
c
c HELAS_CHECK mandatory tree-level relations
c
#ifdef HELAS_CHECK
      idef = 0
#endif
c
c call hdecay or not 
c
      lhdecay=.false.
c
c mass & width catches
c
      if ( hmass.le.Zero ) hmass = mh
      if ( tmass.le.Zero .or. idef.eq.0 ) tmass = mt
      if ( bmass.le.Zero .or. idef.eq.0 ) bmass = mb
      if ( cmass.le.Zero .or. idef.eq.0 ) cmass = mc
      if ( lmass.le.Zero .or. idef.eq.0 ) lmass = ml
      if ( gfermi.le.Zero ) gfermi = hgfermi
c
c option 0 - use predefined MSbar values with tree-level relations
c
      if ( idef.eq.0 ) then
         write(6,*)
         write(6,*) '================================================'
         write(6,*) '| HELAS COUPSM called with option 0 (default): |'
         write(6,*) '| predefined MS-bar parameters will be used,   |'
         write(6,*) '| all parameters calculated at tree level.     |'
         write(6,*) '================================================'
         write(6,*)
         write(6,*) 'Hardwired SM inputs:'
         write(6,17) '  sin^2(W) = ', hss
         write(6,17) '  1/alpha  = ', One/hale
         write(6,18) '  M_Z      = ', mz
         write(6,18) '  Gamma_Z  = ', wz
         write(6,18) '  m_top    = ', mt
         write(6,18) '  m_bot    = ', mb
         write(6,18) '  m_cha    = ', mc
         write(6,18) '  m_tau    = ', ml
         write(6,17) '  alpha_s  = ', hals
         write(6,*)

         sin2w  = hss
         cw     = sqrt( One - sin2w )
         alpha  = hale
         ee2    = alpha * Fourpi
         alpha_s= hals

         zmass  = mz
         zwidth = wz
         wmass  = zmass * cw
c
c option 1 - allow user to select MZ/MW,s2w,alpha,GFermi inputs in OS scheme
c
      else if ( idef.eq.1 ) then
         write(6,*)
         write(6,*) '==============================================='
         write(6,*) '| HELAS COUPSM called with option 1:          |'
         write(6,*) '| user inputs 3 selected basic SM parameters, |'
         write(6,*) '| all others calculated at tree level.        |'
         write(6,*) '==============================================='
         write(6,*)
         write(6,*) 'Inputs are:'

         alpha_s= alfas

         if ( zmass.gt.Zero ) then
            write(6,18) '  M_Z      = ', zmass
            if ( wmass.gt.Zero ) then
               write(6,18) '  M_W      = ', wmass
               if ( sin2w.gt.Zero ) then
                  write(6,*) '*** invalid set of 3 input parameters ***'
                  write(6,*) '*** (or call routine with option 3)   ***'
                  stop
               else
                  cw = wmass/zmass
                  sin2w = One - cw**2
               end if
               if ( alpha.gt.7.246d-3 ) then
                  write(6,17) '  1/alpha  = ', One/alpha
               else
                  write(6,*) '*** invalid input - alpha not set ***'
                  stop
               end if
            else if ( sin2w.gt.Zero .and. alpha.gt.7.246d-3 ) then
               write(6,17) '  1/alpha  = ', One/alpha
               write(6,17) '  sin^2(W) = ', sin2w
               cw = sqrt( One - sin2w )
               wmass = zmass * cw
            else if ( gfermi.gt.Zero .and. alpha.gt.7.246d-3 ) then
               write(6,17) '  1/alpha  = ', One/alpha
               write(6,*) '  G_Fermi  = ', gfermi
               dum = zmass**2 * Half
     &               + sqrt(   zmass**4/Four
     &                       - Pi*alpha*zmass**2/(Rt2*gfermi) )
               wmass = sqrt(dum)
               cw = wmass/zmass
               sin2w = One - cw**2
               sw = sqrt(sin2w)
            else
               write(6,*) '*** invalid set of 3 input parameters ***'
               stop
            end if
         else if ( wmass.gt.Zero .and. sin2w.gt.Zero
     &                           .and. alpha.gt.7.246d-3 ) then
            write(6,18) '  M_W      = ', wmass
            write(6,17) '  1/alpha  = ', One/alpha
            write(6,17) '  sin^2(W) = ', sin2w
            cw = sqrt( One - sin2w )
            zmass = wmass / cw
         else
            write(6,*) '*** ERROR in COUPSM inputs - halting ***'
            stop
         end if
         write(6,*)
         ee2 = alpha * Fourpi
c
c option 2 - MS-bar scheme, uses rho <> 1
c
      else if ( idef.eq.2 ) then
         write(6,*)
         write(6,*) '================================================='
         write(6,*) '| HELAS COUPSM called with option 2:            |'
         write(6,*) '| predefined MS-bar parameters will be used     |'
         write(6,*) '| (except alpha_s, taken as input from user) -  |'
         write(6,*) '| all parameters calculated w/ rho corrections. |'
         write(6,*) '| Gauge invariance not guaranteed!!!            |'
         write(6,*) '================================================='
         write(6,*)

         sin2w  = hss
         gfermi = hgfermi
         cw     = sqrt( One - sin2w )
         alpha_s= alfas
         zwidth = wz

         xt   = gfermi * ( tmass / Pi )**2 / 8.0d0 / sqrt( Two )
         dum  = alphas(tmass,hals,1) / Pi
         rew  = - xt * ( 0.73921d0 + 12.56637d0*hmass/tmass )
         rqcd = - 0.19325d0*dum - 3.9696d0*dum**2
         rho  = One / ( One - Three*xt*( One + rew + rqcd ) )
         write(6,*) 'rho  = ',rho
         write(6,*)
         wmass = sqrt(rho) * zmass * cw
         ee2 = Four*Rt2 * gfermi * sin2w * wmass**2
         alpha = ee2 / Fourpi
c
c option 3 - user must set MZ,MW,sin2w,alpha and know what they're doing!
c
      else if ( idef.eq.3 ) then
         write(6,*)
         write(6,*) '==========================================='
         write(6,*) '| HELAS COUPSM called with option 3:      |'
         write(6,*) '| user inputs M_Z,M_W,sin^2(thetaW),alpha |'
         write(6,*) '| and must know what they are doing.      |'
         write(6,*) '| Gauge invariance not guaranteed!!!      |'
         write(6,*) '==========================================='
         write(6,*)

         if ( sin2w.le.Zero ) then
            write(6,*) '*** sin^2(theta_W) undefined ***'
            stop
         end if
         if ( zmass.le.Zero ) then
            write(6,*) '*** M_Z undefined ***'
            stop
         end if
         if ( wmass.le.Zero ) then
            write(6,*) '*** M_W undefined ***'
            stop
         end if
         if ( alpha.le.Zero ) then
            write(6,*) '*** alpha undefined ***'
            stop
         end if

         cw  = sqrt( One - sin2w )
         ee2 = alpha * Fourpi
c
c else error output and halt
c
      else
         write(6,*) '*** COUPSM called w/ invalid option - halting ***'
         stop
      end if
c
c top quark width
c
      twidth = Gfermi * tmass**3 / (8.0d0*Pi*Rt2)
     &         * ( One - (wmass/tmass)**2 )**2
     &         * ( One + Two*(wmass/tmass)**2 )
     & * ( One-Two*alpha_s*(Two*Pi**2/Three-5.0d0/Two)/Three/Pi )
c
c useful values
c
      sw  = sqrt( sin2w )
      ee  = sqrt( ee2 )
      ez  = ee/(sw*cw)
      ey  = ee*(sw/cw)
      sc2 = sin2w*( One - sin2w )
      v   = Two*zmass*sqrt(sc2)/ee
      lambda = hmass**2 / (Two * v**2)
c
c vector boson couplings
c
      gw   = ee/sw
      gwwa = ee
      gwwz = ee*cw/sw
c
c gauge & higgs boson coupling constants
c
      gwwh  = dcmplx( ee2/sin2w*Half*v, Zero )
      gzzh  = dcmplx( ee2/sc2*Half*v, Zero )
      ghhh  = dcmplx( -hmass**2/v*Three, Zero )
      gwwhh = dcmplx( ee2/sin2w*Half, Zero )
      gzzhh = dcmplx( ee2/sc2*Half, Zero)
      ghhhh = ghhh/v
c
c fermion-fermion-vector couplings
c
      gal(1) = dcmplx(  ee          , Zero )
      gal(2) = dcmplx(  ee          , Zero )
      gau(1) = dcmplx( -ee*Two/Three, Zero )
      gau(2) = dcmplx( -ee*Two/Three, Zero )
      gad(1) = dcmplx(  ee/Three    , Zero )
      gad(2) = dcmplx(  ee/Three    , Zero )

      gwf(1) = dcmplx( -ee/sqrt(Two*sin2w), Zero )
      gwf(2) = dcmplx(  Zero              , Zero )

      gzn(1) = dcmplx( -ez*Half                     , Zero )
      gzn(2) = dcmplx(  Zero                        , Zero )
      gzl(1) = dcmplx( -ez*(-Half + sin2w)          , Zero )
      gzl(2) = dcmplx( -ey                          , Zero )
      gzu(1) = dcmplx( -ez*( Half - sin2w*Two/Three), Zero )
      gzu(2) = dcmplx(  ey*Two/Three                , Zero )
      gzd(1) = dcmplx( -ez*(-Half + sin2w/Three)    , Zero )
      gzd(2) = dcmplx( -ey/Three                    , Zero )
c
c fermion-fermion-Higgs couplings (complex) hff(2)
c
C      ghtop(1) = dcmplx( -tmass/v, Zero )
      write(6,15) 'top mass at Higgs mass = ',
     & mfrun(165.2d0,hmass,hals,1),' GeV'
      ghtop(1) = dcmplx( -mfrun(165.2d0,hmass,hals,1)/v, Zero )
      ghtop(2) = ghtop(1)

C      ghbot(1) = dcmplx( -bmass/v, Zero )
      write(6,15) 'bot mass at Higgs mass = ',
     & mfrun(4.2d0,hmass,hals,1),' GeV'
      ghbot(1) = dcmplx( -mfrun(4.2d0,hmass,hals,1)/v, Zero )
      ghbot(2) = ghbot(1)

C      ghtau(1) = dcmplx( -lmass/v, Zero )
      write(6,15) 'tau mass at Higgs mass = ',
     & mfrun(lmass,hmass,hals,1),' GeV'
      ghtau(1) = dcmplx( -mfrun(lmass,hmass,hals,1)/v, Zero )
      ghtau(2) = ghtau(1)

      write(6,*)
c
c QCD couplings
c
      g = sqrt( Fourpi * alpha_s )
      gg(1) = dcmplx( -g, Zero )
      gg(2) = gg(1)
c
c Z boson width
c
      if ( idef.gt.0 ) then
         decz = zmass / ( 24.0d0 * Pi )
         w_z_nn = decz * ( gzn(1)**2 + gzn(2)**2 )
         w_z_ll = decz * ( gzl(1)**2 + gzl(2)**2 )
         decz = decz * Three * ( One + alpha_s/Pi )
         w_z_uu = decz * ( gzu(1)**2 + gzu(2)**2 )
         w_z_dd = decz * ( gzd(1)**2 + gzd(2)**2 )
         dum = dble( (gzl(2)+gzl(1))/(gzl(2)-gzl(1)) )
         w_z_tau = w_z_ll * heavy( dum, lmass/zmass, lmass/zmass )
         dum = dble( (gzu(2)+gzu(1))/(gzu(2)-gzu(1)) )
         w_z_cc = w_z_uu *  heavy( dum, cmass/zmass, cmass/zmass )
         dum = dble( (gzd(2)+gzd(1))/(gzd(2)-gzd(1)) )
         w_z_bb = w_z_dd *  heavy( dum, bmass/zmass, bmass/zmass )
         zwidth =   Three*w_z_nn + Two*w_z_ll + w_z_tau
     &            + Two*w_z_dd + w_z_uu + w_z_cc + w_z_bb
      end if
c
c W boson width
c
      decw = wmass / ( 24.0d0 * Pi )
      w_w_nl = decw * ( gwf(1)**2 + gwf(2)**2 )
      dum = dble( (gwf(2)+gwf(1))/(gwf(2)-gwf(1)) )
      w_w_tau = w_w_nl * heavy( dum, lmass/wmass, Zero )
      w_w_ud = w_w_nl * Three * ( One + alpha_s/Pi )
      w_w_cs = w_w_ud * heavy( dum, cmass/wmass, Zero )
      wwidth = Two*w_w_nl + w_w_tau + w_w_ud + w_w_cs
c
c interface to HDECAY
c
      ihiggs = 0
      nnlo   = 0
      ipole  = 0

      ionsh   = 0
      ionwz   = 0
      iofsusy = 1

      gf     = gfermi
      alph   = 137.0359895d0
      almass  = lmass
      ammuon = 0.105658389d0
      amz    = zmass
      amw    = wmass

      ams = 0.190d0
      amc = 1.42d0
      amb = 4.7d0
      amt = tmass

      gamw = wwidth
      gamz = zwidth

      hdals   = alpha_s
      hdmhbeg = hmass
      hdmhend = hmass


      if (lhdecay) then
c        call hdecay
         write(*,*) 'Warning need to uncomment hdecay!'
      else
         SMWDTH = 3.7d-3    ! this corresponds to mh=120 GeV
      endif

      hwidth = SMWDTH
c
c Goldstone boson couplings
c
#ifdef HELAS_CHECK
c
c   masses
c
      xzmass = zmass
      xwmass = wmass
c
c   VHS - verified
c
      xzhz  = -ci*Half*gw/cw
      xwhwp = -ci*Half*gw
      xwhwm = dconjg(xwhwp)
c
c   VVS - verified
c
      xwzwp =  ci * gw * zmass * sin2w
      xwzwm = dconjg(xwzwp)
      xwawp = -ci * ee * wmass
      xwawm = dconjg(xwawp)
c
c   VVHS - verified
c
      xwzhwp =  ci * Half * gw**2 * sin2w / cw
      xwzhwm = dconjg(xwzhwp)
      xwahwp = -ci * Half * gw * ee
      xwahwm = dconjg(xwahwp)
c
c   FFS - verified
c
      ghtop(1) = dcmplx( -tmass/v, Zero )
      ghtop(2) = ghtop(1)

      ghbot(1) = dcmplx( -bmass/v, Zero )
      ghbot(2) = ghbot(1)

      ghtau(1) = dcmplx( -lmass/v, Zero )
      ghtau(2) = ghtau(1)

      xzt(1) = -ghtop(1)*ci
      xzt(2) = dconjg(xzt(1))
      xzb(1) =  ghbot(1)*ci
      xzb(2) = dconjg(xzb(1))
      xzl(1) =  ghtau(1)*ci
      xzl(2) = dconjg(xzl(1))

      xwpq(1) =  ci*Rt2*bmass/v
      xwpq(2) = -ci*Rt2*tmass/v
      xwmq(1) = -xwpq(2)
      xwmq(2) = -xwpq(1)
      xwpl(1) =  ci*Rt2*lmass/v
      xwpl(2) = Zero
      xwml(1) = Zero
      xwml(2) = -xwpl(1)
c
c override boson widths for BRS check
c
      hwidth = Zero
      zwidth = Zero
      wwidth = Zero
      twidth = Zero
c
c output BRS check acknowledgement
c
      write(6,*) '*** HELAS COUPSM in BRS check mode ***'
      write(6,*)
#endif
c
c output all info
c
 10   format( 1x,a,f7.3,' GeV        ',a,f7.4,' GeV' )
 11   format( 1x,a,f10.7,2x,f10.7,a,f10.7,2x,f10.7 )
 12   format( 1x,a,f6.2,a )
 13   format( 1x,a,f6.4,a )
 14   format( 1x,2(a,f10.7,', ',f10.7) )
 15   format( 1x,a,f9.5,a )
 16   format( 1x,a,f7.5 )
 17   format( 1x,a,f8.4 )
 18   format( 1x,a,f8.4,' GeV' )
 19   format( 1x,a,f6.4,a,f6.4 )
 20   format( 1x,a,f11.5,1x,f11.5 )

      write(6,*) '*** HELAS Standard Model parameters ***'
      write(6,*)
      write(6,12) '1/alpha = ',One/alpha,
     &            '  from input/MS-bar relations'
      write(6,12) '1/alpha = ',One/alfa(zmass**2),
     &            '  from 1-loop running at Q = M_Z'
      write(6,*)
      write(6,13) 'alpha_s = ',alpha_s,
     &            '  from input'
      write(6,13) 'alpha_s = ',alphas(zmass,hals,1),
     &            '  from 1-loop running at Q = M_Z'
      write(6,*)
      write(6,13) 'sin^2(theta_W) =  ',sin2w
      write(6,*)
      write(6,10) 'Z mass  =  ',zmass, 'Z width  = ',zwidth
      write(6,10) 'W mass  =  ',wmass, 'W width  = ',wwidth
      write(6,10) 'H mass  =  ',hmass, 'H width  = ',hwidth
      write(6,*)
      write(6,10) 'top    mass  =  ', tmass, 'top    width  = ', twidth
      write(6,10) 'bottom mass  =  ', bmass, 'bottom width  = ', Zero
      write(6,10) 'charm  mass  =  ', cmass, 'charm  width  = ', Zero
      write(6,10) 'tau    mass  =  ', lmass, 'tau    width  = ', Zero
      write(6,*) 'all other quark and lepton masses set to zero'
      write(6,*)
      write(6,*) 'Boson couplings:'
      write(6,*)
      write(6,*) 'gw    = ', gw, '  from input/MS-bar relations'
      write(6,*) 'gw    = ', sqrt(Fourpi*alfaw(zmass**2,1)),
     &           '  from 1-loop running at Q = M_Z'
      write(6,*)
      write(6,20) 'gwwa  = ', gwwa
      write(6,20) 'gwwz  = ', gwwz
      write(6,20) 'gwwh  = ', gwwh
      write(6,20) 'gzzh  = ', gzzh
      write(6,20) 'ghhh  = ', ghhh
      write(6,*)
      write(6,20) 'gwwhh = ', gwwhh
      write(6,20) 'gzzhh = ', gzzhh
      write(6,20) 'ghhhh = ', ghhhh
      write(6,*)
      write(6,*) 'FFV couplings:'
      write(6,*)
      write(6,11) 'ggq(L)   =  ',gg(1),  '    ggq(R)   =  ',gg(2)
      write(6,*)
      write(6,11) 'gal(L)   =  ',gal(1), '    gal(R)   =  ',gal(2)
      write(6,11) 'gau(L)   =  ',gau(1), '    gau(R)   =  ',gau(2)
      write(6,11) 'gad(L)   =  ',gad(1), '    gad(R)   =  ',gad(2)
      write(6,*)
      write(6,11) 'gwf(L)   =  ',gwf(1), '    gwf(R)   =  ',gwf(2)
      write(6,*)
      write(6,11) 'gzn(L)   =  ',gzn(1), '    gzn(R)   =  ',gzn(2)
      write(6,11) 'gzl(L)   =  ',gzl(1), '    gzl(R)   =  ',gzl(2)
      write(6,11) 'gzu(L)   =  ',gzu(1), '    gzu(R)   =  ',gzu(2)
      write(6,11) 'gzd(L)   =  ',gzd(1), '    gzd(R)   =  ',gzd(2)
      write(6,*)
      write(6,*) 'FFH couplings:'
      write(6,*)
      write(6,14) 'gHtop(L) =  ',ghtop(1), '    gHtop(R) =  ',ghtop(2)
      write(6,14) 'gHbot(L) =  ',ghbot(1), '    gHbot(R) =  ',ghbot(2)
      write(6,14) 'gHtau(L) =  ',ghtau(1), '    gHtau(R) =  ',ghtau(2)
      write(6,*)
      if ( idef.gt.0 .and. zwidth.gt.Zero .and. wwidth.gt.Zero ) then
         write(6,*) 'Z,W partial widths and effective BRs:'
         write(6,*)
         write(6,19) 'width Z->nn  = ', w_z_nn,
     &               ' GeV,   BR(Z->nn)  = ', Three*w_z_nn/zwidth
         write(6,19) 'width Z->ll  = ', w_z_ll,
     &               ' GeV,   BR(Z->ll)  = ', Two*w_z_ll/zwidth
         write(6,19) 'width Z->TT  = ', w_z_tau,
     &               ' GeV,   BR(Z->TT)  = ', w_z_tau/zwidth
         write(6,19) 'width Z->uu  = ', w_z_uu,
     &               ' GeV,   BR(Z->uu)  = ', w_z_uu/zwidth
         write(6,19) 'width Z->dd  = ', w_z_dd,
     &               ' GeV,   BR(Z->dd)  = ', Two*w_z_dd/zwidth
         write(6,19) 'width Z->cc  = ', w_z_cc,
     &               ' GeV,   BR(Z->cc)  = ', w_z_cc/zwidth
         write(6,19) 'width Z->bb  = ', w_z_bb,
     &               ' GeV,   BR(Z->bb)  = ', w_z_bb/zwidth
         write(6,*)
         write(6,19) 'width W->nl  = ', w_w_nl,
     &               ' GeV,   BR(W->nl)  = ', Two*w_w_nl/wwidth
         write(6,19) 'width W->nl  = ', w_w_nl,
     &               ' GeV,   BR(W->nT)  = ', w_w_tau/wwidth
         write(6,19) 'width W->ud  = ', w_w_ud, 
     &               ' GeV,   BR(W->ud)  = ', w_w_ud/wwidth
         write(6,19) 'width W->cs  = ', w_w_cs, 
     &               ' GeV,   BR(W->cs)  = ', w_w_cs/wwidth
         write(6,*)
      end if
      if (.not. lhdecay) then
         write(6,15) 'Higgs total width =  ', SMWDTH, ' GeV'
         write(6,*)
      else
         write(6,*) 'HDECAY output for Higgs sector:'
         write(6,*)
         write(6,15) 'Higgs total width =  ', SMWDTH, ' GeV'
         write(6,*)
         write(6,16) 'BR (H -> bb~  ) =  ', SMBRB
         write(6,16) 'BR (H -> taus ) =  ', SMBRL
         write(6,16) 'BR (H -> muons) =  ', SMBRM
         write(6,16) 'BR (H -> ss~  ) =  ', SMBRS
         write(6,16) 'BR (H -> cc~  ) =  ', SMBRC
         write(6,16) 'BR (H -> tt~  ) =  ', SMBRT
         write(6,16) 'BR (H -> gg   ) =  ', SMBRG
         write(6,16) 'BR (H -> AA   ) =  ', SMBRGA
         write(6,16) 'BR (H -> ZA   ) =  ', SMBRZGA
         write(6,16) 'BR (H -> W+W- ) =  ', SMBRW
         write(6,16) 'BR (H -> ZZ   ) =  ', SMBRZ
         write(6,*)
      endif
#ifdef HELAS_CHECK
      write(6,*) 'BRS couplings:'
      write(6,*)
      write(6,10) 'xz mass =  ',zmass
      write(6,10) 'xw mass =  ',wmass
      write(6,*)
      write(6,20) 'VHS  ZHxz   = ',xzhz
      write(6,20) 'VHS  WHxwp  = ',xwhwp
      write(6,20) 'VHS  WHxwm  = ',xwhwm
      write(6,*)
      write(6,20) 'VVS  WZxwp  = ',xwzwp
      write(6,20) 'VVS  WZxwm  = ',xwzwm
      write(6,20) 'VVS  WAxwp  = ',xwawp
      write(6,20) 'VVS  WAxwm  = ',xwawm
      write(6,*)
      write(6,20) 'VVHS WZHxwp = ',xwzwp
      write(6,20) 'VVHS WZHxwm = ',xwzwm
      write(6,20) 'VVHS WAHxwp = ',xwawp
      write(6,20) 'VVHS WAHxwm = ',xwawm
      write(6,*)
      write(6,14) 'xztop(L) =  ',xzt(1), '    xztop(R) =  ',xzt(2)
      write(6,14) 'xzbot(L) =  ',xzb(1), '    xzbot(R) =  ',xzb(2)
      write(6,14) 'xztau(L) =  ',xzl(1), '    xztau(R) =  ',xzl(2)
      write(6,*)
      write(6,14) 'xwpq(L)  =  ',xwpq(1), '    xwpq(R)  =  ',xwpq(2)
      write(6,14) 'xwmq(L)  =  ',xwmq(1), '    xwpm(R)  =  ',xwmq(2)
      write(6,14) 'xwpl(L)  =  ',xwpl(1), '    xwpl(R)  =  ',xwpl(2)
      write(6,14) 'xwml(L)  =  ',xwml(1), '    xwml(R)  =  ',xwml(2)
      write(6,*)
#endif
ccc
      return
      end

C
C-----------------------------------------------------------------------------
C
      double precision function alfa( qsq )
C
C-----------------------------------------------------------------------------
C
C	This function returns the 1-loop value of alpha.
C
C	INPUT: 
C		qsq   = Q^2
C
C-----------------------------------------------------------------------------
C
      implicit none
      double precision  qsq
#include "coupl.inc"
c
c constants
c
      double precision  One, Three, Pi
      parameter( One = 1.0d0, Three = 3.0d0 )
      parameter( Pi = 3.14159265358979323846d0 )
cc
      alfa = alpha / ( 1.0d0 - alpha*dlog( qsq/zmass**2 ) /Three /Pi )
ccc
      return
      end

C
C-----------------------------------------------------------------------------
C
      double precision function alfaw( qsq,nh )
C
C-----------------------------------------------------------------------------
C
C	This function returns the 1-loop value of alpha_w.
C
C	INPUT: 
C		qsq = Q^2
C               nh  = # of Higgs doublets
C
C-----------------------------------------------------------------------------
C
      implicit none
      double precision  qsq, alphaw, dum
      integer  nh, nq
#include "coupl.inc"
c
c constants
c
      double precision  Two, Four, Pi, Twpi
      parameter( Two = 2.0d0, Four = 4.0d0 )
      parameter( Pi = 3.14159265358979323846d0 )
      parameter( Twpi = 3.0d0*Four*Pi )
cc
      if ( qsq.ge.tmass**2 ) then
         nq = 6
      else
         nq = 5
      end if
      alphaw = gw**2 / Four / Pi
      dum = ( 22.0d0 - Four*nq - nh/Two ) / Twpi
      alfaw = alphaw / ( 1.0d0 + dum*alphaw*dlog( qsq/zmass**2 ) )
ccc
      return
      end

C-----------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION ALPHAS(Q,AMZ,NLOOP)
c
c     Evaluation of strong coupling constant alpha_S
c     Author: R.K. Ellis
c
c     q -- scale at which alpha_s is to be evaluated
c     amz -- value of alpha_s at the mass of the Z-boson
c     nloop -- the number of loops (1,2, or 3) at which beta 
c     function is evaluated to determine running.
c     the values of the cmass and the bmass should be set
c     in common block qmass.
C-----------------------------------------------------------------------------

      IMPLICIT NONE
      DOUBLE PRECISION Q,T,AMZ,AMZ0,AMB,AMC
      DOUBLE PRECISION AS_OUT
      INTEGER NLOOP,NLOOP0,NF3,NF4,NF5
      PARAMETER(NF5=5,NF4=4,NF3=3)
      double precision     hmass, wmass, zmass, amass,
     &                     tmass, bmass, lmass, cmass
      common /COUPL_MASS/  hmass, wmass, zmass, amass,
     &                     tmass, bmass, lmass, cmass
      SAVE AMZ0,NLOOP0,AMB,AMC
      DATA AMZ0,NLOOP0/0D0,0/
      IF (Q .LE. 0D0) THEN 
         WRITE(6,*) 'q .le. 0 in alphas'
         WRITE(6,*) 'q= ',Q
         PAUSE
      ENDIF
      IF (AMZ .LE. 0D0) THEN 
         WRITE(6,*) 'amz .le. 0 in alphas',AMZ
         WRITE(6,*) 'continue with amz=0.1185?'
         PAUSE
         AMZ=0.1185D0
      ENDIF
      IF (CMASS .LE. 0.3D0) THEN 
         WRITE(6,*) 'cmass .le. 0.3GeV in alphas',CMASS
         WRITE(6,*) 'continue with cmass=1.5GeV?'
         PAUSE
         CMASS=1.42D0
      ENDIF
      IF (BMASS .LE. 0D0) THEN 
         WRITE(6,*) 'bmass .le. 0 in alphas',BMASS
         WRITE(6,*) 'COMMON/QMASS/CMASS,BMASS'
         PAUSE
         BMASS=4.7D0
      ENDIF
c--- establish value of coupling at b- and c-mass and save
      IF ((AMZ .NE. AMZ0) .OR. (NLOOP .NE. NLOOP0)) THEN
         AMZ0=AMZ
         NLOOP0=NLOOP
         T=2D0*DLOG(BMASS/ZMASS)
         CALL NEWTON1(T,AMZ,AMB,NLOOP,NF5)
         T=2D0*DLOG(CMASS/BMASS)
         CALL NEWTON1(T,AMB,AMC,NLOOP,NF4)
      ENDIF

c--- evaluate strong coupling at scale q
      IF (Q  .LT. BMASS) THEN
           IF (Q  .LT. CMASS) THEN
             T=2D0*DLOG(Q/CMASS)
             CALL NEWTON1(T,AMC,AS_OUT,NLOOP,NF3)
           ELSE
             T=2D0*DLOG(Q/BMASS)
             CALL NEWTON1(T,AMB,AS_OUT,NLOOP,NF4)
           ENDIF
      ELSE
      T=2D0*DLOG(Q/ZMASS)
      CALL NEWTON1(T,AMZ,AS_OUT,NLOOP,NF5)
      ENDIF
      ALPHAS=AS_OUT
      RETURN
      END


      SUBROUTINE DIFF(Q,AMZ,NLOOP)
      IMPLICIT NONE
      DOUBLE PRECISION BETA(3:5),B0(3:5),C1(3:5),C2(3:5)
      INTEGER NLOOP,J
      DOUBLE PRECISION Q,QP,QM,AMZ
      DOUBLE PRECISION X1,X2,X3,EP,DIFF1,ALPHAS
      double precision     hmass, wmass, zmass, amass,
     &                     tmass, bmass, lmass, cmass
      common /COUPL_MASS/  hmass, wmass, zmass, amass,
     &                     tmass, bmass, lmass, cmass
C---     B0=(11.-2.*F/3.)/4./PI
      DATA B0/0.716197243913527D0,0.66314559621623D0,0.61009394851893D0/
C---     C1=(102.D0-38.D0/3.D0*F)/4.D0/PI/(11.D0-2.D0/3.D0*F)
      DATA C1/.565884242104515D0,0.49019722472304D0,0.40134724779695D0/
C---     C2=(2857.D0/2.D0-5033*F/18.D0+325*F**2/54)
C---     /16.D0/PI**2/(11.D0-2.D0/3.D0*F)
      DATA C2/0.453013579178645D0,0.30879037953664D0,0.14942733137107D0/
C---     DEL=SQRT(4*C2-C1**2)

      X1=ALPHAS(Q,AMZ,1)
      X2=ALPHAS(Q,AMZ,2)
      X3=ALPHAS(Q,AMZ,3)
      J=3
      IF (Q .GT. CMASS) J=4
      IF (Q .GT. BMASS) J=5
      EP=.001D0
      QP=Q*(1.D0+EP)
      QM=Q*(1.D0-EP)
      IF (NLOOP .EQ.1) THEN 
      BETA(J)=-B0(J)*X1**2
      DIFF1=(ALPHAS(QP,AMZ,1)-ALPHAS(QM,AMZ,1))/4d0/EP/BETA(J)
      ENDIF
      IF (NLOOP .EQ.2) THEN 
      BETA(J)=-B0(J)*X2**2*(1.D0+C1(J)*X2)
      DIFF1=(ALPHAS(QP,AMZ,2)-ALPHAS(QM,AMZ,2))/4d0/EP/BETA(J)
      ENDIF
      IF (NLOOP .EQ.3) THEN 
      BETA(J)=-B0(J)*X3**2*(1.D0+C1(J)*X3+C2(J)*X3**2)
      DIFF1=(ALPHAS(QP,AMZ,3)-ALPHAS(QM,AMZ,3))/4d0/EP/BETA(J)
      ENDIF
      WRITE(6,*) Q,DIFF1,NLOOP
      RETURN
      END

      SUBROUTINE NEWTON1(T,A_IN,A_OUT,NLOOP,NF)
C     Author: R.K. Ellis

c---  calculate a_out using nloop beta-function evolution 
c---  with nf flavours, given starting value as-in
c---  given as_in and logarithmic separation between 
c---  input scale and output scale t.
c---  Evolution is performed using Newton's method,
c---  with a precision given by tol.

      IMPLICIT NONE
      INTEGER NLOOP,NF
      REAL*8 T,A_IN,A_OUT,AS,TOL,F2,F3,F,FP,DELTA
      REAL*8 B0(3:5),C1(3:5),C2(3:5),DEL(3:5)
      PARAMETER(TOL=5.D-4)
C---     B0=(11.-2.*NF/3.)/4./PI
      DATA B0/0.716197243913527D0,0.66314559621623D0,0.61009394851893D0/
C---     C1=(102.D0-38.D0/3.D0*NF)/4.D0/PI/(11.D0-2.D0/3.D0*NF)
      DATA C1/.565884242104515D0,0.49019722472304D0,0.40134724779695D0/
C---     C2=(2857.D0/2.D0-5033*NF/18.D0+325*NF**2/54)
C---     /16.D0/PI**2/(11.D0-2.D0/3.D0*NF)
      DATA C2/0.453013579178645D0,0.30879037953664D0,0.14942733137107D0/
C---     DEL=SQRT(4*C2-C1**2)
      DATA DEL/1.22140465909230D0,0.99743079911360D0,0.66077962451190D0/
      F2(AS)=1D0/AS+C1(NF)*LOG((C1(NF)*AS)/(1D0+C1(NF)*AS))
      F3(AS)=1D0/AS+0.5D0*C1(NF)
     & *LOG((C2(NF)*AS**2)/(1D0+C1(NF)*AS+C2(NF)*AS**2))
     & -(C1(NF)**2-2D0*C2(NF))/DEL(NF)
     & *ATAN((2D0*C2(NF)*AS+C1(NF))/DEL(NF))

           
      A_OUT=A_IN/(1D0+A_IN*B0(NF)*T)
      IF (NLOOP .EQ. 1) RETURN
      A_OUT=A_IN/(1D0+B0(NF)*A_IN*T+C1(NF)*A_IN*LOG(1D0+A_IN*B0(NF)*T))
      IF (A_OUT .LT. 0D0) AS=0.3D0
 30   AS=A_OUT

      IF (NLOOP .EQ. 2) THEN
      F=B0(NF)*T+F2(A_IN)-F2(AS)
      FP=1D0/(AS**2*(1D0+C1(NF)*AS))
      ENDIF
      IF (NLOOP .EQ. 3) THEN
      F=B0(NF)*T+F3(A_IN)-F3(AS)
      FP=1D0/(AS**2*(1D0+C1(NF)*AS+C2(NF)*AS**2))
      ENDIF
      A_OUT=AS-F/FP
      DELTA=ABS(F/FP/AS)
      IF (DELTA .GT. TOL) GO TO 30
      RETURN
      END


C-----------------------------------------------------------------------------
C
      double precision function mfrun(mf,scale,asmz,nloop)
C
C-----------------------------------------------------------------------------
C
C	This function returns the 2-loop value of a MSbar fermion mass
C       at a given scale.
C
C	INPUT: mf    = MSbar mass of fermion at MSbar fermion mass scale 
C	       scale = scale at which the running mass is evaluated
C	       asmz  = AS(MZ) : this is passed to alphas(scale,asmz,nloop)
C              nloop = # of loops in the evolution
C       
C
C
C	EXTERNAL:      double precision alphas(scale,asmz,nloop)
C                      
C-----------------------------------------------------------------------------
C
      implicit none
C
C     ARGUMENTS
C
      double precision  mf,scale,asmz
      integer           nloop
C
C     LOCAL
C
      double precision  beta0, beta1,gamma0,gamma1
      double precision  A1,as,asmf,l2
      integer  nf
C
C     EXTERNAL
C
      double precision  alphas
      external          alphas
C
C     COMMON
C
      double precision     hmass, wmass, zmass, amass,
     &                     tmass, bmass, lmass, cmass
      common /COUPL_MASS/  hmass, wmass, zmass, amass,
     &                     tmass, bmass, lmass, cmass
c
c     CONSTANTS
c
      double precision  One, Two, Three, Pi
      parameter( One = 1.0d0, Two = 2.0d0, Three = 3.0d0 )
      parameter( Pi = 3.14159265358979323846d0) 
cc
C
C
      if ( mf.gt.tmass ) then
         nf = 6
      else
         nf = 5
      end if

      beta0 = ( 11.0d0 - Two/Three *nf )/4d0
      beta1 = ( 102d0  - 38d0/Three*nf )/16d0
      gamma0= 1d0
      gamma1= ( 202d0/3d0  - 20d0/9d0*nf )/16d0
      A1    = -beta1*gamma0/beta0**2+gamma1/beta0
      as    = alphas(scale,asmz,nloop)
      asmf  = alphas(mf   ,asmz,nloop)
      l2    = (1+ A1*as/Pi)/(1+ A1*asmf/Pi)
      
      
      mfrun = mf * (as/asmf)**(gamma0/beta0)

      if(nloop.eq.2) mfrun =mfrun*l2
ccc
      return
      end


