C========================================================================
C
C QQHQQ calculates the matrix elements**2 for light Higgs production by
C electroweak boson fusion in quark quark scattering
C
C    q1 q3 ----> q2 q4 H;   H ---> b5-bar b6 optional
C
C QQHQQI must be called first to initialize some couplings
C
C========================================================================
C
      subroutine qqhqqi( bosdec)  ! >0  includes matrix elements and coupl.
          use globalvars, only: ldoblha
      implicit none                !        for H -> xx (based on bbar)
      integer bosdec               !  0  means the driver does a *|M_dec|^2

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      double precision xmmu
c
c the main task of QQHQQI is to set up the products of coupling constants
c needed in Feynman graphs a, ..., g which are stored in
c
      double precision  fcpl(4,6)
* counterterm couplings for qqV
      double complex fcplCT(4,6)
      double complex fcpla(4,6,2) !hzz,hww
      double complex fcpla1(4,4,2),fcpla2(4,4,2),fcpla3(4,4,2) !hgg,hgz

** fermion qqV couplings
      double complex fcplCTfm(4,6)
      double complex clrCTfm(3:4,2:4,-1:1)
      common /FMclr/ clrCTfm, fcplCTfm
c
c and which are calculated from output of KOPPLN

** clrCT is the counterterm for qqV
      double complex clrCT(3:4,2:4,-1:1)
      common /bkopouCT/ clrCT
      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)

#include "VBFNLO/utilities/koppln.inc"

*sophy: anomalous couplings (resulting from electroweak corrections)
      double complex ahvv(3,4,4), ahvvL(3,4,4)
      common/tensorhvv/ ahvv, ahvvL
c
c here fcpl(sig,i) contains the coupling constant factors. 
c     sig = 1,..,4    identifies the 4 helicity combinations of the quarks
c     i   = 1,..,6    numbers the possible subprocesses (4 NC and 2 CC)
c  
c the various processes are identified by the following
c
      integer  v2ini(6), v3ini(6)
      integer v2ini1(4), v3ini1(4) ! V2 = Z V3 = gamma 
      integer v2ini2(4), v3ini2(4) ! V2 = gamma V3 = Z
      integer v2ini3(4), v3ini3(4) ! V2 = gamma V3 = gamma

      integer  fl1(6), fl3(6)
      integer  sig1(4,2),sig3(4,2)

      data fl1 /3,3,4,4,3,4/
      data fl3 /3,4,3,4,4,3/

      save fl1,fl3
*sophy: added fcplCT to common block cqqhqq
      common /cqqhqq/ fcplCT, fcpl, sig1,sig3, v2ini,v3ini
      common /cqqhqqa/ fcpla,fcpla1,fcpla2,fcpla3,
     $      v2ini1,v3ini1, v2ini2,v3ini2,v2ini3,v3ini3
      data v2ini /2,2,2,2,3,4/, v3ini /2,2,2,2,4,3/
      data v2ini1 /2,2,2,2/, v3ini1 /1,1,1,1/
      data v2ini2 /1,1,1,1/, v3ini2 /2,2,2,2/
      data v2ini3 /1,1,1,1/, v3ini3 /1,1,1,1/
      data sig1 /-1,-1, 1, 1,-1,3*0/
      data sig3 /-1, 1,-1, 1,-1,3*0/

c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"

      real*8 c2,c2o2pi
      parameter (c2=4d0/3d0, c2o2pi=c2/2d0/pi)
      double precision lnmusqp1p2, lnmusqp3p4

c other local variables for QQBQQI
c
      integer  i,k,mu,v2,v3, isig,isig1,isig3
c
c variables for the main part of the program
c
      double precision  pbar(0:3,4+max_v)
      double precision uucc,uuss,ddcc,ddss,udsc,ducs,res(6)
* dummy variable to store 2*Re | M_{born} M^{*}_{loop} |
      double precision loop(6)
      double precision tree(6)
      double precision  p(0:3,4+max_v), p21(0:4), p43(0:4), p65(0:4),
     1                  p87(0:4), ph(0:4)
      integer  sign(4+max_v), nlo
      double complex  psi(2,-1:1,4), j21(0:3,-1:1), j43(0:3,-1:1)
      double complex  j21j43(-1:1,-1:1), zmm(2)
      double complex mma(6,4,2)
      double complex  prop21(1:4), prop43(1:4)
      double precision  lamb, fac
* split 'fac' into 'fac' (prop.) and 'QCDfac' (qcd nlo)
      double precision QCDfac

** FERMION contn to anom matrix element
      double complex mmaFM(6,4,2)

** Flag determining whether ew corns should be included
      logical run_ew

c CPS scheme
      double precision qgammaq
      common /VBFNLO_HIGGS_CPS/ qgammaq

#include "VBFNLO/utilities/mssm.inc"

c
c reset the coupling factors
c
      do i = 1,6
         do isig = 1,4
            fcpl(isig,i) = 0.d0
            fcplCT(isig,i) = 0D0
         enddo
      enddo

      xmmu = MM
      
c
c determine the Yukawa coupling Higgs bb from BR(H--->bb)
c division by 3 takes into account the color factor 3 for H--->qq
c      
      if (bosdec.eq.3) then ! for W+W-jj
            lamb = 2*xmw*clr(1,3,-1)**2*b(6,3,4)!*treefacW
      elseif (bosdec.ne.0) then
         print*," bosdec = ",bosdec," not implemented in qqhqqi "
         stop
      end if

      do i = 1,6
         do isig = 1,4
            isig1 = sig1(isig,(i+3)/4)
            if ( isig1.ne.0 ) then
               isig3 = sig3(isig,(i+3)/4)
               v2 = v2ini(i)
               v3 = v3ini(i)

               fcpl(isig,i) = clr(fl1(i),v2,isig1)*
     &              clr(fl3(i),v3,isig3)*ahvv(1,v2,v3)!b(6,v2,v3)*xmw


** Adding in calculation of counterterm couplings
               if (ewcor_switch .and. (sector .lt. 3)) then
                  fcplCT(isig,i) = (clrCT(fl1(i),v2,isig1)*
     &                 clr(fl3(i),v3,isig3) + clr(fl1(i),v2,isig1)*
     &                 clrCT(fl3(i),v3,isig3))

** Multiplying qqV effective vertex by the appropriate factors
                  fcplCT(isig,i) = fcplCT(isig,i)*ahvv(1,v2,v3)!b(6,v2,v3)*xmw
                  if ((propLEVEL .eq. 1) .and. (model .eq. 2)) then
                     if (HiggsType .eq. 1) then
                        fcplCT(isig,i) = fcplCT(isig,i)*SBA/
     -                       (Hgs(1,1)*SBA + Hgs(1,2)*CBA)
                     else if (HiggsType .eq. 2) then
                        fcplCT(isig,i) = fcplCT(isig,i)*CBA/
     -                       (Hgs(2,1)*SBA + Hgs(2,2)*CBA)
                     else if (HiggsType .eq. 3) then
                        fcplCT(isig,i) = (0D0,0D0)
                     end if
                  end if
               end if


** qqV ct vertex for fermions
** Adding in calculation of counterterm couplings
               if (ewcor_switch) then
                  fcplCTfm(isig,i) = (clrCTfm(fl1(i),v2,isig1)*
     &                 clr(fl3(i),v3,isig3) + clr(fl1(i),v2,isig1)*
     &                 clrCTfm(fl3(i),v3,isig3))
** Multiplying qqV effective vertex by the appropriate factors
                  fcplCTfm(isig,i) = fcplCTfm(isig,i)*ahvv(1,v2,v3)!b(6,v2,v3)*xmw
                  if ((propLEVEL .eq. 1) .and. (model .eq. 2)) then
                     if (HiggsType .eq. 1) then
                        fcplCTfm(isig,i) = fcplCTfm(isig,i)*SBA/
     -                       (Hgs(1,1)*SBA + Hgs(1,2)*CBA)
                     else if (HiggsType .eq. 2) then
                        fcplCTfm(isig,i) = fcplCTfm(isig,i)*CBA/
     -                       (Hgs(2,1)*SBA + Hgs(2,2)*CBA)
                     else if (HiggsType .eq. 3) then
                        fcplCTfm(isig,i) = (0D0,0D0)
                     end if
                  end if
               end if
** end of qqV ct vertex for fermions

               if (bosdec .ne. 0) then
                  fcpl(isig,i) = fcpl(isig,i)*lamb
               end if
*               write(*,*)'fl1(i),v2,isig1 =', fl1(i),v2,isig1
*               write(*,*)'clr(fl1(i),v2,isig1) =', clr(fl1(i),v2,isig1)
*               write(*,*)'clrCT(fl1(i),v2,isig1) =', clrCT(fl1(i),v2,isig1)
            endif
         enddo
      enddo


      if (with_anomHiggs) call qqhqqai(0)
      return
c
c entry for H production
c
      entry qqHqq(pbar,sign,nlo,uucc,uuss,ddcc,ddss,udsc,ducs,
     -     tree)

      run_ew = .false.
      if (ewcor_switch) then
         if (nlo .eq. 2) then
            run_ew = .true.
         end if
      end if

c     
c define the internal momenta
c
      do mu = 0,3
         do i = 1,4+max(n_v,2)
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
         p21(mu) = p(mu,2) - p(mu,1)
         p43(mu) = p(mu,4) - p(mu,3)
         p65(mu) = p(mu,6) - p(mu,5)
         if (n_v.eq.4) then
            p87(mu) = p(mu,8) - p(mu,7)
            ph(mu) = p65(mu) + p87(mu)
         elseif ((n_v.eq.2).or.(n_v.eq.1)) then
            ph(mu) = p65(mu)
         endif
      enddo
      p21(4) = p21(0)**2 - p21(1)**2 - p21(2)**2 - p21(3)**2
      p43(4) = p43(0)**2 - p43(1)**2 - p43(2)**2 - p43(3)**2
      p65(4) = p65(0)**2 - p65(1)**2 - p65(2)**2 - p65(3)**2
      if (n_v.eq.4) p87(4) = 
     1         p87(0)**2 - p87(1)**2 - p87(2)**2 - p87(3)**2
      ph(4) = ph(0)**2 - ph(1)**2 - ph(2)**2 - ph(3)**2
c
c get the vector boson propagator factors
c
      if (p21(4).le.0d0) then
         prop21(2) = 1/(p21(4)-xm2(2))
         prop21(3) = 1/(p21(4)-xm2(3))
      else
         prop21(2) = 1/dcmplx(p21(4)-xm2(2),xmg(2))
         prop21(3) = 1/dcmplx(p21(4)-xm2(3),xmg(3))
      endif
      prop21(4) = prop21(3)
      prop21(1) = 1.0d0/dcmplx(p21(4),0.0d0)

      if (p43(4).le.0d0) then
         prop43(2) = 1/(p43(4)-xm2(2))
         prop43(3) = 1/(p43(4)-xm2(3))
      else
         prop43(2) = 1/dcmplx(p43(4)-xm2(2),xmg(2))
         prop43(3) = 1/dcmplx(p43(4)-xm2(3),xmg(3))
      endif
      prop43(4) = prop43(3)
      prop43(1) = 1.0d0/dcmplx(p43(4),0.0d0)
c
c get the external spinors (including factor sqrt(2E) )
c
      call psi0m(4,pbar(0,1),sign(1),psi)
c
c get the f-fbar currents J21^mu, J43^mu
c
      call curr(1,psi(1,-1,2),psi(1,-1,1),j21)
      call curr(1,psi(1,-1,4),psi(1,-1,3),j43)
c
c get the dot products of the currents for the 4 helicity combinations
c
      do isig1 = -1,1,2
         do isig3 = -1,1,2
            j21j43(isig1,isig3) = j21(0,isig1)*j43(0,isig3) -
     &                            j21(1,isig1)*j43(1,isig3) -
     &                            j21(2,isig1)*j43(2,isig3) -
     &                            j21(3,isig1)*j43(3,isig3)
         enddo
      enddo
      
      if (n_v.gt.1) then
c
c now get the coupling*propagator factors for subprocess k, helicity
c combination isig and sum the contributions from the 4 helicity 
c combinations for all subprocesses. Multiply by the factor corresponding to
c    SUM_pol |D_H(p_6-p_5) * psi-bar(6) psi(5)|^2 * color factor for H->bb
c
        fac = 1.d0 / ( (ph(4)-xm2(6))**2 + xmg(6)**2 ) !propagator^2 only
c
c multiply by 16 pi*M Gamma to compensate for the phase space integral in
c LIPSN
c
c   Int dq^2/2pi * dOmega/(32 pi^2) * 1/(q^2-M*2)^2+(MGamma)^2) 
c
c   = 1/(16 pi MGamma)
c
c (where the propagator factor is fac above)
        fac = fac *16*pi*qgammaq
      else 
c Higgs is final-state - no decays (real or simulated for technical reasons)
        fac = 1d0
      endif

*sophy: took out 'colour factor' and turned into QCDfac
      if (doVirtuals .and. nlo .eq. 1) then
         QCDfac = ( 1d0+(als(1,1)+als(2,1))*c2o2pi*cvirt )
      else
         QCDfac = 1D0           
      endif

      if (nlo .eq. 99) then
         write(*,*)'  '
         write(*,*)'ph(4) =', ph(4)
         write(*,*)'xm2(6) =', xm2(6)
         write(*,*)'xmg(6) =', xmg(6)
         write(*,*)'pi =', pi
         write(*,*)'  '
      end if

*sophy: call to anomalous formfactors caused by ew loop corrections
#ifdef WITH_LT
      if (run_ew) then
         if ((sign(1) .gt. 0) .and. (sign(3) .gt. 0)) then
            call form_factors(p21, p43, ahvvL, 1)
         end if
      end if
#endif

      if (with_anomHiggs .or. run_ew) then
         call qqhqqa(j21,j43,p21,p43,prop21,prop43, !input
     1                       mma,0)  !output
         call qqhqqa(j21,j43,p21,p43,prop21,prop43, !input
     1                       mmaFM,1)  !output
c$$$         write(*,*)'mma =', mma
c$$$         write(*,*)'mma fermions =', mmaFM
      end if


** calling effective qqV coupling.  The final entry is for debugging.
      if ((sector .ge. 3) .and. run_ew) then
#ifdef WITH_LT
         call qqV_eff(pbar,sign,fcplCT,0)
#endif
c         write(*,*)'' 
c         write(*,*)'NEW QUARK PATTERN'
c         do i = 1,6
c            write(*,*)'  '
c            do isig = 1,4
c               write(*,*)'i, isig, fcplCT(isig,i) =', i, isig, fcplCT(isig,i)
c            enddo
c         enddo
c         if ((sign(1) .lt. 0) .and. (sign(3) .lt. 0)) stop

         do i = 1,6
            do isig = 1,4
               isig1 = sig1(isig,(i+3)/4)
               if ( isig1.ne.0 ) then
                  v2 = v2ini(i)
                  v3 = v3ini(i)

** Multiplying qqV effective vertex by the appropriate factors
                  fcplCT(isig,i) = fcplCT(isig,i)*ahvv(1,v2,v3) !b(6,v2,v3)*xmw
                  if ((propLEVEL .eq. 1) .and. (model .eq. 2)) then
                     if (HiggsType .eq. 1) then
                        fcplCT(isig,i) = fcplCT(isig,i)*SBA/
     -                       (Hgs(1,1)*SBA + Hgs(1,2)*CBA)
                     else if (HiggsType .eq. 2) then
                        fcplCT(isig,i) = fcplCT(isig,i)*CBA/
     -                       (Hgs(2,1)*SBA + Hgs(2,2)*CBA)
                     else if (HiggsType .eq. 3) then
                        fcplCT(isig,i) = (0D0,0D0)
                     end if
                  end if

               end if
            enddo
         enddo

      end if


      i=1
      do k = 1,4
* clearing values of res and loop
         res(k) = 0d0
         loop(k) = 0d0
         tree(k) = 0D0
         do isig = 1,4
            isig1 = sig1(isig,1)
            isig3 = sig3(isig,1)
            zmm(i) = prop21(2)*prop43(2)*fcpl(isig,k)*
     1               j21j43(isig1,isig3)
            if (with_anomHiggs) then
               zmm(i) = zmm(i) + mma(k,isig,i)
            endif
            res(k) = res(k) + dreal(zmm(i))**2+dimag(zmm(i))**2

* calculating 2 |M_{loop} M^{*}_{born}| 
            if (run_ew) then

               mma(k,isig,i) = mma(k,isig,i) + prop21(2)*prop43(2)*
     -               fcplCT(isig,k)*j21j43(isig1,isig3)

               loop(k) = loop(k) + 
     -              2D0*DREAL(zmm(i)*conjg(mma(k,isig,i)))

*including |M_{loop}|^2 from fermions
               if (loopSqr_switch .and. run_ew) then
                  mmaFM(k,isig,i) = mmaFM(k,isig,i) + 
     -                 prop21(2)*prop43(2)*
     -                 fcplCTfm(isig,k)*j21j43(isig1,isig3)
                  loop(k) = loop(k) + (dreal(mmaFM(k,isig,i))**2 + 
     -                 dimag(mmaFM(k,isig,i))**2) 
               end if
            end if
         enddo

*sophy: putting colour factor here NB: NOT for ew corns!  Adding 2|M_b M_l|
         tree(k) = res(k) * 9D0 * fac

** If we're running EW with BW distn, a separate result containing ONLY the
** ew corn (i.e. 2*Re{born*loop}) is output
         if (ewcor_switch .and. (nlo .eq. 2)) then
            res(k) = loop(k)
         else
            res(k) = res(k) * QCDfac + loop(k)
         end if

         res(k) = res(k) * 9D0 *fac ! 9 is color factor for color summed |M|^2
      enddo


      do k = 5, 6
*sophy: clearing values of res and loop
         res(k) = 0d0
         loop(k) = 0d0
         tree(k) = 0D0
         zmm(i) = fcpl(1,k)*prop21(3)*prop43(3)*j21j43(-1,-1)
         if (with_anomHiggs) then
            zmm(i) = zmm(i) + mma(k,1,i)
         endif
         res(k) = (dreal(zmm(i))**2+dimag(zmm(i))**2) 

*sophy: calculating 2 |M_{loop} M^{*}_{born}| 
         if (run_ew) then

            mma(k,1,i) = mma(k,1,i) + prop21(3)*prop43(3)*
     -           (fcplCT(1,k)*j21j43(-1,-1) + fcplCT(2,k)*j21j43(-1,1)+
     -           fcplCT(3,k)*j21j43(1,-1) + fcplCT(4,k)*j21j43(1,1))
*sophy: (note: other helicity combinations are included here, so that an 
*        effective qqV coupling can be used properly)
            loop(k) = loop(k) + 2D0*DREAL(zmm(i)*conjg(mma(k,1,i)))

*including |M_{loop}|^2 from fermions
            if (loopSqr_switch .and. run_ew) then
               mmaFM(k,1,i) = mmaFM(k,1,i) + prop21(3)*prop43(3)*
     -              (fcplCTfm(1,k)*j21j43(-1,-1) + 
     -              fcplCTfm(2,k)*j21j43(-1,1)+
     -              fcplCTfm(3,k)*j21j43(1,-1) + fcplCTfm(4,k)*
     -              j21j43(1,1))
               loop(k) = loop(k) + (dreal(mmaFM(k,1,i))**2 + 
     -              dimag(mmaFM(k,1,i))**2) 
            end if

         end if

* putting colour factor here NB: NOT for ew corns!  Ading 2|M_b M_l|
         tree(k) = res(k)* 9D0 * fac
** If we're running EW with BW distn, a separate result containing ONLY the
** ew corn (i.e. 2*Re{born*loop}) is output
         if (ewcor_switch .and. (nlo .eq. 2)) then
            res(k) = loop(k)
         else
            res(k) = res(k) * QCDfac + loop(k)
         end if
         res(k) = res(k) * 9D0 *fac ! 9 is color factor for color summed |M|^2
      enddo

c
c assign the final results for the subprocesses
c
      uucc = res(1) 
      uuss = res(2)
      ddcc = res(3)
      ddss = res(4)
      udsc = res(5)
      ducs = res(6)

      if (ldoblha) then
        if (blha_amptype(blha_curproc).eq.0) then ! tree
          blha_amp(1)=blha_amp(1)+
     &      res(blha_idsubproc(blha_cursubproc,blha_curproc))

        else if (blha_amptype(blha_curproc).eq.1) then ! loop
          lnmusqp1p2 = log(blha_scale**2/(-p21(4)))
          lnmusqp3p4 = log(blha_scale**2/(-p43(4)))

          blha_amp(1)=blha_amp(1)+
     &      blha_CF/(2d0*pi)*2d0*
     &       (als(1,1)+als(2,1))*(-1d0)*
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))
          blha_amp(2)=blha_amp(2)+
     &      1d0/(2d0*pi)*2d0*
     &       ( als(1,1)*(-blha_gammaQuark-blha_CF*lnmusqp1p2)
     &        +als(2,1)*(-blha_gammaQuark-blha_CF*lnmusqp3p4) )*
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))
c this is more complicated than necessary, 
c but should work for other VBF processes as well
          blha_amp(3)= blha_amp(3)+
     &      (res(blha_idsubproc(blha_cursubproc,blha_curproc))     ! tree + loop + I
     &       -tree(blha_idsubproc(blha_cursubproc,blha_curproc)))  ! -tree
     &       * blha_CF/c2                                          ! correct for Nc ?= inf
     &      + (-(als(1,1)+als(2,1))/2d0/pi*                        ! -I
     &          (blha_CF*cvirt+
     &            2d0*(blha_KQuark+blha_gammaQuark+
     &                 blha_tgammaQuark+(pi**2/6d0-1)*blha_CF) )
     &         +als(1,1)/pi*                                       ! +(mu^2/2p1.p2)^epsilon
     &         (-blha_gammaQuark*lnmusqp1p2-blha_CF/2d0*lnmusqp1p2**2)
     &         +als(2,1)/pi*                                       ! +(mu^2/2p3.p4)^epsilon
     &         (-blha_gammaQuark*lnmusqp3p4-blha_CF/2d0*lnmusqp3p4**2)
     &        )*tree(blha_idsubproc(blha_cursubproc,blha_curproc))
          blha_amp(4)= blha_amp(4)+
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))

        else if (blha_amptype(blha_curproc).eq.2) then ! cctree
          ! 1-2
          call BLHA_cctree(1,2,
     &      -blha_CF*
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))
     &    )
          ! 3-4
          call BLHA_cctree(3,4,
     &      -blha_CF*
     &       tree(blha_idsubproc(blha_cursubproc,blha_curproc))
     &    )

        else 
          call BLHA_amptypeerror(blha_amptype(blha_curproc),__FILE__,__LINE__)
        endif
      endif

** sophy: used for debugging
      if (nlo .eq. 99) then
         write(*,*)'  '
         write(*,*)'pre boxes ..'
*         write(*,*)'prop21(3) =', prop21(3)
*         write(*,*)'prop43(3) =', prop43(3)
*         write(*,*)'j21j43(-1,-1) =', j21j43(-1,-1)
*         write(*,*)'fac =', fac
*         write(*,*)'QCDfac =', QCDfac
*         do k = 1, 6
*            write(*,*)'k = ', k
*            write(*,*)'tree() =', tree(k)
*            write(*,*)'loop(k) = ', loop(k)*fac*9D0
*            write(*,*)'res(k) =', res(k)
*         end do
*         do isig = 1,4
*            do k = 1, 6
*               write(*,*)'isig, k =', isig, k
*               write(*,*)'fcpl(isig,k) =', fcpl(isig,k)
*               write(*,*)'fcplCT(isig,k) =', fcplCT(isig,k)
*            end do
*         end do               
         write(*,*)'uucc =', uucc
         write(*,*)'uuss =', uuss
         write(*,*)'ddcc =', ddcc
         write(*,*)'ddss =', ddss
         write(*,*)'udsc =', udsc
         write(*,*)'ducs =', ducs
*         stop
      end if


* calling electroweak boxes and pentagons.  Final entry is debugging label
#ifdef WITH_LT
      if ((sector .ge. 3) .and. run_ew) then
         call boxes(pbar,sign,uucc,uuss,ddcc,ddss,udsc,ducs,
     &           tree,fac,0) 
      end if
#endif


** sophy: used for debugging
      if (nlo .eq. 99) then
         write(*,*)'  '
         write(*,*)'qqhqq, post boxes ..'
         write(*,*)'uucc =', uucc
         write(*,*)'uuss =', uuss
         write(*,*)'ddcc =', ddcc
         write(*,*)'ddss =', ddss
         write(*,*)'udsc =', udsc
         write(*,*)'ducs =', ducs
*         stop
      end if
     

      return
      end


C
C********************  qqhqqg  ***********************************
C
      subroutine qqhqqgi( bosdec ) 
      implicit none                
      integer bosdec             

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/koppln.inc"
#include "VBFNLO/utilities/mssm.inc"

      double precision  fpi
      parameter ( fpi=4d0*pi)
C
C  QQHQQ calculates the matrix elements**2 for light Higgs production by
C  electroweak boson fusion in quark quark scattering
C
C        q1 q3    ---->   q2 q4 g H
C
C  QQHQQI must be called first to initialize some couplings
C
C  The main task of QQHQQI is to set up the products of coupling constants
C  needed in Feynman graphs a, ..., g which are stored in

      double precision  fcpl(4,6) ! hzz,hww
*sophy: counterterm couplings for qqV
      double complex fcplCT(4,6)
      double complex fcpla(4,6,2) !hzz,hww
      double complex fcpla1(4,4,2),fcpla2(4,4,2),fcpla3(4,4,2) !hgg,hgz

** fermion qqV couplings
      double complex fcplCTfm(4,6)
      double complex clrCTfm(3:4,2:4,-1:1)
      common /FMclr/ clrCTfm, fcplCTfm

C  and which are calculated from output of KOPPLN
** clrCT is the counterterm for qqV
      double complex clrCT(3:4,2:4,-1:1)
      common /bkopouCT/ clrCT
      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)
      double complex ahvv(3,4,4), ahvvL(3,4,4)
      common/tensorhvv/ ahvv, ahvvL
C
C  Here fcpl(sig,i) contains the coupling constant factors. 
C     sig = 1,..,4    identifies the 4 different helicity combinations of the
C                     quarks
C     i   = 1,..,6    numbers the possible subprocesses (4 NC and 2 CC)
c
C  
C  The various processes are identified by the following
      integer  v2ini(6), v3ini(6)
      integer v2ini1(4), v3ini1(4) ! V2 = Z V3 = gamma 
      integer v2ini2(4), v3ini2(4) ! V2 = gamma V3 = Z
      integer v2ini3(4), v3ini3(4) ! V2 = gamma V3 = gamma

      integer  fl1(6), fl3(6)
      integer  sig1(4,2),sig3(4,2)

      data fl1 /3,3,4,4,3,4/
      data fl3 /3,4,3,4,4,3/

      save fl1,fl3
*sophy: added fcplCT to common block cqqhqq
      common /cqqhqq/ fcplCT, fcpl, sig1,sig3, v2ini,v3ini
      common /cqqhqqa/ fcpla,fcpla1,fcpla2,fcpla3,
     $      v2ini1,v3ini1, v2ini2,v3ini2,v2ini3,v3ini3

c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"

C  Other local variables for QQBQQI
      integer i,k,mu,v2,v3, isig,isig1,isig3
C
C  Variables for the main part of the program
      real*8 pbar(0:3,4+max_v),qbar(0:4),uucc,uuss,ddcc,ddss,udsc,ducs,
     1       p(0:3,4+max_v),q(0:4),p21(0:4),p43(0:4),p65(0:4),eps(0:3,2),
     2       p21g(0:4),p43g(0:4), res(6,2), l21, l43,
     3       p87(0:4), ph(0:4)
      real*8 uucc_c(3),uuss_c(3),ddcc_c(3),ddss_c(3),udsc_c(3),ducs_c(3)
      integer sign(4+max_v), gsign
      complex*16 psi(2,-1:1,4), braket(2,-1:1,4,2), 
     1        j21(0:3,-1:1), j43(0:3,-1:1), j21e43(-1:1,-1:1,2),
     2        jh1(0:3,-1:1), jh2(0:3,-1:1), e21(0:3,-1:1,2),
     3        e43(0:3,-1:1,2), e21j43(-1:1,-1:1,2),
     4        mm21(6,4,2), mm43(6,4,2),mm21a(6,4,2),mm43a(6,4,2)
      real*8 pk(0:4,4)
      complex*16 prop21(1:4), prop43(1:4), prop21g(1:4), prop43g(1:4)
      real*8 lamb, fac
      logical lcol

      
      logical upline, lowline, initial
      common/uplow/ upline, lowline, initial      
      
c CPS scheme
      double precision qgammaq
      common /VBFNLO_HIGGS_CPS/ qgammaq

C  Reset the coupling factors
      do i = 1,6
         do isig = 1,4
            fcpl(isig,i) = 0.d0
            fcplCT(isig,i) = 0D0
         enddo
      enddo


      if (bosdec.eq.3) then   !   for w+w-jj
            lamb = 2*xmw*clr(1,3,-1)**2*b(6,3,4)!*treefacW
      elseif (bosdec.ne.0) then
         print*," bosdec = ",bosdec," not implemented in qqhqqi "
         stop
      end if

      do i = 1,6
         do isig = 1,4
            isig1 = sig1(isig,(i+3)/4)
            if ( isig1.ne.0 ) then
               isig3 = sig3(isig,(i+3)/4)
               v2 = v2ini(i)
               v3 = v3ini(i)
               if ( bosdec.ne.0 ) then
                  fcpl(isig,i) = clr(fl1(i),v2,isig1)*
     &                           clr(fl3(i),v3,isig3)*
!     &                           b(6,v2,v3)*xmw*lamb
     &                           ahvv(1,v2,v3)*lamb
               else
                  fcpl(isig,i) = clr(fl1(i),v2,isig1)*
     &                           clr(fl3(i),v3,isig3)*
     &                           ahvv(1,v2,v3)
!     &                           b(6,v2,v3)*xmw
*sophy: calculating counterterm couplings
                  if (ewcor_switch) then
                     fcplCT(isig,i) = (clrCT(fl1(i),v2,isig1)*
     &                    clr(fl3(i),v3,isig3) + clr(fl1(i),v2,isig1)*
     &                    clrCT(fl3(i),v3,isig3))*ahvv(1,v2,v3)!b(6,v2,v3)*xmw
* FERMION counterterm couplings                     
                     fcplCTfm(isig,i) = (clrCTfm(fl1(i),v2,isig1)*
     &                    clr(fl3(i),v3,isig3) + clr(fl1(i),v2,isig1)*
     &                    clrCTfm(fl3(i),v3,isig3))*ahvv(1,v2,v3)!b(6,v2,v3)*xmw
                  end if
               end if

            endif
         enddo
      enddo


      if(with_anomHiggs) call qqhqqai(0)

      return
              
C  Entry for H production
      entry qqHqqj(pbar,sign,qbar,gsign,
     1             uucc,uuss,ddcc,ddss,udsc,ducs)
      lcol = .false.
      goto 100

      entry qqHqqj_c(pbar,sign,qbar,gsign,
     1               uucc_c,uuss_c,ddcc_c,ddss_c,udsc_c,ducs_c)
      lcol = .true.
 100  continue
C  Define the internal momenta
      do mu = 0,3
         do i = 1,4+n_v
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
         q(mu) = gsign*qbar(mu)
      enddo
      q(4) = 0
      
C
C  Get the external spinors (including factor sqrt(2E) )
      call psi0m(4,pbar(0,1),sign(1),psi)

C  Get the gluon polarization vector and the gluon emission spinors
      do i = 1,2
         call polvec(qbar,i,eps(0,i))
         do isig = -1,1,2
            call ket2r(psi(1,isig,1),.true.,p(0,1),isig,q,eps(0,i),
     1           braket(1,isig,1,i),pk(0,1))
            call bra2r(psi(1,isig,2),.true.,p(0,2),isig,q,eps(0,i),
     1           braket(1,isig,2,i),pk(0,2))
            call ket2r(psi(1,isig,3),.true.,p(0,3),isig,q,eps(0,i),
     1           braket(1,isig,3,i),pk(0,3))
            call bra2r(psi(1,isig,4),.true.,p(0,4),isig,q,eps(0,i),
     1           braket(1,isig,4,i),pk(0,4))
         enddo
      enddo
C  get the f-fbar currents J21^mu, J43^mu, E21^mu, E43^mu, 
      call curr(1,psi(1,-1,2),psi(1,-1,1),j21)
      call curr(1,psi(1,-1,4),psi(1,-1,3),j43)
      do i = 1,2
         call curr(1,psi(1,-1,2),braket(1,-1,1,i),jh1)
         call curr(1,braket(1,-1,2,i),psi(1,-1,1),jh2)
         do mu = 0,3
            do isig = -1,1,2
               e21(mu,isig,i) = jh1(mu,isig) + jh2(mu,isig)
            enddo
         enddo

         call curr(1,psi(1,-1,4),braket(1,-1,3,i),jh1)
         call curr(1,braket(1,-1,4,i),psi(1,-1,3),jh2)
         do mu = 0,3
            do isig = -1,1,2
               e43(mu,isig,i) = jh1(mu,isig) + jh2(mu,isig)
            enddo
         enddo
      enddo
      
      do mu = 0,3
         p21(mu) = p(mu,2)-p(mu,1)
         p43(mu) = p(mu,4)-p(mu,3)
         p21g(mu) = p21(mu)+q(mu)
         p43g(mu) = p43(mu)+q(mu)
         p65(mu) = p(mu,6)-p(mu,5)
         if (n_v.eq.4) then
            p87(mu) = p(mu,8) - p(mu,7)
            ph(mu) = p65(mu) + p87(mu)
         elseif (n_v.eq.2) then
            ph(mu) = p65(mu)
         endif
      enddo
      p21(4) = p21(0)**2-p21(1)**2-p21(2)**2-p21(3)**2
      p43(4) = p43(0)**2-p43(1)**2-p43(2)**2-p43(3)**2
      p21g(4) = p21g(0)**2-p21g(1)**2-p21g(2)**2-p21g(3)**2
      p43g(4) = p43g(0)**2-p43g(1)**2-p43g(2)**2-p43g(3)**2
      p65(4) = p65(0)**2-p65(1)**2-p65(2)**2-p65(3)**2
      if (n_v.eq.4) p87(4) = 
     1         p87(0)**2 - p87(1)**2 - p87(2)**2 - p87(3)**2
      ph(4) = ph(0)**2 - ph(1)**2 - ph(2)**2 - ph(3)**2
C
C  Get the vector boson propagator factors
      if (p21(4).le.0d0) then
         prop21(2) = 1/(p21(4)-xm2(2))
         prop21(3) = 1/(p21(4)-xm2(3))
      else
         prop21(2) = 1/dcmplx(p21(4)-xm2(2),xmg(2))
         prop21(3) = 1/dcmplx(p21(4)-xm2(3),xmg(3))
      endif
      prop21(4) = prop21(3)
      prop21(1) = 1.0d0/dcmplx(p21(4),0.0d0)

      if (p43(4).le.0d0) then
         prop43(2) = 1/(p43(4)-xm2(2))
         prop43(3) = 1/(p43(4)-xm2(3))
      else
         prop43(2) = 1/dcmplx(p43(4)-xm2(2),xmg(2))
         prop43(3) = 1/dcmplx(p43(4)-xm2(3),xmg(3))
      endif
      prop43(4) = prop43(3)
      prop43(1) = 1.0d0/dcmplx(p43(4),0.0d0)

      if (p21g(4).le.0d0) then
         prop21g(2) = 1/(p21g(4)-xm2(2))
         prop21g(3) = 1/(p21g(4)-xm2(3))
      else
         prop21g(2) = 1/dcmplx(p21g(4)-xm2(2),xmg(2))
         prop21g(3) = 1/dcmplx(p21g(4)-xm2(3),xmg(3))
      endif
      prop21g(4) = prop21g(3)
      prop21g(1) = 1.0d0/dcmplx(p21g(4),0.0d0)

      if (p43g(4).le.0d0) then
         prop43g(2) = 1/(p43g(4)-xm2(2))
         prop43g(3) = 1/(p43g(4)-xm2(3))
      else
         prop43g(2) = 1/dcmplx(p43g(4)-xm2(2),xmg(2))
         prop43g(3) = 1/dcmplx(p43g(4)-xm2(3),xmg(3))
      endif
      prop43g(4) = prop43g(3)
      prop43g(1) = 1.0d0/dcmplx(p43g(4),0.0d0)

C  get the dot products of the currents for the 4 helicity combinations
      do isig1 = -1,1,2
         do isig3 = -1,1,2
            do i = 1,2
               e21j43(isig1,isig3,i) = e21(0,isig1,i)*j43(0,isig3)-
     1                                 e21(1,isig1,i)*j43(1,isig3)-
     2                                 e21(2,isig1,i)*j43(2,isig3)-
     3                                 e21(3,isig1,i)*j43(3,isig3)
               j21e43(isig1,isig3,i) = e43(0,isig3,i)*j21(0,isig1)-
     1                                 e43(1,isig3,i)*j21(1,isig1)-
     2                                 e43(2,isig3,i)*j21(2,isig1)-
     3                                 e43(3,isig3,i)*j21(3,isig1)
            enddo
         enddo
      enddo

C
C  now get the coupling*propagator factors for subprocess k, helicity
C  combination isig
      if(with_anomHiggs) call qqhqqja(j21,j43,e21,e43,p21,p43,
     1          p21g,p43g,prop21,prop43,prop21g,prop43g, !input
     2          mm21a,mm43a)                             !output

      do k = 1,6
         do isig = 1,4
            isig1 = sig1(isig,(k+3)/4)
            if (isig1.ne.0) then
               isig3 = sig3(isig,(k+3)/4)
               v2 = v2ini(k)
               v3 = v3ini(k)
               do i = 1,2
                  mm21(k,isig,i) = fcpl(isig,k)*prop21g(v2)*prop43(v3)*
     1                             e21j43(isig1,isig3,i)
                  if(with_anomHiggs) then
                     mm21(k,isig,i) = mm21(k,isig,i) +
     1                                       mm21a(k,isig,i)
                  endif
                  mm43(k,isig,i) = fcpl(isig,k)*prop21(v2)*prop43g(v3)*
     1                             j21e43(isig1,isig3,i)
                 if(with_anomHiggs) then
                     mm43(k,isig,i) = mm43(k,isig,i) +
     1                                       mm43a(k,isig,i)


                  endif
               enddo
            else
               do i = 1,2
                  mm21(k,isig,i) = 0
                  mm43(k,isig,i) = 0
               enddo
            endif
         enddo
      enddo
C
C  Now sum the contributions from the 4 helicity combinations for all 
C  subprocesses and multiply by the factor corresponding to
C
C    SUM_pol |D_H(p_6-p_5) * psi-bar(6) psi(5)|^2 * color factor
C    * g_s**2 * production color factor
      fac = 1.d0 / ( (ph(4)-xm2(6))**2 + xmg(6)**2 ) !propagator^2 only
c
c multiply by 16 pi*M Gamma to compensate for the phase space integral in
c LIPSN
c
c   Int dq^2/2pi * dOmega/(32 pi^2) * 1/(q^2-M*2)^2+(MGamma)^2) 
c
c   = 1/(16 pi MGamma)
c
c (where the propagator factor is fac above)
      fac = fac *16*pi*qgammaq
      fac = fac * 12                     ! 12 = C_2*9 is the color factor
c
c for initial gluon, i.e. gsign=-1, eliminate diagrams with V-->q qbar decay
      l21 = fpi*als(1,1)
      l43 = fpi*als(2,1)
      
      if(upline.and.(.not.initial)) l43 = 0
      if(lowline.and.(.not.initial)) l21=0

      if (gsign.eq.-1) then
         if (sign(1).eq.-sign(2) .and. 
     1       sign(3).eq.sign(4) ) then
            l43 = 0            !  initial gluon attached to 1-2 line only
         elseif (sign(1).eq.sign(2) .and. 
     1           sign(3).eq.-sign(4) ) then
            l21 = 0            !  initial gluon attached to 3-4 line only
         endif
      endif
      do k = 1,6
         res(k,1) = 0
         res(k,2) = 0
         do isig = 1,4
            do i = 1,2
               res(k,1) = res(k,1) + 
     1      l21*( dreal(mm21(k,isig,i))**2+dimag(mm21(k,isig,i))**2 )
               res(k,2) = res(k,2) + 
     1      l43*( dreal(mm43(k,isig,i))**2+dimag(mm43(k,isig,i))**2 )
            enddo
         enddo
         res(k,1) = res(k,1)*fac
         res(k,2) = res(k,2)*fac              
C   *p65(4)/xm2(6)
      enddo


      if(with_anomHiggs) then
c  set NC processes to zero if photon virtuality falls below cutoff
        if ( abs(p21(4)).lt.qsqAmin ) then
           do k = 1,4
              res(k,1) = 0      ! no graphs with gluon on lower line
           enddo
        endif
        if ( abs(p43(4)).lt.qsqAmin ) then
           do k = 1,4
              res(k,2) = 0      ! no graphs with gluon on upper line
           enddo
        endif
      endif

C  Assign the final results for the subprocesses
      if (lcol) then
         do i = 1,2
            uucc_c(i+1) = res(1,i) 
            uuss_c(i+1) = res(2,i)
            ddcc_c(i+1) = res(3,i)
            ddss_c(i+1) = res(4,i)
            udsc_c(i+1) = res(5,i)
            ducs_c(i+1) = res(6,i)
         enddo
         uucc_c(1) = res(1,1) +  res(1,2)
         uuss_c(1) = res(2,1) +  res(2,2)
         ddcc_c(1) = res(3,1) +  res(3,2)
         ddss_c(1) = res(4,1) +  res(4,2)
         udsc_c(1) = res(5,1) +  res(5,2)
         ducs_c(1) = res(6,1) +  res(6,2)
      else
         uucc = res(1,1) +  res(1,2)
         uuss = res(2,1) +  res(2,2)
         ddcc = res(3,1) +  res(3,2)
         ddss = res(4,1) +  res(4,2)
         udsc = res(5,1) +  res(5,2)
         ducs = res(6,1) +  res(6,2)
      endif
      end

