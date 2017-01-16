c
c	Initial version:  Dieter Zeppenfeld, 2003 February 12
c	Last modified by sophy@particle.uni-karlsruhe.de:   2012 February 24
c       From code by Barbara Jager
c
c  this routine computes |ME|^2 for qq->qq e+ve mu+vm g ("W+W+")
c	(real emission corr. to qq->qqVV, VV->leptons) 
c	 via call of "qqwwppqqj.f" 
c
c	flags can be set such that collin./soft subtraction is tested
c
c ------------------------------------------------------------------------
c
      subroutine m2s_wbfwppg(
     &                   bos,   !in:  Boson identifier, 
                                !     33=W+W+
     &                   nlo,   !in:  NLO=1: create subtraction term; LO = 0
     &                   lok,   !in:  need to calculate m2s(1:3)? T or F 
     &                   xi,	!in:  Feynman x parameters of incoming fermions
     &                   p,	!in:  momenta associated with external fermions
     &                   v,	!in:  momenta associated with WW decay fermions
     &                   rn,    !in:  random number array
     &                   xuz,   !in:  x, u, z for subtraction terms
     &                   m2s    !out: |M|^2*pdf1*pdf2 
     &                        )
          use globalvars, only: ldoblha
      implicit none
c
c declare input/output variables
c
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "tensor.inc"
#include "tensorz.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      integer N_gen_W
      real*8 fac_W, mjj2
      external mjj2

      integer bos
      real*8 xi(nx), p(0:3,max_p,max_kin), v(0:3,max_v), 
     1       rn(1), xuz(2,2:3), m2s(0:3)
      logical lok(3)

c  helicity selection
      integer h
      integer jsig, jsig1, jsig3, jsig5
      common /chelsum/ jsig,jsig1,jsig3,jsig5
c
      complex*16 zero
      parameter (zero=(0d0,0d0))
c
c declare external functions
c
      integer FL_ZHg
      external FL_ZHg
c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
      real*8 x1,x2
      real*8 q12(0:4,3), q34(0:4,3), !qww(0:3), 
     1 	     lnQomu(2:3), omxi(2:3), 
     1       Ax(2:3), Bx(2:3), Cx(2:3), Dxi(2:3), tgs2oqsq(2:3), 
     2       ln1mxi, lnrat, z, lnz
c
c declare local variables
c
c
      real*8 q_sf
      real*8 dotrr
      external dotrr

      integer init/0/, I, J, mu
      save init
c
c declare variables for summation process
c
      INTEGER ires, nmaxold, FSIGN(4+max_v), gsign
      save nmaxold

      integer physToDiag(5), nlo, nmin, nmax, nproc(8)
      save nmin, nmax, nproc
c
c store contributions from subprocess ip in res(ip,ID) where
c  ID = 1  : the real emission |M|^2 * pdf
c     = 2,3: sutraction terms for emision off 12 or 34 line
c     = 0  : subtracted result which drives integration, i.e
c res(*,0) = res(*,1)-res(*,2)-res(*,3)
      real*8 res(maxnumsubproc,0:3)

c
c declare parton distribution variables
c
      real*8 pdf(-6:6,2,3)
c
c define program switches
c
      logical ldebug
c      common/ dbug / ldebug	!store "debugging value"
      parameter (ldebug=.false.)
c      data ldebug /.true./		!output debug information
      logical lcoll,lsoft	! parameters for subtraction tests
      parameter (lcoll=.false.,lsoft=.false.)
c      
      common /hcount / h
      real*8 weight,rnumb,RandomNumber
c
c --------------------------------------------------------------
c
c if first iteration, output acceptance cuts and fix constant input 
c parameters for call of qqZWqqji.f
c
      if ( init .ne. 0) then
         do i = 1,nmax
            do j=0,3
               res(i,j) = 0
            enddo
         enddo
      else if (ldoblha) then
        h = 1
        init = 1
      else
         if (bos.eq.33) then      
            write(6,*) ' '
            if (procID .eq. WPWPjj .or. procID .eq. WPhadWPjj) then
               write(6,*) 'W+W+ jjj amplitude square information:'
            else
               write(6,*) 'W-W- jjj amplitude square information:'
            end if
            write(6,*) '-----------------------------------------------'
            write(6,*) ' '
            call vtowpm_reset
         endif
         call printnfl(.true.)
         print*," "
         print*," creal = ",creal," cvirtual = ",cvirt
         print*," virt factor for alphas = 0.12 is ",
     1            1+0.12/pi*4./3.*cvirt
         h = 1
         init = 1
         do i = 1,maxnumsubproc
            do j = 0,3
               res(i,j) = 0
            enddo
         enddo
      endif   ! end of init
      
c
c -------------------------------------------------------------------
      call Calc_Momentum_Transfer(p, v, q12,q34,3)


      if (bos.eq.33) then   ! W+W+

         if ((lha.or.hepmc).and..not.doNLO) then
            if (procID.eq.WPWPjj) then
               helicity(1)= 1
               helicity(2)=-1
               helicity(3)= 1
               helicity(4)=-1
            else
               helicity(1)=-1
               helicity(2)= 1
               helicity(3)=-1
               helicity(4)= 1
            endif
         endif

         CALL IXXXXX(v(0,1),ZERO , 1,-1,wel) !e+ 
         CALL OXXXXX(v(0,2),ZERO ,-1, 1,wve) !ve 
         CALL IXXXXX(v(0,3),ZERO , 1,-1,wmu) !mu+ 
         CALL OXXXXX(v(0,4),ZERO ,-1, 1,wvu) !vm 
         CALL JIOXXX(wel,wve,GWF ,WMASS,WWIDTH,we) !W+->e+ve
         CALL JIOXXX(wmu,wvu,GWF ,WMASS,WWIDTH,wu) !W+->mu+vm
         
         do mu = 0,3
            qe(mu)  = v(mu,1)+v(mu,2)
            qu(mu)  = v(mu,3)+v(mu,4)
            qww(mu) = qe(mu) + qu(mu)
         enddo
         qe(4)  = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
         qu(4)  = qu(0)**2-qu(1)**2-qu(2)**2-qu(3)**2
         qww(4) = qww(0)**2-qww(1)**2-qww(2)**2-qww(3)**2


         SELECT CASE(procid)
         CASE(WPhadWPjj, WMhadWMjj)

c           using global form factor for all tensors of one phase space point
c           this ensures proper cancellations for anomalous contributions
c           energy scale is invariant WW mass


           do j = 2,3

              call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qu(0))

c leptonic tensors
C for W+W+->e+ve mu+vm
              call wpptowpp_had(q12(0,j),q34(0,j),v,wppwpp(0,0,j))  	 
         
	 
c t-channel processes:	 
C for V1 V2->e+ve
              call WPVtoWP_had(1,1,q34(0,j),v,awen(0,0,1,j),
     -             zwen(0,0,1,j)) !Wu on upper
              call WPVtoWP_had(1,2,q12(0,j),v,awen(0,0,2,j),
     -             zwen(0,0,2,j)) !Wu on lower
c        
              call WPVtoWP_had(1,1,q12(0,j),v,waen(0,0,2,j),
     -             wzen(0,0,2,j)) !Wu on lower
              call WPVtoWP_had(1,2,q34(0,j),v,waen(0,0,1,j),
     -             wzen(0,0,1,j)) !Wu on upper

c for V1 V2->mu+vm
              call WPVtoWP_had(2,1,q34(0,j),v,awmu(0,0,1,j),
     -             zwmu(0,0,1,j)) !We on upper
              call WPVtoWP_had(2,2,q12(0,j),v,awmu(0,0,2,j),
     -             zwmu(0,0,2,j)) !We on lower
c        
              call WPVtoWP_had(2,1,q12(0,j),v,wamu(0,0,2,j),
     -             wzmu(0,0,2,j)) !We on lower
              call WPVtoWP_had(2,2,q34(0,j),v,wamu(0,0,1,j),
     -             wzmu(0,0,1,j)) !We on upper

           enddo                  !j

         CASE DEFAULT

         if (with_anom) then ! anomalous gauge boson couplings
c           using global form factor for all tensors of one phase space point
c           this ensures proper cancellations for anomalous contributions
c           energy scale is invariant WW mass


           do j = 2,3

              call anomal_formfactor(q12(0,j),q34(0,j),qe(0),qu(0))

c leptonic tensors
C for W+W+->e+ve mu+vm
              call wpptowpp_anomal(q12(0,j),q34(0,j),v,wppwpp(0,0,j))  	 
         
	 
c t-channel processes:	 
C for V1 V2->e+ve
              call WPVtoWP_anomal(1,1,q34(0,j),v,awen(0,0,1,j),
     -             zwen(0,0,1,j)) !Wu on upper
              call WPVtoWP_anomal(1,2,q12(0,j),v,awen(0,0,2,j),
     -             zwen(0,0,2,j)) !Wu on lower
c        
              call WPVtoWP_anomal(1,1,q12(0,j),v,waen(0,0,2,j),
     -             wzen(0,0,2,j)) !Wu on lower
              call WPVtoWP_anomal(1,2,q34(0,j),v,waen(0,0,1,j),
     -             wzen(0,0,1,j)) !Wu on upper

c for V1 V2->mu+vm
              call WPVtoWP_anomal(2,1,q34(0,j),v,awmu(0,0,1,j),
     -             zwmu(0,0,1,j)) !We on upper
              call WPVtoWP_anomal(2,2,q12(0,j),v,awmu(0,0,2,j),
     -             zwmu(0,0,2,j)) !We on lower
c        
              call WPVtoWP_anomal(2,1,q12(0,j),v,wamu(0,0,2,j),
     -             wzmu(0,0,2,j)) !We on lower
              call WPVtoWP_anomal(2,2,q34(0,j),v,wamu(0,0,1,j),
     -             wzmu(0,0,1,j)) !We on upper

           enddo                  !j


         else !SM      

c leptonic tensors
C for W+W+->e+ve mu+vm
           do j = 2,3

              call wpptowpp(q12(0,j),q34(0,j),v,wppwpp(0,0,j))  	 
         
	 
c t-channel processes:	 
C for V1 V2->e+ve
              call WPVtoWP(1,1,q34(0,j),v,awen(0,0,1,j),
     -             zwen(0,0,1,j)) !Wu on upper
              call WPVtoWP(1,2,q12(0,j),v,awen(0,0,2,j),
     -             zwen(0,0,2,j)) !Wu on lower
c        
              call WPVtoWP(1,1,q12(0,j),v,waen(0,0,2,j),
     -             wzen(0,0,2,j)) !Wu on lower
              call WPVtoWP(1,2,q34(0,j),v,waen(0,0,1,j),
     -             wzen(0,0,1,j)) !Wu on upper

c for V1 V2->mu+vm
              call WPVtoWP(2,1,q34(0,j),v,awmu(0,0,1,j),
     -             zwmu(0,0,1,j)) !We on upper
              call WPVtoWP(2,2,q12(0,j),v,awmu(0,0,2,j),
     -             zwmu(0,0,2,j)) !We on lower
c        
              call WPVtoWP(2,1,q12(0,j),v,wamu(0,0,2,j),
     -             wzmu(0,0,2,j)) !We on lower
              call WPVtoWP(2,2,q34(0,j),v,wamu(0,0,1,j),
     -             wzmu(0,0,1,j)) !We on upper

           enddo                  !j

         endif                   ! anom or SM

         END SELECT

      endif                     ! end of bos = ...

c
c  ---------------------------------------------------------------------
c
cc  debugging:
cc
cc collinear subtraction:
c	if (lcoll) then 
c      	    ldebug = lok(1).and.lok(2).and.lok(3).and.
c     #      	    dotrr(p(0,5,1),p(0,1,1)).lt.1
c	    if (ldebug) print*,"coll.dot.:",lok(1),lok(2),lok(3),
c     #      		       dotrr(p(0,5,1),p(0,1,1))
cc      
cc soft subtraction:	
c       elseif (lsoft) then		       
c	   ldebug = lok(1).and.lok(2).and.lok(3).and.p(0,5,1).lt.1	
c	   if (ldebug) print*,"soft.mom.:",lok(1),lok(2),lok(3),p(0,5,1)   
c       endif
c
c  ---------------------------------------------------------------------
c
c scales and als:
      if (.false.) then
         write(6,*) " "
         write(6,*) "xi(1), xi(2) =", xi
         write(6,*) "mu_f L=1:",sqrt(mufsq(1,1)),sqrt(mufsq(2,1))
         write(6,*) "mu_f L=2:",sqrt(mufsq(1,2)),sqrt(mufsq(2,2))
         write(6,*) "mu_f L=3:",sqrt(mufsq(1,3)),sqrt(mufsq(2,3))
         write(6,*) " "
         write(6,*) "xi(1), xi(2) =", xi
         write(6,*) "alphas L=1:",als(1,1), als(2,1)
         write(6,*) "alphas L=2:",als(1,2), als(2,2)
         write(6,*) "alphas L=3:",als(1,3), als(2,3)
      end if
      if (lwarn) then           ! check connection between scales
         if(abs(mufsq(1,2)/mufsq(1,1)-1).gt.1d-8) then
            print*," check muf 1,2 vs 1,1 ",mufsq(1,2),mufsq(1,1)
         endif
         if(abs(mufsq(2,3)/mufsq(2,1)-1).gt.1d-8) then
            print*," check muf 2,3 vs 2,1 ",mufsq(2,3),mufsq(2,1)
         endif
      endif

      if ( ldebug ) then

 10      format( " p(", i1, ") = ", 4(f10.3, 2x) )
 20      format( " v(", i1, ") = ", 4(f10.3, 2x) )
c
         write(6,*) " "
         print*," vector boson decay products momenta "
         do i = 1, n_v
            write(6,20) i, v(0,i), v(1,i), v(2,i), v(3,i)
         end do
         print*," parton momenta "
         do i = 1, n_p
            write(6,10) i, p(0,i,1), p(1,i,1), p(2,i,1), p(3,i,1)
         end do
c
         write(6,*) " "
         write(6,*) "xi(1), xi(2) =", xi
         write(6,*) "mu_f =",sqrt(mufsq(1,1)),sqrt(mufsq(2,1))
         read(*,*)

      end if
c
c  end debugging
c
c  ---------------------------------------------------------------------

       if (ldoblha) then

         if (procid .eq. WmWmjj) then
           physToDiag(1)=blha_physdiagmap(1,blha_cursubproc,blha_curproc) 
           physToDiag(2)=blha_physdiagmap(2,blha_cursubproc,blha_curproc)
           physToDiag(3)=blha_physdiagmap(3,blha_cursubproc,blha_curproc)
           physToDiag(4)=blha_physdiagmap(4,blha_cursubproc,blha_curproc)
           physToDiag(5)=blha_physdiagmap(5,blha_cursubproc,blha_curproc)
c flip both fermion lines
           do i=1,5
             if (physToDiag(i).eq.1) then 
               physToDiag(i)=2 
             elseif (physToDiag(i).eq.2) then 
               physToDiag(i)=1 
             elseif (physToDiag(i).eq.3) then 
               physToDiag(i)=4 
             elseif (physToDiag(i).eq.4) then 
               physToDiag(i)=3 
             endif
           enddo
c first argument of blha_fsign is diagtophys(i)!
           fsign(1) = -blha_fsign(2,blha_cursubproc,blha_curproc)
           fsign(2) = -blha_fsign(1,blha_cursubproc,blha_curproc)
           fsign(3) = -blha_fsign(4,blha_cursubproc,blha_curproc)
           fsign(4) = -blha_fsign(3,blha_cursubproc,blha_curproc)
           gsign    =  blha_fsign(5,blha_cursubproc,blha_curproc) ! initial-/final-state is conserved
         else
           physToDiag(1)=blha_physdiagmap(1,blha_cursubproc,blha_curproc)
           physToDiag(2)=blha_physdiagmap(2,blha_cursubproc,blha_curproc)
           physToDiag(3)=blha_physdiagmap(3,blha_cursubproc,blha_curproc)
           physToDiag(4)=blha_physdiagmap(4,blha_cursubproc,blha_curproc)
           physToDiag(5)=blha_physdiagmap(5,blha_cursubproc,blha_curproc)
           fsign(1) = blha_fsign(1,blha_cursubproc,blha_curproc)
           fsign(2) = blha_fsign(2,blha_cursubproc,blha_curproc)
           fsign(3) = blha_fsign(3,blha_cursubproc,blha_curproc)
           fsign(4) = blha_fsign(4,blha_cursubproc,blha_curproc)
           gsign    = blha_fsign(5,blha_cursubproc,blha_curproc)
         endif

         call wbf_wpp3j(xi,p,v,physToDiag,fsign,gsign,bos,  nlo,lok,xuz,
     1            pdf,res,nmin,nmax)

         return
       endif

c
c call PDF subroutine in order to determine parton
c distributions in the incoming (anti)protons.
c for the NLO contributions x1 = x*y in my notes with x=z=xuz(1,j), y=xi(j-1)
c
      x1 = xi(1)*xuz(1,2)
      q_sf = sqrt(mufsq(1,1))
      call pdfproton( xi(1), q_sf, pdf(-6,1,1) )      ! f_a(y)=f_a(x1/z)

      if (nlo.gt.0) then
      q_sf = sqrt(mufsq(1,2))                         ! f_a(x1) for upper line
      call pdfproton( x1, q_sf, pdf(-6,1,2) )         !   NLO correction
      endif

      if (mufsq(1,3).ne.mufsq(1,1)) then              ! f_a(x1) for lower line
         q_sf = sqrt(mufsq(1,3))                      !   NLO correction
         call pdfproton( xi(1), q_sf, pdf(-6,1,3) )
      else
         do i = -6,6
            pdf(i,1,3) = pdf(i,1,1)
         enddo
      endif

      x2 = xi(2)*xuz(1,3)
      q_sf = sqrt(mufsq(2,1))
      call pdfproton( xi(2), q_sf, pdf(-6,2,1) )      ! f_b(y)=f_a(x2/z)

      if (nlo.gt.0) then
      q_sf = sqrt(mufsq(2,3))                         ! f_b(x2) for lower line
      call pdfproton( x2, q_sf, pdf(-6,2,3) )         !   NLO correction
      endif

      if (mufsq(2,2).ne.mufsq(2,1)) then              ! f_b(x2) for upper line
         q_sf = sqrt(mufsq(2,2))                      !   NLO correction
         call pdfproton( xi(2), q_sf, pdf(-6,2,2) )
      else
         do i = -6,6
            pdf(i,2,2) = pdf(i,2,1)
         enddo
      endif

c
c and fill the coefficient functions for the finite subtraction terms
C Note that the color factors TR and C2 are NOT!!! included here
      lnQomu(2) = log(q12(4,2)/mufsq(1,2))
      lnQomu(3) = log(q34(4,3)/mufsq(2,3))
      omxi(2) = 1d0-x1
      omxi(3) = 1d0-x2
      tgs2oqsq(2) = 8d0*pi*als(1,2)/q12(4,2)         !2g_s^2/Q^2
      tgs2oqsq(3) = 8d0*pi*als(2,3)/q34(4,3)         !2g_s^2/Q^2

      do j = 2,3
         z = xuz(1,j)
         ln1mxi = log(omxi(j))
c         tgs2oqsq(j) = tgs2oqsq(j) * 6*xuz(2,j)*(1d0-xuz(2,j))
C cut off Int_x1^1 dz log(1-z) at z<1-5E-6. This assures that the relative
C error of Int_x1^1 dz log(1-z) is less than 1E-4/(1-x1). see notes p 28.4
         if (z.lt.0.999995) then
            lnrat = lnQomu(j) + log((1d0-z)/z)
            lnz = log(z)
            Ax(j) = (z**2+(1-z)**2)*lnrat + 2*z*(1-z)
            Bx(j) = ( 2d0*(lnrat+lnz)-1.5d0 )/(1d0-z)
            Cx(j) = 1-z - 2d0*lnz/(1d0-z) - (1+z)*lnrat
         else
            Ax(j) = 0d0
            Bx(j) = 0d0
            Cx(j) = 0d0
         endif
         Dxi(j) = 1.5d0*lnQomu(j) + 2d0*ln1mxi*lnQomu(j)
     1            + ln1mxi**2 - 1.5d0*ln1mxi 
     2            + creal
         do i = 1,5
            pdf(-i,j-1,j)=(pdf(-i,j-1,j)*(Dxi(j)/omxi(j)-Bx(j)) +
     1                     pdf(-i,j-1,1)*(Bx(j)+Cx(j)))*tgs2oqsq(j)
            pdf(i,j-1,j)= (pdf(i,j-1,j)*(Dxi(j)/omxi(j)-Bx(j)) +
     1                     pdf(i,j-1,1)*(Bx(j)+Cx(j)))*tgs2oqsq(j)
         enddo
         pdf(0,j-1,j)=pdf(0,j-1,1)*Ax(j)*tgs2oqsq(j)
      enddo
c select helicity
      jsig = 0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Now sum over the subprocesses contributing to H production

      nmax = FL_ZHg(fsign,-1)            !reset counter for subprocesses to 0

C*******************  q1 q3 ---> q2 q4 g V V   **********************
         
c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
      physToDiag(5)=5   ! gluon
C NOTE: for call of wbf_wpp3j it is important that p(*,1,*) and p(*,3,*)
c correspond to 1-2 fermion line ALWAYS, i.e physToDiag(1/2)={1,3} and 
c similarly physToDiag(3/4)={2,4} for the 3-4 fermion line
      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = 1
      fsign(4) = 1
      gsign    = 1

      call wbf_wpp3j(xi,p,v,physToDiag,fsign,gsign,bos,  nlo,lok,xuz,
     1  	    pdf,res,nmin,nmax)

c      if (ldebug) call resprint(nmin,nmax,res)
      if (init.eq.1) nproc(1) = nmax

C*******************  q1 qb4 ---> q2 qb3 g V V   **********************

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3
      physToDiag(5)=5   ! gluon

      fsign(3) = -1
      fsign(4) = -1

      call wbf_wpp3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	   pdf,res,nmin,nmax)

c      if (ldebug) call resprint(nmin,nmax,res)
      if (init.eq.1) nproc(2) = nmax

C*******************  qbar2 q3 ---> qbar1 q4 g V V   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
      physToDiag(5)=5   ! gluon

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = 1
      fsign(4) = 1
      
      call wbf_wpp3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	   pdf,res,nmin,nmax)

c      if (ldebug) call resprint(nmin,nmax,res)
      if (init.eq.1) nproc(3) = nmax

C*******************  qbar2 qb4 ---> qbar1 qb3 g V V  ******************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3
      physToDiag(5)=5   ! gluon

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1 
      fsign(4) = -1
         
      call wbf_wpp3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	   pdf,res,nmin,nmax)

c      if (ldebug) call resprint(nmin,nmax,res)
      if (init.eq.1) nproc(4) = nmax

C*******************  g q3 ---> qb1 q2 q4 V V  **********************

      physToDiag(1)=5
      physToDiag(2)=3
      physToDiag(3)=1
      physToDiag(4)=4
      physToDiag(5)=2

      fsign(1) = -1
      fsign(2) =  1
      fsign(3) =  1
      fsign(4) =  1
      gsign    = -1

      call wbf_wpp3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	   pdf,res,nmin,nmax)

c      if (ldebug) call resprint(nmin,nmax,res)
      if (init.eq.1) nproc(5) = nmax

C*******************  q1 g ---> q2 qb3 q4 V V   **********************

      physToDiag(1)=1
      physToDiag(2)=5
      physToDiag(3)=2
      physToDiag(4)=3
      physToDiag(5)=4

      fsign(1) =  1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1

      call wbf_wpp3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	   pdf,res,nmin,nmax)

c      if (ldebug) call resprint(nmin,nmax,res)
      if (init.eq.1) nproc(6) = nmax

C*******************  g qbar4 ---> qbar3 qb1 q2 V V   **********************

      physToDiag(1)=5
      physToDiag(2)=4
      physToDiag(3)=1
      physToDiag(4)=3
      physToDiag(5)=2

      fsign(1) = -1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) = -1

      call wbf_wpp3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	   pdf,res,nmin,nmax)

c      if (ldebug) call resprint(nmin,nmax,res)
      if (init.eq.1) nproc(7) = nmax

C*******************  qbar2 g ---> qbar1 qb3 q4 V V  **********************

      physToDiag(1)=2
      physToDiag(2)=5
      physToDiag(3)=1
      physToDiag(4)=3
      physToDiag(5)=4
      
      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) =  1

      call wbf_wpp3j(xi,p,v,physToDiag,fsign,gsign,bos, nlo,lok,xuz,
     1  	   pdf,res,nmin,nmax)

c      if (ldebug) call resprint(nmin,nmax,res)
      if (init.eq.1) nproc(8) = nmax

c*****************  end of process evaluation  *********************

      if (init.eq.1) then
         init = init+1
         if (lwarn) print 199," proc #s for Z/Hjjj are ",nproc
 199     format(a,8i5)
      endif


      SELECT CASE(procid)
      CASE(WPhadWPjj, WMhadWMjj)
       call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,1),v(0,2))), N_gen_W, fac_W)
       do j=0,3
        DO IRES = 1,NMAX
          res(IRES,j) = res(IRES,j) * fac_W            ! factor fac_W for hadronic decay
        ENDDO
       enddo
      END SELECT


      do j=0,3
         m2s(j) = 0      
         DO IRES = 1,NMAX
            m2s(j) = m2s(j) + RES(IRES,j)
         ENDDO

         if(j.eq.0)then
c...Les Houches interface - the most propable subprocess 3jets at LO  
            if ((lha.or.hepmc).and..not.doNLO) then
               i=0
               weight=0.d0
               rnumb=RandomNumber()
               do while((i.le.nmax).and.(weight.le.rnumb*m2s(0)))
                  i=i+1
                  weight=weight+res(i,0)
                  iprocess=i
               enddo
               SELECT CASE(procid)
               CASE(WPhadWPjj)
                  if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
                     rnumb=RandomNumber()
                     finalquarks_psp(1) =  2 + 2* INT(rnumb*2)
                     finalquarks_psp(2) = -1 - 2* INT(rnumb*2)
                  endif
               CASE(WMhadWMjj)
                  if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
                     rnumb=RandomNumber()
                     finalquarks_psp(1) =  1 + 2* INT(rnumb*2)
                     finalquarks_psp(2) = -2 - 2* INT(rnumb*2)
                  endif
               END SELECT
            endif
         endif

c         DO IRES = nproc(4)+1,nproc(8)
c            m2s(j) = m2s(j) + RES(IRES,j)
c         ENDDO

      enddo  ! j loop


c
c ---------------------------------------------------------------
c debugging:
      if (ldebug) then
         do j=1,3
            if (lok(j) .and. m2s(j).eq.0 ) then
               print*,j," lok(j) = ",lok(j)," m2s = ",m2s(j)
            endif
         enddo
         if (m2s(0).ne.0) then
c            if (abs((m2s(1)+m2s(2)+m2s(3))/m2s(0)-1).gt.1d-5) then
               print*," m2s(0) = ",m2s(0),(m2s(1)+m2s(2)+m2s(3))/m2s(0)
               print*," ratio subtr/real = ",(m2s(2)+m2s(3))/m2s(1)
               print*," m2s(0,1,2,3) = ",m2s
	       print*
c            endif
         else
            print*," m2s = ",m2s
         endif
         if (nmax.ne.nmaxold) print*," nmax: ",nmaxold,nmax
         nmaxold = nmax
      endif
c      if (ldebug) print*," lokt = ",lok
c      if (ldebug) read(*,113) text
c 113  format(a)

      RETURN
      END


c******************************************************************
c
c   begin subroutine wbf_wpp3j
c
c*****************************************************************
      subroutine wbf_wpp3j(xi,p,v,physToDiag,fsign,gsign,bos,
     1                    nlo,lok,xuz,pdf,res,nmin,nmax)
          use globalvars, only: ldoblha

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"

      real*8 p(0:3,max_p,max_kin), v(0:3,max_v), xi(nx), xuz(2,2:3)
      real*8 pdf(-6:6,2,3), res(maxnumsubproc,0:3)
      integer physToDiag(5), fsign(4+max_v), 
     1        gsign, h, bos, nlo, nmin, nmax
      logical lok(3)
c      
      common /hcount / h
c
c wbf_wpp3j calls the amplitude square routines 
c             qqwppqq    for qq-->qqW+W+   for bos=33
c for the subtraction terms for the NLO cross section calculation
c
c  INPUT:  p(0:3,5,3)      external physical parton momenta
c          v(0:3,nv)       Z decay momenta, for Higgs production only
c                          the sum q(mu) = v(mu,1)+...+v(mu,nv) is needed
c          physToDiag(5)   physToDiag(ext.mom. label) = Feynman diagram label
c          fsign,gsign     sign factors for amplitude calls; see qqZqq(j)
c	   h		   specify lepton helicity combination
c          nlo             nlo = 0: calculate LO only (i.e. no subtraction)
c                          nlo = 1: calculate full NLO subtraction
c          lok(3)          lok(ID)=T means momenta set ID passes acceptance
c                          cuts, i.e res(k,ID) needs to be calculated
c  OUTPUT:
c          uucc(ID)   output in format as in qqZqq(j), but corresponding
c          etc.            to res(*,ID) in m2s_qqZqq
c                          ID = 1  : the real emission |M|^2 * pdf
c                             = 2,3: sutraction terms for emision off 12 or 
c     or   res(k,ID)                 34 line
c                             = 0  : subtracted result which drives 
c                                    integration, i.e 
c                                    res(*,0) = res(*,1)+res(*,2)+res(*,3)
c
c  In and Output
c    nmin, nmax            range of process numbers for this call
c
      real*8 C2, TR, N                              ! color factors
      parameter(N=3d0,TR=0.5d0,C2=TR*(N**2-1d0)/N)

c Note: Factor 9 from color sum included in qqbqq. I am using color summed
c amplitudes here, while CS use color averages. Thus I later divide by
c 8*3 for an initial qg state and by 3*3 for an initial qq state. So
c
c  2*gs2*C2*9/(8*3) = 2 gs^2 * TR *8/3*9/(8*3) = 2 gs^2 * TR
c
c consistent with Catani 
c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
c
c  helicity selection
      integer jsig, jsig1, jsig3, jsig5
      common /chelsum/ jsig,jsig1,jsig3,jsig5
c
      real*8 pbar(0:3,4+max_v), qbar(0:4),q2,v2,q3,v3,sub(2:3), 
     1       e_in(2), dotrr, gs2(2:3), polcolq, polcolg

      integer iflav(5), diagToPhys(5), FL_ZHg
      external dotrr, FL_ZHg

      logical loldp, ldebug
      data e_in/0d0,0d0/
      save e_in,loldp, gs2,polcolq,polcolg
      integer i,if1,if2,j,k,mu, fc(6), fcb(6)
      
      parameter (ldebug=.false.)

      double precision matr3sq(2:3)
      double precision matr2sq(2:3)
      real*8 ucdsb(-1:1,-1:1,2:3),ucds(0:3)
      real*8 dsucb(-1:1,-1:1,2:3),dsuc(0:3)
      integer if1min,if2min

      save ucdsb, dsucb

      real*8 ucds_temp(0:3)  ! for debugging


      do j = 0,3
         ucds(j) = 0d0
         dsuc(j) = 0d0
      enddo
      do j = 2,3
         matr3sq(j) = 0d0
         matr2sq(j) = 0d0
      enddo


      nmin = nmax+1
      do i = 1,5
         diagToPhys(physToDiag(i)) = i
      enddo

      loldp = e_in(1).eq.p(0,3,1) .and. e_in(2).eq.p(0,4,1)
c  reset the LO amplitude to 0 to avoid wrong subtraction for initial gluon
      if (.not.loldp.or.ldoblha) then            ! this is a new phase space point
         e_in(1) = p(0,3,1)
         e_in(2) = p(0,4,1)
         do i = -1,1
            do j = -1,1
               do k = 2,3
                  ucdsb(i,j,k) = 0d0
                  dsucb(i,j,k) = 0d0
               enddo
            enddo
         enddo
c  determine strong coupling gs for the two quark lines and factor for 
c  polarization and spin average
         gs2(2) = 4d0*pi*als(1,1)
         gs2(3) = 4d0*pi*als(2,1)
         polcolq = 1d0/(4d0*N**2*xi(1)*xi(2))
         polcolg = 1d0/(4d0*N*(N**2-1)*xi(1)*xi(2))
      endif
         
c get the real emission amplitude squared, store it in uucc(-1,1) etc.
      do mu = 0,3
         do i = 1,5
            pbar(mu,physToDiag(i))=p(mu,i,1)
         enddo
         do i = 6,4+max_v
            pbar(mu,i)=0d0
         enddo
         qbar(mu) = pbar(mu,5)
      enddo
      qbar(4)=0d0


      if (bos .eq. 33) then  !  W+W+jj
c
         do mu = 0,3            ! kinematics for 4-lepton final state
            pbar(mu,5) = v(mu,1) ! l+
            pbar(mu,6) = v(mu,2) ! nu
            pbar(mu,7) = v(mu,3) ! l'+
            pbar(mu,8) = v(mu,4) ! nu'
         enddo
         fsign(5) = -1
         fsign(6) = 1
         fsign(7) = -1
         fsign(8) = 1         
         
	 if (lok(1).or.nlo.eq.0) then	 
c            call qqwpmqqj_comp(pbar,fsign,qbar,gsign,bos,ucds(1)) 
            call qqwppqqj(pbar,fsign,qbar,gsign,bos,ucds(1)) 
	 endif
          		 
      else
         print*," Invalid entry BOS = ",bos
         print*," in call of wbf_wpp3j. Must be 33"
         stop
      endif   ! end of if bos = ..

      if (ldoblha) then
        return
      endif


      if (bos.eq.33) then  ! W+W+jj
         do j = 2,3
            matr3sq(j) = ucds(j) 	    
         enddo
      endif                     !bos

c for the NLO case get the subtraction terms; 
c first the case with a final state gluon
      if (gsign.eq.1) then
         if (ldebug) then
            print*," final state gluon section in wbf_wpp3j "
            print*," jsig = ",jsig," jsig1,3 = ",jsig1,jsig3
            print 101," fsign = ",fsign
 101        format(a,6i5,a,i5)
         endif
         if (nlo.eq.1) then
            do j = 2,3               ! j=2: emission off 1,2 line
               do mu = 0,3           ! j=3: emission off 3,4 line
                  do i = 1,4
                     pbar(mu,physToDiag(i))=p(mu,i,j)
                  enddo
               enddo

               if (bos.eq.33.and.lok(j)) then  ! W+W+jj
                     call qqwppqq(pbar,fsign,0,j,bos,ucds(j))
               else
                  ucds(j) = 0d0
                  dsuc(j) = 0d0 
               endif
	       	       
	       if (bos.eq.33) then   ! W+W+jj
                  matr2sq(j) = ucds(j)
c save matrix elements for later use with initial gluons
                  ucdsb(fsign(1),fsign(3),j) = ucds(j)
               end if
c
               q2 = 2*dotrr(qbar,p(0,j-1,1))*xuz(1,j) !p(mu,j-1,1) is inc.quark
               v2 = 2d0*gs2(j)*                      !i.e  j-1 = 1,2
     $            ( 2d0/(1-xuz(1,j)+xuz(2,j))-(1+xuz(1,j)))
               q3 = 2*dotrr(qbar,p(0,j+1,1))*xuz(1,j) !p(mu,j+1,1) is out.quark
               v3 = 2d0*gs2(j)*                      !i.e. j+1 = 3,4
     $            ( 2d0/(1-xuz(1,j)+xuz(2,j))-(2-xuz(2,j)))
               sub(j) = v2/q2+v3/q3
            enddo  ! loop j
         endif   ! nlo=1
         iflav(5) = 0           ! final state gluon id

	
         IF (bos.eq.33) THEN    ! W+W+
		 
        ! ifmin = 2 for up and charm quark, 1 for down and strange anti-quark
            if1min = (fsign(physToDiag(1))+3)/2
            if2min = (fsign(physToDiag(2))+3)/2

            do if1 = if1min,int(nfl/2)*2,2 
               do if2 = if2min,int(nfl/2)*2,2
                  iflav(1)=if1*fsign(physToDiag(1)) !in
                  iflav(2)=if2*fsign(physToDiag(2)) !in
	    
                  iflav(3)=iflav(1)-1 !out
                  iflav(4)=iflav(2)-1 !out  
		  		  
                  do j = 2,3
                     
                     k=FL_ZHg(iflav,j)
                     if (lok(1)) then
                        res(k,1)=pdf(sign1*iflav(1),1,2*j-3) !1 for j=2;
                                                             !3 for j=3
     &                       *pdf(sign2*iflav(2),2,4-j  ) !2 for j=2;1 for j=3
     &                       *matr3sq(j)*polcolq
                     else
                        res(k,1) = 0
                     endif
                     if (nlo.eq.1 .and. lok(j)) then
                        res(k,j)=(pdf(sign1*iflav(1),1,j)
     &                           *pdf(sign2*iflav(2),2,j)
     &                           -pdf(sign1*iflav(1),1,2*j-3)
     &                           *pdf(sign2*iflav(2),2,4-j  )*sub(j))*C2
     &                           *matr2sq(j)*polcolq
                        res(k,5-j)=0
                        res(k,0) = res(k,1)+res(k,j)
                     else
                        res(k,0) = res(k,1)
                        res(k,2) = 0
                        res(k,3) = 0
                     endif
                  enddo         ! j = 2,3

c s-channel exchange diagrams:     	  
c	   if (fsign(physToDiag(1))*fsign(physToDiag(2)).eq.-1) then
c	   	! one quark plus one anti-quark
c	     
c	     if ((abs(if1*if2).eq.2).or.(abs(if1*if2).eq.12)) then
c	     	! same generation only
c	   
c	   	fs3 = fsign(physToDiag(1))
c	   	fs4 = fsign(physToDiag(2))
c	   
c	      	iflav(3)=fs3*5-iflav(1)
c	      	iflav(4)=fs4*5-iflav(2)
c
c                  do j = 2,3
c	          k=FL_ZHg(iflav,1)
c	   
c                     if (lok(1)) then
c                        res(k,1)=pdf(iflav(1),1,2*j-3)  !1 for j=2;3 for j=3
c     &                          *pdf(iflav(2),2,4-j  )  !2 for j=2;1 for j=3
c     &                          *matr3sq(j)*polcolq
c                     else
c                        res(k,1) = 0
c                     endif
c                     if (nlo.eq.1 .and. lok(j)) then
c                        res(k,j)=(pdf(iflav(1),1,j)
c     &                           *pdf(iflav(2),2,j)
c     &                           -pdf(iflav(1),1,2*j-3)
c     &                           *pdf(iflav(2),2,4-j  )*sub(j))*C2
c     &                           *matr2sq(j)*polcolq
c                        res(k,5-j)=0
c                        res(k,0) = res(k,1)+res(k,j)
c                     else
c                        res(k,0) = res(k,1)
c                        res(k,2) = 0
c                        res(k,3) = 0
c                     endif
c                  enddo    ! j = 2,3
c     
c     	    endif  !abs ...     
c	  endif  !fsign		  
		  
               enddo            ! if2
            enddo               ! if1

         end if  ! end W+W+ or oTw

      elseif (gsign.eq.-1) then           !initial gluon section
         j = 0
         if (nlo.eq.0) then
            if (fsign(1).eq.-fsign(2)) then 
               j=2              ! j=2: emission off 1,2 line
            elseif(fsign(3).eq.-fsign(4)) then 
               j=3              ! j=3: emission off 3,4 line
            endif
         elseif (nlo.eq.1) then
            do i=1,6
               fc(i)=fsign(i)
               fcb(i)=fsign(i)
            enddo
            if (fsign(1).eq.-fsign(2)) then 
               j=2              ! j=2: emission off 1,2 line
               fc(2)=-fsign(2)
               q2 = 2*dotrr(qbar,pbar(0,2))*xuz(1,j)
               fcb(1)=-fsign(1)
               q3 = 2*dotrr(qbar,pbar(0,1))*xuz(1,j)
            elseif(fsign(3).eq.-fsign(4)) then 
               j=3              ! j=3: emission off 3,4 line
               fc(4)=-fsign(4)
               q2 = 2*dotrr(qbar,pbar(0,4))*xuz(1,j)
               fcb(3)=-fsign(3)
               q3 = 2*dotrr(qbar,pbar(0,3))*xuz(1,j)
            endif
            v2 = 2d0*gs2(j) * ( (1-xuz(1,j))**2 + xuz(1,j)**2 )
            v3 = v2


            if (bos.eq.33) then  ! W+W+
               matr2sq(2) = ucdsb(fc(1),fc(3),j)
               matr2sq(3) = ucdsb(fcb(1),fcb(3),j)
            endif
         end if  ! nlo

         iflav(j-1) = 0 !gluon


         IF (bos.eq.33) THEN
 
       ! ifmin = 2 for up quark, 1 for down anti-quark
            if1min = (fsign(physToDiag(4-j))+3)/2 ! for incoming (anti-)quark
            if2min = (-fsign(physToDiag(j+1))+3)/2 ! for first outcoming 
                                                   ! (anti-)quark
            
            do if1 = if1min,int(nfl/2)*2,2
               do if2=if2min,int(nfl/2)*2,2
                  
                  iflav(4-j) = if1*fsign(physToDiag(4-j)) !2 or 
                  iflav(6-j) = if1*fsign(physToDiag(6-j))-1 !4 or 
                  iflav(j+1) = if2*fsign(physToDiag(j+1)) !3 or
                  iflav(5)   = if2*fsign(physToDiag(5))-1
		
                  k=FL_ZHg(iflav,j)
                  if (lok(1)) then
                     res(k,1)=pdf(sign1*iflav(1),1,2*j-3) !1 for j=2;3 for j=3
     &                    *pdf(sign2*iflav(2),2,4-j) !2 for j=2;1 for j=3
     &                    *matr3sq(j)*polcolg
                  else
                     res(k,1) = 0
                  endif
                  if ( nlo.eq.1 .and. lok(j) ) then
                     res(k,j) = ( pdf(sign1*iflav(1),1,j)
     &                    *pdf(sign2*iflav(2),2,j)
     &                    * (matr2sq(2) + matr2sq(3))
     &                    - pdf(sign1*iflav(1),1,2*j-3)
     &                    * pdf(sign2*iflav(2),2,4-j) *
     &                    (matr2sq(2)*v2/q2 + matr2sq(3)*v3/q3) )
     &                    * C2*polcolg
                     res(k,0) = res(k,1)+res(k,j)
                     res(k,5-j)=0
                  else
                     res(k,0) = res(k,1)
                     res(k,2) = 0
                     res(k,3) = 0
                  endif
                  
               enddo            !if2
            enddo               !if1
            
         end if  !W+W+ or not

      end if

      nmax = FL_ZHg(iflav,-2)

      end
c******************************************************************
c
c   end subroutine wbf_wpp3j
c
c*****************************************************************

