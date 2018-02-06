c     Modified by Terrance Figy
c     2 Oct 2008
C
      subroutine init_qqhqqg_virt 
          use globalvars, only: ldoblha
c >0  includes matrix elements and coupl
      implicit none             !        for H -> xx (based on bbar)
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
#include "VBFNLO/utilities/mssm.inc"
      double precision  fpi
      parameter ( fpi=4d0*pi)
        
c    This subroutine computes the born graphs and finite virtual graphs.
c
C     nlo =  0  |M_born|^2
c     nlo =  1  |M_born|^2 + 2Re[M_born conj(M_virt)]
c     nlo =  2  |M_born|^2 + 2Re[M_born conj(M_virt-M_virt,boxline)]
c     nlo = -2  2Re[M_born conj(M_virt,boxline)]
c
C  QQHQQ calculates the matrix elements**2 for light Higgs production by
C  electroweak boson fusion in quark quark scattering
C
C        q1 q3    ---->   q2 q4 g H
C
C  QQHQQI must be called first to initialize some couplings
C
C  The main task of QQHQQI is to set up the products of coupling constants
C  needed in Feynman graphs a, ..., g which are stored in
      real*8 fcpl(4,6)
C  and which are calculated from output of KOPPLN
      real*8 clr,xm2,xmg,b,v,a
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)
#include "VBFNLO/utilities/koppln.inc"  
c CPS scheme
      double precision qgammaq
      common /VBFNLO_HIGGS_CPS/ qgammaq
C
C  Here fcpl(sig,i) contains the coupling constant factors. 
C     sig = 1,..,4    identifies the 4 different helicity combinations of the
C                     quarks
C     i   = 1,..,6    numbers the possible subprocesses (4 NC and 2 CC)
c
C  
C  The various processes are identified by the following
      integer v2ini(6), v3ini(6)
      integer fl1(6), fl3(6)
      integer sig1(4,2),sig3(4,2)

      data fl1 /3,3,4,4,3,4/
      data fl3 /3,4,3,4,4,3/

      save fl1,fl3
      common /cqqhqqj/ fcpl, sig1,sig3, v2ini,v3ini
      data v2ini /2,2,2,2,3,4/, v3ini /2,2,2,2,4,3/
      data sig1 /-1,-1, 1, 1,-1,3*0/
      data sig3 /-1, 1,-1, 1,-1,3*0/
c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
#include "color.inc"
C  Other local variables for QQBQQI
      integer i,k,mu,v2,v3, isig,isig1,isig3
      integer imin,nlo,nd
C
C  Variables for the main part of the program
      real*8 pbar(0:3,4+max_v),qbar(0:4),uucc,uuss,ddcc,ddss,udsc,ducs,
     1       p(0:3,4+max_v),q(0:4),p21(0:4),p43(0:4),p65(0:4),
     2       eps(0:3,0:2),p21g(0:4),p43g(0:4), res(6,2), resv(6,2),
     3       l21, l43
      double complex scres(6,2)
      real*8 p87(0:4),ph(0:4)
      real*8 uucc_c(3),uuss_c(3),ddcc_c(3),ddss_c(3),udsc_c(3),ducs_c(3)
      integer sign(8), gsign
      complex*16 psi(2,-1:1,4), braket(2,-1:1,4,0:2), 
     1        j21(0:3,-1:1), j43(0:3,-1:1), j21e43(-1:1,-1:1,0:2),
     2        jh1(0:3,-1:1), jh2(0:3,-1:1), e21(0:3,-1:1,0:2),
     3        e43(0:3,-1:1,0:2), e21j43(-1:1,-1:1,0:2)
      complex*16 mm21(6,4,0:2), mm43(6,4,0:2),mv21(6,4,0:2),
     1     mv43(6,4,0:2)
      complex*16 virt(4),born(2),mbox21(-1:1,-1:1,0:2),
     1     mbox43(-1:1,-1:1,0:2)
      complex*16 eps1(0:3,0:2)
      real*8 pk(0:4,4)
      complex*16 prop21(2:4), prop43(2:4), prop21g(2:4), prop43g(2:4)
      real*8 fac,colfac,ratio
      double precision lnmusqp1p2, lnmusqp3p4, lnmusqp1pg, lnmusqp2pg,
     &                 lnmusqp3pg, lnmusqp4pg

      real*8 qvec(0:3)          !q_mu from dips subr

      logical lcol
      logical lspcor            !controls spin correlation
      logical lborn,lbox,lbox21,lbox43,lvirt,ltri
      logical ldebug 
      parameter(ldebug = .false.)
     
      double precision dotrr
      external dotrr

C  Reset the coupling factors
      do i = 1,6
         do isig = 1,4
            fcpl(isig,i) = 0d0
         enddo
      enddo


c determine the Yukawa coupling Higgs bb from BR(H--->bb)
c division by 3 takes into account the color factor 3 for H--->qq
c      
c
c
      do i = 1,6
         do isig = 1,4
            isig1 = sig1(isig,(i+3)/4)
            if ( isig1.ne.0 ) then
               isig3 = sig3(isig,(i+3)/4)
               v2 = v2ini(i)
               v3 = v3ini(i)
               fcpl(isig,i) = clr(fl1(i),v2,isig1)*
     &              clr(fl3(i),v3,isig3)*
     &              b(6,v2,v3)*xmw
            endif
         enddo
      enddo

      return
              
C  Entry for H production
! put entry with nd as argument first to fix ifort compilation error
      entry qqHqqj_c_virt(nlo,pbar,sign,qbar,gsign,nd,
     1               uucc_c,uuss_c,ddcc_c,ddss_c,udsc_c,ducs_c)
      imin = 1
      if(nlo.eq.0) then
         lborn = .true.
         lvirt = .false.
         lbox = .false.
         ltri = .false. 
      elseif(nlo.eq.1) then
         lborn = .true.
         lvirt = .true.
         lbox  = .true.
         ltri = .true. 
         imin = 0
      elseif(nlo.eq.2) then
         lborn = .true.
         lvirt = .true.
         lbox = .false. 
         ltri = .true. 
      else
         lborn = .false.
         lvirt = .true.
         lbox = .true.  
         ltri = .false.
         imin = 0
      endif
      lcol = .true.
      lspcor = .false.
      goto 100


      entry qqHqqj_virt(nlo,pbar,sign,qbar,gsign,
     1             uucc,uuss,ddcc,ddss,udsc,ducs)  
      imin = 1
      nd=0
      if(nlo.eq.0) then         !born
         lborn = .true.
         lvirt = .false.
         lbox = .false.
         ltri = .false. 
      elseif(nlo.eq.1) then     !born+box+tri
         lborn = .true.
         lvirt = .true.
         lbox  = .true.
         ltri = .true.
         imin = 0               !do gauge inv.chk.
c         print*,'new event'
      elseif(nlo.eq.2) then     !born+tri
         lborn = .true.
         lvirt = .true.
         lbox = .false.   
         ltri = .true. 
      else                      !box
         lborn = .false.
         lvirt = .true.
         lbox = .true.
         ltri = .false.
         imin = 0               !do gauge inv.chk.
      endif
      lcol = .false.
      lspcor =.false.
      goto 100

      entry qqHqqj_spcor(qvec,pbar,sign,qbar,gsign,nd,
     1               uucc_c,uuss_c,ddcc_c,ddss_c,udsc_c,ducs_c)
c      print*,'inside qqhqqj_spcor'
      imin = 0
      lcol = .true.
      lspcor = .true.
      lborn = .true.            !just compute born
      lvirt = .false.
      lbox  = .false.
      ltri  =.false.
c      print*,'inside qqhqqj_spcor'
 100  continue
c      print*,'new event'
c
C  Define the internal momenta
      do mu = 0,3
         do i = 1,4+max(n_v,2)
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
         q(mu) = gsign*qbar(mu)

         p21(mu) = p(mu,1)-p(mu,2)
         p43(mu) = p(mu,3)-p(mu,4)
         p21g(mu) = p21(mu)-q(mu)
         p43g(mu) = p43(mu)-q(mu)
         p65(mu) = p(mu,5)-p(mu,6)
         if (n_v.eq.4) then
            p87(mu) = p(mu,8) - p(mu,7)
            ph(mu) = p65(mu) + p87(mu)
         elseif ((n_v.eq.2).or.(n_v.eq.1)) then
            ph(mu) = p65(mu)
         endif
   
      enddo
c
      q(4) = 0
      p21(4) = p21(0)**2-p21(1)**2-p21(2)**2-p21(3)**2
      p43(4) = p43(0)**2-p43(1)**2-p43(2)**2-p43(3)**2
      p21g(4) = p21g(0)**2-p21g(1)**2-p21g(2)**2-p21g(3)**2
      p43g(4) = p43g(0)**2-p43g(1)**2-p43g(2)**2-p43g(3)**2
      p65(4) = p65(0)**2-p65(1)**2-p65(2)**2-p65(3)**2
      if (n_v.eq.4) p87(4) = 
     1         p87(0)**2 - p87(1)**2 - p87(2)**2 - p87(3)**2
      ph(4) = ph(0)**2 - ph(1)**2 - ph(2)**2 - ph(3)**2 

c         
C  Get the vector boson propagator factors
c
      if (p21(4).le.0d0) then
         prop21(2) = 1/(p21(4)-xm2(2))
         prop21(3) = 1/(p21(4)-xm2(3))
c      
      else
         prop21(2) = 1/dcmplx(p21(4)-xm2(2),xmg(2))
         prop21(3) = 1/dcmplx(p21(4)-xm2(3),xmg(3))
c
      endif
      prop21(4) = prop21(3)
c
      if (p43(4).le.0d0) then
         prop43(2) = 1/(p43(4)-xm2(2))
         prop43(3) = 1/(p43(4)-xm2(3))
      else
         prop43(2) = 1/dcmplx(p43(4)-xm2(2),xmg(2))
         prop43(3) = 1/dcmplx(p43(4)-xm2(3),xmg(3))
      endif
      prop43(4) = prop43(3)
      if (p21g(4).le.0d0) then
         prop21g(2) = 1/(p21g(4)-xm2(2))
         prop21g(3) = 1/(p21g(4)-xm2(3))
      else
         prop21g(2) = 1/dcmplx(p21g(4)-xm2(2),xmg(2))
         prop21g(3) = 1/dcmplx(p21g(4)-xm2(3),xmg(3))
      endif
      prop21g(4) = prop21g(3)
      if (p43g(4).le.0d0) then
         prop43g(2) = 1/(p43g(4)-xm2(2))
         prop43g(3) = 1/(p43g(4)-xm2(3))
      else
         prop43g(2) = 1/dcmplx(p43g(4)-xm2(2),xmg(2))
         prop43g(3) = 1/dcmplx(p43g(4)-xm2(3),xmg(3))
   
      endif
      prop43g(4) = prop43g(3)
C
C  Get the external spinors (including factor sqrt(2E) )
      call psi0m(4,pbar(0,1),sign(1),psi)

C  Get the gluon polarization vector and the gluon emission spinors
c  or qvec(mu) 
      do i = imin,2
         if(lspcor) then
            do mu = 0,3
               eps(mu,i) = qvec(mu) ! eps_mu is replaced by q_mu
            enddo
         else
            if(i.eq.0) then
               do mu =0,3
                  eps(mu,i) = q(mu)/q(0) !for gauge inv.check
               enddo
            else 
               call polvec(qbar,i,eps(0,i)) 
            endif 
         endif
c
c     fill complex*16 eps1(0,i) needed for boxline
c     Is it better to just hard code this into boxline?
         do mu=0,3
            eps1(mu,i) = dcmplx(eps(mu,i),0d0)
         enddo
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
c  The calculation of one-loop qcd corrections will be done here.
      call curr(1,psi(1,-1,2),psi(1,-1,1),j21)
      call curr(1,psi(1,-1,4),psi(1,-1,3),j43)
      do i = imin,2
         call curr(1,psi(1,-1,2),braket(1,-1,1,i),jh1)
         call curr(1,braket(1,-1,2,i),psi(1,-1,1),jh2)
         do isig=-1,1,2
            do mu = 0,3
               e21(mu,isig,i) = jh1(mu,isig) + jh2(mu,isig)
            enddo
         enddo

         call curr(1,psi(1,-1,4),braket(1,-1,3,i),jh1)
         call curr(1,braket(1,-1,4,i),psi(1,-1,3),jh2)
         do isig = -1,1,2
            do mu = 0,3
               e43(mu,isig,i) = jh1(mu,isig) + jh2(mu,isig)
            enddo
         enddo
      enddo
c
c     initialize
c
      do i=0,2
         do isig1=-1,1,2
            do isig3=-1,1,2
               mbox43(isig1,isig3,i)=(0d0,0d0)
               mbox21(isig1,isig3,i)=(0d0,0d0)
            enddo
         enddo
      enddo

#ifdef WITH_NLO
      ! set scale
      if (ldoblha) mursq(1,1) = blha_scale**2
      if(lbox) then
         lbox21 = .true.
         lbox43 = .true.
         if (gsign.eq.-1) then
            if (sign(1).eq.-sign(2) .and. 
     1           sign(3).eq.sign(4) ) then
               lbox43 = .false. !  initial gluon attached to 1-2 line only
            elseif (sign(1).eq.sign(2) .and. 
     1              sign(3).eq.-sign(4) ) then
               lbox21 = .false.  !  initial gluon attached to 3-4 line only
            endif
         endif
c     boxline contributions 
c    
c     1-loop corrections to 21 line with gluon emission
c     First, get the Dij's
c
         if(lbox21) then
            if(ldebug) print*,'box on 21 line'
c     print*,'----------------------------'
            call bcd_fill_v(p(0,1),p(0,2),q,p21g)
            if(ldebug) print*,'after bcd_fill_v'
c     
            do i=0,2            ! gluon polarization, 0 is for gauge inv.check
               call curr(1,psi(1,-1,2),braket(1,-1,1,i),jh1) !mborn1
               call curr(1,braket(1,-1,2,i),psi(1,-1,1),jh2) !mborn2
               do isig1 = -1,1,2 ! fermion 1 chirality
                  do isig3 = -1,1,2 ! fermion 3 chirality
c     born-like amplitudes
                     born(1) = jh1(0,isig1)*j43(0,isig3)-
     1                    jh1(1,isig1)*j43(1,isig3)-
     2                    jh1(2,isig1)*j43(2,isig3)-
     3                    jh1(3,isig1)*j43(3,isig3)
                     born(2) = jh2(0,isig1)*j43(0,isig3)-
     1                    jh2(1,isig1)*j43(1,isig3)-
     2                    jh2(2,isig1)*j43(2,isig3)-
     3                    jh2(3,isig1)*j43(3,isig3)
c     
                     if(ldebug) print*,'before boxline'
                     call boxline_vg(psi(1,isig1,1),psi(1,isig1,2),
     1                    p(0,1),p(0,2),isig1,
     2                    eps1(0,i),j43(0,isig3),q,p21g,
     3                    born,virt)
c
                     if(ldebug) print*,'afterboxline'
                     mbox21(isig1,isig3,i) = virt(4)
cvirt(1)+virt(2)
cvirt(4) !(cf-CAc/2)mvirt(1+2) + CAc/2 mvirt(3)
c                     print*,'virt(4)=',virt(4)
                     if(ldebug) print*,'isig1,isig3,i=',isig1,isig3,i
c     do gauge invariance check
c                     if(i.eq.0) then
c                        print*,'isig1,isig3,i=',isig1,isig3,i
c                        print*,'mbox21=', mbox21(isig1,isig3,i)
c                     endif
                  enddo !isig3
               enddo !isig1
            enddo ! gluon pol (i)
         endif !lbox21

c     1-loop corrections to 43 line with gluon emission
c     First, get the Dij's
         if(lbox43) then
            if(ldebug) print*,'box on 43 line'
            call bcd_fill_v(p(0,3),p(0,4),q,p43g)

            do i=0,2               ! gluon polarizations
               call curr(1,psi(1,-1,4),braket(1,-1,3,i),jh1) !mborn1
               call curr(1,braket(1,-1,4,i),psi(1,-1,3),jh2) !mborn2
               do isig1 = -1,1,2   ! fermion 1 chirality
                  do isig3 = -1,1,2 ! fermion 3 chirality
c     born-like amplitudes
                     born(1) = jh1(0,isig3)*j21(0,isig1)-
     1                 jh1(1,isig3)*j21(1,isig1)-jh1(2,isig3)*j21(2,isig1)-
     2                 jh1(3,isig3)*j21(3,isig1)
                     born(2) = jh2(0,isig3)*j21(0,isig1)-
     1                 jh2(1,isig3)*j21(1,isig1)-jh2(2,isig3)*j21(2,isig1)-
     2                 jh2(3,isig3)*j21(3,isig1)

                     call boxline_vg(psi(1,isig3,3),psi(1,isig3,4),
     1                    p(0,3),p(0,4),isig3,
     2                    eps1(0,i),j21(0,isig1),q,p43g,
     3                    born,virt)
                     mbox43(isig1,isig3,i) = virt(4)

                     if(ldebug) print*,'isig1,isig3,i=',isig1,isig3,i
                  enddo !isig3
               enddo !isig1
            enddo !gluon pol (i)
         endif !lbox43

      endif !lbox                     !do not compute c0,b0,d0's
#endif /* WITH_NLO */

C  get the dot products of the currents for the 4 helicity combinations
      do isig1 = -1,1,2
         do isig3 = -1,1,2
            do i = imin,2
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
c     Done for both born and virtual graphs
c
      colfac = cf               !color factors 
c      print*,'colfac=',colfac
c      print*,'cvirt=',cvirt
      do k = 1,6
         do isig = 1,4
            isig1 = sig1(isig,(k+3)/4)
            if (isig1.ne.0) then
               isig3 = sig3(isig,(k+3)/4)
               v2 = v2ini(k)
               v3 = v3ini(k)
               do i = imin,2
c
                  mv21(k,isig,i)=0.0d0
                  mv43(k,isig,i)=0.0d0
c     born 
                  mm21(k,isig,i) = fcpl(isig,k)*prop21g(v2)*
     $                 prop43(v3)*e21j43(isig1,isig3,i) !born term 21 line
c     
                  mm43(k,isig,i) = fcpl(isig,k)*prop21(v2)*
     $                 prop43g(v3)*j21e43(isig1,isig3,i) !born term 43 line
c     
c     boxes
#ifdef WITH_NLO
                  if(lbox21) then !check that the correct factors are here
                     mv21(k,isig,i) = fcpl(isig,k)*prop21g(v2)*
     1                    prop43(v3)*mbox21(isig1,isig3,i)+ !21box
     2                    colfac*cvirtH3j*mm21(k,isig,i)

                     if(i.gt.0) then !gauge invariance testing
                        if(ldebug) then
                           print*,'-----------------------'
                           print*,'mv21q/mm21(',k,isig,0,')=',
     $                          abs(mv21(k,isig,0)/mm21(k,isig,i))
                           print*,'mm21(',k,isig,i,')=',mm21(k,isig,i)
                        endif
                        ratio = abs(mv21(k,isig,0)/mm21(k,isig,i))
                        if(ratio.gt.0.1d0) then ! use constant factor
                           mv21(k,isig,i) = colfac*cvirtH3j*mm21(k,isig,i)
                        endif
                     endif
                  endif
                  if(lbox43) then
                     mv43(k,isig,i) = fcpl(isig,k)*prop21(v2)*
     1                    prop43g(v3)*mbox43(isig1,isig3,i)+ !43box
     2                    colfac*cvirtH3j*mm43(k,isig,i)

                     if(i.gt.0) then ! gauge invariance testing
                        if(ldebug) then
                           print*,'-----------------------'
                           print*,'mv43q/mm43(',k,isig,0,')=',
     $                          abs(mv43(k,isig,0)/mm43(k,isig,i))
                           print*,'mm43(',k,isig,i,')=',mm43(k,isig,i)
                        endif
                        ratio = abs(mv43(k,isig,0)/mm43(k,isig,i))
                        if(ratio.gt.0.1d0) then ! use constant factor           
                           mv43(k,isig,i) = colfac*cvirtH3j*mm43(k,isig,i)
                        endif
                     endif
                  endif
c     triangles
                  if(ltri) then
                     mv21(k,isig,i) = mv21(k,isig,i) + 
     1                    cf*cvirtH3j*mm21(k,isig,i) !21 line 43 tri
                     mv43(k,isig,i) = mv43(k,isig,i) +
     1                    cf*cvirtH3j*mm43(k,isig,i) !41 line 21 tri
                  endif
#endif /* WITH_NLO */

               enddo
            else
               do i = 1,2
                  mm21(k,isig,i) = 0
                  mm43(k,isig,i) = 0
                  mv21(k,isig,i) = 0
                  mv43(k,isig,i) = 0
               enddo
            endif
         enddo
      enddo
      
      if (n_v.gt.1) then
C
C  Now sum the contributions from the 4 helicity combinations for all 
C  subprocesses and multiply by the factor corresponding to
C
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

c     needed for spin correlation
      if(lspcor) then
         fac = fac * 6          ! 12 = C_2*9 is the color factor times 1/2 
                                ! to compensate for over computing q_mu M^mu
c     Its probably more efficient just to compute i = 1
      else
        if (ldoblha) then
          fac = fac * blha_CF*blha_CA**2
        else
          fac = fac * 12         ! 12 = C_2*9 is the color factor
        endif
      endif
c     
c for initial gluon, i.e. gsign=-1, eliminate diagrams with V-->q qbar decay
c     All alpha_s have the same renormalization scale!
c
      l21 = fpi*als(1,1)                 !*als(gnlo) is in dipolesub.f
      l43 = fpi*als(2,1)                 !*als(gnlo) is in dipolesub.f
      if(nd.gt.0) then
          l21 = fpi*als(1,nd)                 !*als(gnlo) is in dipolesub.f
          l43 = fpi*als(2,nd)                 !*als(gnlo) is in dipolesub.f       
      endif
      
      
c      if(LOisON) then           !needed LO cross section 
c         l21 = fpi*als(0)
c         l43 = fpi*als(0)
c       endif
c     
c      print*,'l21,l43',l21,l43
      if (gsign.eq.-1) then
         if (sign(1).eq.-sign(2) .and. 
     1       sign(3).eq.sign(4) ) then
            l43 = 0            !  initial gluon attached to 1-2 line only
         elseif (sign(1).eq.sign(2) .and. 
     1           sign(3).eq.-sign(4) ) then
            l21 = 0            !  initial gluon attached to 3-4 line only
         endif
      endif
c
c     sum over helicity combinations
c     
c
      do k = 1,6
         res(k,1) = 0
         res(k,2) = 0
         resv(k,1) = 0
         resv(k,2) = 0
         scres(k,1) = 0
         scres(k,2) = 0
         do isig = 1,4
c undo changes !
            do i = 1,2 
               if(lborn) then
                  res(k,1) = res(k,1) + 
     $                 l21*( dreal(mm21(k,isig,i))**2+
     $                 dimag(mm21(k,isig,i))**2) !|mm21|^2     
                  res(k,2) = res(k,2) + 
     $                 l43*( dreal(mm43(k,isig,i))**2+
     $                 dimag(mm43(k,isig,i))**2) !|mm43|^2
                  
               endif

#ifdef WITH_NLO
               if(lvirt) then
                  resv(k,1) = resv(k,1)+
     1                 l21*2d0*dreal(mm21(k,isig,i)*
     $                 conjg(mv21(k,isig,i))) !2Re[mm21 conj(mv21)]  
                  resv(k,2) = resv(k,2) + 
     1                 l43*2d0*dreal(mm43(k,isig,i)*
     $                 conjg(mv43(k,isig,i))) !2Re[mm43 conj(mv43)]
               endif
#endif /* WITH_NLO */
            enddo
            ! sctree amplitudes - convert polvec to BLHA-helvec
            ! <M,-|M,+> = (<M,1|M,1>-<M,2|M,2>-2i*gsign*Re(<M,1|M,2>))/2
            if (ldoblha) then
              scres(k,1) = scres(k,1) + l21*fac*(
     &            abs(mm21(k,isig,1))**2
     &          - abs(mm21(k,isig,2))**2
     &          - 2d0*dcmplx(0d0,1d0)*gsign*dreal(
     &            dconjg(mm21(k,isig,1))*mm21(k,isig,2)))/2d0
              scres(k,2) = scres(k,2) + l43*fac*(
     &            abs(mm43(k,isig,1))**2
     &          - abs(mm43(k,isig,2))**2
     &          - 2d0*dcmplx(0d0,1d0)*gsign*dreal(
     &            dconjg(mm43(k,isig,1))*mm43(k,isig,2)))/2d0
            endif
         enddo
      enddo
      
      do k =1,6
c     born**2
         if(lborn) then
            res(k,1) = res(k,1)*fac !21 line
            res(k,2) = res(k,2)*fac !43 line
         endif
c     2Re[mborn mvirt*]
#ifdef WITH_NLO
         if(lvirt) then
            if(nd.gt.0) then
            resv(k,1) = resv(k,1)*fac*als(1,nd)/fpi !21 line
            resv(k,2) = resv(k,2)*fac*als(2,nd)/fpi !43 line 
            else
            resv(k,1) = resv(k,1)*fac*als(1,1)/fpi !21 line
            resv(k,2) = resv(k,2)*fac*als(2,1)/fpi !43 line 
            endif
         endif
#endif /* WITH_NLO */

      enddo
      
      if (ldoblha) then
        if (blha_amptype(blha_curproc).eq.0) then ! tree
          if ( (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 1) .or.
     &         (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 2) ) then
            blha_amp(1)=blha_amp(1)+
     &       res(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),1)
          endif
          if ( (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 1) .or.
     &         (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 3) ) then
            blha_amp(1)=blha_amp(1)+
     &       res(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),2)
          endif

        else if (blha_amptype(blha_curproc).eq.1) then ! loop
          lnmusqp1p2 = log(blha_scale**2/(2d0*dotrr(pbar(0,1),pbar(0,2))))
          lnmusqp3p4 = log(blha_scale**2/(2d0*dotrr(pbar(0,3),pbar(0,4))))
          lnmusqp1pg = log(blha_scale**2/(2d0*dotrr(pbar(0,1),qbar(0))))
          lnmusqp2pg = log(blha_scale**2/(2d0*dotrr(pbar(0,2),qbar(0))))
          lnmusqp3pg = log(blha_scale**2/(2d0*dotrr(pbar(0,3),qbar(0))))
          lnmusqp4pg = log(blha_scale**2/(2d0*dotrr(pbar(0,4),qbar(0))))

          if ( (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 1) .or.
     &         (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 2) ) then
            blha_amp(1)=blha_amp(1)+
     &       ( als(1,1)/(2d0*pi)*(-blha_CA-2*blha_CF)
     &        +als(2,1)/(2d0*pi)*(-2*blha_CF) )*
     &       res(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),1)
            blha_amp(2)=blha_amp(2)+
     &       ( als(1,1)/(2d0*pi)*
     &         (-2*(blha_gammaQuark+blha_CF*lnmusqp1p2)-blha_gammaGluon
     &          +blha_CA*(lnmusqp1p2-lnmusqp1pg-lnmusqp2pg) )
     &        +als(2,1)/(2d0*pi)*
     &         (-2*(blha_gammaQuark+blha_CF*lnmusqp3p4)) )*
     &       res(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),1)
             blha_amp(3)= blha_amp(3)+
     &         resv(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),1)
     &         +( als(1,1)/(2d0*pi)*(
     &              blha_CF*(-cvirtH3j+pi**2/3d0-8d0)-2*blha_tgammaQuark
     &             +(-blha_CA-2*blha_CF)*pi**2/6d0
     &            )
     &           +als(2,1)/(2d0*pi)*(
     &              blha_CF*(-cvirtH3j+pi**2/3d0-8d0)-2*blha_tgammaQuark
     &             +(-2*blha_CF)*pi**2/6d0
     &             +(-2)*blha_gammaQuark*lnmusqp3p4
     &             +(-blha_CF)*lnmusqp3p4**2
     &            )
     &          )*res(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),1)
            blha_amp(4)=blha_amp(4)+
     &       res(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),1)
          endif
          if ( (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 1) .or.
     &         (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 3) ) then
            blha_amp(1)=blha_amp(1)+
     &       ( als(1,1)/(2d0*pi)*(-2*blha_CF)
     &        +als(2,1)/(2d0*pi)*(-blha_CA-2*blha_CF) )*
     &       res(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),2)
            blha_amp(2)=blha_amp(2)+
     &       ( als(1,1)/(2d0*pi)*
     &         (-2*(blha_gammaQuark+blha_CF*lnmusqp1p2))
     &        +als(2,1)/(2d0*pi)*
     &         (-2*(blha_gammaQuark+blha_CF*lnmusqp3p4)-blha_gammaGluon
     &          +blha_CA*(lnmusqp3p4-lnmusqp3pg-lnmusqp4pg) ) )*
     &       res(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),2)
             blha_amp(3)= blha_amp(3)+
     &         resv(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),2)
     &         +( als(1,1)/(2d0*pi)*(
     &              blha_CF*(-cvirtH3j+pi**2/3d0-8d0)-2*blha_tgammaQuark
     &             +(-2*blha_CF)*pi**2/6d0
     &             +(-2)*blha_gammaQuark*lnmusqp1p2
     &             +(-blha_CF)*lnmusqp1p2**2
     &            )
     &           +als(2,1)/(2d0*pi)*(
     &              blha_CF*(-cvirtH3j+pi**2/3d0-8d0)-2*blha_tgammaQuark
     &             +(-blha_CA-2*blha_CF)*pi**2/6d0
     &            )
     &          )*res(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),2)
            blha_amp(4)=blha_amp(4)+
     &       res(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),2)
          endif

        else if (blha_amptype(blha_curproc).eq.2) then ! cctree
          if ( (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 1) .or.
     &         (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 2) ) then
          ! gluon attached to 1-2 line
            ! 1-2
            call BLHA_cctree(1,2,
     &         (blha_CA-2*blha_CF)/2d0*
     &         res(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),1)
     &      )
            ! 1-5
            call BLHA_cctree(1,5,
     &         -blha_CA/2d0*
     &         res(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),1)
     &      )
            ! 2-5
            call BLHA_cctree(2,5,
     &         -blha_CA/2d0*
     &         res(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),1)
     &      )
            ! 3-4
            call BLHA_cctree(3,4,
     &        -blha_CF*
     &         res(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),1)
     &      )
          endif
          if ( (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 1) .or.
     &         (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 3) ) then
          ! gluon attached to 3-4 line
            ! 1-2
            call BLHA_cctree(1,2,
     &        -blha_CF*
     &         res(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),2)
     &      )
            ! 3-4
            call BLHA_cctree(3,4,
     &         (blha_CA-2*blha_CF)/2d0*
     &         res(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),2)
     &      )
            ! 3-5
            call BLHA_cctree(3,5,
     &         -blha_CA/2d0*
     &         res(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),2)
     &      )
            ! 4-5
            call BLHA_cctree(4,5,
     &         -blha_CA/2d0*
     &         res(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),2)
     &      )
          endif

        else if (blha_amptype(blha_curproc).eq.3) then ! sctree
          if ( (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 1) .or.
     &         (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 2) ) then
          ! gluon attached to 1-2 line
            ! 5-1 = <5|15|5>
            call BLHA_sctree(5,1,
     &         -blha_CA/2d0*
     &         scres(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),1)
     &      )
            ! 5-2 = <5|25|5>
            call BLHA_sctree(5,2,
     &         -blha_CA/2d0*
     &         scres(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),1)
     &      )
            ! 5-5 = CA<5|5>
            call BLHA_sctree(5,5,
     &         blha_CA*
     &         scres(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),1)
     &      )
          endif
          if ( (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 1) .or.
     &         (blha_idsubproc(blha_cursubproc,blha_curproc)/10 .eq. 3) ) then
          ! gluon attached to 3-4 line
            ! 5-3 = <5|35|5>
            call BLHA_sctree(5,3,
     &         -blha_CA/2d0*
     &         scres(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),2)
     &      )
            ! 5-4 = <5|45|5>
            call BLHA_sctree(5,4,
     &         -blha_CA/2d0*
     &         scres(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),2)
     &      )
            ! 5-5 = CA<5|5>
            call BLHA_sctree(5,5,
     &         blha_CA*
     &         scres(mod(blha_idsubproc(blha_cursubproc,blha_curproc),10),2)
     &      )
          endif

        else 
          call BLHA_amptypeerror(blha_amptype(blha_curproc),__FILE__,__LINE__)
        endif
      endif
c     
c    
C     Assign the final results for the subprocesses
      if (lcol) then
         do i = 1,2
            uucc_c(i+1) = res(1,i) + resv(1,i)
            uuss_c(i+1) = res(2,i) + resv(2,i)
            ddcc_c(i+1) = res(3,i) + resv(3,i)
            ddss_c(i+1) = res(4,i) + resv(4,i)
            udsc_c(i+1) = res(5,i) + resv(5,i)
            ducs_c(i+1) = res(6,i) + resv(6,i)
         enddo
         uucc_c(1) = res(1,1) +  res(1,2) + resv(1,1) + resv(1,2)
         uuss_c(1) = res(2,1) +  res(2,2) + resv(2,1) + resv(2,2)
         ddcc_c(1) = res(3,1) +  res(3,2) + resv(3,1) + resv(3,2)
         ddss_c(1) = res(4,1) +  res(4,2) + resv(4,1) + resv(4,2)
         udsc_c(1) = res(5,1) +  res(5,2) + resv(5,1) + resv(5,2)
         ducs_c(1) = res(6,1) +  res(6,2) + resv(6,1) + resv(6,2)
      else
         uucc = res(1,1) +  res(1,2) + resv(1,1) +  resv(1,2)
         uuss = res(2,1) +  res(2,2) + resv(2,1) +  resv(2,2)
         ddcc = res(3,1) +  res(3,2) + resv(3,1) +  resv(3,2)
         ddss = res(4,1) +  res(4,2) + resv(4,1) +  resv(4,2)
         udsc = res(5,1) +  res(5,2) + resv(5,1) +  resv(5,2)
         ducs = res(6,1) +  res(6,2) + resv(6,1) +  resv(6,2)
      endif
      
      return
      end
c------------------------------------------------------------------------

