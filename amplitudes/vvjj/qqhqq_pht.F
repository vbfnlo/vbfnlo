*******************************************************************************
*******************************************************************************
*** sophy: 27-10-2009:  NEW FILE                                            ***
*** code for real photon emissions, originally written by Terrance Figy     ***
*** All alterations/comments marked with *S or *sophy                       ***
*******************************************************************************
*******************************************************************************

C********************  qqhqqa ***********************************
C
      subroutine qqhqqPHTi 
      implicit none               

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/mssm.inc"

      double precision  fpi, xmb
      parameter ( fpi=4d0*Pi, xmb=4.5d0 )

      double precision p87(0:4), ph(0:4)
C
C  QQHQQ calculates the matrix elements**2 for light Higgs production by
C  electroweak boson fusion in quark quark scattering
C
C        q1 q3    ---->   q2 q4 A H
C
C  QQHQQI must be called first to initialize some couplings
C
C  The main task of QQHQQI is to set up the products of coupling constants
C  needed in Feynman graphs a, ..., g which are stored in
      real*8 fcpl(4,6)
C  and which are calculated from output of KOPPLN
      real*8 clr,xm2,xmg,b,v,a
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)
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
      double complex fcplCT(4,6)

      data fl1 /3,3,4,4,3,4/
      data fl3 /3,4,3,4,4,3/

      save fl1,fl3
      common /cqqhqq/ fcplCT, fcpl, sig1,sig3, v2ini,v3ini
c note: initialization for v2ini, v3ini, sig1,sig3 in block data
C
c      data v2ini /2,2,2,2,3,4/, v3ini /2,2,2,2,4,3/
c      data sig1 /-1,-1, 1, 1,-1,3*0/
c      data sig3 /-1, 1,-1, 1,-1,3*0/
c     
c     relative charge for fermions in mutilplies of 3
c
      integer Qf1(6),Qf2(6),Qf3(6),Qf4(6)
      common/icharge/Qf1,Qf2,Qf3,Qf4
      data Qf1 /2,2,-1,-1,2,-1/
      data Qf2 /2,2,-1,-1,-1,2/
      data Qf3 /2,-1,2,-1,-1,2/
      data Qf4 /2,-1,2,-1,2,-1/

c      save Qf1,Qf2,Qf3,Qf4
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"

C  Other local variables for QQBQQI
      integer i,k,mu,v2,v3, isig,isig1,isig3
C
C  Variables for the main part of the program
      real*8 pbar(0:3,4+max_v),qbar(0:4),uucc,uuss,ddcc,ddss,udsc,ducs,
     1       p(0:3,4+max_v),q(0:4),p21(0:4),p43(0:4),p65(0:4),eps(0:3,2),
     2       p21g(0:4),p43g(0:4), res(6)
      real*8 uucc_c(3),uuss_c(3),ddcc_c(3),ddss_c(3),udsc_c(3),ducs_c(3)
      integer fsign(4+max_v),gsign
      complex*16 psi(2,-1:1,4), braket(2,-1:1,4,2), 
     1        j21(0:3,-1:1), j43(0:3,-1:1), j21e43(-1:1,-1:1,2),
     2        jh1(0:3,-1:1), jh2(0:3,-1:1), e21(0:3,-1:1,2),
     3        e43(0:3,-1:1,2), e21j43(-1:1,-1:1,2),
     4        mm21(6,4,2), mm43(6,4,2)
      complex*16 f21(0:3,-1:1,2),f43(0:3,-1:1,2),f21j43(-1:1,-1:1,2),
     $     j21f43(-1:1,-1:1,2),j21j43a(-1:1,-1:1,2,2),amp(6,4,2)

      real*8 pk(0:4,4)
      complex*16 prop21(1:4), prop43(1:4), prop21g(1:4), prop43g(1:4)


*S sophy: new output of call to formfactor calculation


      real*8 xmw, fac
      logical lcol
      double complex  czero
      parameter ( czero = (0.0d0,0.0d0) )
c     PS slicing parameters
      logical lPSslice
      parameter (lPSslice=.false.)
      double precision deltaE,deltaTH
      common/PSslice/deltaE,deltaTH
      double precision theta(4)
      double precision angle
      external angle
     
c CPS scheme
      double precision qgammaq
      common /VBFNLO_HIGGS_CPS/ qgammaq

      if(lPSslice) then
c     PS slicing parameters
      deltaE = 1.0d0            ! 1 GeV
      deltaTH = 1.0d-3          ! 
      endif
      
C  Reset the coupling factors
      do i = 1,6
         do isig = 1,4
            fcpl(isig,i) = 0.d0
         enddo
      enddo

c
c determine the Yukawa coupling Higgs bb from BR(H--->bb)
c division by 3 takes into account the color factor 3 for H--->qq
c
      xmw = sqrt(xm2(3))

      do i = 1,6
         do isig = 1,4
            isig1 = sig1(isig,(i+3)/4)
            if ( isig1.ne.0 ) then
               isig3 = sig3(isig,(i+3)/4)
               v2 = v2ini(i)
               v3 = v3ini(i)
                  fcpl(isig,i) = clr(fl1(i),v2,isig1)*
     &                           clr(fl3(i),v3,isig3)

c                      b(6,v2,v3)
            endif
         enddo
      enddo

      print*,'Qf1=',Qf1
      print*,'Qf2=',Qf2
      print*,'Qf3=',Qf3
      print*,'Qf4=',Qf4
c      read(*,*)
      return
C  Entry for H production
      entry qqHqqPHT(pbar,fsign,qbar,gsign,
     1             uucc,uuss,ddcc,ddss,udsc,ducs)
c      print*,'nlosusy=',nlosusy
      lcol = .false.
      goto 100

      entry qqHqqPHT_c(pbar,fsign,qbar,gsign,
     1               uucc_c,uuss_c,ddcc_c,ddss_c,udsc_c,ducs_c)
c      print*,'nlosusy=',nlosusy
      lcol = .true.
 100  continue
             

c     print*,'fcpl=',fcpl
c
c     Initialize res to zero 
      do k = 1,6
         do isig=1,4
            do i=1,2
               mm21(k,isig,i) = czero
               mm43(k,isig,i) = czero
               amp(k,isig,i) = czero 
            enddo
         enddo
         res(k)=0.0d0 
      enddo
      if(lPSslice) then
c
c     Apply phase space slicing 
c     Energy of photon

      if(qbar(0).lt.deltaE) then
c         print*,'event fails PS slice'
         goto 200
      endif
c     The angles
      theta(1)=angle(pbar(0,1),qbar)
      theta(2)=angle(pbar(0,2),qbar)
      theta(3)=angle(pbar(0,3),qbar)
      theta(4)=angle(pbar(0,4),qbar)
      
      do i = 1,4
         if(theta(i).gt.1.0d0-deltaTH) then
c            print*,'event fails PS slice'
            goto 200
         endif
      enddo
c     
c      print*,'event passes PS slice'
      endif

C  Define the internal momenta
      do mu = 0,3
         do i = 1,4+n_v
            p(mu,i) = pbar(mu,i)*fsign(i)
         enddo
         q(mu) = gsign*qbar(mu)
      enddo
      q(4) = 0

C
C  Get the external spinors (including factor sqrt(2E) )
      call psi0m(4,pbar(0,1),fsign(1),psi)

C  Get the photon polarization vector and the photon emission spinors
      do i = 1,2
         call polvec(qbar,i,eps(0,i))
c     gauge invariance check:
c         do mu = 0,3
c            eps(mu,i) = qbar(mu)/qbar(0)
c         enddo
c
         do isig = -1,1,2
            call ket2r(psi(1,isig,1),.true.,p(0,1),isig,q,eps(0,i),
     1           braket(1,isig,1,i),pk(0,1))
            call bra2r(psi(1,isig,2),.true.,p(0,2),isig,q,eps(0,i),
     1           braket(1,isig,2,i),pk(0,2))
            call ket2r(psi(1,isig,3),.true.,p(0,3),isig,q,eps(0,i),
     1           braket(1,isig,3,i),pk(0,3))
            call bra2r(psi(1,isig,4),.true.,p(0,4),isig,q,eps(0,i),
     1           braket(1,isig,4,i),pk(0,4))
*            write(*,*)'  '
*            write(*,*)'i, isig =', i, isig
*            write(*,*)'psi(1,isig,4) =', psi(1,isig,4)
*            write(*,*)'eps(0,i) =', eps(0,i)
*            write(*,*)'p(0,4) =', p(0,4)
*            write(*,*)'pk(0,4) =', pk(0,4)
*            write(*,*)'braket(1,isig,4,i) =', braket(1,isig,4,i)
*            write(*,*)'q =', q
         enddo
      enddo

C  get the f-fbar currents J21^mu, J43^mu, E21^mu, E43^mu, F21^mu, F43^mu
      call curr(1,psi(1,-1,2),psi(1,-1,1),j21)
      call curr(1,psi(1,-1,4),psi(1,-1,3),j43)
      do i = 1,2
         call curr(1,psi(1,-1,2),braket(1,-1,1,i),jh1)
         call curr(1,braket(1,-1,2,i),psi(1,-1,1),jh2)
*         write(*,*)'  '
*         write(*,*)'i, braket(1,-1,1,i) =', i, braket(1,-1,1,i)
*         write(*,*)'i, braket(1,-1,2,i) =', i, braket(1,-1,2,i)
*         write(*,*)'i, psi(1,-1,1) =', i, psi(1,-1,1)
*         write(*,*)'i, psi(1,-1,2) =', i, psi(1,-1,2)
         do mu = 0,3
            do isig = -1,1,2
               e21(mu,isig,i) = jh1(mu,isig) + jh2(mu,isig)
               f21(mu,isig,i) = jh1(mu,isig) - jh2(mu,isig)
*               write(*,*)' '
*               write(*,*)'i, mu, isig =', i, mu, isig
*               write(*,*)'e21 =', e21(mu,isig,i)
*               write(*,*)'f21 =', f21(mu,isig,i)
*               write(*,*)'jh1 =', jh1(mu,isig)
*               write(*,*)'jh2 =', jh1(mu,isig)
            enddo
         enddo

         call curr(1,psi(1,-1,4),braket(1,-1,3,i),jh1)
         call curr(1,braket(1,-1,4,i),psi(1,-1,3),jh2)
         do mu = 0,3
            do isig = -1,1,2
               e43(mu,isig,i) = jh1(mu,isig) + jh2(mu,isig)
               f43(mu,isig,i) = jh1(mu,isig) - jh2(mu,isig)
*               write(*,*)' '
*               write(*,*)'i, mu, isig =', i, mu, isig
*               write(*,*)'e43 =', e43(mu,isig,i)
*               write(*,*)'f43 =', f43(mu,isig,i)
            enddo
         enddo
      enddo
      
      do mu = 0,3
         p21(mu) = p(mu,1)-p(mu,2)
         p43(mu) = p(mu,3)-p(mu,4)
         p21g(mu) = p21(mu)-q(mu)
         p43g(mu) = p43(mu)-q(mu)
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
c      prop21(1) = 1.0d0/dcmplx(p21(4),0.0d0)
c      print*,'p21(1)=',p21(1)

      if (p43(4).le.0d0) then
         prop43(2) = 1/(p43(4)-xm2(2))
         prop43(3) = 1/(p43(4)-xm2(3))

      else
         prop43(2) = 1/dcmplx(p43(4)-xm2(2),xmg(2))
         prop43(3) = 1/dcmplx(p43(4)-xm2(3),xmg(3))
 
      endif
      prop43(4) = prop43(3)
c      prop43(1) = 1.0d0/dcmplx(p43(4),0.0d0)
c      print*,'p43(1)=',p43(1)

      if (p21g(4).le.0d0) then
         prop21g(2) = 1/(p21g(4)-xm2(2))
         prop21g(3) = 1/(p21g(4)-xm2(3))
    
      else
         prop21g(2) = 1/dcmplx(p21g(4)-xm2(2),xmg(2))
         prop21g(3) = 1/dcmplx(p21g(4)-xm2(3),xmg(3))
     
      endif
      prop21g(4) = prop21g(3)
c      prop21g(1) = 1.0d0/dcmplx(p21g(4),0.0d0)
c      print*,'p21g(1)=',p21g(1)
      if (p43g(4).le.0d0) then
         prop43g(2) = 1/(p43g(4)-xm2(2))
         prop43g(3) = 1/(p43g(4)-xm2(3))
     
      else
         prop43g(2) = 1/dcmplx(p43g(4)-xm2(2),xmg(2))
         prop43g(3) = 1/dcmplx(p43g(4)-xm2(3),xmg(3))
   
      endif
      prop43g(4) = prop43g(3)
c      prop43g(1) = 1.0d0/dcmplx(p43g(4),0.0d0)
c      print*,'p43g(1)=',p43g(1) 


c We always have ZZ and WW fusion graphs
c Build the T^mu^nu J43^mu J21^nu 
c      

c$$$*sophy: I am fairly certain these calls to form_factors are now totally
c$$$*       irrelevant as we take our couplings from koppln.F .  I'll leave them
c$$$*       in for now just to be safe 
c$$$*S     call form_factors(born21(3),a1factor21(3),a2factor21(3),
c$$$*S    $     a3factor21(3),born21(2),a1factor21(2),a2factor21(2),
c$$$*S    $     a3factor21(2),p21g,p43,0,0,HiggsType) 
c$$$#ifdef WITH_LT
c$$$      call form_factors(p21g,p43,factor,0)
c$$$#endif
c$$$
c$$$      a1factor21(3) = factor(1,3,4)  !a1W
c$$$      a2factor21(3) = factor(2,3,4)
c$$$      a3factor21(3) = factor(3,3,4)
c$$$
c$$$      a1factor21(2) = factor(1,2,2)  !a1Z
c$$$      a2factor21(2) = factor(2,2,2)
c$$$      a3factor21(2) = factor(3,2,2)
c$$$      
c$$$      a1factor21(4) = a1factor21(3)
c$$$      a2factor21(4) = a2factor21(3)
c$$$      a3factor21(4) = a3factor21(3)
c$$$
c$$$*S      born21(4) = born21(3)
c$$$
c$$$*S     call form_factors(born43(3),a1factor43(3),a2factor43(3),
c$$$*S     $     a3factor43(3),born43(2),a1factor43(2),a2factor43(2),
c$$$*S     $     a3factor43(2),p21,p43g,0,0,HiggsType) 
c$$$#ifdef WITH_LT
c$$$      call form_factors(p21,p43g,factor,0)
c$$$#endif
c$$$
c$$$      a1factor43(3) = factor(1,3,4)  !a1W
c$$$      a2factor43(3) = factor(2,3,4)
c$$$      a3factor43(3) = factor(3,3,4)
c$$$
c$$$      a1factor43(2) = factor(1,2,2)  !a1Z
c$$$      a2factor43(2) = factor(2,2,2)
c$$$      a3factor43(2) = factor(3,2,2)
c$$$
c$$$      a1factor43(4) = a1factor43(3)
c$$$      a2factor43(4) = a2factor43(3)
c$$$      a3factor43(4) = a3factor43(3)
c$$$*S      born43(4) = born43(3)
c
      do isig1 = -1,1,2
         do isig3 = -1,1,2
            do i=1,2
               e21j43(isig1,isig3,i) = e21(0,isig1,i)*j43(0,isig3) -
     &              e21(1,isig1,i)*j43(1,isig3) -
     &              e21(2,isig1,i)*j43(2,isig3) -
     &              e21(3,isig1,i)*j43(3,isig3)
c     
               j21e43(isig1,isig3,i) = j21(0,isig1)*e43(0,isig3,i) -
     &              j21(1,isig1)*e43(1,isig3,i) -
     &              j21(2,isig1)*e43(2,isig3,i) -
     &              j21(3,isig1)*e43(3,isig3,i)
               
               f21j43(isig1,isig3,i) = f21(0,isig1,i)*j43(0,isig3) -
     &              f21(1,isig1,i)*j43(1,isig3) -
     &              f21(2,isig1,i)*j43(2,isig3) -
     &              f21(3,isig1,i)*j43(3,isig3)
c     
               j21f43(isig1,isig3,i) = j21(0,isig1)*f43(0,isig3,i) -
     &              j21(1,isig1)*f43(1,isig3,i) -
     &              j21(2,isig1)*f43(2,isig3,i) -
     &              j21(3,isig1)*f43(3,isig3,i)

               call calc_tens_wwah(p21(0),p43(0),q(0),j21(0,isig1),
     $              j43(0,isig3),eps(0,i),xm2(3),
     $              j21j43a(isig1,isig3,i,1)) !piece for udsc
c               
               call calc_tens_wwah(p43(0),p21(0),q(0),j43(0,isig3),
     $              j21(0,isig1),eps(0,i),xm2(3),
     $              j21j43a(isig1,isig3,i,2)) !piece for ducs

*               write(*,*)'  '
*               write(*,*)'isig1, isig3, i =', isig1, isig3, i
*               write(*,*)'e21j43 =', e21j43(isig1,isig3,i)
*               write(*,*)'j21e43 =', j21e43(isig1,isig3,i)
*               write(*,*)'f21j43 =', f21j43(isig1,isig3,i)
*               write(*,*)'j21f43 =', j21f43(isig1,isig3,i)
cc                       
            enddo
         enddo
      enddo


C
C  now get the coupling*propagator factors for subprocess k, helicity
C  combination isig

       
      do k = 1,6
c
c         print*,'k=',k
c         print*,'Qf1(k)-Qf2(k)',Qf1(k)-Qf2(k)
c         print*,'Qf1(k)+Qf2(k)',Qf1(k)+Qf2(k)
c         print*,'Qf3(k)-Qf4(k)',Qf3(k)-Qf4(k)
c         print*,'Qf3(k)+Qf4(k)',Qf3(k)+Qf4(k)

c
         do isig = 1,4
            isig1 = sig1(isig,(k+3)/4)
            if (isig1.ne.0) then
               isig3 = sig3(isig,(k+3)/4)
c               print*,'isig1=',isig1
c               print*,'isig3=',isig3
               v2 = v2ini(k)
               v3 = v3ini(k)
               do i = 1,2     
c     higgs coupling is an overall factor here
c     we therefore factor it out
                  mm21(k,isig,i) = fcpl(isig,k)*
     $                 (dble(Qf1(k) + Qf2(k))*e21j43(isig1,isig3,i)+
     $                 dble(Qf1(k) - Qf2(k))*f21j43(isig1,isig3,i))* 
     $                 prop21g(v2)*prop43(v3)/6.0d0
*S                  mm21(k,isig,i) = fcpl(isig,k)*born21(v2)*
*S     $                 (dble(Qf1(k) + Qf2(k))*e21j43(isig1,isig3,i)+
*S     $                 dble(Qf1(k) - Qf2(k))*f21j43(isig1,isig3,i))* 
*S     $                 prop21g(v2)*prop43(v3)/6.0d0     
c
                  mm43(k,isig,i) = fcpl(isig,k)*
     $                 (dble(Qf3(k)+Qf4(k))*j21e43(isig1,isig3,i)+
     $                  dble(Qf3(k)-Qf4(k))*j21f43(isig1,isig3,i))* 
     $                 prop21(v2)*prop43g(v3)/6.0d0
*S                  mm43(k,isig,i) = fcpl(isig,k)*born43(v2)*
*S     $                 (dble(Qf3(k)+Qf4(k))*j21e43(isig1,isig3,i)+
*S     $                  dble(Qf3(k)-Qf4(k))*j21f43(isig1,isig3,i))* 
*S     $                 prop21(v2)*prop43g(v3)/6.0d0

c                  amp(k,isig,i) =  mm21(k,isig,i) + mm43(k,isig,i)
                  
                  if(k.ge.5) then ! additional graphs included for charged current processes
                     amp(k,isig,i) =  mm21(k,isig,i) + mm43(k,isig,i) + 
     $                    fcpl(isig,k)*
     $                    j21j43a(isig1,isig3,i,k-4)*prop21(v2)*
     $                    prop43(v3) 
*S                     amp(k,isig,i) =  mm21(k,isig,i) + mm43(k,isig,i) + 
*S     $                    fcpl(isig,k)*born21(v2)*
*S     $                    j21j43a(isig1,isig3,i,k-4)*prop21(v2)*
*S     $                    prop43(v3) 
                  else
                     amp(k,isig,i) =  mm21(k,isig,i) + mm43(k,isig,i) 
    
                  endif
               enddo
            else
               do i = 1,2
                  amp(k,isig,i) = czero
               enddo
            endif
         enddo
      enddo
c
c     also beware of averaging factors: qq > qq + H
c                                       photon + q > qqbar q + H 
c
c     This is the color summed, spin summed matrix element squared.
c
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
c   Int dq^2/2pi * dOmega/(32 Pi^2) * 1/(q^2-M*2)^2+(MGamma)^2) 
c
c   = 1/(16 Pi MGamma)
c
c (where the propagator factor is fac above)
      fac = fac *16*pi*qgammaq

      fac = fac * EL**2 * 9.0D0    ! 9 is the color factor * positron charge

*      write(*,*)'fac =', fac

      if(gsign.eq.-1 .and. lcol) fac = 2.0D0 * fac
c
      do k = 1,6
         res(k) = 0
         do isig = 1,4
            do i = 1,2
               res(k) = res(k) + 
     1              ( dreal(amp(k,isig,i))**2 + 
     2              dimag(amp(k,isig,i))**2 )
               
            enddo
         enddo
         res(k) = res(k)*fac          
C   *p65(4)/xm2(6)
      enddo

C  Assign the final results for the subprocesses
     
 200  continue

      if (lcol) then
         do i = 1,2
            uucc_c(i+1) = res(1)/2.0d0
            uuss_c(i+1) = res(2)/2.0d0
            ddcc_c(i+1) = res(3)/2.0d0
            ddss_c(i+1) = res(4)/2.0d0
            udsc_c(i+1) = res(5)/2.0d0
            ducs_c(i+1) = res(6)/2.0d0
         enddo
         uucc_c(1) = res(1)
         uuss_c(1) = res(2)
         ddcc_c(1) = res(3)
         ddss_c(1) = res(4)
         udsc_c(1) = res(5)
         ducs_c(1) = res(6)
      else
         uucc = res(1)
         uuss = res(2)
         ddcc = res(3)
         ddss = res(4)
         udsc = res(5)
         ducs = res(6)
      endif
     
*      write(*,*)'uucc =', uucc
*      write(*,*)'uuss =', uuss
*      write(*,*)'ddcc =', ddcc
*      write(*,*)'ddss =', ddss
*      write(*,*)'udsc =', udsc
*      write(*,*)'ducs =', ducs
*      stop


c     do k = 5,6
c     if(res(k).gt. 1.0d-5) then
c         print*,'res=',res
c         read(*,*) 
c        endif
c     enddo
      end
ccccccccccccccc
c This computes the tensor T(k1,k2,k3,eps1,eps2,eps3)
c
c W+(k1,eps1) W-(k2,eps2) -> A(k3,esp3) + h(k4)
c
cccccccccccccccc
      subroutine calc_tens_wwah(k1,k2,k3,e1,e2,e3,mw2,Ctree)
      implicit none
c
      double complex e1(0:3),e2(0:3) ! effective pol vectors
      double precision e3(0:3)         ! photon pol vector
      double precision k1(0:3),k2(0:3),k3(0:3)
      double precision mw2
      integer mu
      double complex Pair6, Pair9, Pair4, Pair12, Pair5, Pair10
      double complex Pair11, Pair8, Pair15, Pair1, Pair13, Pair2
      double complex Pair7, Pair14, Pair3, Abb1, Abb2, Abb7, Abb3
      double complex Abb4, Abb5, Abb6, AbbSum1, AbbSum2
      double complex ctree
      double precision k4(0:3)
c
c    
      double precision dotrr
      external dotrr
c
      double complex dotrc,dotcc
      external dotrc
      external dotcc
      
c
      do mu = 0,3
         k4(mu) = k1(mu) + k2(mu) - k3(mu)
      enddo
c
      Pair13 = dotrr(k1,k3)
      Pair2 = dotrr(k1,k4)
      Pair14 = dotrr(k3,k3)
      Pair3 = dotrr(k4,k4)
      Pair6 = dotcc(e1,e2)
      Pair9 = dotrc(e3,e1)
      Pair4 = dotrc(k1,e1)
      Pair12 = dotrc(k3,e1)
      Pair5 = dotrc(e3,e2)
      Pair10 = dotrc(k2,e2)
      Pair11 = dotrc(k3,e2)
      Pair8 = dotrr(e3,k3)
      Pair15 = dotrr(k1,e3)
      Pair7 = dotrr(k2,e3)
      Abb1 = Pair12*Pair5
      Abb2 = Pair4*Pair5
      Abb7 = Pair15*Pair6
      Abb3 = Pair6*Pair7

      Abb4 = Pair6*Pair8
      Abb5 = Pair10*Pair9
      Abb6 = Pair11*Pair9
      AbbSum1 = -4*Abb1 + Abb2 + 4*Abb3 - 2*Abb4 - 2*Abb5 + 4*Abb6
      AbbSum2 = -2*Abb1 + Abb2 + Abb4 + 2*Abb6 - 2*Abb7 - 1/2.D0*Abb5

      Pair1 = dotrr(k1,k1)

      Ctree = - AbbSum2/
     -     ((-MW2 + Pair1 - 2*Pair13 + Pair14)) - 
     -     AbbSum1/(2.0*(-MW2 + Pair1 - 2*Pair2 + Pair3))

   
      return 
      end 
ccccccccccccccccccccccccccc
c     angle between to 4-vectors
c
c
      double precision function angle(p1,p2)
      implicit none
c     input
      double precision p1(0:4),p2(0:4)
      double precision p1sq,p2sq,p1p2
      double precision costheta12
c
      p1sq = p1(1)**2 + p1(2)**2 + p1(3)**2
      p2sq = p2(1)**2 + p2(2)**2 + p2(3)**2

      p1p2 = p1(1)*p2(1) + p1(2)*p2(2) + p1(3)*p2(3)

      costheta12 = p1p2/sqrt(p1sq)/sqrt(p2sq)


c      print*,'costheta12',costheta12

      angle = costheta12

      end
