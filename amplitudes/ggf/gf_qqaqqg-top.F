***********************************************************************
      subroutine AQQ_top (PBAR,SIGN,KBAR,GSIGN,phbar,hsign,
     &     ipolext,isig1,isig3,ampAtop3j21,ampAtop3j41) 
*----------------------------------------------------------------------
c... author   : Michael Kubocz
c... date     : 20.12.2008
c... modified : 10.04.2013 (MK)
*----------------------------------------------------------------------
c  modified: Francisco Campanario 
c  add:  Ward identities for pentagons and hexagons
c       + control of instabilities in Dble and QUAD

      IMPLICIT NONE

#include "ggf_process.inc"
c... for checks with MadGraph
#include "ggfloMG.inc" 

      logical SingularC(2)
      common/singluonver/SingularC
      
      character*17 wr1,wr2
      character*50 wr3,wr4
      PARAMETER(wr4='(A10,I1,A1,I2,A2,E28.20,E28.20)',
     &          wr1='(E28.20,E28.20)',
     &          wr2='(A10,E28.20)',wr3='(A5,I1,A2,E28.20)')

      real*8 ALFAS,XMT,ALFA,XMZ,XMW,SIN2W,XMH
      real*8 clr, xm2, xmg, b, v, a
      COMMON /BKOPIN/ ALFAS,XMT,ALFA,XMZ,XMW,SIN2W,XMH
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)

      real*8 col0,col1,col2,pi
      PARAMETER (col0=2.6666666666666667d0, col1=-0.3333333333333333d0,
     &           col2=0.1111111111111111d0)
      PARAMETER (pi=3.14159 26535 89793238d0)

      real*8 p(0:3,4),pbar(0:3,4),kbar(0:4),K(0:4),phbar(0:4),PH(0:3)


      real*8 p21(0:4),p43(0:4),p41(0:4),p23(0:4),
     &       p21K(0:4),p43K(0:4),p41K(0:4),p23K(0:4),
     &       p21KH(0:4),p43KH(0:4),p41KH(0:4),p23KH(0:4),
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!FC ATTENTION THIS ARE 0:3
     &       p21p43(0:3),p41p23(0:3),
!FC END ATTENTION
     &       K_H(0:4),DUMMY(0:4),DUMMY2(0:4)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cFC
c new arrays from 0:3, append a "p"
      real*8 Kp(0:3)

      real*8 p21p(0:3),p43p(0:3),p41p(0:3),p23p(0:3),
     &       p21Kp(0:3),p43Kp(0:3),p41Kp(0:3),p23Kp(0:3),
     &       p21KHp(0:3),p43KHp(0:3),p41KHp(0:3),p23KHp(0:3),
     &       p21p43p(0:3),p41p23p(0:3),
     &       K_Hp(0:3),DUMMYp(0:3),DUMMY2p(0:3)
cFC
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      complex*16 epsgC(0:3,2),GGG(0:3,-1:1,2,2)
      real*8 epsg(0:3,2)
      complex*16 epsAtop(0:3,2)

      complex*16 j21(0:3,-1:1),j21f(0:3,-1:1,2),j21i(0:3,-1:1,2),
     &           j43(0:3,-1:1),j43f(0:3,-1:1,2),j43i(0:3,-1:1,2),
     &           j41(0:3,-1:1),j41f(0:3,-1:1,2),j41i(0:3,-1:1,2),
     &           j23(0:3,-1:1),j23f(0:3,-1:1,2),j23i(0:3,-1:1,2),
     &           psi(2,-1:1,4),braketg(2,-1:1,8,2),
     &           j21Atopf(0:3,-1:1,2),j21Atopi(0:3,-1:1,2),
     &           j43Atopf(0:3,-1:1,2),j43Atopi(0:3,-1:1,2),
     &           j41Atopf(0:3,-1:1,2),j41Atopi(0:3,-1:1,2),
     &           j23Atopf(0:3,-1:1,2),j23Atopi(0:3,-1:1,2)

      integer ipol,ipolext,gsign,hsign,sign(4)

      complex*16 ampAtop3j21(12),ampAtop3j41(12)
 
      complex*16 amp21i,amp21f,amp43i,amp43f,
     &           amp41i,amp41f,amp23i,amp23f

c... Other local variables 
      integer i,j,mu,isig,isig1,isig3,icall
      integer ipolmin,isig1min,isig3min,ipolmax,isig1max,isig3max
      
      data icall /0/

      real*8 fac, vev
      save icall,fac

c... Variables for the main part of the program
      complex*16 aux234,aux34,aux5678(4)
      real*8 res12,res14,resI(3),rh(0:2)

      real*8 aux
c... Vertex functions
      complex*16 A2G_Loop,A3G_Loop

c... Other external functions
      complex*16 dotcc,dotrc,GGG_vertex
      real*8 dotrr

      EXTERNAL dotrr,dotcc,dotrc,GGG_vertex,
     &     A2G_Loop,A3G_Loop

cfc      verdbleSin=.false.
cfc      verquadSin=.false.
cfc      boxdbleSin=.false.
cfc      boxquadSin=.false.
cfc      boxSingular(1)=.false.
cfc      boxSingular(2)=.false.
cfc      verSingular(1)=.false.
cfc      verSingular(2)=.false.
cfc      SingularC(1)=.false.
cfc      SingularC(2)=.false.
      
      if (icall.eq.0) then
          vev = 2d0*xmw/(v(1,3)*dsqrt(8d0))
c... common overall factor
          fac=64d0*pi*alfas**5/(vev*vev)
          icall = 1
      endif

c... Define the internal momenta
      do mu = 0,3
         do i = 1,4
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
         K(mu)   = gsign*kbar(mu)
         PH(mu)  = hsign*phbar(mu)
c new arrays from 0:3, append a "p"
         Kp(mu) = K(mu)
       enddo

c... assign momenta combinations
      do mu = 0,3
         p21(mu) = p(mu,2)-p(mu,1) ! virtual gluon -q1
         p43(mu) = p(mu,4)-p(mu,3) ! virtual gluon -q2
         p41(mu) = p(mu,4)-p(mu,1) ! virtual gluon -q1 for res14
         p23(mu) = p(mu,2)-p(mu,3) ! virtual gluon -q2 for res14
         
         p21K(mu) = p21(mu)+K(mu)
         p43K(mu) = p43(mu)+K(mu)
         p41K(mu) = p41(mu)+K(mu)
         p23K(mu) = p23(mu)+K(mu)
         
         p21KH(mu) = p21K(mu)+PH(mu)
         p43KH(mu) = p43K(mu)+PH(mu)
         p41KH(mu) = p41K(mu)+PH(mu)
         p23KH(mu) = p23K(mu)+PH(mu)
         
         p21p43(mu)=-(p21(mu)+p43(mu))
         p41p23(mu)=-(p41(mu)+p23(mu))
c new arrays from 0:3, append a "p"
         p21p(mu) =          p21(mu)  
         p43p(mu) =          p43(mu)  
         p41p(mu) =          p41(mu)  
         p23p(mu) =          p23(mu)  
                                     
         p21Kp(mu) =         p21K(mu) 
         p43Kp(mu) =         p43K(mu) 
         p41Kp(mu) =         p41K(mu) 
         p23Kp(mu) =         p23K(mu) 
                                       
         p21KHp(mu)  =        p21KH(mu)
         p43KHp(mu)  =        p43KH(mu)
         p41KHp(mu)  =        p41KH(mu)
         p23KHp(mu)  =        p23KH(mu)
                                      
         p21p43p(mu) =        p21p43(mu)
         p41p23p(mu) =        p41p23(mu)
      enddo
      
c... assign m^2 to particles
      p21(4) = dotrr(p21p,p21p)
      p43(4) = dotrr(p43p,p43p)
      p41(4) = dotrr(p41p,p41p)
      p23(4) = dotrr(p23p,p23p)

      K(4)=0.0d0

      p21K(4) = dotrr(p21Kp,p21Kp) 
      p43K(4) = dotrr(p43Kp,p43Kp) 
      p41K(4) = dotrr(p41Kp,p41Kp) 
      p23K(4) = dotrr(p23Kp,p23Kp) 
      
      p21KH(4)=dotrr(p21KHp,p21KHp) 
      p43KH(4)=dotrr(p43KHp,p43KHp) 
      p41KH(4)=dotrr(p41KHp,p41KHp) 
      p23KH(4)=dotrr(p23KHp,p23KHp) 

c... get the external spinors (including factor sqrt(2E) )
      CALL psi0m(4,pbar(0,1),sign(1),psi)

c... get the f-fbar currents J21^mu, J43^mu
      CALL curr(1,psi(1,-1,2),psi(1,-1,1),j21(0,-1))
      CALL curr(1,psi(1,-1,4),psi(1,-1,3),j43(0,-1))
      CALL curr(1,psi(1,-1,4),psi(1,-1,1),j41(0,-1))
      CALL curr(1,psi(1,-1,2),psi(1,-1,3),j23(0,-1))

c... get the gluon polarization vector and the gluon emission spinors
      do ipol = 1,2	! 2 gluon polarization
         CALL polvec(kbar,ipol,epsg(0,ipol)) ! get gluon pol.vectors

      enddo
      
      if(MGdebug) then          ! MadGraph CHECK !!!!!!!!!!!!
         do mu=0,3              
            epsgC(mu,1)=1d0/sqrt(2d0)*( epsg(mu,1)+(0d0,1d0)*epsg(mu,2))
            epsgC(mu,2)=1d0/sqrt(2d0)*(-epsg(mu,1)+(0d0,1d0)*epsg(mu,2))
         enddo
      else
         do ipol = 1,2          ! 2 gluon polarizations
            do mu=0,3
               epsgC(mu,ipol)= epsg(mu,ipol) ! complex pol. vector
            enddo
         enddo
      endif
      
      
c...  for gauge check:
      if (lgauge2) then   ! contract amplitude with q rather than eps(q)
         do ipol = 1,2   ! 2 gluon polarizations
            do mu = 0,3
               epsgC(mu,ipol)= K(mu)/abs(K(0))
            enddo	
         enddo
      endif

      do ipol = 1,2
         do isig = -1,1,2  ! quark helicities 
            CALL ket2c(psi(1,isig,1),.false.,p(0,1),isig,K,epsgC(0,ipol)
     &           ,braketg(1,isig,1,ipol),DUMMY)       ! |K,1>_ipol,isig1
            CALL bra2c(psi(1,isig,2),.false.,p(0,2),isig,K,epsgC(0,ipol)
     &           ,braketg(1,isig,2,ipol),DUMMY)       ! <2,K|_ipol,isig2
            CALL ket2c(psi(1,isig,3),.false.,p(0,3),isig,K,epsgC(0,ipol)
     &           ,braketg(1,isig,3,ipol),DUMMY)       ! |K,3>_ipol,isig3
            CALL bra2c(psi(1,isig,4),.false.,p(0,4),isig,K,epsgC(0,ipol)
     &           ,braketg(1,isig,4,ipol),DUMMY)       ! <4,K|_ipol,isig4
         enddo
      enddo

c... gluon from upper line:
      do ipol = 1,2 ! gluon polarizations
         CALL curr(1,psi(1,-1,2),braketg(1,-1,1,ipol),j21i(0,-1,ipol))	
c...                                        =  <2|vertex|q,1>_ipol,isig1
         CALL curr(1,braketg(1,-1,2,ipol),psi(1,-1,1),j21f(0,-1,ipol))	
c...                                        =  <2,q|vertex|1>_ipol,isig1

         CALL curr(1,psi(1,-1,4),braketg(1,-1,1,ipol),j41i(0,-1,ipol))	
c...                                        =  <4|vertex|q,1>_ipol,isig1
         CALL curr(1,braketg(1,-1,4,ipol),psi(1,-1,1),j41f(0,-1,ipol))	
c...                                        =  <4,q|vertex|1>_ipol,isig1

c... gluon from lower line:
         CALL curr(1,psi(1,-1,4),braketg(1,-1,3,ipol),j43i(0,-1,ipol))	
c...                                        =  <4|gam.mu|q,3>_ipol,isig3
         CALL curr(1,braketg(1,-1,4,ipol),psi(1,-1,3),j43f(0,-1,ipol))	
c...                                        =  <4,q|gam.mu|3>_ipol,isig3

         CALL curr(1,psi(1,-1,2),braketg(1,-1,3,ipol),j23i(0,-1,ipol))	
c...                                        =  <2|gam.mu|q,3>_ipol,isig3
         CALL curr(1,braketg(1,-1,2,ipol),psi(1,-1,3),j23f(0,-1,ipol))	
c...                                        =  <2,q|gam.mu|3>_ipol,isig3
      enddo
      do ipol=1,2 ! gluon polarisations
         if(AloopTop) then
            CALL EpsA_Loop(Kp,epsgC(0,ipol),PH,K_Hp,
     &           epsAtop(0,ipol),m_top)
         endif
c... 4-dim momenta variables
         do mu = 0,3
            K_H(mu)=K_Hp(mu)
         enddo
         K_H(4)=dotrr(K_Hp,K_Hp)
         do isig=-1,1,2 ! quark helicities
            
         if(AloopTop) then
c... with effective polarization vector 
c... |A,K,1>_ipol,isig1
            CALL ket2c(psi(1,isig,1),.false.,p(0,1),isig,K_H,
     &           epsAtop(0,ipol),braketg(1,isig,5,ipol),DUMMY)
c... <2,A,K|_ipol,isig2
            CALL bra2c(psi(1,isig,2),.false.,p(0,2),isig,K_H,
     &           epsAtop(0,ipol),braketg(1,isig,6,ipol),DUMMY)
c... |A,K,3>_ipol,isig3
            CALL ket2c(psi(1,isig,3),.false.,p(0,3),isig,K_H,
     &           epsAtop(0,ipol),braketg(1,isig,7,ipol),DUMMY)
c... <4,A,K|_ipol,isig4
            CALL bra2c(psi(1,isig,4),.false.,p(0,4),isig,K_H,
     &           epsAtop(0,ipol),braketg(1,isig,8,ipol),DUMMY)

c... gluon from upper line:
c... <2|vertex|A,q,1>_ipol,isig1
            CALL curr(1,psi(1,-1,2),braketg(1,-1,5,ipol),
     &           j21Atopi(0,-1,ipol))	
c... <2,A,q|vertex|1>_ipol,isig1
            CALL curr(1,braketg(1,-1,6,ipol),psi(1,-1,1),
     &           j21Atopf(0,-1,ipol))	
c... <4|vertex|A,q,1>_ipol,isig1
            CALL curr(1,psi(1,-1,4),braketg(1,-1,5,ipol),
     &           j41Atopi(0,-1,ipol))	
c... <4,A,q|vertex|1>_ipol,isig1
            CALL curr(1,braketg(1,-1,8,ipol),psi(1,-1,1),
     &           j41Atopf(0,-1,ipol))	

c... gluon from lower line:
c... <4|gam.mu|A,q,3>_ipol,isig3
            CALL curr(1,psi(1,-1,4),braketg(1,-1,7,ipol),
     &           j43Atopi(0,-1,ipol))	
c... <4,A,q|gam.mu|3>_ipol,isig3
            CALL curr(1,braketg(1,-1,8,ipol),psi(1,-1,3),
     &           j43Atopf(0,-1,ipol))	
c... <2|gam.mu|A,q,3>_ipol,isig3

            CALL curr(1,psi(1,-1,2),braketg(1,-1,7,ipol),
     &           j23Atopi(0,-1,ipol))	
c... <2,A,q|gam.mu|3>_ipol,isig3
            CALL curr(1,braketg(1,-1,6,ipol),psi(1,-1,3),
     &           j23Atopf(0,-1,ipol))	
         endif
      enddo
      enddo

      do i=1,12
         ampAtop3j21(i)= DCMPLX(0d0,0d0)       
         ampAtop3j41(i)= DCMPLX(0d0,0d0)     
      enddo
      amp21i = DCMPLX(0d0,0d0)
      amp21f = DCMPLX(0d0,0d0)
      amp43i = DCMPLX(0d0,0d0)
      amp43f = DCMPLX(0d0,0d0)
      amp41i = DCMPLX(0d0,0d0)
      amp41f = DCMPLX(0d0,0d0)
      amp23i = DCMPLX(0d0,0d0)
      amp23f = DCMPLX(0d0,0d0)
      aux234 = DCMPLX(0d0,0d0)
      aux34  = DCMPLX(0d0,0d0)
      do i=1,4 
         aux5678(i) = DCMPLX(0d0,0d0)
      enddo
      res12 = 0.0d0
      res14 = 0.0d0
      do i=1,3
         resI(i) = 0.0d0
      enddo
      do i=0,2
         rh(i)=0.0d0
      enddo

c$$$      if(MGdebug) then
c$$$        ipolmin = 1
c$$$        isig1min= 1
c$$$        isig3min= 1
c$$$        ipolmax = 1
c$$$        isig1max= 1
c$$$        isig3max= 1
c$$$      else
c$$$        ipolmin = 1
c$$$        isig1min= -1
c$$$        isig3min= -1
c$$$        ipolmax = 2
c$$$        isig1max= 1
c$$$        isig3max= 1
c$$$      endif


c     if (MGdebug.and.ipol.ne.1.and.isig1.ne.1
c     &                .and.isig3.ne.1) cycle
      if(AloopTop) then
         ampAtop3j21(1)=cplModAtop
     &        *dotcc(J21Atopi(0,isig1,ipolext),J43(0,isig3))/p21KH(4)
         ampAtop3j21(2)=cplModAtop
     &        *A2G_Loop(p21K,p43,J21i(0,isig1,ipolext),J43(0,isig3),
     &        m_top)/(p21K(4)*p43(4))
         ampAtop3j21(3)=cplModAtop
     &        *dotcc(J21(0,isig1),J43Atopi(0,isig3,ipolext))/p43KH(4)
         ampAtop3j21(4)=cplModAtop
     &        *A2G_Loop(p21,p43K,J21(0,isig1),J43i(0,isig3,ipolext),
     &        m_top)/(p21(4)*p43K(4))
         ampAtop3j21(5)=cplModAtop
     &        *dotcc(J21Atopf(0,isig1,ipolext),J43(0,isig3))/p21KH(4)
         ampAtop3j21(6)=cplModAtop
     &        *A2G_Loop(p21K,p43,J21f(0,isig1,ipolext),J43(0,isig3),
     &        m_top)/(p21K(4)*p43(4))
         ampAtop3j21(7)=cplModAtop
     &        *dotcc(J21(0,isig1),J43Atopf(0,isig3,ipolext))/p43KH(4)
         ampAtop3j21(8)=cplModAtop
     &        *A2G_Loop(p21,p43K,J21(0,isig1),J43f(0,isig3,ipolext),
     &        m_top)/(p21(4)*p43K(4))
         ampAtop3j21(9)=-cplModAtop
     &        *GGG_vertex(p21p,p43p,p21p43,J21(0,isig1),
     &        J43(0,isig3),epsAtop(0,ipolext))/(p21(4)*p43(4))
         ampAtop3j21(10)=cplModAtop
     &        *A3G_Loop(K,p21,p43,epsgC(0,ipolext),J21(0,isig1),
     &        J43(0,isig3),m_top)/(p21(4)*p43(4))
         CALL GGG_current(p43p,J43(0,isig3),Kp,
     &        epsgC(0,ipolext),DUMMY2p,GGG(0,isig3,ipolext,1))
         CALL GGG_current(p21p,J21(0,isig1),Kp,
     &        epsgC(0,ipolext),DUMMY2p,GGG(0,isig1,ipolext,2))
         ampAtop3j21(11)=cplModAtop
     &        *A2G_Loop(p21,p43K,J21(0,isig1),GGG(0,isig3,ipolext,1),
     &        m_top)/(p21(4)*p43(4))
         ampAtop3j21(12)=cplModAtop
     &        *A2G_Loop(p21K,p43,GGG(0,isig1,ipolext,2),J43(0,isig3),
     &        m_top)/(p21(4)*p43(4))
c...  interchanged final state quarks
         ampAtop3j41(1)=cplModAtop
     &        *dotcc(J41Atopi(0,isig1,ipolext),J23(0,isig3))/p41KH(4)
         ampAtop3j41(2)=cplModAtop
     &        *A2G_Loop(p41K,p23,J41i(0,isig1,ipolext),J23(0,isig3),
     &        m_top)/(p41K(4)*p23(4))
         ampAtop3j41(3)=cplModAtop
     &        *dotcc(J41(0,isig1),J23Atopi(0,isig3,ipolext))/p23KH(4)
         ampAtop3j41(4)=cplModAtop
     &        *A2G_Loop(p41,p23K,J41(0,isig1),J23i(0,isig3,ipolext),
     &        m_top)/(p41(4)*p23K(4))
         ampAtop3j41(5)=cplModAtop
     &        *dotcc(J41(0,isig1),J23Atopf(0,isig3,ipolext))/p23KH(4)
         ampAtop3j41(6)=cplModAtop
     &        *A2G_Loop(p41,p23K,J41(0,isig1),J23f(0,isig3,ipolext),
     &        m_top)/(p41(4)*p23K(4))
         ampAtop3j41(7)=cplModAtop
     &        *dotcc(J41Atopf(0,isig1,ipolext),J23(0,isig3))/p41KH(4)
         ampAtop3j41(8)=cplModAtop
     &        *A2G_Loop(p41K,p23,J41f(0,isig1,ipolext),J23(0,isig3),
     &        m_top)/(p41K(4)*p23(4))
         ampAtop3j41(9)=cplModAtop
     &        *GGG_vertex(p41p,p23p,p41p23,J41(0,isig1),
     &        J23(0,isig3),epsAtop(0,ipolext))/(p41(4)*p23(4))
         ampAtop3j41(10)=-cplModAtop 
     &        *A3G_Loop(K,p41,p23,epsgC(0,ipolext),J41(0,isig1),
     &        J23(0,isig3),m_top)/(p41(4)*p23(4))
         CALL GGG_current(p23p,J23(0,isig3),Kp,
     &        epsgC(0,ipolext),DUMMY2p,GGG(0,isig3,ipolext,1))
         CALL GGG_current(p41p,J41(0,isig1),Kp,
     &        epsgC(0,ipolext),DUMMY2p,GGG(0,isig1,ipolext,2))
         ampAtop3j41(11)=cplModAtop
     &        *A2G_Loop(p41K,p23,GGG(0,isig1,ipolext,2),J23(0,isig3),
     &        m_top)/(p41(4)*p23(4))
         ampAtop3j41(12)=cplModAtop
     &        *A2G_Loop(p41,p23K,J41(0,isig1),GGG(0,isig3,ipolext,1),
     &        m_top)/(p41(4)*p23(4))
      endif


c... store amplitudes for check with MadGraph
      if(MGdebug) then
         do i=1,12
            if(AloopTop) then
               GGFLOamp1(isig1,isig3,ipolext,i)
     &              =ampAtop3j21(i)*dsqrt(fac)
               GGFLOamp1(isig1,isig3,ipolext,i+12)
     &              =ampAtop3j41(i)*dsqrt(fac)
            endif
         enddo
c     print*,'ipolext,isig1,isig3',ipolext,isig1,isig3
c     do i=1,24
c     print*,'GGFLOamp1(',i,')',GGFLOamp1(isig1,
c     &                    isig3,ipolext,i)
c     enddo
c     print*,''
      endif
               
c... gauge check
      if (lgauge2) then
         aux=0.0d0
         if(AloopTop) then
            do i=1,12
               aux=aux+abs(ampAtop3j21(i))
            enddo
            do i=1,12
               aux=aux+abs(ampAtop3j41(i))
            enddo
         endif
         print*,'ipolext',ipolext
         print*,'aux',aux
         print*,'norm. amp21i',amp21i/aux
         print*,'norm. amp21f',amp21f/aux
         print*,'norm. amp43i',amp43i/aux
         print*,'norm. amp43f',amp43f/aux
      endif
      
cfcccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
cfcc 
cfcc         Check whether there is instabilities
cfcc
cfcc!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cfcc!                                                       !! 
cfcc!         IF SO ::::-->   Rh(i) -> 0                     !!
cfcc!                                                       !!
cfcc!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
cfcc
cfcc         Count instabilites for hex,pen,box,vertices
cfcc         both in quad and dble precision
cfcc           
cfccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cfcc 
cfcc   Check Vertices
cfcc     
cfcc        if(boxeva) then
cfcc   Dble precision
cfc           If(verdbleSin) then
cfcc          print*,"verdbleSi",verdbleSin
cfcc          stop
cfcc bad dble
cfc        verinst(gf_iteration,1,1)=verinst(gf_iteration,1,1)+1
cfcc        print*, " verinst(gf_iteration,1,1)", verinst(gf_iteration,1,1)
cfc       else
cfcc good dble 
cfc        verinst(gf_iteration,2,1)=verinst(gf_iteration,2,1)+1
cfcc        print*, " verinst(gf_iteration,2,1)", verinst(gf_iteration,2,1)
cfc        endif
cfc
cfc       If(verQUADSin) then
cfcc bad QUAD
cfc        verinst(gf_iteration,1,2)=verinst(gf_iteration,1,2)+1
cfccccccccccccccccccccccccccccccccccccccccccccccccccccc
cfcc
cfc        rh(1)=0d0
cfc        rh(2)=0d0
cfc        rh(0)=0d0
cfcc
cfccccccccccccccccccccccccccccccccccccccccccccccccccccc
cfc       else
cfcc good QUAD 
cfc        verinst(gf_iteration,2,2)=verinst(gf_iteration,2,2)+1
cfc        endif
cfcc        endif
cfcc 
cfcc   Check Boxes
cfcc     
cfc        if(boxeva) then
cfcc   Dble precision
cfc       If(boxdbleSin) then
cfcc bad dble
cfc        boxinst(gf_iteration,1,1)=boxinst(gf_iteration,1,1)+1
cfc       else
cfcc good dble 
cfc        boxinst(gf_iteration,2,1)=boxinst(gf_iteration,2,1)+1
cfc        endif
cfc
cfc       If(boxQUADSin) then
cfcc bad QUAD
cfc        boxinst(gf_iteration,1,2)=boxinst(gf_iteration,1,2)+1
cfccccccccccccccccccccccccccccccccccccccccccccccccccccc
cfcc
cfc        rh(1)=0d0
cfc        rh(2)=0d0
cfc        rh(0)=0d0        
cfcc
cfccccccccccccccccccccccccccccccccccccccccccccccccccccc
cfc       else
cfcc good QUAD 
cfc        boxinst(gf_iteration,2,2)=boxinst(gf_iteration,2,2)+1
cfc        endif
cfc        endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c... switch between (non-)identical quarks for check with MadGraph
      if(MGdebug) then
         if(mguu) then
            GGFLOampSq1= rh(0)  ! for u u
         else
            GGFLOampSq1= rh(1)  ! for u d 
         endif

      endif

c... gauge check
      if (lgauge2) then
         print*,'rh(1)',rh(1)/aux*aux
         print*,'rh(2)',rh(2)/aux*aux
         print*,'rh(0)',rh(0)/aux*aux
         STOP
      endif


      END
