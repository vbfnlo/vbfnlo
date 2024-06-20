***********************************************************************
      subroutine Higgs3juuxccxgH(PBAR,SIGN,KBAR,GSIGN,phbar,hsign,rh) 

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
      complex*16 epsHeff(0:3,2), epsHtop(0:3,2)

      complex*16 j21(0:3,-1:1),j21f(0:3,-1:1,2),j21i(0:3,-1:1,2),
     &           j43(0:3,-1:1),j43f(0:3,-1:1,2),j43i(0:3,-1:1,2),
     &           j41(0:3,-1:1),j41f(0:3,-1:1,2),j41i(0:3,-1:1,2),
     &           j23(0:3,-1:1),j23f(0:3,-1:1,2),j23i(0:3,-1:1,2),
     &           psi(2,-1:1,4),braketg(2,-1:1,8,2),
     &           j21Hefff(0:3,-1:1,2),j21Heffi(0:3,-1:1,2),
     &           j43Hefff(0:3,-1:1,2),j43Heffi(0:3,-1:1,2),
     &           j41Hefff(0:3,-1:1,2),j41Heffi(0:3,-1:1,2),
     &           j23Hefff(0:3,-1:1,2),j23Heffi(0:3,-1:1,2),
     &           j21Htopf(0:3,-1:1,2),j21Htopi(0:3,-1:1,2),
     &           j43Htopf(0:3,-1:1,2),j43Htopi(0:3,-1:1,2),
     &           j41Htopf(0:3,-1:1,2),j41Htopi(0:3,-1:1,2),
     &           j23Htopf(0:3,-1:1,2),j23Htopi(0:3,-1:1,2)

      integer ipol,gsign,hsign,sign(4)

      complex*16 ampHeff3j21(12),ampHeff3j41(12),
     &  ampAeff3j21(12),ampAeff3j41(12),
     &  ampHtop3j21(12),ampHtop3j41(12),ampHbottom3j21(12),
     &  ampHbottom3j41(12),ampAtop3j21(12),ampAtop3j41(12),
     &  ampAbottom3j21(12),ampAbottom3j41(12)   
 
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
      complex*16 A2G_Loop,A3G_Loop,H2G_Loop,H3G_Loop,
     &           H2G_effective,H3G_effective,
     &           A2G_effective,A3G_effective
      complex*16 H1G_SqLoop,H2G_SqLoop

c... MG specific variables for color structures and amplitudes
      real*8 RGB(4,4)
      complex*16 color(4),colorT,ampQQ(12)

c... Other external functions
      complex*16 dotcc,dotrc,GGG_vertex
      real*8 dotrr

      EXTERNAL dotrr,dotcc,dotrc,GGG_vertex,
     &     H2G_effective,H3G_effective,A2G_effective,A3G_effective,
     &     H2G_Loop,H3G_Loop,A2G_Loop,A3G_Loop


      GGFLOampSq1=0d0

      verdbleSin=.false.
      verquadSin=.false.
      boxdbleSin=.false.
      boxquadSin=.false.
      boxSingular(1)=.false.
      boxSingular(2)=.false.
      verSingular(1)=.false.
      verSingular(2)=.false.
      
      SingularC(1)=.false.
      SingularC(2)=.false.
      
      if (icall.eq.0) then
          vev = 2d0*xmw/(v(1,3)*dsqrt(8d0))
c... common overall factor
          fac=64d0*pi*alfas**5/(vev*vev)
          icall = 1
      endif

      RGB(1,2)=-1d0/3d0  !-2d0/3d0  2 --> 3
      RGB(1,3)=7d0/3d0   !-1d0/3d0  3 --> 4
      RGB(1,4)=-2d0/3d0  ! 7d0/3d0  4 --> 2
      RGB(2,3)=-2d0/3d0  ! 7d0/3d0
      RGB(2,4)=7d0/3d0   !-1d0/3d0
      RGB(3,4)=-1d0/3d0  !-2d0/3d0
      do i=1,4
         RGB(i,i)=8d0/3d0
         do j=1,4
            if(i.gt.j) then
               RGB(i,j)=RGB(j,i)
            endif
         enddo
      enddo
cmk      RGB(1,2)=0d0
cmk      RGB(1,3)=0d0
cmk      RGB(1,4)=0d0
cmk      RGB(2,3)=0d0
cmk      RGB(2,4)=0d0
cmk      RGB(3,4)=0d0
cmk      do i=1,4
cmk         RGB(i,i)=0d0
cmk         do j=1,4
cmk            if(i.gt.j) then
cmk               RGB(i,j)=RGB(j,i)
cmk            endif
cmk         enddo
cmk      enddo

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

cmk         Print*, "psi1 hel 1",  psi(:,1,1)
cmk         Print*, "psi1 hel -1", psi(:,-1,1)
cmk         print*
cmk         Print*, "psi2 hel 1",  psi(:,1,2)
cmk         Print*, "psi2 hel -1", psi(:,-1,2)
cmk         print*
cmk         Print*, "psi3 hel 1",  psi(:,1,3)
cmk         Print*, "psi3 hel -1", psi(:,-1,4)
cmk         print*
cmk         Print*, "psi4 hel 1",  psi(:,1,4)
cmk         Print*, "psi4 hel -1", psi(:,-1,4)
cmk         print*

cfc         Print*, "J21 current hel 1",  j21(0:4,1)
cfc         Print*, "J21 current hel -1", j21(0:4,-1)
cfc
cfc         Print*, "J43 current hel 1",  j43(0:4,1)
cfc         Print*, "J43 current hel -1", j43(0:4,-1)
cfc
cfc         Print*, "J41 current hel 1",  j41(0:4,1)
cfc         Print*, "J41 current hel -1", j41(0:4,-1)
cfc
cfc         Print*, "J23 current hel 1",  j23(0:4,1)
cfc         Print*, "J23 current hel -1", j23(0:4,-1)

         do mu=0,3              
            epsgC(mu,1)=1d0/sqrt(2d0)*( epsg(mu,1)+(0d0,1d0)*epsg(mu,2))
            epsgC(mu,2)=1d0/sqrt(2d0)*(-epsg(mu,1)+(0d0,1d0)*epsg(mu,2))
         enddo
cmk         Print*,"MG_check pol 1 ", epsgC(0:3,1)
cmk         print*
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
         if(Heff) then
            CALL EpsH_effective(Kp,.false.,epsgC(0,ipol),PH,K_Hp
     &           ,epsHeff(0,ipol))
          endif
         if(HloopTop) then
            CALL EpsH_Loop(K,.false.,epsgC(0,ipol),PH,K_Hp,
     &           epsHtop(0,ipol),m_top)
         endif

c... 4-dim momenta variables
         do mu = 0,3
            K_H(mu)=K_Hp(mu)
         enddo
         K_H(4)=dotrr(K_Hp,K_Hp)
         do isig=-1,1,2 ! quark helicities
         if(Heff) then
c... with effective polarization vector 
c... |H,K,1>_ipol,isig1
            CALL ket2c(psi(1,isig,1),.false.,p(0,1),isig,K_H,
     &           epsHeff(0,ipol),braketg(1,isig,5,ipol),DUMMY) 
c... <2,H,K|_ipol,isig2
            CALL bra2c(psi(1,isig,2),.false.,p(0,2),isig,K_H,
     &           epsHeff(0,ipol),braketg(1,isig,6,ipol),DUMMY) 
c... |H,K,3>_ipol,isig3
            CALL ket2c(psi(1,isig,3),.false.,p(0,3),isig,K_H,
     &           epsHeff(0,ipol),braketg(1,isig,7,ipol),DUMMY) 
c... <4,H,K|_ipol,isig4
            CALL bra2c(psi(1,isig,4),.false.,p(0,4),isig,K_H,
     &           epsHeff(0,ipol),braketg(1,isig,8,ipol),DUMMY) 

c... gluon from upper line:
c... <2|vertex|H,q,1>_ipol,isig1
            CALL curr(1,psi(1,-1,2),braketg(1,-1,5,ipol),
     &           j21Heffi(0,-1,ipol))	
c... <2,H,q|vertex|1>_ipol,isig1
            CALL curr(1,braketg(1,-1,6,ipol),psi(1,-1,1),
     &           j21Hefff(0,-1,ipol))	
c... <4|vertex|H,q,1>_ipol,isig1
            CALL curr(1,psi(1,-1,4),braketg(1,-1,5,ipol),
     &           j41Heffi(0,-1,ipol))	
c... <4,H,q|vertex|1>_ipol,isig1
            CALL curr(1,braketg(1,-1,8,ipol),psi(1,-1,1),
     &           j41Hefff(0,-1,ipol))	

c... gluon from lower line:
c... <4|gam.mu|H,q,3>_ipol,isig3
            CALL curr(1,psi(1,-1,4),braketg(1,-1,7,ipol),
     &           j43Heffi(0,-1,ipol))	
c... <4,H,q|gam.mu|3>_ipol,isig3
            CALL curr(1,braketg(1,-1,8,ipol),psi(1,-1,3),
     &           j43Hefff(0,-1,ipol))	
c... <2|gam.mu|H,q,3>_ipol,isig3
            CALL curr(1,psi(1,-1,2),braketg(1,-1,7,ipol),
     &           j23Heffi(0,-1,ipol))	
c...  <2,H,q|gam.mu|3>_ipol,isig3
            CALL curr(1,braketg(1,-1,6,ipol),psi(1,-1,3),
     &           j23Hefff(0,-1,ipol))	
         endif
         if(HloopTop) then
c... with effective polarization vector 
c... |H,K,1>_ipol,isig1
            CALL ket2c(psi(1,isig,1),.false.,p(0,1),isig,K_H,
     &           epsHtop(0,ipol),braketg(1,isig,5,ipol),DUMMY) 
c... <2,H,K|_ipol,isig2
            CALL bra2c(psi(1,isig,2),.false.,p(0,2),isig,K_H,
     &           epsHtop(0,ipol),braketg(1,isig,6,ipol),DUMMY) 
c... |H,K,3>_ipol,isig3
            CALL ket2c(psi(1,isig,3),.false.,p(0,3),isig,K_H,
     &           epsHtop(0,ipol),braketg(1,isig,7,ipol),DUMMY)
c... <4,H,K|_ipol,isig4 
            CALL bra2c(psi(1,isig,4),.false.,p(0,4),isig,K_H,
     &           epsHtop(0,ipol),braketg(1,isig,8,ipol),DUMMY)
c... gluon from upper line:
c... <2|vertex|H,q,1>_ipol,isig1
            CALL curr(1,psi(1,-1,2),braketg(1,-1,5,ipol),
     &           j21Htopi(0,-1,ipol))	
c... <2,H,q|vertex|1>_ipol,isig1
            CALL curr(1,braketg(1,-1,6,ipol),psi(1,-1,1),
     &           j21Htopf(0,-1,ipol))	
c... <4|vertex|H,q,1>_ipol,isig1
            CALL curr(1,psi(1,-1,4),braketg(1,-1,5,ipol),
     &           j41Htopi(0,-1,ipol))	
c... <4,H,q|vertex|1>_ipol,isig1
            CALL curr(1,braketg(1,-1,8,ipol),psi(1,-1,1),
     &           j41Htopf(0,-1,ipol))	

c... gluon from lower line:
c... <4|gam.mu|H,q,3>_ipol,isig3
            CALL curr(1,psi(1,-1,4),braketg(1,-1,7,ipol),
     &           j43Htopi(0,-1,ipol))	
c... <4,H,q|gam.mu|3>_ipol,isig3
            CALL curr(1,braketg(1,-1,8,ipol),psi(1,-1,3),
     &           j43Htopf(0,-1,ipol))	
c... <2|gam.mu|H,q,3>_ipol,isig3
            CALL curr(1,psi(1,-1,2),braketg(1,-1,7,ipol),
     &           j23Htopi(0,-1,ipol))	
c... <2,H,q|gam.mu|3>_ipol,isig3
            CALL curr(1,braketg(1,-1,6,ipol),psi(1,-1,3),
     &           j23Htopf(0,-1,ipol))	
         endif
      enddo
      enddo

      do i=1,12
         ampHeff3j21(i)= DCMPLX(0d0,0d0)       
         ampHtop3j21(i)= DCMPLX(0d0,0d0)       
         ampQQ(i)=DCMPLX(0d0,0d0)
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
      do i=1,4
         color(i) =DCMPLX(0d0,0d0)
      enddo
      colorT=DCMPLX(0d0,0d0)

      if(MGdebug) then
         
        ipolmin = 1
        isig1min= 1
        isig3min= 1
        ipolmax = 1
        isig1max= 1
        isig3max= 1
      else
        ipolmin = 1
        ipolmax = 2
        isig1min= -1
        isig3min= -1
        isig1max= 1
        isig3max= 1
      endif

**** loop over the quark and gluon polarizations ****
      do ipol=ipolmin,ipolmax
        do isig1 = isig1min,isig1max,2
           do isig3 = isig3min,isig3max,2
c                  if (MGdebug.and.ipol.ne.1.and.isig1.ne.1
c     &                .and.isig3.ne.1) cycle
               if(Heff) then
                ampHeff3j21(1)=cplModHeff
     &           *dotcc(J21Heffi(0,isig1,ipol),J43(0,isig3))/p21KH(4)
                ampHeff3j21(2)=cplModHeff
     &           *H2G_effective(p21Kp,p43p,J21i(0,isig1,ipol),
     &               J43(0,isig3))/(p21K(4)*p43(4))
                ampHeff3j21(3)=cplModHeff
     &           *dotcc(J21(0,isig1),J43Heffi(0,isig3,ipol))/p43KH(4)
                ampHeff3j21(4)=cplModHeff
     &           *H2G_effective(p21p,p43Kp,J21(0,isig1),
     &               J43i(0,isig3,ipol))/(p21(4)*p43K(4))
                ampHeff3j21(5)=cplModHeff
     &           *dotcc(J21Hefff(0,isig1,ipol),J43(0,isig3))/p21KH(4)
                ampHeff3j21(6)=cplModHeff
     &           *H2G_effective(p21Kp,p43p,J21f(0,isig1,ipol),
     &               J43(0,isig3))/(p21K(4)*p43(4))
                ampHeff3j21(7)=cplModHeff
     &           *dotcc(J21(0,isig1),J43Hefff(0,isig3,ipol))/p43KH(4)
                ampHeff3j21(8)=cplModHeff
     &           *H2G_effective(p21p,p43Kp,J21(0,isig1),
     &               J43f(0,isig3,ipol))/(p21(4)*p43K(4))
                ampHeff3j21(9)=-cplModHeff
     &           *GGG_vertex(p21p,p43p,p21p43,J21(0,isig1),
     &               J43(0,isig3),epsHeff(0,ipol))/(p21(4)*p43(4))
                ampHeff3j21(10)=cplModHeff
     &           *H3G_effective(Kp,p21p,p43p,epsgC(0,ipol),
     &               J21(0,isig1),J43(0,isig3))/(p21(4)*p43(4))
                CALL GGG_current(p43p,J43(0,isig3),Kp,
     &               epsgC(0,ipol),DUMMY2p,GGG(0,isig3,ipol,1))
                CALL GGG_current(p21p,J21(0,isig1),Kp,
     &               epsgC(0,ipol),DUMMY2p,GGG(0,isig1,ipol,2))
                ampHeff3j21(11)=cplModHeff
     &               *H2G_effective(p21p,p43Kp,J21(0,isig1),
     &               GGG(0,isig3,ipol,1))/(p21(4)*p43(4))
                ampHeff3j21(12)=cplModHeff  
     &            *H2G_effective(p21Kp,p43p,GGG(0,isig1,ipol,2),
     &            J43(0,isig3))/(p21(4)*p43(4))
               endif
*-----------------------------------------------------------------------
               if(HloopTop) then
                ampHtop3j21(1)=cplModHtop
     &           *dotcc(J21Htopi(0,isig1,ipol),J43(0,isig3))/p21KH(4)
                ampHtop3j21(2)=cplModHtop
     &           *H2G_Loop(p21K,p43,J21i(0,isig1,ipol),
     &               J43(0,isig3),m_top)/(p21K(4)*p43(4))
                ampHtop3j21(3)=cplModHtop
     &           *dotcc(J21(0,isig1),J43Htopi(0,isig3,ipol))/p43KH(4)
                ampHtop3j21(4)=cplModHtop
     &           *H2G_Loop(p21,p43K,J21(0,isig1),J43i(0,isig3,ipol),
     &               m_top)/(p21(4)*p43K(4))
                ampHtop3j21(5)=cplModHtop
     &           *dotcc(J21Htopf(0,isig1,ipol),J43(0,isig3))/p21KH(4)
                ampHtop3j21(6)=cplModHtop
     &           *H2G_Loop(p21K,p43,J21f(0,isig1,ipol),J43(0,isig3),
     &               m_top)/(p21K(4)*p43(4))
                ampHtop3j21(7)=cplModHtop
     &           *dotcc(J21(0,isig1),J43Htopf(0,isig3,ipol))/p43KH(4)
                ampHtop3j21(8)=cplModHtop
     &           *H2G_Loop(p21,p43K,J21(0,isig1),J43f(0,isig3,ipol),
     &               m_top)/(p21(4)*p43K(4))
                ampHtop3j21(9)=-cplModHtop
     &           *GGG_vertex(p21p,p43p,p21p43,J21(0,isig1),
     &               J43(0,isig3),epsHtop(0,ipol))/(p21(4)*p43(4))
                ampHtop3j21(10)=cplModHtop
     &           *H3G_Loop(K,p21,p43,epsgC(0,ipol),J21(0,isig1),
     &                     J43(0,isig3),m_top)/(p21(4)*p43(4))
                CALL GGG_current(p43p,J43(0,isig3),Kp,
     &               epsgC(0,ipol),DUMMY2p,GGG(0,isig3,ipol,1))
                CALL GGG_current(p21p,J21(0,isig1),Kp,
     &               epsgC(0,ipol),DUMMY2p,GGG(0,isig1,ipol,2))
                ampHtop3j21(11)=cplModHtop
     &           *H2G_Loop(p21,p43K,J21(0,isig1),GGG(0,isig3,ipol,1),
     &               m_top)/(p21(4)*p43(4))
                ampHtop3j21(12)=cplModHtop
     &           *H2G_Loop(p21K,p43,GGG(0,isig1,ipol,2),J43(0,isig3),
     &               m_top)/(p21(4)*p43(4))
               endif
               
               do i=1,12
                  ampQQ(i)=ampHeff3j21(i)+ampHtop3j21(i)
               enddo
               ampQQ(3)=ampHeff3j21(5)+ampHtop3j21(5)
               ampQQ(4)=ampHeff3j21(6)+ampHtop3j21(6)
               ampQQ(5)=ampHeff3j21(7)+ampHtop3j21(7)
               ampQQ(6)=ampHeff3j21(8)+ampHtop3j21(8)
               ampQQ(7)=ampHeff3j21(3)+ampHtop3j21(3)
               ampQQ(8)=ampHeff3j21(4)+ampHtop3j21(4)

               color(1)=-ampQQ(1)-ampQQ(2)+ampQQ(9)-ampQQ(10)-ampQQ(12)
cmk            color(1)=-ampQQ(1)-ampQQ(2)+ampQQ(9)-ampQQ(10)-ampQQ(12)
               color(2)=-ampQQ(3)-ampQQ(4)-ampQQ(9)+ampQQ(10)+ampQQ(12)
cmk            color(2)=-ampQQ(3)-ampQQ(4)-ampQQ(11)
               color(3)=-ampQQ(5)-ampQQ(6)+ampQQ(11)
cmk            color(3)=-ampQQ(5)-ampQQ(6)-ampQQ(9)+ampQQ(10)+ampQQ(12)
               color(4)=-ampQQ(7)-ampQQ(8)-ampQQ(11)
cmk            color(4)=-ampQQ(7)-ampQQ(8)+ampQQ(11)

               do i=1,4
                  colorT=DCMPLX(0d0,0d0)
                  do j=1,4
                     colorT=colorT+RGB(i,j)*color(j)
                  enddo
                  res12=res12+(dreal(colorT)*dreal(color(i))
     &                 + dimag(colorT)*dimag(color(i)))   
               enddo

c... store amplitudes for check with MadGraph
               if(MGdebug) then
                  do i=1,4
                     GGFLOcolor(i)=color(i)*dsqrt(fac)
cmk                     print*,'color(',i,')',color(i)*dsqrt(fac)
cmk                     STOP
                  enddo
                  do i=1,12
                     GGFLOamp1(isig1,isig3,ipol,i)=ampQQ(i)*dsqrt(fac)
                     GGFLOamp1(isig1,isig3,ipol,i+12)=DCMPLX(0d0,0d0)
                  enddo
c                  print*,'ipol,isig1,isig3',ipol,isig1,isig3
c                  do i=1,24
c                     print*,'GGFLOamp1(',i,')',GGFLOamp1(isig1,
c     &                    isig3,ipol,i)
c                  enddo
c                  print*,''
               endif
               
c... gauge check
               if (lgauge2) then
                  aux=0.0d0
                  if(Heff) then
                     do i=1,12
                        aux=aux+abs(ampHeff3j21(i))
                     enddo
                     do i=1,12
                        aux=aux+abs(ampHeff3j41(i))
                     enddo
                  endif
                  if(HloopTop) then
                     do i=1,12
                        aux=aux+abs(ampHtop3j21(i))
                     enddo
                     do i=1,12
                        aux=aux+abs(ampHtop3j41(i))
                     enddo
                  endif
                  print*,'ipol',ipol
                  print*,'aux',aux
                  print*,'norm. amp21i',amp21i/aux
                  print*,'norm. amp21f',amp21f/aux
                  print*,'norm. amp43i',amp43i/aux
                  print*,'norm. amp43f',amp43f/aux
               endif
           enddo
        enddo
      enddo
      
c... multiply by the coupling factor
      rh(2)= res12*fac

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
c 
c         Check whether there is instabilities
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!                                                       !! 
c!         IF SO ::::-->   Rh(i) -> 0                     !!
c!                                                       !!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
c
c         Count instabilites for hex,pen,box,vertices
c         both in quad and dble precision
c           
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c   Check Vertices
c     
c        if(boxeva) then
c   Dble precision
           If(verdbleSin) then
c          print*,"verdbleSi",verdbleSin
c          stop
c bad dble
        verinst(gf_iteration,1,1)=verinst(gf_iteration,1,1)+1
c        print*, " verinst(gf_iteration,1,1)", verinst(gf_iteration,1,1)
       else
c good dble 
        verinst(gf_iteration,2,1)=verinst(gf_iteration,2,1)+1
c        print*, " verinst(gf_iteration,2,1)", verinst(gf_iteration,2,1)
        endif

       If(verQUADSin) then
c bad QUAD
        verinst(gf_iteration,1,2)=verinst(gf_iteration,1,2)+1
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        rh(1)=0d0
        rh(2)=0d0
        rh(0)=0d0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccc
       else
c good QUAD 
        verinst(gf_iteration,2,2)=verinst(gf_iteration,2,2)+1
        endif
c        endif
c 
c   Check Boxes
c     
        if(boxeva) then
c   Dble precision
       If(boxdbleSin) then
c bad dble
        boxinst(gf_iteration,1,1)=boxinst(gf_iteration,1,1)+1
       else
c good dble 
        boxinst(gf_iteration,2,1)=boxinst(gf_iteration,2,1)+1
        endif

       If(boxQUADSin) then
c bad QUAD
        boxinst(gf_iteration,1,2)=boxinst(gf_iteration,1,2)+1
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        rh(1)=0d0
        rh(2)=0d0
        rh(0)=0d0        
c
cccccccccccccccccccccccccccccccccccccccccccccccccccc
       else
c good QUAD 
        boxinst(gf_iteration,2,2)=boxinst(gf_iteration,2,2)+1
        endif
        endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c... switch between (non-)identical quarks for check with MadGraph
        if(MGdebug) then
           if(mguu) then
              if(equal_f) then
                 GGFLOampSq1= rh(0) ! for u u
              else
                 GGFLOampSq1= rh(1) ! for u d 
              endif
           endif
           if(mguub) then
              if(equal_f) then
                 GGFLOampSq1= rh(0) ! for u u
              else
                 GGFLOampSq1= rh(2) ! for u d 
              endif
           endif
           if(mggu) then
              if(equal_f) then
                 GGFLOampSq1= rh(0) ! for u u
              else
                 GGFLOampSq1= rh(2) ! for u d 
              endif
           endif
        endif

c...  gauge check
        if (lgauge2) then
           print*,'rh(1)',rh(1)/aux*aux
           print*,'rh(2)',rh(2)/aux*aux
           print*,'rh(0)',rh(0)/aux*aux
           STOP
        endif


        END
