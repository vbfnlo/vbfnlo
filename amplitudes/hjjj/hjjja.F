C********************  qqhqq ***********************************
      subroutine hjjjai 
         implicit none                
#include "VBFNLO/utilities/global.inc"
         double precision  fpi, xmb
         parameter ( fpi=4d0*pi, xmb=4.5d0 )
C
C  QQHQQ calculates the matrix elements**2 for light Higgs production by
C  electroweak boson fusion in quark quark scattering
C
C        q1 q3    ---->   q2 q4 g H
C
C  QQHQQI must be called first to initialize some couplings
C
C  The main task of QQHQQgaI is to set up the products of coupling constants
C  needed in Feynman graphs a, ..., g which are stored in

         double precision  fcpl(4,6) ! hzz,hww
         double complex fcpla(4,6,2) !hzz,hww
         double complex fcpla1(4,4,2),fcpla2(4,4,2),fcpla3(4,4,2) !hgg,hgz
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
c     For additional graphs in NC processes that arise 
c     from loop-induced couplings of the Higgs to v2 and v3 bosons
c
         integer v2ini1(4), v3ini1(4) ! V2 = Z V3 = gamma 
         integer v2ini2(4), v3ini2(4) ! V2 = gamma V3 = Z
         integer v2ini3(4), v3ini3(4) ! V2 = gamma V3 = gamma
         integer fl1(6), fl3(6)
         integer sig1(4,2),sig3(4,2)
         
         data fl1 /3,3,4,4,3,4/
         data fl3 /3,4,3,4,4,3/

         save fl1,fl3
      
         common /cqqhqqj/ fcpl,sig1,sig3, v2ini,v3ini
      !   common /cqqhqqja/ fcpla,fcpla1,fcpla2,fcpla3

c     ahvv(i,v1,v2) i = 1,2,3 are coefficients of eq 1 of 
c     Physics Letters B 591,297
         double complex ahvv(3,4,4), ahvvL(3,4,4)
         common/tensorhvv/ ahvv, ahvvL
c note: initialization for v2ini, v3ini, sig1,sig3 in block data
C
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"

C  Other local variables for QQBQQI
         integer i,k,v2,v3, isig,isig1,isig3
C         
C  Variables for the main part of the program
         real*8 
     1       p21(0:4),p43(0:4),
     2       p21g(0:4),p43g(0:4)


         complex*16 
     1        j21(0:3,-1:1), j43(0:3,-1:1),
     2        e21(0:3,-1:1,0:2),
     3        e43(0:3,-1:1,0:2),
     4        mm21(6,4,0:2), mm43(6,4,0:2)
         complex*16 mm211(6,4,0:2), mm431(6,4,0:2), mm212(6,4,0:2), mm432(6,4,0:2),
     $      mm213(6,4,0:2), mm433(6,4,0:2)
         complex*16 j21e43a(-1:1,-1:1,0:2,2),e21j43a(-1:1,-1:1,0:2,2)

         complex*16 prop21(1:4), prop43(1:4), prop21g(1:4), prop43g(1:4)
         real*8 xmw
         integer ffac
         double complex formfactor
         double precision mff, m2ff, psq
         logical lff
         common/ formfacmass/ m2ff,mff,ffac,lff
         double complex  czero
         parameter ( czero = (0.0d0,0.0d0) )
         real*8 two
         parameter(two = 2.0d0)


         double complex I3point
         external I3point

c set constant to identify the boson v2,v3
         v2ini1=(/2,2,2,2/)
         v3ini1= (/1,1,1,1/)
         v2ini2=(/1,1,1,1/) 
         v3ini2=(/2,2,2,2/)
         v2ini3=(/1,1,1,1/)
         v3ini3=(/1,1,1,1/)
c---------------------------------------------------------------

C  Reset the coupling factors
c-----------------------------------------------
c  set the SM coupling equal to ahvv(1,v2,v3)
      do i = 1,6
         do isig = 1,4
            isig1 = sig1(isig,(i+3)/4)
            if ( isig1.ne.0 ) then
               isig3 = sig3(isig,(i+3)/4)
               v2 = v2ini(i)
               v3 = v3ini(i)
               fcpl(isig,i) = clr(fl1(i),v2,isig1)*
     &              clr(fl3(i),v3,isig3)*
     &              ahvv(1,v2,v3)
            endif
         enddo
      enddo 
c----------------------------------------------------------------
         do i = 1,6
            do isig = 1,4
               fcpla(isig,i,1)=czero
               fcpla(isig,i,2)=czero
            enddo
         enddo
         do i = 1,4
            do isig = 1,4
         
               fcpla1(isig,i,1) = czero ! CP even
               fcpla1(isig,i,2) = czero ! CP odd
         
               fcpla2(isig,i,1) = czero ! CP even
               fcpla2(isig,i,2) = czero ! CP odd
         
               fcpla3(isig,i,1) = czero ! CP even
               fcpla3(isig,i,2) = czero ! CP odd     
            enddo
         enddo
c
         xmw = sqrt(xm2(3))
         do i = 1,6
            do isig = 1,4
               isig1 = sig1(isig,(i+3)/4)
               if ( isig1.ne.0 ) then
                  isig3 = sig3(isig,(i+3)/4)
                  v2 = v2ini(i)
                  v3 = v3ini(i)
                  
                  fcpla(isig,i,1) = clr(fl1(i),v2,isig1)*
     &              clr(fl3(i),v3,isig3)*ahvv(2,v2,v3)
                  fcpla(isig,i,2) = clr(fl1(i),v2,isig1)*
     &              clr(fl3(i),v3,isig3)*ahvv(3,v2,v3)
               endif         
            enddo
         enddo
c     initialize additional graphs
c     initialize additional NC graphs

         do i = 1,4
            do isig = 1,4
               isig1 = sig1(isig,(i+3)/4)
               if ( isig1.ne.0 ) then
                  isig3 = sig3(isig,(i+3)/4)
                  v2 = v2ini1(i)
                  v3 = v3ini1(i)
                  fcpla1(isig,i,1) = clr(fl1(i),v2,isig1)*
     &              clr(fl3(i),v3,isig3)*ahvv(2,v2,v3)
                  fcpla1(isig,i,2) = clr(fl1(i),v2,isig1)*
     &              clr(fl3(i),v3,isig3)*ahvv(3,v2,v3)
                  v2 = v2ini2(i)
                  v3 = v3ini2(i)
                  fcpla2(isig,i,1) = clr(fl1(i),v2,isig1)*
     &              clr(fl3(i),v3,isig3)*ahvv(2,v2,v3)
                  fcpla2(isig,i,2) = clr(fl1(i),v2,isig1)*
     &              clr(fl3(i),v3,isig3)*ahvv(3,v2,v3)
                  v2 = v2ini3(i)
                  v3 = v3ini3(i)
                  fcpla3(isig,i,1) = clr(fl1(i),v2,isig1)*
     &              clr(fl3(i),v3,isig3)*ahvv(2,v2,v3)
                  fcpla3(isig,i,2) = clr(fl1(i),v2,isig1)*
     &              clr(fl3(i),v3,isig3)*ahvv(3,v2,v3)
               else
                  fcpla1(isig,i,1) = czero ! CP even
                  fcpla1(isig,i,2) = czero ! CP odd
            
                  fcpla2(isig,i,1) = czero ! CP even
                  fcpla2(isig,i,2) = czero ! CP odd
            
                  fcpla3(isig,i,1) = czero ! CP even
                  fcpla3(isig,i,2) = czero ! CP odd 
               endif
            enddo
         enddo
         
!       do i = 1,4
!          do isig = 1,4
!  !           print*,"fcpl",fcpl(isig,i)
!             print*,"fcpla1",fcpla1(isig,i,1) ! CP even
!             print*,"fcpla2",fcpla2(isig,i,1) ! CP odd
!             print*,"fcpla3",fcpla3(isig,i,1)       
!          enddo
!       enddo
       return
               
C  Entry for anomalous H production
c 
      entry hjjja(j21,j43,e21,e43,p21,p43,p21g,
     $     p43g,prop21,prop43,
     $     prop21g,prop43g,                       ! input
     $     mm21,mm43)                             !output
   
c  initialization
         do isig1 = -1,1,2
            do isig3 = -1,1,2
               do i =1,2
                  e21j43a(isig1,isig3,i,1) = czero !CP even
                  e21j43a(isig1,isig3,i,2) = czero !CP odd
                  j21e43a(isig1,isig3,i,1) = czero
                  j21e43a(isig1,isig3,i,2) = czero
               enddo  
            enddo
         enddo    
c
C  get the contractions of the currents with cp even and odd 
c  tensor structures for the 4 helicity combinations
c     e21j43a(isig1,isig3,i,icp)
c     e43j21a(isig1,isig3,i,icp)
c     i=1..2 is the polarization of the gluon
c     icp=1..2 is contraction of currents with cp odd tensor structure
c

         do isig1 = -1,1,2
            do isig3 = -1,1,2
               do i = 1,2
            call contract_CPO(e21j43a(isig1,isig3,i,2),e21(0,isig1,i),
     $              j43(0,isig3),p21g,p43)
                  if(lff) then
                     if(ffac.eq.1) then
                        formfactor =M2ff/(p21g(4)-M2ff)*M2ff/(p43(4)-M2ff)
                     else if(ffac.eq.2) then
                        psq = (p21g(0)+p43(0))**2-(p21g(1)+p43(1))**2
     &                    -(p21g(2)+p43(2))**2-(p21g(3)+p43(3))**2
                        formfactor = -2*m2ff*I3point(mff,p21g(4),p43(4),psq)
                     else
                        write(*,*)"Please select ffac=1or2 in anom_HVV.dat"
                        stop
                     endif                   
                     e21j43a(isig1,isig3,i,2)=
     1                 e21j43a(isig1,isig3,i,2)*formfactor
                  endif
c     
                  call contract_CPO(j21e43a(isig1,isig3,i,2),j21(0,isig1),
     $              e43(0,isig3,i),p21,p43g)
                  if(lff) then
                     if(ffac.eq.1) then
                        formfactor =M2ff/(p21(4)-M2ff)*M2ff/(p43g(4)-M2ff)
                     else if(ffac.eq.2) then
                        psq = (p21(0)+p43g(0))**2-(p21(1)+p43g(1))**2
     &                    -(p21(2)+p43g(2))**2-(p21(3)+p43g(3))**2
                        formfactor = -2*m2ff*I3point(mff,p21(4),p43g(4),psq)
                     endif
                     j21e43a(isig1,isig3,i,2) = 
     1                 formfactor*j21e43a(isig1,isig3,i,2)
                  endif
c     
               call contract_CPE(e21j43a(isig1,isig3,i,1),e21(0,isig1,i),
     $              j43(0,isig3),p21g,p43)
                  if(lff) then
                     if(ffac.eq.1) then
                        formfactor =M2ff/(p21g(4)-M2ff)*M2ff/(p43(4)-M2ff)
                     else if(ffac.eq.2) then
                        psq = (p21g(0)+p43(0))**2-(p21g(1)+p43(1))**2
     &                    -(p21g(2)+p43(2))**2-(p21g(3)+p43(3))**2
                        formfactor = -2*m2ff*I3point(mff,p21g(4),p43(4),psq)
                     endif                     
                     e21j43a(isig1,isig3,i,1) = 
     1                 formfactor*e21j43a(isig1,isig3,i,1)
                  endif
c     
                  call contract_CPE(j21e43a(isig1,isig3,i,1),j21(0,isig1),
     $              e43(0,isig3,i),p21,p43g)
                  if(lff) then
                     if(ffac.eq.1) then
                        formfactor =M2ff/(p21(4)-M2ff)*M2ff/(p43g(4)-M2ff)
                     else if(ffac.eq.2) then
                        psq = (p21(0)+p43g(0))**2-(p21(1)+p43g(1))**2
     &                    -(p21(2)+p43g(2))**2-(p21(3)+p43g(3))**2
                        formfactor = -2*m2ff*I3point(mff,p21(4),p43g(4),psq)
                     endif     
                     j21e43a(isig1,isig3,i,1) = 
     1                 formfactor*j21e43a(isig1,isig3,i,1)
                  endif
               enddo
            enddo
         enddo

!       do isig1 = -1,1,2
!          do isig3 = -1,1,2
!             do i =1,2
! !               print*,e21j43(isig1,isig3,i) 
! !               print*,j21e43(isig1,isig3,i) 
!                print*,e21j43a(isig1,isig3,i,1)  !CP even
!                print*,e21j43a(isig1,isig3,i,2)  !CP odd
!                print*,j21e43a(isig1,isig3,i,1) 
!                print*,j21e43a(isig1,isig3,i,2) 
!             enddo  
!          enddo
!       enddo    
C
C  now get the coupling*propagator factors for subprocess k, helicity
C  combination isig, gluon polarization i
c     mm21(k,isig,i) gluon emission off line 2-1
c     mm43           gluon emission off line 4-3 for zz and ww fusion
c     mm211..3 
c     mm431..3 are z-gamma and gamma-gamma fusion contributions to the
c     amplitude

         do k = 1,6
            do isig = 1,4
               isig1 = sig1(isig,(k+3)/4)
               if (isig1.ne.0) then
                  isig3 = sig3(isig,(k+3)/4)
                  v2 = v2ini(k)
                  v3 = v3ini(k)
                  do i = 1,2
                     mm21(k,isig,i) = 
     1                 fcpla(isig,k,1) * e21j43a(isig1,isig3,i,1) +
     2                 fcpla(isig,k,2) * e21j43a(isig1,isig3,i,2)
                  mm21(k,isig,i) = prop21g(v2)*prop43(v3)*mm21(k,isig,i)
c
                     mm43(k,isig,i) = 
     1                 fcpla(isig,k,1)* j21e43a(isig1,isig3,i,1) +
     2                 fcpla(isig,k,2)* j21e43a(isig1,isig3,i,2)
c
                  mm43(k,isig,i) = prop21(v2)*prop43g(v3)*mm43(k,isig,i)
                  enddo
               else
                  do i = 1,2
                     mm21(k,isig,i) = czero
                     mm43(k,isig,i) = czero
                  enddo
               endif
            enddo
         enddo
c     additional NC process that arise from loop-induced couplings of the higgs
c     
         do k = 1,4
            do isig = 1,4
               isig1 = sig1(isig,(k+3)/4)
               if (isig1.ne.0) then
                  isig3 = sig3(isig,(k+3)/4)
                  v2 = v2ini1(k)
                  v3 = v3ini1(k)
                  do i = 1,2
                     mm211(k,isig,i) = 
     2                 fcpla1(isig,k,1) * e21j43a(isig1,isig3,i,1)+
     3                 fcpla1(isig,k,2) * e21j43a(isig1,isig3,i,2)
                     mm211(k,isig,i) = 
     1                 prop21g(v2)*prop43(v3)*mm211(k,isig,i)
c
                     mm431(k,isig,i) = 
     2                 fcpla1(isig,k,1)* j21e43a(isig1,isig3,i,1)+
     3                 fcpla1(isig,k,2)* j21e43a(isig1,isig3,i,2)
c
                     mm431(k,isig,i) = 
     1                 prop21(v2)*prop43g(v3)*mm431(k,isig,i)
                  enddo
                  v2 = v2ini2(k)
                  v3 = v3ini2(k)
                  do i = 1,2
                     mm212(k,isig,i) = 
     2                 fcpla2(isig,k,1) * e21j43a(isig1,isig3,i,1)+
     3                 fcpla2(isig,k,2) * e21j43a(isig1,isig3,i,2)
                     mm212(k,isig,i) = 
     1                 prop21g(v2)*prop43(v3)*mm212(k,isig,i)
c
                     mm432(k,isig,i) = 
     2                 fcpla2(isig,k,1)* j21e43a(isig1,isig3,i,1)+
     3                 fcpla2(isig,k,2)* j21e43a(isig1,isig3,i,2)
c
                     mm432(k,isig,i) = 
     1                 prop21(v2)*prop43g(v3)*mm432(k,isig,i)
                  enddo
                  v2 = v2ini3(k)
                  v3 = v3ini3(k)
                  do i = 1,2
                     mm213(k,isig,i) = 
     2                 fcpla3(isig,k,1) * e21j43a(isig1,isig3,i,1)+
     3                 fcpla3(isig,k,2) * e21j43a(isig1,isig3,i,2)
                     mm213(k,isig,i) = 
     1                 prop21g(v2)*prop43(v3)*mm213(k,isig,i)
c
                     mm433(k,isig,i) = 
     2                 fcpla3(isig,k,1)* j21e43a(isig1,isig3,i,1)+
     3                 fcpla3(isig,k,2)* j21e43a(isig1,isig3,i,2)
c
                     mm433(k,isig,i) = 
     1                 prop21(v2)*prop43g(v3)*mm433(k,isig,i)
                  enddo
               else
                  do i = 1,2
                     mm211(k,isig,i) = czero
                     mm431(k,isig,i) = czero
                     mm212(k,isig,i) = czero
                     mm432(k,isig,i) = czero
                     mm213(k,isig,i) = czero
                     mm433(k,isig,i) = czero
                  enddo
               endif
            enddo
         enddo      
C
C  add additional n.c.contributions

         do k = 1,6
            do isig = 1,4
               do i = 1,2
                  if(k.lt.5) then  !NC additional graphs added in
                      mm21(k,isig,i) = mm21(k,isig,i)+mm211(k,isig,i)+
     $                 mm212(k,isig,i)+mm213(k,isig,i)
                      mm43(k,isig,i) = mm43(k,isig,i)+mm431(k,isig,i)+
     $                 mm432(k,isig,i)+mm433(k,isig,i)
                  endif
               enddo
            enddo            
         enddo
          
         return
         end
