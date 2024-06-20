C========================================================================
C
C QQHQQa calculates the matrix elements for light Higgs production by
C electroweak boson fusion in quark quark scattering
C for anomalous higgs coulpings
C    q1 q3 ----> q2 q4 H;   H ---> b5-bar b6 optional
C
C QQHQQaI must be called first to initialize some couplings
C
C========================================================================
C
      subroutine qqhqqai
      implicit none
      
#include "VBFNLO/utilities/global.inc"
c      double precision  xmb
c      parameter ( xmb = 4.5d0 )
c
c the main task of QQHQQaI is to set up the products of coupling constants
c needed in Feynman graphs:
c     fcpla(isig,kgraph,icp) are couplings for ww and zz fusion
c     fcpla1..3(isig,kgraph,icp) are coupling for the 3
c     additional n.c. graphs. 
c     isig = 1,4 helicity configuration
c     kgraph = 1,6 is the graph id
c     icp = 1 for cp even coupling
c         = 2 for cp odd coulping
c
c       p1 ---------------------------p2
c                   |Z
c                   |------------
c                   |Z
c       p3 ---------------------------p4
c
      double precision  fcpl(4,6) ! hzz,hww
      double precision fcpl3(4,4) !hgg,hgz
      double complex fcpla(4,6,2) !hzz,hww 
      double complex fcpla1(4,4,2),fcpla2(4,4,2),fcpla3(4,4,2) !hgg,hgz
*sophy: fcplLa contains the loop induced value of a1: hzz/ww, hzg, hgz, hgg
      double complex fcplLa(4,6),fcplLa1(4,4),fcplLa2(4,4),fcplLa3(4,4)
c
c and which are calculated from output of KOPPLN
c
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)
#include "VBFNLO/utilities/koppln.inc"
c
c here fcpl(sig,i) contains the coupling constant factors. 
c     sig = 1,..,4    identifies the 4 helicity combinations of the quarks
c     i   = 1,..,6    numbers the possible subprocesses (4 NC and 2 CC)
c  
c the various processes are identified by the following
c
      integer  v2ini(6), v3ini(6)
c     For additional graphs in NC processes that arise 
c     from loop-induced couplings of the Higgs to v2 and v3 bosons
c
      integer v2ini1(4), v3ini1(4) ! V2 = Z V3 = gamma 
      integer v2ini2(4), v3ini2(4) ! V2 = gamma V3 = Z
      integer v2ini3(4), v3ini3(4) ! V2 = gamma V3 = gamma
      integer  fl1(6), fl3(6)
      integer  sig1(4,2),sig3(4,2)
      double complex fcplCT(4,6)
      
      data fl1 /3,3,4,4,3,4/
      data fl3 /3,4,3,4,4,3/

      save fl1,fl3
      common /cqqhqq/ fcplCT, fcpl, sig1,sig3, v2ini,v3ini
      common /cqqhqqa/ fcpla,fcpla1,fcpla2,fcpla3,
     $      v2ini1,v3ini1, v2ini2,v3ini2,v2ini3,v3ini3

c     ahvv(i,v1,v2) i = 1,2,3 are coefficients of eq 1 of 
c     Physics Letters B 591,297
      double complex ahvv(3,4,4), ahvvL(3,4,4), ahvvUSED(3,4,4)
      common/tensorhvv/ ahvv, ahvvL
c note: initialization for v2ini, v3ini, sig1,sig3 in block data
c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"

      real*8 c2,c2o2pi
      parameter (c2=4d0/3d0, c2o2pi=c2/2d0/pi)
c other local variables for QQBQQI
c
      integer  i,k,v2,v3, isig,isig1,isig3
c
c variables for the main part of the program
c

      double precision  p21(0:4), p43(0:4)

      double complex  j21(0:3,-1:1), j43(0:3,-1:1)
      double complex  j21j43(-1:1,-1:1), zmm(6,4,2),zmm1(6,4,2),
     $     zmm0(6,4,2),zmm2(6,4,2),zmm3(6,4,2)

      double complex j21j43a(-1:1,-1:1,2)
      double complex  prop21(1:4), prop43(1:4)

      integer ffac
      double precision m2ff, mff, psq
      double complex formfactor
      logical lff
      common/ formfacmass/ m2ff,mff,ffac,lff
      double complex  czero
      parameter ( czero = (0.0d0,0.0d0) )
      real*8 two
      parameter(two = 2.0d0)
 
      double complex I3point
      external I3point


** FERMION contribution: fm = 0 actual loop contn, = 1 fermion contn
      integer fm
** formfactors for FERMIONS.  Labelling: (i,2,2) are the z formfactors, 
** (i,3,4) and (i,4,3) are the w formfactors.
**    i     =     1      a1
**          =     2      a2
**          =     3      a3
      double complex ahvvFM(3,4,4)  
      common /fmFormFac/ ahvvFM


** setting couplings used in production: ew or anomalous
      do i = 1, 3
         do isig = 1, 4
            do isig1 = 1, 4
               if (ewcor_switch) then
                  ahvvUSED(i,isig,isig1) = ahvvL(i,isig,isig1)
               else
                  ahvvUSED(i,isig,isig1) = ahvv(i,isig,isig1)
               end if
            end do
         end do
      end do
     
c
c reset the coupling factors
c
      do i = 1,6
         do isig = 1,4
            fcpla(isig,i,1) = czero ! CP even
            fcpla(isig,i,2) = czero ! CP odd     
         enddo
      enddo
c
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


      do i = 1,6
         do isig = 1,4
            isig1 = sig1(isig,(i+3)/4)
            if ( isig1.ne.0 ) then
               isig3 = sig3(isig,(i+3)/4)
               v2 = v2ini(i)
               v3 = v3ini(i)
               fcpla(isig,i,1) = clr(fl1(i),v2,isig1)*
     &              clr(fl3(i),v3,isig3)*ahvvUSED(2,v2,v3)
               fcpla(isig,i,2) = clr(fl1(i),v2,isig1)*
     &              clr(fl3(i),v3,isig3)*ahvvUSED(3,v2,v3)
c     
            endif      
         enddo
      enddo
c     initialize additional NC graphs

      do i = 1,4
         do isig = 1,4
            isig1 = sig1(isig,(i+3)/4)
            if ( isig1.ne.0 ) then
               isig3 = sig3(isig,(i+3)/4)
               v2 = v2ini1(i)
               v3 = v3ini1(i)
               fcpla1(isig,i,1) = clr(fl1(i),v2,isig1)*
     &              clr(fl3(i),v3,isig3)*ahvvUSED(2,v2,v3)
               fcpla1(isig,i,2) = clr(fl1(i),v2,isig1)*
     &              clr(fl3(i),v3,isig3)*ahvvUSED(3,v2,v3)
               v2 = v2ini2(i)
               v3 = v3ini2(i)
               fcpla2(isig,i,1) = clr(fl1(i),v2,isig1)*
     &              clr(fl3(i),v3,isig3)*ahvvUSED(2,v2,v3)
               fcpla2(isig,i,2) = clr(fl1(i),v2,isig1)*
     &              clr(fl3(i),v3,isig3)*ahvvUSED(3,v2,v3)
               v2 = v2ini3(i)
               v3 = v3ini3(i)
               fcpla3(isig,i,1) = clr(fl1(i),v2,isig1)*
     &              clr(fl3(i),v3,isig3)*ahvvUSED(2,v2,v3)
               fcpla3(isig,i,2) = clr(fl1(i),v2,isig1)*
     &              clr(fl3(i),v3,isig3)*ahvvUSED(3,v3,v3)
               
            endif
         enddo
      enddo
      
      return
c
c entry for anomalous H production
c
      entry qqhqqa(j21,j43,p21,p43,prop21,prop43,  !input
     1             zmm,fm)                            ! output


** setting couplings used in production: ew or anomalous
      do i = 1, 3
         do isig = 1, 4
            do isig1 = 1, 4
               if (ewcor_switch) then
                  ahvvUSED(i,isig,isig1) = ahvvL(i,isig,isig1)
               else
                  ahvvUSED(i,isig,isig1) = ahvv(i,isig,isig1)
               end if
            end do
         end do
      end do


*sophy: If we're running the ew corns, clear and reset the couplings
      if (ewcor_switch) then
         do i = 1,6
            do isig = 1,4
               fcpla(isig,i,1) = czero ! CP even
               fcpla(isig,i,2) = czero ! CP odd    
*sophy: clearing loop induced a1
               fcplLa(isig,i) = czero
               if (i .lt. 5) then
                  fcplLa1(isig,i) = czero
                  fcplLa2(isig,i) = czero
                  fcplLa3(isig,i) = czero
               end if                  
            enddo
         enddo
c
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

         
         if (fm .eq. 0) then

*sophy: setting hww and hzz
            do i = 1,6
               do isig = 1,4
                  isig1 = sig1(isig,(i+3)/4)
                  if ( isig1.ne.0 ) then
                     isig3 = sig3(isig,(i+3)/4)
                     v2 = v2ini(i)
                     v3 = v3ini(i)
*sophy: anomalous a1
                     fcplLa(isig,i) = clr(fl1(i),v2,isig1)*
     &                    clr(fl3(i),v3,isig3)*ahvvUSED(1,v2,v3)
*sophy: anomalous a2
                     fcpla(isig,i,1) = clr(fl1(i),v2,isig1)*
     &                    clr(fl3(i),v3,isig3)*ahvvUSED(2,v2,v3)
*sophy: anomalous a3
                     fcpla(isig,i,2) = clr(fl1(i),v2,isig1)*
     &                    clr(fl3(i),v3,isig3)*ahvvUSED(3,v2,v3)
c     
c$$$                     write(*,*)'fm,isig,i,fcplLa(isig,i) =',
c$$$     -                    fm,isig,i,fcplLa(isig,i)
c$$$                     write(*,*)'fm,isig,i,fcpla(isig,i,1) =',
c$$$     -                    fm,isig,i,fcpla(isig,i,1)
c$$$                     write(*,*)'fm,isig,i,fcpla(isig,i,2) =',
c$$$     -                    fm,isig,i,fcpla(isig,i,2)
c$$$                     write(*,*)'ahvvUSED(1,v2,v3) =', ahvvUSED(1,v2,v3)
c$$$                     write(*,*)'ahvvUSED(2,v2,v3) =', ahvvUSED(2,v2,v3)
c$$$                     write(*,*)'ahvvUSED(3,v2,v3) =', ahvvUSED(3,v2,v3)
                  endif   
   
               enddo
            enddo


*sophy: setting additional NC graphs
            do i = 1,4
               do isig = 1,4
                  isig1 = sig1(isig,(i+3)/4)
                  if ( isig1.ne.0 ) then
                     isig3 = sig3(isig,(i+3)/4)
*sophy: h-z-photon
                     v2 = v2ini1(i)
                     v3 = v3ini1(i)
                     fcplLa1(isig,i) = clr(fl1(i),v2,isig1)*
     &                    clr(fl3(i),v3,isig3)*ahvvUSED(1,v2,v3)
                     fcpla1(isig,i,1) = clr(fl1(i),v2,isig1)*
     &                    clr(fl3(i),v3,isig3)*ahvvUSED(2,v2,v3)
                     fcpla1(isig,i,2) = clr(fl1(i),v2,isig1)*
     &                    clr(fl3(i),v3,isig3)*ahvvUSED(3,v2,v3)
                     
*sophy: h-photon-z
                     v2 = v2ini2(i)
                     v3 = v3ini2(i)
                     fcplLa2(isig,i) = clr(fl1(i),v2,isig1)*
     &                    clr(fl3(i),v3,isig3)*ahvvUSED(1,v2,v3)
                     fcpla2(isig,i,1) = clr(fl1(i),v2,isig1)*
     &                    clr(fl3(i),v3,isig3)*ahvvUSED(2,v2,v3)
                     fcpla2(isig,i,2) = clr(fl1(i),v2,isig1)*
     &                    clr(fl3(i),v3,isig3)*ahvvUSED(3,v2,v3)
                        
*sophy: h-photon-photon
                     v2 = v2ini3(i)
                     v3 = v3ini3(i)
                     fcplLa3(isig,i) = clr(fl1(i),v2,isig1)*
     &                    clr(fl3(i),v3,isig3)*ahvvUSED(1,v2,v3)
                     fcpla3(isig,i,1) = clr(fl1(i),v2,isig1)*
     &                    clr(fl3(i),v3,isig3)*ahvvUSED(2,v2,v3)
                     fcpla3(isig,i,2) = clr(fl1(i),v2,isig1)*
     &                    clr(fl3(i),v3,isig3)*ahvvUSED(3,v3,v3)
                  
                  endif
               enddo
            enddo
*sophy: end of reinitialisation of couplings for ew corns

         else if (fm .eq. 1) then

* START of FM initialisation!
*sophy: setting hww and hzz
            do i = 1,6
               do isig = 1,4
                  isig1 = sig1(isig,(i+3)/4)
                  if ( isig1.ne.0 ) then
                     isig3 = sig3(isig,(i+3)/4)
                     v2 = v2ini(i)
                     v3 = v3ini(i)
*sophy: anomalous a1
                     fcplLa(isig,i) = clr(fl1(i),v2,isig1)*
     &                    clr(fl3(i),v3,isig3)*ahvvFM(1,v2,v3)
*sophy: anomalous a2
                     fcpla(isig,i,1) = clr(fl1(i),v2,isig1)*
     &                    clr(fl3(i),v3,isig3)*ahvvFM(2,v2,v3)
*sophy: anomalous a3
                     fcpla(isig,i,2) = clr(fl1(i),v2,isig1)*
     &                    clr(fl3(i),v3,isig3)*ahvvFM(3,v2,v3)

c$$$                     write(*,*)'fm,isig,i,fcplLa(isig,i) =',
c$$$     -                    fm,isig,i,fcplLa(isig,i)
c$$$                     write(*,*)'fm,isig,i,fcpla(isig,i,1) =',
c$$$     -                    fm,isig,i,fcpla(isig,i,1)
c$$$                     write(*,*)'fm,isig,i,fcpla(isig,i,2) =',
c$$$     -                    fm,isig,i,fcpla(isig,i,2)
c$$$                     write(*,*)'ahvvFM(1,v2,v3) =', ahvvFM(1,v2,v3)
c$$$                     write(*,*)'ahvvFM(2,v2,v3) =', ahvvFM(2,v2,v3)
c$$$                     write(*,*)'ahvvFM(3,v2,v3) =', ahvvFM(3,v2,v3)

c     
                  endif      
               enddo
            enddo


*sophy: setting additional NC graphs
            do i = 1,4
               do isig = 1,4
                  isig1 = sig1(isig,(i+3)/4)
                  if ( isig1.ne.0 ) then
                     isig3 = sig3(isig,(i+3)/4)
*sophy: h-z-photon
                     v2 = v2ini1(i)
                     v3 = v3ini1(i)
                     fcplLa1(isig,i) = clr(fl1(i),v2,isig1)*
     &                    clr(fl3(i),v3,isig3)*ahvvFM(1,v2,v3)
                     fcpla1(isig,i,1) = clr(fl1(i),v2,isig1)*
     &                    clr(fl3(i),v3,isig3)*ahvvFM(2,v2,v3)
                     fcpla1(isig,i,2) = clr(fl1(i),v2,isig1)*
     &                    clr(fl3(i),v3,isig3)*ahvvFM(3,v2,v3)
                     
*sophy: h-photon-z
                     v2 = v2ini2(i)
                     v3 = v3ini2(i)
                     fcplLa2(isig,i) = clr(fl1(i),v2,isig1)*
     &                    clr(fl3(i),v3,isig3)*ahvvFM(1,v2,v3)
                     fcpla2(isig,i,1) = clr(fl1(i),v2,isig1)*
     &                    clr(fl3(i),v3,isig3)*ahvvFM(2,v2,v3)
                     fcpla2(isig,i,2) = clr(fl1(i),v2,isig1)*
     &                    clr(fl3(i),v3,isig3)*ahvvFM(3,v2,v3)

*sophy: h-photon-photon
                     v2 = v2ini3(i)
                     v3 = v3ini3(i)
                     fcplLa3(isig,i) = clr(fl1(i),v2,isig1)*
     &                    clr(fl3(i),v3,isig3)*ahvvFM(1,v2,v3)
                     fcpla3(isig,i,1) = clr(fl1(i),v2,isig1)*
     &                    clr(fl3(i),v3,isig3)*ahvvFM(2,v2,v3)
                     fcpla3(isig,i,2) = clr(fl1(i),v2,isig1)*
     &                    clr(fl3(i),v3,isig3)*ahvvFM(3,v3,v3)
                     
                  endif
               enddo
            enddo
            
         end if  ! end of fm = 0/1


      end if
*sophy: end of FM reinitialisation of couplings for ew corns


c initialization
      do isig1 = -1,1,2
         do isig3 = -1,1,2
            j21j43a(isig1,isig3,1) = czero
            j21j43a(isig1,isig3,2) = czero
         end do
      end do

c   
      if(lff) then
         if(ffac.eq.1) then
            formfactor = M2ff/(-p21(4)+M2ff)*M2ff/(-p43(4)+M2ff)
         else if(ffac.eq.2) then
            psq = (p21(0)+p43(0))**2-(p21(1)+p43(1))**2-(p21(2)+p43(2))**2
     &           -(p21(3)+p43(3))**2
            formfactor = -2*m2ff*I3point(mff,p21(4),p43(4),psq)
         else
            write(*,*) "Please select ffac=1 or ffac=2 in anom_HVV.dat"
            stop
         endif
      endif

*sophy: get the dot products of the currents for a1
      do isig1 = -1,1,2
         do isig3 = -1,1,2
            j21j43(isig1,isig3) = j21(0,isig1)*j43(0,isig3) -
     &                            j21(1,isig1)*j43(1,isig3) -
     &                            j21(2,isig1)*j43(2,isig3) -
     &                            j21(3,isig1)*j43(3,isig3)
         enddo
      enddo

      do isig1 = -1,1,2
         do isig3 = -1,1,2
            call contract_CPO(j21j43a(isig1,isig3,2),j21(0,isig1),
     $           j43(0,isig3),p21,p43)
            if(lff) then
               j21j43a(isig1,isig3,2) = 
     1           formfactor *j21j43a(isig1,isig3,2)
            endif
c     
c      print*,j21j43a(isig1,isig3,2)
            call contract_CPE(j21j43a(isig1,isig3,1),j21(0,isig1),
     $           j43(0,isig3),p21,p43)
            if(lff) then
               j21j43a(isig1,isig3,1)=formfactor*j21j43a(isig1,isig3,1)
            endif
c      print*,j21j43a(isig1,isig3,1)
         enddo
      enddo
c

c
c now get the coupling*propagator factors for subprocess k, helicity
c combination isig
c     zmm0..3(k,isig,i)
c     0 is zz fusion
c     1 and 2 are z - gamma fusion
c     3 is gamma-gamma fusion
c     k=1..4 refer to the flavor configuration
c     isig is the helicity configuration
c     i = 1 (no gluon emission here)

      i=1
      do k = 1,4
     
         do isig = 1,4
            isig1 = sig1(isig,1)
            isig3 = sig3(isig,1)
            zmm0(k,isig,i) = fcpla(isig,k,1)*j21j43a(isig1,isig3,1)+
     $           fcpla(isig,k,2)*j21j43a(isig1,isig3,2)         ! hzz
c
            zmm1(k,isig,i) = fcpla1(isig,k,1)*j21j43a(isig1,isig3,1)+
     $                fcpla1(isig,k,2)*j21j43a(isig1,isig3,2)   ! hzg
c
            zmm2(k,isig,i) = fcpla2(isig,k,1)*j21j43a(isig1,isig3,1)+
     $               fcpla2(isig,k,2)*j21j43a(isig1,isig3,2)    ! hgz
c
            zmm3(k,isig,i) = fcpla3(isig,k,1)*j21j43a(isig1,isig3,1)+
     $               fcpla3(isig,k,2)*j21j43a(isig1,isig3,2)    ! hgg

*sophy: adding in anomalous (i.e. loop-induced) a1 terms
            if (ewcor_switch) then
               zmm0(k,isig,i) = zmm0(k,isig,i) + 
     $              fcplLa(isig,k)*j21j43(isig1,isig3)          ! hzz
               zmm1(k,isig,i) = zmm1(k,isig,i) + 
     $              fcplLa1(isig,k)*j21j43(isig1,isig3)         ! hzg
               zmm2(k,isig,i) = zmm2(k,isig,i) + 
     $              fcplLa2(isig,k)*j21j43(isig1,isig3)         ! hgz
               zmm3(k,isig,i) = zmm3(k,isig,i) + 
     $              fcplLa3(isig,k)*j21j43(isig1,isig3)         ! hgg
            end if

            zmm0(k,isig,i) = prop21(2)*prop43(2)*zmm0(k,isig,i) ! hzz
            zmm1(k,isig,i) = prop21(2)*prop43(1)*zmm1(k,isig,i) ! hzg
            zmm2(k,isig,i) = prop21(1)*prop43(2)*zmm2(k,isig,i) ! hgz
            zmm3(k,isig,i) = prop21(1)*prop43(1)*zmm3(k,isig,i) ! hgg
c
            zmm(k,isig,i) = zmm0(k,isig,i) + zmm1(k,isig,i)+ 
     $             zmm2(k,isig,i) + zmm3(k,isig,i)
 
         enddo
      enddo

      do k = 5,6

         zmm(k,1,i) = fcpla(1,k,1)*j21j43a(-1,-1,1)+
     $        fcpla(1,k,2)*j21j43a(-1,-1,2)
*sophy: adding in anomalous (i.e. loop-induced) a1 term
         if (ewcor_switch) then
            zmm(k,1,i) = zmm(k,1,i) + fcplLa(1,k)*j21j43(-1,-1)
         end if

         zmm(k,1,i) =  prop21(3)*prop43(3)*zmm(k,1,i)

      enddo

c
      return
      end

C
C********************  qqhqq ***********************************
C
      subroutine qqhqqgai 
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
c      real*8 fcpl(4,6)
c      double complex fcpla(4,6,2)
      double precision  fcpl(4,6) ! hzz,hww
      double precision fcpl3(4,4) !hgg,hgz
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
      double complex fcplCT(4,6)
      
      data fl1 /3,3,4,4,3,4/
      data fl3 /3,4,3,4,4,3/

      save fl1,fl3
     
      common /cqqhqq/ fcplCT, fcpl,sig1,sig3, v2ini,v3ini
      common /cqqhqqa/ fcpla,fcpla1,fcpla2,fcpla3,
     $      v2ini1,v3ini1, v2ini2,v3ini2,v2ini3,v3ini3

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
     2        e21(0:3,-1:1,2),
     3        e43(0:3,-1:1,2),
     4        mm21(6,4,2), mm43(6,4,2)
      complex*16 mm211(4,4,2), mm431(4,4,2), mm212(4,4,2), mm432(4,4,2),
     $      mm213(4,4,2), mm433(4,4,2)
      complex*16 j21e43a(-1:1,-1:1,2,2),e21j43a(-1:1,-1:1,2,2)

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

c
C  Reset the coupling factors
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
            endif
         enddo
      enddo
      
c      do i = 1,6
c         do isig = 1,4
c            print*,"fcpl",fcpl(isig,i)
c            print*,"fcpla",fcpla(isig,i,1) ! CP even
c            print*,"fcpla",fcpla(isig,i,2) ! CP odd     
c         enddo
c      enddo
      return
              
C  Entry for anomalous H production
c 
      entry qqHqqja(j21,j43,e21,e43,p21,p43,p21g,
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

c      do isig1 = -1,1,2
c         do isig3 = -1,1,2
c            do i =1,2
c               print*,e21j43(isig1,isig3,i) 
c               print*,j21e43(isig1,isig3,i) 
c               print*,e21j43a(isig1,isig3,i,1)  !CP even
c               print*,e21j43a(isig1,isig3,i,2)  !CP odd
c               print*,j21e43a(isig1,isig3,i,1) 
c               print*,j21e43a(isig1,isig3,i,2) 
c            enddo  
c         enddo
c      enddo    
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
                  mm21(k,isig,i) = 0
                  mm43(k,isig,i) = 0
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
                  mm211(k,isig,i) = 0
                  mm431(k,isig,i) = 0
                  mm212(k,isig,i) = 0
                  mm432(k,isig,i) = 0
                  mm213(k,isig,i) = 0
                  mm433(k,isig,i) = 0
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
