cfc       complex*16 FUNCTION FD3gTri(i1,i2,i3,i4,i5)               
cfc       complex*16 FUNCTION FD1ePV3g(i1,i2,i3,i4,i5,virtual)    
cfc       complex*16 FUNCTION FD2gBox(i1,i2,i3,i4,i5,virtual)
cfc       complex*16 FUNCTION FDePV4g(i1,i2,i3,i4,i5,gp)
cfc       complex*16 FUNCTION FD4gBox(i1,i2,i3,i4,i5)
cfc       complex*16 FUNCTION FDeffPen(i1,i2,i3,i4,i5)
cfc       complex*16 FUNCTION FDeffPenLoop(i1,i2,i3,i4)
cfc       complex*16 FUNCTION FDHexLoop(i1,i2,i3,i4,i5)              
***********************************************************************
      complex*16 FUNCTION FD3gTri(i1,i2,i3,i4,i5)
*----------------------------------------------------------------------
c     3 three-gluon-vertices and effective polarization vector
***********************************************************************
      
      IMPLICIT NONE

#include "ggf_process.inc"

      integer i1,i2,i3,i4,i5

      integer i,mu
      complex*16 aux(26)
      
      real*8 K(0:4,5),Kij(0:4,5,5),Kijk(0:4,5,5,5),PH(0:4),KPH(0:4,5),
     &     epsg(0:3,5),
     &     KPHeffij(0:4,5,5),KPAeffij(0:4,5,5),KPHtopij(0:4,5,5),
     &     KPHbottomij(0:4,5,5),KPAtopij(0:4,5,5),KPAbottompij(0:4,5,5),
     &     Kpen(0:4,4)
c... 4-dim momentum variables
      real*8 Kp(0:3,5),Kijp(0:3,5,5),Kijkp(0:3,5,5,5),PHp(0:3),
     &     KPHp(0:3,5),K123p(0:3),KPHeffijp(0:3,5,5),KPAeffijp(0:3,5,5),
     &     KPHtopijp(0:3,5,5),KPHbottomijp(0:3,5,5),KPAtopijp(0:3,5,5),
     &     KPAbottomijp(0:3,5,5)
      complex*16 epsgC(0:3,5),Jij(0:3,5,5),Jijk(0:3,5,5,5),
     &           epsgPen(0:3,4)
      complex*16 epsHeff(0:3,5),epsAeff(0:3,5),epsHtop(0:3,5),
     &     epsHbottom(0:3,5),epsAtop(0:3,5),epsAbottom(0:3,5),
     &     JPHeffij(0:3,5,5),JPAeffij(0:3,5,5),JPHtopij(0:3,5,5),
     &     JPHbottomij(0:3,5,5),JPAtopij(0:3,5,5),JPAbottomij(0:3,5,5)
      COMMON /GG3JET/ K,Kij,Kijk,epsg,epsgC,Jij,Jijk,KPH,PH,
     &     KPHeffij,KPAeffij,KPHtopij,KPHbottomij,KPAtopij,
     &     KPAbottompij,
     &     epsHeff,epsAeff,epsHtop,epsHbottom,epsAtop,epsAbottom,
     &     JPHeffij,JPAeffij,JPHtopij,JPHbottomij,JPAtopij,
     &     JPAbottomij,
c... 4-dim momentum variables
     &     Kp,Kijp,Kijkp,PHp,KPHp,KPHeffijp,KPAeffijp,K123p,
     &     KPHtopijp,KPHbottomijp,KPAtopijp,KPAbottomijp

      Common/Forpent/Kpen,epsgPen
     
      complex*16 H2G_effective,A2G_effective,H2G_Loop,A2G_Loop
      EXTERNAL H2G_effective,A2G_effective,H2G_Loop,A2G_Loop

      FD3gTri=DCMPLX(0d0,0d0)
      do i=1,26
         aux(i)=DCMPLX(0d0,0d0)
      enddo

      if (Heff) then
         aux(1)=cplModHeff
     &     *H2G_effective(Kijkp(0,i1,i2,i3),Kijp(0,i4,i5),
     &                    Jijk(0,i1,i2,i3),Jij(0,i4,i5))
      endif
*----------------------------------------------------------------------
      if (Aeff) then
         aux(2)=cplModAeff
     &     *A2G_effective(Kijkp(0,i1,i2,i3),Kijp(0,i4,i5),
     &                    Jijk(0,i1,i2,i3),Jij(0,i4,i5))
      endif
*----------------------------------------------------------------------
      if (HloopTop) then
         aux(3)=cplModHtop
     &     *H2G_Loop(Kijk(0,i1,i2,i3),Kij(0,i4,i5),
     &                    Jijk(0,i1,i2,i3),Jij(0,i4,i5),m_top)
      endif
*----------------------------------------------------------------------
      if (HloopBottom) then
         aux(4)=cplModHbottom
     &     *H2G_Loop(Kijk(0,i1,i2,i3),Kij(0,i4,i5),
     &                    Jijk(0,i1,i2,i3),Jij(0,i4,i5),m_bottom)
      endif
*----------------------------------------------------------------------
      if (AloopTop) then
         aux(5)=cplModAtop
     &     *A2G_Loop(Kijk(0,i1,i2,i3),Kij(0,i4,i5),
     &                    Jijk(0,i1,i2,i3),Jij(0,i4,i5),m_top)
      endif
*----------------------------------------------------------------------
      if (AloopBottom) then
         aux(6)=cplModAbottom
     &     *A2G_Loop(Kijk(0,i1,i2,i3),Kij(0,i4,i5),
     &                    Jijk(0,i1,i2,i3),Jij(0,i4,i5),m_bottom)
      endif
*----------------------------------------------------------------------
      do i=1,26
         FD3gTri=FD3gTri+aux(i)
      enddo

      RETURN
      END

***********************************************************************
      complex*16 FUNCTION FD1ePV3g(i1,i2,i3,i4,i5,virtual)
*----------------------------------------------------------------------
c     3 three-gluon-vertices and effective polarization vector
***********************************************************************
      
      IMPLICIT NONE

#include "ggf_process.inc"
      
      integer i1,i2,i3,i4,i5
      logical virtual
      
      integer i,mu
      complex*16 aux(26)
      
      real*8 K(0:4,5),Kij(0:4,5,5),Kijk(0:4,5,5,5),PH(0:4),KPH(0:4,5),
     &     epsg(0:3,5),
     &     KPHeffij(0:4,5,5),KPAeffij(0:4,5,5),KPHtopij(0:4,5,5),
     &     KPHbottomij(0:4,5,5),KPAtopij(0:4,5,5),KPAbottompij(0:4,5,5),
     &     Kpen(0:4,4)
c... 4-dim momentum variables
      real*8 Kp(0:3,5),Kijp(0:3,5,5),Kijkp(0:3,5,5,5),PHp(0:3),
     &     KPHp(0:3,5),K123p(0:3),KPHeffijp(0:3,5,5),KPAeffijp(0:3,5,5),
     &     KPHtopijp(0:3,5,5),KPHbottomijp(0:3,5,5),KPAtopijp(0:3,5,5),
     &     KPAbottomijp(0:3,5,5)
      complex*16 epsgC(0:3,5),Jij(0:3,5,5),Jijk(0:3,5,5,5),
     &           epsgPen(0:3,4)
      complex*16 epsHeff(0:3,5),epsAeff(0:3,5),epsHtop(0:3,5),
     &     epsHbottom(0:3,5),epsAtop(0:3,5),epsAbottom(0:3,5),
     &     JPHeffij(0:3,5,5),JPAeffij(0:3,5,5),JPHtopij(0:3,5,5),
     &     JPHbottomij(0:3,5,5),JPAtopij(0:3,5,5),JPAbottomij(0:3,5,5)
      COMMON /GG3JET/ K,Kij,Kijk,epsg,epsgC,Jij,Jijk,KPH,PH,
     &     KPHeffij,KPAeffij,KPHtopij,KPHbottomij,KPAtopij,
     &     KPAbottompij,
     &     epsHeff,epsAeff,epsHtop,epsHbottom,epsAtop,epsAbottom,
     &     JPHeffij,JPAeffij,JPHtopij,JPHbottomij,JPAtopij,
     &     JPAbottomij,
c... 4-dim momentum variables
     &     Kp,Kijp,Kijkp,PHp,KPHp,KPHeffijp,KPAeffijp,K123p,
     &     KPHtopijp,KPHbottomijp,KPAtopijp,KPAbottomijp

      Common/Forpent/Kpen,epsgPen

      complex*16 GGG_vertex
      EXTERNAL GGG_vertex
      
      FD1ePV3g=DCMPLX(0d0,0d0)
      do i=1,26
         aux(i)=DCMPLX(0d0,0d0)
      enddo
      
      if (Heff) then
         if(virtual) then
            aux(1)=cplModHeff
     &            *GGG_vertex(Kijp(0,i1,i2),KPHp(0,i3),Kijp(0,i4,i5),
     &           Jij(0,i1,i2),epsHeff(0,i3),Jij(0,i4,i5))
         else
            aux(1)=cplModHeff
     &           *GGG_vertex(KPHeffijp(0,i1,i2),Kp(0,i3),Kijp(0,i4,i5),
     &           JPHeffij(0:3,i1,i2),epsgC(0,i3),Jij(0,i4,i5))
         endif
      endif
*----------------------------------------------------------------------
      if (Aeff) then
         if(virtual) then
            aux(2)=cplModAeff
     &            *GGG_vertex(Kijp(0,i1,i2),KPHp(0,i3),Kijp(0,i4,i5),
     &           Jij(0,i1,i2),epsAeff(0,i3),Jij(0,i4,i5))
         else
            aux(2)=cplModAeff
     &           *GGG_vertex(KPAeffijp(0,i1,i2),Kp(0,i3),Kijp(0,i4,i5),
     &           JPAeffij(0:3,i1,i2),epsgC(0,i3),Jij(0,i4,i5))
         endif
      endif
*----------------------------------------------------------------------
      if (HloopTop) then
         if(virtual) then
            aux(3)=cplModHtop
     &           *GGG_vertex(Kijp(0,i1,i2),KPHp(0,i3),Kijp(0,i4,i5),
     &           Jij(0,i1,i2),epsHtop(0,i3),Jij(0,i4,i5))
         else
            aux(3)=cplModHtop
     &           *GGG_vertex(KPHtopijp(0,i1,i2),Kp(0,i3),Kijp(0,i4,i5),
     &           JPHtopij(0,i1,i2),epsgC(0,i3),Jij(0,i4,i5))
         endif
      endif
*----------------------------------------------------------------------
      if (HloopBottom) then
         if(virtual) then
            aux(4)=cplModHbottom
     &           *GGG_vertex(Kijp(0,i1,i2),KPHp(0,i3),Kijp(0,i4,i5),
     &           Jij(0,i1,i2),epsHbottom(0,i3),Jij(0,i4,i5))
         else
            aux(4)=cplModHbottom
     &           *GGG_vertex(KPHbottomijp(0,i1,i2),Kp(0,i3),
     &           Kijp(0,i4,i5),JPHbottomij(0,i1,i2),epsgC(0,i3),
     &           Jij(0,i4,i5))
         endif
      endif
*----------------------------------------------------------------------
      if (AloopTop) then
         if(virtual) then
            aux(5)=cplModAtop
     &           *GGG_vertex(Kijp(0,i1,i2),KPHp(0,i3),Kijp(0,i4,i5),
     &           Jij(0,i1,i2),epsAtop(0,i3),Jij(0,i4,i5))
         else
            aux(5)=cplModAtop
     &           *GGG_vertex(KPAtopijp(0,i1,i2),Kp(0,i3),Kijp(0,i4,i5),
     &           JPAtopij(0,i1,i2),epsgC(0,i3),Jij(0,i4,i5))
         endif
      endif
*----------------------------------------------------------------------
      if (AloopBottom) then
         if(virtual) then
            aux(6)=cplModAbottom
     &           *GGG_vertex(Kijp(0,i1,i2),KPHp(0,i3),Kijp(0,i4,i5),
     &           Jij(0,i1,i2),epsAbottom(0,i3),Jij(0,i4,i5))
         else
            aux(6)=cplModAbottom
     &           *GGG_vertex(KPAbottomijp(0,i1,i2),Kp(0,i3),
     &           Kijp(0,i4,i5),JPAbottomij(0,i1,i2),epsgC(0,i3),
     &           Jij(0,i4,i5))
         endif
      endif
*----------------------------------------------------------------------
      do i=1,26
         FD1ePV3g=FD1ePV3g+aux(i)
      enddo

      RETURN
      END
 
***********************************************************************
      complex*16 FUNCTION FD2gBox(i1,i2,i3,i4,i5,virtual)
*----------------------------------------------------------------------
c     2 three-gluon-vertices and effective polarization vector
***********************************************************************
      
      IMPLICIT NONE

#include "ggf_process.inc"

      integer i1,i2,i3,i4,i5
      logical virtual

      integer i,mu
      complex*16 aux(26)
      
      real*8 K(0:4,5),Kij(0:4,5,5),Kijk(0:4,5,5,5),PH(0:4),KPH(0:4,5),
     &     epsg(0:3,5),
     &     KPHeffij(0:4,5,5),KPAeffij(0:4,5,5),KPHtopij(0:4,5,5),
     &     KPHbottomij(0:4,5,5),KPAtopij(0:4,5,5),KPAbottompij(0:4,5,5),
     &     Kpen(0:4,4)
c... 4-dim momentum variables
      real*8 Kp(0:3,5),Kijp(0:3,5,5),Kijkp(0:3,5,5,5),PHp(0:3),
     &     KPHp(0:3,5),K123p(0:3),KPHeffijp(0:3,5,5),KPAeffijp(0:3,5,5),
     &     KPHtopijp(0:3,5,5),KPHbottomijp(0:3,5,5),KPAtopijp(0:3,5,5),
     &     KPAbottomijp(0:3,5,5)
      complex*16 epsgC(0:3,5),Jij(0:3,5,5),Jijk(0:3,5,5,5),
     &           epsgPen(0:3,4)
      complex*16 epsHeff(0:3,5),epsAeff(0:3,5),epsHtop(0:3,5),
     &     epsHbottom(0:3,5),epsAtop(0:3,5),epsAbottom(0:3,5),
     &     JPHeffij(0:3,5,5),JPAeffij(0:3,5,5),JPHtopij(0:3,5,5),
     &     JPHbottomij(0:3,5,5),JPAtopij(0:3,5,5),JPAbottomij(0:3,5,5)
      COMMON /GG3JET/ K,Kij,Kijk,epsg,epsgC,Jij,Jijk,KPH,PH,
     &     KPHeffij,KPAeffij,KPHtopij,KPHbottomij,KPAtopij,
     &     KPAbottompij,
     &     epsHeff,epsAeff,epsHtop,epsHbottom,epsAtop,epsAbottom,
     &     JPHeffij,JPAeffij,JPHtopij,JPHbottomij,JPAtopij,
     &     JPAbottomij,
c... 4-dim momentum variables
     &     Kp,Kijp,Kijkp,PHp,KPHp,KPHeffijp,KPAeffijp,K123p,
     &     KPHtopijp,KPHbottomijp,KPAtopijp,KPAbottomijp

      Common/Forpent/Kpen,epsgPen
      
      complex*16 GGG_vertex,H3G_effective,A3G_effective,H3G_Loop,
     &     A3G_Loop
      EXTERNAL GGG_vertex,H3G_effective,A3G_effective,H3G_Loop,
     &     A3G_Loop

      FD2gBox=DCMPLX(0d0,0d0)
      do i=1,26
         aux(i)=DCMPLX(0d0,0d0)
      enddo

      if (Heff) then
      if(virtual) then
         aux(1)=cplModHeff*
     &       H3G_effective(Kp(0,i3),Kijp(0,i1,i2),Kijp(0,i4,i5),
     &        epsgC(0,i3),Jij(0,i1,i2),Jij(0,i4,i5))
      else
         aux(1)=cplModHeff*
     &       H3G_effective(Kp(0,i4),Kijkp(0,i1,i2,i3),Kp(0,i5),
     &        epsgC(0,i4),Jijk(0,i1,i2,i3),epsgC(0,i5))
      endif
      endif
*----------------------------------------------------------------------
      if (Aeff) then
      if(virtual) then
         aux(2)=cplModAeff*
     &       A3G_effective(PHp,epsgC(0,i3),Jij(0,i1,i2),Jij(0,i4,i5))
      else
         aux(2)=cplModAeff*
     &       A3G_effective(PHp,epsgC(0,i4),Jijk(0,i1,i2,i3),epsgC(0,i5))
      endif
      endif
*----------------------------------------------------------------------
      if (HloopTop) then
      if(virtual) then
         aux(3)=cplModHtop*
     &       H3G_Loop(K(0,i3),Kij(0,i1,i2),Kij(0,i4,i5),
     &        epsgC(0,i3),Jij(0,i1,i2),Jij(0,i4,i5),m_top)
      else
         aux(3)=cplModHtop*
     &       H3G_Loop(K(0,i4),Kijk(0,i1,i2,i3),K(0,i5),
     &        epsgC(0,i4),Jijk(0,i1,i2,i3),epsgC(0,i5),m_top)
      endif
      endif
*----------------------------------------------------------------------
      if (HloopBottom) then
      if(virtual) then
         aux(4)=cplModHbottom*
     &       H3G_Loop(K(0,i3),Kij(0,i1,i2),Kij(0,i4,i5),
     &        epsgC(0,i3),Jij(0,i1,i2),Jij(0,i4,i5),m_bottom)
      else
         aux(4)=cplModHbottom*
     &       H3G_Loop(K(0,i4),Kijk(0,i1,i2,i3),K(0,i5),
     &        epsgC(0,i4),Jijk(0,i1,i2,i3),epsgC(0,i5),m_bottom)
      endif
      endif
*----------------------------------------------------------------------
      if (AloopTop) then
      if(virtual) then
c         print*,'------- in FD2gBox virtual -------'
         aux(5)=cplModAtop*
     &       A3G_Loop(K(0,i3),Kij(0,i1,i2),Kij(0,i4,i5),
     &        epsgC(0,i3),Jij(0,i1,i2),Jij(0,i4,i5),m_top)
      else
c         print*,'------- in FD2gBox  -------'
         aux(5)=cplModAtop*
     &       A3G_Loop(K(0,i4),Kijk(0,i1,i2,i3),K(0,i5),
     &        epsgC(0,i4),Jijk(0,i1,i2,i3),epsgC(0,i5),m_top)
      endif
      endif
*----------------------------------------------------------------------
      if (AloopBottom) then
      if(virtual) then
         aux(6)=cplModAbottom*
     &       A3G_Loop(K(0,i3),Kij(0,i1,i2),Kij(0,i4,i5),
     &        epsgC(0,i3),Jij(0,i1,i2),Jij(0,i4,i5),m_bottom)
      else
         aux(6)=cplModAbottom*
     &       A3G_Loop(K(0,i4),Kijk(0,i1,i2,i3),K(0,i5),
     &        epsgC(0,i4),Jijk(0,i1,i2,i3),epsgC(0,i5),m_bottom)
      endif
      endif
*----------------------------------------------------------------------

      do i=1,26
         FD2gBox=FD2gBox+aux(i)
      enddo
      
      RETURN
      END

***********************************************************************
      complex*16 FUNCTION FDePV4g(i1,i2,i3,i4,i5,gp)
*----------------------------------------------------------------------

***********************************************************************
      
      IMPLICIT NONE

#include "ggf_process.inc"

      integer i1,i2,i3,i4,i5,gp
      logical ePV

      integer i,mu
      complex*16 Vert4G(0:3),aux(26)
      real*8 K123(0:4)

      real*8 K(0:4,5),Kij(0:4,5,5),Kijk(0:4,5,5,5),PH(0:4),KPH(0:4,5),
     &     epsg(0:3,5),
     &     KPHeffij(0:4,5,5),KPAeffij(0:4,5,5),KPHtopij(0:4,5,5),
     &     KPHbottomij(0:4,5,5),KPAtopij(0:4,5,5),KPAbottompij(0:4,5,5),
     &     Kpen(0:4,4)
c... 4-dim momentum variables
      real*8 Kp(0:3,5),Kijp(0:3,5,5),Kijkp(0:3,5,5,5),PHp(0:3),
     &     KPHp(0:3,5),K123p(0:3),KPHeffijp(0:3,5,5),KPAeffijp(0:3,5,5),
     &     KPHtopijp(0:3,5,5),KPHbottomijp(0:3,5,5),KPAtopijp(0:3,5,5),
     &     KPAbottomijp(0:3,5,5)
      complex*16 epsgC(0:3,5),Jij(0:3,5,5),Jijk(0:3,5,5,5),
     &           epsgPen(0:3,4)
      complex*16 epsHeff(0:3,5),epsAeff(0:3,5),epsHtop(0:3,5),
     &     epsHbottom(0:3,5),epsAtop(0:3,5),epsAbottom(0:3,5),
     &     JPHeffij(0:3,5,5),JPAeffij(0:3,5,5),JPHtopij(0:3,5,5),
     &     JPHbottomij(0:3,5,5),JPAtopij(0:3,5,5),JPAbottomij(0:3,5,5)
      COMMON /GG3JET/ K,Kij,Kijk,epsg,epsgC,Jij,Jijk,KPH,PH,
     &     KPHeffij,KPAeffij,KPHtopij,KPHbottomij,KPAtopij,
     &     KPAbottompij,
     &     epsHeff,epsAeff,epsHtop,epsHbottom,epsAtop,epsAbottom,
     &     JPHeffij,JPAeffij,JPHtopij,JPHbottomij,JPAtopij,
     &     JPAbottomij,
c... 4-dim momentum variables
     &     Kp,Kijp,Kijkp,PHp,KPHp,KPHeffijp,KPAeffijp,K123p,
     &     KPHtopijp,KPHbottomijp,KPAtopijp,KPAbottomijp

      Common/Forpent/Kpen,epsgPen

      real*8 dotrr
      complex*16 H4G_effective2,H2G_effective,A2G_effective,H2G_Loop,
     &     A2G_Loop
      EXTERNAL dotrr,H4G_effective2,H2G_effective,A2G_effective,
     &     H2G_Loop,A2G_Loop

      FDePV4g=DCMPLX(0d0,0d0)
      do i=1,26
         aux(i)=DCMPLX(0d0,0d0)
      enddo

      if (Heff) then
      if (gp.eq.0) then
         aux(1)=cplModHeff*H4G_effective2(epsgC(0,i1),epsgC(0,i2),
     &                       epsHeff(0,i3),Jij(0,i4,i5))
      endif
      if (gp.eq.1) then
         aux(1)=cplModHeff*H4G_effective2(epsgC(0,i1),epsHeff(0,i2),
     &                     epsgC(0,i3),Jij(0,i4,i5))
      endif
      if (gp.eq.2) then
         aux(1)=cplModHeff*H4G_effective2(JPHeffij(0,i1,i2),epsgC(0,i3),
     &                     epsgC(0,i4),epsgC(0,i5))
      endif
      if (gp.eq.3) then
         CALL Vertex4g(epsgC(0,i1),epsgC(0,i2),epsgC(0,i3),Vert4G)
         do mu=0,3
            K123p(mu)=(Kp(mu,i1)+Kp(mu,i2)+Kp(mu,i3))
            K123(mu)=K123p(mu)
         enddo
         K123(4)=dotrr(K123p,K123p)
         aux(1)=cplModHeff*
     &    H2G_effective(K123p,Kijp(0,i4,i5),Vert4G,Jij(0,i4,i5))/K123(4)
      endif
      endif
*----------------------------------------------------------------------
      if (Aeff) then
      if (gp.eq.0) then
         aux(2)=cplModAeff*H4G_effective2(epsgC(0,i1),epsgC(0,i2),
     &                       epsAeff(0,i3),Jij(0,i4,i5))
      endif
      if (gp.eq.1) then
         aux(2)=cplModAeff*H4G_effective2(epsgC(0,i1),epsAeff(0,i2),
     &                     epsgC(0,i3),Jij(0,i4,i5))
      endif
      if (gp.eq.2) then
         aux(2)=cplModAeff*H4G_effective2(JPAeffij(0,i1,i2),epsgC(0,i3),
     &                     epsgC(0,i4),epsgC(0,i5))
      endif
      if (gp.eq.3) then
         CALL Vertex4g(epsgC(0,i1),epsgC(0,i2),epsgC(0,i3),Vert4G)
         do mu=0,3
            K123p(mu)=(Kp(mu,i1)+Kp(mu,i2)+Kp(mu,i3))
            K123(mu)=K123p(mu)
         enddo
         K123(4)=dotrr(K123p,K123p)
         aux(2)=cplModAeff*
     &    A2G_effective(K123p,Kijp(0,i4,i5),Vert4G,Jij(0,i4,i5))/K123(4)
      endif
      endif
*----------------------------------------------------------------------
      if (HloopTop) then
      if (gp.eq.0) then
         aux(3)=cplModHtop*H4G_effective2(epsgC(0,i1),epsgC(0,i2),
     &                       epsHtop(0,i3),Jij(0,i4,i5))
      endif
      if (gp.eq.1) then
         aux(3)=cplModHtop*H4G_effective2(epsgC(0,i1),epsHtop(0,i2),
     &                     epsgC(0,i3),Jij(0,i4,i5))
      endif
      if (gp.eq.2) then
         aux(3)=cplModHtop*H4G_effective2(JPHtopij(0,i1,i2),epsgC(0,i3),
     &                     epsgC(0,i4),epsgC(0,i5))
      endif
      if (gp.eq.3) then
         CALL Vertex4g(epsgC(0,i1),epsgC(0,i2),epsgC(0,i3),Vert4G)
         do mu=0,3
            K123p(mu)=(Kp(mu,i1)+Kp(mu,i2)+Kp(mu,i3))
            K123(mu)=K123p(mu)
         enddo
         K123(4)=dotrr(K123p,K123p) 
         aux(3)=cplModHtop*
     &     H2G_Loop(K123,Kij(0,i4,i5),Vert4G,Jij(0,i4,i5),m_top)/K123(4)
      endif
      endif
*----------------------------------------------------------------------
      if (HloopBottom) then
      if (gp.eq.0) then
         aux(4)=cplModHbottom*H4G_effective2(epsgC(0,i1),epsgC(0,i2),
     &                       epsHbottom(0,i3),Jij(0,i4,i5))
      endif
      if (gp.eq.1) then
         aux(4)=cplModHbottom*H4G_effective2(epsgC(0,i1),
     &        epsHbottom(0,i2),epsgC(0,i3),Jij(0,i4,i5))
      endif
      if (gp.eq.2) then
         aux(4)=cplModHbottom*H4G_effective2(JPHbottomij(0,i1,i2),
     &        epsgC(0,i3),epsgC(0,i4),epsgC(0,i5))
      endif
      if (gp.eq.3) then
         CALL Vertex4g(epsgC(0,i1),epsgC(0,i2),epsgC(0,i3),Vert4G)
         do mu=0,3
            K123p(mu)=(Kp(mu,i1)+Kp(mu,i2)+Kp(mu,i3))
            K123(mu)=K123p(mu)
         enddo
         K123(4)=dotrr(K123p,K123p) 
         aux(4)=cplModHbottom*
     &     H2G_Loop(K123,Kij(0,i4,i5),Vert4G,Jij(0,i4,i5),m_bottom)/
     &        K123(4)
      endif
      endif
*----------------------------------------------------------------------
      if (AloopTop) then
      if (gp.eq.0) then
         aux(5)=cplModAtop*H4G_effective2(epsgC(0,i1),epsgC(0,i2),
     &                       epsAtop(0,i3),Jij(0,i4,i5))
      endif
      if (gp.eq.1) then
         aux(5)=cplModAtop*H4G_effective2(epsgC(0,i1),epsAtop(0,i2),
     &                     epsgC(0,i3),Jij(0,i4,i5))
      endif
      if (gp.eq.2) then
         aux(5)=cplModAtop*H4G_effective2(JPAtopij(0,i1,i2),epsgC(0,i3),
     &                     epsgC(0,i4),epsgC(0,i5))
      endif
      if (gp.eq.3) then
         CALL Vertex4g(epsgC(0,i1),epsgC(0,i2),epsgC(0,i3),Vert4G)
         do mu=0,3
            K123p(mu)=(Kp(mu,i1)+Kp(mu,i2)+Kp(mu,i3))
            K123(mu)=K123p(mu)
         enddo
         K123(4)=dotrr(K123p,K123p) 
         aux(5)=cplModAtop*
     &     A2G_Loop(K123,Kij(0,i4,i5),Vert4G,Jij(0,i4,i5),m_top)/K123(4)
      endif
      endif
*----------------------------------------------------------------------
      if (AloopBottom) then
      if (gp.eq.0) then
         aux(6)=cplModAbottom*H4G_effective2(epsgC(0,i1),epsgC(0,i2),
     &                       epsAbottom(0,i3),Jij(0,i4,i5))
      endif
      if (gp.eq.1) then
         aux(6)=cplModAbottom*H4G_effective2(epsgC(0,i1),
     &        epsAbottom(0,i2),epsgC(0,i3),Jij(0,i4,i5))
      endif
      if (gp.eq.2) then
         aux(6)=cplModAbottom*H4G_effective2(JPAbottomij(0,i1,i2),
     &        epsgC(0,i3),epsgC(0,i4),epsgC(0,i5))
      endif
      if (gp.eq.3) then
         CALL Vertex4g(epsgC(0,i1),epsgC(0,i2),epsgC(0,i3),Vert4G)
         do mu=0,3
            K123p(mu)=(Kp(mu,i1)+Kp(mu,i2)+Kp(mu,i3))
            K123(mu)=K123p(mu)
         enddo
         K123(4)=dotrr(K123p,K123p) 
         aux(6)=cplModAbottom*
     &     A2G_Loop(K123,Kij(0,i4,i5),Vert4G,Jij(0,i4,i5),m_bottom)/
     &        K123(4)
      endif
      endif
*----------------------------------------------------------------------
      do i=1,26
         FDePV4g=FDePV4g+aux(i)
      enddo
      
      RETURN
      END


***********************************************************************
      complex*16 FUNCTION FD4gBox(i1,i2,i3,i4,i5)
*----------------------------------------------------------------------

***********************************************************************
      
      IMPLICIT NONE

#include "ggf_process.inc"

      integer i1,i2,i3,i4,i5
      logical box

      integer i,mu
      complex*16 Vert4G(0:3),aux(26)
      real*8 K123(0:4)

      real*8 K(0:4,5),Kij(0:4,5,5),Kijk(0:4,5,5,5),PH(0:4),KPH(0:4,5),
     &     epsg(0:3,5),
     &     KPHeffij(0:4,5,5),KPAeffij(0:4,5,5),KPHtopij(0:4,5,5),
     &     KPHbottomij(0:4,5,5),KPAtopij(0:4,5,5),KPAbottompij(0:4,5,5),
     &     Kpen(0:4,4)
c... 4-dim momentum variables
      real*8 Kp(0:3,5),Kijp(0:3,5,5),Kijkp(0:3,5,5,5),PHp(0:3),
     &     KPHp(0:3,5),K123p(0:3),KPHeffijp(0:3,5,5),KPAeffijp(0:3,5,5),
     &     KPHtopijp(0:3,5,5),KPHbottomijp(0:3,5,5),KPAtopijp(0:3,5,5),
     &     KPAbottomijp(0:3,5,5)
      complex*16 epsgC(0:3,5),Jij(0:3,5,5),Jijk(0:3,5,5,5),
     &           epsgPen(0:3,4)
      complex*16 epsHeff(0:3,5),epsAeff(0:3,5),epsHtop(0:3,5),
     &     epsHbottom(0:3,5),epsAtop(0:3,5),epsAbottom(0:3,5),
     &     JPHeffij(0:3,5,5),JPAeffij(0:3,5,5),JPHtopij(0:3,5,5),
     &     JPHbottomij(0:3,5,5),JPAtopij(0:3,5,5),JPAbottomij(0:3,5,5)
      COMMON /GG3JET/ K,Kij,Kijk,epsg,epsgC,Jij,Jijk,KPH,PH,
     &     KPHeffij,KPAeffij,KPHtopij,KPHbottomij,KPAtopij,
     &     KPAbottompij,
     &     epsHeff,epsAeff,epsHtop,epsHbottom,epsAtop,epsAbottom,
     &     JPHeffij,JPAeffij,JPHtopij,JPHbottomij,JPAtopij,
     &     JPAbottomij,
c... 4-dim momentum variables
     &     Kp,Kijp,Kijkp,PHp,KPHp,KPHeffijp,KPAeffijp,K123p,
     &     KPHtopijp,KPHbottomijp,KPAtopijp,KPAbottomijp

      Common/Forpent/Kpen,epsgPen

      real*8 dotrr
      complex*16 H3G_effective,A3G_effective,H3G_Loop,A3G_Loop
      EXTERNAL dotrr,H3G_effective,A3G_effective,H3G_Loop,A3G_Loop

      FD4gBox=DCMPLX(0d0,0d0)
      do i=1,26
         aux(i)=DCMPLX(0d0,0d0)
      enddo

      if (Heff) then
         CALL Vertex4g(epsgC(0,i1),epsgC(0,i2),epsgC(0,i3),Vert4G)
         do mu=0,3
            K123p(mu)=(Kp(mu,i1)+Kp(mu,i2)+Kp(mu,i3))
            K123(mu)=K123p(mu)
         enddo
         K123(4)=dotrr(K123p,K123p)
         aux(1)=cplModHeff*
     &       H3G_effective(Kp(0,i4),Kp(0,i5),K123p,
     &                     epsgC(0,i4),epsgC(0,i5),Vert4G)/K123(4)
      endif
*----------------------------------------------------------------------
      if (Aeff) then
         CALL Vertex4g(epsgC(0,i1),epsgC(0,i2),epsgC(0,i3),Vert4G)
         do mu=0,3
            K123p(mu)=(Kp(mu,i1)+Kp(mu,i2)+Kp(mu,i3))
            K123(mu)=K123p(mu)
         enddo
         K123(4)=dotrr(K123p,K123p)
         aux(2)=cplModAeff*
     &       A3G_effective(PHp,epsgC(0,i4),epsgC(0,i5),Vert4G)/K123(4)
      endif
*----------------------------------------------------------------------
      if (HloopTop) then
         CALL Vertex4g(epsgC(0,i1),epsgC(0,i2),epsgC(0,i3),Vert4G)
         do mu=0,3
            K123p(mu)=(Kp(mu,i1)+Kp(mu,i2)+Kp(mu,i3))
            K123(mu)=K123p(mu)
         enddo
         K123(4)=dotrr(K123p,K123p)
         aux(3)=cplModHtop*
     &       H3G_Loop(K(0,i4),K(0,i5),K123,epsgC(0,i4),epsgC(0,i5),
     &                Vert4G,m_top)/K123(4)
      endif
*----------------------------------------------------------------------
      if (HloopBottom) then
         CALL Vertex4g(epsgC(0,i1),epsgC(0,i2),epsgC(0,i3),Vert4G)
         do mu=0,3
            K123p(mu)=(Kp(mu,i1)+Kp(mu,i2)+Kp(mu,i3))
            K123(mu)=K123p(mu)
         enddo
         K123(4)=dotrr(K123p,K123p)
         aux(4)=cplModHbottom*
     &       H3G_Loop(K(0,i4),K(0,i5),K123,epsgC(0,i4),epsgC(0,i5),
     &                Vert4G,m_bottom)/K123(4)
      endif
*----------------------------------------------------------------------
      if (AloopTop) then
c         print*,'------- in FD4gBox  -------'
         CALL Vertex4g(epsgC(0,i1),epsgC(0,i2),epsgC(0,i3),Vert4G)
         do mu=0,3
            K123p(mu)=(Kp(mu,i1)+Kp(mu,i2)+Kp(mu,i3))
            K123(mu)=K123p(mu)
         enddo
         K123(4)=dotrr(K123p,K123p)
         aux(5)=cplModAtop*
     &       A3G_Loop(K(0,i4),K(0,i5),K123,epsgC(0,i4),epsgC(0,i5),
     &                Vert4G,m_top)/K123(4)
      endif
*----------------------------------------------------------------------
      if (AloopBottom) then
         CALL Vertex4g(epsgC(0,i1),epsgC(0,i2),epsgC(0,i3),Vert4G)
         do mu=0,3
            K123p(mu)=(Kp(mu,i1)+Kp(mu,i2)+Kp(mu,i3))
            K123(mu)=K123p(mu)
         enddo
         K123(4)=dotrr(K123p,K123p)
         aux(6)=cplModAbottom*
     &       A3G_Loop(K(0,i4),K(0,i5),K123,epsgC(0,i4),epsgC(0,i5),
     &                Vert4G,m_bottom)/K123(4)
      endif
*----------------------------------------------------------------------

      do i=1,26
         FD4gBox=FD4gBox+aux(i)
      enddo
      
      RETURN
      END

***********************************************************************
      complex*16 FUNCTION FDeffPen(i1,i2,i3,i4,i5)
*----------------------------------------------------------------------

***********************************************************************
      
      IMPLICIT NONE

#include "ggf_process.inc"

      integer i1,i2,i3,i4,i5

      integer i,mu
      complex*16 Vert4G(0:3),aux(26)
      real*8 K123(0:4)
     
      real*8 K(0:4,5),Kij(0:4,5,5),Kijk(0:4,5,5,5),PH(0:4),KPH(0:4,5),
     &     epsg(0:3,5),
     &     KPHeffij(0:4,5,5),KPAeffij(0:4,5,5),KPHtopij(0:4,5,5),
     &     KPHbottomij(0:4,5,5),KPAtopij(0:4,5,5),KPAbottompij(0:4,5,5),
     &     Kpen(0:4,4)
c... 4-dim momentum variables
      real*8 Kp(0:3,5),Kijp(0:3,5,5),Kijkp(0:3,5,5,5),PHp(0:3),
     &     KPHp(0:3,5),K123p(0:3),KPHeffijp(0:3,5,5),KPAeffijp(0:3,5,5),
     &     KPHtopijp(0:3,5,5),KPHbottomijp(0:3,5,5),KPAtopijp(0:3,5,5),
     &     KPAbottomijp(0:3,5,5)
      complex*16 epsgC(0:3,5),Jij(0:3,5,5),Jijk(0:3,5,5,5),
     &           epsgPen(0:3,4)
      complex*16 epsHeff(0:3,5),epsAeff(0:3,5),epsHtop(0:3,5),
     &     epsHbottom(0:3,5),epsAtop(0:3,5),epsAbottom(0:3,5),
     &     JPHeffij(0:3,5,5),JPAeffij(0:3,5,5),JPHtopij(0:3,5,5),
     &     JPHbottomij(0:3,5,5),JPAtopij(0:3,5,5),JPAbottomij(0:3,5,5)
      COMMON /GG3JET/ K,Kij,Kijk,epsg,epsgC,Jij,Jijk,KPH,PH,
     &     KPHeffij,KPAeffij,KPHtopij,KPHbottomij,KPAtopij,
     &     KPAbottompij,
     &     epsHeff,epsAeff,epsHtop,epsHbottom,epsAtop,epsAbottom,
     &     JPHeffij,JPAeffij,JPHtopij,JPHbottomij,JPAtopij,
     &     JPAbottomij,
c... 4-dim momentum variables
     &     Kp,Kijp,Kijkp,PHp,KPHp,KPHeffijp,KPAeffijp,K123p,
     &     KPHtopijp,KPHbottomijp,KPAtopijp,KPAbottomijp

      Common/Forpent/Kpen,epsgPen
      
      complex*16 H4G_effective2
      EXTERNAL H4G_effective2
      logical MGdebug1,MGPen,MGHex,PRINT_AMP
       common/MGdebugC1/MGdebug1,MGPen,MGHex,PRINT_AMP

      FDeffPen=DCMPLX(0d0,0d0)
      do i=1,26
         aux(i)=DCMPLX(0d0,0d0)
      enddo

      if (Heff) then
         aux(1)=cplModHeff*H4G_effective2(Jij(0,i1,i2),epsgC(0,i3),
     &        epsgC(0,i4),epsgC(0,i5))
c         print*, "here",aux(1)
      endif
*----------------------------------------------------------------------
      if (Aeff) then
         aux(2)=DCMPLX(0d0,0d0)
      endif
*----------------------------------------------------------------------
      if (HloopTop.or.HloopBottom.or.AloopTop.or.AloopBottom) then
         aux(3)=DCMPLX(0d0,0d0)
      endif
*----------------------------------------------------------------------

      if(MGpen) then
        do i=1,26
         aux(i)=DCMPLX(0d0,0d0)
       enddo
      endif 

      do i=1,26
         FDeffPen=FDeffPen+aux(i)
      enddo
c      print*,"here",FDeffPen
      RETURN
      END

***********************************************************************
      complex*16 FUNCTION FDeffPenLoop(i1,i2,i3,i4)
*----------------------------------------------------------------------

***********************************************************************
      
      IMPLICIT NONE

#include "ggf_process.inc"

      integer i1,i2,i3,i4,i5

      integer i,mu
      complex*16 Vert4G(0:3),aux(26),qpentH_top(4)
      real*8 K123(0:4)
     
      real*8 K(0:4,5),Kij(0:4,5,5),Kijk(0:4,5,5,5),PH(0:4),KPH(0:4,5),
     &     epsg(0:3,5),
     &     KPHeffij(0:4,5,5),KPAeffij(0:4,5,5),KPHtopij(0:4,5,5),
     &     KPHbottomij(0:4,5,5),KPAtopij(0:4,5,5),KPAbottompij(0:4,5,5),
     &     Kpen(0:4,4)
c... 4-dim momentum variables
      real*8 Kp(0:3,5),Kijp(0:3,5,5),Kijkp(0:3,5,5,5),PHp(0:3),
     &     KPHp(0:3,5),K123p(0:3),KPHeffijp(0:3,5,5),KPAeffijp(0:3,5,5),
     &     KPHtopijp(0:3,5,5),KPHbottomijp(0:3,5,5),KPAtopijp(0:3,5,5),
     &     KPAbottomijp(0:3,5,5)
      complex*16 epsgC(0:3,5),Jij(0:3,5,5),Jijk(0:3,5,5,5),
     &           epsgPen(0:3,4)
      complex*16 epsHeff(0:3,5),epsAeff(0:3,5),epsHtop(0:3,5),
     &     epsHbottom(0:3,5),epsAtop(0:3,5),epsAbottom(0:3,5),
     &     JPHeffij(0:3,5,5),JPAeffij(0:3,5,5),JPHtopij(0:3,5,5),
     &     JPHbottomij(0:3,5,5),JPAtopij(0:3,5,5),JPAbottomij(0:3,5,5)
      COMMON /GG3JET/ K,Kij,Kijk,epsg,epsgC,Jij,Jijk,KPH,PH,
     &     KPHeffij,KPAeffij,KPHtopij,KPHbottomij,KPAtopij,
     &     KPAbottompij,
     &     epsHeff,epsAeff,epsHtop,epsHbottom,epsAtop,epsAbottom,
     &     JPHeffij,JPAeffij,JPHtopij,JPHbottomij,JPAtopij,
     &     JPAbottomij,
c... 4-dim momentum variables
     &     Kp,Kijp,Kijkp,PHp,KPHp,KPHeffijp,KPAeffijp,K123p,
     &     KPHtopijp,KPHbottomijp,KPAtopijp,KPAbottomijp

      Common/Forpent/Kpen,epsgPen
    
      complex*16 H4G_effective2
      EXTERNAL H4G_effective2
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      complex*16 mup1temp(0:3),mup2temp(0:3),mup3temp(0:3),mup4temp(0:3)
      complex*16 pen_Ward1,pen_Ward2,pen_Ward3,pen_Ward4
      real*8 pp12(0:3),pp23(0:3),pp34(0:3),pp45(0:3),pp15(0:3)
      real*8 gauge1,gauge2,gauge3,gauge4
      real*8 gauge
      complex*16 box1_for_Ward1, box1_for_Ward2
      complex*16 box2_for_ward2, box2_for_ward3
      complex*16 box3_for_ward3, box3_for_ward4
      complex*16 box4_for_ward4, box4_for_ward1


      logical QUAD,Dble1
      real*8 Kpen1(0:3,4),PH1(0:3)
      complex*16 pen1
      integer jtemp
       real*8 gaugelimit 
       common/gauge/gaugelimit
      common/dbleswitch/dble1,quad
c
c     
      logical eva
      real*8 phc(0:3)
      common/phc1/phc,eva
      
      integer initdo,finaldo,jmass
      real*8 massj(8),cplmassj(8)
      integer shift

      logical SelectABottom,SelectATop,SelectHBottom,SelectHTop
      common/select/SelectABottom,SelectATop,SelectHBottom,SelectHTop

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
c       penSingular(1)=.true.
c       penSingular(2)=.true.


      do jtemp=1,4
      do i=0,3
         Kpen1(i,jtemp)=Kpen(i,jtemp)
c         PH1(i)=PH(i)
      enddo
      enddo

      do i=0,3
        PH1(i)=PH(i)
      enddo

c      print*,"eva",eva
      if(eva) then
      do i=0,3
         PH1(i)=PHc(i)
c         print*,"PH1",PH1(i)
      enddo
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FDeffPenLoop=DCMPLX(0d0,0d0)
      do i=1,26
         aux(i)=DCMPLX(0d0,0d0)
      enddo

      if(HloopTop.or.HloopBottom) then

      if (HloopTop) then
       initdo=1
       else
       initdo=2
      endif
      if(HloopBottom) then
         finaldo=2
      else
        finaldo=1
      endif


      If(SelectHbottom) then
         initdo=2
         finaldo=2
      endif
      If(SelectHtop) then
         initdo=1
         finaldo=1
      endif

      If(SelectATop.or.SelectABottom) then
         go to 100
      endif
         
      
      massj(1)=m_top
      massj(2)=m_bottom
      cplmassj(1)=cplModHtop
      cplmassj(2)=cplModHbottom
      

      do jmass=initdo,finaldo


        If(dble1) then

c         CALL H4G_LoopDD(i4,i3,i2,i1,m_top,qpentH_top(1))
c         CALL H4G_LoopDD(i2,i3,i4,i1,m_top,qpentH_top(2))
c         CALL H4G_LoopDD(i2,i1,i4,i3,m_top,qpentH_top(3))
c         CALL H4G_LoopDD(i3,i2,i1,i4,m_top,qpentH_top(4))
c         aux(1)=cplModHtop
c     &    *(qpentH_top(1)+qpentH_top(2)+qpentH_top(3)+qpentH_top(4))
c         CALL H4G_LoopDD(i1,i2,i3,i4,m_top,qpentH_top(1))
c         aux(1)=cplModHtop*qpentH_top(1)
cc start check
c         aux(1)=-cplModHtop*qpentH_top(1)
c
         CALL gluonpenHiggsEven(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),
     &        Kpen1(0,i3),Kpen1(0,i4),PH1,epsgPen(0,i1),epsgPen(0,i2),
     &        epsgPen(0,i3),epsgPen(0,i4),pen1,1)
         pen1=-cplmassj(jmass)*pen1/(12d0*massj(jmass)*massj(jmass))
c         print*, "pen1",abs(pen1/aux(jmass)-1d0)
         aux(jmass)= pen1
c start Ward identity
c The momenta of the Hexline is going to be the same. 
c                       --->  comp=-1 to avoid redundant calculation
c First Gauge test Eps_i -->   p_i
c +  auxiliary momenta
c         print*, "i1--i5", i1,i2,i3,i4,i5
         do jtemp=0,3
c            print*, "epsgC(0,i1)",epsgC(jtemp,i1)
            mup1temp(jtemp)= Kpen1(jtemp,i1)
            mup2temp(jtemp)= Kpen1(jtemp,i2)
            mup3temp(jtemp)= Kpen1(jtemp,i3)
            mup4temp(jtemp)= Kpen1(jtemp,i4)
            pp12(jtemp)=Kpen1(jtemp,i1)+Kpen1(jtemp,i2)
            pp23(jtemp)=Kpen1(jtemp,i2)+Kpen1(jtemp,i3)
            pp34(jtemp)=Kpen1(jtemp,i3)+Kpen1(jtemp,i4)
            pp45(jtemp)=-(Kpen1(jtemp,i1)+Kpen1(jtemp,i2)+Kpen1(jtemp,i3))
            pp15(jtemp)=-(Kpen1(jtemp,i2)+Kpen1(jtemp,i3)+Kpen1(jtemp,i4))
         enddo   

c Ward 1
         CALL gluonpenHiggsEven(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),
     &        Kpen1(0,i3),Kpen1(0,i4),PH1,mup1temp,epsgPen(0,i2),
     &        epsgPen(0,i3),epsgPen(0,i4),pen_Ward1,-1)
c Ward 2
         CALL gluonpenHiggsEven(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),
     &        Kpen1(0,i3),Kpen1(0,i4),PH1,epsgPen(0,i1),mup2temp,
     &        epsgPen(0,i3),epsgPen(0,i4),pen_Ward2,-1)
c Ward 3
         CALL gluonpenHiggsEven(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),
     &        Kpen1(0,i3),Kpen1(0,i4),PH1,epsgPen(0,i1),
     &        epsgPen(0,i2),mup3temp,epsgPen(0,i4),pen_Ward3,-1)
c Ward 4
         CALL gluonpenHiggsEven(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),
     &        Kpen1(0,i3),Kpen1(0,i4),PH1,epsgPen(0,i1),
     &        epsgPen(0,i2),epsgPen(0,i3),mup4temp,pen_Ward4,-1)
c         print*, "aux(jmass)",aux(jmass)
c         print*, "pen1Ward1",pen_Ward1
c         print*, "pen1Ward2",pen_Ward2
c         print*, "pen1Ward3",pen_Ward3
c         print*, "pen1Ward4",pen_Ward4
c         stop
c Call of pentagons
c  Each pentagon will appear twice. eg., pen1:  pen2 -pen1,  pen1 -pen0
c wich different sign depending for which contraction is used.
c      print*, "epsgC(0,i5)_2",epsgC(0,i5),epsgC(1,i5),epsgC(2,i5),epsgC(3,i5)
c  pp12
       CALL gluonboxHiggsEven(massj(jmass),pp12,Kpen1(0,i3),Kpen1(0,i4),PH1,
     &  epsgPen(0,i2),epsgPen(0,i3),epsgPen(0,i4),box1_for_Ward1,1)
       CALL gluonboxHiggsEven(massj(jmass),pp12,Kpen1(0,i3),Kpen1(0,i4),PH1,
     &  epsgPen(0,i1),epsgPen(0,i3),epsgPen(0,i4),box1_for_Ward2,-1)
c  pp23
       CALL gluonboxHiggsEven(massj(jmass),Kpen1(0,i1),pp23,Kpen1(0,i4),PH1,
     &  epsgPen(0,i1),epsgPen(0,i3),epsgPen(0,i4),box2_for_Ward2,1)
       CALL gluonboxHiggsEven(massj(jmass),Kpen1(0,i1),pp23,Kpen1(0,i4),PH1,
     &  epsgPen(0,i1),epsgPen(0,i2),epsgPen(0,i4),box2_for_Ward3,-1)
c  pp34
       CALL gluonboxHiggsEven(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),pp34,PH1,
     &  epsgPen(0,i1),epsgPen(0,i2),epsgPen(0,i4),box3_for_Ward3,1)
       CALL gluonboxHiggsEven(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),pp34,PH1,
     &  epsgPen(0,i1),epsgPen(0,i2),epsgPen(0,i3),box3_for_Ward4,-1)
c  pp45
       CALL gluonboxHiggsEven(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),Kpen1(0,i3),pp45,
     &  epsgPen(0,i1),epsgPen(0,i2),epsgPen(0,i3),box4_for_Ward4,1)
c pp15
       CALL gluonboxHiggsEven(massj(jmass),Kpen1(0,i2),Kpen1(0,i3),Kpen1(0,i4),pp15,
     &  epsgPen(0,i2),epsgPen(0,i3),epsgPen(0,i4),box4_for_Ward1,1)

c gauge1
         if(abs(Pen_Ward1).gt.1d-2) then
         gauge1=abs((box1_for_Ward1-box4_for_Ward1)/(Pen_Ward1)-1d0)
         else
         gauge1=abs((box1_for_Ward1-box4_for_Ward1) - (Pen_Ward1))
         endif
c gauge2
         if(abs(Pen_Ward2).gt.1d-2) then
         gauge2=abs((box2_for_Ward2-box1_for_Ward2)/(Pen_Ward2)-1d0)
         else
         gauge2=abs((box2_for_Ward2-box1_for_Ward2) - (Pen_Ward2))
         endif
c gauge3
         if(abs(Pen_Ward3).gt.1d-2) then
         gauge3=abs((box3_for_Ward3-box2_for_Ward3)/(Pen_Ward3)-1d0)
         else
         gauge3=abs((box3_for_Ward3-box2_for_Ward3) - (Pen_Ward3))
         endif
c gauge4
         if(abs(Pen_Ward4).gt.1d-2) then
         gauge4=abs((box4_for_Ward4-box3_for_Ward4)/(Pen_Ward4)-1d0)
         else
         gauge4=abs((box4_for_Ward4-box3_for_Ward4) - (Pen_Ward4))
         endif


         if(gauge1+1.eq.gauge1) then
            gauge1=1d99
         endif
         if(gauge1.ne.gauge1) then
            gauge1=1d99
         endif

         if(gauge2+1.eq.gauge2) then
            gauge2=1d99
         endif
         if(gauge2.ne.gauge2) then
            gauge2=1d99
         endif

         if(gauge3+1.eq.gauge3) then
            gauge3=1d99
         endif
         if(gauge3.ne.gauge3) then
            gauge3=1d99
         endif

         if(gauge4+1.eq.gauge4) then
            gauge4=1d99
         endif
         if(gauge4.ne.gauge4) then
            gauge4=1d99
         endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         gauge=Max(gauge1,gauge2,gauge3,gauge4)
         
         if(gauge+1.eq.gauge) then
            gauge=1d99
         endif
         if(gauge.ne.gauge) then
            gauge=1d99
         endif


         if (gauge.lt.gaugelimit) then
c             PenSingular(2)=.false.
c             print*,"gauge:true:P", aux(jmass),gauge,PenSingular
         else
            
                PenSingular(1)=.true.
                aux(jmass)=0d0
                
c             print*,"gauge:false:P",aux(jmass), gauge,PenSingular(1)
          endif
cfc         print*, "box4_for_ward4",box4_for_ward4/(12d0*massj(jmass)*massj(jmass))  
cfc         print*, "box3_for_ward4",box3_for_ward4/(12d0*massj(jmass)*massj(jmass))  
cfc         print*, "pen_ward4",pen_ward4/(3d0*massj(jmass)*massj(jmass)) 
cfc
cfc         print*, "pen1_for_ward1",pen1_for_ward1/(12d0*massj(jmass)*massj(jmass))  
cfc         print*, "pen5_for_ward1",pen5_for_ward1/(12d0*massj(jmass)*massj(jmass))  
cfc         print*, "hex_ward1",Hex_ward1/(3d0*massj(jmass)*massj(jmass)) 
cfc         print*, "gauge1",gauge1
cfc         print*, "gauge2",gauge2
cfc         print*, "gauge3",gauge3
cfc         print*, "gauge4",gauge4
c         print*, "gauge5",gauge5
c         print*, "aux(jmass)", aux(jmass)
cfc         print*, "pen_Ward1",pen_Ward1
cfc         print*, "pen_Ward2",pen_Ward2
cfc         print*, "pen_Ward3",pen_Ward3
cfc         print*, "pen_Ward4",pen_Ward4
cc end check
          endif ! dble1


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    START QUAD precision
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef WITH_QUAD
          
          If(QUAD.or.pensingular(1)) then

 

         CALL gluonpenHiggsEven_QUAD(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),
     &        Kpen1(0,i3),Kpen1(0,i4),PH1,epsgPen(0,i1),epsgPen(0,i2),
     &        epsgPen(0,i3),epsgPen(0,i4),pen1,1)
         pen1=-cplmassj(jmass)*pen1/(12d0*massj(jmass)*massj(jmass))
c         print*, "pen1",abs(pen1/aux(jmass)-1d0)
         aux(jmass)= pen1
c start Ward identity
c The momenta of the Hexline is going to be the same. 
c                       --->  comp=-1 to avoid redundant calculation
c First Gauge test Eps_i -->   p_i
c +  auxiliary momenta
c         print*, "i1--i5", i1,i2,i3,i4,i5
         do jtemp=0,3
c            print*, "epsgC(0,i1)",epsgC(jtemp,i1)
            mup1temp(jtemp)= Kpen1(jtemp,i1)
            mup2temp(jtemp)= Kpen1(jtemp,i2)
            mup3temp(jtemp)= Kpen1(jtemp,i3)
            mup4temp(jtemp)= Kpen1(jtemp,i4)
            pp12(jtemp)=Kpen1(jtemp,i1)+Kpen1(jtemp,i2)
            pp23(jtemp)=Kpen1(jtemp,i2)+Kpen1(jtemp,i3)
            pp34(jtemp)=Kpen1(jtemp,i3)+Kpen1(jtemp,i4)
            pp45(jtemp)=-(Kpen1(jtemp,i1)+Kpen1(jtemp,i2)+Kpen1(jtemp,i3))
            pp15(jtemp)=-(Kpen1(jtemp,i2)+Kpen1(jtemp,i3)+Kpen1(jtemp,i4))
         enddo   

c Ward 1
         CALL gluonpenHiggsEven_QUAD(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),
     &        Kpen1(0,i3),Kpen1(0,i4),PH1,mup1temp,epsgPen(0,i2),
     &        epsgPen(0,i3),epsgPen(0,i4),pen_Ward1,-1)
c Ward 2
         CALL gluonpenHiggsEven_QUAD(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),
     &        Kpen1(0,i3),Kpen1(0,i4),PH1,epsgPen(0,i1),mup2temp,
     &        epsgPen(0,i3),epsgPen(0,i4),pen_Ward2,-1)
c Ward 3
         CALL gluonpenHiggsEven_QUAD(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),
     &        Kpen1(0,i3),Kpen1(0,i4),PH1,epsgPen(0,i1),
     &        epsgPen(0,i2),mup3temp,epsgPen(0,i4),pen_Ward3,-1)
c Ward 4
         CALL gluonpenHiggsEven_QUAD(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),
     &        Kpen1(0,i3),Kpen1(0,i4),PH1,epsgPen(0,i1),
     &        epsgPen(0,i2),epsgPen(0,i3),mup4temp,pen_Ward4,-1)
c         print*, "aux(jmass)",aux(jmass)
c         print*, "pen1Ward1",pen_Ward1
c         print*, "pen1Ward2",pen_Ward2
c         print*, "pen1Ward3",pen_Ward3
c         print*, "pen1Ward4",pen_Ward4
c         stop
c Call of pentagons
c  Each pentagon will appear twice. eg., pen1:  pen2 -pen1,  pen1 -pen0
c wich different sign depending for which contraction is used.
c      print*, "epsgC(0,i5)_2",epsgC(0,i5),epsgC(1,i5),epsgC(2,i5),epsgC(3,i5)
c  pp12
       CALL gluonboxHiggsEven_QUAD(massj(jmass),pp12,Kpen1(0,i3),Kpen1(0,i4),PH1,
     &  epsgPen(0,i2),epsgPen(0,i3),epsgPen(0,i4),box1_for_Ward1,1)
       CALL gluonboxHiggsEven_QUAD(massj(jmass),pp12,Kpen1(0,i3),Kpen1(0,i4),PH1,
     &  epsgPen(0,i1),epsgPen(0,i3),epsgPen(0,i4),box1_for_Ward2,-1)
c  pp23
       CALL gluonboxHiggsEven_QUAD(massj(jmass),Kpen1(0,i1),pp23,Kpen1(0,i4),PH1,
     &  epsgPen(0,i1),epsgPen(0,i3),epsgPen(0,i4),box2_for_Ward2,1)
       CALL gluonboxHiggsEven_QUAD(massj(jmass),Kpen1(0,i1),pp23,Kpen1(0,i4),PH1,
     &  epsgPen(0,i1),epsgPen(0,i2),epsgPen(0,i4),box2_for_Ward3,-1)
c  pp34
       CALL gluonboxHiggsEven_QUAD(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),pp34,PH1,
     &  epsgPen(0,i1),epsgPen(0,i2),epsgPen(0,i4),box3_for_Ward3,1)
       CALL gluonboxHiggsEven_QUAD(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),pp34,PH1,
     &  epsgPen(0,i1),epsgPen(0,i2),epsgPen(0,i3),box3_for_Ward4,-1)
c  pp45
       CALL gluonboxHiggsEven_QUAD(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),Kpen1(0,i3),pp45,
     &  epsgPen(0,i1),epsgPen(0,i2),epsgPen(0,i3),box4_for_Ward4,1)
c pp15
       CALL gluonboxHiggsEven_QUAD(massj(jmass),Kpen1(0,i2),Kpen1(0,i3),Kpen1(0,i4),pp15,
     &  epsgPen(0,i2),epsgPen(0,i3),epsgPen(0,i4),box4_for_Ward1,1)

c gauge1
         if(abs(Pen_Ward1).gt.1d-2) then
         gauge1=abs((box1_for_Ward1-box4_for_Ward1)/(Pen_Ward1)-1d0)
         else
         gauge1=abs((box1_for_Ward1-box4_for_Ward1) - (Pen_Ward1))
         endif
c gauge2
         if(abs(Pen_Ward2).gt.1d-2) then
         gauge2=abs((box2_for_Ward2-box1_for_Ward2)/(Pen_Ward2)-1d0)
         else
         gauge2=abs((box2_for_Ward2-box1_for_Ward2) - (Pen_Ward2))
         endif
c gauge3
         if(abs(Pen_Ward3).gt.1d-2) then
         gauge3=abs((box3_for_Ward3-box2_for_Ward3)/(Pen_Ward3)-1d0)
         else
         gauge3=abs((box3_for_Ward3-box2_for_Ward3) - (Pen_Ward3))
         endif
c gauge4
         if(abs(Pen_Ward4).gt.1d-2) then
         gauge4=abs((box4_for_Ward4-box3_for_Ward4)/(Pen_Ward4)-1d0)
         else
         gauge4=abs((box4_for_Ward4-box3_for_Ward4) - (Pen_Ward4))
         endif

        if(gauge1+1.eq.gauge1) then
            gauge1=1d99
         endif
         if(gauge1.ne.gauge1) then
            gauge1=1d99
         endif

         if(gauge2+1.eq.gauge2) then
            gauge2=1d99
         endif
         if(gauge2.ne.gauge2) then
            gauge2=1d99
         endif

         if(gauge3+1.eq.gauge3) then
            gauge3=1d99
         endif
         if(gauge3.ne.gauge3) then
            gauge3=1d99
         endif

         if(gauge4+1.eq.gauge4) then
            gauge4=1d99
         endif
         if(gauge4.ne.gauge4) then
            gauge4=1d99
         endif


         gauge=Max(gauge1,gauge2,gauge3,gauge4)


         if(gauge+1.eq.gauge) then
            gauge=1d99
         endif
         if(gauge.ne.gauge) then
            gauge=1d99
         endif

         if (gauge.lt.gaugelimit) then
c             PenSingular(2)=.false.
c            if(PenSingular(1))then
c             print*,"gauge:true:PQ",aux(jmass),gauge, PenSingular
c            endif
             else
                PenSingular(2)=.true.
                aux(jmass)= 0d0
              print*,"gauge:false:PQ",aux(jmass),gauge, PenSingular(2)
          endif



       endif                    ! QUAD

#endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    END  QUAD precision
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


          
          enddo
cc start check
c      print*,'qpentH_top(1)',qpentH_top(1)
c      print*,'qpentH_top(2)',qpentH_top(2)
c      print*,'qpentH_top(3)',qpentH_top(3)
c      print*,'qpentH_top(4)',qpentH_top(4)
c      STOP
cc end check 
c
C
C
C
C
C
          endif
C
C
C
C
C
C
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
C
C
C
C               CP ODD CASE
C
C
C
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
 100      shift=2
       if(AloopTop.or.AloopBottom) then

      if (AloopTop) then
       initdo=1
       else
       initdo=2
      endif
      if(AloopBottom) then
         finaldo=2
      else
        finaldo=1
      endif

      If(SelectAbottom) then
         initdo=2
         finaldo=2
      endif
      If(SelectAtop) then
         initdo=1
         finaldo=1
      endif

      If(SelectHTop.or.SelectHBottom) then
         go to 200
      endif


      
      massj(shift+1)=m_top
      massj(shift+2)=m_bottom
      cplmassj(shift+1)=cplModAtop
      cplmassj(shift+2)=cplModAbottom
      

      do jmass=shift+initdo,shift+finaldo

c         print*, "jamss", jmass
        If(dble1) then

c         CALL H4G_LoopDD(i4,i3,i2,i1,m_top,qpentH_top(1))
c         CALL H4G_LoopDD(i2,i3,i4,i1,m_top,qpentH_top(2))
c         CALL H4G_LoopDD(i2,i1,i4,i3,m_top,qpentH_top(3))
c         CALL H4G_LoopDD(i3,i2,i1,i4,m_top,qpentH_top(4))
c         aux(1)=cplModHtop
c     &    *(qpentH_top(1)+qpentH_top(2)+qpentH_top(3)+qpentH_top(4))
c         CALL H4G_LoopDD(i1,i2,i3,i4,m_top,qpentH_top(1))
c         aux(1)=cplModHtop*qpentH_top(1)
cc start check
c         aux(1)=-cplModHtop*qpentH_top(1)
c
         CALL gluonpenHiggsOdd(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),
     &        Kpen1(0,i3),Kpen1(0,i4),PH1,epsgPen(0,i1),epsgPen(0,i2),
     &        epsgPen(0,i3),epsgPen(0,i4),pen1,1)

c         print*,"pen:mass", massj(jmass),cplmassj(jmass)/massj(jmass)**2
c         pause
         pen1=-cplmassj(jmass)*pen1/(12d0*massj(jmass)*massj(jmass))
c         print*, "pen1",abs(pen1/aux(jmass)-1d0)
         aux(jmass)= pen1
c start Ward identity
c The momenta of the Hexline is going to be the same. 
c                       --->  comp=-1 to avoid redundant calculation
c First Gauge test Eps_i -->   p_i
c +  auxiliary momenta
c         print*, "i1--i5", i1,i2,i3,i4,i5
         do jtemp=0,3
c            print*, "epsgC(0,i1)",epsgC(jtemp,i1)
            mup1temp(jtemp)= Kpen1(jtemp,i1)
            mup2temp(jtemp)= Kpen1(jtemp,i2)
            mup3temp(jtemp)= Kpen1(jtemp,i3)
            mup4temp(jtemp)= Kpen1(jtemp,i4)
            pp12(jtemp)=Kpen1(jtemp,i1)+Kpen1(jtemp,i2)
            pp23(jtemp)=Kpen1(jtemp,i2)+Kpen1(jtemp,i3)
            pp34(jtemp)=Kpen1(jtemp,i3)+Kpen1(jtemp,i4)
            pp45(jtemp)=-(Kpen1(jtemp,i1)+Kpen1(jtemp,i2)+Kpen1(jtemp,i3))
            pp15(jtemp)=-(Kpen1(jtemp,i2)+Kpen1(jtemp,i3)+Kpen1(jtemp,i4))
         enddo   

c Ward 1
         CALL gluonpenHiggsOdd(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),
     &        Kpen1(0,i3),Kpen1(0,i4),PH1,mup1temp,epsgPen(0,i2),
     &        epsgPen(0,i3),epsgPen(0,i4),pen_Ward1,-1)
c Ward 2
         CALL gluonpenHiggsOdd(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),
     &        Kpen1(0,i3),Kpen1(0,i4),PH1,epsgPen(0,i1),mup2temp,
     &        epsgPen(0,i3),epsgPen(0,i4),pen_Ward2,-1)
c Ward 3
         CALL gluonpenHiggsOdd(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),
     &        Kpen1(0,i3),Kpen1(0,i4),PH1,epsgPen(0,i1),
     &        epsgPen(0,i2),mup3temp,epsgPen(0,i4),pen_Ward3,-1)
c Ward 4
         CALL gluonpenHiggsOdd(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),
     &        Kpen1(0,i3),Kpen1(0,i4),PH1,epsgPen(0,i1),
     &        epsgPen(0,i2),epsgPen(0,i3),mup4temp,pen_Ward4,-1)
c         print*, "aux(jmass)",aux(jmass)
c         print*, "pen1Ward1",pen_Ward1
c         print*, "pen1Ward2",pen_Ward2
c         print*, "pen1Ward3",pen_Ward3
c         print*, "pen1Ward4",pen_Ward4
c         stop
c Call of pentagons
c  Each pentagon will appear twice. eg., pen1:  pen2 -pen1,  pen1 -pen0
c wich different sign depending for which contraction is used.
c      print*, "epsgC(0,i5)_2",epsgC(0,i5),epsgC(1,i5),epsgC(2,i5),epsgC(3,i5)
c  pp12
       CALL gluonboxHiggsOdd(massj(jmass),pp12,Kpen1(0,i3),Kpen1(0,i4),PH1,
     &  epsgPen(0,i2),epsgPen(0,i3),epsgPen(0,i4),box1_for_Ward1,1)
       CALL gluonboxHiggsOdd(massj(jmass),pp12,Kpen1(0,i3),Kpen1(0,i4),PH1,
     &  epsgPen(0,i1),epsgPen(0,i3),epsgPen(0,i4),box1_for_Ward2,-1)
c  pp23
       CALL gluonboxHiggsOdd(massj(jmass),Kpen1(0,i1),pp23,Kpen1(0,i4),PH1,
     &  epsgPen(0,i1),epsgPen(0,i3),epsgPen(0,i4),box2_for_Ward2,1)
       CALL gluonboxHiggsOdd(massj(jmass),Kpen1(0,i1),pp23,Kpen1(0,i4),PH1,
     &  epsgPen(0,i1),epsgPen(0,i2),epsgPen(0,i4),box2_for_Ward3,-1)
c  pp34
       CALL gluonboxHiggsOdd(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),pp34,PH1,
     &  epsgPen(0,i1),epsgPen(0,i2),epsgPen(0,i4),box3_for_Ward3,1)
       CALL gluonboxHiggsOdd(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),pp34,PH1,
     &  epsgPen(0,i1),epsgPen(0,i2),epsgPen(0,i3),box3_for_Ward4,-1)
c  pp45
       CALL gluonboxHiggsOdd(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),Kpen1(0,i3),pp45,
     &  epsgPen(0,i1),epsgPen(0,i2),epsgPen(0,i3),box4_for_Ward4,1)
c pp15
       CALL gluonboxHiggsOdd(massj(jmass),Kpen1(0,i2),Kpen1(0,i3),Kpen1(0,i4),pp15,
     &  epsgPen(0,i2),epsgPen(0,i3),epsgPen(0,i4),box4_for_Ward1,1)

c gauge1
         if(abs(Pen_Ward1).gt.1d-2) then
         gauge1=abs((box1_for_Ward1-box4_for_Ward1)/(Pen_Ward1)-1d0)
         else
         gauge1=abs((box1_for_Ward1-box4_for_Ward1) - (Pen_Ward1))
         endif
c gauge2
         if(abs(Pen_Ward2).gt.1d-2) then
         gauge2=abs((box2_for_Ward2-box1_for_Ward2)/(Pen_Ward2)-1d0)
         else
         gauge2=abs((box2_for_Ward2-box1_for_Ward2) - (Pen_Ward2))
         endif
c gauge3
         if(abs(Pen_Ward3).gt.1d-2) then
         gauge3=abs((box3_for_Ward3-box2_for_Ward3)/(Pen_Ward3)-1d0)
         else
         gauge3=abs((box3_for_Ward3-box2_for_Ward3) - (Pen_Ward3))
         endif
c gauge4
         if(abs(Pen_Ward4).gt.1d-2) then
         gauge4=abs((box4_for_Ward4-box3_for_Ward4)/(Pen_Ward4)-1d0)
         else
         gauge4=abs((box4_for_Ward4-box3_for_Ward4) - (Pen_Ward4))
         endif


         if(gauge1+1.eq.gauge1) then
            gauge1=1d99
         endif
         if(gauge1.ne.gauge1) then
            gauge1=1d99
         endif


         if(gauge2+1.eq.gauge2) then
            gauge2=1d99
         endif
         if(gauge2.ne.gauge2) then
            gauge2=1d99
         endif


         if(gauge3+1.eq.gauge3) then
            gauge3=1d99
         endif
         if(gauge3.ne.gauge3) then
            gauge3=1d99
         endif

         if(gauge4+1.eq.gauge4) then
            gauge4=1d99
         endif
         if(gauge4.ne.gauge4) then
            gauge4=1d99
         endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
         gauge=Max(gauge1,gauge2,gauge3,gauge4)

         
         if(gauge+1.eq.gauge) then
            gauge=1d99
         endif
        if(gauge.ne.gauge) then
            gauge=1d99
         endif

         if (gauge.lt.gaugelimit) then
c             PenSingular(2)=.false.
c             print*,"gauge:true:P", aux(jmass),gauge,PenSingular
             else
                PenSingular(1)=.true.
                aux(jmass)=0d0
c             print*,"gauge:false:P",aux(jmass), gauge,PenSingular(1)
          endif
cfc         print*, "box4_for_ward4",box4_for_ward4/(12d0*massj(jmass)*massj(jmass))  
cfc         print*, "box3_for_ward4",box3_for_ward4/(12d0*massj(jmass)*massj(jmass))  
cfc         print*, "pen_ward4",pen_ward4/(3d0*massj(jmass)*massj(jmass)) 
cfc
cfc         print*, "pen1_for_ward1",pen1_for_ward1/(12d0*massj(jmass)*massj(jmass))  
cfc         print*, "pen5_for_ward1",pen5_for_ward1/(12d0*massj(jmass)*massj(jmass))  
cfc         print*, "hex_ward1",Hex_ward1/(3d0*massj(jmass)*massj(jmass)) 
cfc         print*, "gauge1",gauge1
cfc         print*, "gauge2",gauge2
cfc         print*, "gauge3",gauge3
cfc         print*, "gauge4",gauge4
c         print*, "gauge5",gauge5
c         print*, "aux(jmass)", aux(jmass)
cfc         print*, "pen_Ward1",pen_Ward1
cfc         print*, "pen_Ward2",pen_Ward2
cfc         print*, "pen_Ward3",pen_Ward3
cfc         print*, "pen_Ward4",pen_Ward4
cc end check
          endif ! dble1


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    START QUAD precision
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef WITH_QUAD
          If(QUAD.or.pensingular(1)) then
 

         CALL gluonpenHiggsOdd_QUAD(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),
     &        Kpen1(0,i3),Kpen1(0,i4),PH1,epsgPen(0,i1),epsgPen(0,i2),
     &        epsgPen(0,i3),epsgPen(0,i4),pen1,1)
         pen1=-cplmassj(jmass)*pen1/(12d0*massj(jmass)*massj(jmass))
c         print*, "pen1",abs(pen1/aux(jmass)-1d0)
         aux(jmass)= pen1
c start Ward identity
c The momenta of the Hexline is going to be the same. 
c                       --->  comp=-1 to avoid redundant calculation
c First Gauge test Eps_i -->   p_i
c +  auxiliary momenta
c         print*, "i1--i5", i1,i2,i3,i4,i5
         do jtemp=0,3
c            print*, "epsgC(0,i1)",epsgC(jtemp,i1)
            mup1temp(jtemp)= Kpen1(jtemp,i1)
            mup2temp(jtemp)= Kpen1(jtemp,i2)
            mup3temp(jtemp)= Kpen1(jtemp,i3)
            mup4temp(jtemp)= Kpen1(jtemp,i4)
            pp12(jtemp)=Kpen1(jtemp,i1)+Kpen1(jtemp,i2)
            pp23(jtemp)=Kpen1(jtemp,i2)+Kpen1(jtemp,i3)
            pp34(jtemp)=Kpen1(jtemp,i3)+Kpen1(jtemp,i4)
            pp45(jtemp)=-(Kpen1(jtemp,i1)+Kpen1(jtemp,i2)+Kpen1(jtemp,i3))
            pp15(jtemp)=-(Kpen1(jtemp,i2)+Kpen1(jtemp,i3)+Kpen1(jtemp,i4))
         enddo   

c Ward 1
         CALL gluonpenHiggsOdd_QUAD(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),
     &        Kpen1(0,i3),Kpen1(0,i4),PH1,mup1temp,epsgPen(0,i2),
     &        epsgPen(0,i3),epsgPen(0,i4),pen_Ward1,-1)
c Ward 2
         CALL gluonpenHiggsOdd_QUAD(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),
     &        Kpen1(0,i3),Kpen1(0,i4),PH1,epsgPen(0,i1),mup2temp,
     &        epsgPen(0,i3),epsgPen(0,i4),pen_Ward2,-1)
c Ward 3
         CALL gluonpenHiggsOdd_QUAD(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),
     &        Kpen1(0,i3),Kpen1(0,i4),PH1,epsgPen(0,i1),
     &        epsgPen(0,i2),mup3temp,epsgPen(0,i4),pen_Ward3,-1)
c Ward 4
         CALL gluonpenHiggsOdd_QUAD(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),
     &        Kpen1(0,i3),Kpen1(0,i4),PH1,epsgPen(0,i1),
     &        epsgPen(0,i2),epsgPen(0,i3),mup4temp,pen_Ward4,-1)
c         print*, "aux(jmass)",aux(jmass)
c         print*, "pen1Ward1",pen_Ward1
c         print*, "pen1Ward2",pen_Ward2
c         print*, "pen1Ward3",pen_Ward3
c         print*, "pen1Ward4",pen_Ward4
c         stop
c Call of pentagons
c  Each pentagon will appear twice. eg., pen1:  pen2 -pen1,  pen1 -pen0
c wich different sign depending for which contraction is used.
c      print*, "epsgC(0,i5)_2",epsgC(0,i5),epsgC(1,i5),epsgC(2,i5),epsgC(3,i5)
c  pp12
       CALL gluonboxHiggsOdd_QUAD(massj(jmass),pp12,Kpen1(0,i3),Kpen1(0,i4),PH1,
     &  epsgPen(0,i2),epsgPen(0,i3),epsgPen(0,i4),box1_for_Ward1,1)
       CALL gluonboxHiggsOdd_QUAD(massj(jmass),pp12,Kpen1(0,i3),Kpen1(0,i4),PH1,
     &  epsgPen(0,i1),epsgPen(0,i3),epsgPen(0,i4),box1_for_Ward2,-1)
c  pp23
       CALL gluonboxHiggsOdd_QUAD(massj(jmass),Kpen1(0,i1),pp23,Kpen1(0,i4),PH1,
     &  epsgPen(0,i1),epsgPen(0,i3),epsgPen(0,i4),box2_for_Ward2,1)
       CALL gluonboxHiggsOdd_QUAD(massj(jmass),Kpen1(0,i1),pp23,Kpen1(0,i4),PH1,
     &  epsgPen(0,i1),epsgPen(0,i2),epsgPen(0,i4),box2_for_Ward3,-1)
c  pp34
       CALL gluonboxHiggsOdd_QUAD(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),pp34,PH1,
     &  epsgPen(0,i1),epsgPen(0,i2),epsgPen(0,i4),box3_for_Ward3,1)
       CALL gluonboxHiggsOdd_QUAD(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),pp34,PH1,
     &  epsgPen(0,i1),epsgPen(0,i2),epsgPen(0,i3),box3_for_Ward4,-1)
c  pp45
       CALL gluonboxHiggsOdd_QUAD(massj(jmass),Kpen1(0,i1),Kpen1(0,i2),Kpen1(0,i3),pp45,
     &  epsgPen(0,i1),epsgPen(0,i2),epsgPen(0,i3),box4_for_Ward4,1)
c pp15
       CALL gluonboxHiggsOdd_QUAD(massj(jmass),Kpen1(0,i2),Kpen1(0,i3),Kpen1(0,i4),pp15,
     &  epsgPen(0,i2),epsgPen(0,i3),epsgPen(0,i4),box4_for_Ward1,1)

c gauge1
         if(abs(Pen_Ward1).gt.1d-2) then
         gauge1=abs((box1_for_Ward1-box4_for_Ward1)/(Pen_Ward1)-1d0)
         else
         gauge1=abs((box1_for_Ward1-box4_for_Ward1) - (Pen_Ward1))
         endif
c gauge2
         if(abs(Pen_Ward2).gt.1d-2) then
         gauge2=abs((box2_for_Ward2-box1_for_Ward2)/(Pen_Ward2)-1d0)
         else
         gauge2=abs((box2_for_Ward2-box1_for_Ward2) - (Pen_Ward2))
         endif
c gauge3
         if(abs(Pen_Ward3).gt.1d-2) then
         gauge3=abs((box3_for_Ward3-box2_for_Ward3)/(Pen_Ward3)-1d0)
         else
         gauge3=abs((box3_for_Ward3-box2_for_Ward3) - (Pen_Ward3))
         endif
c gauge4
         if(abs(Pen_Ward4).gt.1d-2) then
         gauge4=abs((box4_for_Ward4-box3_for_Ward4)/(Pen_Ward4)-1d0)
         else
         gauge4=abs((box4_for_Ward4-box3_for_Ward4) - (Pen_Ward4))
         endif

        if(gauge1+1.eq.gauge1) then
            gauge1=1d99
         endif
         if(gauge1.ne.gauge1) then
            gauge1=1d99
         endif


         if(gauge2+1.eq.gauge2) then
            gauge2=1d99
         endif
         if(gauge2.ne.gauge2) then
            gauge2=1d99
         endif


         if(gauge3+1.eq.gauge3) then
            gauge3=1d99
         endif
         if(gauge3.ne.gauge3) then
            gauge3=1d99
         endif

         if(gauge4+1.eq.gauge4) then
            gauge4=1d99
         endif
         if(gauge4.ne.gauge4) then
            gauge4=1d99
         endif


         gauge=Max(gauge1,gauge2,gauge3,gauge4)

c         print*, "gauge:pen:quad",gauge
c         pause


         if(gauge+1.eq.gauge) then
            gauge=1d99
         endif
         if(gauge.ne.gauge) then
            gauge=1d99
         endif


         if (gauge.lt.gaugelimit) then
c             PenSingular(2)=.false.
c            if(PenSingular(1))then
c             print*,"gauge:true:PQ",aux(jmass),gauge, PenSingular
c            endif
             else
                PenSingular(2)=.true.
                aux(jmass)=0d0
              print*,"gauge:false:PQ",aux(jmass),gauge, PenSingular(2)
          endif



       endif                    ! QUAD

#endif       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    END  QUAD precision
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          enddo




          endif







 

 200  do i=1,26
          FDeffPenLoop=FDeffPenLoop+aux(i)
      enddo
      
c      print*, "FDeffPenLoop",FDeffPenLoop,aux(3)
      RETURN
      END


***********************************************************************
      complex*16 FUNCTION FDHexLoop(i1,i2,i3,i4,i5)
*----------------------------------------------------------------------
c... Date: 09.07.2012      
***********************************************************************
      
      IMPLICIT NONE

#include "ggf_process.inc"

      integer i1,i2,i3,i4,i5,i,comp
      complex*16 aux(26),qhexH_top,qhexH_bottom
     
      real*8 K(0:4,5),Kij(0:4,5,5),Kijk(0:4,5,5,5),PH(0:4),KPH(0:4,5),
     &     epsg(0:3,5),
     &     KPHeffij(0:4,5,5),KPAeffij(0:4,5,5),KPHtopij(0:4,5,5),
     &     KPHbottomij(0:4,5,5),KPAtopij(0:4,5,5),KPAbottompij(0:4,5,5),
     &     Kpen(0:4,4)
c... 4-dim momentum variables
      real*8 Kp(0:3,5),Kijp(0:3,5,5),Kijkp(0:3,5,5,5),PHp(0:3),
     &     KPHp(0:3,5),K123p(0:3),KPHeffijp(0:3,5,5),KPAeffijp(0:3,5,5),
     &     KPHtopijp(0:3,5,5),KPHbottomijp(0:3,5,5),KPAtopijp(0:3,5,5),
     &     KPAbottomijp(0:3,5,5)
      complex*16 epsgC(0:3,5),Jij(0:3,5,5),Jijk(0:3,5,5,5),
     &           epsgPen(0:3,4)
      complex*16 epsHeff(0:3,5),epsAeff(0:3,5),epsHtop(0:3,5),
     &     epsHbottom(0:3,5),epsAtop(0:3,5),epsAbottom(0:3,5),
     &     JPHeffij(0:3,5,5),JPAeffij(0:3,5,5),JPHtopij(0:3,5,5),
     &     JPHbottomij(0:3,5,5),JPAtopij(0:3,5,5),JPAbottomij(0:3,5,5)
      COMMON /GG3JET/ K,Kij,Kijk,epsg,epsgC,Jij,Jijk,KPH,PH,
     &     KPHeffij,KPAeffij,KPHtopij,KPHbottomij,KPAtopij,
     &     KPAbottompij,
     &     epsHeff,epsAeff,epsHtop,epsHbottom,epsAtop,epsAbottom,
     &     JPHeffij,JPAeffij,JPHtopij,JPHbottomij,JPAtopij,
     &     JPAbottomij,
c... 4-dim momentum variables
     &     Kp,Kijp,Kijkp,PHp,KPHp,KPHeffijp,KPAeffijp,K123p,
     &     KPHtopijp,KPHbottomijp,KPAtopijp,KPAbottomijp

      Common/Forpent/Kpen,epsgPen

      complex*16 mup1temp(0:3),mup2temp(0:3),mup3temp(0:3),mup4temp(0:3)
      complex*16 mup5temp(0:3)
      INTEGER jtemp
      complex*16 Hex_Ward1,Hex_Ward2,Hex_Ward3,Hex_Ward4,Hex_Ward5
      real*8 pp12(0:3),pp23(0:3),pp34(0:3),pp45(0:3),pp56(0:3),pp16(0:3)
      real*8 gauge1,gauge2,gauge3,gauge4,gauge5
      real*8 gauge
      complex*16 pen1_for_Ward1, pen1_for_Ward2
      complex*16 pen2_for_ward2, pen2_for_ward3
      complex*16 pen3_for_ward3, pen3_for_ward4
      complex*16 pen4_for_ward4, pen4_for_ward5
      complex*16 pen5_for_ward5, pen5_for_ward1
      logical QUAD,Dble1
      real*8 K1(0:3,5),PH1(0:3)
       real*8 gaugelimit 
       common/gauge/gaugelimit
      common/dbleswitch/dble1,quad

      integer initdo,finaldo,jmass
      real*8 massj(8),cplmassj(8)
      integer shift
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      HexSingular(1)=.true.
c      HexSingular(2)=.true.

      do jtemp=1,5
      do i=0,3
         K1(i,jtemp)=K(i,jtemp)
         PH1(i)=PH(i)
      enddo
      enddo
      do i=0,3
        PH1(i)=PH(i)
      enddo

      comp=1

      FDHexLoop=DCMPLX(0d0,0d0)
      do i=1,26
         aux(i)=DCMPLX(0d0,0d0)
      enddo
   



      if (HloopTop) then
       initdo=1
       else
       initdo=2
      endif
      if(HloopBottom) then
         finaldo=2
      else
        finaldo=1
      endif
      massj(1)=m_top
      massj(2)=m_bottom
      cplmassj(1)=cplModHtop
      cplmassj(2)=cplModHbottom
      

      do jmass=initdo,finaldo




         if(Dble1) then
         
         CALL gluonHexHiggsEven(massj(jmass),K1(0,i1),K1(0,i2),K1(0,i3),K1(0,i4),
     &        K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i2),epsgC(0,i3),
     &        epsgC(0,i4),epsgC(0,i5),1,qhexH_top)
         aux(jmass)=-cplmassj(jmass)*qhexH_top/(4d0*massj(jmass)*massj(jmass))  
 

cc start check
c         aux(jmass)=-cplModHtop*qhexH_top/(12d0*massj(jmass)*massj(jmass))  
cc end check
c start Ward identity
c The momenta of the Hexline is going to be the same. 
c                       --->  comp=-1 to avoid redundant calculation
c
c
c First Gauge test Eps_i -->   p_i
c +  auxiliary momenta
c         print*, "i1--i5", i1,i2,i3,i4,i5
         do jtemp=0,3
c            print*, "epsgC(0,i1)",epsgC(jtemp,i1)
            mup1temp(jtemp)= K1(jtemp,i1)
            mup2temp(jtemp)= K1(jtemp,i2)
            mup3temp(jtemp)= K1(jtemp,i3)
            mup4temp(jtemp)= K1(jtemp,i4)
            mup5temp(jtemp)= K1(jtemp,i5)
            pp12(jtemp)=K1(jtemp,i1)+K1(jtemp,i2)
            pp23(jtemp)=K1(jtemp,i2)+K1(jtemp,i3)
            pp34(jtemp)=K1(jtemp,i3)+K1(jtemp,i4)
            pp45(jtemp)=K1(jtemp,i4)+K1(jtemp,i5)
            pp56(jtemp)=-(K1(jtemp,i1)+K1(jtemp,i2)+K1(jtemp,i3)+K1(jtemp,i4))
            pp16(jtemp)=-(K1(jtemp,i2)+K1(jtemp,i3)+K1(jtemp,i4)+K1(jtemp,i5))
           enddo   
c   Ward 1
         CALL gluonHexHiggsEven(massj(jmass),K1(0,i1),K1(0,i2),K1(0,i3),K1(0,i4),
     &        K1(0,i5),PH1,mup1temp,epsgC(0:3,i2),epsgC(0:3,i3),
     &        epsgC(0:3,i4),epsgC(0:3,i5),-1,Hex_Ward1)
c   Ward 2
         CALL gluonHexHiggsEven(massj(jmass),K1(0,i1),K1(0,i2),K1(0,i3),K1(0,i4),
     &        K1(0,i5),PH1,epsgC(0,i1),mup2temp,epsgC(0,i3),
     &        epsgC(0,i4),epsgC(0,i5),-1,Hex_Ward2)
c   Ward 3
         CALL gluonHexHiggsEven(massj(jmass),K1(0,i1),K1(0,i2),K1(0,i3),K1(0,i4),
     &        K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i2),mup3temp,
     &        epsgC(0,i4),epsgC(0,i5),-1,Hex_Ward3)
c   Ward 4
         CALL gluonHexHiggsEven(massj(jmass),K1(0,i1),K1(0,i2),K1(0,i3),K1(0,i4),
     &        K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i2),epsgC(0,i3),
     &        mup4temp,epsgC(0,i5),-1,Hex_Ward4)
c   Ward 5
         CALL gluonHexHiggsEven(massj(jmass),K1(0,i1),K1(0,i2),K1(0,i3),K1(0,i4),
     &        K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i2),epsgC(0,i3),
     &        epsgC(0,i4),mup5temp,-1,Hex_Ward5)

c Call of pentagons
c  Each pentagon will appear twice. eg., pen1:  pen2 -pen1,  pen1 -pen0
c wich different sign depending for which contraction is used.
c      print*, "epsgC(0,i5)_2",epsgC(0,i5),epsgC(1,i5),epsgC(2,i5),epsgC(3,i5)
c  pp12
         CALL gluonpenHiggsEven(massj(jmass),pp12,K1(0,i3),
     &        K1(0,i4),K1(0,i5),PH1,epsgC(0,i2),epsgC(0,i3),
     &        epsgC(0,i4),epsgC(0,i5),pen1_for_Ward1,1)

         CALL gluonpenHiggsEven(massj(jmass),pp12,K1(0,i3),
     &        K1(0,i4),K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i3),
     &        epsgC(0,i4),epsgC(0,i5),pen1_for_Ward2,-1)

c  pp23
         CALL gluonpenHiggsEven(massj(jmass),K1(0,i1),pp23,
     &        K1(0,i4),K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i3),
     &        epsgC(0,i4),epsgC(0,i5),pen2_for_Ward2,1)

         CALL gluonpenHiggsEven(massj(jmass),K1(0,i1),pp23,
     &        K1(0,i4),K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i2),
     &        epsgC(0,i4),epsgC(0,i5),pen2_for_Ward3,-1)

c pp34
         CALL gluonpenHiggsEven(massj(jmass),K1(0,i1),K1(0,i2),
     &        pp34,K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i2),
     &        epsgC(0,i4),epsgC(0,i5),pen3_for_Ward3,1)

         CALL gluonpenHiggsEven(massj(jmass),K1(0,i1),K1(0,i2),
     &        pp34,K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i2),
     &        epsgC(0,i3),epsgC(0,i5),pen3_for_Ward4,-1)

c pp45
         CALL gluonpenHiggsEven(massj(jmass),K1(0,i1),K1(0,i2),
     &        K1(0,i3),pp45,PH1,epsgC(0,i1),epsgC(0,i2),
     &        epsgC(0,i3),epsgC(0,i5),pen4_for_Ward4,1)

         CALL gluonpenHiggsEven(massj(jmass),K1(0,i1),K1(0,i2),
     &        K1(0,i3),pp45,PH1,epsgC(0,i1),epsgC(0,i2),
     &        epsgC(0,i3),epsgC(0,i4),pen4_for_Ward5,-1)

c pp56
         CALL gluonpenHiggsEven(massj(jmass),K1(0,i1),K1(0,i2),
     &        K1(0,i3),K1(0,i4),pp56,epsgC(0,i1),epsgC(0,i2),
     &        epsgC(0,i3),epsgC(0,i4),pen5_for_Ward5,1)

c         CALL gluonpenHiggsEven(massj(jmass),K1(0,i2),K1(0,i3),
c     &        K1(0,i4),K1(0,i5),pp56,epsgC(0,i2),
c     &        epsgC(0,i3),epsgC(0,i4),epsgC(0,i5),pen5_for_Ward1)
c      print*, "epsgC(0,i5)_3",epsgC(0,i5),epsgC(1,i5),epsgC(2,i5),epsgC(3,i5)
c pp16 
         CALL gluonpenHiggsEven(massj(jmass),K1(0,i2),K1(0,i3),
     &        K1(0,i4),K1(0,i5),pp16,epsgC(0,i2),epsgC(0,i3),
     &        epsgC(0,i4),epsgC(0,i5),pen5_for_Ward1,1)

c      print*, "epsgC(0,i5)",epsgC(0,i5),epsgC(1,i5),epsgC(2,i5),epsgC(3,i5)
c gauge1
         if(abs(Hex_Ward1).gt.1d-2) then
         gauge1=abs((pen1_for_Ward1-pen5_for_Ward1)/(3d0*Hex_Ward1)-1d0)
         else
         gauge1=abs((pen1_for_Ward1-pen5_for_Ward1) - (3d0*Hex_Ward1))
         endif
c gauge2
         if(abs(Hex_Ward2).gt.1d-2) then
         gauge2=abs((pen2_for_Ward2-pen1_for_Ward2)/(3d0*Hex_Ward2)-1d0)
         else
         gauge2=abs((pen2_for_Ward2-pen1_for_Ward2) - (3d0*Hex_Ward2))
         endif
c gauge3
         if(abs(Hex_Ward3).gt.1d-2) then
         gauge3=abs((pen3_for_Ward3-pen2_for_Ward3)/(3d0*Hex_Ward3)-1d0)
         else
         gauge3=abs((pen3_for_Ward3-pen2_for_Ward3) - (3d0*Hex_Ward3))
         endif
c gauge4
         if(abs(Hex_Ward4).gt.1d-2) then
         gauge4=abs((pen4_for_Ward4-pen3_for_Ward4)/(3d0*Hex_Ward4)-1d0)
         else
         gauge4=abs((pen4_for_Ward4-pen3_for_Ward4) - (3d0*Hex_Ward4))
         endif
c gauge5
         if(abs(Hex_Ward5).gt.1d-2) then
         gauge5=abs((pen5_for_Ward5-pen4_for_Ward5)/(3d0*Hex_Ward5)-1d0)
         else
         gauge5=abs((pen5_for_Ward5-pen4_for_Ward5) - (3d0*Hex_Ward5))
         endif

cfc         print*, "pen5_for_ward5",pen5_for_ward5/(12d0*massj(jmass)*massj(jmass))  
cfc         print*, "pen4_for_ward5",pen4_for_ward5/(12d0*massj(jmass)*massj(jmass))  
cfc         print*, "hex_ward5",Hex_ward5/(3d0*massj(jmass)*massj(jmass)) 
cfc
cfc         print*, "pen1_for_ward1",pen1_for_ward1/(12d0*massj(jmass)*massj(jmass))  
cfc         print*, "pen5_for_ward1",pen5_for_ward1/(12d0*massj(jmass)*massj(jmass))  
cfc         print*, "hex_ward1",Hex_ward1/(3d0*massj(jmass)*massj(jmass)) 
cfc         print*, "gauge1",gauge1
cfc         print*, "gauge2",gauge2
cfc         print*, "gauge3",gauge3
cfc         print*, "gauge4",gauge4
cfc         print*, "gauge5",gauge5
cfc         print*, "aux(jmass)", aux(jmass)
cfc         print*, "Hex_Ward1",Hex_Ward1
cfc         print*, "Hex_Ward2",Hex_Ward2
cfc         print*, "Hex_Ward3",Hex_Ward3
cfc         print*, "Hex_Ward4",Hex_Ward4
cfc         print*, "Hex_Ward5",Hex_Ward5
 

check for NANS
         if(gauge1+1.eq.gauge1) then
            gauge1=1d99
         endif
         if(gauge1.ne.gauge1) then
            gauge1=1d99
         endif

         if(gauge2+1.eq.gauge2) then
            gauge2=1d99
         endif
         if(gauge2.ne.gauge2) then
            gauge2=1d99
         endif

         if(gauge3+1.eq.gauge3) then
            gauge3=1d99
         endif
         if(gauge3.ne.gauge3) then
            gauge3=1d99
         endif

         if(gauge4+1.eq.gauge4) then
            gauge4=1d99
         endif
         if(gauge4.ne.gauge4) then
            gauge4=1d99
         endif

         if(gauge5+1.eq.gauge5) then
            gauge5=1d99
         endif
         if(gauge5.ne.gauge5) then
            gauge5=1d99
         endif



          gauge=Max(gauge1,gauge2,gauge3,gauge4,gauge5)

         if(gauge+1.eq.gauge) then
            gauge=1d99
         endif
         if(gauge.ne.gauge) then
            gauge=1d99
         endif



          if (gauge.lt.gaugelimit) then
c             print*,"gauge:true:H",aux(jmass),gauge, HexSingular
             else
                HexSingular(1)=.true.
                aux(jmass)=0d0
c             print*,"gauge:false:Hd",aux(jmass),gauge, HexSingular(1)
          endif
          endif ! dble

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    START QUAD precision
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef WITH_QUAD          
          If(QUAD.or.hexSingular(1)) then
        
         comp=1

         CALL gluonHexHiggsEven_QUAD(massj(jmass),K1(0,i1),K1(0,i2),K1(0,i3),K1(0,i4),
     &        K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i2),epsgC(0,i3),
     &        epsgC(0,i4),epsgC(0,i5),comp,qhexH_top)
         aux(jmass)=-cplmassj(jmass)*qhexH_top/(4d0*massj(jmass)*massj(jmass))  
     
c start Ward identity
c The momenta of the Hexline is going to be the same. 
c                       --->  comp=-1 to avoid redundant calculation
c First Gauge test Eps_i -->   p_i
c +  auxiliary momenta
c         print*, "i1--i5", i1,i2,i3,i4,i5
         do jtemp=0,3
c            print*, "epsgC(0,i1)",epsgC(jtemp,i1)
            mup1temp(jtemp)= K1(jtemp,i1)
            mup2temp(jtemp)= K1(jtemp,i2)
            mup3temp(jtemp)= K1(jtemp,i3)
            mup4temp(jtemp)= K1(jtemp,i4)
            mup5temp(jtemp)= K1(jtemp,i5)
            pp12(jtemp)=K1(jtemp,i1)+K1(jtemp,i2)
            pp23(jtemp)=K1(jtemp,i2)+K1(jtemp,i3)
            pp34(jtemp)=K1(jtemp,i3)+K1(jtemp,i4)
            pp45(jtemp)=K1(jtemp,i4)+K1(jtemp,i5)
            pp56(jtemp)=-(K1(jtemp,i1)+K1(jtemp,i2)+K1(jtemp,i3)+K1(jtemp,i4))
            pp16(jtemp)=-(K1(jtemp,i2)+K1(jtemp,i3)+K1(jtemp,i4)+K1(jtemp,i5))
           enddo   

c         print*, "epsgC(0,i5)_1",epsgC(0,i5),epsgC(1,i5),epsgC(2,i5),epsgC(3,i5)

c   Ward 1
         CALL gluonHexHiggsEven_QUAD(massj(jmass),K1(0,i1),K1(0,i2),K1(0,i3),K1(0,i4),
     &        K1(0,i5),PH1,mup1temp,epsgC(0:3,i2),epsgC(0:3,i3),
     &        epsgC(0:3,i4),epsgC(0:3,i5),-1,Hex_Ward1)
c         print*, "Hex_Ward1/res",Hex_Ward1/qhexH_top/K1(0,i1)
c   Ward 2
         CALL gluonHexHiggsEven_QUAD(massj(jmass),K1(0,i1),K1(0,i2),K1(0,i3),K1(0,i4),
     &        K1(0,i5),PH1,epsgC(0,i1),mup2temp,epsgC(0,i3),
     &        epsgC(0,i4),epsgC(0,i5),-1,Hex_Ward2)
c   Ward 3
         CALL gluonHexHiggsEven_QUAD(massj(jmass),K1(0,i1),K1(0,i2),K1(0,i3),K1(0,i4),
     &        K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i2),mup3temp,
     &        epsgC(0,i4),epsgC(0,i5),-1,Hex_Ward3)
c   Ward 4
         CALL gluonHexHiggsEven_QUAD(massj(jmass),K1(0,i1),K1(0,i2),K1(0,i3),K1(0,i4),
     &        K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i2),epsgC(0,i3),
     &        mup4temp,epsgC(0,i5),-1,Hex_Ward4)
c   Ward 5
         CALL gluonHexHiggsEven_QUAD(massj(jmass),K1(0,i1),K1(0,i2),K1(0,i3),K1(0,i4),
     &        K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i2),epsgC(0,i3),
     &        epsgC(0,i4),mup5temp,-1,Hex_Ward5)

c Call of pentagons
c  Each pentagon will appear twice. eg., pen1:  pen2 -pen1,  pen1 -pen0
c wich different sign depending for which contraction is used.
c      print*, "epsgC(0,i5)_2",epsgC(0,i5),epsgC(1,i5),epsgC(2,i5),epsgC(3,i5)
c  pp12
         CALL gluonpenHiggsEven_QUAD(massj(jmass),pp12,K1(0,i3),
     &        K1(0,i4),K1(0,i5),PH1,epsgC(0,i2),epsgC(0,i3),
     &        epsgC(0,i4),epsgC(0,i5),pen1_for_Ward1,1)

         CALL gluonpenHiggsEven_QUAD(massj(jmass),pp12,K1(0,i3),
     &        K1(0,i4),K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i3),
     &        epsgC(0,i4),epsgC(0,i5),pen1_for_Ward2,-1)

c  pp23
         CALL gluonpenHiggsEven_QUAD(massj(jmass),K1(0,i1),pp23,
     &        K1(0,i4),K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i3),
     &        epsgC(0,i4),epsgC(0,i5),pen2_for_Ward2,1)

         CALL gluonpenHiggsEven_QUAD(massj(jmass),K1(0,i1),pp23,
     &        K1(0,i4),K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i2),
     &        epsgC(0,i4),epsgC(0,i5),pen2_for_Ward3,-1)

c pp34
         CALL gluonpenHiggsEven_QUAD(massj(jmass),K1(0,i1),K1(0,i2),
     &        pp34,K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i2),
     &        epsgC(0,i4),epsgC(0,i5),pen3_for_Ward3,1)

         CALL gluonpenHiggsEven_QUAD(massj(jmass),K1(0,i1),K1(0,i2),
     &        pp34,K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i2),
     &        epsgC(0,i3),epsgC(0,i5),pen3_for_Ward4,-1)

c pp45
         CALL gluonpenHiggsEven_QUAD(massj(jmass),K1(0,i1),K1(0,i2),
     &        K1(0,i3),pp45,PH1,epsgC(0,i1),epsgC(0,i2),
     &        epsgC(0,i3),epsgC(0,i5),pen4_for_Ward4,1)

         CALL gluonpenHiggsEven_QUAD(massj(jmass),K1(0,i1),K1(0,i2),
     &        K1(0,i3),pp45,PH1,epsgC(0,i1),epsgC(0,i2),
     &        epsgC(0,i3),epsgC(0,i4),pen4_for_Ward5,-1)

c pp56
         CALL gluonpenHiggsEven_QUAD(massj(jmass),K1(0,i1),K1(0,i2),
     &        K1(0,i3),K1(0,i4),pp56,epsgC(0,i1),epsgC(0,i2),
     &        epsgC(0,i3),epsgC(0,i4),pen5_for_Ward5,1)

c         CALL gluonpenHiggsEven(massj(jmass),K1(0,i2),K1(0,i3),
c     &        K1(0,i4),K1(0,i5),pp56,epsgC(0,i2),
c     &        epsgC(0,i3),epsgC(0,i4),epsgC(0,i5),pen5_for_Ward1)
c      print*, "epsgC(0,i5)_3",epsgC(0,i5),epsgC(1,i5),epsgC(2,i5),epsgC(3,i5)
c pp16 
         CALL gluonpenHiggsEven_QUAD(massj(jmass),K1(0,i2),K1(0,i3),
     &        K1(0,i4),K1(0,i5),pp16,epsgC(0,i2),epsgC(0,i3),
     &        epsgC(0,i4),epsgC(0,i5),pen5_for_Ward1,1)

c      print*, "epsgC(0,i5)",epsgC(0,i5),epsgC(1,i5),epsgC(2,i5),epsgC(3,i5)
c gauge1
         if(abs(Hex_Ward1).gt.1d-2) then
         gauge1=abs((pen1_for_Ward1-pen5_for_Ward1)/(3d0*Hex_Ward1)-1d0)
         else
         gauge1=abs((pen1_for_Ward1-pen5_for_Ward1) - (3d0*Hex_Ward1))
         endif
c gauge2
         if(abs(Hex_Ward2).gt.1d-2) then
         gauge2=abs((pen2_for_Ward2-pen1_for_Ward2)/(3d0*Hex_Ward2)-1d0)
         else
         gauge2=abs((pen2_for_Ward2-pen1_for_Ward2) - (3d0*Hex_Ward2))
         endif
c gauge3
         if(abs(Hex_Ward3).gt.1d-2) then
         gauge3=abs((pen3_for_Ward3-pen2_for_Ward3)/(3d0*Hex_Ward3)-1d0)
         else
         gauge3=abs((pen3_for_Ward3-pen2_for_Ward3) - (3d0*Hex_Ward3))
         endif
c gauge4
         if(abs(Hex_Ward4).gt.1d-2) then
         gauge4=abs((pen4_for_Ward4-pen3_for_Ward4)/(3d0*Hex_Ward4)-1d0)
         else
         gauge4=abs((pen4_for_Ward4-pen3_for_Ward4) - (3d0*Hex_Ward4))
         endif
c gauge5
         if(abs(Hex_Ward5).gt.1d-2) then
         gauge5=abs((pen5_for_Ward5-pen4_for_Ward5)/(3d0*Hex_Ward5)-1d0)
         else
         gauge5=abs((pen5_for_Ward5-pen4_for_Ward5) - (3d0*Hex_Ward5))
         endif
cfc         print*, "aux(jmass)", aux(jmass)
cfc         print*, "Hex_Ward1",Hex_Ward1
cfc         print*, "Hex_Ward2",Hex_Ward2
cfc         print*, "Hex_Ward3",Hex_Ward3
cfc         print*, "Hex_Ward4",Hex_Ward4
cfc         print*, "Hex_Ward5",Hex_Ward5
cfc         print*, "pen5_for_ward5",pen5_for_ward5/(12d0*massj(jmass)*massj(jmass))  
cfc         print*, "pen4_for_ward5",pen4_for_ward5/(12d0*massj(jmass)*massj(jmass))  
cfc         print*, "hex_ward5",Hex_ward5/(3d0*massj(jmass)*massj(jmass)) 
cfc
cfc         print*, "pen1_for_ward1",pen1_for_ward1/(12d0*massj(jmass)*massj(jmass))  
cfc         print*, "pen5_for_ward1",pen5_for_ward1/(12d0*massj(jmass)*massj(jmass))  
cfc         print*, "hex_ward1",Hex_ward1/(3d0*massj(jmass)*massj(jmass)) 
cfc         print*, "gauge1",gauge1
cfc         print*, "gauge2",gauge2
cfc         print*, "gauge3",gauge3
cfc         print*, "gauge4",gauge4
cfc         print*, "gauge5",gauge5

         if(gauge1+1.eq.gauge1) then
            gauge1=1d99
         endif
         if(gauge1.ne.gauge1) then
            gauge1=1d99
         endif

         if(gauge2+1.eq.gauge2) then
            gauge2=1d99
         endif
         if(gauge2.ne.gauge2) then
            gauge2=1d99
         endif


         if(gauge3+1.eq.gauge3) then
            gauge3=1d99
         endif
         if(gauge3.ne.gauge3) then
            gauge3=1d99
         endif


         if(gauge4+1.eq.gauge4) then
            gauge4=1d99
         endif
         if(gauge4.ne.gauge4) then
            gauge4=1d99
         endif


         if(gauge5+1.eq.gauge5) then
            gauge5=1d99
         endif
         if(gauge5.ne.gauge5) then
            gauge5=1d99
         endif


          gauge=Max(gauge1,gauge2,gauge3,gauge4,gauge5)


         if(gauge+1.eq.gauge) then
            gauge=1d99
         endif
         if(gauge.ne.gauge) then
            gauge=1d99
         endif


          if (gauge.lt.gaugelimit) then
c             if(.not.QUAD) then
c             print*,"gauge:true:H",aux(jmass),gauge, HexSingular
c              endif
             else
                HexSingular(2)=.true.
                aux(jmass)=0d0
             print*,"gauge:false:HQ",aux(jmass),gauge, HexSingular(2)
          endif




          endif ! QUAD

#endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c           END QUAD PRECISION
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         stop
c
          enddo
      

c
C
C
C
C
C
C
C
C
C
C
C
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
C
C
C
C               CP ODD CASE
C
C
C
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C

       shift=2

      if (AloopTop) then
       initdo=1
       else
       initdo=2
      endif
      if(AloopBottom) then
         finaldo=2
      else
        finaldo=1
      endif

      massj(shift+1)=m_top
      massj(shift+2)=m_bottom
      cplmassj(shift+1)=cplModAtop
      cplmassj(shift+2)=cplModAbottom
      

      do jmass=shift+initdo,shift+finaldo


         if(Dble1) then

c         CALL gluonHexHiggsOdd_optimize(massj(jmass),K1(0,i1),K1(0,i2),K1(0,i3),K1(0,i4),
c     &        K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i2),epsgC(0,i3),
c     &        epsgC(0,i4),epsgC(0,i5),1,aux(jmass))

            
         CALL gluonHexHiggsOdd_optimize(massj(jmass),K1(0,i1),K1(0,i2),K1(0,i3),K1(0,i4),
     &        K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i2),epsgC(0,i3),
     &        epsgC(0,i4),epsgC(0,i5),1,qhexH_top)

cFC         if (abs(AUX(jmass)/qhexH_top-1d0).gt.1d-10)then
cFC         print*, "zero", AUX(jmass)/qhexH_top-1d0
cFC         endif
c         print*,"hex:mass", massj(jmass),cplmassj(jmass)/massj(jmass)**2

        
         aux(jmass)=-cplmassj(jmass)*qhexH_top/(4d0*massj(jmass)*massj(jmass))  

c     print*, "aux", aux(1),aux(2),aux(3),aux(jmass)
c         print*, "aux", aux(jmass)
c         pause

cc start check
c         aux(jmass)=-cplModHtop*qhexH_top/(12d0*massj(jmass)*massj(jmass))  
cc end check
c start Ward identity
c The momenta of the Hexline is going to be the same. 
c                       --->  comp=-1 to avoid redundant calculation
c
c
c First Gauge test Eps_i -->   p_i
c +  auxiliary momenta
c         print*, "i1--i5", i1,i2,i3,i4,i5
         do jtemp=0,3
c            print*, "epsgC(0,i1)",epsgC(jtemp,i1)
            mup1temp(jtemp)= K1(jtemp,i1)
            mup2temp(jtemp)= K1(jtemp,i2)
            mup3temp(jtemp)= K1(jtemp,i3)
            mup4temp(jtemp)= K1(jtemp,i4)
            mup5temp(jtemp)= K1(jtemp,i5)
            pp12(jtemp)=K1(jtemp,i1)+K1(jtemp,i2)
            pp23(jtemp)=K1(jtemp,i2)+K1(jtemp,i3)
            pp34(jtemp)=K1(jtemp,i3)+K1(jtemp,i4)
            pp45(jtemp)=K1(jtemp,i4)+K1(jtemp,i5)
            pp56(jtemp)=-(K1(jtemp,i1)+K1(jtemp,i2)+K1(jtemp,i3)+K1(jtemp,i4))
            pp16(jtemp)=-(K1(jtemp,i2)+K1(jtemp,i3)+K1(jtemp,i4)+K1(jtemp,i5))
           enddo   
c   Ward 1
         CALL gluonHexHiggsOdd_optimize(massj(jmass),K1(0,i1),K1(0,i2),K1(0,i3),K1(0,i4),
     &        K1(0,i5),PH1,mup1temp,epsgC(0:3,i2),epsgC(0:3,i3),
     &        epsgC(0:3,i4),epsgC(0:3,i5),-1,Hex_Ward1)
c   Ward 2
         CALL gluonHexHiggsOdd_optimize(massj(jmass),K1(0,i1),K1(0,i2),K1(0,i3),K1(0,i4),
     &        K1(0,i5),PH1,epsgC(0,i1),mup2temp,epsgC(0,i3),
     &        epsgC(0,i4),epsgC(0,i5),-1,Hex_Ward2)
c   Ward 3
         CALL gluonHexHiggsOdd_optimize(massj(jmass),K1(0,i1),K1(0,i2),K1(0,i3),K1(0,i4),
     &        K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i2),mup3temp,
     &        epsgC(0,i4),epsgC(0,i5),-1,Hex_Ward3)
c   Ward 4
         CALL gluonHexHiggsOdd_optimize(massj(jmass),K1(0,i1),K1(0,i2),K1(0,i3),K1(0,i4),
     &        K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i2),epsgC(0,i3),
     &        mup4temp,epsgC(0,i5),-1,Hex_Ward4)
c   Ward 5
         CALL gluonHexHiggsOdd_optimize(massj(jmass),K1(0,i1),K1(0,i2),K1(0,i3),K1(0,i4),
     &        K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i2),epsgC(0,i3),
     &        epsgC(0,i4),mup5temp,-1,Hex_Ward5)

c Call of pentagons
c  Each pentagon will appear twice. eg., pen1:  pen2 -pen1,  pen1 -pen0
c wich different sign depending for which contraction is used.
c      print*, "epsgC(0,i5)_2",epsgC(0,i5),epsgC(1,i5),epsgC(2,i5),epsgC(3,i5)
c  pp12
         CALL gluonpenHiggsOdd(massj(jmass),pp12,K1(0,i3),
     &        K1(0,i4),K1(0,i5),PH1,epsgC(0,i2),epsgC(0,i3),
     &        epsgC(0,i4),epsgC(0,i5),pen1_for_Ward1,1)

         CALL gluonpenHiggsOdd(massj(jmass),pp12,K1(0,i3),
     &        K1(0,i4),K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i3),
     &        epsgC(0,i4),epsgC(0,i5),pen1_for_Ward2,-1)

c  pp23
         CALL gluonpenHiggsOdd(massj(jmass),K1(0,i1),pp23,
     &        K1(0,i4),K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i3),
     &        epsgC(0,i4),epsgC(0,i5),pen2_for_Ward2,1)

         CALL gluonpenHiggsOdd(massj(jmass),K1(0,i1),pp23,
     &        K1(0,i4),K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i2),
     &        epsgC(0,i4),epsgC(0,i5),pen2_for_Ward3,-1)

c pp34
         CALL gluonpenHiggsOdd(massj(jmass),K1(0,i1),K1(0,i2),
     &        pp34,K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i2),
     &        epsgC(0,i4),epsgC(0,i5),pen3_for_Ward3,1)

         CALL gluonpenHiggsOdd(massj(jmass),K1(0,i1),K1(0,i2),
     &        pp34,K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i2),
     &        epsgC(0,i3),epsgC(0,i5),pen3_for_Ward4,-1)

c pp45
         CALL gluonpenHiggsOdd(massj(jmass),K1(0,i1),K1(0,i2),
     &        K1(0,i3),pp45,PH1,epsgC(0,i1),epsgC(0,i2),
     &        epsgC(0,i3),epsgC(0,i5),pen4_for_Ward4,1)

         CALL gluonpenHiggsOdd(massj(jmass),K1(0,i1),K1(0,i2),
     &        K1(0,i3),pp45,PH1,epsgC(0,i1),epsgC(0,i2),
     &        epsgC(0,i3),epsgC(0,i4),pen4_for_Ward5,-1)

c pp56
         CALL gluonpenHiggsOdd(massj(jmass),K1(0,i1),K1(0,i2),
     &        K1(0,i3),K1(0,i4),pp56,epsgC(0,i1),epsgC(0,i2),
     &        epsgC(0,i3),epsgC(0,i4),pen5_for_Ward5,1)

c         CALL gluonpenHiggsOdd(massj(jmass),K1(0,i2),K1(0,i3),
c     &        K1(0,i4),K1(0,i5),pp56,epsgC(0,i2),
c     &        epsgC(0,i3),epsgC(0,i4),epsgC(0,i5),pen5_for_Ward1)
c      print*, "epsgC(0,i5)_3",epsgC(0,i5),epsgC(1,i5),epsgC(2,i5),epsgC(3,i5)
c pp16 
         CALL gluonpenHiggsOdd(massj(jmass),K1(0,i2),K1(0,i3),
     &        K1(0,i4),K1(0,i5),pp16,epsgC(0,i2),epsgC(0,i3),
     &        epsgC(0,i4),epsgC(0,i5),pen5_for_Ward1,1)

c      print*, "epsgC(0,i5)",epsgC(0,i5),epsgC(1,i5),epsgC(2,i5),epsgC(3,i5)
c gauge1
         if(abs(Hex_Ward1).gt.1d-2) then
         gauge1=abs((pen1_for_Ward1-pen5_for_Ward1)/(3d0*Hex_Ward1)-1d0)
         else
         gauge1=abs((pen1_for_Ward1-pen5_for_Ward1) - (3d0*Hex_Ward1))
         endif
c gauge2
         if(abs(Hex_Ward2).gt.1d-2) then
         gauge2=abs((pen2_for_Ward2-pen1_for_Ward2)/(3d0*Hex_Ward2)-1d0)
         else
         gauge2=abs((pen2_for_Ward2-pen1_for_Ward2) - (3d0*Hex_Ward2))
         endif
c gauge3
         if(abs(Hex_Ward3).gt.1d-2) then
         gauge3=abs((pen3_for_Ward3-pen2_for_Ward3)/(3d0*Hex_Ward3)-1d0)
         else
         gauge3=abs((pen3_for_Ward3-pen2_for_Ward3) - (3d0*Hex_Ward3))
         endif
c gauge4
         if(abs(Hex_Ward4).gt.1d-2) then
         gauge4=abs((pen4_for_Ward4-pen3_for_Ward4)/(3d0*Hex_Ward4)-1d0)
         else
         gauge4=abs((pen4_for_Ward4-pen3_for_Ward4) - (3d0*Hex_Ward4))
         endif
c gauge5
         if(abs(Hex_Ward5).gt.1d-2) then
         gauge5=abs((pen5_for_Ward5-pen4_for_Ward5)/(3d0*Hex_Ward5)-1d0)
         else
         gauge5=abs((pen5_for_Ward5-pen4_for_Ward5) - (3d0*Hex_Ward5))
         endif

cfc         print*, "pen5_for_ward5",pen5_for_ward5/(12d0*massj(jmass)*massj(jmass))  
cfc         print*, "pen4_for_ward5",pen4_for_ward5/(12d0*massj(jmass)*massj(jmass))  
cfc         print*, "hex_ward5",Hex_ward5/(3d0*massj(jmass)*massj(jmass)) 
cfc
cfc         print*, "pen1_for_ward1",pen1_for_ward1/(12d0*massj(jmass)*massj(jmass))  
cfc         print*, "pen5_for_ward1",pen5_for_ward1/(12d0*massj(jmass)*massj(jmass))  
cfc         print*, "hex_ward1",Hex_ward1/(3d0*massj(jmass)*massj(jmass)) 
CFC           print*, "gauge1",gauge1,(pen1_for_Ward1-pen5_for_Ward1)/(3d0*Hex_Ward1)
cfc         print*, "gauge2",gauge2
cfc         print*, "gauge3",gauge3
cfc         print*, "gauge4",gauge4
cfc         print*, "gauge5",gauge5
cfc         print*, "aux(jmass)", aux(jmass)
cfc         print*, "Hex_Ward1",Hex_Ward1
cfc         print*, "Hex_Ward2",Hex_Ward2
cfc         print*, "Hex_Ward3",Hex_Ward3
cfc         print*, "Hex_Ward4",Hex_Ward4
cfc         print*, "Hex_Ward5",Hex_Ward5
 

check for NANS


         if(gauge1+1.eq.gauge1) then
            gauge1=1d99
         endif
         if(gauge1.ne.gauge1) then
            gauge1=1d99
         endif

         if(gauge2+1.eq.gauge2) then
            gauge2=1d99
         endif
         if(gauge2.ne.gauge2) then
            gauge2=1d99
         endif


         if(gauge3+1.eq.gauge3) then
            gauge3=1d99
         endif
         if(gauge3.ne.gauge3) then
            gauge3=1d99
         endif


         if(gauge4+1.eq.gauge4) then
            gauge4=1d99
         endif
         if(gauge4.ne.gauge4) then
            gauge4=1d99
         endif


         if(gauge5+1.eq.gauge5) then
            gauge5=1d99
         endif
         if(gauge5.ne.gauge5) then
            gauge5=1d99
         endif

          gauge=Max(gauge1,gauge2,gauge3,gauge4,gauge5)

c           print*, "gauge:Ho",gauge

         if(gauge+1.eq.gauge) then
            gauge=1d99
         endif
         if(gauge.ne.gauge) then
            gauge=1d99
         endif



          if (gauge.lt.gaugelimit) then
c             print*,"gauge:true:H",aux(jmass),gauge, HexSingular
             else
                HexSingular(1)=.true.
                aux(jmass)=0d0
c             print*,"gauge:false:Hd",aux(jmass),gauge, HexSingular(1)
          endif
          endif ! dble

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    START QUAD precision
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef WITH_QUAD          
          If(QUAD.or.hexSingular(1)) then
        
         comp=1

         CALL gluonHexHiggsOdd_QUAD_optimize(massj(jmass),K1(0,i1),K1(0,i2),K1(0,i3),K1(0,i4),
     &        K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i2),epsgC(0,i3),
     &        epsgC(0,i4),epsgC(0,i5),comp,qhexH_top)
         aux(jmass)=-cplmassj(jmass)*qhexH_top/(4d0*massj(jmass)*massj(jmass))  
     
c start Ward identity
c The momenta of the Hexline is going to be the same. 
c                       --->  comp=-1 to avoid redundant calculation
c First Gauge test Eps_i -->   p_i
c +  auxiliary momenta
c         print*, "i1--i5", i1,i2,i3,i4,i5
         do jtemp=0,3
c            print*, "epsgC(0,i1)",epsgC(jtemp,i1)
            mup1temp(jtemp)= K1(jtemp,i1)
            mup2temp(jtemp)= K1(jtemp,i2)
            mup3temp(jtemp)= K1(jtemp,i3)
            mup4temp(jtemp)= K1(jtemp,i4)
            mup5temp(jtemp)= K1(jtemp,i5)
            pp12(jtemp)=K1(jtemp,i1)+K1(jtemp,i2)
            pp23(jtemp)=K1(jtemp,i2)+K1(jtemp,i3)
            pp34(jtemp)=K1(jtemp,i3)+K1(jtemp,i4)
            pp45(jtemp)=K1(jtemp,i4)+K1(jtemp,i5)
            pp56(jtemp)=-(K1(jtemp,i1)+K1(jtemp,i2)+K1(jtemp,i3)+K1(jtemp,i4))
            pp16(jtemp)=-(K1(jtemp,i2)+K1(jtemp,i3)+K1(jtemp,i4)+K1(jtemp,i5))
           enddo   

c         print*, "epsgC(0,i5)_1",epsgC(0,i5),epsgC(1,i5),epsgC(2,i5),epsgC(3,i5)

c   Ward 1
         CALL gluonHexHiggsOdd_QUAD_optimize(massj(jmass),K1(0,i1),K1(0,i2),K1(0,i3),K1(0,i4),
     &        K1(0,i5),PH1,mup1temp,epsgC(0:3,i2),epsgC(0:3,i3),
     &        epsgC(0:3,i4),epsgC(0:3,i5),-1,Hex_Ward1)
c         print*, "Hex_Ward1/res",Hex_Ward1/qhexH_top/K1(0,i1)
c   Ward 2
         CALL gluonHexHiggsOdd_QUAD_optimize(massj(jmass),K1(0,i1),K1(0,i2),K1(0,i3),K1(0,i4),
     &        K1(0,i5),PH1,epsgC(0,i1),mup2temp,epsgC(0,i3),
     &        epsgC(0,i4),epsgC(0,i5),-1,Hex_Ward2)
c   Ward 3
         CALL gluonHexHiggsOdd_QUAD_optimize(massj(jmass),K1(0,i1),K1(0,i2),K1(0,i3),K1(0,i4),
     &        K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i2),mup3temp,
     &        epsgC(0,i4),epsgC(0,i5),-1,Hex_Ward3)
c   Ward 4
         CALL gluonHexHiggsOdd_QUAD_optimize(massj(jmass),K1(0,i1),K1(0,i2),K1(0,i3),K1(0,i4),
     &        K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i2),epsgC(0,i3),
     &        mup4temp,epsgC(0,i5),-1,Hex_Ward4)
c   Ward 5
         CALL gluonHexHiggsOdd_QUAD_optimize(massj(jmass),K1(0,i1),K1(0,i2),K1(0,i3),K1(0,i4),
     &        K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i2),epsgC(0,i3),
     &        epsgC(0,i4),mup5temp,-1,Hex_Ward5)

c Call of pentagons
c  Each pentagon will appear twice. eg., pen1:  pen2 -pen1,  pen1 -pen0
c wich different sign depending for which contraction is used.
c      print*, "epsgC(0,i5)_2",epsgC(0,i5),epsgC(1,i5),epsgC(2,i5),epsgC(3,i5)
c  pp12
         CALL gluonpenHiggsOdd_QUAD(massj(jmass),pp12,K1(0,i3),
     &        K1(0,i4),K1(0,i5),PH1,epsgC(0,i2),epsgC(0,i3),
     &        epsgC(0,i4),epsgC(0,i5),pen1_for_Ward1,1)

         CALL gluonpenHiggsOdd_QUAD(massj(jmass),pp12,K1(0,i3),
     &        K1(0,i4),K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i3),
     &        epsgC(0,i4),epsgC(0,i5),pen1_for_Ward2,-1)

c  pp23
         CALL gluonpenHiggsOdd_QUAD(massj(jmass),K1(0,i1),pp23,
     &        K1(0,i4),K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i3),
     &        epsgC(0,i4),epsgC(0,i5),pen2_for_Ward2,1)

         CALL gluonpenHiggsOdd_QUAD(massj(jmass),K1(0,i1),pp23,
     &        K1(0,i4),K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i2),
     &        epsgC(0,i4),epsgC(0,i5),pen2_for_Ward3,-1)

c pp34
         CALL gluonpenHiggsOdd_QUAD(massj(jmass),K1(0,i1),K1(0,i2),
     &        pp34,K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i2),
     &        epsgC(0,i4),epsgC(0,i5),pen3_for_Ward3,1)

         CALL gluonpenHiggsOdd_QUAD(massj(jmass),K1(0,i1),K1(0,i2),
     &        pp34,K1(0,i5),PH1,epsgC(0,i1),epsgC(0,i2),
     &        epsgC(0,i3),epsgC(0,i5),pen3_for_Ward4,-1)

c pp45
         CALL gluonpenHiggsOdd_QUAD(massj(jmass),K1(0,i1),K1(0,i2),
     &        K1(0,i3),pp45,PH1,epsgC(0,i1),epsgC(0,i2),
     &        epsgC(0,i3),epsgC(0,i5),pen4_for_Ward4,1)

         CALL gluonpenHiggsOdd_QUAD(massj(jmass),K1(0,i1),K1(0,i2),
     &        K1(0,i3),pp45,PH1,epsgC(0,i1),epsgC(0,i2),
     &        epsgC(0,i3),epsgC(0,i4),pen4_for_Ward5,-1)

c pp56
         CALL gluonpenHiggsOdd_QUAD(massj(jmass),K1(0,i1),K1(0,i2),
     &        K1(0,i3),K1(0,i4),pp56,epsgC(0,i1),epsgC(0,i2),
     &        epsgC(0,i3),epsgC(0,i4),pen5_for_Ward5,1)

c         CALL gluonpenHiggsOdd(massj(jmass),K1(0,i2),K1(0,i3),
c     &        K1(0,i4),K1(0,i5),pp56,epsgC(0,i2),
c     &        epsgC(0,i3),epsgC(0,i4),epsgC(0,i5),pen5_for_Ward1)
c      print*, "epsgC(0,i5)_3",epsgC(0,i5),epsgC(1,i5),epsgC(2,i5),epsgC(3,i5)
c pp16 
         CALL gluonpenHiggsOdd_QUAD(massj(jmass),K1(0,i2),K1(0,i3),
     &        K1(0,i4),K1(0,i5),pp16,epsgC(0,i2),epsgC(0,i3),
     &        epsgC(0,i4),epsgC(0,i5),pen5_for_Ward1,1)

c      print*, "epsgC(0,i5)",epsgC(0,i5),epsgC(1,i5),epsgC(2,i5),epsgC(3,i5)
c gauge1
         if(abs(Hex_Ward1).gt.1d-2) then
         gauge1=abs((pen1_for_Ward1-pen5_for_Ward1)/(3d0*Hex_Ward1)-1d0)
         else
         gauge1=abs((pen1_for_Ward1-pen5_for_Ward1) - (3d0*Hex_Ward1))
         endif
c gauge2
         if(abs(Hex_Ward2).gt.1d-2) then
         gauge2=abs((pen2_for_Ward2-pen1_for_Ward2)/(3d0*Hex_Ward2)-1d0)
         else
         gauge2=abs((pen2_for_Ward2-pen1_for_Ward2) - (3d0*Hex_Ward2))
         endif
c gauge3
         if(abs(Hex_Ward3).gt.1d-2) then
         gauge3=abs((pen3_for_Ward3-pen2_for_Ward3)/(3d0*Hex_Ward3)-1d0)
         else
         gauge3=abs((pen3_for_Ward3-pen2_for_Ward3) - (3d0*Hex_Ward3))
         endif
c gauge4
         if(abs(Hex_Ward4).gt.1d-2) then
         gauge4=abs((pen4_for_Ward4-pen3_for_Ward4)/(3d0*Hex_Ward4)-1d0)
         else
         gauge4=abs((pen4_for_Ward4-pen3_for_Ward4) - (3d0*Hex_Ward4))
         endif
c gauge5
         if(abs(Hex_Ward5).gt.1d-2) then
         gauge5=abs((pen5_for_Ward5-pen4_for_Ward5)/(3d0*Hex_Ward5)-1d0)
         else
         gauge5=abs((pen5_for_Ward5-pen4_for_Ward5) - (3d0*Hex_Ward5))
         endif
cfc         print*, "aux(jmass)", aux(jmass)
cfc         print*, "Hex_Ward1",Hex_Ward1
cfc         print*, "Hex_Ward2",Hex_Ward2
cfc         print*, "Hex_Ward3",Hex_Ward3
cfc         print*, "Hex_Ward4",Hex_Ward4
cfc         print*, "Hex_Ward5",Hex_Ward5
cfc         print*, "pen5_for_ward5",pen5_for_ward5/(12d0*massj(jmass)*massj(jmass))  
cfc         print*, "pen4_for_ward5",pen4_for_ward5/(12d0*massj(jmass)*massj(jmass))  
cfc         print*, "hex_ward5",Hex_ward5/(3d0*massj(jmass)*massj(jmass)) 
cfc
cfc         print*, "pen1_for_ward1",pen1_for_ward1/(12d0*massj(jmass)*massj(jmass))  
cfc         print*, "pen5_for_ward1",pen5_for_ward1/(12d0*massj(jmass)*massj(jmass))  
cfc         print*, "hex_ward1",Hex_ward1/(3d0*massj(jmass)*massj(jmass)) 
c            print*, "gauge1:quad",gauge1
cfc         print*, "gauge2",gauge2
cfc         print*, "gauge3",gauge3
cfc         print*, "gauge4",gauge4
cfc         print*, "gauge5",gauge5


         if(gauge1+1.eq.gauge1) then
            gauge1=1d99
         endif
         if(gauge1.ne.gauge1) then
            gauge1=1d99
         endif

         if(gauge2+1.eq.gauge2) then
            gauge2=1d99
         endif
         if(gauge2.ne.gauge2) then
            gauge2=1d99
         endif


         if(gauge3+1.eq.gauge3) then
            gauge3=1d99
         endif
         if(gauge3.ne.gauge3) then
            gauge3=1d99
         endif


         if(gauge4+1.eq.gauge4) then
            gauge4=1d99
         endif
         if(gauge4.ne.gauge4) then
            gauge4=1d99
         endif


         if(gauge5+1.eq.gauge5) then
            gauge5=1d99
         endif
         if(gauge5.ne.gauge5) then
            gauge5=1d99
         endif


          gauge=Max(gauge1,gauge2,gauge3,gauge4,gauge5)

c           print*, "gauge:Ho:quad",gauge
c           pause
         if(gauge+1.eq.gauge) then
            gauge=1d99
         endif

         if(gauge.ne.gauge) then
            gauge=1d99
         endif



          if (gauge.lt.gaugelimit) then
c             if(.not.QUAD) then
c             print*,"gauge:true:H",aux(jmass),gauge, HexSingular
c              endif
             else
             HexSingular(2)=.true.
             print*,"gauge:false:HQ",aux(jmass),gauge, HexSingular(2)
             aux(jmass)= 0d0
          endif


          endif ! QUAD

#endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c           END QUAD PRECISION
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         stop
c
          enddo
      


      do i=1,26
         FDHexLoop=FDHexLoop+aux(i)
      enddo
      
      RETURN
      END
