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
              print*,"gauge:false:PQ",aux(jmass),gauge, PenSingular(2)
          endif



          endif ! QUAD
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
              print*,"gauge:false:PQ",aux(jmass),gauge, PenSingular(2)
          endif



          endif ! QUAD
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
c             print*,"gauge:false:Hd",aux(jmass),gauge, HexSingular(1)
          endif
          endif ! dble

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    START QUAD precision
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
             print*,"gauge:false:HQ",aux(jmass),gauge, HexSingular(2)
          endif




          endif ! QUAD


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
c             print*,"gauge:false:Hd",aux(jmass),gauge, HexSingular(1)
          endif
          endif ! dble

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    START QUAD precision
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
          endif


          endif ! QUAD


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






      
      
*****************
      SUBROUTINE Higgs3jGG(KBAR,GPOL,GSIGN,RHsum)
*----------------------------------------------------------------------
c...   author  : Michael Kubocz
c...   date    : 15.04.2009
c...   modified: 26.03.2012
c
c...  modified: Francisco Campanario 
c...  add:  Ward identities for pentagons and hexagons
c...        + control of instabilities in Dble and QUAD
***********************************************************************
      IMPLICIT NONE

#include "ggf_process.inc"
#include "ggf_dotproducts.inc"
#include "ggfloMG.inc" 

      logical SingularC(2)
      common/singluonver/SingularC

      logical eva
      real*8 phc(0:3)
      common/phc1/phc,eva

c... information on coupling constants from KOPPELN
      real*8 ALFAS,XMT,ALFA,XMZ,XMW,SIN2W,XMH
      real*8 clr, xm2, xmg, b, v, a
      COMMON /BKOPIN/ ALFAS,XMT,ALFA,XMZ,XMW,SIN2W,XMH
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)

      real*8 KBAR(0:4,5),PHBAR(0:4),epsgMG(0:3,2,5)

      real*8 K(0:4,5),Kij(0:4,5,5),Kijk(0:4,5,5,5),PH(0:4),KPH(0:4,5),
     &     epsg(0:3,5),
     &     Kpen(0:4,4),KPHeffij(0:4,5,5),KPAeffij(0:4,5,5),
     &     KPHtopij(0:4,5,5),KPHbottomij(0:4,5,5),KPAtopij(0:4,5,5),
     &     KPAbottomij(0:4,5,5)
c... 4-dim momentum variables
      real*8 Kp(0:3,5),Kijp(0:3,5,5),Kijkp(0:3,5,5,5),PHp(0:3),
     &     KPHp(0:3,5),KPHeffijp(0:3,5,5),KPAeffijp(0:3,5,5),K123p(0:3),
     &     KPHtopijp(0:3,5,5),KPHbottomijp(0:3,5,5),KPAtopijp(0:3,5,5),
     &     KPAbottomijp(0:3,5,5)

      complex*16 epsgC(0:3,5),Jij(0:3,5,5),Jijk(0:3,5,5,5),
     &           epsgPen(0:3,4)
      complex*16 epsHeff(0:3,5),epsAeff(0:3,5),epsHtop(0:3,5),
     &           epsHbottom(0:3,5),epsAtop(0:3,5),epsAbottom(0:3,5),
     &           JPHeffij(0:3,5,5),JPAeffij(0:3,5,5),JPHtopij(0:3,5,5),
     &           JPHbottomij(0:3,5,5),JPAtopij(0:3,5,5),JPAbottomij(0:3,5,5)
      COMMON /GG3JET/ K,Kij,Kijk,epsg,epsgC,Jij,Jijk,KPH,PH,
     &        KPHeffij,KPAeffij,KPHtopij,KPHbottomij,KPAtopij,
     &        KPAbottomij,
     &        epsHeff,epsAeff,epsHtop,epsHbottom,epsAtop,epsAbottom,
     &        JPHeffij,JPAeffij,JPHtopij,JPHbottomij,JPAtopij,
     &        JPAbottomij,
c... 4-dim momentum variables
     &        Kp,Kijp,Kijkp,PHp,KPHp,KPHeffijp,KPAeffijp,K123p,
     &        KPHtopijp,KPHbottomijp,KPAtopijp,KPAbottomijp

      Common/Forpent/Kpen,epsgPen
      
      complex*16 epsChelas(1:6,-1:1,5)

      integer GPOL(5),GSIGN(5),HSIGN,ipol

c... different momenta combinations
      real*8 p21(0:4),p21K1(0:4),p21K2(0:4),p21K3(0:4),
     &       p21K12(0:4),p21K13(0:4),p21K23(0:4),
     &       p21K1H(0:4),p21K2H(0:4),p21K3H(0:4),
     &       K12(0:4),K13(0:4),K23(0:4),K123(0:4),
     &       PHK1(0:4),PHK2(0:4),PHK3(0:4),PHp21(0:4),
     &       PHK13(0:4),PHK23(0:4),PHK12(0:4)

c... auxiliary "current" variables
      complex*16 GGG12(0:4,2,2),GGG13(0:4,2,2),GGG23(0:4,2,2),
     &           GGG21(0:4,2,2),GGG31(0:4,2,2),GGG32(0:4,2,2),
     &           GGGj211(0:4,2),GGG2j21(0:4,2),GGG3j21(0:4,2),
     &           GGG123(0:4,2,2,2),
     &           GGG132a(0:4,2,2,2),GGG132b(0:4,2,2,2)

      complex*16 GGGPH12(0:4,2,2),GGGPH13(0:4,2,2),GGGPH23(0:4,2,2),
     &           GGGPH32(0:4,2,2),
     &           GGG1PH2(0:4,2,2),GGG1PH3(0:4,2,2)

c... effective vertices and loop functions    
      complex*16 FD3gTri,FD1ePV3g,FD2gBox,FDePV4g,FD4gBox,FDeffPen,
     &           FDeffPenLoop,FDHexLoop
      complex*16 H2G_Loop,A2G_Loop,H3G_Loop,A3G_Loop
      complex*16 H1G_SqLoop,H2G_SqLoop,H3G_SqLoop
      complex*16 H2G_effective,A2G_effective,H3G_effective,A3G_effective,
     &           GGG_vertex,H4G_effective2

      complex*16 pent(15,15,15,15),
     &           qhex12345,qhex12354,qhex12453,qhex12543,qhex13425,
     &           qhex13254,qhex13524,qhex13245,qhex12435,qhex12534,
     &           qhex14325,qhex14235

c... further local variables
      complex*16 amp3j(390),colorT,color(24),
     &     colorPentagons(24),colorHexagons(24)
     &    ,colorHexagons1(24),colorpentagons1(24)

      integer ienter,mu,i,j,l
      real*8 pi,fac,vev,dummy(0:4),RH,RHsum
      parameter (pi=3.14159 26535 89793d0)
      save fac
      data ienter /0/
      save ienter
      real*8 auxmom(0:3,2)
      real*8 ratio,shat
      complex*16 colorN1(24)
      common/colorna/colorN1
      complex*16 Box1,Box2,Pen1,Pen2,Hex

      integer i1,i2,i3,i4,i5

c... external functions
      complex*16 dotrc,dotcc
      real*8 dotrr
      EXTERNAL dotrr,dotrc,dotcc,
     &         H2G_Loop,A2G_Loop,H3G_Loop,A3G_Loop,
     &         H1G_SqLoop,H2G_SqLoop,H3G_SqLoop,
     &         H2G_effective,A2G_effective,H3G_effective,A3G_effective,
     &         GGG_vertex

      real*8 RGB(24,24)
c... color matrix with 24 different color structures (taken from MadGraph)
c... color(1) --> abcde
      DATA (RGB(i,1  ),i=1  ,6  ) /   455,  -58,   14,   68,  -58,   -4/    
      DATA (RGB(i,1  ),i=7  ,12 ) /    68,  -40,   14,   -4,    5,    5/    
      DATA (RGB(i,1  ),i=13 ,18 ) /    14,   -4,  -58,   68,  -58,   14/    
      DATA (RGB(i,1  ),i=19 ,24 ) /    -4,   68,   68,  -58,   14,   -4/    
c... color(2) --> acbde
      DATA (RGB(i,2  ),i=1  ,6  ) /   -58,  455,   68,   14,   -4,  -58/    
      DATA (RGB(i,2  ),i=7  ,12 ) /   -40,   68,    5,  -58,   14,   -4/    
      DATA (RGB(i,2  ),i=13 ,18 ) /    68,    5,   -4,   14,   14,  -58/    
      DATA (RGB(i,2  ),i=19 ,24 ) /    68,   -4,   -4,   14,  -58,   68/    
c... color(3) --> adbce
      DATA (RGB(i,3  ),i=1  ,6  ) /    14,   68,  455,  -58,   68,  -40/    
      DATA (RGB(i,3  ),i=7  ,12 ) /   -58,   -4,   -4,   -4,   68,  -58/    
      DATA (RGB(i,3  ),i=13 ,18 ) /    14,   14,  -58,   68,  -58,   14/    
      DATA (RGB(i,3  ),i=19 ,24 ) /    -4,   68,    5,    5,   -4,   14/    
c... color(4) --> adcbe
      DATA (RGB(i,4  ),i=1  ,6  ) /    68,   14,  -58,  455,  -40,   68/    
      DATA (RGB(i,4  ),i=7  ,12 ) /    -4,  -58,   68,  -58,   -4,   14/    
      DATA (RGB(i,4  ),i=13 ,18 ) /    68,  -58,   -4,   14,   14,  -58/    
      DATA (RGB(i,4  ),i=19 ,24 ) /    68,   -4,   14,   -4,    5,    5/    
c... color(5) --> aebcd
      DATA (RGB(i,5  ),i=1  ,6  ) /   -58,   -4,   68,  -40,  455,  -58/    
      DATA (RGB(i,5  ),i=7  ,12 ) /    14,   68,  -58,   68,   14,   -4/    
      DATA (RGB(i,5  ),i=13 ,18 ) /   -58,   68,   14,   -4,   -4,   68/    
      DATA (RGB(i,5  ),i=19 ,24 ) /   -58,   14,   -4,   14,    5,    5/    
c... color(6) --> aecbd
      DATA (RGB(i,6  ),i=1  ,6  ) /    -4,  -58,  -40,   68,  -58,  455/    
      DATA (RGB(i,6  ),i=7  ,12 ) /    68,   14,   14,   14,  -58,   68/    
      DATA (RGB(i,6  ),i=13 ,18 ) /    -4,   -4,   68,  -58,   68,   -4/    
      DATA (RGB(i,6  ),i=19 ,24 ) /    14,  -58,    5,    5,   14,   -4/    
c... color(7) --> aedbc
      DATA (RGB(i,7  ),i=1  ,6  ) /    68,  -40,  -58,   -4,   14,   68/    
      DATA (RGB(i,7  ),i=7  ,12 ) /   455,  -58,    5,   68,   -4,   14/    
      DATA (RGB(i,7  ),i=13 ,18 ) /   -58,    5,   14,   -4,   -4,   68/    
      DATA (RGB(i,7  ),i=19 ,24 ) /   -58,   14,   14,   -4,   68,  -58/    
c... color(8) --> aedcb
      DATA (RGB(i,8  ),i=1  ,6  ) /   -40,   68,   -4,  -58,   68,   14/    
      DATA (RGB(i,8  ),i=7  ,12 ) /   -58,  455,   -4,   14,    5,    5/    
      DATA (RGB(i,8  ),i=13 ,18 ) /    -4,   14,   68,  -58,   68,   -4/    
      DATA (RGB(i,8  ),i=19 ,24 ) /    14,  -58,  -58,   68,   -4,   14/    
c... color(9) --> abecd
      DATA (RGB(i,9  ),i=1  ,6  ) /    14,    5,   -4,   68,  -58,   14/    
      DATA (RGB(i,9  ),i=7  ,12 ) /     5,   -4,  455,   14,   68,  -58/    
      DATA (RGB(i,9  ),i=13 ,18 ) /    -4,  -40,  -58,   68,   68,   14/    
      DATA (RGB(i,9  ),i=19 ,24 ) /    -4,  -58,  -58,   68,   14,   -4/    
c... color(10) --> acbed
      DATA (RGB(i,10 ),i=1  ,6  ) /    -4,  -58,   -4,  -58,   68,   14/    
      DATA (RGB(i,10 ),i=7  ,12 ) /    68,   14,   14,  455,  -58,   68/    
      DATA (RGB(i,10 ),i=13 ,18 ) /   -40,   -4,  -58,   68,    5,   14/    
      DATA (RGB(i,10 ),i=19 ,24 ) /    -4,    5,  -58,   68,   -4,   14/    
c... color(11) --> acebd
      DATA (RGB(i,11 ),i=1  ,6  ) /     5,   14,   68,   -4,   14,  -58/    
      DATA (RGB(i,11 ),i=7  ,12 ) /    -4,    5,   68,  -58,  455,  -40/    
      DATA (RGB(i,11 ),i=13 ,18 ) /    68,  -58,   14,   -4,   14,   68/    
      DATA (RGB(i,11 ),i=19 ,24 ) /   -58,   -4,   -4,   14,   68,  -58/    
c... color(12) --> adbec
      DATA (RGB(i,12 ),i=1  ,6  ) /     5,   -4,  -58,   14,   -4,   68/    
      DATA (RGB(i,12 ),i=7  ,12 ) /    14,    5,  -58,   68,  -40,  455/    
      DATA (RGB(i,12 ),i=13 ,18 ) /   -58,   68,   -4,   14,   -4,  -58/    
      DATA (RGB(i,12 ),i=19 ,24 ) /    68,   14,   14,   -4,  -58,   68/    
c... color(13) --> adebc
      DATA (RGB(i,13 ),i=1  ,6  ) /    14,   68,   14,   68,  -58,   -4/    
      DATA (RGB(i,13 ),i=7  ,12 ) /   -58,   -4,   -4,  -40,   68,  -58/    
      DATA (RGB(i,13 ),i=13 ,18 ) /   455,   14,   68,  -58,    5,   -4/    
      DATA (RGB(i,13 ),i=19 ,24 ) /    14,    5,   68,  -58,   14,   -4/    
c... color(14) --> adceb
      DATA (RGB(i,14 ),i=1  ,6  ) /    -4,    5,   14,  -58,   68,   -4/    
      DATA (RGB(i,14 ),i=7  ,12 ) /     5,   14,  -40,   -4,  -58,   68/    
      DATA (RGB(i,14 ),i=13 ,18 ) /    14,  455,   68,  -58,  -58,   -4/    
      DATA (RGB(i,14 ),i=19 ,24 ) /    14,   68,   68,  -58,   -4,   14/    
c... color(15) --> abced
      DATA (RGB(i,15 ),i=1  ,6  ) /   -58,   -4,  -58,   -4,   14,   68/    
      DATA (RGB(i,15 ),i=7  ,12 ) /    14,   68,  -58,  -58,   14,   -4/    
      DATA (RGB(i,15 ),i=13 ,18 ) /    68,   68,  455,  -40,   14,    5/    
      DATA (RGB(i,15 ),i=19 ,24 ) /     5,   -4,   14,   -4,   68,  -58/    
c... color(16) --> adecb
      DATA (RGB(i,16 ),i=1  ,6  ) /    68,   14,   68,   14,   -4,  -58/    
      DATA (RGB(i,16 ),i=7  ,12 ) /    -4,  -58,   68,   68,   -4,   14/    
      DATA (RGB(i,16 ),i=13 ,18 ) /   -58,  -58,  -40,  455,   -4,    5/    
      DATA (RGB(i,16 ),i=19 ,24 ) /     5,   14,   -4,   14,  -58,   68/    
c... color(17) --> abdce
      DATA (RGB(i,17 ),i=1  ,6  ) /   -58,   14,  -58,   14,   -4,   68/    
      DATA (RGB(i,17 ),i=7  ,12 ) /    -4,   68,   68,    5,   14,   -4/    
      DATA (RGB(i,17 ),i=13 ,18 ) /     5,  -58,   14,   -4,  455,   68/    
      DATA (RGB(i,17 ),i=19 ,24 ) /   -58,  -40,   14,   -4,  -58,   68/    
c... color(18) --> acdbe
      DATA (RGB(i,18 ),i=1  ,6  ) /    14,  -58,   14,  -58,   68,   -4/    
      DATA (RGB(i,18 ),i=7  ,12 ) /    68,   -4,   14,   14,   68,  -58/    
      DATA (RGB(i,18 ),i=13 ,18 ) /    -4,   -4,    5,    5,   68,  455/    
      DATA (RGB(i,18 ),i=19 ,24 ) /   -40,  -58,   68,  -58,   -4,   14/    
c... color(19) --> aebdc
      DATA (RGB(i,19 ),i=1  ,6  ) /    -4,   68,   -4,   68,  -58,   14/    
      DATA (RGB(i,19 ),i=7  ,12 ) /   -58,   14,   -4,   -4,  -58,   68/    
      DATA (RGB(i,19 ),i=13 ,18 ) /    14,   14,    5,    5,  -58,  -40/    
      DATA (RGB(i,19 ),i=19 ,24 ) /   455,   68,  -58,   68,   14,   -4/    
c... color(20) --> aecdb
      DATA (RGB(i,20 ),i=1  ,6  ) /    68,   -4,   68,   -4,   14,  -58/    
      DATA (RGB(i,20 ),i=7  ,12 ) /    14,  -58,  -58,    5,   -4,   14/    
      DATA (RGB(i,20 ),i=13 ,18 ) /     5,   68,   -4,   14,  -40,  -58/    
      DATA (RGB(i,20 ),i=19 ,24 ) /    68,  455,   -4,   14,   68,  -58/    
c... color(21) --> abedc
      DATA (RGB(i,21 ),i=1  ,6  ) /    68,   -4,    5,   14,   -4,    5/    
      DATA (RGB(i,21 ),i=7  ,12 ) /    14,  -58,  -58,  -58,   -4,   14/    
      DATA (RGB(i,21 ),i=13 ,18 ) /    68,   68,   14,   -4,   14,   68/    
      DATA (RGB(i,21 ),i=19 ,24 ) /   -58,   -4,  455,  -40,  -58,   68/    
c... color(22) --> acdeb
      DATA (RGB(i,22 ),i=1  ,6  ) /   -58,   14,    5,   -4,   14,    5/    
      DATA (RGB(i,22 ),i=7  ,12 ) /    -4,   68,   68,   68,   14,   -4/    
      DATA (RGB(i,22 ),i=13 ,18 ) /   -58,  -58,   -4,   14,   -4,  -58/    
      DATA (RGB(i,22 ),i=19 ,24 ) /    68,   14,  -40,  455,   68,  -58/    
c... color(23) --> abdec
      DATA (RGB(i,23 ),i=1  ,6  ) /    14,  -58,   -4,    5,    5,   14/    
      DATA (RGB(i,23 ),i=7  ,12 ) /    68,   -4,   14,   -4,   68,  -58/    
      DATA (RGB(i,23 ),i=13 ,18 ) /    14,   -4,   68,  -58,  -58,   -4/    
      DATA (RGB(i,23 ),i=19 ,24 ) /    14,   68,  -58,   68,  455,  -40/    
c... color(24) --> acedb
      DATA (RGB(i,24 ),i=1  ,6  ) /    -4,   68,   14,    5,    5,   -4/    
      DATA (RGB(i,24 ),i=7  ,12 ) /   -58,   14,   -4,   14,  -58,   68/    
      DATA (RGB(i,24 ),i=13 ,18 ) /    -4,   14,  -58,   68,   68,   14/    
      DATA (RGB(i,24 ),i=19 ,24 ) /    -4,  -58,   68,  -58,  -40,  455/

      logical SelectABottom,SelectATop,SelectHBottom,SelectHTop
      common/select/SelectABottom,SelectATop,SelectHBottom,SelectHTop

      SelectAbottom=   .false.
      SelectATop=      .false.
      SelectHbottom=   .false.
      SelectHTop=      .false.     

      eva=.false.
      
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
      
      if (ienter.eq.0) then
         ienter = ienter+1
         vev = 2d0*xmw/(v(1,3)*dsqrt(8d0))
c... common overall factor
         fac=64d0*pi*alfas**5/(vev*vev)
      endif
      
      do i=1,5
         do mu = 0,3
            K(mu,i) = -GSIGN(i)*KBAR(mu,i)
            Kp(mu,i)=K(mu,i)
         enddo
      enddo
      do i=1,5
         K(4,i)=0d0
      enddo
      do mu=0,3
        PH(mu)=-(K(mu,1)+K(mu,2)+K(mu,3)+K(mu,4)+K(mu,5))
        PHp(mu)=PH(mu)
      enddo

c... get the gluon polarization vector and the gluon emission spinors
      do i = 1,5
         CALL POLVEC(KBAR(0,i),gpol(i),epsg(0,i))
         if(MGdebug) then
            do mu=0,3
               CALL vxxxxx(KBAR(0,1),0d0,1,-1,epsChelas(1,1,1))
               epsgC(mu,1)=epsChelas(mu+1, 1,1)
               CALL vxxxxx(KBAR(0,2),0d0,1,-1,epsChelas(1,1,2))
               epsgC(mu,2)=epsChelas(mu+1, 1,2)
               CALL vxxxxx(KBAR(0,3),0d0,-1,1,epsChelas(1,1,3))
               epsgC(mu,3)=epsChelas(mu+1, 1,3)
               CALL vxxxxx(KBAR(0,4),0d0,-1,1,epsChelas(1,1,4))
               epsgC(mu,4)=epsChelas(mu+1, 1,4)
               CALL vxxxxx(KBAR(0,5),0d0,-1,1,epsChelas(1,1,5))
               epsgC(mu,5)=epsChelas(mu+1, 1,5)
            enddo
         else
            do mu=0,3
               epsgC(mu,i) = epsg(mu,i)
            enddo
         endif
      enddo
      
      if(lgauge2) then
         i=1
         do mu=0,3
            epsgC(mu,i) = K(mu,i)/K(0,i)
         enddo
      endif
 
c... calculation of gluon "currents"
      do i=1,5
         do j=1,5
c            if(i.eq.j) cycle
            CALL  GGG_current(Kp(0,i),epsgC(0,i),Kp(0,j),epsgC(0,j),
     &                        Kijp(0,i,j),Jij(0,i,j))
c... 4/5-dim momentum variables
            do mu=0,3
               Kij(mu,i,j)=Kijp(mu,i,j)
            enddo
            Kij(4,i,j)=dotrr(Kijp(0,i,j),Kijp(0,i,j))
         enddo
      enddo
      do i=1,5
         do j=1,5
            do l=1,5
               CALL GGG_current(Kijp(0,i,j),Jij(0,i,j),Kp(0,l),
     &                epsgC(0,l),Kijkp(0,i,j,l),Jijk(0,i,j,l))
c... 4/5-dim momentum variables
               do mu=0,3
                  Kijk(mu,i,j,l)=Kijkp(mu,i,j,l)
               enddo
               Kijk(4,i,j,l)=dotrr(Kijkp(0,i,j,l),Kijkp(0,i,j,l))
            enddo
         enddo
      enddo

*********************************************************************
******************** EFFECTIVE CP-EVEN HIGGS ************************
*********************************************************************
c... calculation of effective polarization vectors
      if(Heff) then
         do i=1,5
            CALL EpsH_effective(Kp(0,i),.false.,epsgC(0,i),PHp,KPHp(0,i),
     &           epsHeff(0,i))
c... 4/5-dim momentum variables
            do mu=0,3
               KPH(mu,i)=KPHp(mu,i)
            enddo
         enddo         
c... calculation of gluon "currents" with additional Higgs emission
         do i=1,5
            do j=1,5
               if (i.ne.j) then
                  CALL GGG_current(Kp(0,i),epsgC(0,i),KPHp(0,j),
     &                 epsHeff(0,j),KPHeffijp(0,i,j),JPHeffij(0,i,j))
c... 4/5-dim momentum variables
                  do mu=0,3
                     KPHeffij(mu,i,j)=KPHeffijp(mu,i,j)
                  enddo
                  KPHeffij(4,i,j)=dotrr(KPHeffijp(0,i,j),
     &                                  KPHeffijp(0,i,j))
               else
                  do mu=0,4
                     KPHeffij(mu,i,j)=0d0
                  enddo
               endif
            enddo
         enddo
      endif
*********************************************************************
******************** EFFECTIVE CP-ODD HIGGS *************************
*********************************************************************
c... calculation of effective polarization vectors
      if(Aeff) then
         do i=1,5
            CALL EpsA_effective(Kp(0,i),epsgC(0,i),PHp,KPHp(0,i),
     &           epsAeff(0,i))
c... 4/5-dim momentum variables
            do mu=0,3
               KPH(mu,i)=KPHp(mu,i)
            enddo
         enddo         
c... calculation of gluon "currents" with additional Higgs emission
         do i=1,5
            do j=1,5
               if (i.ne.j) then
                  CALL GGG_current(Kp(0,i),epsgC(0,i),KPHp(0,j),
     &                 epsAeff(0,j),KPAeffijp(0,i,j),JPAeffij(0,i,j))
c... 4/5-dim momentum variables
                  do mu=0,3
                     KPAeffij(mu,i,j)=KPAeffijp(mu,i,j)
                  enddo
                  KPAeffij(4,i,j)=dotrr(KPAeffijp(0,i,j),
     &                                  KPAeffijp(0,i,j))
               else
                  do mu=0,4
                     KPAeffij(mu,i,j)=0d0
                  enddo
               endif
            enddo
         enddo
      endif
*********************************************************************
******************** TOP-LOOP CP-EVEN HIGGS *************************
*********************************************************************
c... calculation of effective polarization vectors
      if(HloopTop) then
         do i=1,5
            CALL EpsH_Loop(K(0,i),.false.,epsgC(0,i),PHp,KPHp(0,i),
     &                     epsHtop(0,i),m_top)
c... 4/5-dim momentum variables
            do mu=0,3
               KPH(mu,i)=KPHp(mu,i)
            enddo
         enddo
c... calculation of gluon "currents" with additional Higgs emission
         do i=1,5
            do j=1,5
               if (i.ne.j) then
                  CALL GGG_current(Kp(0,i),epsgC(0,i),KPHp(0,j),
     &                 epsHtop(0,j),KPHtopijp(0,i,j),JPHtopij(0,i,j))
c... 4/5-dim momentum variables
                  do mu=0,3
                     KPHtopij(mu,i,j)=KPHtopijp(mu,i,j)
                  enddo
                  KPHtopij(4,i,j)=dotrr(KPHtopijp(0,i,j),
     &                                  KPHtopijp(0,i,j))
               else
                  do mu=0,4
                     KPHtopij(mu,i,j)=0d0
                  enddo
               endif
            enddo
         enddo
      endif
*********************************************************************
****************** BOTTOM-LOOP CP-EVEN HIGGS ************************
*********************************************************************
c... calculation of effective polarization vectors
      if(HloopBottom) then
         do i=1,5
            CALL EpsH_Loop(K(0,i),.false.,epsgC(0,i),PHp,KPHp(0,i),
     &                     epsHbottom(0,i),m_bottom)
c... 4/5-dim momentum variables
            do mu=0,3
               KPH(mu,i)=KPHp(mu,i)
            enddo
         enddo
c... calculation of gluon "currents" with additional Higgs emission
         do i=1,5
            do j=1,5
               if (i.ne.j) then
                  CALL GGG_current(Kp(0,i),epsgC(0,i),KPHp(0,j),
     &                 epsHbottom(0,j),KPHbottomijp(0,i,j),
     &                 JPHbottomij(0,i,j))
c... 4/5-dim momentum variables
                  do mu=0,3
                     KPHbottomij(mu,i,j)=KPHbottomijp(mu,i,j)
                  enddo
                  KPHbottomij(4,i,j)=dotrr(KPHbottomijp(0,i,j),
     &                                  KPHbottomijp(0,i,j))
               else
                  do mu=0,4
                     KPHbottomij(mu,i,j)=0d0
                  enddo
               endif
            enddo
         enddo
      endif
*********************************************************************
******************** TOP-LOOP CP-ODD HIGGS **************************
*********************************************************************
c... calculation of effective polarization vectors
      if(AloopTop) then
         do i=1,5
            CALL EpsA_Loop(Kp(0,i),epsgC(0,i),PHp,KPHp(0,i),
     &                     epsAtop(0,i),m_top)
c... 4/5-dim momentum variables
            do mu=0,3
               KPH(mu,i)=KPHp(mu,i)
            enddo
         enddo
c... calculation of gluon "currents" with additional Higgs emission
         do i=1,5
            do j=1,5
               if (i.ne.j) then
                  CALL GGG_current(Kp(0,i),epsgC(0,i),KPHp(0,j),
     &                 epsAtop(0,j),KPAtopijp(0,i,j),JPAtopij(0,i,j))
c... 4/5-dim momentum variables
                  do mu=0,3
                     KPAtopij(mu,i,j)=KPAtopijp(mu,i,j)
                  enddo
                  KPAtopij(4,i,j)=dotrr(KPAtopijp(0,i,j),
     &                                  KPAtopijp(0,i,j))
               else
                  do mu=0,4
                     KPAtopij(mu,i,j)=0d0
                  enddo
               endif
            enddo
         enddo
      endif
      do i=1,5
         KPH(4,i)=dotrr(KPHp(0,i),KPHp(0,i))
      enddo
*********************************************************************
***************** BOTTOM-LOOP CP-ODD HIGGS **************************
*********************************************************************
c... calculation of effective polarization vectors
      if(AloopBottom) then
         do i=1,5
            CALL EpsA_Loop(Kp(0,i),epsgC(0,i),PHp,KPHp(0,i),
     &                     epsAbottom(0,i),m_bottom)
c... 4/5-dim momentum variables
            do mu=0,3
               KPH(mu,i)=KPHp(mu,i)
            enddo
         enddo
c... calculation of gluon "currents" with additional Higgs emission
         do i=1,5
            do j=1,5
               if (i.ne.j) then
                  CALL GGG_current(Kp(0,i),epsgC(0,i),KPHp(0,j),
     &                 epsAbottom(0,j),KPAbottomijp(0,i,j),
     &                 JPAbottomij(0,i,j))
c... 4/5-dim momentum variables
                  do mu=0,3
                     KPAbottomij(mu,i,j)=KPAbottomijp(mu,i,j)
                  enddo
                  KPAbottomij(4,i,j)=dotrr(KPAbottomijp(0,i,j),
     &                                  KPAbottomijp(0,i,j))
               else
                  do mu=0,4
                     KPAbottomij(mu,i,j)=0d0
                  enddo
               endif
            enddo
         enddo
      endif
      do i=1,5
         KPH(4,i)=dotrr(KPHp(0,i),KPHp(0,i))
      enddo
*********************************************************************
*********************************************************************
      amp3j(1)  = FD1ePV3g(5,1,4,2,3,.false.)    ! 1
      amp3j(2)  = FDePV4g(5,1,4,2,3,1)           ! 2
      amp3j(3)  = FDePV4g(5,1,4,2,3,2)           ! 3
      amp3j(4)  = FD1ePV3g(3,2,4,5,1,.false.)    ! 4  
      amp3j(5)  =-FDePV4g(3,2,4,5,1,1)           ! 5
      amp3j(6)  = FDePV4g(3,2,4,5,1,2)           ! 6
      amp3j(7)  = FD1ePV3g(2,3,4,5,1,.false.)    ! 7
      amp3j(8)  =-FDePV4g(2,3,4,5,1,1)           ! 8
      amp3j(9)  = FDePV4g(2,3,4,5,1,2)           ! 9
      amp3j(10) =-FD1ePV3g(5,1,4,2,3,.true.)     ! 10
      amp3j(11) =-FDePV4g(5,1,4,2,3,0)           ! 11
      amp3j(12) =-FDePV4g(2,3,4,5,1,0)           ! 12
      amp3j(13) =-FD1ePV3g(1,5,4,2,3,.false.)    ! 13 
      amp3j(14) = FDePV4g(1,5,4,2,3,1)           ! 14
      amp3j(15) =-FDePV4g(1,5,4,2,3,2)           ! 15
      amp3j(16) = FD2gBox(5,1,4,2,3,.false.)     ! 16
      amp3j(17) = FD4gBox(5,1,4,2,3)             ! 17
      amp3j(18) =-FD2gBox(5,1,4,2,3,.true.)      ! 18
      amp3j(19) = FD2gBox(2,3,4,1,5,.false.)     ! 19
      amp3j(20) = FD4gBox(2,3,4,1,5)             ! 20
      amp3j(21) =-FD3gTri(5,1,4,2,3)             ! 21
      amp3j(22) =-FDePV4g(1,5,4,2,3,3)           ! 22
      amp3j(23) = FDeffPen(5,1,4,2,3)            ! 23
cc      print*,"here:amp3j(23),",amp3j(23)
      amp3j(24) =-FD3gTri(2,3,4,5,1)             ! 24
      amp3j(25) = FDeffPen(2,3,4,5,1)            ! 25
      amp3j(26) = FDePV4g(2,3,4,5,1,3)           ! 26
      amp3j(27) = FD1ePV3g(4,1,3,2,5,.false.)    ! 27
      amp3j(28) = FDePV4g(4,1,3,2,5,1)           ! 28
      amp3j(29) = FDePV4g(4,1,3,2,5,2)           ! 29
      amp3j(30) = FD1ePV3g(5,2,3,4,1,.false.)    ! 30
      amp3j(31) = FDePV4g(5,2,3,4,1,1)           ! 31
      amp3j(32) =-FDePV4g(5,2,3,4,1,2)           ! 32
      amp3j(33) = FD1ePV3g(4,1,3,2,5,.true.)     ! 33
      amp3j(34) = FDePV4g(5,2,3,4,1,0)           ! 34
      amp3j(35) =-FDePV4g(4,1,3,2,5,0)           ! 35
      amp3j(36) =-FD1ePV3g(1,4,3,2,5,.false.)    ! 36
      amp3j(37) = FDePV4g(1,4,3,2,5,1)           ! 37
      amp3j(38) =-FDePV4g(1,4,3,2,5,2)           ! 38
      amp3j(39) =-FD1ePV3g(2,5,3,4,1,.false.)    ! 39
      amp3j(40) = FDePV4g(2,5,3,4,1,1)           ! 40
      amp3j(41) = FDePV4g(2,5,3,4,1,2)           ! 41
      amp3j(42) =-FD2gBox(4,1,3,2,5,.true.)      ! 42
      amp3j(43) = FD2gBox(2,5,3,1,4,.false.)     ! 43
      amp3j(44) =-FD4gBox(2,5,3,4,1)             ! 44
      amp3j(45) =-FD2gBox(4,1,3,5,2,.false.)     ! 45
      amp3j(46) = FD4gBox(4,1,3,5,2)             ! 46
      amp3j(47) =-FD3gTri(2,5,3,4,1)             ! 47
      amp3j(48) = FDePV4g(2,5,3,4,1,3)           ! 48
      amp3j(49) = FDeffPen(2,5,3,4,1)            ! 49   
      amp3j(50) = FD3gTri(4,1,3,2,5)             ! 50
      amp3j(51) =-FDeffPen(4,1,3,2,5)            ! 51
      amp3j(52) = FDePV4g(4,1,3,2,5,3)           ! 52
      amp3j(53) = FD1ePV3g(4,1,2,5,3,.false.)    ! 53
      amp3j(54) = FDePV4g(4,1,2,5,3,1)           ! 54
      amp3j(55) =-FDePV4g(4,1,2,5,3,2)           ! 55
      amp3j(56) =-FD1ePV3g(4,1,2,5,3,.true.)     ! 56
      amp3j(57) = FDePV4g(5,3,2,4,1,0)           ! 57
      amp3j(58) = FDePV4g(4,1,2,5,3,0)           ! 58
      amp3j(59) = FD1ePV3g(5,3,2,4,1,.false.)    ! 59
      amp3j(60) = FDePV4g(5,3,2,4,1,1)           ! 60
      amp3j(61) =-FDePV4g(5,3,2,4,1,2)           ! 61
      amp3j(62) =-FD1ePV3g(1,4,2,5,3,.false.)    ! 62
      amp3j(63) = FDePV4g(1,4,2,5,3,1)           ! 63
      amp3j(64) = FDePV4g(1,4,2,5,3,2)           ! 64
      amp3j(65) =-FD1ePV3g(3,5,2,4,1,.false.)    ! 65
      amp3j(66) = FDePV4g(3,5,2,4,1,1)           ! 66
      amp3j(67) = FDePV4g(3,5,2,4,1,2)           ! 67
      amp3j(68) =-FD2gBox(4,1,2,5,3,.true.)      ! 68
      amp3j(69) =-FD2gBox(5,3,2,4,1,.false.)     ! 69
      amp3j(70) = FD4gBox(5,3,2,4,1)             ! 70
      amp3j(71) =-FD2gBox(4,1,2,5,3,.false.)     ! 71
      amp3j(72) = FD4gBox(4,1,2,5,3)             ! 72
      amp3j(73) =-FD3gTri(5,3,2,4,1)             ! 73
      amp3j(74) =-FDePV4g(5,3,2,4,1,3)           ! 74
      amp3j(75) = FDeffPen(5,3,2,4,1)            ! 75
      amp3j(76) =-FD3gTri(4,1,2,5,3)             ! 76
      amp3j(77) = FDeffPen(4,1,2,5,3)            ! 77
      amp3j(78) =-FDePV4g(4,1,2,5,3,3)           ! 78
      amp3j(79) =-FD1ePV3g(4,5,1,2,3,.true.)     ! 79 
      amp3j(80) =-FDePV4g(5,4,1,2,3,0)           ! 80
      amp3j(81) =-FDePV4g(2,3,1,5,4,0)           ! 81
      amp3j(82) = FD1ePV3g(3,2,1,4,5,.false.)    ! 82
      amp3j(83) =-FDePV4g(3,2,1,4,5,1)           ! 83
      amp3j(84) =-FDePV4g(3,2,1,4,5,2)           ! 84
      amp3j(85) = FD1ePV3g(2,3,1,4,5,.false.)    ! 85
      amp3j(86) =-FDePV4g(2,3,1,4,5,1)           ! 86
      amp3j(87) =-FDePV4g(2,3,1,4,5,2)           ! 87
      amp3j(88) =-FD1ePV3g(5,4,1,2,3,.false.)    ! 88
      amp3j(89) =-FDePV4g(5,4,1,2,3,1)           ! 89
      amp3j(90) = FDePV4g(5,4,1,2,3,2)           ! 90
      amp3j(91) = FD1ePV3g(4,5,1,2,3,.false.)    ! 91
      amp3j(92) =-FDePV4g(4,5,1,2,3,1)           ! 92
      amp3j(93) =-FDePV4g(4,5,1,2,3,2)           ! 93
      amp3j(94) =-FD2gBox(4,5,1,3,2,.false.)     ! 94
      amp3j(95) = FD4gBox(4,5,1,3,2)             ! 95
      amp3j(96) =-FD2gBox(4,5,1,2,3,.true.)      ! 96 
      amp3j(97) =-FD2gBox(2,3,1,4,5,.false.)     ! 97
      amp3j(98) =-FD4gBox(2,3,1,5,4)             ! 98
      amp3j(99) =-FD3gTri(4,5,1,2,3)             ! 99
      amp3j(100)=-FDePV4g(4,5,1,2,3,3)           ! 100
      amp3j(101)= FDeffPen(4,5,1,2,3)            ! 101
      amp3j(102)=-FD3gTri(2,3,1,4,5)             ! 102
      amp3j(103)= FDeffPen(2,3,1,4,5)            ! 103
      amp3j(104)=-FDePV4g(2,3,1,4,5,3)           ! 104
      amp3j(105)= FD1ePV3g(4,1,5,2,3,.false.)    ! 105
      amp3j(106)= FDePV4g(4,1,5,2,3,1)           ! 106
      amp3j(107)= FDePV4g(4,1,5,2,3,2)           ! 107
      amp3j(108)= FD1ePV3g(3,2,5,4,1,.false.)    ! 108
      amp3j(109)= FDePV4g(3,2,5,1,4,1)           ! 109
      amp3j(110)=-FDePV4g(3,2,5,1,4,2)           ! 110
      amp3j(111)= FD1ePV3g(2,3,5,4,1,.false.)    ! 111
      amp3j(112)=-FDePV4g(2,3,5,4,1,1)           ! 112
      amp3j(113)= FDePV4g(2,3,5,4,1,2)           ! 113
      amp3j(114)=-FD1ePV3g(1,4,5,2,3,.false.)    ! 114 
      amp3j(115)= FDePV4g(1,4,5,2,3,1)           ! 115
      amp3j(116)=-FDePV4g(1,4,5,2,3,2)           ! 116
      amp3j(117)= FD1ePV3g(4,1,5,2,3,.true.)     ! 117 
      amp3j(118)=-FDePV4g(4,1,5,2,3,0)           ! 118
      amp3j(119)=-FDePV4g(2,3,5,4,1,0)           ! 119
      amp3j(120)=-FD2gBox(4,1,5,3,2,.false.)     ! 120
      amp3j(121)=-FD4gBox(4,1,5,3,2)             ! 121
      amp3j(122)= FD2gBox(2,3,5,1,4,.false.)     ! 122
      amp3j(123)= FD4gBox(2,3,5,1,4)             ! 123
      amp3j(124)= FD2gBox(4,1,5,2,3,.true.)      ! 124
      amp3j(125)=-FD3gTri(4,1,5,2,3)             ! 125
      amp3j(126)= FDePV4g(4,1,5,2,3,3)           ! 126
      amp3j(127)= FDeffPen(4,1,5,2,3)            ! 127
      amp3j(128)=-FD3gTri(2,3,5,4,1)             ! 128
      amp3j(129)= FDeffPen(2,3,5,4,1)            ! 129
      amp3j(130)= FDePV4g(2,3,5,4,1,3)           ! 130
      amp3j(131)= FD1ePV3g(5,1,3,2,4,.false.)    ! 131
      amp3j(132)= FDePV4g(5,1,3,2,4,1)           ! 132
      amp3j(133)= FDePV4g(5,1,3,2,4,2)           ! 133
      amp3j(134)= FD1ePV3g(4,2,3,5,1,.false.)    ! 134
      amp3j(135)= FDePV4g(4,2,3,5,1,1)           ! 135
      amp3j(136)=-FDePV4g(4,2,3,5,1,2)           ! 136
      amp3j(137)= FD1ePV3g(5,1,3,2,4,.true.)     ! 137
      amp3j(138)=-FDePV4g(5,1,3,2,4,0)           ! 138
      amp3j(139)=-FDePV4g(2,4,3,5,1,0)           ! 139
      amp3j(140)=-FD1ePV3g(2,4,3,5,1,.false.)    ! 140
      amp3j(141)= FDePV4g(2,4,3,5,1,1)           ! 141
      amp3j(142)= FDePV4g(2,4,3,5,1,2)           ! 142
      amp3j(143)=-FD1ePV3g(1,5,3,2,4,.false.)    ! 143
      amp3j(144)= FDePV4g(1,5,3,2,4,1)           ! 144
      amp3j(145)=-FDePV4g(1,5,3,2,4,2)           ! 145
      amp3j(146)=-FD2gBox(5,1,3,2,4,.true.)      ! 146 
      amp3j(147)=-FD2gBox(5,1,3,4,2,.false.)     ! 147
      amp3j(148)= FD4gBox(5,1,3,4,2)             ! 148
      amp3j(149)=-FD2gBox(2,4,3,5,1,.false.)     ! 149
      amp3j(150)=-FD4gBox(2,4,3,5,1)             ! 150
      amp3j(151)= FD3gTri(5,1,3,2,4)             ! 151 
      amp3j(152)= FDePV4g(5,1,3,2,4,3)           ! 152
      amp3j(153)=-FDeffPen(5,1,3,2,4)            ! 153
      amp3j(154)=-FD3gTri(2,4,3,5,1)             ! 154
      amp3j(155)= FDeffPen(2,4,3,5,1)            ! 155
      amp3j(156)= FDePV4g(2,4,3,5,1,3)           ! 156
      amp3j(157)= FD1ePV3g(3,1,4,2,5,.false.)    ! 157
      amp3j(158)= FDePV4g(3,1,4,2,5,1)           ! 158
      amp3j(159)= FDePV4g(3,1,4,2,5,2)           ! 159
      amp3j(160)= FD1ePV3g(5,2,4,3,1,.false.)    ! 160
      amp3j(161)= FDePV4g(5,2,4,3,1,1)           ! 161
      amp3j(162)= FDePV4g(5,2,4,3,1,2)           ! 162
      amp3j(163)= FD1ePV3g(1,3,4,2,5,.false.)    ! 163
      amp3j(164)= FDePV4g(1,3,4,2,5,1)           ! 164
      amp3j(165)= FDePV4g(1,3,4,2,5,2)           ! 165
      amp3j(166)=-FD1ePV3g(3,1,4,2,5,.true.)     ! 166
      amp3j(167)= FDePV4g(2,5,4,3,1,0)           ! 167
      amp3j(168)= FDePV4g(3,1,4,2,5,0)           ! 168
      amp3j(169)=-FD1ePV3g(2,5,4,3,1,.false.)    ! 169
      amp3j(170)= FDePV4g(2,5,4,3,1,1)           ! 170
      amp3j(171)=-FDePV4g(2,5,4,3,1,2)           ! 171
      amp3j(172)=-FD2gBox(2,5,4,3,1,.false.)     ! 172
      amp3j(173)=-FD4gBox(2,5,4,3,1)             ! 173
      amp3j(174)=-FD2gBox(3,1,4,2,5,.true.)      ! 174
      amp3j(175)=-FD2gBox(3,1,4,5,2,.false.)     ! 175
      amp3j(176)=-FD4gBox(3,1,4,5,2)             ! 176
      amp3j(177)= FD3gTri(2,5,4,3,1)             ! 177
      amp3j(178)=-FDePV4g(2,5,4,3,1,3)           ! 178
      amp3j(179)= FDeffPen(5,2,4,3,1)            ! 179
      amp3j(180)= FD3gTri(3,1,4,2,5)             ! 180
      amp3j(181)=-FDeffPen(3,1,4,2,5)            ! 181
      amp3j(182)=-FDePV4g(3,1,4,2,5,3)           ! 182
      amp3j(183)= FD1ePV3g(5,3,1,2,4,.true.)     ! 183
      amp3j(184)=-FDePV4g(5,3,1,2,4,0)           ! 184
      amp3j(185)=-FDePV4g(2,4,1,5,3,0)           ! 185
      amp3j(186)=-FD1ePV3g(4,2,1,5,3,.false.)    ! 186
      amp3j(187)=-FDePV4g(4,2,1,5,3,1)           ! 187
      amp3j(188)=-FDePV4g(4,2,1,5,3,2)           ! 188
      amp3j(189)= FD1ePV3g(5,3,1,2,4,.false.)    ! 189
      amp3j(190)= FDePV4g(5,3,1,2,4,1)           ! 190
      amp3j(191)= FDePV4g(5,3,1,2,4,2)           ! 191
      amp3j(192)= FD1ePV3g(2,4,1,5,3,.false.)    ! 192
      amp3j(193)=-FDePV4g(2,4,1,5,3,1)           ! 193
      amp3j(194)= FDePV4g(2,4,1,5,3,2)           ! 194
      amp3j(195)=-FD1ePV3g(3,5,1,2,4,.false.)    ! 195
      amp3j(196)= FDePV4g(3,5,1,2,4,1)           ! 196
      amp3j(197)=-FDePV4g(3,5,1,2,4,2)           ! 197
      amp3j(198)=-FD2gBox(3,5,1,2,4,.true.)      ! 198
      amp3j(199)= FD2gBox(5,3,1,4,2,.false.)     ! 199
      amp3j(200)= FD4gBox(5,3,1,4,2)             ! 200
      amp3j(201)=-FD2gBox(4,2,1,3,5,.false.)     ! 201
      amp3j(202)=-FD4gBox(2,4,1,5,3)             ! 202
      amp3j(203)=-FD3gTri(5,3,1,2,4)             ! 203
      amp3j(204)= FDePV4g(5,3,1,2,4,3)           ! 204
      amp3j(205)= FDeffPen(5,3,1,2,4)            ! 205
      amp3j(206)=-FD3gTri(2,4,1,5,3)             ! 206
      amp3j(207)= FDeffPen(2,4,1,5,3)            ! 207
      amp3j(208)= FDePV4g(2,4,1,5,3,3)           ! 208
      amp3j(209)= FD1ePV3g(3,1,2,5,4,.false.)    ! 209
      amp3j(210)= FDePV4g(3,1,2,5,4,1)           ! 210
      amp3j(211)=-FDePV4g(3,1,2,5,4,2)           ! 211
      amp3j(212)= FD1ePV3g(5,4,2,3,1,.true.)     ! 212
      amp3j(213)= FDePV4g(4,5,2,3,1,0)           ! 213
      amp3j(214)=-FDePV4g(3,1,2,5,4,0)           ! 214
      amp3j(215)= FD1ePV3g(1,3,2,5,4,.false.)    ! 215
      amp3j(216)= FDePV4g(1,3,2,5,4,1)           ! 216
      amp3j(217)=-FDePV4g(1,3,2,5,4,2)           ! 217
      amp3j(218)=-FD1ePV3g(5,4,2,3,1,.false.)    ! 218
      amp3j(219)=-FDePV4g(5,4,2,3,1,1)           ! 219
      amp3j(220)= FDePV4g(5,4,2,3,1,2)           ! 220
      amp3j(221)= FD1ePV3g(4,5,2,3,1,.false.)    ! 221
      amp3j(222)=-FDePV4g(4,5,2,3,1,1)           ! 222
      amp3j(223)=-FDePV4g(4,5,2,3,1,2)           ! 223
      amp3j(224)=-FD2gBox(5,4,2,3,1,.false.)     ! 224
      amp3j(225)= FD4gBox(5,4,2,3,1)             ! 225
      amp3j(226)=-FD2gBox(3,1,2,5,4,.true.)      ! 226
      amp3j(227)= FD2gBox(3,1,2,5,4,.false.)     ! 227
      amp3j(228)=-FD4gBox(3,1,2,5,4)             ! 228
      amp3j(229)=-FD3gTri(4,5,2,3,1)             ! 229
      amp3j(230)= FDePV4g(5,4,2,3,1,3)           ! 230
      amp3j(231)=-FDeffPen(5,4,2,3,1)            ! 231
      amp3j(232)= FD3gTri(3,1,2,5,4)             ! 232
      amp3j(233)=-FDeffPen(3,1,2,5,4)            ! 233
      amp3j(234)= FDePV4g(3,1,2,5,4,3)           ! 234
      amp3j(235)= FD1ePV3g(3,1,5,2,4,.false.)    ! 235
      amp3j(236)= FDePV4g(3,1,5,2,4,1)           ! 236
      amp3j(237)= FDePV4g(3,1,5,2,4,2)           ! 237
      amp3j(238)= FD1ePV3g(4,2,5,3,1,.false.)    ! 238
      amp3j(239)= FDePV4g(4,2,5,3,1,1)           ! 239
      amp3j(240)= FDePV4g(4,2,5,3,1,2)           ! 240
      amp3j(241)= FD1ePV3g(1,3,5,2,4,.false.)    ! 241
      amp3j(242)= FDePV4g(1,3,5,2,4,1)           ! 242
      amp3j(243)= FDePV4g(1,3,5,2,4,2)           ! 243
      amp3j(244)=-FD1ePV3g(2,4,5,3,1,.false.)    ! 244
      amp3j(245)= FDePV4g(2,4,5,3,1,1)           ! 245
      amp3j(246)=-FDePV4g(2,4,5,3,1,2)           ! 246
      amp3j(247)= FD1ePV3g(3,1,5,2,4,.true.)     ! 247
      amp3j(248)= FDePV4g(2,4,5,3,1,0)           ! 248
      amp3j(249)= FDePV4g(3,1,5,2,4,0)           ! 249
      amp3j(250)= FD2gBox(4,2,5,3,1,.false.)     ! 250
      amp3j(251)=-FD4gBox(2,4,5,3,1)             ! 251
      amp3j(252)=-FD2gBox(3,1,5,4,2,.false.)     ! 252
      amp3j(253)=-FD4gBox(3,1,5,4,2)             ! 253
      amp3j(254)= FD2gBox(3,1,5,2,4,.true.)      ! 254
      amp3j(255)= FD3gTri(2,4,5,3,1)             ! 255
      amp3j(256)=-FDePV4g(2,4,5,3,1,3)           ! 256
      amp3j(257)= FDeffPen(4,2,5,3,1)            ! 257
      amp3j(258)= FD3gTri(3,1,5,2,4)             ! 258
      amp3j(259)=-FDeffPen(3,1,5,2,4)            ! 259
      amp3j(260)=-FDePV4g(3,1,5,2,4,3)           ! 260
      amp3j(261)= FD1ePV3g(5,1,2,4,3,.false.)    ! 261
      amp3j(262)= FDePV4g(5,1,2,4,3,1)           ! 262
      amp3j(263)=-FDePV4g(5,1,2,4,3,2)           ! 263
      amp3j(264)= FD1ePV3g(4,3,2,5,1,.true.)     ! 264
      amp3j(265)= FDePV4g(5,1,2,4,3,0)           ! 265
      amp3j(266)= FDePV4g(4,3,2,5,1,0)           ! 266
      amp3j(267)= FD1ePV3g(4,3,2,5,1,.false.)    ! 267
      amp3j(268)= FDePV4g(4,3,2,5,1,1)           ! 268
      amp3j(269)=-FDePV4g(4,3,2,5,1,2)           ! 269
      amp3j(270)=-FD1ePV3g(3,4,2,5,1,.false.)    ! 270
      amp3j(271)= FDePV4g(3,4,2,5,1,1)           ! 271
      amp3j(272)= FDePV4g(3,4,2,5,1,2)           ! 272
      amp3j(273)=-FD1ePV3g(1,5,2,4,3,.false.)    ! 273
      amp3j(274)= FDePV4g(1,5,2,4,3,1)           ! 274
      amp3j(275)= FDePV4g(1,5,2,4,3,2)           ! 275
      amp3j(276)=-FD2gBox(5,1,2,4,3,.true.)      ! 276
      amp3j(277)= FD2gBox(5,1,2,3,4,.false.)     ! 277
      amp3j(278)=-FD4gBox(5,1,2,3,4)             ! 278
      amp3j(279)= FD2gBox(3,4,2,5,1,.false.)     ! 279
      amp3j(280)= FD4gBox(4,3,2,5,1)             ! 280
      amp3j(281)=-FD3gTri(5,1,2,4,3)             ! 281
      amp3j(282)=-FDePV4g(5,1,2,4,3,3)           ! 282
      amp3j(283)= FDeffPen(5,1,2,4,3)            ! 283
      amp3j(284)=-FD3gTri(4,3,2,5,1)             ! 284
      amp3j(285)= FDeffPen(4,3,2,5,1)            ! 285
      amp3j(286)=-FDePV4g(4,3,2,5,1,3)           ! 286
      amp3j(287)= FD1ePV3g(4,3,1,2,5,.true.)     ! 287
      amp3j(288)=-FDePV4g(2,5,1,4,3,0)           ! 288
      amp3j(289)=-FDePV4g(4,3,1,2,5,0)           ! 289
      amp3j(290)=-FD1ePV3g(5,2,1,4,3,.false.)    ! 290
      amp3j(291)=-FDePV4g(5,2,1,4,3,1)           ! 291
      amp3j(292)=-FDePV4g(5,2,1,4,3,2)           ! 292
      amp3j(293)= FD1ePV3g(4,3,1,2,5,.false.)    ! 293
      amp3j(294)= FDePV4g(4,3,1,2,5,1)           ! 294
      amp3j(295)= FDePV4g(4,3,1,2,5,2)           ! 295
      amp3j(296)=-FD1ePV3g(3,4,1,2,5,.false.)    ! 296
      amp3j(297)= FDePV4g(3,4,1,2,5,1)           ! 297
      amp3j(298)=-FDePV4g(3,4,1,2,5,2)           ! 298
      amp3j(299)= FD1ePV3g(2,5,1,4,3,.false.)    ! 299
      amp3j(300)=-FDePV4g(2,5,1,4,3,1)           ! 300
      amp3j(301)= FDePV4g(2,5,1,4,3,2)           ! 301
      amp3j(302)= FD2gBox(4,3,1,2,5,.true.)      ! 302
      amp3j(303)=-FD2gBox(2,5,1,4,3,.false.)     ! 303
      amp3j(304)=-FD4gBox(2,5,1,4,3)             ! 304
      amp3j(305)= FD2gBox(4,3,1,5,2,.false.)     ! 305
      amp3j(306)= FD4gBox(4,3,1,5,2)             ! 306
      amp3j(307)=-FD3gTri(2,5,1,4,3)             ! 307
      amp3j(308)= FDePV4g(2,5,1,4,3,3)           ! 308
      amp3j(309)= FDeffPen(2,5,1,4,3)            ! 309
      amp3j(310)=-FD3gTri(4,3,1,2,5)             ! 310
      amp3j(311)= FDeffPen(4,3,1,2,5)            ! 311
      amp3j(312)= FDePV4g(4,3,1,2,5,3)           ! 312
      amp3j(313)= FD1ePV3g(2,1,4,5,3,.false.)    ! 313
      amp3j(314)=-FDePV4g(2,1,4,5,3,1)           ! 314
      amp3j(315)= FDePV4g(2,1,4,5,3,2)           ! 315
      amp3j(316)=-FD1ePV3g(1,2,4,5,3,.false.)    ! 316
      amp3j(317)=-FDePV4g(1,2,4,5,3,1)           ! 317
      amp3j(318)=-FDePV4g(1,2,4,5,3,2)           ! 318
      amp3j(319)=-FD1ePV3g(5,3,4,1,2,.false.)    ! 319
      amp3j(320)=-FDePV4g(5,3,4,1,2,1)           ! 320
      amp3j(321)=-FDePV4g(5,3,4,1,2,2)           ! 321
      amp3j(322)=-FD1ePV3g(1,2,4,5,3,.true.)     ! 322
      amp3j(323)= FDePV4g(5,3,4,1,2,0)           ! 323
      amp3j(324)= FDePV4g(1,2,4,5,3,0)           ! 324
      amp3j(325)= FD1ePV3g(3,5,4,1,2,.false.)    ! 325
      amp3j(326)=-FDePV4g(3,5,4,1,2,1)           ! 326
      amp3j(327)= FDePV4g(3,5,4,1,2,2)           ! 327
      amp3j(328)= FD2gBox(5,3,4,1,2,.false.)     ! 328
      amp3j(329)=-FD4gBox(3,5,4,1,2)             ! 329
      amp3j(330)= FD2gBox(1,2,4,3,5,.true.)      ! 330
      amp3j(331)=-FD2gBox(1,2,4,3,5,.false.)     ! 331
      amp3j(332)=-FD4gBox(1,2,4,3,5)             ! 332
      amp3j(333)= FD3gTri(5,3,4,1,2)             ! 333
      amp3j(334)=-FDePV4g(5,3,4,1,2,3)           ! 334
      amp3j(335)=-FDeffPen(5,3,4,1,2)            ! 335
      amp3j(336)=-FD3gTri(1,2,4,3,5)             ! 336
      amp3j(337)= FDeffPen(1,2,4,3,5)            ! 337
      amp3j(338)= FDePV4g(1,2,4,3,5,3)           ! 338
      amp3j(339)=-FD1ePV3g(2,1,3,5,4,.false.)    ! 339
      amp3j(340)= FDePV4g(2,1,3,5,4,1)           ! 340
      amp3j(341)= FDePV4g(2,1,3,5,4,2)           ! 341
      amp3j(342)= FD1ePV3g(1,2,3,5,4,.false.)    ! 342
      amp3j(343)=-FDePV4g(1,2,3,4,5,1)           ! 343
      amp3j(344)= FDePV4g(1,2,3,4,5,2)           ! 344
      amp3j(345)= FD1ePV3g(1,2,3,5,4,.true.)     ! 345
      amp3j(346)= FDePV4g(5,4,3,1,2,0)           ! 346
      amp3j(347)= FDePV4g(1,2,3,5,4,0)           ! 347
      amp3j(348)= FD1ePV3g(5,4,3,1,2,.false.)    ! 348
      amp3j(349)= FDePV4g(5,4,3,1,2,1)           ! 349
      amp3j(350)=-FDePV4g(5,4,3,1,2,2)           ! 350
      amp3j(351)=-FD1ePV3g(4,5,3,1,2,.false.)    ! 351
      amp3j(352)= FDePV4g(4,5,3,1,2,1)           ! 352
      amp3j(353)= FDePV4g(4,5,3,1,2,2)           ! 353
      amp3j(354)=-FD2gBox(5,4,3,1,2,.false.)     ! 354
      amp3j(355)= FD4gBox(5,4,3,1,2)             ! 355
      amp3j(356)= FD2gBox(1,2,3,5,4,.true.)      ! 356
      amp3j(357)=-FD2gBox(1,2,3,5,4,.false.)     ! 357
      amp3j(358)= FD4gBox(1,2,3,5,4)             ! 358
      amp3j(359)=-FD3gTri(5,4,3,1,2)             ! 359
      amp3j(360)=-FDePV4g(5,4,3,1,2,3)           ! 360
      amp3j(361)= FDeffPen(5,4,3,1,2)            ! 361
      amp3j(362)=-FD3gTri(1,2,3,5,4)             ! 362
      amp3j(363)= FDeffPen(1,2,3,5,4)            ! 363
      amp3j(364)=-FDePV4g(1,2,3,5,4,3)           ! 364
      amp3j(365)= FD1ePV3g(2,1,5,4,3,.false.)    ! 365
      amp3j(366)=-FDePV4g(2,1,5,4,3,1)           ! 366
      amp3j(367)= FDePV4g(2,1,5,4,3,2)           ! 367
      amp3j(368)=-FD1ePV3g(1,2,5,4,3,.false.)    ! 368
      amp3j(369)=-FDePV4g(1,2,5,4,3,1)           ! 369
      amp3j(370)=-FDePV4g(1,2,5,4,3,2)           ! 370
      amp3j(371)=-FD1ePV3g(4,3,5,1,2,.false.)    ! 371
      amp3j(372)=-FDePV4g(4,3,5,1,2,1)           ! 372
      amp3j(373)=-FDePV4g(4,3,5,1,2,2)           ! 373
      amp3j(374)= FD1ePV3g(3,4,5,1,2,.false.)    ! 374
      amp3j(375)=-FDePV4g(3,4,5,1,2,1)           ! 375
      amp3j(376)= FDePV4g(3,4,5,1,2,2)           ! 376
      amp3j(377)= FD1ePV3g(1,2,5,4,3,.true.)     ! 377
      amp3j(378)= FDePV4g(4,3,5,1,2,0)           ! 378
      amp3j(379)= FDePV4g(1,2,5,4,3,0)           ! 379
      amp3j(380)= FD2gBox(4,3,5,1,2,.false.)     ! 380
      amp3j(381)= FD4gBox(4,3,5,1,2)             ! 381
      amp3j(382)=-FD2gBox(1,2,5,3,4,.false.)     ! 382
      amp3j(383)=-FD4gBox(1,2,5,3,4)             ! 383
      amp3j(384)=-FD2gBox(1,2,5,3,4,.true.)      ! 384
      amp3j(385)= FD3gTri(4,3,5,1,2)             ! 385
      amp3j(386)=-FDePV4g(4,3,5,1,2,3)           ! 386
      amp3j(387)=-FDeffPen(4,3,5,1,2)            ! 387
      amp3j(388)=-FD3gTri(1,2,5,3,4)             ! 388
      amp3j(389)= FDeffPen(1,2,5,3,4)            ! 389
      amp3j(390)= FDePV4g(1,2,5,3,4,3)           ! 390

c      print*, "boxdbleSin",boxdbleSin



      do i=1,24
         color(i) =DCMPLX(0d0,0d0)
         colorPentagons(i)=DCMPLX(0d0,0d0)
         colorHexagons(i)=DCMPLX(0d0,0d0)
      enddo

*------------- Further contributions from full theory ------------
      if(HloopTop.or.HloopBottom.or.AloopTop.or.AloopBottom) then
         If(peneva) then
         pendbleSin=.false.
         penquadSin=.false.
c     print*, "peneva", peneva
c$$$*------------ [ 1,2,3,(4,5=15) ] ------------
      do mu=0,4
         Kpen(mu,1)=K(mu,1)
         Kpen(mu,2)=K(mu,2)
         Kpen(mu,3)=K(mu,3)
         Kpen(mu,4)=Kij(mu,4,5)
         if(mu.eq.4) exit
         epsgPen(mu,1)=epsgC(mu,1)
         epsgPen(mu,2)=epsgC(mu,2)
         epsgPen(mu,3)=epsgC(mu,3)
         epsgPen(mu,4)=Jij(mu,4,5)
      enddo
      PenSingular(1)=.false.
      PenSingular(2)=.false.
c      CALL auxiliary_calculation(epsgPen(0,1),Kpen(0,1))
c      CALL calculate_E_functions() !depend on q_dot_q(i,j) only
      pent(1,2,3,15)=-(FDeffPenLoop(4,3,2,1)+FDeffPenLoop(2,3,4,1)
     &     +FDeffPenLoop(2,1,4,3)+FDeffPenLoop(3,2,1,4))
      pent(1,3,15,2)=-(FDeffPenLoop(2,1,3,4)+FDeffPenLoop(3,4,2,1)
     &     +FDeffPenLoop(4,2,1,3)+FDeffPenLoop(2,4,3,1))
      pent(1,15,2,3)=-(FDeffPenLoop(2,3,1,4)+FDeffPenLoop(3,1,4,2)
     &     +FDeffPenLoop(4,2,3,1)+FDeffPenLoop(3,2,4,1))
      IF(PenSingular(1)) then
         PendbleSin=.true.
      endif
      IF(PenSingular(2)) then
         PenquadSin=.true.
      endif
cc start check
c      print*,'pent(1,2,3,15)',pent(1,2,3,15)*dsqrt(fac)!*(-1)
c      print*,'pent(1,3,15,2)',pent(1,3,15,2)*dsqrt(fac)*(-1)
c      print*,'pent(1,15,2,3)',pent(1,15,2,3)*dsqrt(fac)*(-1)
c      STOP
cc end check
*------------ [ 3,4,5,(1,2=6) ] ------------
      do mu=0,4
         Kpen(mu,1)=K(mu,3)
         Kpen(mu,2)=K(mu,4)
         Kpen(mu,3)=K(mu,5)
         Kpen(mu,4)=Kij(mu,1,2)
         if(mu.eq.4) exit
         epsgPen(mu,1)=epsgC(mu,3)
         epsgPen(mu,2)=epsgC(mu,4)
         epsgPen(mu,3)=epsgC(mu,5)
         epsgPen(mu,4)=Jij(mu,1,2)
      enddo
      PenSingular(1)=.false.
      PenSingular(2)=.false.
c      CALL auxiliary_calculation(epsgPen(0,1),Kpen(0,1))
c      CALL calculate_E_functions() !depend on q_dot_q(i,j) only
      pent(3,4,5,6)=-(FDeffPenLoop(4,3,2,1)+FDeffPenLoop(2,3,4,1)
     &     +FDeffPenLoop(2,1,4,3)+FDeffPenLoop(3,2,1,4))
      pent(3,5,6,4)=-(FDeffPenLoop(2,1,3,4)+FDeffPenLoop(3,4,2,1)
     &     +FDeffPenLoop(4,2,1,3)+FDeffPenLoop(2,4,3,1))
      pent(3,6,4,5)=-(FDeffPenLoop(2,3,1,4)+FDeffPenLoop(3,1,4,2)
     &     +FDeffPenLoop(4,2,3,1)+FDeffPenLoop(3,2,4,1))
      IF(PenSingular(1)) then
         PendbleSin=.true.
      endif
      IF(PenSingular(2)) then
         PenquadSin=.true.
      endif
cc start check
c      print*,'pent(3,4,5,6)',pent(3,4,5,6)*dsqrt(fac)*(-1)
c      print*,'pent(3,5,6,4)',pent(3,5,6,4)*dsqrt(fac)*(-1)
c      print*,'pent(3,6,4,5)',pent(3,6,4,5)*dsqrt(fac)*(-1)
c      print*,''
      if(lgauge2) then
         print*,'sum^2 [ 3,4,5,(1,2=6) ]',
     &        dreal(pent(3,4,5,6)+pent(3,5,6,4)+pent(3,6,4,5))
     &       *dreal(pent(3,4,5,6)+pent(3,5,6,4)+pent(3,6,4,5))
     &       +dimag(pent(3,4,5,6)+pent(3,5,6,4)+pent(3,6,4,5))
     &       *dimag(pent(3,4,5,6)+pent(3,5,6,4)+pent(3,6,4,5))
         print*,''
      endif
c      STOP
cc end check
*------------ [ 2,4,5,(1,3=7) ] ------------
      do mu=0,4
         Kpen(mu,1)=K(mu,2)
         Kpen(mu,2)=K(mu,4)
         Kpen(mu,3)=K(mu,5)
         Kpen(mu,4)=Kij(mu,1,3)
         if(mu.eq.4) exit
         epsgPen(mu,1)=epsgC(mu,2)
         epsgPen(mu,2)=epsgC(mu,4)
         epsgPen(mu,3)=epsgC(mu,5)
         epsgPen(mu,4)=Jij(mu,1,3)
      enddo
c      CALL auxiliary_calculation(epsgPen(0,1),Kpen(0,1))
c      CALL calculate_E_functions() !depend on q_dot_q(i,j) only
      PenSingular(1)=.false.
      PenSingular(2)=.false.
      pent(2,4,5,7)=-(FDeffPenLoop(4,3,2,1)+FDeffPenLoop(2,3,4,1)
     &     +FDeffPenLoop(2,1,4,3)+FDeffPenLoop(3,2,1,4))
      pent(2,5,7,4)=-(FDeffPenLoop(2,1,3,4)+FDeffPenLoop(3,4,2,1)
     &     +FDeffPenLoop(4,2,1,3)+FDeffPenLoop(2,4,3,1))
      pent(2,7,4,5)=-(FDeffPenLoop(2,3,1,4)+FDeffPenLoop(3,1,4,2)
     &     +FDeffPenLoop(4,2,3,1)+FDeffPenLoop(3,2,4,1))
      IF(PenSingular(1)) then
         PendbleSin=.true.
      endif
      IF(PenSingular(2)) then
         PenquadSin=.true.
      endif
cc start check
c      print*,'pent(2,4,5,7)',pent(2,4,5,7)*dsqrt(fac)!*(-1)
c      print*,'pent(2,5,7,4)',pent(2,5,7,4)*dsqrt(fac)*(-1)
c      print*,'pent(2,7,4,5)',pent(2,7,4,5)*dsqrt(fac)!*(-1)
c      print*,''
      if(lgauge2) then
         print*,'sum^2 [ 2,4,5,(1,3=7) ]',
     &        dreal(pent(2,4,5,7)+pent(2,5,7,4)+pent(2,7,4,5))
     &       *dreal(pent(2,4,5,7)+pent(2,5,7,4)+pent(2,7,4,5))
     &       +dimag(pent(2,4,5,7)+pent(2,5,7,4)+pent(2,7,4,5))
     &       *dimag(pent(2,4,5,7)+pent(2,5,7,4)+pent(2,7,4,5))
         print*,''
      endif
cc end check
*------------ [ 2,3,5,(1,4=8) ] ------------
      do mu=0,4
         Kpen(mu,1)=K(mu,2)
         Kpen(mu,2)=K(mu,3)
         Kpen(mu,3)=K(mu,5)
         Kpen(mu,4)=Kij(mu,1,4)
         if(mu.eq.4) exit
         epsgPen(mu,1)=epsgC(mu,2)
         epsgPen(mu,2)=epsgC(mu,3)
         epsgPen(mu,3)=epsgC(mu,5)
         epsgPen(mu,4)=Jij(mu,1,4)
      enddo
c      CALL auxiliary_calculation(epsgPen(0,1),Kpen(0,1))
c      CALL calculate_E_functions() !depend on q_dot_q(i,j) only
      PenSingular(1)=.false.
      PenSingular(2)=.false.
      pent(2,3,5,8)=-(FDeffPenLoop(4,3,2,1)+FDeffPenLoop(2,3,4,1)
     &     +FDeffPenLoop(2,1,4,3)+FDeffPenLoop(3,2,1,4))
      pent(2,5,8,3)=-(FDeffPenLoop(2,1,3,4)+FDeffPenLoop(3,4,2,1)
     &     +FDeffPenLoop(4,2,1,3)+FDeffPenLoop(2,4,3,1))
      pent(2,8,3,5)=-(FDeffPenLoop(2,3,1,4)+FDeffPenLoop(3,1,4,2)
     &     +FDeffPenLoop(4,2,3,1)+FDeffPenLoop(3,2,4,1))
      IF(PenSingular(1)) then
         PendbleSin=.true.
      endif
      IF(PenSingular(2)) then
         PenquadSin=.true.
      endif
cc start check
c      print*,'pent(2,3,5,8)',pent(2,3,5,8)*dsqrt(fac)!*(-1)
c      print*,'pent(2,5,8,3)',pent(2,5,8,3)*dsqrt(fac)*(-1)
c      print*,'pent(2,8,3,5)',pent(2,8,3,5)*dsqrt(fac)*(-1)
c      print*,''
      if(lgauge2) then
         print*,'sum^2 [ 2,3,5,(1,4=8) ]',
     &        dreal(pent(2,3,5,8)+pent(2,5,8,3)+pent(2,8,3,5))
     &       *dreal(pent(2,3,5,8)+pent(2,5,8,3)+pent(2,8,3,5))
     &       +dimag(pent(2,3,5,8)+pent(2,5,8,3)+pent(2,8,3,5))
     &       *dimag(pent(2,3,5,8)+pent(2,5,8,3)+pent(2,8,3,5))
         print*,''
      endif
cc end check
*------------ [ 2,3,4,(1,5=9) ] ------------
      do mu=0,4
         Kpen(mu,1)=K(mu,2)
         Kpen(mu,2)=K(mu,3)
         Kpen(mu,3)=K(mu,4)
         Kpen(mu,4)=Kij(mu,1,5)
         if(mu.eq.4) exit
         epsgPen(mu,1)=epsgC(mu,2)
         epsgPen(mu,2)=epsgC(mu,3)
         epsgPen(mu,3)=epsgC(mu,4)
         epsgPen(mu,4)=Jij(mu,1,5)
      enddo
c      CALL auxiliary_calculation(epsgPen(0,1),Kpen(0,1))
c      CALL calculate_E_functions() !depend on q_dot_q(i,j) only
      PenSingular(1)=.false.
      PenSingular(2)=.false.
      pent(2,3,4,9)=-(FDeffPenLoop(4,3,2,1)+FDeffPenLoop(2,3,4,1)
     &     +FDeffPenLoop(2,1,4,3)+FDeffPenLoop(3,2,1,4))
      pent(2,4,9,3)=-(FDeffPenLoop(2,1,3,4)+FDeffPenLoop(3,4,2,1)
     &     +FDeffPenLoop(4,2,1,3)+FDeffPenLoop(2,4,3,1))
      pent(2,9,3,4)=-(FDeffPenLoop(2,3,1,4)+FDeffPenLoop(3,1,4,2)
     &     +FDeffPenLoop(4,2,3,1)+FDeffPenLoop(3,2,4,1))
      IF(PenSingular(1)) then
         PendbleSin=.true.
      endif
      IF(PenSingular(2)) then
         PenquadSin=.true.
      endif
cc start check
c      print*,'pent(2,3,4,9)',pent(2,3,4,9)*dsqrt(fac)*(-1)
c      print*,'pent(2,4,9,3)',pent(2,4,9,3)*dsqrt(fac)*(-1)
c      print*,'pent(2,9,3,4)',pent(2,9,3,4)*dsqrt(fac)!*(-1)
c      print*,''

      if(lgauge2) then
         print*,'sum^2 [ 2,3,4,(1,5=9) ]',
     &        dreal(pent(2,3,4,9)+pent(2,4,9,3)+pent(2,9,3,4))
     &       *dreal(pent(2,3,4,9)+pent(2,4,9,3)+pent(2,9,3,4))
     &       +dimag(pent(2,3,4,9)+pent(2,4,9,3)+pent(2,9,3,4))
     &       *dimag(pent(2,3,4,9)+pent(2,4,9,3)+pent(2,9,3,4))
         print*,''
      endif
cc end check
*------------ [ 1,4,5,(2,3=10) ] ------------
      do mu=0,4
         Kpen(mu,1)=K(mu,1)
         Kpen(mu,2)=K(mu,4)
         Kpen(mu,3)=K(mu,5)
         Kpen(mu,4)=Kij(mu,2,3)
         if(mu.eq.4) exit
         epsgPen(mu,1)=epsgC(mu,1)
         epsgPen(mu,2)=epsgC(mu,4)
         epsgPen(mu,3)=epsgC(mu,5)
         epsgPen(mu,4)=Jij(mu,2,3)
      enddo
c      CALL auxiliary_calculation(epsgPen(0,1),Kpen(0,1))
c      CALL calculate_E_functions() !depend on q_dot_q(i,j) only
      PenSingular(1)=.false.
      PenSingular(2)=.false.
      pent(1,4,5,10)=-(FDeffPenLoop(4,3,2,1)+FDeffPenLoop(2,3,4,1)
     &     +FDeffPenLoop(2,1,4,3)+FDeffPenLoop(3,2,1,4))
      pent(1,5,10,4)=-(FDeffPenLoop(2,1,3,4)+FDeffPenLoop(3,4,2,1)
     &     +FDeffPenLoop(4,2,1,3)+FDeffPenLoop(2,4,3,1))
      pent(1,10,4,5)=-(FDeffPenLoop(2,3,1,4)+FDeffPenLoop(3,1,4,2)
     &     +FDeffPenLoop(4,2,3,1)+FDeffPenLoop(3,2,4,1))
      IF(PenSingular(1)) then
         PendbleSin=.true.
      endif
      IF(PenSingular(2)) then
         PenquadSin=.true.
      endif
cc start check
c      print*,'pent(1,4,5,10)',pent(1,4,5,10)*dsqrt(fac)*(-1)
c      print*,'pent(1,5,10,4)',pent(1,5,10,4)*dsqrt(fac)*(-1)
c      print*,'pent(1,10,4,5)',pent(1,10,4,5)*dsqrt(fac)*(-1)
c$$$      print*,''
      if(lgauge2) then
         print*,'sum^2 [ 1,4,5,(2,3=10) ]',
     &        dreal(pent(1,4,5,10)+pent(1,5,10,4)+pent(1,10,4,5))
     &       *dreal(pent(1,4,5,10)+pent(1,5,10,4)+pent(1,10,4,5))
     &       +dimag(pent(1,4,5,10)+pent(1,5,10,4)+pent(1,10,4,5))
     &       *dimag(pent(1,4,5,10)+pent(1,5,10,4)+pent(1,10,4,5))
         print*,''
      endif
c      STOP
cc end check
*------------ [ 1,3,5,(2,4=11) ] ------------
      do mu=0,4
         Kpen(mu,1)=K(mu,1)
         Kpen(mu,2)=K(mu,3)
         Kpen(mu,3)=K(mu,5)
         Kpen(mu,4)=Kij(mu,2,4)
         if(mu.eq.4) exit
         epsgPen(mu,1)=epsgC(mu,1)
         epsgPen(mu,2)=epsgC(mu,3)
         epsgPen(mu,3)=epsgC(mu,5)
         epsgPen(mu,4)=Jij(mu,2,4)
      enddo
c      CALL auxiliary_calculation(epsgPen(0,1),Kpen(0,1))
c      CALL calculate_E_functions() !depend on q_dot_q(i,j) only
      PenSingular(1)=.false.
      PenSingular(2)=.false.
      pent(1,3,5,11)=-(FDeffPenLoop(4,3,2,1)+FDeffPenLoop(2,3,4,1)
     &     +FDeffPenLoop(2,1,4,3)+FDeffPenLoop(3,2,1,4))
      pent(1,5,11,3)=-(FDeffPenLoop(2,1,3,4)+FDeffPenLoop(3,4,2,1)
     &     +FDeffPenLoop(4,2,1,3)+FDeffPenLoop(2,4,3,1))
      pent(1,11,3,5)=-(FDeffPenLoop(2,3,1,4)+FDeffPenLoop(3,1,4,2)
     &     +FDeffPenLoop(4,2,3,1)+FDeffPenLoop(3,2,4,1))
      IF(PenSingular(1)) then
         PendbleSin=.true.
      endif
      IF(PenSingular(2)) then
         PenquadSin=.true.
      endif
cc start check
c      print*,'pent(1,3,5,11)',pent(1,3,5,11)*dsqrt(fac)*(-1)
c      print*,'pent(1,5,11,3)',pent(1,5,11,3)*dsqrt(fac)!*(-1)
c      print*,'pent(1,11,3,5)',pent(1,11,3,5)*dsqrt(fac)*(-1)
c      print*,''
      if(lgauge2) then
         print*,'sum^2 [ 1,3,5,(2,4=11) ]',
     &        dreal(pent(1,3,5,11)+pent(1,5,11,3)+pent(1,11,3,5))
     &       *dreal(pent(1,3,5,11)+pent(1,5,11,3)+pent(1,11,3,5))
     &       +dimag(pent(1,3,5,11)+pent(1,5,11,3)+pent(1,11,3,5))
     &       *dimag(pent(1,3,5,11)+pent(1,5,11,3)+pent(1,11,3,5))
         print*,''
      endif
cc end check
*------------ [ 1,3,4,(2,5=12) ] ------------
      do mu=0,4
         Kpen(mu,1)=K(mu,1)
         Kpen(mu,2)=K(mu,3)
         Kpen(mu,3)=K(mu,4)
         Kpen(mu,4)=Kij(mu,2,5)
         if(mu.eq.4) exit
         epsgPen(mu,1)=epsgC(mu,1)
         epsgPen(mu,2)=epsgC(mu,3)
         epsgPen(mu,3)=epsgC(mu,4)
         epsgPen(mu,4)=Jij(mu,2,5)
      enddo
c      CALL auxiliary_calculation(epsgPen(0,1),Kpen(0,1))
c      CALL calculate_E_functions() !depend on q_dot_q(i,j) only
      PenSingular(1)=.false.
      PenSingular(2)=.false.
      pent(1,3,4,12)=-(FDeffPenLoop(4,3,2,1)+FDeffPenLoop(2,3,4,1)
     &     +FDeffPenLoop(2,1,4,3)+FDeffPenLoop(3,2,1,4))
      pent(1,4,12,3)=-(FDeffPenLoop(2,1,3,4)+FDeffPenLoop(3,4,2,1)
     &     +FDeffPenLoop(4,2,1,3)+FDeffPenLoop(2,4,3,1))
      pent(1,12,3,4)=-(FDeffPenLoop(2,3,1,4)+FDeffPenLoop(3,1,4,2)
     &     +FDeffPenLoop(4,2,3,1)+FDeffPenLoop(3,2,4,1))
      IF(PenSingular(1)) then
         PendbleSin=.true.
      endif
      IF(PenSingular(2)) then
         PenquadSin=.true.
      endif
cc start check
c      print*,'pent(1,3,4,12)',pent(1,3,4,12)*dsqrt(fac)*(-1)
c      print*,'pent(1,4,12,3)',pent(1,4,12,3)*dsqrt(fac)!*(-1)
c      print*,'pent(1,12,3,4)',pent(1,12,3,4)*dsqrt(fac)*(-1)
c      print*,''
      if(lgauge2) then
         print*,'sum^2 [ 1,3,4,(2,5=12) ]',
     &        dreal(pent(1,3,4,12)+pent(1,4,12,3)+pent(1,12,3,4))
     &       *dreal(pent(1,3,4,12)+pent(1,4,12,3)+pent(1,12,3,4))
     &       +dimag(pent(1,3,4,12)+pent(1,4,12,3)+pent(1,12,3,4))
     &       *dimag(pent(1,3,4,12)+pent(1,4,12,3)+pent(1,12,3,4))
         print*,''
      endif
cc end check
*------------ [ 1,2,5,(3,4=13) ] ------------
      do mu=0,4
         Kpen(mu,1)=K(mu,1)
         Kpen(mu,2)=K(mu,2)
         Kpen(mu,3)=K(mu,5)
         Kpen(mu,4)=Kij(mu,3,4)
         if(mu.eq.4) exit
         epsgPen(mu,1)=epsgC(mu,1)
         epsgPen(mu,2)=epsgC(mu,2)
         epsgPen(mu,3)=epsgC(mu,5)
         epsgPen(mu,4)=Jij(mu,3,4)
      enddo
c      CALL auxiliary_calculation(epsgPen(0,1),Kpen(0,1))
c      CALL calculate_E_functions() !depend on q_dot_q(i,j) only
      PenSingular(1)=.false.
      PenSingular(2)=.false.
      pent(1,2,5,13)=-(FDeffPenLoop(4,3,2,1)+FDeffPenLoop(2,3,4,1)
     &     +FDeffPenLoop(2,1,4,3)+FDeffPenLoop(3,2,1,4))
      pent(1,5,13,2)=-(FDeffPenLoop(2,1,3,4)+FDeffPenLoop(3,4,2,1)
     &     +FDeffPenLoop(4,2,1,3)+FDeffPenLoop(2,4,3,1))
      pent(1,13,2,5)=-(FDeffPenLoop(2,3,1,4)+FDeffPenLoop(3,1,4,2)
     &     +FDeffPenLoop(4,2,3,1)+FDeffPenLoop(3,2,4,1))
      IF(PenSingular(1)) then
         PendbleSin=.true.
      endif
      IF(PenSingular(2)) then
         PenquadSin=.true.
      endif
cc start check
c      print*,'pent(1,2,5,13)',pent(1,2,5,13)*dsqrt(fac)!*(-1)
c      print*,'pent(1,5,13,2)',pent(1,5,13,2)*dsqrt(fac)!*(-1)
c      print*,'pent(1,13,2,5)',pent(1,13,2,5)*dsqrt(fac)*(-1)
c      print*,''
      if(lgauge2) then
         print*,'sum^2 [ 1,2,5,(3,4=13) ]',
     &        dreal(pent(1,2,5,13)+pent(1,5,13,2)+pent(1,13,2,5))
     &       *dreal(pent(1,2,5,13)+pent(1,5,13,2)+pent(1,13,2,5))
     &       +dimag(pent(1,2,5,13)+pent(1,5,13,2)+pent(1,13,2,5))
     &       *dimag(pent(1,2,5,13)+pent(1,5,13,2)+pent(1,13,2,5))
         print*,''
      endif
cc end check
*------------ [ 1,2,4,(3,5=14) ] ------------
      do mu=0,4
         Kpen(mu,1)=K(mu,1)
         Kpen(mu,2)=K(mu,2)
         Kpen(mu,3)=K(mu,4)
         Kpen(mu,4)=Kij(mu,3,5)
         if(mu.eq.4) exit
         epsgPen(mu,1)=epsgC(mu,1)
         epsgPen(mu,2)=epsgC(mu,2)
         epsgPen(mu,3)=epsgC(mu,4)
         epsgPen(mu,4)=Jij(mu,3,5)
      enddo
c      CALL auxiliary_calculation(epsgPen(0,1),Kpen(0,1))
c      CALL calculate_E_functions() !depend on q_dot_q(i,j) only
      PenSingular(1)=.false.
      PenSingular(2)=.false.
      pent(1,2,4,14)=-(FDeffPenLoop(4,3,2,1)+FDeffPenLoop(2,3,4,1)
     &     +FDeffPenLoop(2,1,4,3)+FDeffPenLoop(3,2,1,4))
      pent(1,4,14,2)=-(FDeffPenLoop(2,1,3,4)+FDeffPenLoop(3,4,2,1)
     &     +FDeffPenLoop(4,2,1,3)+FDeffPenLoop(2,4,3,1))
      pent(1,14,2,4)=-(FDeffPenLoop(2,3,1,4)+FDeffPenLoop(3,1,4,2)
     &     +FDeffPenLoop(4,2,3,1)+FDeffPenLoop(3,2,4,1))
      IF(PenSingular(1)) then
         PendbleSin=.true.
      endif
      IF(PenSingular(2)) then
         PenquadSin=.true.
      endif
cc start check
c      print*,'pent(1,2,4,14)',pent(1,2,4,14)*dsqrt(fac)!*(-1)
c      print*,'pent(1,4,14,2)',pent(1,4,14,2)*dsqrt(fac)*(-1)
c      print*,'pent(1,14,2,4)',pent(1,14,2,4)*dsqrt(fac)*(-1)
c      print*,''
      if(lgauge2) then
         print*,'sum^2 [ 1,2,4,(3,5=14) ]',
     &        dreal(pent(1,2,4,14)+pent(1,4,14,2)+pent(1,14,2,4))
     &       *dreal(pent(1,2,4,14)+pent(1,4,14,2)+pent(1,14,2,4))
     &       +dimag(pent(1,2,4,14)+pent(1,4,14,2)+pent(1,14,2,4))
     &       *dimag(pent(1,2,4,14)+pent(1,4,14,2)+pent(1,14,2,4))
         print*,''
      endif
c      STOP
cc end check
      colorPentagons(1) =(+pent(1,2,3,15)+pent(3,4,5,6)-pent(2,3,4,9)
     &            +pent(1,10,4,5)+pent(1,5,13,2)) 
      colorPentagons(2)=(+pent(1,15,2,3)+pent(2,4,5,7)-pent(2,4,9,3)
     &            -pent(1,10,4,5)+pent(1,5,11,3))
      colorPentagons(3)=(+pent(2,3,5,8)-pent(2,4,9,3)+pent(1,5,10,4)
     &            -pent(1,11,3,5)+pent(1,14,2,4))
      colorPentagons(4)=(+pent(2,5,8,3)-pent(2,3,4,9)-pent(1,5,10,4)
     &            +pent(1,12,3,4)-pent(1,13,2,5))
      colorPentagons(5)=(-pent(2,5,8,3)+pent(2,3,4,9)+pent(1,5,10,4)
     &            -pent(1,12,3,4)+pent(1,13,2,5))
      colorPentagons(6)=(-pent(2,3,5,8)+pent(2,4,9,3)-pent(1,5,10,4)
     &            +pent(1,11,3,5)-pent(1,14,2,4))
      colorPentagons(7)=(-pent(1,15,2,3)-pent(2,4,5,7)+pent(2,4,9,3)
     &            +pent(1,10,4,5)-pent(1,5,11,3))
      colorPentagons(8) =(-pent(1,2,3,15)-pent(3,4,5,6)+pent(2,3,4,9)
     &            -pent(1,10,4,5)-pent(1,5,13,2))
      colorPentagons(9)=(+pent(3,5,6,4)-pent(2,8,3,5)+pent(1,12,3,4)
     &            +pent(1,2,5,13)-pent(1,4,14,2))
      colorPentagons(10)=(-pent(1,15,2,3)+pent(2,7,4,5)-pent(2,5,8,3)
     &            -pent(1,4,5,10)+pent(1,4,12,3))
      colorPentagons(11)=(+pent(2,5,7,4)-pent(2,8,3,5)+pent(1,3,5,11)
     &            -pent(1,4,12,3)+pent(1,14,2,4))
      colorPentagons(12)=(-pent(2,5,7,4)+pent(2,8,3,5)-pent(1,3,5,11)
     &            +pent(1,4,12,3)-pent(1,14,2,4))
      colorPentagons(13)=(+pent(1,15,2,3)-pent(2,7,4,5)+pent(2,5,8,3)
     &            +pent(1,4,5,10)-pent(1,4,12,3))
      colorPentagons(14)=(-pent(3,5,6,4)+pent(2,8,3,5)-pent(1,12,3,4)
     &            -pent(1,2,5,13)+pent(1,4,14,2))
      colorPentagons(15) =(-pent(1,2,3,15)+pent(3,6,4,5)-pent(2,3,5,8)
     &            +pent(1,4,5,10)+pent(1,4,14,2))
      colorPentagons(16) =(+pent(1,2,3,15)-pent(3,6,4,5)+pent(2,3,5,8)
     &            -pent(1,4,5,10)-pent(1,4,14,2))
      colorPentagons(17)=(+pent(3,5,6,4)-pent(2,9,3,4)+pent(1,11,3,5)
     &            -pent(1,5,13,2)+pent(1,2,4,14))
      colorPentagons(18) =(+pent(2,5,7,4)-pent(2,9,3,4)-pent(1,5,11,3)
     &            +pent(1,3,4,12)+pent(1,13,2,5))
      colorPentagons(19)=(-pent(2,5,7,4)+pent(2,9,3,4)+pent(1,5,11,3)
     &            -pent(1,3,4,12)-pent(1,13,2,5))
      colorPentagons(20)=(-pent(3,5,6,4)+pent(2,9,3,4)-pent(1,11,3,5)
     &            +pent(1,5,13,2)-pent(1,2,4,14))
      colorPentagons(21) =(-pent(1,3,15,2)+pent(3,4,5,6)-pent(2,7,4,5)
     &            +pent(1,3,4,12)-pent(1,2,5,13))
      colorPentagons(22) =(+pent(1,3,15,2)-pent(3,4,5,6)+pent(2,7,4,5)
     &            -pent(1,3,4,12)+pent(1,2,5,13))
      colorPentagons(23) =(+pent(1,3,15,2)+pent(3,6,4,5)-pent(2,4,5,7)
     &            +pent(1,3,5,11)-pent(1,2,4,14))
      colorPentagons(24) =(-pent(1,3,15,2)-pent(3,6,4,5)+pent(2,4,5,7)
     &            -pent(1,3,5,11)+pent(1,2,4,14))

cc start check
      do i=1,24
         colorN1(i)=colorPentagons(i)*dsqrt(fac)
      enddo
      endif

      if(hexeva) then
c       print*, "hexeva", hexeva
cc end check
        hexdbleSin=.false.
        hexquadSin=.false.
c... Hexagon contributions
c abcde         
      HexSingular(1)=.false.
      HexSingular(2)=.false.
      colorHexagons(1)=FDHexLoop(1,2,3,4,5)+FDHexLoop(2,3,4,5,1)
     &     +FDHexLoop(3,4,5,1,2)+FDHexLoop(4,5,1,2,3)
     &     +FDHexLoop(5,1,2,3,4)
      IF(HexSingular(1)) then
         hexdbleSin=.true.
      endif
      IF(HexSingular(2)) then
         hexquadSin=.true.
      endif
c acbde
      HexSingular(1)=.false.
      HexSingular(2)=.false.
      colorHexagons(2)=FDHexLoop(1,3,2,4,5)+FDHexLoop(3,2,4,5,1)
     &     +FDHexLoop(2,4,5,1,3)+FDHexLoop(4,5,1,3,2)
     &     +FDHexLoop(5,1,3,2,4)
      IF(HexSingular(1)) then
         hexdbleSin=.true.
      endif
      IF(HexSingular(2)) then
         hexquadSin=.true.
      endif
c adbce
      HexSingular(1)=.false.
      HexSingular(2)=.false.
      colorHexagons(3)=FDHexLoop(1,4,2,3,5)+FDHexLoop(4,2,3,5,1)
     &     +FDHexLoop(2,3,5,1,4)+FDHexLoop(3,5,1,4,2)
     &     +FDHexLoop(5,1,4,2,3)
      IF(HexSingular(1)) then
         hexdbleSin=.true.
      endif
      IF(HexSingular(2)) then
         hexquadSin=.true.
      endif
c adcbe
      HexSingular(1)=.false.
      HexSingular(2)=.false.
      colorHexagons(4)=FDHexLoop(1,4,3,2,5)+FDHexLoop(4,3,2,5,1)
     &     +FDHexLoop(3,2,5,1,4)+FDHexLoop(2,5,1,4,3)
     &     +FDHexLoop(5,1,4,3,2)
      IF(HexSingular(1)) then
         hexdbleSin=.true.
      endif
      IF(HexSingular(2)) then
         hexquadSin=.true.
      endif
c aebcd
      colorHexagons(5)=-colorHexagons(4)
c aecbd
      colorHexagons(6)=-colorHexagons(3)
c aedbc
      colorHexagons(7)=-colorHexagons(2)
c aedcb
      colorHexagons(8)=-colorHexagons(1)
c abecd
      HexSingular(1)=.false.
      HexSingular(2)=.false.
      colorHexagons(9)=FDHexLoop(1,2,5,3,4)+FDHexLoop(2,5,3,4,1)
     &     +FDHexLoop(5,3,4,1,2)+FDHexLoop(3,4,1,2,5)
     &     +FDHexLoop(4,1,2,5,3)
      IF(HexSingular(1)) then
         hexdbleSin=.true.
      endif
      IF(HexSingular(2)) then
         hexquadSin=.true.
      endif
c acbed
      HexSingular(1)=.false.
      HexSingular(2)=.false.
      colorHexagons(10)=FDHexLoop(1,3,2,5,4)+FDHexLoop(3,2,5,4,1)
     &     +FDHexLoop(2,5,4,1,3)+FDHexLoop(5,4,1,3,2)
     &     +FDHexLoop(4,1,3,2,5)
      IF(HexSingular(1)) then
         hexdbleSin=.true.
      endif
      IF(HexSingular(2)) then
         hexquadSin=.true.
      endif
c acebd
      HexSingular(1)=.false.
      HexSingular(2)=.false.
      colorHexagons(11)=FDHexLoop(1,3,5,2,4)+FDHexLoop(3,5,2,4,1)
     &     +FDHexLoop(5,2,4,1,3)+FDHexLoop(2,4,1,3,5)
     &     +FDHexLoop(4,1,3,5,2)
      IF(HexSingular(1)) then
         hexdbleSin=.true.
      endif
      IF(HexSingular(2)) then
         hexquadSin=.true.
      endif
c adbec
      colorHexagons(12)=-colorHexagons(11)
c adebc
      colorHexagons(13)=-colorHexagons(10)
c adceb
      colorHexagons(14)=-colorHexagons(9)
c abced
      HexSingular(1)=.false.
      HexSingular(2)=.false.
      colorHexagons(15)=FDHexLoop(1,2,3,5,4)+FDHexLoop(2,3,5,4,1)
     &     +FDHexLoop(3,5,4,1,2)+FDHexLoop(5,4,1,2,3)
     &     +FDHexLoop(4,1,2,3,5)
      IF(HexSingular(1)) then
         hexdbleSin=.true.
      endif
      IF(HexSingular(2)) then
         hexquadSin=.true.
      endif
c adecb
      colorHexagons(16)=-colorHexagons(15)
c abdce
      HexSingular(1)=.false.
      HexSingular(2)=.false.
      colorHexagons(17)=FDHexLoop(1,2,4,3,5)+FDHexLoop(2,4,3,5,1)
     &     +FDHexLoop(4,3,5,1,2)+FDHexLoop(3,5,1,2,4)
     &     +FDHexLoop(5,1,2,4,3)
      IF(HexSingular(1)) then
         hexdbleSin=.true.
      endif
      IF(HexSingular(2)) then
         hexquadSin=.true.
      endif
c acdbe
      HexSingular(1)=.false.
      HexSingular(2)=.false.
      colorHexagons(18)=FDHexLoop(1,3,4,2,5)+FDHexLoop(3,4,2,5,1)
     &     +FDHexLoop(4,2,5,1,3)+FDHexLoop(2,5,1,3,4)
     &     +FDHexLoop(5,1,3,4,2)
      IF(HexSingular(1)) then
         hexdbleSin=.true.
      endif
      IF(HexSingular(2)) then
         hexquadSin=.true.
      endif
c aebdc
      colorHexagons(19)=-colorHexagons(18)
c aecdb
      colorHexagons(20)=-colorHexagons(17)
c abedc
      HexSingular(1)=.false.
      HexSingular(2)=.false.
      colorHexagons(21)=FDHexLoop(1,2,5,4,3)+FDHexLoop(2,5,4,3,1)
     &     +FDHexLoop(5,4,3,1,2)+FDHexLoop(4,3,1,2,5)
     &     +FDHexLoop(3,1,2,5,4)
      IF(HexSingular(1)) then
         hexdbleSin=.true.
      endif
      IF(HexSingular(2)) then
         hexquadSin=.true.
      endif
c acdeb
      colorHexagons(22)=-colorHexagons(21)
c abdec
      HexSingular(1)=.false.
      HexSingular(2)=.false.
      colorHexagons(23)=FDHexLoop(1,2,4,5,3)+FDHexLoop(2,4,5,3,1)
     &     +FDHexLoop(4,5,3,1,2)+FDHexLoop(5,3,1,2,4)
     &     +FDHexLoop(3,1,2,4,5)
      IF(HexSingular(1)) then
         hexdbleSin=.true.
      endif
      IF(HexSingular(2)) then
         hexquadSin=.true.
      endif
c acedb
      colorHexagons(24)=-colorHexagons(23)
      endif

      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
c start check
      if(Mgdebug1) then
         do i=1,24
            print*, "ColorPentagons",i,colorPentagons(i)!/sha
         enddo
         do i=1,24
            print*, "ColorHexagons",i,colorHexagons(i)!/sha
         enddo
      if(Mgpen) then
         do i=1,24
c            print*, "ColorPentagons",i,colorPentagons(i)!/sha
           colorpentagons1(i)=colorpentagons(i)!/shat
           colorPentagons(i)=0d0
c           colorHexagons(i)=0d0
         enddo
      endif

      if(Mghex) then
         shat=dotrr(2d0*K(0,1),2d0*K(0,2))
         do i=1,24
c            print*, "ColorHexagons",i,colorHexagons(i)!/shat
           colorHexagons1(i)=colorHexagons(i)!/shat
           colorHexagons(i)=0d0
         enddo
      endif
      endif
c end check
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      color(1) =2d0*(-amp3j(1)-amp3j(2)-amp3j(3)-amp3j(4)
     &    +amp3j(5)-amp3j(6)+amp3j(7)-amp3j(8)+amp3j(9)
     &    +amp3j(10)+amp3j(11)-amp3j(12)-amp3j(13)+amp3j(14)
     &    -amp3j(15)+amp3j(16)+amp3j(17)-amp3j(18)+amp3j(19)
     &    +amp3j(20)-amp3j(21)-amp3j(22)-amp3j(23)+amp3j(24)
     &    +amp3j(25)+amp3j(26)-amp3j(79)+amp3j(80)-amp3j(81)
     &    +amp3j(82)-amp3j(83)-amp3j(84)-amp3j(85)+amp3j(86)
     &    +amp3j(87)+amp3j(88)+amp3j(89)-amp3j(90)+amp3j(91)
     &    -amp3j(92)-amp3j(93)-amp3j(94)+amp3j(95)+amp3j(96)
     &    -amp3j(97)+amp3j(98)+amp3j(99)-amp3j(100)+amp3j(101)
     &    -amp3j(102)-amp3j(103)+amp3j(104)-amp3j(261)-amp3j(262)
     &    +amp3j(263)+amp3j(264)-amp3j(265)+amp3j(266)+amp3j(267)
     &    +amp3j(268)-amp3j(269)+amp3j(270)-amp3j(271)-amp3j(272)
     &    -amp3j(273)+amp3j(274)+amp3j(275)-amp3j(276)-amp3j(277)
     &    +amp3j(278)+amp3j(279)-amp3j(280)-amp3j(281)+amp3j(282)
     &    -amp3j(283)+amp3j(284)+amp3j(285)-amp3j(286)-amp3j(339)
     &    +amp3j(340)+amp3j(341)-amp3j(342)-amp3j(343)+amp3j(344)
     &    -amp3j(345)+amp3j(346)-amp3j(347)+amp3j(348)+amp3j(349)
     &    -amp3j(350)+amp3j(351)-amp3j(352)-amp3j(353)+amp3j(354)
     &    -amp3j(355)+amp3j(356)-amp3j(357)+amp3j(358)+amp3j(359)
     &    -amp3j(360)+amp3j(361)-amp3j(362)-amp3j(363)+amp3j(364)
     &    -amp3j(365)+amp3j(366)-amp3j(367)-amp3j(368)-amp3j(369)
     &    -amp3j(370)+amp3j(371)+amp3j(372)+amp3j(373)+amp3j(374)
     &    -amp3j(375)+amp3j(376)+amp3j(377)-amp3j(378)+amp3j(379)
     &    +amp3j(380)+amp3j(381)-amp3j(382)-amp3j(383)-amp3j(384)
     &    +amp3j(385)+amp3j(386)+amp3j(387)-amp3j(388)-amp3j(389)
     &    -amp3j(390))+colorPentagons(1)+colorHexagons(1)
      color(2) =2d0*(+amp3j(1)+amp3j(2)+amp3j(3)+amp3j(4)
     &    -amp3j(5)+amp3j(6)-amp3j(7)+amp3j(8)-amp3j(9)
     &    -amp3j(10)-amp3j(11)+amp3j(12)+amp3j(13)-amp3j(14)
     &    +amp3j(15)-amp3j(16)-amp3j(17)+amp3j(18)-amp3j(19)
     &    -amp3j(20)+amp3j(21)+amp3j(22)+amp3j(23)-amp3j(24)
     &    -amp3j(25)-amp3j(26)+amp3j(79)-amp3j(80)+amp3j(81)
     &    -amp3j(82)+amp3j(83)+amp3j(84)+amp3j(85)-amp3j(86)
     &    -amp3j(87)-amp3j(88)-amp3j(89)+amp3j(90)-amp3j(91)
     &    +amp3j(92)+amp3j(93)+amp3j(94)-amp3j(95)-amp3j(96)
     &    +amp3j(97)-amp3j(98)-amp3j(99)+amp3j(100)-amp3j(101)
     &    +amp3j(102)+amp3j(103)-amp3j(104)+amp3j(131)+amp3j(132)
     &    +amp3j(133)+amp3j(134)+amp3j(135)-amp3j(136)+amp3j(137)
     &    -amp3j(138)+amp3j(139)+amp3j(140)-amp3j(141)-amp3j(142)
     &    +amp3j(143)-amp3j(144)+amp3j(145)+amp3j(146)-amp3j(147)
     &    +amp3j(148)-amp3j(149)-amp3j(150)-amp3j(151)+amp3j(152)
     &    -amp3j(153)-amp3j(154)-amp3j(155)-amp3j(156)+amp3j(209)
     &    +amp3j(210)-amp3j(211)-amp3j(212)+amp3j(213)-amp3j(214)
     &    -amp3j(215)-amp3j(216)+amp3j(217)+amp3j(218)+amp3j(219)
     &    -amp3j(220)+amp3j(221)-amp3j(222)-amp3j(223)-amp3j(224)
     &    +amp3j(225)+amp3j(226)-amp3j(227)+amp3j(228)+amp3j(229)
     &    -amp3j(230)+amp3j(231)-amp3j(232)-amp3j(233)+amp3j(234)
     &    +amp3j(235)+amp3j(236)+amp3j(237)+amp3j(238)+amp3j(239)
     &    +amp3j(240)-amp3j(241)-amp3j(242)-amp3j(243)+amp3j(244)
     &    -amp3j(245)+amp3j(246)+amp3j(247)-amp3j(248)+amp3j(249)
     &    -amp3j(250)-amp3j(251)-amp3j(252)-amp3j(253)-amp3j(254)
     &    +amp3j(255)+amp3j(256)+amp3j(257)-amp3j(258)-amp3j(259)
     &    -amp3j(260))+colorPentagons(2)+colorHexagons(2)
      color(3) =2d0*(+amp3j(1)+amp3j(2)+amp3j(3)+amp3j(4)
     &    -amp3j(5)+amp3j(6)-amp3j(7)+amp3j(8)-amp3j(9)
     &    -amp3j(10)-amp3j(11)+amp3j(12)+amp3j(13)-amp3j(14)
     &    +amp3j(15)-amp3j(16)-amp3j(17)+amp3j(18)-amp3j(19)
     &    -amp3j(20)+amp3j(21)+amp3j(22)+amp3j(23)-amp3j(24)
     &    -amp3j(25)-amp3j(26)+amp3j(53)+amp3j(54)-amp3j(55)
     &    -amp3j(56)-amp3j(57)+amp3j(58)-amp3j(59)-amp3j(60)
     &    +amp3j(61)+amp3j(62)-amp3j(63)-amp3j(64)-amp3j(65)
     &    +amp3j(66)+amp3j(67)+amp3j(68)-amp3j(69)+amp3j(70)
     &    +amp3j(71)-amp3j(72)-amp3j(73)+amp3j(74)-amp3j(75)
     &    +amp3j(76)+amp3j(77)-amp3j(78)+amp3j(105)+amp3j(106)
     &    +amp3j(107)+amp3j(108)-amp3j(109)+amp3j(110)-amp3j(111)
     &    +amp3j(112)-amp3j(113)+amp3j(114)-amp3j(115)+amp3j(116)
     &    +amp3j(117)-amp3j(118)+amp3j(119)-amp3j(120)-amp3j(121)
     &    -amp3j(122)-amp3j(123)-amp3j(124)+amp3j(125)+amp3j(126)
     &    +amp3j(127)-amp3j(128)-amp3j(129)-amp3j(130)+amp3j(131)
     &    +amp3j(132)+amp3j(133)+amp3j(134)+amp3j(135)-amp3j(136)
     &    +amp3j(137)-amp3j(138)+amp3j(139)+amp3j(140)-amp3j(141)
     &    -amp3j(142)+amp3j(143)-amp3j(144)+amp3j(145)+amp3j(146)
     &    -amp3j(147)+amp3j(148)-amp3j(149)-amp3j(150)-amp3j(151)
     &    +amp3j(152)-amp3j(153)-amp3j(154)-amp3j(155)-amp3j(156)
     &    +amp3j(183)-amp3j(184)+amp3j(185)-amp3j(186)-amp3j(187)
     &    -amp3j(188)+amp3j(189)+amp3j(190)+amp3j(191)-amp3j(192)
     &    +amp3j(193)-amp3j(194)+amp3j(195)-amp3j(196)+amp3j(197)
     &    -amp3j(198)+amp3j(199)+amp3j(200)-amp3j(201)-amp3j(202)
     &    +amp3j(203)+amp3j(204)+amp3j(205)-amp3j(206)-amp3j(207)
     &    -amp3j(208))+colorPentagons(3)+colorHexagons(3)
      color(4) =2d0*(-amp3j(1)-amp3j(2)-amp3j(3)-amp3j(4)
     &    +amp3j(5)-amp3j(6)+amp3j(7)-amp3j(8)+amp3j(9)
     &    +amp3j(10)+amp3j(11)-amp3j(12)-amp3j(13)+amp3j(14)
     &    -amp3j(15)+amp3j(16)+amp3j(17)-amp3j(18)+amp3j(19)
     &    +amp3j(20)-amp3j(21)-amp3j(22)-amp3j(23)+amp3j(24)
     &    +amp3j(25)+amp3j(26)-amp3j(27)-amp3j(28)-amp3j(29)
     &    -amp3j(30)-amp3j(31)+amp3j(32)-amp3j(33)-amp3j(34)
     &    +amp3j(35)-amp3j(36)+amp3j(37)-amp3j(38)-amp3j(39)
     &    +amp3j(40)+amp3j(41)-amp3j(42)+amp3j(43)+amp3j(44)
     &    +amp3j(45)-amp3j(46)+amp3j(47)+amp3j(48)+amp3j(49)
     &    +amp3j(50)+amp3j(51)-amp3j(52)-amp3j(105)-amp3j(106)
     &    -amp3j(107)-amp3j(108)+amp3j(109)-amp3j(110)+amp3j(111)
     &    -amp3j(112)+amp3j(113)-amp3j(114)+amp3j(115)-amp3j(116)
     &    -amp3j(117)+amp3j(118)-amp3j(119)+amp3j(120)+amp3j(121)
     &    +amp3j(122)+amp3j(123)+amp3j(124)-amp3j(125)-amp3j(126)
     &    -amp3j(127)+amp3j(128)+amp3j(129)+amp3j(130)-amp3j(261)
     &    -amp3j(262)+amp3j(263)+amp3j(264)-amp3j(265)+amp3j(266)
     &    +amp3j(267)+amp3j(268)-amp3j(269)+amp3j(270)-amp3j(271)
     &    -amp3j(272)-amp3j(273)+amp3j(274)+amp3j(275)-amp3j(276)
     &    -amp3j(277)+amp3j(278)+amp3j(279)-amp3j(280)-amp3j(281)
     &    +amp3j(282)-amp3j(283)+amp3j(284)+amp3j(285)-amp3j(286)
     &    -amp3j(287)-amp3j(288)+amp3j(289)+amp3j(290)+amp3j(291)
     &    +amp3j(292)-amp3j(293)-amp3j(294)-amp3j(295)-amp3j(296)
     &    +amp3j(297)-amp3j(298)+amp3j(299)-amp3j(300)+amp3j(301)
     &    +amp3j(302)+amp3j(303)+amp3j(304)-amp3j(305)-amp3j(306)
     &    +amp3j(307)+amp3j(308)+amp3j(309)-amp3j(310)-amp3j(311)
     &    -amp3j(312))+colorPentagons(4)+colorHexagons(4)
      color(5) =2d0*(+amp3j(1)+amp3j(2)+amp3j(3)+amp3j(4)
     &    -amp3j(5)+amp3j(6)-amp3j(7)+amp3j(8)-amp3j(9)
     &    -amp3j(10)-amp3j(11)+amp3j(12)+amp3j(13)-amp3j(14)
     &    +amp3j(15)-amp3j(16)-amp3j(17)+amp3j(18)-amp3j(19)
     &    -amp3j(20)+amp3j(21)+amp3j(22)+amp3j(23)-amp3j(24)
     &    -amp3j(25)-amp3j(26)+amp3j(27)+amp3j(28)+amp3j(29)
     &    +amp3j(30)+amp3j(31)-amp3j(32)+amp3j(33)+amp3j(34)
     &    -amp3j(35)+amp3j(36)-amp3j(37)+amp3j(38)+amp3j(39)
     &    -amp3j(40)-amp3j(41)+amp3j(42)-amp3j(43)-amp3j(44)
     &    -amp3j(45)+amp3j(46)-amp3j(47)-amp3j(48)-amp3j(49)
     &    -amp3j(50)-amp3j(51)+amp3j(52)+amp3j(105)+amp3j(106)
     &    +amp3j(107)+amp3j(108)-amp3j(109)+amp3j(110)-amp3j(111)
     &    +amp3j(112)-amp3j(113)+amp3j(114)-amp3j(115)+amp3j(116)
     &    +amp3j(117)-amp3j(118)+amp3j(119)-amp3j(120)-amp3j(121)
     &    -amp3j(122)-amp3j(123)-amp3j(124)+amp3j(125)+amp3j(126)
     &    +amp3j(127)-amp3j(128)-amp3j(129)-amp3j(130)+amp3j(261)
     &    +amp3j(262)-amp3j(263)-amp3j(264)+amp3j(265)-amp3j(266)
     &    -amp3j(267)-amp3j(268)+amp3j(269)-amp3j(270)+amp3j(271)
     &    +amp3j(272)+amp3j(273)-amp3j(274)-amp3j(275)+amp3j(276)
     &    +amp3j(277)-amp3j(278)-amp3j(279)+amp3j(280)+amp3j(281)
     &    -amp3j(282)+amp3j(283)-amp3j(284)-amp3j(285)+amp3j(286)
     &    +amp3j(287)+amp3j(288)-amp3j(289)-amp3j(290)-amp3j(291)
     &    -amp3j(292)+amp3j(293)+amp3j(294)+amp3j(295)+amp3j(296)
     &    -amp3j(297)+amp3j(298)-amp3j(299)+amp3j(300)-amp3j(301)
     &    -amp3j(302)-amp3j(303)-amp3j(304)+amp3j(305)+amp3j(306)
     &    -amp3j(307)-amp3j(308)-amp3j(309)+amp3j(310)+amp3j(311)
     &    +amp3j(312))+colorPentagons(5)+colorHexagons(5)
      color(6) =2d0*(-amp3j(1)-amp3j(2)-amp3j(3)-amp3j(4)
     &    +amp3j(5)-amp3j(6)+amp3j(7)-amp3j(8)+amp3j(9)
     &    +amp3j(10)+amp3j(11)-amp3j(12)-amp3j(13)+amp3j(14)
     &    -amp3j(15)+amp3j(16)+amp3j(17)-amp3j(18)+amp3j(19)
     &    +amp3j(20)-amp3j(21)-amp3j(22)-amp3j(23)+amp3j(24)
     &    +amp3j(25)+amp3j(26)-amp3j(53)-amp3j(54)+amp3j(55)
     &    +amp3j(56)+amp3j(57)-amp3j(58)+amp3j(59)+amp3j(60)
     &    -amp3j(61)-amp3j(62)+amp3j(63)+amp3j(64)+amp3j(65)
     &    -amp3j(66)-amp3j(67)-amp3j(68)+amp3j(69)-amp3j(70)
     &    -amp3j(71)+amp3j(72)+amp3j(73)-amp3j(74)+amp3j(75)
     &    -amp3j(76)-amp3j(77)+amp3j(78)-amp3j(105)-amp3j(106)
     &    -amp3j(107)-amp3j(108)+amp3j(109)-amp3j(110)+amp3j(111)
     &    -amp3j(112)+amp3j(113)-amp3j(114)+amp3j(115)-amp3j(116)
     &    -amp3j(117)+amp3j(118)-amp3j(119)+amp3j(120)+amp3j(121)
     &    +amp3j(122)+amp3j(123)+amp3j(124)-amp3j(125)-amp3j(126)
     &    -amp3j(127)+amp3j(128)+amp3j(129)+amp3j(130)-amp3j(131)
     &    -amp3j(132)-amp3j(133)-amp3j(134)-amp3j(135)+amp3j(136)
     &    -amp3j(137)+amp3j(138)-amp3j(139)-amp3j(140)+amp3j(141)
     &    +amp3j(142)-amp3j(143)+amp3j(144)-amp3j(145)-amp3j(146)
     &    +amp3j(147)-amp3j(148)+amp3j(149)+amp3j(150)+amp3j(151)
     &    -amp3j(152)+amp3j(153)+amp3j(154)+amp3j(155)+amp3j(156)
     &    -amp3j(183)+amp3j(184)-amp3j(185)+amp3j(186)+amp3j(187)
     &    +amp3j(188)-amp3j(189)-amp3j(190)-amp3j(191)+amp3j(192)
     &    -amp3j(193)+amp3j(194)-amp3j(195)+amp3j(196)-amp3j(197)
     &    +amp3j(198)-amp3j(199)-amp3j(200)+amp3j(201)+amp3j(202)
     &    -amp3j(203)-amp3j(204)-amp3j(205)+amp3j(206)+amp3j(207)
     &    +amp3j(208))+colorPentagons(6)+colorHexagons(6)
      color(7) =2d0*(-amp3j(1)-amp3j(2)-amp3j(3)-amp3j(4)
     &    +amp3j(5)-amp3j(6)+amp3j(7)-amp3j(8)+amp3j(9)
     &    +amp3j(10)+amp3j(11)-amp3j(12)-amp3j(13)+amp3j(14)
     &    -amp3j(15)+amp3j(16)+amp3j(17)-amp3j(18)+amp3j(19)
     &    +amp3j(20)-amp3j(21)-amp3j(22)-amp3j(23)+amp3j(24)
     &    +amp3j(25)+amp3j(26)-amp3j(79)+amp3j(80)-amp3j(81)
     &    +amp3j(82)-amp3j(83)-amp3j(84)-amp3j(85)+amp3j(86)
     &    +amp3j(87)+amp3j(88)+amp3j(89)-amp3j(90)+amp3j(91)
     &    -amp3j(92)-amp3j(93)-amp3j(94)+amp3j(95)+amp3j(96)
     &    -amp3j(97)+amp3j(98)+amp3j(99)-amp3j(100)+amp3j(101)
     &    -amp3j(102)-amp3j(103)+amp3j(104)-amp3j(131)-amp3j(132)
     &    -amp3j(133)-amp3j(134)-amp3j(135)+amp3j(136)-amp3j(137)
     &    +amp3j(138)-amp3j(139)-amp3j(140)+amp3j(141)+amp3j(142)
     &    -amp3j(143)+amp3j(144)-amp3j(145)-amp3j(146)+amp3j(147)
     &    -amp3j(148)+amp3j(149)+amp3j(150)+amp3j(151)-amp3j(152)
     &    +amp3j(153)+amp3j(154)+amp3j(155)+amp3j(156)-amp3j(209)
     &    -amp3j(210)+amp3j(211)+amp3j(212)-amp3j(213)+amp3j(214)
     &    +amp3j(215)+amp3j(216)-amp3j(217)-amp3j(218)-amp3j(219)
     &    +amp3j(220)-amp3j(221)+amp3j(222)+amp3j(223)+amp3j(224)
     &    -amp3j(225)-amp3j(226)+amp3j(227)-amp3j(228)-amp3j(229)
     &    +amp3j(230)-amp3j(231)+amp3j(232)+amp3j(233)-amp3j(234)
     &    -amp3j(235)-amp3j(236)-amp3j(237)-amp3j(238)-amp3j(239)
     &    -amp3j(240)+amp3j(241)+amp3j(242)+amp3j(243)-amp3j(244)
     &    +amp3j(245)-amp3j(246)-amp3j(247)+amp3j(248)-amp3j(249)
     &    +amp3j(250)+amp3j(251)+amp3j(252)+amp3j(253)+amp3j(254)
     &    -amp3j(255)-amp3j(256)-amp3j(257)+amp3j(258)+amp3j(259)
     &    +amp3j(260))+colorPentagons(7)+colorHexagons(7)
      color(8) =2d0*(+amp3j(1)+amp3j(2)+amp3j(3)+amp3j(4)
     &    -amp3j(5)+amp3j(6)-amp3j(7)+amp3j(8)-amp3j(9)
     &    -amp3j(10)-amp3j(11)+amp3j(12)+amp3j(13)-amp3j(14)
     &    +amp3j(15)-amp3j(16)-amp3j(17)+amp3j(18)-amp3j(19)
     &    -amp3j(20)+amp3j(21)+amp3j(22)+amp3j(23)-amp3j(24)
     &    -amp3j(25)-amp3j(26)+amp3j(79)-amp3j(80)+amp3j(81)
     &    -amp3j(82)+amp3j(83)+amp3j(84)+amp3j(85)-amp3j(86)
     &    -amp3j(87)-amp3j(88)-amp3j(89)+amp3j(90)-amp3j(91)
     &    +amp3j(92)+amp3j(93)+amp3j(94)-amp3j(95)-amp3j(96)
     &    +amp3j(97)-amp3j(98)-amp3j(99)+amp3j(100)-amp3j(101)
     &    +amp3j(102)+amp3j(103)-amp3j(104)+amp3j(261)+amp3j(262)
     &    -amp3j(263)-amp3j(264)+amp3j(265)-amp3j(266)-amp3j(267)
     &    -amp3j(268)+amp3j(269)-amp3j(270)+amp3j(271)+amp3j(272)
     &    +amp3j(273)-amp3j(274)-amp3j(275)+amp3j(276)+amp3j(277)
     &    -amp3j(278)-amp3j(279)+amp3j(280)+amp3j(281)-amp3j(282)
     &    +amp3j(283)-amp3j(284)-amp3j(285)+amp3j(286)+amp3j(339)
     &    -amp3j(340)-amp3j(341)+amp3j(342)+amp3j(343)-amp3j(344)
     &    +amp3j(345)-amp3j(346)+amp3j(347)-amp3j(348)-amp3j(349)
     &    +amp3j(350)-amp3j(351)+amp3j(352)+amp3j(353)-amp3j(354)
     &    +amp3j(355)-amp3j(356)+amp3j(357)-amp3j(358)-amp3j(359)
     &    +amp3j(360)-amp3j(361)+amp3j(362)+amp3j(363)-amp3j(364)
     &    +amp3j(365)-amp3j(366)+amp3j(367)+amp3j(368)+amp3j(369)
     &    +amp3j(370)-amp3j(371)-amp3j(372)-amp3j(373)-amp3j(374)
     &    +amp3j(375)-amp3j(376)-amp3j(377)+amp3j(378)-amp3j(379)
     &    -amp3j(380)-amp3j(381)+amp3j(382)+amp3j(383)+amp3j(384)
     &    -amp3j(385)-amp3j(386)-amp3j(387)+amp3j(388)+amp3j(389)
     &    +amp3j(390))+colorPentagons(8)+colorHexagons(8)
      color(9) =2d0*(-amp3j(27)-amp3j(28)-amp3j(29)-amp3j(30)
     &    -amp3j(31)+amp3j(32)-amp3j(33)-amp3j(34)+amp3j(35)
     &    -amp3j(36)+amp3j(37)-amp3j(38)-amp3j(39)+amp3j(40)
     &    +amp3j(41)-amp3j(42)+amp3j(43)+amp3j(44)+amp3j(45)
     &    -amp3j(46)+amp3j(47)+amp3j(48)+amp3j(49)+amp3j(50)
     &    +amp3j(51)-amp3j(52)+amp3j(53)+amp3j(54)-amp3j(55)
     &    -amp3j(56)-amp3j(57)+amp3j(58)-amp3j(59)-amp3j(60)
     &    +amp3j(61)+amp3j(62)-amp3j(63)-amp3j(64)-amp3j(65)
     &    +amp3j(66)+amp3j(67)+amp3j(68)-amp3j(69)+amp3j(70)
     &    +amp3j(71)-amp3j(72)-amp3j(73)+amp3j(74)-amp3j(75)
     &    +amp3j(76)+amp3j(77)-amp3j(78)-amp3j(287)-amp3j(288)
     &    +amp3j(289)+amp3j(290)+amp3j(291)+amp3j(292)-amp3j(293)
     &    -amp3j(294)-amp3j(295)-amp3j(296)+amp3j(297)-amp3j(298)
     &    +amp3j(299)-amp3j(300)+amp3j(301)+amp3j(302)+amp3j(303)
     &    +amp3j(304)-amp3j(305)-amp3j(306)+amp3j(307)+amp3j(308)
     &    +amp3j(309)-amp3j(310)-amp3j(311)-amp3j(312)+amp3j(313)
     &    -amp3j(314)+amp3j(315)+amp3j(316)+amp3j(317)+amp3j(318)
     &    -amp3j(319)-amp3j(320)-amp3j(321)+amp3j(322)+amp3j(323)
     &    -amp3j(324)-amp3j(325)+amp3j(326)-amp3j(327)-amp3j(328)
     &    -amp3j(329)-amp3j(330)+amp3j(331)+amp3j(332)-amp3j(333)
     &    -amp3j(334)-amp3j(335)+amp3j(336)+amp3j(337)+amp3j(338)
     &    +amp3j(365)-amp3j(366)+amp3j(367)+amp3j(368)+amp3j(369)
     &    +amp3j(370)-amp3j(371)-amp3j(372)-amp3j(373)-amp3j(374)
     &    +amp3j(375)-amp3j(376)-amp3j(377)+amp3j(378)-amp3j(379)
     &    -amp3j(380)-amp3j(381)+amp3j(382)+amp3j(383)+amp3j(384)
     &    -amp3j(385)-amp3j(386)-amp3j(387)+amp3j(388)+amp3j(389)
     &    +amp3j(390))+colorPentagons(9)+colorHexagons(9)
      color(10) =2d0*(+amp3j(27)+amp3j(28)+amp3j(29)+amp3j(30)
     &    +amp3j(31)-amp3j(32)+amp3j(33)+amp3j(34)-amp3j(35)
     &    +amp3j(36)-amp3j(37)+amp3j(38)+amp3j(39)-amp3j(40)
     &    -amp3j(41)+amp3j(42)-amp3j(43)-amp3j(44)-amp3j(45)
     &    +amp3j(46)-amp3j(47)-amp3j(48)-amp3j(49)-amp3j(50)
     &    -amp3j(51)+amp3j(52)-amp3j(79)+amp3j(80)-amp3j(81)
     &    +amp3j(82)-amp3j(83)-amp3j(84)-amp3j(85)+amp3j(86)
     &    +amp3j(87)+amp3j(88)+amp3j(89)-amp3j(90)+amp3j(91)
     &    -amp3j(92)-amp3j(93)-amp3j(94)+amp3j(95)+amp3j(96)
     &    -amp3j(97)+amp3j(98)+amp3j(99)-amp3j(100)+amp3j(101)
     &    -amp3j(102)-amp3j(103)+amp3j(104)+amp3j(105)+amp3j(106)
     &    +amp3j(107)+amp3j(108)-amp3j(109)+amp3j(110)-amp3j(111)
     &    +amp3j(112)-amp3j(113)+amp3j(114)-amp3j(115)+amp3j(116)
     &    +amp3j(117)-amp3j(118)+amp3j(119)-amp3j(120)-amp3j(121)
     &    -amp3j(122)-amp3j(123)-amp3j(124)+amp3j(125)+amp3j(126)
     &    +amp3j(127)-amp3j(128)-amp3j(129)-amp3j(130)+amp3j(157)
     &    +amp3j(158)+amp3j(159)+amp3j(160)+amp3j(161)+amp3j(162)
     &    -amp3j(163)-amp3j(164)-amp3j(165)-amp3j(166)-amp3j(167)
     &    +amp3j(168)+amp3j(169)-amp3j(170)+amp3j(171)-amp3j(172)
     &    -amp3j(173)+amp3j(174)-amp3j(175)-amp3j(176)+amp3j(177)
     &    +amp3j(178)+amp3j(179)-amp3j(180)-amp3j(181)-amp3j(182)
     &    -amp3j(209)-amp3j(210)+amp3j(211)+amp3j(212)-amp3j(213)
     &    +amp3j(214)+amp3j(215)+amp3j(216)-amp3j(217)-amp3j(218)
     &    -amp3j(219)+amp3j(220)-amp3j(221)+amp3j(222)+amp3j(223)
     &    +amp3j(224)-amp3j(225)-amp3j(226)+amp3j(227)-amp3j(228)
     &    -amp3j(229)+amp3j(230)-amp3j(231)+amp3j(232)+amp3j(233)
     &    -amp3j(234))+colorPentagons(10)+colorHexagons(10)
      color(11) =2d0*(-amp3j(27)-amp3j(28)-amp3j(29)-amp3j(30)
     &    -amp3j(31)+amp3j(32)-amp3j(33)-amp3j(34)+amp3j(35)
     &    -amp3j(36)+amp3j(37)-amp3j(38)-amp3j(39)+amp3j(40)
     &    +amp3j(41)-amp3j(42)+amp3j(43)+amp3j(44)+amp3j(45)
     &    -amp3j(46)+amp3j(47)+amp3j(48)+amp3j(49)+amp3j(50)
     &    +amp3j(51)-amp3j(52)+amp3j(53)+amp3j(54)-amp3j(55)
     &    -amp3j(56)-amp3j(57)+amp3j(58)-amp3j(59)-amp3j(60)
     &    +amp3j(61)+amp3j(62)-amp3j(63)-amp3j(64)-amp3j(65)
     &    +amp3j(66)+amp3j(67)+amp3j(68)-amp3j(69)+amp3j(70)
     &    +amp3j(71)-amp3j(72)-amp3j(73)+amp3j(74)-amp3j(75)
     &    +amp3j(76)+amp3j(77)-amp3j(78)-amp3j(157)-amp3j(158)
     &    -amp3j(159)-amp3j(160)-amp3j(161)-amp3j(162)+amp3j(163)
     &    +amp3j(164)+amp3j(165)+amp3j(166)+amp3j(167)-amp3j(168)
     &    -amp3j(169)+amp3j(170)-amp3j(171)+amp3j(172)+amp3j(173)
     &    -amp3j(174)+amp3j(175)+amp3j(176)-amp3j(177)-amp3j(178)
     &    -amp3j(179)+amp3j(180)+amp3j(181)+amp3j(182)+amp3j(183)
     &    -amp3j(184)+amp3j(185)-amp3j(186)-amp3j(187)-amp3j(188)
     &    +amp3j(189)+amp3j(190)+amp3j(191)-amp3j(192)+amp3j(193)
     &    -amp3j(194)+amp3j(195)-amp3j(196)+amp3j(197)-amp3j(198)
     &    +amp3j(199)+amp3j(200)-amp3j(201)-amp3j(202)+amp3j(203)
     &    +amp3j(204)+amp3j(205)-amp3j(206)-amp3j(207)-amp3j(208)
     &    -amp3j(235)-amp3j(236)-amp3j(237)-amp3j(238)-amp3j(239)
     &    -amp3j(240)+amp3j(241)+amp3j(242)+amp3j(243)-amp3j(244)
     &    +amp3j(245)-amp3j(246)-amp3j(247)+amp3j(248)-amp3j(249)
     &    +amp3j(250)+amp3j(251)+amp3j(252)+amp3j(253)+amp3j(254)
     &    -amp3j(255)-amp3j(256)-amp3j(257)+amp3j(258)+amp3j(259)
     &    +amp3j(260))+colorPentagons(11)+colorHexagons(11)
      color(12) =2d0*(+amp3j(27)+amp3j(28)+amp3j(29)+amp3j(30)
     &    +amp3j(31)-amp3j(32)+amp3j(33)+amp3j(34)-amp3j(35)
     &    +amp3j(36)-amp3j(37)+amp3j(38)+amp3j(39)-amp3j(40)
     &    -amp3j(41)+amp3j(42)-amp3j(43)-amp3j(44)-amp3j(45)
     &    +amp3j(46)-amp3j(47)-amp3j(48)-amp3j(49)-amp3j(50)
     &    -amp3j(51)+amp3j(52)-amp3j(53)-amp3j(54)+amp3j(55)
     &    +amp3j(56)+amp3j(57)-amp3j(58)+amp3j(59)+amp3j(60)
     &    -amp3j(61)-amp3j(62)+amp3j(63)+amp3j(64)+amp3j(65)
     &    -amp3j(66)-amp3j(67)-amp3j(68)+amp3j(69)-amp3j(70)
     &    -amp3j(71)+amp3j(72)+amp3j(73)-amp3j(74)+amp3j(75)
     &    -amp3j(76)-amp3j(77)+amp3j(78)+amp3j(157)+amp3j(158)
     &    +amp3j(159)+amp3j(160)+amp3j(161)+amp3j(162)-amp3j(163)
     &    -amp3j(164)-amp3j(165)-amp3j(166)-amp3j(167)+amp3j(168)
     &    +amp3j(169)-amp3j(170)+amp3j(171)-amp3j(172)-amp3j(173)
     &    +amp3j(174)-amp3j(175)-amp3j(176)+amp3j(177)+amp3j(178)
     &    +amp3j(179)-amp3j(180)-amp3j(181)-amp3j(182)-amp3j(183)
     &    +amp3j(184)-amp3j(185)+amp3j(186)+amp3j(187)+amp3j(188)
     &    -amp3j(189)-amp3j(190)-amp3j(191)+amp3j(192)-amp3j(193)
     &    +amp3j(194)-amp3j(195)+amp3j(196)-amp3j(197)+amp3j(198)
     &    -amp3j(199)-amp3j(200)+amp3j(201)+amp3j(202)-amp3j(203)
     &    -amp3j(204)-amp3j(205)+amp3j(206)+amp3j(207)+amp3j(208)
     &    +amp3j(235)+amp3j(236)+amp3j(237)+amp3j(238)+amp3j(239)
     &    +amp3j(240)-amp3j(241)-amp3j(242)-amp3j(243)+amp3j(244)
     &    -amp3j(245)+amp3j(246)+amp3j(247)-amp3j(248)+amp3j(249)
     &    -amp3j(250)-amp3j(251)-amp3j(252)-amp3j(253)-amp3j(254)
     &    +amp3j(255)+amp3j(256)+amp3j(257)-amp3j(258)-amp3j(259)
     &    -amp3j(260))+colorPentagons(12)+colorHexagons(12)
      color(13) =2d0*(-amp3j(27)-amp3j(28)-amp3j(29)-amp3j(30)
     &    -amp3j(31)+amp3j(32)-amp3j(33)-amp3j(34)+amp3j(35)
     &    -amp3j(36)+amp3j(37)-amp3j(38)-amp3j(39)+amp3j(40)
     &    +amp3j(41)-amp3j(42)+amp3j(43)+amp3j(44)+amp3j(45)
     &    -amp3j(46)+amp3j(47)+amp3j(48)+amp3j(49)+amp3j(50)
     &    +amp3j(51)-amp3j(52)+amp3j(79)-amp3j(80)+amp3j(81)
     &    -amp3j(82)+amp3j(83)+amp3j(84)+amp3j(85)-amp3j(86)
     &    -amp3j(87)-amp3j(88)-amp3j(89)+amp3j(90)-amp3j(91)
     &    +amp3j(92)+amp3j(93)+amp3j(94)-amp3j(95)-amp3j(96)
     &    +amp3j(97)-amp3j(98)-amp3j(99)+amp3j(100)-amp3j(101)
     &    +amp3j(102)+amp3j(103)-amp3j(104)-amp3j(105)-amp3j(106)
     &    -amp3j(107)-amp3j(108)+amp3j(109)-amp3j(110)+amp3j(111)
     &    -amp3j(112)+amp3j(113)-amp3j(114)+amp3j(115)-amp3j(116)
     &    -amp3j(117)+amp3j(118)-amp3j(119)+amp3j(120)+amp3j(121)
     &    +amp3j(122)+amp3j(123)+amp3j(124)-amp3j(125)-amp3j(126)
     &    -amp3j(127)+amp3j(128)+amp3j(129)+amp3j(130)-amp3j(157)
     &    -amp3j(158)-amp3j(159)-amp3j(160)-amp3j(161)-amp3j(162)
     &    +amp3j(163)+amp3j(164)+amp3j(165)+amp3j(166)+amp3j(167)
     &    -amp3j(168)-amp3j(169)+amp3j(170)-amp3j(171)+amp3j(172)
     &    +amp3j(173)-amp3j(174)+amp3j(175)+amp3j(176)-amp3j(177)
     &    -amp3j(178)-amp3j(179)+amp3j(180)+amp3j(181)+amp3j(182)
     &    +amp3j(209)+amp3j(210)-amp3j(211)-amp3j(212)+amp3j(213)
     &    -amp3j(214)-amp3j(215)-amp3j(216)+amp3j(217)+amp3j(218)
     &    +amp3j(219)-amp3j(220)+amp3j(221)-amp3j(222)-amp3j(223)
     &    -amp3j(224)+amp3j(225)+amp3j(226)-amp3j(227)+amp3j(228)
     &    +amp3j(229)-amp3j(230)+amp3j(231)-amp3j(232)-amp3j(233)
     &    +amp3j(234))+colorPentagons(13)+colorHexagons(13)
      color(14) =2d0*(+amp3j(27)+amp3j(28)+amp3j(29)+amp3j(30)
     &    +amp3j(31)-amp3j(32)+amp3j(33)+amp3j(34)-amp3j(35)
     &    +amp3j(36)-amp3j(37)+amp3j(38)+amp3j(39)-amp3j(40)
     &    -amp3j(41)+amp3j(42)-amp3j(43)-amp3j(44)-amp3j(45)
     &    +amp3j(46)-amp3j(47)-amp3j(48)-amp3j(49)-amp3j(50)
     &    -amp3j(51)+amp3j(52)-amp3j(53)-amp3j(54)+amp3j(55)
     &    +amp3j(56)+amp3j(57)-amp3j(58)+amp3j(59)+amp3j(60)
     &    -amp3j(61)-amp3j(62)+amp3j(63)+amp3j(64)+amp3j(65)
     &    -amp3j(66)-amp3j(67)-amp3j(68)+amp3j(69)-amp3j(70)
     &    -amp3j(71)+amp3j(72)+amp3j(73)-amp3j(74)+amp3j(75)
     &    -amp3j(76)-amp3j(77)+amp3j(78)+amp3j(287)+amp3j(288)
     &    -amp3j(289)-amp3j(290)-amp3j(291)-amp3j(292)+amp3j(293)
     &    +amp3j(294)+amp3j(295)+amp3j(296)-amp3j(297)+amp3j(298)
     &    -amp3j(299)+amp3j(300)-amp3j(301)-amp3j(302)-amp3j(303)
     &    -amp3j(304)+amp3j(305)+amp3j(306)-amp3j(307)-amp3j(308)
     &    -amp3j(309)+amp3j(310)+amp3j(311)+amp3j(312)-amp3j(313)
     &    +amp3j(314)-amp3j(315)-amp3j(316)-amp3j(317)-amp3j(318)
     &    +amp3j(319)+amp3j(320)+amp3j(321)-amp3j(322)-amp3j(323)
     &    +amp3j(324)+amp3j(325)-amp3j(326)+amp3j(327)+amp3j(328)
     &    +amp3j(329)+amp3j(330)-amp3j(331)-amp3j(332)+amp3j(333)
     &    +amp3j(334)+amp3j(335)-amp3j(336)-amp3j(337)-amp3j(338)
     &    -amp3j(365)+amp3j(366)-amp3j(367)-amp3j(368)-amp3j(369)
     &    -amp3j(370)+amp3j(371)+amp3j(372)+amp3j(373)+amp3j(374)
     &    -amp3j(375)+amp3j(376)+amp3j(377)-amp3j(378)+amp3j(379)
     &    +amp3j(380)+amp3j(381)-amp3j(382)-amp3j(383)-amp3j(384)
     &    +amp3j(385)+amp3j(386)+amp3j(387)-amp3j(388)-amp3j(389)
     &    -amp3j(390))+colorPentagons(14)+colorHexagons(14)
      color(15) =2d0*(-amp3j(53)-amp3j(54)+amp3j(55)+amp3j(56)
     &    +amp3j(57)-amp3j(58)+amp3j(59)+amp3j(60)-amp3j(61)
     &    -amp3j(62)+amp3j(63)+amp3j(64)+amp3j(65)-amp3j(66)
     &    -amp3j(67)-amp3j(68)+amp3j(69)-amp3j(70)-amp3j(71)
     &    +amp3j(72)+amp3j(73)-amp3j(74)+amp3j(75)-amp3j(76)
     &    -amp3j(77)+amp3j(78)+amp3j(79)-amp3j(80)+amp3j(81)
     &    -amp3j(82)+amp3j(83)+amp3j(84)+amp3j(85)-amp3j(86)
     &    -amp3j(87)-amp3j(88)-amp3j(89)+amp3j(90)-amp3j(91)
     &    +amp3j(92)+amp3j(93)+amp3j(94)-amp3j(95)-amp3j(96)
     &    +amp3j(97)-amp3j(98)-amp3j(99)+amp3j(100)-amp3j(101)
     &    +amp3j(102)+amp3j(103)-amp3j(104)-amp3j(105)-amp3j(106)
     &    -amp3j(107)-amp3j(108)+amp3j(109)-amp3j(110)+amp3j(111)
     &    -amp3j(112)+amp3j(113)-amp3j(114)+amp3j(115)-amp3j(116)
     &    -amp3j(117)+amp3j(118)-amp3j(119)+amp3j(120)+amp3j(121)
     &    +amp3j(122)+amp3j(123)+amp3j(124)-amp3j(125)-amp3j(126)
     &    -amp3j(127)+amp3j(128)+amp3j(129)+amp3j(130)-amp3j(313)
     &    +amp3j(314)-amp3j(315)-amp3j(316)-amp3j(317)-amp3j(318)
     &    +amp3j(319)+amp3j(320)+amp3j(321)-amp3j(322)-amp3j(323)
     &    +amp3j(324)+amp3j(325)-amp3j(326)+amp3j(327)+amp3j(328)
     &    +amp3j(329)+amp3j(330)-amp3j(331)-amp3j(332)+amp3j(333)
     &    +amp3j(334)+amp3j(335)-amp3j(336)-amp3j(337)-amp3j(338)
     &    +amp3j(339)-amp3j(340)-amp3j(341)+amp3j(342)+amp3j(343)
     &    -amp3j(344)+amp3j(345)-amp3j(346)+amp3j(347)-amp3j(348)
     &    -amp3j(349)+amp3j(350)-amp3j(351)+amp3j(352)+amp3j(353)
     &    -amp3j(354)+amp3j(355)-amp3j(356)+amp3j(357)-amp3j(358)
     &    -amp3j(359)+amp3j(360)-amp3j(361)+amp3j(362)+amp3j(363)
     &    -amp3j(364))+colorPentagons(15)+colorHexagons(15)
      color(16) =2d0*(+amp3j(53)+amp3j(54)-amp3j(55)-amp3j(56)
     &    -amp3j(57)+amp3j(58)-amp3j(59)-amp3j(60)+amp3j(61)
     &    +amp3j(62)-amp3j(63)-amp3j(64)-amp3j(65)+amp3j(66)
     &    +amp3j(67)+amp3j(68)-amp3j(69)+amp3j(70)+amp3j(71)
     &    -amp3j(72)-amp3j(73)+amp3j(74)-amp3j(75)+amp3j(76)
     &    +amp3j(77)-amp3j(78)-amp3j(79)+amp3j(80)-amp3j(81)
     &    +amp3j(82)-amp3j(83)-amp3j(84)-amp3j(85)+amp3j(86)
     &    +amp3j(87)+amp3j(88)+amp3j(89)-amp3j(90)+amp3j(91)
     &    -amp3j(92)-amp3j(93)-amp3j(94)+amp3j(95)+amp3j(96)
     &    -amp3j(97)+amp3j(98)+amp3j(99)-amp3j(100)+amp3j(101)
     &    -amp3j(102)-amp3j(103)+amp3j(104)+amp3j(105)+amp3j(106)
     &    +amp3j(107)+amp3j(108)-amp3j(109)+amp3j(110)-amp3j(111)
     &    +amp3j(112)-amp3j(113)+amp3j(114)-amp3j(115)+amp3j(116)
     &    +amp3j(117)-amp3j(118)+amp3j(119)-amp3j(120)-amp3j(121)
     &    -amp3j(122)-amp3j(123)-amp3j(124)+amp3j(125)+amp3j(126)
     &    +amp3j(127)-amp3j(128)-amp3j(129)-amp3j(130)+amp3j(313)
     &    -amp3j(314)+amp3j(315)+amp3j(316)+amp3j(317)+amp3j(318)
     &    -amp3j(319)-amp3j(320)-amp3j(321)+amp3j(322)+amp3j(323)
     &    -amp3j(324)-amp3j(325)+amp3j(326)-amp3j(327)-amp3j(328)
     &    -amp3j(329)-amp3j(330)+amp3j(331)+amp3j(332)-amp3j(333)
     &    -amp3j(334)-amp3j(335)+amp3j(336)+amp3j(337)+amp3j(338)
     &    -amp3j(339)+amp3j(340)+amp3j(341)-amp3j(342)-amp3j(343)
     &    +amp3j(344)-amp3j(345)+amp3j(346)-amp3j(347)+amp3j(348)
     &    +amp3j(349)-amp3j(350)+amp3j(351)-amp3j(352)-amp3j(353)
     &    +amp3j(354)-amp3j(355)+amp3j(356)-amp3j(357)+amp3j(358)
     &    +amp3j(359)-amp3j(360)+amp3j(361)-amp3j(362)-amp3j(363)
     &    +amp3j(364))+colorPentagons(16)+colorHexagons(16)
      color(17) =2d0*(-amp3j(131)-amp3j(132)-amp3j(133)-amp3j(134)
     &    -amp3j(135)+amp3j(136)-amp3j(137)+amp3j(138)-amp3j(139)
     &    -amp3j(140)+amp3j(141)+amp3j(142)-amp3j(143)+amp3j(144)
     &    -amp3j(145)-amp3j(146)+amp3j(147)-amp3j(148)+amp3j(149)
     &    +amp3j(150)+amp3j(151)-amp3j(152)+amp3j(153)+amp3j(154)
     &    +amp3j(155)+amp3j(156)-amp3j(183)+amp3j(184)-amp3j(185)
     &    +amp3j(186)+amp3j(187)+amp3j(188)-amp3j(189)-amp3j(190)
     &    -amp3j(191)+amp3j(192)-amp3j(193)+amp3j(194)-amp3j(195)
     &    +amp3j(196)-amp3j(197)+amp3j(198)-amp3j(199)-amp3j(200)
     &    +amp3j(201)+amp3j(202)-amp3j(203)-amp3j(204)-amp3j(205)
     &    +amp3j(206)+amp3j(207)+amp3j(208)+amp3j(261)+amp3j(262)
     &    -amp3j(263)-amp3j(264)+amp3j(265)-amp3j(266)-amp3j(267)
     &    -amp3j(268)+amp3j(269)-amp3j(270)+amp3j(271)+amp3j(272)
     &    +amp3j(273)-amp3j(274)-amp3j(275)+amp3j(276)+amp3j(277)
     &    -amp3j(278)-amp3j(279)+amp3j(280)+amp3j(281)-amp3j(282)
     &    +amp3j(283)-amp3j(284)-amp3j(285)+amp3j(286)+amp3j(313)
     &    -amp3j(314)+amp3j(315)+amp3j(316)+amp3j(317)+amp3j(318)
     &    -amp3j(319)-amp3j(320)-amp3j(321)+amp3j(322)+amp3j(323)
     &    -amp3j(324)-amp3j(325)+amp3j(326)-amp3j(327)-amp3j(328)
     &    -amp3j(329)-amp3j(330)+amp3j(331)+amp3j(332)-amp3j(333)
     &    -amp3j(334)-amp3j(335)+amp3j(336)+amp3j(337)+amp3j(338)
     &    +amp3j(365)-amp3j(366)+amp3j(367)+amp3j(368)+amp3j(369)
     &    +amp3j(370)-amp3j(371)-amp3j(372)-amp3j(373)-amp3j(374)
     &    +amp3j(375)-amp3j(376)-amp3j(377)+amp3j(378)-amp3j(379)
     &    -amp3j(380)-amp3j(381)+amp3j(382)+amp3j(383)+amp3j(384)
     &    -amp3j(385)-amp3j(386)-amp3j(387)+amp3j(388)+amp3j(389)
     &    +amp3j(390))+colorPentagons(17)+colorHexagons(17)
      color(18) =2d0*(-amp3j(131)-amp3j(132)-amp3j(133)-amp3j(134)
     &    -amp3j(135)+amp3j(136)-amp3j(137)+amp3j(138)-amp3j(139)
     &    -amp3j(140)+amp3j(141)+amp3j(142)-amp3j(143)+amp3j(144)
     &    -amp3j(145)-amp3j(146)+amp3j(147)-amp3j(148)+amp3j(149)
     &    +amp3j(150)+amp3j(151)-amp3j(152)+amp3j(153)+amp3j(154)
     &    +amp3j(155)+amp3j(156)-amp3j(157)-amp3j(158)-amp3j(159)
     &    -amp3j(160)-amp3j(161)-amp3j(162)+amp3j(163)+amp3j(164)
     &    +amp3j(165)+amp3j(166)+amp3j(167)-amp3j(168)-amp3j(169)
     &    +amp3j(170)-amp3j(171)+amp3j(172)+amp3j(173)-amp3j(174)
     &    +amp3j(175)+amp3j(176)-amp3j(177)-amp3j(178)-amp3j(179)
     &    +amp3j(180)+amp3j(181)+amp3j(182)-amp3j(235)-amp3j(236)
     &    -amp3j(237)-amp3j(238)-amp3j(239)-amp3j(240)+amp3j(241)
     &    +amp3j(242)+amp3j(243)-amp3j(244)+amp3j(245)-amp3j(246)
     &    -amp3j(247)+amp3j(248)-amp3j(249)+amp3j(250)+amp3j(251)
     &    +amp3j(252)+amp3j(253)+amp3j(254)-amp3j(255)-amp3j(256)
     &    -amp3j(257)+amp3j(258)+amp3j(259)+amp3j(260)+amp3j(261)
     &    +amp3j(262)-amp3j(263)-amp3j(264)+amp3j(265)-amp3j(266)
     &    -amp3j(267)-amp3j(268)+amp3j(269)-amp3j(270)+amp3j(271)
     &    +amp3j(272)+amp3j(273)-amp3j(274)-amp3j(275)+amp3j(276)
     &    +amp3j(277)-amp3j(278)-amp3j(279)+amp3j(280)+amp3j(281)
     &    -amp3j(282)+amp3j(283)-amp3j(284)-amp3j(285)+amp3j(286)
     &    +amp3j(287)+amp3j(288)-amp3j(289)-amp3j(290)-amp3j(291)
     &    -amp3j(292)+amp3j(293)+amp3j(294)+amp3j(295)+amp3j(296)
     &    -amp3j(297)+amp3j(298)-amp3j(299)+amp3j(300)-amp3j(301)
     &    -amp3j(302)-amp3j(303)-amp3j(304)+amp3j(305)+amp3j(306)
     &    -amp3j(307)-amp3j(308)-amp3j(309)+amp3j(310)+amp3j(311)
     &    +amp3j(312))+colorPentagons(18)+colorHexagons(18)
      color(19) =2d0*(+amp3j(131)+amp3j(132)+amp3j(133)+amp3j(134)
     &    +amp3j(135)-amp3j(136)+amp3j(137)-amp3j(138)+amp3j(139)
     &    +amp3j(140)-amp3j(141)-amp3j(142)+amp3j(143)-amp3j(144)
     &    +amp3j(145)+amp3j(146)-amp3j(147)+amp3j(148)-amp3j(149)
     &    -amp3j(150)-amp3j(151)+amp3j(152)-amp3j(153)-amp3j(154)
     &    -amp3j(155)-amp3j(156)+amp3j(157)+amp3j(158)+amp3j(159)
     &    +amp3j(160)+amp3j(161)+amp3j(162)-amp3j(163)-amp3j(164)
     &    -amp3j(165)-amp3j(166)-amp3j(167)+amp3j(168)+amp3j(169)
     &    -amp3j(170)+amp3j(171)-amp3j(172)-amp3j(173)+amp3j(174)
     &    -amp3j(175)-amp3j(176)+amp3j(177)+amp3j(178)+amp3j(179)
     &    -amp3j(180)-amp3j(181)-amp3j(182)+amp3j(235)+amp3j(236)
     &    +amp3j(237)+amp3j(238)+amp3j(239)+amp3j(240)-amp3j(241)
     &    -amp3j(242)-amp3j(243)+amp3j(244)-amp3j(245)+amp3j(246)
     &    +amp3j(247)-amp3j(248)+amp3j(249)-amp3j(250)-amp3j(251)
     &    -amp3j(252)-amp3j(253)-amp3j(254)+amp3j(255)+amp3j(256)
     &    +amp3j(257)-amp3j(258)-amp3j(259)-amp3j(260)-amp3j(261)
     &    -amp3j(262)+amp3j(263)+amp3j(264)-amp3j(265)+amp3j(266)
     &    +amp3j(267)+amp3j(268)-amp3j(269)+amp3j(270)-amp3j(271)
     &    -amp3j(272)-amp3j(273)+amp3j(274)+amp3j(275)-amp3j(276)
     &    -amp3j(277)+amp3j(278)+amp3j(279)-amp3j(280)-amp3j(281)
     &    +amp3j(282)-amp3j(283)+amp3j(284)+amp3j(285)-amp3j(286)
     &    -amp3j(287)-amp3j(288)+amp3j(289)+amp3j(290)+amp3j(291)
     &    +amp3j(292)-amp3j(293)-amp3j(294)-amp3j(295)-amp3j(296)
     &    +amp3j(297)-amp3j(298)+amp3j(299)-amp3j(300)+amp3j(301)
     &    +amp3j(302)+amp3j(303)+amp3j(304)-amp3j(305)-amp3j(306)
     &    +amp3j(307)+amp3j(308)+amp3j(309)-amp3j(310)-amp3j(311)
     &    -amp3j(312))+colorPentagons(19)+colorHexagons(19)
      color(20) =2d0*(+amp3j(131)+amp3j(132)+amp3j(133)+amp3j(134)
     &    +amp3j(135)-amp3j(136)+amp3j(137)-amp3j(138)+amp3j(139)
     &    +amp3j(140)-amp3j(141)-amp3j(142)+amp3j(143)-amp3j(144)
     &    +amp3j(145)+amp3j(146)-amp3j(147)+amp3j(148)-amp3j(149)
     &    -amp3j(150)-amp3j(151)+amp3j(152)-amp3j(153)-amp3j(154)
     &    -amp3j(155)-amp3j(156)+amp3j(183)-amp3j(184)+amp3j(185)
     &    -amp3j(186)-amp3j(187)-amp3j(188)+amp3j(189)+amp3j(190)
     &    +amp3j(191)-amp3j(192)+amp3j(193)-amp3j(194)+amp3j(195)
     &    -amp3j(196)+amp3j(197)-amp3j(198)+amp3j(199)+amp3j(200)
     &    -amp3j(201)-amp3j(202)+amp3j(203)+amp3j(204)+amp3j(205)
     &    -amp3j(206)-amp3j(207)-amp3j(208)-amp3j(261)-amp3j(262)
     &    +amp3j(263)+amp3j(264)-amp3j(265)+amp3j(266)+amp3j(267)
     &    +amp3j(268)-amp3j(269)+amp3j(270)-amp3j(271)-amp3j(272)
     &    -amp3j(273)+amp3j(274)+amp3j(275)-amp3j(276)-amp3j(277)
     &    +amp3j(278)+amp3j(279)-amp3j(280)-amp3j(281)+amp3j(282)
     &    -amp3j(283)+amp3j(284)+amp3j(285)-amp3j(286)-amp3j(313)
     &    +amp3j(314)-amp3j(315)-amp3j(316)-amp3j(317)-amp3j(318)
     &    +amp3j(319)+amp3j(320)+amp3j(321)-amp3j(322)-amp3j(323)
     &    +amp3j(324)+amp3j(325)-amp3j(326)+amp3j(327)+amp3j(328)
     &    +amp3j(329)+amp3j(330)-amp3j(331)-amp3j(332)+amp3j(333)
     &    +amp3j(334)+amp3j(335)-amp3j(336)-amp3j(337)-amp3j(338)
     &    -amp3j(365)+amp3j(366)-amp3j(367)-amp3j(368)-amp3j(369)
     &    -amp3j(370)+amp3j(371)+amp3j(372)+amp3j(373)+amp3j(374)
     &    -amp3j(375)+amp3j(376)+amp3j(377)-amp3j(378)+amp3j(379)
     &    +amp3j(380)+amp3j(381)-amp3j(382)-amp3j(383)-amp3j(384)
     &    +amp3j(385)+amp3j(386)+amp3j(387)-amp3j(388)-amp3j(389)
     &    -amp3j(390))+colorPentagons(20)+colorHexagons(20)
      color(21) =2d0*(-amp3j(157)-amp3j(158)-amp3j(159)-amp3j(160)
     &    -amp3j(161)-amp3j(162)+amp3j(163)+amp3j(164)+amp3j(165)
     &    +amp3j(166)+amp3j(167)-amp3j(168)-amp3j(169)+amp3j(170)
     &    -amp3j(171)+amp3j(172)+amp3j(173)-amp3j(174)+amp3j(175)
     &    +amp3j(176)-amp3j(177)-amp3j(178)-amp3j(179)+amp3j(180)
     &    +amp3j(181)+amp3j(182)+amp3j(209)+amp3j(210)-amp3j(211)
     &    -amp3j(212)+amp3j(213)-amp3j(214)-amp3j(215)-amp3j(216)
     &    +amp3j(217)+amp3j(218)+amp3j(219)-amp3j(220)+amp3j(221)
     &    -amp3j(222)-amp3j(223)-amp3j(224)+amp3j(225)+amp3j(226)
     &    -amp3j(227)+amp3j(228)+amp3j(229)-amp3j(230)+amp3j(231)
     &    -amp3j(232)-amp3j(233)+amp3j(234)+amp3j(287)+amp3j(288)
     &    -amp3j(289)-amp3j(290)-amp3j(291)-amp3j(292)+amp3j(293)
     &    +amp3j(294)+amp3j(295)+amp3j(296)-amp3j(297)+amp3j(298)
     &    -amp3j(299)+amp3j(300)-amp3j(301)-amp3j(302)-amp3j(303)
     &    -amp3j(304)+amp3j(305)+amp3j(306)-amp3j(307)-amp3j(308)
     &    -amp3j(309)+amp3j(310)+amp3j(311)+amp3j(312)-amp3j(339)
     &    +amp3j(340)+amp3j(341)-amp3j(342)-amp3j(343)+amp3j(344)
     &    -amp3j(345)+amp3j(346)-amp3j(347)+amp3j(348)+amp3j(349)
     &    -amp3j(350)+amp3j(351)-amp3j(352)-amp3j(353)+amp3j(354)
     &    -amp3j(355)+amp3j(356)-amp3j(357)+amp3j(358)+amp3j(359)
     &    -amp3j(360)+amp3j(361)-amp3j(362)-amp3j(363)+amp3j(364)
     &    -amp3j(365)+amp3j(366)-amp3j(367)-amp3j(368)-amp3j(369)
     &    -amp3j(370)+amp3j(371)+amp3j(372)+amp3j(373)+amp3j(374)
     &    -amp3j(375)+amp3j(376)+amp3j(377)-amp3j(378)+amp3j(379)
     &    +amp3j(380)+amp3j(381)-amp3j(382)-amp3j(383)-amp3j(384)
     &    +amp3j(385)+amp3j(386)+amp3j(387)-amp3j(388)-amp3j(389)
     &    -amp3j(390))+colorPentagons(21)+colorHexagons(21)
      color(22) =2d0*(+amp3j(157)+amp3j(158)+amp3j(159)+amp3j(160)
     &    +amp3j(161)+amp3j(162)-amp3j(163)-amp3j(164)-amp3j(165)
     &    -amp3j(166)-amp3j(167)+amp3j(168)+amp3j(169)-amp3j(170)
     &    +amp3j(171)-amp3j(172)-amp3j(173)+amp3j(174)-amp3j(175)
     &    -amp3j(176)+amp3j(177)+amp3j(178)+amp3j(179)-amp3j(180)
     &    -amp3j(181)-amp3j(182)-amp3j(209)-amp3j(210)+amp3j(211)
     &    +amp3j(212)-amp3j(213)+amp3j(214)+amp3j(215)+amp3j(216)
     &    -amp3j(217)-amp3j(218)-amp3j(219)+amp3j(220)-amp3j(221)
     &    +amp3j(222)+amp3j(223)+amp3j(224)-amp3j(225)-amp3j(226)
     &    +amp3j(227)-amp3j(228)-amp3j(229)+amp3j(230)-amp3j(231)
     &    +amp3j(232)+amp3j(233)-amp3j(234)-amp3j(287)-amp3j(288)
     &    +amp3j(289)+amp3j(290)+amp3j(291)+amp3j(292)-amp3j(293)
     &    -amp3j(294)-amp3j(295)-amp3j(296)+amp3j(297)-amp3j(298)
     &    +amp3j(299)-amp3j(300)+amp3j(301)+amp3j(302)+amp3j(303)
     &    +amp3j(304)-amp3j(305)-amp3j(306)+amp3j(307)+amp3j(308)
     &    +amp3j(309)-amp3j(310)-amp3j(311)-amp3j(312)+amp3j(339)
     &    -amp3j(340)-amp3j(341)+amp3j(342)+amp3j(343)-amp3j(344)
     &    +amp3j(345)-amp3j(346)+amp3j(347)-amp3j(348)-amp3j(349)
     &    +amp3j(350)-amp3j(351)+amp3j(352)+amp3j(353)-amp3j(354)
     &    +amp3j(355)-amp3j(356)+amp3j(357)-amp3j(358)-amp3j(359)
     &    +amp3j(360)-amp3j(361)+amp3j(362)+amp3j(363)-amp3j(364)
     &    +amp3j(365)-amp3j(366)+amp3j(367)+amp3j(368)+amp3j(369)
     &    +amp3j(370)-amp3j(371)-amp3j(372)-amp3j(373)-amp3j(374)
     &    +amp3j(375)-amp3j(376)-amp3j(377)+amp3j(378)-amp3j(379)
     &    -amp3j(380)-amp3j(381)+amp3j(382)+amp3j(383)+amp3j(384)
     &    -amp3j(385)-amp3j(386)-amp3j(387)+amp3j(388)+amp3j(389)
     &    +amp3j(390))+colorPentagons(22)+colorHexagons(22)
      color(23) =2d0*(+amp3j(183)-amp3j(184)+amp3j(185)-amp3j(186)
     &    -amp3j(187)-amp3j(188)+amp3j(189)+amp3j(190)+amp3j(191)
     &    -amp3j(192)+amp3j(193)-amp3j(194)+amp3j(195)-amp3j(196)
     &    +amp3j(197)-amp3j(198)+amp3j(199)+amp3j(200)-amp3j(201)
     &    -amp3j(202)+amp3j(203)+amp3j(204)+amp3j(205)-amp3j(206)
     &    -amp3j(207)-amp3j(208)-amp3j(209)-amp3j(210)+amp3j(211)
     &    +amp3j(212)-amp3j(213)+amp3j(214)+amp3j(215)+amp3j(216)
     &    -amp3j(217)-amp3j(218)-amp3j(219)+amp3j(220)-amp3j(221)
     &    +amp3j(222)+amp3j(223)+amp3j(224)-amp3j(225)-amp3j(226)
     &    +amp3j(227)-amp3j(228)-amp3j(229)+amp3j(230)-amp3j(231)
     &    +amp3j(232)+amp3j(233)-amp3j(234)-amp3j(235)-amp3j(236)
     &    -amp3j(237)-amp3j(238)-amp3j(239)-amp3j(240)+amp3j(241)
     &    +amp3j(242)+amp3j(243)-amp3j(244)+amp3j(245)-amp3j(246)
     &    -amp3j(247)+amp3j(248)-amp3j(249)+amp3j(250)+amp3j(251)
     &    +amp3j(252)+amp3j(253)+amp3j(254)-amp3j(255)-amp3j(256)
     &    -amp3j(257)+amp3j(258)+amp3j(259)+amp3j(260)-amp3j(313)
     &    +amp3j(314)-amp3j(315)-amp3j(316)-amp3j(317)-amp3j(318)
     &    +amp3j(319)+amp3j(320)+amp3j(321)-amp3j(322)-amp3j(323)
     &    +amp3j(324)+amp3j(325)-amp3j(326)+amp3j(327)+amp3j(328)
     &    +amp3j(329)+amp3j(330)-amp3j(331)-amp3j(332)+amp3j(333)
     &    +amp3j(334)+amp3j(335)-amp3j(336)-amp3j(337)-amp3j(338)
     &    +amp3j(339)-amp3j(340)-amp3j(341)+amp3j(342)+amp3j(343)
     &    -amp3j(344)+amp3j(345)-amp3j(346)+amp3j(347)-amp3j(348)
     &    -amp3j(349)+amp3j(350)-amp3j(351)+amp3j(352)+amp3j(353)
     &    -amp3j(354)+amp3j(355)-amp3j(356)+amp3j(357)-amp3j(358)
     &    -amp3j(359)+amp3j(360)-amp3j(361)+amp3j(362)+amp3j(363)
     &    -amp3j(364))+colorPentagons(23)+colorHexagons(23)
      color(24) =2d0*(-amp3j(183)+amp3j(184)-amp3j(185)+amp3j(186)
     &    +amp3j(187)+amp3j(188)-amp3j(189)-amp3j(190)-amp3j(191)
     &    +amp3j(192)-amp3j(193)+amp3j(194)-amp3j(195)+amp3j(196)
     &    -amp3j(197)+amp3j(198)-amp3j(199)-amp3j(200)+amp3j(201)
     &    +amp3j(202)-amp3j(203)-amp3j(204)-amp3j(205)+amp3j(206)
     &    +amp3j(207)+amp3j(208)+amp3j(209)+amp3j(210)-amp3j(211)
     &    -amp3j(212)+amp3j(213)-amp3j(214)-amp3j(215)-amp3j(216)
     &    +amp3j(217)+amp3j(218)+amp3j(219)-amp3j(220)+amp3j(221)
     &    -amp3j(222)-amp3j(223)-amp3j(224)+amp3j(225)+amp3j(226)
     &    -amp3j(227)+amp3j(228)+amp3j(229)-amp3j(230)+amp3j(231)
     &    -amp3j(232)-amp3j(233)+amp3j(234)+amp3j(235)+amp3j(236)
     &    +amp3j(237)+amp3j(238)+amp3j(239)+amp3j(240)-amp3j(241)
     &    -amp3j(242)-amp3j(243)+amp3j(244)-amp3j(245)+amp3j(246)
     &    +amp3j(247)-amp3j(248)+amp3j(249)-amp3j(250)-amp3j(251)
     &    -amp3j(252)-amp3j(253)-amp3j(254)+amp3j(255)+amp3j(256)
     &    +amp3j(257)-amp3j(258)-amp3j(259)-amp3j(260)+amp3j(313)
     &    -amp3j(314)+amp3j(315)+amp3j(316)+amp3j(317)+amp3j(318)
     &    -amp3j(319)-amp3j(320)-amp3j(321)+amp3j(322)+amp3j(323)
     &    -amp3j(324)-amp3j(325)+amp3j(326)-amp3j(327)-amp3j(328)
     &    -amp3j(329)-amp3j(330)+amp3j(331)+amp3j(332)-amp3j(333)
     &    -amp3j(334)-amp3j(335)+amp3j(336)+amp3j(337)+amp3j(338)
     &    -amp3j(339)+amp3j(340)+amp3j(341)-amp3j(342)-amp3j(343)
     &    +amp3j(344)-amp3j(345)+amp3j(346)-amp3j(347)+amp3j(348)
     &    +amp3j(349)-amp3j(350)+amp3j(351)-amp3j(352)-amp3j(353)
     &    +amp3j(354)-amp3j(355)+amp3j(356)-amp3j(357)+amp3j(358)
     &    +amp3j(359)-amp3j(360)+amp3j(361)-amp3j(362)-amp3j(363)
     &    +amp3j(364))+colorPentagons(24)+colorHexagons(24)

      RH=0d0
      RHsum =0d0

      do i=1,24
         colorT=DCMPLX(0d0,0d0)
         do j=1,24
            colorT=colorT+RGB(i,j)*color(j)
         enddo
         RH=RH+(dreal(colorT)*dreal(color(i))
     &        + dimag(colorT)*dimag(color(i)))/108d0   
      enddo
      RHsum=RH*fac

 

c      pause
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
c 
c         Check whether there is instabilities
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!                                                       !! 
c!         IF SO ::::-->   RSUM -> 0                     !!
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
        RHsum=0d0
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
        RHsum=0d0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccc
       else
c good QUAD 
        boxinst(gf_iteration,2,2)=boxinst(gf_iteration,2,2)+1
        endif
        endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c   Check Penagons
c     
        if(peneva) then
c   Dble precision
       If(pendbleSin) then
c bad dble
        peninst(gf_iteration,1,1)=peninst(gf_iteration,1,1)+1
       else
c good dble 
        peninst(gf_iteration,2,1)=peninst(gf_iteration,2,1)+1
        endif

       If(penQUADSin) then
c bad QUAD
        peninst(gf_iteration,1,2)=peninst(gf_iteration,1,2)+1
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        RHsum=0d0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccc
       else
c good QUAD 
        peninst(gf_iteration,2,2)=peninst(gf_iteration,2,2)+1
        endif
        endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c   Check Hexagons
c     
c   Dble precision
        if(hexeva) then
       If(hexdbleSin) then
c bad dble
        hexinst(gf_iteration,1,1)=hexinst(gf_iteration,1,1)+1
       else
c good dble 
        hexinst(gf_iteration,2,1)=hexinst(gf_iteration,2,1)+1
        endif

       If(hexQUADSin) then
c bad QUAD
        hexinst(gf_iteration,1,2)=hexinst(gf_iteration,1,2)+1
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        RHsum=0d0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccc
       else
c good QUAD 
        hexinst(gf_iteration,2,2)=hexinst(gf_iteration,2,2)+1
        endif
        endif
        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc start check
c... gauge check
        if(mgdebug1) then
           shat=dotrr(2d0*K(0,1),2d0*K(0,2))
           print*,''
           do i=1,24
cfc comment out the shat to compare with the hexagons 
              print*,'color mod(',i,')',color(i) !/shat
cfc end fc 
           enddo

           do i=1,24
cfc comment out the shat to compare with the hexagons 
              print*,'color -penta -hex (',i,')',color(i)- colorpentagons(i)- colorhexagons(i) !/shat
cfc end fc 
           enddo
         
           if(Mgpen) then
              do i=1,24
cfc comment out the shat to compare with the hexagons 
                 print*,'color/ pentagons  mod(',i,')',color(i)/colorpentagons1(i)
cfc end fc 
              enddo
           endif


           if(Mghex) then
              do i=1,24
cfc comment out the shat to compare with the hexagons 
                 print*,'color/ hexaxons  mod(',i,')',color(i)/colorhexagons1(i)
cfc end fc 
              enddo
           endif
        endif
        if(mgdebug) then
           print*,'shat',shat
           print*,'RHsum    ',RHsum
c         print*,'RHsum mod',RHsum/(abs(amp3j(4))**2*fac)
           print*,'RHsum mod',RHsum/shat
        endif
        
        
        if(MGdebug1) then
           print*, "FUll"
           do i =1,390
              print*, "amp",i,dsqrt(fac)*Amp3j(i)
           enddo   
        endif
           
        if(mgdebug) then
*--------------- CHECK with MadGraph ---------------- 
 
           do i=1,390
              GGFLOamp3(i)=dsqrt(fac)*amp3j(i)
c            print*,'GGFLOamp3(',i,')',GGFLOamp3(i)
           enddo
           do i=1,24
              GGFLOcolor(i)=color(i)*dsqrt(fac)
c            print*,'color(',i,')',color(i)*dsqrt(fac)
           enddo
           do i=1,24
c            print*,'colorPentagons(',i,')',colorPentagons(i)*dsqrt(fac)
           enddo
c         print*,'RHsum,fac',RHsum,fac
           GGFLOampSq3=RHsum
c         STOP
        endif
*----------------------------------------------------  

       if(print_amp) then
         print*, "AMP",RHSUM
      endif

cc start check
c      print*, "AMP",RHSUM
c      STOP
cc end check

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      END



