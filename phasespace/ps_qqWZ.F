c
c phase space generator for pp --> V1 + V2 + 1 Jet
c                                  V1 --> l1 l2, V2 --> l3 l4
c
c       Franziska Schissler: schissler@particle.uni-karlruhe.de
c       Initial version:  April 2011
c


      subroutine def_ps_qqWZ(bos1,bos2)
          use globalvars, only: ldoblha


      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/cuts.inc"

      integer bos(2), i, bos1, bos2

      REAL*8 RM2,RMG,RM2MIN,RM2MAX, S,YJMIN,YJMAX,PTJMIN
      REAL*8 RM2MINX2, RM2MAXX2, M2MIN, RM2MINX1, RM2MAXX1
      COMMON /BLIPSIWZ/ RM2(2),RMG(2),RM2MIN(2),RM2MAX(2),
     1                RM2MINX1, RM2MAXX1, M2MIN, RM2MINX2, RM2MAXX2,
     2                S,YJMIN,YJMAX,PTJMIN

      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)

      s = ecm**2

      bos(1) = bos1
      bos(2) = bos2

      do i = 1,2
         rm2(i) = xm2(bos(i))        ! mass squared
         rmg(i) = xmg(bos(i))        ! mass * decay width
         if (bos(i).eq.2) then
            rm2min(i) = 1d-3
            rm2max(i) = s/2d0
         else
            rm2min(i) = 1d-3           ! minimal/maximal value of mass squared
            rm2max(i) = s/2d0          ! for tan mapping
          endif
      enddo
      
      
      rm2minx1 = ((rm2(1)- 15d0*rmg(1))/sqrt(rm2(1)) )**2
      rm2maxx1 = ((rm2(1)+15d0*rmg(1))/sqrt(rm2(1)) )**2    
      
      rm2minx2 = ((rm2(2)- 15d0*rmg(2))/sqrt(rm2(2)) )**2
      rm2maxx2 = ((rm2(2)+15d0*rmg(2))/sqrt(rm2(2)) )**2          
!      print*, rm2minx

c cutoff for numerical purposes

      m2min = 0.01d0**2

      ptjmin = max(ptmin_def(1),10.0d0)


 60   format ( " machine energy = ", f10.1, "  GeV" )
 70   format ( " integration range for boson ",i3," = ",f10.3,
     1         "  to ",f12.3," GeV" )
 72   format ( " boson has mass ",f8.2," GeV and width ",
     1         f10.4," GeV" )
      if (.not.ldoblha) then
        write(6,*)
        write(6,*) "--------------------PHASE SPACE--------------------"
        write(6,*)
        do i = 1,2
           write(6,72) sqrt(rm2(i)),rmg(i)/sqrt(rm2(i))
           write(6,70) bos(i),sqrt(rm2min(i)),sqrt(rm2max(i))
        enddo
        write(6,*)
        write(6,60) sqrt(s)
        write(6,*)
        write(6,*) "---------------------------------------------------"
        write(6,*)
      endif

      end


C**********************  LIPSVVJ   ****************************************
C
      SUBROUTINE LIPSWZ(N,RD,RN,
     1                   K1,K2,V,Pjet,
     2                   X1,X2,W)
      IMPLICIT NONE
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/kkcoupl.inc"
c
      INTEGER N,NRANDL
      REAL*8 TPI
      PARAMETER (nrandl=(1+3)*3, TPI=2D0*PI )
      real*8 RD(*),RN
      real*8 K1(0:3),K2(0:3),V(0:3,4),X1,X2,W
      real*8 PJET(0:3)

C Phase space generator for production of 2 bosons + 1 jet at hadron colliders
C
C  INPUT:   RD        random number array (between 0 and 1)
C           RN        additional random number which will not be driven by
C                     VEGAS
C
C  OUTPUT:  K1,K2     4-vectors for 2 incoming partons
C           X1,X2     corresponding values of Feynman x
C           Q(0:4,2)  4-vectors for the vector boson resonances V1 and V2,
C                     the 4th component contains the mass**2
C           V(0:3,4)  4-vectors for the 4 decay leptons, assumed massless
C                     V1(Q(*,1)) --> D(*,1)+D(*,2)
C                     V2(Q(*,2)) --> D(*,3)+D(*,4)
C           P(0:3)    4-vector for the jet, assumed massless
C           W         the weight of the phase space point


      REAL*8 RM2,RMG,RM2MIN,RM2MAX, S,YJMIN,YJMAX,PTJMIN
      REAL*8 RM2MINX2, RM2MAXX2, M2MIN, RM2MINX1, RM2MAXX1
      COMMON /BLIPSIWZ/ RM2(2),RMG(2),RM2MIN(2),RM2MAX(2),
     1                RM2MINX1, RM2MAXX1, M2MIN,RM2MINX2, RM2MAXX2,
     2                S,YJMIN,YJMAX,PTJMIN
C
C   RM2            resonance mass squared for resonance
C   RMG            resonance total width * mass for resonance
C   RM2_MIN;_MAX   min & max values of the inv. mass**2 of resonance
C   YJMIN;YJMAX    rapidity range of the individual jets:
C                  YJMIN(I) < ABS(Y(I)) < YJMAX(I)
C   PTJMIN         minimal PT values      REAL*8 Y,PT,PHI,PNJ


C  OUTPUT INFORMATION is transferred via the common block /DIBOSLIPSO/


C  local variables

      INTEGER MU


      double precision xxwidth, xxm
      parameter ( xxwidth = 10d1, xxm = 50d1 )
      real*8 test
      real*8 y, tau, q(0:4,2), ql(0:4,2), y2


      double precision dotrr
      external dotrr
      LOGICAL RESONANCE
      EXTERNAL RESONANCE


      logical ldebug
      parameter (ldebug=.false.)


      w=1d0


!        if ( rd(1) .le. 0.05d0 ) then
!           if (.not. Resonance(xxm**2, xxm*xxwidth, rm2min(1), rm2minx1,
!      1                            rd(1)/0.05d0, W, ql(4,1))) return
!           W = W/0.05d0
!	  print*,'1'
!       else if ( rd(1) .le.  0.95d0 ) then
          if (.not. Resonance(rm2(1), rmg(1), rm2min(1), rm2max(1),
     1                    rd(1), W, ql(4,1))) return
!           W = W/0.9d0
! !	  print*,'2'  
!        else
!           if (.not. Resonance(xxm**2, xxm*xxwidth, rm2maxx1, rm2max(1),
!      1                    (rd(1)-0.95d0)/0.05d0, W, ql(4,1))) return
!           W = W/0.05d0
!	  print*,'3'	  
	  
!       endif
!        print*, w

         if ( rd(2) .le. 0.05d0 ) then
           if (.not. Resonance(xxm**2, xxm*xxwidth, rm2min(2), rm2minx2,
     1                             rd(2)/0.05, W, ql(4,2))) return
           W = W/0.05d0

         else if ( rd(2) .le.  0.98d0 ) then
           if (.not. Resonance(rm2(2), rmg(2), rm2minx2, rm2maxx2,
     1                     (rd(2)-0.05)/0.93d0, W, ql(4,2))) return
           W = W/0.93d0
 
        else
           if (.not. Resonance(xxm**2, xxm*xxwidth, rm2maxx2, rm2max(2),
     1                     (rd(2)-0.98d0)/0.02d0, W, ql(4,2))) return
           W = W/0.02d0
 	endif
	  
!        else
!         if ( rd(2) .le. 0.05d0 ) then
!            if (.not. Resonance(xxm**2, xxm*xxwidth, rm2min(2), rm2minx2,
!      1                             rd(2)/0.05d0, W, ql(4,2))) return
!            W = W/0.05d0
! 
!         else if ( rd(2) .le.  0.95d0 ) then
!            if (.not. Resonance(rm2(2), rmg(2), rm2minx2, rm2maxx2,
!      1                     (rd(2)-0.05d0)/0.9d0, W, ql(4,2))) return
!            W = W/0.9d0
!  
!         else
!            if (.not. Resonance(xxm**2, xxm*xxwidth, rm2maxx2, rm2max(2),
!      1                     (rd(2)-0.95d0)/0.05d0, W, ql(4,2))) return
!            W = W/0.05d0       
!        
!        endif
	  
	  
 !      endif
!       print*, sqrt(ql(4,1)), sqrt(ql(4,2))

!	  rm2max(2)=rm2max(1)-ql(4,1)
! 	  print*, rm2max(1)
! 	  print*, rm2max(2)	
! 	  print*, rn
!       if(.not. Resonance(rm2(2), rmg(2), rm2min(2), rm2max(2), rd(2), w, ql(4,2))) return	 
! 	  print*, '1', w
!   	  print*, ql(4,1)
!   	  print*, ql(4,2)

         call TwoToTwo(rd(3),rn,ECM,ql(4,1),ql(4,2),K1,K2,X1,X2,ql(0,1),
     -       ql(0,2),w)     
! !      ELSE ! 2 -> Y
! c Generate y, tau -> x1, x2, shat from RD(1) and RD(2) first **************
!          tau = qsq/s  !rd(1)
!          W = W / s
!          y_border = log(1/tau)
!          y = y_border*( rn - 0.5d0 )
!          W = W * y_border !(prop density)^-1 of y
! 	 
! !	  print*, '2', w
!          x1 = sqrt(tau) * exp(y)
!          x2 = sqrt(tau) * exp(-y)
! 
! c generate incoming momenta in lab frame out of x1 x2
!          k1(0) = 0.5d0 * x1 * ecm
!          k1(1) = 0d0
!          k1(2) = 0d0
!          k1(3) = 0.5d0 * x1 * ecm
! 
!          k2(0) = 0.5d0 * x2 * ecm
!          k2(1) = 0d0
!          k2(2) = 0d0
!          k2(3) = - 0.5d0 * x2 * ecm
! 
!


!      call TwoBodyDecay(rd(4), rd(5), K(0), K(4),ql(4,1),ql(4,2),ql(0,1),ql(0,2),W)

!	  print*, '3', w

c generate phase space for decay of V1 -> l1 + l2 out of rd(10), rd(11)
      call TwoBodyDecay0(rd(6),rd(7),ql(0,1),ql(4,1),
     1      v(0,1),v(0,2),W)
!	  print*, 'qz', sqrt(ql(0,1)**2-ql(1,1)**2-ql(2,1)**2-ql(3,1)**2)
c generate phase space for decay of V2 -> l3 + l4 out of rd(12), rd(13)
      call TwoBodyDecay0(rd(8),rd(9),ql(0,2),ql(4,2),
     1              v(0,3),v(0,4),W)
!	  print*, '5', w
!       if (.not. TwoBodyDecay(rd(3), rd(4), K(0), K(4),
!      1     ql(0,1),ql(0,2),W)) return
! 
! 
! c generate phase space for decay of V1 -> l1 + l2 out of rd(10), rd(11)
!       if (.not. TwoBodyDecay0(rd(5),rd(6),q(0,1),q(4,1),
!      1      v(0,1),v(0,2),W)) return
! 
! c generate phase space for decay of V2 -> l3 + l4 out of rd(12), rd(13)
!       if (.not. TwoBodyDecay0(rd(7),rd(8),q(0,2),q(4,2),
!      1              v(0,3),v(0,4),W)) return
     
     
!         W = W * 1.d0/2d0/K(4)        ! Flux
!         W = W * TPI                  ! (2*pi)^(4-3) from delta- and ps-int
!         W = W * 3.89379304d11        ! conversion into fb
!          	  print*, '6', w
! 		  stop


! 
!       if (ldebug) print*," enter LIPSVV "
! 
!       W = 1d0             !init of weightfactor
! 
! c---------------------------------------------------------------------
! c pX^2 of intermediate particle X -> V1 + V2 out of rd(10)
! 
!       if (N.eq.1) then
! 
! c---------- Breit-Wigner-distributed-X -----------------
! 
!       xminx = atan( (rm2minx-xxm**2)/xxm/xxwidth )
!       xmaxx = atan( (rm2maxx-xxm**2)/xxm/xxwidth )
!       delxx = xmaxx-xminx
! 
!       xx = xminx + rd(nrandl-2) * delxx
! 
!       tanxx = tan(xx)
! 
!       pX(4) = xxm**2 + xxwidth * xxm * tanxx
! 
! c  ... this factor to the weight W corresponds to dq**2/2pi
! 
!       W = W * xxwidth * xxm * (1 + tanxx**2 )*delxx/TPI
! 
! c--------- plain X ----------------------------------------
! 
! c         pX(4) = rm2minX + rd(nrandl-2)*(rm2maxX - rm2minX)
! c         rm2max(1) = pX(4)
! c         rm2max(2) = pX(4)
! 
! c         W = W * (rm2maxX-rm2minX)/TPI
! 
!        if (pX(4).ge.s) then
! 
!          print*, "cccccccccc WARNING ccccccccccccccc"
!          print*, "c                                c"
!          print*, "c         wrong X-mass           c"
!          print*, "c                                c"
!          print*, "cccccccccccccccccccccccccccccccccc"
! 
!          pX(4) = s - 1d-3
!       endif
! 
!       endif
! 
! c end of pX^2 generation of X
! c--------------------------------------------------------------------
! c q^2 of V1 is generated via tan mapping out of rd(1)
! 
!       xmin1 = atan((rm2min(1)-rm2(1))/rmg(1))
!       if (N.eq.1) rm2max(1) = pX(4)
!       xmax1 = atan((rm2max(1)-rm2(1))/rmg(1))
!       delx1 = xmax1-xmin1
! 
!       xval1 = xmin1 + rd(1) * delx1
! 
!       tanx1 = tan(xval1)
! 
!       q(4,1) = rm2(1)+rmg(1)*tanx1
! 
!       if (N.eq.1.and.(q(4,1).ge.pX(4))) then
! 
!          print*, "cccccccccc WARNING ccccccccccccccc"
!          print*, "c                                c"
!          print*, "c         tan-mapping            c"
!          print*, "c   leads to: W-mass > X-mass    c"
!          print*, "c                                c"
!          print*, "cccccccccccccccccccccccccccccccccc"
!          print*," "
!          print*, " W-mass = ",sqrt(q(4,1))
!          print*, " X-mass = ",sqrt(pX(4))    
!          print*," "
! 
!          q(4,1) = pX(4) - 1d-3
!       endif
! 
! c  ... this factor to the weight W corresponds to dq**2/2pi
! 
!       W = W * rmg(1)*(1 + tanx1**2 )*delx1/TPI
! 
! c end of q^2 generation of V1
! c-----------------------------------------------------------------------
! c q^2 of V2 is generated via tan mapping out of rd(2)
! 
!       xmin2 = atan( (rm2min(2)-rm2(2)) /rmg(2) )
!       if (N.eq.1) rm2max(2) = pX(4)-q(4,1)
!       xmax2 = atan( (rm2max(2)-rm2(2)) /rmg(2) )
!       delx2 = xmax2-xmin2
! 
!       xval2 = xmin2 + rd(2) * delx2
! 
!       tanx2 = tan(xval2)
! 
!       q(4,2)  = rm2(2) + rmg(2)*tanx2
! 
!       if (N.eq.1.and.q(4,2).ge.pX(4)-q(4,1)) then
! 
!          print*, "cccccccccc WARNING ccccccccccccccc"
!          print*, "c                                c"
!          print*, "c         tan-mapping            c"
!          print*, "c leads to: W-mass > (X-W)-mass  c"
!          print*, "c                                c"
!          print*, "cccccccccccccccccccccccccccccccccc"
!          print*," "
!          print*, " W-mass = ",sqrt(q(4,1))
!          print*, " X-mass = ",sqrt(pX(4))    
!          print*," "
! 
!          q(4,2) = pX(4) - 1d-3
!       endif
! 
! c  ... this factor to the weight W corresponds to dq**2/2pi
! 
!       W = W * rmg(2)*( 1 + tanx2**2 )*delx2 /TPI
! 
! c end of q^2 generation of V2
! c------------------------------------------------------------------------
! c set ql^2 for particle 1 and 2 and pTmin
! 
!       if (N.eq.1) then
!          ql(4,1) = 0.0d0
!          ql(4,2) = pX(4)
!          pTmin = 0.0d0
!       else if (N.eq.0) then
!          ql(4,1) = q(4,1)
!          ql(4,2) = q(4,2)
!          pTmin = 0.0d0
!       endif
! 
! c end of ql^2
! c------------------------------------------------------------------------
! c generate y, p_T, y*, phi from plain distributed rd(3),rd(4),rd(5), rn
! 
!       phi = rn * tpi - tpi/2d0
! 
!       if (rd(3).le.0.2) then
!          pTmax = ptjmin
!          pT = rd(3)*(pTmax - pTmin)/0.2d0 + pTmin
!          W = W/0.2d0
!       else
!          pTmin = 2.0d0/ecm
!          pTmax = 1/ptjmin
!          pT = pTmin + (rd(3)-0.2)/0.8d0*(pTmax-pTmin)
!          pT = 1/pT
!          W = W/0.8d0*pT**2
!       endif
! 
!       ET = sqrt(pT**2+ql(4,1))
!       temp = (s+ql(4,1)-ql(4,2))/(2*ecm*ET)
! 
!       y2_border = log( temp+sqrt(temp**2-1) )
! 
!       y2 = y2_border*2d0*( rd(4) - 0.5d0 )          ! y*
! 
!       sinhy2 = sinh(y2)
!       coshy2 = sqrt(sinhy2**2+1)
!       abs_p = sqrt(pT**2*coshy2**2+ql(4,1)*sinhy2**2)
!       pL = sinhy2*ET
!       tau=((coshy2*ET+sqrt(coshy2**2*ET**2-ql(4,1)+ql(4,2)))**2)/s
! 
!       y_border = log(1/tau)
! 
!       y = y_border*( rd(5) - 0.5d0 )                !y
! 
! c-----------------------------------------------------------------------
! c generate phase space for decay of X -> V1 + V2 out of rd(11), rd(12)
! 
!       if (N.eq.1) then
!          phiX = rd(nrandl-1) * tpi - tpi/2d0
!          costhetaX = 2d0 * rd(nrandl) - 1
!          sinthetaX = sqrt(1 - costhetaX**2)
! 
!          abs_pX = 1/(2*sqrt(pX(4)))
!          abs_pX=abs_pX*sqrt((pX(4)-(q(4,1)+q(4,2)))**2-4*q(4,1)*q(4,2))
! 
!          W = W * abs_pX / tpi / 2d0 / sqrt(pX(4))
!       end if
! 
! c-----------------------------------------------------------------------
! c generate phase space for decay of V1 -> l1 + l2 out of rd(6), rd(7)
! 
!       phi2 = rd(6) * tpi - tpi/2d0
!       costheta2 = 2d0 * rd(7) - 1
! 
!       w = w / tpi / 4d0
! 
! c generate phase space for decay of V2 -> l3 + l4 out of rd(8), rd(9)
! 
!       phi3 = rd(8) * tpi - tpi/2d0
!       costheta3 = 2d0 * rd(9) - 1
! 
!       w = w / tpi / 4d0
! 
! c      if ( w .le. 0d0 ) print*, "dps<=0, dps =",w
! 
! c------------------------------------------------------------------------
! c generate incoming momenta in lab frame + x1 x2
! 
!       x1 = sqrt(tau) * exp(y)
!       x2 = sqrt(tau) * exp(-y)
! 
!       k1(0) = 0.5d0 * x1 * sqrt(s)
!       k1(1) = 0d0
!       k1(2) = 0d0
!       k1(3) = 0.5d0 * x1 * sqrt(s)
! 
!       k2(0) = 0.5d0 * x2 * sqrt(s)
!       k2(1) = 0d0
!       k2(2) = 0d0
!       k2(3) = - 0.5d0 * x2 * sqrt(s)
! 
!       do mu = 0, 3
!          K(mu) = k1(mu)+k2(mu)
!       enddo
! 
!       if (N.eq.1) then
! 
! c Momenta of Jet
! 
!          pjet(0) = sqrt( pT**2 + pL**2 + ql(4,1) )
!          pjet(1) = pT * cos(phi)
!          pjet(2) = pT * sin(phi)
!          pjet(3) = pL
! 
!          if (2*min( pjet(0)*k1(0)-pjet(3)*k1(3),
!      &               pjet(0)*k2(0)-pjet(3)*k2(3) ).lt.m2min) then
!             W = 0
!             return
!          endif
! 
! c Momenta of particle X
! 
!          pX(0) = sqrt( pT**2 + pL**2 + ql(4,2) )
!          pX(1) = -pT * cos(phi)
!          pX(2) = -pT * sin(phi)
!          pX(3) = -pL
! 
! c boost in lab frame
! 
!          call boostn( pjet, K, pjet )
!          call boostn( pX, K, pX )
! 
! c Momenta of V1
! 
!          q(0,1) = sqrt( (abs_pX)**2 + q(4,1) )
!          q(1,1) = abs_pX * sinthetaX * cos(phiX)
!          q(2,1) = abs_pX * sinthetaX * sin(phiX)
!          q(3,1) = abs_pX * costhetaX
! 
!          call boostn(q(0,1), pX, q(0,1))
! 
! c determine momenta of V2 from momentum conservation
! 
!          do mu = 0, 3
!             q(mu,2) = pX(mu) - q(mu,1)
!          enddo
! 
!       else if (N.eq.0) then
!             
!          pjet(0) = 0d0
!          pjet(1) = 0d0
!          pjet(2) = 0d0
!          pjet(3) = 0d0
! 
! c Momenta of V1
! 
!          q(0,1) = sqrt( pT**2 + pL**2 + q(4,1) )
!          q(1,1) = pT * cos(phi)
!          q(2,1) = pT * sin(phi)
!          q(3,1) = pL
! 
! 
! c Momenta of V2
! 
!          q(0,2) = sqrt( pT**2 + pL**2 + q(4,2) )
!          q(1,2) = - pT * cos(phi)
!          q(2,2) = - pT * sin(phi)
!          q(3,2) = - pL
! 
! c boost in lab frame
! 
!          call boostn( q(0,1), K, q(0,1) )
!          call boostn( q(0,2), K, q(0,2) )
! 
!       end if
! 
! c Momenta of l1 and l2
! 
!       abs_p2 = 1d0/2d0 * sqrt(q(4,1))
!       sintheta2 = sqrt(1 - costheta2**2)
! 
!       v(0,1) = abs_p2
!       v(1,1) = abs_p2 * sintheta2 * cos(phi2)
!       v(2,1) = abs_p2 * sintheta2 * sin(phi2)
!       v(3,1) = abs_p2 * costheta2
! 
! c boost ve-vector in lab-frame (the frame of W+ [q(0,1)] )
! 
!       call boostn( v(0,1), q(0,1), v(0,1) )
! 
! c determine e+ momenta [v(mu,2)] from momentum conservation
! c of W+ -> ve e+
! 
!       do mu = 0, 3
!          v(mu,2) = q(mu,1) - v(mu,1)
!       enddo
! 
! c Momenta of l3 and l4
! 
!       abs_p3 = 1d0/2d0 * sqrt(q(4,2))
!       sintheta3 = sqrt(1 - costheta3**2)
! 
!       v(0,3) = abs_p3
!       v(1,3) = abs_p3 * sintheta3 * cos(phi3)
!       v(2,3) = abs_p3 * sintheta3 * sin(phi3)
!       v(3,3) = abs_p3 * costheta3
! 
! c boost mu- -vector in lab-frame (the frame of W- [q(0,2)] )]
! 
!       call boostn( v(0,3), q(0,2), v(0,3) )
! 
! c determine vm~ momenta [v(mu,4)] from momentum conservation
! c of W- -> mu- vm~
! 
!       do mu = 0,3
!          v(mu,4) = q(mu,2) - v(mu,3)
!       enddo
! 
! c end of momenta generation
! c-----------------------------------------------------------------------
! c additional phase-space-factors, probability densities,...
! 
!       w = w * 1d0/tpi/4d0/s
!       w = w * pT/(sqrt(tau*s))
!       w = w/(sqrt(coshy2**2*ET**2-ql(4,1)+ql(4,2)))
!       w = w * (pTmax - pTmin)      !(prop density)^-1 of pT
!       w = w * y_border             !(prop density)^-1 of y
!       w = w * y2_border * 2d0      !(prop density)^-1 of y2
!       w = w * 3.89379304d11         !conversion into fb
! 
! c end -----------------------------------------------------
! 
      if (ldebug) then
         print*, "--------- new event -----------"
         print*, ""
         print*, "random numbers"
         print*, rn
         print*, rd(1), rd(2), rd(3)
         print*, rd(4), rd(5), rd(6)
         print*, rd(7), rd(8), rd(9)
         print*, rd(10), rd(11), rd(12)
         print*, "y     = ", y
         print*, "tau   = ", tau
         print*, "x1    = ", x1
         print*, "x2    = ", x2
         print*, "y2    = ", y2
         print*, "wgt   = ", w
         print*, "---------------------------"

         print*, "--- momenta conservation ? ---"
         do mu=0,3
            test = k1(mu)+k2(mu)-q(mu,1)-q(mu,2)!-pjet(mu)
            if (test.gt.1E-9) print*, test
            test = q(mu,1)-v(mu,1)-v(mu,2)
            if (test.gt.1E-9) print*, test
            test = q(mu,2)-v(mu,3)-v(mu,4)
             if (test.gt.1E-9) print*, test
            test = k1(mu)+k2(mu)-v(mu,1)-v(mu,2)!-pjet(mu)
     &               -v(mu,3)-v(mu,4)
             if (test.gt.1E-9) print*, test
         enddo
         print*,"---tau*s = (q1 + q2)^2 ?---"
         test = tau*s-dotrr(q(0,1),q(0,1))-dotrr(q(0,2),q(0,2))
     &                -2*dotrr(q(0,1),q(0,2))
         print*, test
         print*, "--- mass equal to zero ? ---"
         print*,"p1^2 = ",dotrr(k1(0),k1(0))
         print*,"p2^2 = ",dotrr(k2(0),k2(0))
         print*,"pjet^2 = ",dotrr(pjet(0),pjet(0))
         print*,"v1^2 = ",dotrr(v(0,1),v(0,1))
         print*,"v2^2 = ",dotrr(v(0,2),v(0,2))
         print*,"v3^2 = ",dotrr(v(0,3),v(0,3))
         print*,"v4^2 = ",dotrr(v(0,4),v(0,4))
         print*, "-----------------------"

      endif

      return
      end subroutine

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                           c
c     more genereral useful functions, added by             c
c     Manuel Baehr: mb@particle.uni-karlsruhe.de            c
c     2005 March 23                                         c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c      subroutine boostl(y, p)

c      implicit none

c      double precision y, p(0:3), tmp(0:3)
c      double precision coshy, sinhy

c      sinhy = sinh(y)
c      coshy = sqrt(1+sinhy**2)

c      tmp(0) = coshy*p(0) - sinhy*p(3)
c      tmp(3) = coshy*p(3) - sinhy*p(0)

c      p(0) = tmp(0)
c      p(3) = tmp(3)

c      end



