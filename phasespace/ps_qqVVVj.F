c
c phase space generator for pp --> V1 + V2 + V3 + 1 Jet
c                                  V1 --> l1 l2, V2 --> l3 l4, V3 --> l5 l6
c
c       Vera Hankele, <vera@particle.uni-karlsruhe.de>
c       Initial version: June 2006
c       Last modified: September 2012

      SUBROUTINE DEF_PS_QQVVVj(bos1,bos2,bos3)
          use globalvars, only: ldoblha

      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/cuts.inc"
      INTEGER bos(3), i, bos1, bos2, bos3
      common /pswwzbos/bos

      DOUBLE PRECISION RM2, RMG, RM2MIN, RM2MAX, S, PTJMIN
      DOUBLE PRECISION RM2MINX, RM2MAXX, M2MIN, RM2MINY, RM2MAXY
      DOUBLE PRECISION RM2MAX1, RM2MAX2
      COMMON /TRIBOSBLIPSI/ RM2(3), RMG(3), RM2MIN(3), RM2MAX(3),
     1                RM2MINX, RM2MAXX, RM2MINY, RM2MAXY, M2MIN,
     2                S, PTJMIN, RM2MAX1, RM2MAX2

      DOUBLE PRECISION xxwidth, xxm
      common /psVVVx/ xxwidth, xxm

      DOUBLE PRECISION CLR, XM2, XMG, B, V, A
      COMMON /BKOPOU/ CLR(4,5,-1:1), XM2(6), XMG(6),
     1                B(6,6,6), V(4,5), A(4,5)

      double precision M2_i(max_res),MG_i(max_res),F_i(max_res)
      integer N_res

      s = ecm**2

      bos(1) = bos1        ! W+ = 3
      bos(2) = bos2        ! W- = 4
      bos(3) = bos3        ! Z  = 2

      if ( (bos1.ne.bos2) .and. abs(bos1-bos2).ne.1 .and. (bos1+bos2).ne.7 ) then
            print*, "Error:"
            print*, "This phase space is designed only for bos1 and bos2 being the same bosons!"
            stop
      endif

      do i = 1,3
         rm2(i) = xm2(bos(i))        ! mass squared
         rmg(i) = xmg(bos(i))        ! mass * decay width
         rm2min(i) = 1d-3            ! minimal/maximal value of mass squared
         rm2max(i) = s/2d0           ! for tan mapping
      enddo

      ! for intermediate particle Y in NLO
      rm2miny = 1.0d-3
      rm2maxy = s/2d0

      ! for intermediate particle X -> V1 V2
      rm2minx = 1d-3
      rm2maxx = s/2d0   
!      xxm = (sqrt(xm2(bos1)) + sqrt(xm2(bos2)) + sqrt(xm2(bos3)))*0.75
      xxm = 180d0
      xxwidth = 300d0

      ! select relevant resonances for VV invariant mass
      N_res = 0
      if (xm2(6) .lt. 600d0**2) then
        N_res = N_res + 1
        M2_i(N_res) = xm2(6)
        MG_i(N_res) = xmg(6)
        if (bos1.eq.2) then
          F_i(N_res)  = 0.4
        else
          F_i(N_res)  = 0.6
        endif
      endif
      if (bos1 .eq. 2) then         ! Z -> 4l
        N_res = N_res + 1
        M2_i(N_res) = xm2(2)
        MG_i(N_res) = xmg(2)
        F_i(N_res)  = 0.25
      endif
      call MultiResonancesInit(1, xxm**2, xxm*xxwidth, 1d0, M2_i, MG_i, F_i, N_res)


c cutoff for numerical purposes

      m2min = 0.01d0**2
      ptjmin = max(ptmin_def(1),10d0)

 60   format ( " machine energy = ", f10.1, "  GeV" )
 70   format ( " integration range for boson ",i3," = ",f10.3,
     1         "  to ",f12.3," GeV" )
 72   format ( " boson has mass ",f8.2," GeV and width ",
     1         f10.4," GeV" )
      if (.not.ldoblha) then
        write(6,*)
        write(6,*) "--------------------PHASE SPACE--------------------"
        write(6,*)
        do i = 1,3
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

C**********************  LIPSVVVJ   ****************************************
C
      SUBROUTINE LIPSVVVJ(N,RD,RN,
     1                   K1,K2,V,Pjet,
     2                   X1,X2,W, NW)

      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
c
      INTEGER N,NRANDL, NW
      DOUBLE PRECISION TPI
      PARAMETER (nrandl=18, TPI=2D0*PI )
      DOUBLE PRECISION RD(nrandl),RN
      DOUBLE PRECISION K1(0:3), K2(0:3), V(0:3,6), X1, X2, W
      DOUBLE PRECISION Pjet(0:3)

C Phase space generator for the production of 3 bosons (+ 1 jet)
C at hadron colliders
C
C  INPUT:   RD         random number array (between 0 and 1)
C           RN         additional random number which will not be driven by
C                      VEGAS
C           NW         NW = 1 corresponds to Narrow Width Approximation,
C                      NW = 0 means Full calculation
C
C  OUTPUT:  K1,K2      4-vectors for 2 incoming partons
C           X1,X2      corresponding values of Feynman x
C           V(0:3,6)   4-vectors for the 4 decay leptons, assumed massless
C                      Q(*,1) = V(*,1) + V(*,2)
C                      Q(*,2) = V(*,3) + V(*,4)
C                      Q(*,3) = V(*,5) + V(*,6)
C           Pjet(0:3)  4-vector for the jet, assumed massless
C           W          the weight of the phase space point
C


      INTEGER bos(3)
      common /pswwzbos/bos

      DOUBLE PRECISION RM2, RMG, RM2MIN, RM2MAX, S, PTC
      DOUBLE PRECISION RM2MINX, RM2MAXX, RM2MINY, RM2MAXY, M2MIN
      DOUBLE PRECISION RM2MAX1, RM2MAX2
      COMMON /TRIBOSBLIPSI/ RM2(3), RMG(3), RM2MIN(3), RM2MAX(3),
     1                RM2MINX, RM2MAXX, RM2MINY, RM2MAXY, M2MIN,
     2                S, PTC, RM2MAX1, RM2MAX2
C
C   RM2            resonance mass squared
C   RMG            resonance total width * mass
C   RM2_MIN;_MAX   min & max values of the inv. mass**2 of resonance
C   PTJMIN         minimal PT value of the jet

      DOUBLE PRECISION CLR, XM2, XMG, B, VE, A
      COMMON /BKOPOU/ CLR(4,5,-1:1), XM2(6), XMG(6),
     1                B(6,6,6), VE(4,5), A(4,5)

C  local variables

      INTEGER MU

      DOUBLE PRECISION Y
      DOUBLE PRECISION test, q(0:4,3), pX(0:4), pY(0:4)


      DOUBLE PRECISION xxwidth, xxm, xywidth, xym
      DOUBLE PRECISION dotrr

      PARAMETER ( xywidth = 10d1, xym = 45d1 )

      common /psVVVx/ xxwidth, xxm

      EXTERNAL dotrr

      LOGICAL ldebug
      PARAMETER (ldebug=.false.)

      LOGICAL NONRESONANCE, RESONANCE, TWOBODYDECAY, RESONANCE2, ResonanceQsq, multiresonance
      EXTERNAL NONRESONANCE, RESONANCE, TWOBODYDECAY, RESONANCE2, ResonanceQsq, multiresonance

c symmetrisation of q1 and q3
      INTEGER symm
      DOUBLE PRECISION Qtemp(0:4)
      data symm /1/
      save symm    

c--------------------------------------------------------------------

c---------------------------------------------------------------------

      if (ldebug) print*," enter LIPSWWZ "

      W = 1d0             !init of weightfactor


      if (NW.eq.1) then

c q^2_1 in narrow width approximation
         q(4,1) = rm2(1)
         W = W * rmg(1)/2.0d0

c q^2_2 in narrow width approximation
         q(4,2) = rm2(2)
         W = W * rmg(2)/2.0d0   

c q^2_3 in narrow width approximation
         q(4,3) = rm2(3)
         W = W * rmg(3)/2.0d0

c pX^2 of intermediate particle X -> W + W out of rd(1)
         rm2minx = (sqrt(q(4,1)) + sqrt(q(4,2)))**2
         rm2maxx = (sqrt(s) - sqrt(q(4,3)))**2
         if (.not. MultiResonance(1, rd(1), rm2minx, rm2maxx, W, pX(4))) return

c q^2_Y generated out of rd(nrandl-2)
         if (N.eq.1) then
            rm2miny = (sqrt(pX(4)) + sqrt(q(4,3)))**2
            if (.not. Resonance(xym**2, xym*xywidth, rm2miny, rm2maxy,
     1                              rd(nrandl-2), W, pY(4))) return
         endif

      else !NW.ne.1

c q^2_3 generated out of rd(4)
            if (rd(4).lt.0.3d0) then
               W=W/0.3d0
               if (bos(3).eq.2) then
                  if (.not.ResonanceQsq(0d0,rm2min(3), rm2(3)-16d0*rmg(3),
     1                               rd(4)/0.3d0, W, q(4,3))) RETURN 
               else
                  if (.not.Resonance(rm2(3), rmg(3)*2d0, rm2min(3), rm2(3)-16d0*rmg(3),
     1                               rd(4)/0.3d0, W, q(4,3))) RETURN 
               endif
            else
               W=W/0.7d0
               if (.not.Resonance2(rm2(3), rmg(3), rm2(3)-16d0*rmg(3),
     1              rm2max(3), 16d0, (rd(4)-0.3d0)/0.7d0, W, q(4,3))) RETURN 
            endif

c pX^2 of intermediate particle X -> W+ + W- out of rd(1)
         rm2maxx = (sqrt(s) - sqrt(q(4,3)))**2
         if (.not. MultiResonance(1, rd(1), rm2minx, rm2maxx, W, pX(4))) return

c q^2_Y generated out of rd(nrandl-2)
         if (N.eq.1) then
            rm2miny = (sqrt(pX(4)) + sqrt(q(4,3)))**2
            if (.not. Resonance(xym**2, xym*xywidth, rm2miny, rm2maxy,
     1                              rd(nrandl-2), W, pY(4))) return
         endif

c q^2 of W1 is generated via tan mapping out of rd(2)
         rm2max(1) = pX(4)
         if (rd(2).lt.0.3d0) then
            W=W/0.3d0
            if (bos(1).eq.2) then
               if (.not.NonResonance(rm2min(1),
     1               rm2(1)-4d0*rmg(1), rd(2)/0.3d0, W, q(4,1))) RETURN 
            else
               if (.not.Resonance(rm2(1), rmg(1)*2d0, rm2min(1),
     1               rm2(1)-4d0*rmg(1), rd(2)/0.3d0, W, q(4,1))) RETURN 
            endif
         else
            W=W/0.7d0
            if (.not.Resonance2(rm2(1), rmg(1), rm2(1)-4d0*rmg(1), rm2max(1), 16d0,
     1                 (rd(2)-0.3d0)/0.7d0, W, q(4,1))) RETURN 
         endif

c q^2 of W2 is generated via tan mapping out of rd(3)
         rm2max(2) = (sqrt(pX(4)) - sqrt(q(4,1)))**2  ! pX(4) > q(4,1) is assured!
         if (rm2min(1).eq.rm2min(2)) then             ! m_W1 > m_W2
            rm2max(2) = min( rm2max(2), q(4,1) )
            w = w*2d0
         endif
         if (rd(3).lt.0.5d0) then
            W=W/0.5d0
            if (.not.NonResonance(rm2min(2), min(rm2(2)-16d0*rmg(2), rm2max(2)),
     1                               rd(3)/0.5d0, W, q(4,2))) RETURN 
         else
            W=W/0.5d0
            if (.not.Resonance(rm2(2), rmg(2), rm2(2)-16d0*rmg(2),
     1                 rm2max(2), (rd(3)-0.5d0)/0.5d0, W, q(4,2))) RETURN 
         endif
      endif !NW

c-------------------------------------------------------------------
c-------------------------------------------------------------------

      if (N.eq.1) then
c generate 2 -> 2 process for k1 k2 -> Pjet q_Y
         call TwoToJetPlusX(rd(5), rn, ecm, pTc, pY(4),
     1              k1(0), k2(0), x1, x2, pY(0), pJet(0), W)

         if (w.eq.0) return

c generate phase space for decay of Y -> X + Z out of rd(17), rd(18)
         if (.not. TwoBodyDecay(rd(nrandl), rd(nrandl-1), pY(0), pY(4),
     1                 pX(4), q(4,3), pX(0), q(0,3), W)) return
      else if (N.eq.0) then
c generate 2 -> 2 process for k1 k2 -> q_Z q_X
         call TwoToTwo(rd(5), rn, ecm, pX(4), q(4,3),
     1              k1(0), k2(0), x1, x2, pX(0), q(0,3), W)
         pjet(0) = 0.0d0
         pjet(1) = 0.0d0
         pjet(2) = 0.0d0
         pjet(3) = 0.0d0
      endif


      if (W .eq. 0D0) return


c generate phase space for decay of X -> W+ + W- out of rd(8), rd(9)
      if (.not. TwoBodyDecay(rd(9), rd(8), pX(0), pX(4), q(4,1), q(4,2),
     1                                    q(0,1),q(0,2),W)) return

c symmetrisation of W1 and W2
      symm = mod(symm,2) + 1
      if (symm.eq.1) then !interchange W1 <-> W2
         do mu = 0,4
            Qtemp(mu) = q(mu,1)
            q(mu,1)   = q(mu,2)
            q(mu,2)   = Qtemp(mu)
         enddo
      endif

c generate phase space for decay of W+ -> l1 + l2 out of rd(10), rd(11)
      call TwoBodyDecay0(rd(11),rd(10),q(0,1),q(4,1),v(0,1),v(0,2),W)

c generate phase space for decay of W- -> l3 + l4 out of rd(12), rd(13)
      call TwoBodyDecay0(rd(13),rd(12),q(0,2),q(4,2),v(0,3),v(0,4),W)

c generate phase space for decay of Z -> l5 + l6 out of rd(14), rd(15)
      call TwoBodyDecay0(rd(15),rd(14),q(0,3),q(4,3),v(0,5),v(0,6),W)

c ----------------------------------------------------

      if (ldebug) then
         print*, "--------- new event -----------"
         print*, ""
         print*, "random numbers"
         print*, rn
         print*, rd(1), rd(2), rd(3)
         print*, rd(4), rd(5), rd(6)
         print*, rd(7), rd(8), rd(9)
         print*, rd(10), rd(11), rd(12)
         print*, rd(13), rd(14), rd(15)
         if (N.eq.1) print*, rd(16), rd(17), rd(18)
         print*, "y     = ", y
         print*, "x1    = ", x1
         print*, "x2    = ", x2
         print*, "wgt   = ", w
         print*, "---------------------------"

         print*, "--- momenta conservation ? ---"
         do mu=0,3
            test = k1(mu)+k2(mu)-Pjet(mu)-q(mu,1)-q(mu,2)-q(mu,3)
            if (test.gt.1E-9) print*, test
            test = q(mu,1)-v(mu,1)-v(mu,2)
            if (test.gt.1E-9) print*, test
            test = q(mu,2)-v(mu,3)-v(mu,4)
             if (test.gt.1E-9) print*, test
            test = q(mu,3)-v(mu,5)-v(mu,6)
             if (test.gt.1E-9) print*, test
            test = k1(mu)+k2(mu)-Pjet(mu)-v(mu,1)-v(mu,2)
     &             -v(mu,3)-v(mu,4)-v(mu,5)-v(mu,6)
             if (test.gt.1E-9) print*, test
         enddo
         print*, "--- mass equal to zero ? ---"
         print*,"p1^2 = ",dotrr(k1(0),k1(0))
         print*,"p2^2 = ",dotrr(k2(0),k2(0))
         print*,"Pjet^2 = ",dotrr(Pjet(0),Pjet(0))
         print*,"v1^2 = ",dotrr(v(0,1),v(0,1))
         print*,"v2^2 = ",dotrr(v(0,2),v(0,2))
         print*,"v3^2 = ",dotrr(v(0,3),v(0,3))
         print*,"v4^2 = ",dotrr(v(0,4),v(0,4))
         print*,"v5^2 = ",dotrr(v(0,5),v(0,5))
         print*,"v6^2 = ",dotrr(v(0,6),v(0,6))
         print*, "-----------------------"

      endif

      return
      end

